# adapted from http://wiki.cython.org/PackageHierarchy
from __future__ import print_function

import sys, os, shutil, site
import multiprocessing
import subprocess as sb
import tempfile as tmp

import mpi4py

from distutils.core import setup
from distutils.extension import Extension
from distutils import log
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# Use parallel compilation on this number of cores.
nthreads = int(os.getenv('COMPILE_NTHREADS', multiprocessing.cpu_count() ))
isRoot = os.geteuid() == 0  # Do we have root privileges?

class inTempFolder:
    """Context manager for working in temporary folder.

    Creates a temporary folder enters in it and removes it when done.
    """
    def __enter__(self):
        self.originalDir = os.getcwd()
        self.tmpdir = tmp.mkdtemp()
        os.chdir(self.tmpdir)
        return self.tmpdir, self.originalDir
    def __exit__(self, type, value, traceback):
        shutil.rmtree(self.tmpdir)
        os.chdir(self.originalDir)


def installJDFTx(isRoot=False, enableGPU=False):
    """Search and return the path to libjdfx.so and compile jdftx if not found.

    Check user folder for jdftx/libjdftx.so. Unless jdftx library is found
    compile the source code, copy the library to the user directory and
    return the path to the library. If GPU is enabled, look for libjdftx_gpu.so
    too.
    """
    # is there a valid installation in user folders:
    if os.path.exists(os.path.join(site.USER_BASE, "jdftx/libjdftx.so")):
        if (not enableGPU) or \
           os.path.exists(os.path.join(site.USER_BASE, "jdftx/libjdftx_gpu.so")):
           return os.path.join(site.USER_BASE, "jdftx")

    with inTempFolder() as (jdftxDir, pythonJDFTxDir):
        log.info("JDFTx Compilation:")
        log.info("Running cmake...")
        jdftxCodeDir = os.path.join(pythonJDFTxDir, "jdftx")
        if enableGPU:
            sb.check_call(["cmake", "-D", "EnableCUDA=yes",
                           jdftxCodeDir])
        else:
            sb.check_call(["cmake", "-D", jdftxCodeDir])
        log.info("Running make. This takes a few minutes.")
        sb.check_call(["make", "-j%d" % nthreads])
        if isRoot:
            jdftxLibDir = "/usr/local/jdftx"
        else:
            jdftxLibDir = os.path.join(site.USER_BASE, "jdftx")
        if not os.path.exists(jdftxLibDir):
            os.mkdir(jdftxLibDir)

        shutil.move("libjdftx.so", jdftxLibDir)
        if enableGPU:
            shutil.move("libjdftx_gpu.so", jdftxLibDir)
        if not os.path.exists(os.path.join(site.USER_BASE, "jdftx/pseudopotentials")):
            shutil.move("pseudopotentials", jdftxLibDir)
        try:
            os.symlink("/usr/local/jdftx/libjdftx.so",
                       "/usr/lib/libjdftx.so")
            if enableGPU:
                os.symlink("/usr/local/jdftx/libjdftx_gpu.so",
                           "/usr/lib/libjdftx_gpu.so")
            return ""
        except OSError:
            return jdftxLibDir

def make_extension(ext_name, ext_libraries=(), is_directory=False, enableGPU=0):
    try:
        sb.check_call(["ld", "-ljdftx"], stderr=open("/dev/null"))
        if enableGPU:
            sb.check_call(["ld", "-ljdftx_gpu"])
        jdftxLibDir = ""
    except sb.CalledProcessError:
        jdftxLibDir = installJDFTx(isRoot, enableGPU)

    jdftxIncDirs = ["jdftx", ".", mpi4py.get_include()]

    if enableGPU: # adds only the default directory for now
        jdftxIncDirs.append("/usr/local/cuda/include")

    ext_path = ext_name
    if is_directory:
        ext_path += ".__init__"
    return Extension(
        ext_name,
        [ext_path.replace(".", os.path.sep) + ".pyx"],
        include_dirs=jdftxIncDirs,
        language="c++",
        libraries=ext_libraries,
        library_dirs=[jdftxLibDir],
        runtime_library_dirs=[jdftxLibDir],
        extra_compile_args=['-std=c++0x', '-O3', '-DMPI_ENABLED'] +
                            ['-DGPU_ENABLED'] * enableGPU,
        #depends=["jdftx/libjdftx.so"],
    )

def writeTargetToSourceCode(target, baseFName="JDFTxCalcBase.pyx"):
    targetFName = "JDFTxCalc"+target+".pyx"
    with open(baseFName) as baseFile, open(targetFName, 'w') as t:
        t.write("# Generated by the setup script from " + baseFName + "\n")

        # Do not copy the first comment block that is specific to JDFTxCalcBase.pyx.
        while True:
            if baseFile.readline()[0] != '#': break
        t.write(baseFile.read().format(TARGET = target))

def main():
    enableGPU = (len(sys.argv) >= 2 and "--GPU" in sys.argv)
    if enableGPU: sys.argv.pop(sys.argv.index("--GPU"))

    extensions = []
    if enableGPU:
        writeTargetToSourceCode(target = 'GPU')
        extensions.append(
            make_extension("JDFTxCalcGPU", ["jdftx_gpu"], enableGPU=True,
                            ),
    	)

    writeTargetToSourceCode(target = 'CPU')
    extensions.append(
        make_extension("JDFTxCalcCPU", ["jdftx"]),
    )

    for e in extensions:
        e.cython_directives = {"boundscheck": False,
                               "wraparound": False,
                               "infer_types": True}

    mpiCompilers = mpi4py.get_config()
    os.environ['CC'] = mpiCompilers['mpicc']
    os.environ['CXX'] = mpiCompilers['mpicxx']

    pyVersion = sys.version_info[0]

    extensions = cythonize(extensions, nthreads=nthreads,
                            compiler_directives = {'language_level': pyVersion})

    setup(**{
        "name": "pythonJDFTx",
        # "packages": [
        #     "core",
        #     "electronic",
        #     "includes",
        #     "fluid",
        # ],
        "py_modules":["ElectronicMinimize"],
        "ext_modules": extensions,
        "cmdclass": {'build_ext': build_ext},
    })

if __name__ == "__main__":
    main()

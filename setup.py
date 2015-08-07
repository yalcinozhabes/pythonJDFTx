# adapted from http://wiki.cython.org/PackageHierarchy

import sys
import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# Use parallel compilation on this number of cores.
nthreads = int(os.environ.get('COMPILE_NTHREADS', 0))
nthreads = 2

def make_extension(ext_name, ext_libraries=(), is_directory=False):
    ext_path = ext_name
    if is_directory:
        ext_path += ".__init__"
    return Extension(
        ext_name,
        [ext_path.replace(".", os.path.sep) + ".pyx"],
        include_dirs=(["../../jdftx", "."]),
        language="c++",
        libraries=ext_libraries,
        library_dirs=["../"],
        extra_compile_args=['-std=c++11'],
        #depends=["jdftx/libjdftx.so"],
    )

extensions = [
    make_extension("JDFTCalculator", ["jdftx"]),
]

setup(**{
    "name": "pythonJDFTx",
    "packages": [
        "core",
        "electronic",
        "includes",
        "fluid",
    ],
    "ext_modules": cythonize(extensions, nthreads=nthreads,compiler_directives = {'language_level':3}),
    "cmdclass": {'build_ext': build_ext},
})

pythonJDFTx
===========================
This is a computational quantum chemistry tool, written in [Cython](http://cython.org),
which brings a fast C++ package (`JDFTx`)
to the Python environment where lots of fun stuff happen.

JDFT stands for Joint Density Functional Theory and `JDFTx` is the implementation of this
theory which lives at http://sourceforge.net/p/jdftx/wiki/Home/ . Feel free to visit and
investigate the code and appreciate its great capabilities such as GPU parallelization.

`pythonJDFTx` is a Python wrapper that provides a calculator interface for [ASE
(Atomic Simulation Environment)](https://wiki.fysik.dtu.dk/ase/). ASE makes it possible
to use many other calculators, aside from JDFTx, using the same syntax. This way, a PhD
student such as myself, does not need to read through the documentation of many other
legacy codes, learn their input file syntax and write scripts in order to convert
the parameter files from one to another. Letting ASE do the hard work, one can use not
only Joint Density Functional Theory but many other theories to understand and explain
physical phenomena much easily.

Technical details
===========================

Installation
---------------------------
`pythonJDFTx` needs both ASE and JDFTx installed on your computer. Right now, it has not
been tested in Python 2.7 (works in Python 3.4). Also small changes in `JDFTx` is needed
to make some class variables and methods public and accessible in Python. You can find a
`diff` file that you can patch to JDFTx source tree before compilation.

After compiling `JDFTx` following the steps described in JDFTx wiki, you should have
a library named `libjdftx.so`. You should modify `setup.py` to show the cython compiler
where the library is. I also needed to set an environmental variable for `ld` as
`export LD_LIBRARY_PATH=<path-to-libjdftx.so>`. The setup file should work without any
modification if you cloned this folder in your `build` directory.

    python3 setup.py build_ext --inplace

Pseudopotential Files
---------------------------
Support for Ultrasoft Pseudopotential format has been implemented. Support for
norm conserving pseudo-potentials is under development. `pythonJDFTx` looks for
an environmental variable, `PSEUDOPOT_HOME`, at runtime. Make sure this variable shows
the path to the pseudo-potential files and the files are named all lowercase with
extension `.uspp`.

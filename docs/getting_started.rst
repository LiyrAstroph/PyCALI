.. _getting_started:

***************
Getting Started
***************

.. _installing-docdir:

Requirements
============
pyCALI requires the following third-party packages

* **CMake**: https://cmake.org/
  
  a software for package compilation.


* **Pybind11**: https://github.com/pybind/pybind11
  
  a lightweight header-only library for python wrapper.


* **LAPACKE**: http://performance.netlib.org/lapack/
  
  C/C++ interface to LAPACK.


* **CBLAS**: https://www.netlib.org/blas/

  C/C++ interface to BLAS.


* **GSL**: https://www.gnu.org/software/gsl/
  
  GNU Scientific Library.

Installation
============

pyCALI uses CMake to do building and compilation. 

.. code-block:: bash

  cmake .
  make

This will create an exectuable binary file ``cali`` and a Python module ``pycali``.

If one wants to create only ``cali``, use the command 

.. code-block:: bash 

  cmake .
  make cali 

Similarly, if one wants to create only ``pycali``, use the command 

.. code-block:: bash 

  cmake .
  make pycali 

The installation presumes that LAPACKE and CBLAS are installed in the default paths, namely, for LAPACKE, headers placed 
at /usr/include/lapacke and libraries at /usr/lib or /usr/lib64; for CBLAS, headers placed 
at /usr/include/cblas and libraries at /usr/lib or /usr/lib64.  If this is not the case, use the CMake GUI to 
make editing

.. code-block:: bash 
  
  ccmake .

The triggered GUI generally looks like 

.. code-block:: bash 

  CBLAS_INCLUDE_DIR                /usr/include/cblas
  CBLAS_LIB                        /usr/lib64/libcblas.so
  CMAKE_BUILD_TYPE
  CMAKE_INSTALL_PREFIX             /usr/local
  LAPACKE_INCLUDE_DIR              /usr/include/lapacke
  LAPACKE_LIB                      /usr/lib64/liblapacke.so
  PYBIND11_CPP_STANDARD            -std=c++14
  PYBIND11_PYTHON_VERSION
  pybind11_DIR                     /usr/share/cmake/pybind11


Usage
=====

Either ``cali`` or ``pycali`` can be used to do intercalibrating 

.. code-block:: bash
  
  ./cali param.txt 

in which ``param.txt`` specifies the configurations passed to ``cali``.

For the Python module ``pycali``, a Python script ``example.py`` shows
an example regarding the usage.

A python script ``plot_results.py`` in the subdirtory ``data/`` shows how to plot 
the merged light curves and the posterior distributions of parameters. 

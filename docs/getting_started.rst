.. _getting_started:

***************
Getting Started
***************

.. _installing-docdir:

Requirements
============
pyCALI requires the following third-party packages

* **CMake** 
  
  a software for package compilation, available at https://cmake.org/.


* **Pybind11**
  
  a lightweight header-only library for python wrapper, available at https://github.com/pybind/pybind11.


* **LAPACKE**
  
  C/C++ interface to LAPACK.


* **CBLAS**

  C/C++ interface to BLAS.


* **GSL**
  
  GNU Scientific Library.

Installation
============

pyCALI uses CMake to do building and compilation. 

.. code-block:: bash

  cmake .
  make

This will create an exectuable binary file ``cali`` and a Python module ``pycali``.

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

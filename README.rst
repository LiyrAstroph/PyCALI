pyCALI
======

A Python package for intercalibrating light curves. 

References: `Li, Y.-R., et al. 2014, ApJL, 786, 6 <https://ui.adsabs.harvard.edu/abs/2014ApJ...786L...6L/abstract>`_

Third-party packages
--------------------

pyCALI requires the following third-party packages

* **CMake** 
  
  a software for package compilation, available at https://cmake.org/.

* **Pybind11**
  
  a lightweight header-only library for python wrapper, available at https://github.com/pybind/pybind11.

* **Lapacke**
  
  C/C++ interface to Lapack.

* **Cblas**

  C/C++ interface to Blas.

* **GSL**
  
  GNU Scientific Library.

Installation
------------

pyCALI uses CMake to do building and compilation. 

.. code-block:: bash

  cmake .
  make

This will create an exectuable binary file ``cali`` and a Python module ``pycali``.

Usage
-------

Either ``cali`` or ``pycali`` can be used to do intercalibrating 

.. code-block:: bash
  
  ./cali param.txt 

in which ``param.txt`` specifies the configurations passed to ``cali``.

For the Python module ``pycali``, a Python script ``example.py`` shows
an example regarding the usage.

A python script ``plot_results.py`` in the subdirtory ``data/`` shows how to plot 
the merged light curves and the posterior distributions of parameters. 
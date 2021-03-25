.. _getting_started:

***************
Getting Started
***************

.. _installing-docdir:

Requirements
============
PyCALI requires the following third-party packages. 
(Thanks to Yong-Jie Chen and Yu-Yang Songsheng for help) 

* **CMake**: https://cmake.org/
  
  a software for package compilation.

  In Fedora/Redhat distributions, use the following command to install CMake

  .. code-block:: bash
  
    sudo dnf install cmake cmake-gui
  
  In Debian/Ubuntu distributions, use the command 

  .. code-block:: bash
    
    sudo apt install cmake cmake-curses-gui


* **LAPACKE** (optional): http://performance.netlib.org/lapack/
  
  C/C++ interface to LAPACK.

  In Fedora/Redhat distributions, use the following command to install LAPACKE

  .. code-block:: bash
  
    sudo dnf install lapack lapack-devel
  
  In Debian/Ubuntu distributions, use the command 

  .. code-block:: bash 

    sudo apt install liblapacke-dev
  
  .. note::

    * There is a LAPACK variant written purely in C, called CLAPACK. Do not confuse it with LAPACKE. 
  
    * **If the system doest not have LAPACKE libraries, the code will automatically compile the lapacke source 
      codes (verion 3.9.0) packaged along with PyCALI.** 

* **CBLAS** (optional): https://www.netlib.org/blas/

  C/C++ interface to BLAS.

  In Fedora/Redhat distributions, use the following command to install CBLAS

  .. code-block:: bash
  
    sudo dnf install blas blas-devel

  In Debian/Ubuntu distributions, use the command 

  .. code-block:: bash 

    sudo apt install libblas-dev
  
  .. note::
    
     **If the system doest not have CBLAS libraries, the code will automatically compile the cblas source 
     codes (provided by LPACKE verion 3.9.0) packaged along with PyCALI.**

* **GSL**: https://www.gnu.org/software/gsl/
  
  GNU Scientific Library.

  In Fedora/Redhat distributions, use the following command to install GSL

  .. code-block:: bash
  
    sudo dnf install gsl gsl-devel
  
  In Debian/Ubuntu distributions, use the command 

  .. code-block:: bash 

    sudo apt install libgsl-dev

* **Pybind11**: https://github.com/pybind/pybind11
  
  a lightweight header-only library for python wrapper.

  Use the following command to install Pybind11

  .. code-block:: bash

    pip install pybind11
  
  Refer to `Installing Pybind11 <https://pybind11.readthedocs.io/en/stable/installing.html#>`_ for details.

  .. note::

    For Python provided by anaconda, ``pip install pybind11`` will put configuration file ``pybind11Config.cmake`` into  
    ``<install-dir-of-pybind11>/share/cmake/pybind11``, which can not be found by CMake unless you specify it via
    ``cmake -D pybind11_DIR=<install-dir-of-pybind11>/share/cmake/pybind11`` explicitly when using CMake.
    However, ``pip install "pybind11[global]"`` will put the configuration file into ``<install-dir-of-anaconda>
    share/cmake/pybind11``, which can be found by CMake as long as the path ``<install-dir-of-anaconda>`` is included in the
    $PATH. 

* **cmaketools**: https://pypi.org/project/cmaketools/
  
  An integration of Cmake build system to Python setuptools/distutils.
  Only used for Python wrapper.

  Use the following command to install camketoolds

  .. code-block:: bash

    pip install cmaketools

Installation
============
PyCALI uses CMake to do building and compilation. 

The following installations presume that LAPACKE and CBLAS are installed in the default paths, namely, for LAPACKE, headers placed 
at /usr/include/lapacke and libraries at /usr/lib or /usr/lib64; for CBLAS, headers placed 
at /usr/include/cblas and libraries at /usr/lib or /usr/lib64. (Note that this generally works in Fedora/Redhat distributions.
See below for Ubuntu/Debian distributions.) 

If the above libraries are not installed in the default paths, use the CMake GUI to 
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


.. note::

  * If using **clang** compiler, one may explicitly add **-std=c++11** or something like in **CMakeLists.txt**
    that to support the C++ standards, see https://clang.llvm.org/cxx_status.html.  

  * Debian/Ubuntu science team maintainers have merged the CBLAS ABI into **libblas.so**. 
    Everything one needs from **libcblas.so** are included in **libblas.so**. So for Debian/Ubuntu systems, 
    one shoud refer **CBLAS_LIB** to **libblas.so** instead of **libcblas.so**.
  
  * For Debian/Ubuntu systems, if one insists on using **libcblas.so**,  install **libatlas3-base (/-dev)**, 
    which is the only provider in archives. That **libcblas.so** provided by **libatlas3-base** is quite
    slow in terms of performance if not re-compiled locally. In this case,  the header file **cblas.h**
    (usually in /usr/include/x86_64-linux-gnu/ for amd64 architecture) is indeed a soft link to
    **cblas-atlas.h**. A problem with **cblas_atlas.h** is that it can not be called from C++ program. 
    To amend it, one should modify cblas-atlas.h as the following: 
    add
    
    .. code-block:: C
      
      #ifdef __cplusplus
      extern "C" { /* Assume C declarations for C++ */ 
      #endif /* __cplusplus */ 

    after the first line
    
    .. code-block:: C
      
      #ifndef CBLAS_H

    and add 
    
    .. code-block:: C
      
      #ifdef __cplusplus
      } 
      #endif 

    before the last line

    .. code-block:: C
    
      #endif 
  
  * When installing **pycali**, one may encounter errors like::
    
      fatal error: Python.h: No such file or directory

      #include <Python.h>
    
    This error can be solved by installing the header file of Python, e.g.,

    .. code-block:: Python 

      dnf install python-devel


C/C++ executable binary: cali
-----------------------------

If one wants to create executable binary file ``cali``, use the command 

.. code-block:: bash 

  cmake .
  make cali 

Python module: pycali
---------------------

If one wants to create Python module ``pycali``, use the command 

.. code-block:: bash
  
  cmake .
  python setup.py install --user 

This will install pycali module to a path that can be reconginzed by the Python interpretor.
Usually this path is located at, e.g., .local/lib/python3.9/site-packages. 


Basic Usage
===========

Either ``cali`` or ``pycali`` can be used to do intercalibrating.  ``cali`` is an executable binary file 
and can directly executed in a Linux terminal as

.. code-block:: bash
  
  ./cali param.txt 

in which ``param.txt`` specifies the configurations passed to ``cali``.

For the Python module ``pycali``, a Python script ``example.py`` shows
an example regarding the usage.

.. note::

  A directory "data/" in the present working directory is needed to place ouput files. ``cali`` and ``pycali``
  automatically check whether the directory exists. If not, it will be created.

A python script ``plot_results.py`` in the subdirtory ``data/`` shows how to plot 
the merged light curves and the posterior distributions of parameters. 

The final intercalibrated light curves are output to files with a name by adding a postfix "_cali" 
to the input file name. For example, if your intput file name is "exmaple.txt", the output 
file name is "example.txt_cali".

Please also refer to :ref:`faq` for more details not covered here.
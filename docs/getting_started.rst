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
  
  a lightweight header-only library for python wrapper. **Note that Pybind11 is only needed when installing with Python.
  If only installing the C version, no need to install Pybin11.**

  Use the following command to install Pybind11

  .. code-block:: bash

    pip install pybind11
    sudo dnf install pybind11-devel #(Fedora/Redhat)
    sudo apt install pybind11-dev   #(Debian/Ununtu)
  
  Refer to `Installing Pybind11 <https://pybind11.readthedocs.io/en/stable/installing.html#>`_ for details.

  .. note::

    For Python provided by anaconda, ``pip install pybind11`` will put configuration file ``pybind11Config.cmake`` into  
    ``<install-dir-of-pybind11>/share/cmake/pybind11``, which can not be found by CMake unless you specify it via
    ``cmake -D pybind11_DIR=<install-dir-of-pybind11>/share/cmake/pybind11`` explicitly when using CMake.
    However, ``pip install "pybind11[global]"`` will put the configuration file into ``<install-dir-of-anaconda>
    share/cmake/pybind11``, which can be found by CMake as long as the path ``<install-dir-of-anaconda>`` is included in the
    $PATH. 

Installation with CMake
=======================
This only installs executable binary ``cali``.

A common error occuring frequently is LAPACKE and CBLAS libraries not found. PyCALI also packages the source codes 
of LAPACKE and CBLAS. One can use these source codes if encountering problems with installing LAPACKE and CBLAS.
If so, one usually do no need to edit CMake configurations described below and keep things unchanged.

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

  BLAS_LIB                         /usr/lib64/libblas.so
  CBLAS_INCLUDE_DIR                /usr/include
  CBLAS_LIB                        /usr/lib64/libcblas.so
  CMAKE_BUILD_TYPE
  CMAKE_INSTALL_PREFIX             /home/liyropt/Projects/GIT/PyCALI/dist
  LAPACKE_INCLUDE_DIR              /usr/include
  LAPACKE_LIB                      /usr/lib64/liblapacke.so
  LAPACK_LIB                       /usr/lib64/liblapack.so


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

Then use the command 

.. code-block:: bash 

  cmake .
  make cali 

Installation with Python
========================
This only installs Python version ``pycali``.

Use the command 

.. code-block:: bash
  
  cmake .
  python setup.py install --user 

This will install pycali module to a path that can be reconginzed by the Python interpretor.
Usually this path is located at, e.g., .local/lib/python3.9/site-packages. 

Installation with Makefile
==========================
If your system does not have latest CMake or the installation with Python does not work, you 
may take a try with Makefile. 

First copy the file "Makefile_old" in the package to "Makefile", i.e., 

.. code-block:: bash

  cp Makefile_old Makefile

Then edit some configuration options in Makefile according your system's settings. After that, 
execute the command 

.. code-block:: bash

  make 

This will crate a executable binary ``cali``.

Basic Usage
===========

Either ``cali`` or ``pycali`` can be used to do intercalibrating.  

Usage in a terminal
-------------------
``cali`` is an executable binary file and can directly executed in a Linux terminal as

.. code-block:: bash
  
  ./cali param.txt 

in which ``param.txt`` specifies the configurations passed to ``cali``.

For the Python module ``pycali``, a Python script ``example.py`` shows
an example regarding the usage.

.. note::

  A directory "data/" in the present working directory is needed to place ouput files. ``cali`` and ``pycali``
  automatically check whether the directory exists. If not, it will be created.

Usage in Python 
---------------
A python script ``plot_results.py`` in the subdirtory ``data/`` shows how to plot 
the merged light curves and the posterior distributions of parameters. 

The final intercalibrated light curves are output to files with a name by adding a postfix "_cali" 
to the input file name. For example, if your intput file name is "exmaple.txt", the output 
file name is "example.txt_cali".

Please also refer to :ref:`faq` for more details not covered here.

Format of Input Data files
===========================

``cali`` or ``pycali`` reads input data files with the following format::

  # code1 120     
  7517.0   1.98   0.08
  7534.0   2.06   0.08
  ...
  7719.0   2.03   0.08
  7725.0   1.97   0.08
  7778.0   2.02   0.08
  # code2 45
  7573.0   2.73   0.11
  7584.0   2.73   0.11
  ...
  7644.0   3.45   0.14
  7661.0   3.26   0.13
  # code3 33
  7509.0   1.92   0.08
  7530.0   1.97   0.08
  ...
  7556.0   2.21   0.09
  7614.0   2.31   0.09

In the above file, there are three codes (code1, code2 and code3) with 120, 45 and 33 points respectively. 

``pycali`` provides a function to generate input formatted data file as 

.. code-block:: python

  import pycali 

  pycali.format(fname, data)
  # "fname" is the file name to generate
  # "data" is a python dict that stores the data, in which the keys represent the codes

Besides, ``pycali`` provides functions to convert ASAS-SN and ZTF data

.. code-block:: python 

  import pycali 
  
  ztf = pycali.convert_ztf("ZTF.csv", rebin=True, errlimit=0.079)   
  # rebin:  whether rebin the points within one day
  # errlimit: discard these points with errors larger than this limit
  # return a dict
  
  asassn = pycali.convert_asassn("asas.csv", rebin=True, errlimit=0.079, diffcamera=False)
  # diffcamera: whether treat different cameras as different datasets
  
  data = ztf | asassn  # combine the two dicts
  
  # write to a file named "test.txt"
  # trange is the range of time to use
  pycali.format("test.txt", data, trange=[t1, t2]) 
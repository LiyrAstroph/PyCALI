from glob import glob
from setuptools import setup
import pybind11
from pybind11.setup_helpers import Pybind11Extension
import os, sys
import pkgconfig

# force to use g++
os.environ["CC"] = "g++" #"clang++"
os.environ["CXX"] = "g++" #"clang++"

basedir = os.path.dirname(os.path.abspath(__file__))

def configure_gsl():
  """
  get configuration of gsl
  """
  if pkgconfig.exists('gsl'):
    gslconf = pkgconfig.parse('gsl')
  else:
    raise SystemError("Not found GSL installed. \n" 
                      + "Please install GSL (GNU Scientific Library) first. \n"
                      + "You can install it via your package manager, e.g., \n"
                      + "'apt install libgsl-dev' on Debian/Ubuntu, \n"
                      + "'yum install gsl-devel' on CentOS/RHEL, \n"
                      + "'dnf install gsl-devel' on Fedora, \n"
                      + "'brew install gsl' on macOS. \n"
                      )

  return gslconf

def configure_lapack():
  """
  get configuration of lapack/lapacke
  """
  if pkgconfig.exists('lapacke'):
    lapackconf = pkgconfig.parse('lapacke')
  else:
    raise SystemError("Not found LAPACK installed. \n" 
                      + "Please install LAPACK (Linear Algebra PACKage) first. \n"
                      + "You can install it via your package manager, e.g., \n"
                      + "'apt install liblapacke-dev' on Debian/Ubuntu, \n"
                      + "'yum install lapack-devel' on CentOS/RHEL, \n"
                      + "'dnf install lapack-devel' on Fedora, \n"
                      + "'brew install lapack' on macOS. \n"
                      )

  return lapackconf

gslconf = configure_gsl()
lapackconf = configure_lapack()

# get path to pybind11/pybind11.h
pybind11_include_dir = pybind11.get_include()
pybind11_include_dir = os.path.abspath(os.path.join(pybind11_include_dir, "..", ".."))
if not os.path.exists(pybind11_include_dir):
    raise SystemError("pybind11 include directory not found. Please install pybind11.")

# libraries
libraries = ['m', 'c', 'gsl', 'gslcblas', 'lapack', 'lapacke']
compiler_args = ['-O3', '-ffast-math', '-fcommon', '-fpermissive','-std=c++11']

# if gsl, lapack are not in the standard path, put their paths here.
include_dirs=[basedir] + gslconf['include_dirs'] + lapackconf['include_dirs'] + [pybind11_include_dir]
library_dirs=[basedir] + gslconf['library_dirs'] + lapackconf['library_dirs'] 

# source files
src = glob(os.path.join(basedir, "src/pycali/pycali", "*.cpp")) + glob(os.path.join(basedir, "src/pycali/pycali", "*.c")) \
    + glob(os.path.join(basedir,"src/pycali/cdnest", "*.c"))

# headers
headerfiles = glob(os.path.join(basedir, "src/pycali/pycali", "*.hpp")) + glob(os.path.join(basedir, "src/pycali/pycali", "*.h")) \
            + glob(os.path.join(basedir, "src/pycali/cdnest", "*.h"))

# python module
ext_modules = [
    Pybind11Extension(
        "pycali.pycali",
        sources=src,
        depends=headerfiles,
        extra_compile_args=compiler_args,
        libraries=libraries,
        include_dirs=include_dirs,
        library_dirs=library_dirs
    )
    ]

setup(
    name="pycali",
    version="0.2.3",
    author="Yan-Rong Li",
    author_email="liyanrong@mail.ihep.ac.cn",
    description="A Bayesian approach to intercalibrate light curves",
    url="https://github.com/LiyrAstroph/PyCALI",
    license="MIT License",
    packages=["pycali"],
    package_dir={"":"src"},
    ext_modules=ext_modules,
    setup_requires=["pybind11"],
    install_requires=["pybind11", "numpy", "matplotlib", "corner"],
)

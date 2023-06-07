from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
import os, sys

basedir = os.path.dirname(os.path.abspath(__file__))

# libraries
libraries = ['m', 'c', 'gsl', 'gslcblas', 'lapack', 'lapacke','blas']
compiler_args = ['-O3', '-ffast-math', '-fcommon'] 

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
    )
    ]

setup(
    name="pycali",
    version="0.2.2",
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

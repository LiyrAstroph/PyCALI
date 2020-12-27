from cmaketools import setup
import os, sys

setup(
    name="pycali",
    version="0.0",
    author="Yan-Rong Li",
    author_email="liyanrong@mail.ihep.ac.cn",
    description="A Bayesian approach to intercalibrate light curves",
    url="https://github.com/LiyrAstroph/pyCALI",
    license="MIT License",
    src_dir="src",
    ext_module_hint=r"pybind11_add_module",
    has_package_data=False,
    install_requires=["cmaketools"]
)

"""Build the compiled term loop:  sage -python scripts/setup_cy.py build_ext --inplace"""
from setuptools import Extension, setup
from Cython.Build import cythonize
from sage.env import sage_include_directories

setup(
    name="symfn_cy",
    ext_modules=cythonize(
        [Extension("symfn_cy", ["scripts/symfn_cy.pyx"],
                   include_dirs=sage_include_directories())],
        language_level=3,
    ),
)

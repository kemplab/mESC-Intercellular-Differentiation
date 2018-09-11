try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:

    from distutils.core import setup
    from distutils.extension import Extension
import numpy

from Cython.Build import cythonize
"""
This script is meant to compile the intercellular diffusion module to C-code. 

To compile go to a command console and:
1. Change to this directory (type cd, press space then drag/drop this folder to the console. Press enter)
2. Type out: "python setup.py build_ext --inplace" without the quotation marks. Press enter.

"""
setup(
    ext_modules=cythonize('IntC_FP.pyx'),
    include_dirs=[numpy.get_include()]
    )

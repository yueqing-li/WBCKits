import os
import numpy as np

from setuptools import setup, Extension
from Cython.Build import cythonize

BASEDIR = os.path.dirname(os.path.abspath(__file__))
BASEDIR = os.path.dirname(BASEDIR)
BASEDIR = os.path.dirname(BASEDIR)

ext = Extension(name="pyWBCKits",
                sources=["WBCKits.pyx"],
                include_dirs=["/usr/local/include", np.get_include()],
                # -L 库路径
                library_dirs=["/usr/local/lib", os.path.join(BASEDIR, 'build')],
                # -l 库名称
                libraries=["wbckits"]
                )
setup(
    name='pyWBCKits',
    ext_modules=cythonize(ext, language_level=3)
)
import setuptools
from distutils.core import setup,Extension

from Cython.Distutils import build_ext
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

with open ("README.md","r") as f:
    long_description = f.read()

setup(
    name = 'TAPpy',
    version = '0.0.1',
    author = 'Riley Vickers',
    author_email = 'rvickers93@gmail.com',
    description = 'Trajectory Analysis of Particles in python.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/rileyvickers/TAPpy',
    packages=setuptools.find_packages(),
    install_requires=[
    'numpy>=1.26.0',
    'scipy>=1.11.3',
    'pandas>=2.1.2',
    'scikit-learn>=1.3.2',
    'Cython>=3.0.5',
    ],
    ext_modules=cythonize(Extension("TAPpy.MSD",['TAPpy/MSD.pyx'],include_dirs=['TAPpy']))
)
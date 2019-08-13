#!/usr/bin/env python3

from setuptools import setup

requirements = [ 'Biopython', 'sys', 'os', 'random', 'string', 'shutil', 'argparse', 'functions']

readme = open('README.md').read()
setup(
    name='Pysbi',
    version='1.0',
    author='G. A. Laura, K. LL. Olfat, M. R. Marina',
    author_email=' laurag.antiga@hotmail.com, olfatk96@gmail.com, marinamrbiotech@gmail.com ',
    description='Installation package to install the Pysbi constructor',
    long_description=readme,
    requires=requirements,
    scripts=['maincompl.py', 'functions.py', 'modules_arg.py']
   
)
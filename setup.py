#!/usr/bin/env python
import numpy
from distutils.core import setup
from Cython.Build import cythonize
#from distutils.extension import Extension
#from Cython.Distutils import build_ext

#ext_modules = [
#    Extension("c_haplotypes", ["c_haplotypes.pyx"]),
#    Extension("c_reads", ["c_reads.pyx"]),
#    Extension("c_branchprune", ["c_branchprune.pyx"]),
#    Extension("HapTreePop", ["HapTreePop.pyx"])
#]

#setup(
  #name = 'HapTreePop',
#  ext_modules = ext_modules,
#  cmdclass = {'build_ext': build_ext}
#)

setup(
  ext_modules = cythonize("logprob2.pyx"),
include_dirs=[numpy.get_include()]
)
#python setup.py build_ext --inplace

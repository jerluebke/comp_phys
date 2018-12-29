# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

ext = Extension(
    "visualise",
    sources = ["./visualise.pyx", "./cvisualise.c", "../src/extern.c",
               "../src/morton.c", "../src/quadtree.c"],
    include_dirs = ["../include", ".", numpy.get_include()],
    #  libraries = ["m"]
)

setup(cmdclass = {"build_ext" : build_ext},
      ext_modules = cythonize(ext))

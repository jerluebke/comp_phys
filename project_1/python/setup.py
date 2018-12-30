# -*- coding: utf-8 -*-

import os
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

exclude = ["print_tree.c"]

ext = Extension(
    "visualise",
    sources = ["./visualise.pyx", *[os.path.join("..", "src", d) for d in
                                    os.listdir("../src") if d not in exclude]],
    include_dirs = ["../include", numpy.get_include()],
    #  libraries = ["m"],
)

setup(cmdclass = {"build_ext" : build_ext},
      ext_modules = cythonize(ext))

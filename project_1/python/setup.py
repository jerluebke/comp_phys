# -*- coding: utf-8 -*-

import os
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

# when compiling with MSVS on windows
#  os.environ["CFLAGS"] = "-std=c11"

exclude_1 = ["csearch.c"]
exclude_2 = ["cvisualise.c"]

ext_modules = [
    Extension(
        "visualise",
        sources = ["./visualise.pyx", *[os.path.join("..", "src", d) for d in
                                        os.listdir("../src") if d not in
                                        exclude_1]],
        include_dirs = ["../include", numpy.get_include(),
                        "/usr/include/graphviz"],
        #  libraries = ["m", "gvc", "cgraph", "cdt"],
        #  extra_compile_args = ["-std=c11"]
    ),
    Extension(
        "search" ,
        sources = ["./search.pyx", *[os.path.join("..", "src", d) for d in
                                     os.listdir("../src") if d not in
                                     exclude_2]],
        include_dirs = ["../include", numpy.get_include()],
        #  libraries = ["m"],
        extra_compile_args = ["-DDO_PRINT"]
    )
]

setup(
    cmdclass = {"build_ext" : build_ext},
    ext_modules = cythonize(ext_modules[0]) #, gdb_debug=True
)

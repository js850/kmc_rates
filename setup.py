import os

from numpy.distutils.core import setup
from numpy.distutils.core import Extension

depends = ["source/graph.hpp", "source/ngt.hpp"]
extra_compile_args = ["-Wall", "-Wextra", "-g", "-O2", '-funroll-loops', "-mtune=native"]
extra_link_args = []
include_dirs = ["source/"]
#extra_link_args = ["-lprofiler"]

setup(
      ext_modules=[Extension("kmc_rates/ngt", ["kmc_rates/ngt.cpp"], 
                             extra_compile_args=extra_compile_args,
                             extra_link_args=extra_link_args,
                             include_dirs=include_dirs,
                             language="c++", depends=depends,
                             ),
          ]
      )

#from distutils.core import setup
#from Cython.Build import cythonize
#
#ext_modules=cythonize(
#           "ngt.pyx",                 # our Cython source
##           sources=["source/graph.hpp", "source/ngt.hpp"],  # additional source file(s)
#           language="c++",             # generate C++ code
#           )
#
#setup(ext_modules=ext_modules)

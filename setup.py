import os

from numpy.distutils.core import setup
from numpy.distutils.core import Extension

depends = [ f for f in os.listdir(".") 
           if f.endswith(".cpp") or f.endswith(".hpp")]
extra_compile_args = ["-Wall", "-Wextra", "-g", "-O2", '-funroll-loops', "-mtune=native"]
extra_link_args = []
#extra_link_args = ["-lprofiler"]

setup(
      ext_modules=[Extension("ngt_wrapper", ["ngt_wrapper.cpp"], 
                             extra_compile_args=extra_compile_args,
                             extra_link_args=extra_link_args,
                             language="c++", depends=depends,
                             ),
          ]
      )

#!/usr/bin/env py
"""
setup.py file for gtmfe (GTFold)
"""

from distutils.core import setup, Extension

gtmfe_module = Extension('_gtmfe',
                         sources=['gtmfe-without-openmp_wrap.cxx',
                                  'src/algorithms.cc', 'src/algorithms-partition.c', 'src/energy.cc', 'src/partition-dangle.c', 'src/partition-func.c',
                                  'src/traceback.cc', 'src/AdvancedDouble.cc', 'src/constraints.cc',  'src/global.cc', 'src/key.cc', 'src/loader.cc', 'src/mfe_main.cc', 'src/options.cc', 'src/partition-func-d2.cc', 'src/shapereader.cc','src/utils.cc'],
                         swig_opts=['-c++', '-verbose'],
                         include_dirs=['include'],
                     )

setup (name = 'gtmfe',
       version = '0.1',
       author      = "Andrew Gainer-Dewar",
       description = """SWIG plumbing for gtmfe""",
       ext_modules = [gtmfe_module],
       py_modules = ["gtmfe"],
       )

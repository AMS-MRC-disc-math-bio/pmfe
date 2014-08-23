#!/usr/bin/env py
"""
setup.py file for gtmfe (GTFold)
"""

from distutils.core import setup, Extension

gtmfe_module = Extension('_gtmfe',
                         sources=['gtmfe-with-openmp_wrap.cxx',
                                  'src/algorithms.cc',
                                  'src/energy.cc',
                                  'src/traceback.cc',
                                  'src/global.cc',
                                  'src/key.cc',
                                  'src/loader.cc',
                                  'src/mfe_main.cc',
                                  'src/constraints.cc',
                                  'src/utils.cc',
                                  '../parametrizer-types/parametrizer-types.cc',],
                         swig_opts=['-c++', '-verbose'],
                         include_dirs=['include','../parametrizer-types/'],
                         libraries = ['gomp','gmp'],
                         extra_compile_args = ['-fopenmp']
                     )

setup (name = 'gtmfe',
       version = '0.1',
       author      = "Andrew Gainer-Dewar",
       description = """SWIG plumbing for gtmfe""",
       ext_modules = [gtmfe_module],
       py_modules = ["gtmfe"],
       )

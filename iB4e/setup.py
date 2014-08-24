#!/usr/bin/env py
"""
setup.py file for iB4e
"""

from distutils.core import setup, Extension

iB4e_module = Extension('_iB4e',
                        sources=['iB4e_wrap.cxx',
                                 'BBpolytope.cc',
                                 'euclideanvector.cc',
                                 'faces.cc',
                                 'linalg.cc',
                                 'stack.cc',
                                 '../parametrizer_types/parametrizer_types.cc'],
                        swig_opts=['-c++', '-verbose'],
                        include_dirs = ['../parametrizer_types/'],
                        libraries = ['gmp'],
                    )

setup (name = 'iB4e',
       version = '0.1',
       author      = "Andrew Gainer-Dewar",
       description = """SWIG plumbing for iB4e""",
       ext_modules = [iB4e_module],
       py_modules = ["iB4e"],
)

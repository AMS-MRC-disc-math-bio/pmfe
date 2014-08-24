#!/usr/bin/env py
"""
setup.py file for parametrizer_types
"""

from distutils.core import setup, Extension

parametrizer_types_module = Extension('_parametrizer_types',
                        sources=['parametrizer_types_wrap.cxx',
                                 'parametrizer_types.cc'],
                        swig_opts=['-c++', '-verbose'],
                        libraries = ['gmp'],
                                      extra_compile_args = ['-fpic'],
                    )

setup (name = 'parametrizer_types',
       version = '0.1',
       author      = "Andrew Gainer-Dewar",
       description = """SWIG plumbing for parametrizer_types library files""",
       ext_modules = [parametrizer_types_module],
       py_modules = ["parametrizer_types"],
)

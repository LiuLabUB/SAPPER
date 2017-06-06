#!/usr/bin/env python
# Time-stamp: <2017-06-01 16:04:45 Tao Liu>

"""Description: 

Setup script for SAPPER

Use this when you need Cython regenerate .c files.

Copyright (c) 2017 Tao Liu

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  beta
@version: $Revision$
@author:  Tao Liu
@contact: tliu4@buffalo.edu
"""

import os
import sys
from setuptools import setup, Extension

# Use build_ext from Cython if found
command_classes = {}
try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext
    has_cython = True
except:
    has_cython = False

try: 
    from numpy import get_include as numpy_get_include 
    numpy_include_dir = [numpy_get_include()] 
except: 
    numpy_include_dir = [] 

def main():
    if float(sys.version[:3])<3.6 or float(sys.version[:3])>=3.7:
        sys.stderr.write("CRITICAL: Python version must be 3.6!\n")
        sys.exit(1)

    # I intend to use -Ofast, however if gcc version < 4.6, this option is unavailable so...
    extra_c_args = ["-w","-O3","-ffast-math"] # for C, -Ofast implies -O3 and -ffast-math

    ext_modules = [Extension("SAPPER.PeakIO",["SAPPER/PeakIO.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("SAPPER.ReadAlignment",["SAPPER/ReadAlignment.pyx",],libraries=["m"],include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("SAPPER.RACollection",["SAPPER/RACollection.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("SAPPER.PosReadsInfo",["SAPPER/PosReadsInfo.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("SAPPER.Stat",["SAPPER/Stat.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("SAPPER.BAM",["SAPPER/BAM.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   #Extension("SAPPER.io",["SAPPER/io.pyx","SAPPER/cIO.cpp"],libraries=["m"], include_dirs=["SAPPER/",], extra_compile_args=extra_c_args),
                   #Extension("SAPPER.math",["SAPPER/math.pyx","SAPPER/cMath.cpp"],libraries=["m"], include_dirs=["SAPPER/",], extra_compile_args=extra_c_args),
                   #Extension("SAPPER.analysis",["SAPPER/analysis.pyx","SAPPER/cAnalysis.cpp"],libraries=["m"], include_dirs=["SAPPER/",], extra_compile_args=extra_c_args),
                   #Extension("SAPPER.swalign",["SAPPER/swalign.pyx","SAPPER/cSwalign.cpp"],libraries=["m"], include_dirs=["SAPPER/",], extra_compile_args=extra_c_args),
                   ]

    setup(name="SAPPER",
          version="1.0.0.20170520",
          description="de novo Variant caller for ChIP-Seq",
          author='Tao Liu',
          author_email='tliu4@buffalo.edu',
          url='http://github.com/taoliu/SAPPER/',
          package_dir={'SAPPER' : 'SAPPER'},
          packages=['SAPPER',],
          scripts=['bin/sapper',
                   ],
          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',              
              'License :: OSI Approved :: BSD License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python',
              ],
          install_requires=[
              'cython>=0.25',
              #'scipy',
              ],
          cmdclass = command_classes,
          ext_modules = ext_modules
          )

if __name__ == '__main__':
    main()

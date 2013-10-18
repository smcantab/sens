import os
import numpy as np
from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext

## Python include 
#py_include = get_python_inc() 

gsl_include = '/usr/local/include' 

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

## Gslpy include 
gslpy_include = [gsl_include, numpy_include] 



setup(
    name="sens",
    version='0.1', 
    description="Implementation of the superposition enhanced nested sampling algorithm",
    url="https://github.com/smcantab/sens",
    packages=["sens",
              "sens.src",
              "sens.models",
              "sens.tests",
             ],
    #cmdclass = {'build_ext': build_ext},
    ext_modules = 
        [
            Extension("sens.src.weighted_pick", ["sens/src/weighted_pick.c"],
                      extra_compile_args = ['-Wextra','-pedantic','-funroll-loops','-O3','-march=native','-mtune=native'],
                      include_dirs=[numpy_include],
                       ),
            Extension("sens.src.runmc", ["sens/src/runmc.c", "sens/src/mc.c", "sens/src/lj.c"],
                      include_dirs=gslpy_include,
                      extra_compile_args = ['-Wextra','-pedantic','-funroll-loops','-O3','-march=native','-mtune=native'],
                      libraries=['gsl', 'gslcblas', 'm']
                        ),
            Extension("sens.src.run_ising_mc", ["sens/src/run_ising_mc.c", "sens/src/mcising.c"],
                      include_dirs=gslpy_include,
                      extra_compile_args = ['-Wextra','-pedantic','-funroll-loops','-O3','-march=native','-mtune=native'],
                      libraries=['gsl', 'gslcblas', 'm']
                        ), 
        ]
)     

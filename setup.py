import numpy
import os
import sys
from setuptools import setup, find_packages, Extension



# Setup C module include directories
include_dirs = [numpy.get_include()]

# Setup C module macros
define_macros = [('NUMPY', '1')]

# Handle MSVC `wcsset` redefinition
if sys.platform == 'win32':
    define_macros += [
        ('_CRT_SECURE_NO_WARNING', None),
        ('__STDC__', 1)
    ]

setup(
    name="stsciutils",
    use_scm_version=True,
    url="https://github.com/spacetelescope/stsciutils",
    maintainer="J. Hunkeler",
    maintainer_email="jhunk@stsci.edu",
    setup_requires=['setuptools_scm'],
    packages=find_packages(),
    package_data={},
    ext_modules=[
        Extension('stsciutils.image._combine',
                  ['stsciutils/image/src/_combinemodule.c'],
                  include_dirs=include_dirs,
                  define_macros=define_macros),
        Extension('stsciutils.imagestats.buildHistogram',
                  ['stsciutils/imagestats/src/buildHistogram.c'],
                  include_dirs=include_dirs,
                  define_macros=define_macros),
        Extension('stsciutils.imagestats.computeMean',
                  ['stsciutils/imagestats/src/computeMean.c'],
                  include_dirs=include_dirs,
                  define_macros=define_macros),
        Extension('stsciutils.stimage._stimage',
                  ['stsciutils/stimage/src/immatch/geomap.c',
                   'stsciutils/stimage/src/immatch/xyxymatch.c',
                   'stsciutils/stimage/src/immatch/lib/tolerance.c',
                   'stsciutils/stimage/src/immatch/lib/triangles.c',
                   'stsciutils/stimage/src/immatch/lib/triangles_vote.c',
                   'stsciutils/stimage/src/lib/error.c',
                   'stsciutils/stimage/src/lib/lintransform.c',
                   'stsciutils/stimage/src/lib/polynomial.c',
                   'stsciutils/stimage/src/lib/util.c',
                   'stsciutils/stimage/src/lib/xybbox.c',
                   'stsciutils/stimage/src/lib/xycoincide.c',
                   'stsciutils/stimage/src/lib/xysort.c',
                   'stsciutils/stimage/src/surface/cholesky.c',
                   'stsciutils/stimage/src/surface/fit.c',
                   'stsciutils/stimage/src/surface/surface.c',
                   'stsciutils/stimage/src/surface/vector.c',
                   'stsciutils/stimage/src_wrap/stimage_module.c',
                   'stsciutils/stimage/src_wrap/wrap_util.c',
                   'stsciutils/stimage/src_wrap/immatch/py_xyxymatch.c',
                   'stsciutils/stimage/src_wrap/immatch/py_geomap.c'],
                   include_dirs=include_dirs + ['stsciutils/stimage/include', 'stsciutils/stimage/src_wrap'],
                   define_macros=define_macros,
        )
    ],
)

import distutils.core
import Cython.Build
distutils.core.setup(
    ext_modules=Cython.Build.cythonize(["evolution.pyx",
                                        "bodysolver.pyx",
                                        "util.pyx",
                                        "modulecore.pyx",
                                        "canon_units.pyx"]))

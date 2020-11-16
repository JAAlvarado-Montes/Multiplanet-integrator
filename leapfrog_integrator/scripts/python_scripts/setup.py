import distutils.core
import Cython.Build
distutils.core.setup(
    ext_modules=Cython.Build.cythonize(["example.pyx",
                                        "util.pyx",
                                        "functions.pyx",
                                        "canon_units.pyx"]))

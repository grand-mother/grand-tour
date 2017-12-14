import os
from distutils.core import setup, Extension

# Check the build environment.
deps_dir = os.getenv("DEPS_DIR")
if not deps_dir: deps_dir = "deps"
lib_dir = os.getenv("LIB_DIR")
if lib_dir:
    if lib_dir.endswith("python"): lib_dir = "/".join(lib_dir.split("/")[:-1])
else: lib_dir = "deps/turtle/lib"

# Configure the extension module.
module = Extension("grand_tour",
                    include_dirs = [os.path.join(deps_dir, "turtle",
                                                 "include")],
                    libraries = ["turtle"],
                    library_dirs = [lib_dir,],
                    runtime_library_dirs = [lib_dir,],
                    sources = ["src/grand-tour.c"])

setup(name = "grand_tour",
      version = "1.0",
      description = "GRAND TOpography with tURtle",
      author = "Valentin Niess",
      ext_modules = [module])

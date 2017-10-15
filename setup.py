from distutils.core import setup, Extension

module = Extension("grand_tour",
                    include_dirs = ["deps/turtle/include"],
                    libraries = ["turtle"],
                    library_dirs = ["deps/turtle/lib"],
                    sources = ["src/grand-tour.c"])

setup(name = "grand_tour",
      version = "1.0",
      description = "GRAND TOpography with tURtle",
      author = "Valentin Niess",
      ext_modules = [module])

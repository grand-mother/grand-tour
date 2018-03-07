# GRAND-TOUR
( **GRAND** **TO**pography with t**UR**tle )

## Description

GRAND-TOUR is an extension of the [TURTLE][TURTLE] library providing topographic
utilities for GRAND. It can be built as a Python C-extension module or directly
embedded in a C99 project.

[TURTLE]: https://github.com/niess/turtle

## Installation

GRAND-TOUR is built automatically as a component of
[RETRO](https://github.com/grand-mother/retro). For a standalone Python
installation you can run the following:
```bash
make deps && make
```

**Note** that you'll need ASTER-GDME2 tiles downloaded to `share/topography`
for the [example](examples/example.py) to work.

The C API is defined in [include/grand-tour.h](include/grand-tour.h). The source
code is in [src/grand-tour.c](src/grand-tour.c).

## License

The GRAND-TOUR module is under the **GNU LGPLv3** license. See the provided
`LICENSE` and `COPYING.LESSER` files.

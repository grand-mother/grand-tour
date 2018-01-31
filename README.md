# GRAND-TOUR
( **GRAND** **TO**pography with t**UR**tle )

## Description

GRAND-TOUR is a Python C-extension module providing topographic utilities for
GRAND. It is built over the [TURTLE][TURTLE] library.

[TURTLE]: https://github.com/niess/turtle

## Installation

GRAND-TOUR is built automatically as a component of
[RETRO](https://github.com/grand-mother/retro). For a standalone Installation
you can run the following:
```bash
make deps && make
. setup.sh
```

**Note** that you'll need ASTER-GDME2 tiles downloaded to `share/topography`
for the [example](examples/example.py) to work.

## License

The GRAND-TOUR module is under the **GNU LGPLv3** license. See the provided
`LICENSE` and `COPYING.LESSER` files.

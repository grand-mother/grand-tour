LIB_DIR := lib

.PHONY: all clean

all: lib/grand_tour.so

clean:
	@rm -rf lib build

lib/grand_tour.so: src/grand-tour.c
	@mkdir -p $(LIB_DIR)
	@python setup.py --quiet build --build-lib=$(LIB_DIR)
	@rm -rf build

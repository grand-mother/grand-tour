LIB_DIR := lib

.PHONY: all clean deps

all: lib/grand_tour.so

clean:
	@rm -rf lib build

lib/grand_tour.so: src/grand-tour.c
	@mkdir -p $(LIB_DIR)
	@python setup.py --quiet build --build-lib=$(LIB_DIR)
	@rm -rf build

deps: $(LIB_DIR)/libturtle.so

$(LIB_DIR)/libturtle.so: deps/turtle
	@git submodule update --init --recursive
	@$(MAKE) -C deps/turtle TURTLE_USE_PNG=0

## ------ language="Make" file="Makefile"
# allow spaces for indenting
.RECIPEPREFIX +=

# find sources
build_dir = ./build
cc_files = $(shell find ./src -name *.cc)
obj_files = $(cc_files:%.cc=$(build_dir)/%.o)
dep_files = $(obj_files:%.o=%.d)

# set compiler
compile = g++
link = g++

# libfmt
fmtlib_lflags = -lfmt

# eigen3
eigen_lflags = $(shell pkg-config --libs eigen3)
eigen_cflags = $(shell pkg-config --cflags eigen3)

# compile and link flags
compile_flags = -O3 -std=c++17 -Wall -Werror $(eigen_cflags)
link_flags = $(fmt_lflags) $(eigen_lflags)

# rules
.PHONY: clean build

build: parareal

clean:
    rm -rf $(build_dir)

-include $(dep_files)

$(build_dir)/%.o: %.cc Makefile
    @mkdir -p $(@D)
    $(compile) $(compile_flags) -MMD -c $< -o $@

parareal: $(obj_files)
    @mkdir -p $(@D)
    $(link) $^ $(link_flags) -o $@
## ------ end

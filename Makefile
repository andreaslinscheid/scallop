# To make "debug" the default configuration if invoked with just "make" except we want to build libraries.
#
ifneq ($(MAKECMDGOALS),lib)
ifeq ($(CFG),)
CFG=release
endif
else
ifeq ($(CFG),)
CFG=debug
endif
endif

# The source files: We fetch them from the source tree. We exclude the .../test/tests.cpp files
# since they contain small subprograms that excecute the tests in the given module individually.
SOURCE_FILES := $(shell find -L ./scallop/ -type f | grep .*\.cpp | grep -vP "(?:\./)?tests?\/*" )

# Build a Dependency list and an Object list, by replacing the .cpp
# extension to .d for dependency files, and .o for object files.
DEP = $(patsubst %.cpp, deps.$(CFG)/%.d, ${SOURCE_FILES})
OBJ = $(patsubst %.cpp, objs.$(CFG)/%.o, ${SOURCE_FILES})

# A list of libraries that can be compiled. We take the detour via first fetching the
# cpp files since for modules that do not have cpp files (i.e. are header only) we do not need to create a lib
BUILDLIBLIST := $(shell find -L ./scallop/ -type f | grep .*\.cpp|grep -v main.cpp | grep -vP "(?:\./)?tests?\/*"|sed -e 's/\.\/scallop\/\(\w*\).*/\1/'|uniq)


# The final binary
TARGET=scallop.x

# What compiler to use for generating dependencies: 
# it will be invoked with -MM -MP
CXXDEP = g++-5
CXX = g++-5

# What include flags to pass to the compiler
INCLUDES = -I./

# Separate compile options per configuration
ifeq ($(CFG),debug)
-include ./make.sys.debug
else
-include ./make.sys.release
endif
# A common link flag for all configurations, assables the
-include ./make.sys.common

all:	inform bin.$(CFG)/${TARGET} 

inform:
ifneq ($(CFG),release)
ifneq ($(CFG),debug)
	@echo "Invalid configuration "$(CFG)" specified."
	@echo "You must specify a configuration when running make, e.g."
	@echo  "make CFG=debug"
	@echo  
	@echo  "Possible choices for configuration are 'release' and 'debug'"
	@exit 1
endif
endif
	@echo "Configuration "$(CFG)
	@echo "------------------------"

bin.$(CFG)/${TARGET}: ${OBJ} | inform
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) $(INCLUDES) -o $@ $^ $(LDLIBS) 

objs.$(CFG)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -o $@ $<

deps.$(CFG)/%.d: %.cpp
	@mkdir -p $(dir $@)
	@echo Generating dependencies for $<
	@set -e ; $(CXXDEP) $(CXXFLAGS) $(INCLUDES) -MM -MP $< > $@.$$$$; \
	sed 's,\(.*\)\.o[ :]*,objs.$(CFG)\/$(shell echo "$@"| sed -e 's/^[^\/]*\/\(.*\)\.d$$/\1/').o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$


.PHONY: clean lib
clean:
	@rm -rf \
	deps.debug objs.debug bin.debug \
	deps.release objs.release bin.release \
	lib debuglib

lib : ${OBJ} | inform
ifeq ($(CFG),release)
	@mkdir -p lib
	@$(foreach lib, $(BUILDLIBLIST), echo "Building library $(lib)";a="";for i in $(OBJ); do a=$$a" `echo $$i|grep $(lib)`";done;ar rs lib/lib$(lib).a $$a; echo "";) 
else
	@mkdir -p debuglib
	@$(foreach lib, $(BUILDLIBLIST), echo "Building debugging library $(lib)";a="";for i in $(OBJ); do a=$$a" `echo $$i|grep $(lib)`";done;ar rs debuglib/lib$(lib).a $$a; echo "";) 
endif

# Generate the list of objects that belong to this library
		  


# Unless "make clean" is called, include the dependency files
# which are auto-generated. Don't fail if they are missing
# (-include), since they will be missing in the first invocation!
ifneq ($(MAKECMDGOALS),clean)
-include ${DEP}
endif



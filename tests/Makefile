# To make "debug" the default configuration if invoked with just "make":
#
ifeq ($(CFG),)
CFG=debug
endif

# The source files: We fetch them from the source tree
SOURCE_FILES := $(shell find -L . -type f | grep .*\.cpp )
#SOURCE_FILES = ../gslpp/

# Build a Dependency list and an Object list, by replacing the .cpp
# extension to .d for dependency files, and .o for object files.
DEP = $(patsubst %.cpp, deps.$(CFG)/%.d, ${SOURCE_FILES})
OBJ = $(patsubst %.cpp, objs.$(CFG)/%.o, ${SOURCE_FILES})

# The final binary
TARGET=tests.x

# What compiler to use for generating dependencies: 
# it will be invoked with -MM -MP
CXXDEP = g++

# What include flags to pass to the compiler
INCLUDES = -I./

# Separate compile options per configuration
ifeq ($(CFG),debug)
-include ../make.sys.debug
else
-include ../make.sys.release
endif
# A common link flag for all configurations, assables the
-include ../make.sys.common

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

.PHONY: clean
clean:
	@rm -rf \
	deps.debug objs.debug bin.debug \
	deps.release objs.release bin.release


# Unless "make clean" is called, include the dependency files
# which are auto-generated. Don't fail if they are missing
# (-include), since they will be missing in the first invocation!
ifneq ($(MAKECMDGOALS),clean)
-include ${DEP}
endif



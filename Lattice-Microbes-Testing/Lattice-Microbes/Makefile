
# Include local settings.
ifdef LOCALMK
-include $(LOCALMK)
else
-include local.mk
endif

# Define the object lists filled in by the subdir files.
MAIN :=
MAIN_OBJS :=
OBJS :=
CPP_DEPS :=
PROTOBUF_GEN_SRC :=

# Global build settings.
INCLUDE_DIRS := -I./src/c $(HDF5_INCLUDE_DIR) $(PROTOBUF_INCLUDE_DIR) -I./profile -I./$(BUILD_DIR)/src/protobuf
LIB_DIRS := $(HDF5_LIB_DIR) $(PROTOBUF_LIB_DIR)
LIBS := $(HDF5_LIB) $(PROTOBUF_LIB)
INSTALL_TARGETS := installlm

# Add definitions for utils.
UTIL_BINS :=
UTIL_OBJS :=
UTIL_INCLUDE_DIRS := -I./utils/c

# Main binary libs
LIBS += $(PYTHON_LIB)
LM_LIBS := $(LIBS)
LM_LIB_DIRS := $(LIB_DIRS)

# Add definitions for CUDA.
ifeq ($(USE_CUDA),1)
CCFLAGS += 
CXXFLAGS += 
CUDA_FLAGS += 
INCLUDE_DIRS += -I./src/cuda $(CUDA_INCLUDE_DIR)
LIB_DIRS += $(CUDA_LIB_DIR)
LIBS += $(CUDA_LIB)
LM_LIB_DIRS += $(CUDA_LIB_DIR)
LM_LIBS += $(CUDA_LIB)
endif

# Add definitions for Python.
ifeq ($(USE_PYTHON),1)
INCLUDE_DIRS += $(PYTHON_INCLUDE_DIR)
LIB_DIRS += $(PYTHON_LIB_DIR)
SWIG_GEN_SRC := 
INSTALL_LIBS += ./$(BUILD_DIR)/lib/_lm.so
endif

# Add definitions for MPI.
ifeq ($(USE_MPI),1)
CCFLAGS += $(MPI_COMPILE_FLAGS)
CXXFLAGS += $(MPI_COMPILE_FLAGS)
CUDA_FLAGS += $(MPI_COMPILE_FLAGS)
LDFLAGS += $(MPI_LINK_FLAGS)
endif

# Add definitions for SBML.
ifeq ($(USE_SBML),1)
CCFLAGS += 
CXXFLAGS += 
INCLUDE_DIRS += $(SBML_INCLUDE_DIR)
LIB_DIRS += $(SBML_LIB_DIR)
LIBS += $(SBML_LIB)
endif

# Add definitions for VMD.
ifeq ($(USE_VMD),1)
VMD_PLUGIN :=
VMD_PLUGIN_OBJS :=
VMD_PLUGIN_DEP_OBJS :=
INCLUDE_DIRS += -I./vmd/c $(VMD_INCLUDE_DIR)
INSTALL_TARGETS += installvmd
endif


# Add definitions for debugging.
CCFLAGS += 
CXXFLAGS += 
ifeq ($(USE_VERBOSITY_LEVEL),10)
CCFLAGS += 
CXXFLAGS += 
CUDA_FLAGS += 
endif

# Timing specific build settings.
TIMING_BINS :=
TIMING_OBJS :=
TIMING_INCLUDE_DIRS := -I./timing/c -I./timing/cuda
TIMING_CCFLAGS= $(CCFLAGS)
TIMING_CXXFLAGS := $(CXXFLAGS) 
TIMING_CUDA_FLAGS := $(CUDA_FLAGS)

# Testing specific build settings.
TESTING_BINS :=
TESTING_OBJS :=
TESTING_CXXFLAGS := $(BOOST_TEST_CXXFLAGS)
TESTING_INCLUDE_DIRS := -I./testing/c -I./testing/cuda $(BOOST_TEST_INCLUDE_DIR)
TESTING_LIB_DIRS := $(BOOST_TEST_LIB_DIR)
TESTING_LIBS := $(BOOST_TEST_LIB)

# Add definitions for profiling.
ifeq ($(USE_PROF),1)
ifeq ($(USE_PROF_NVTX),1)
CCFLAGS += 
CXXFLAGS += 
CUDA_FLAGS += 
LIBS += -lnvToolsExt
else
CCFLAGS +=
CXXFLAGS +=
CUDA_FLAGS +=
endif
endif

# All of the sources participating in the build are defined here
-include src/c/subdir.mk
-include src/protobuf/subdir.mk
-include testing/c/subdir.mk
-include timing/c/subdir.mk
-include utils/c/subdir.mk

ifeq ($(USE_CUDA),1)
-include src/cuda/subdir.mk
#-include testing/cuda/subdir.mk
-include timing/cuda/subdir.mk
endif

ifeq ($(USE_PYTHON),1)
-include src/swig/subdir.mk
-include src/python/subdir.mk
endif

ifeq ($(USE_VMD),1)
-include vmd/c/subdir.mk
endif

# Set the build info.
BUILD_INFO := "\"by $(shell whoami) on $(shell hostname -s) at $(shell date "+%Y-%m-%d %H:%M:%S")\""
CCFLAGS +=
CXXFLAGS +=

.SECONDARY:

default: all
.DEFAULT_GOAL := all

# All target.
ifeq ($(USE_PYTHON),1)
ifeq ($(USE_MPI),1)
all: $(MAIN) util vmdplugin pylm pylmExamples pylmInstall ./$(BUILD_DIR)/lm-mpmpd
else
all: $(MAIN) util vmdplugin pylm pylmExamples pylmInstall
endif
else
ifeq ($(USE_MPI),1)
all: $(MAIN) util vmdplugin ./$(BUILD_DIR)/lm-mpmpd
else
all: $(MAIN) util vmdplugin
endif
endif

# Build main LK program.
$(MAIN): protobuf $(MAIN_OBJS) $(OBJS)
	mkdir -p $(BUILD_DIR)
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	$(LD) $(LDFLAGS) $(MAIN_OBJS) $(OBJS) $(LM_LIB_DIRS) $(LM_LIBS) -o "$@" 
	mkdir -p ./$(BUILD_DIR)/bin
	cp "$@" ./$(BUILD_DIR)/bin
	@echo 'Finished building target: $@'
	@echo ' '

# Build MPMPD Binary
./$(BUILD_DIR)/lm-mpmpd: protobuf $(MAIN_MPMPD_OBJS) $(OBJS)
	mkdir -p $(BUILD_DIR)
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	$(LD) $(LDFLAGS) $(MAIN_MPMPD_OBJS) $(OBJS) $(LM_LIB_DIRS) $(LM_LIBS) -o "$@" 
	mkdir -p ./$(BUILD_DIR)/bin
	cp "$@" ./$(BUILD_DIR)/bin
	@echo 'Finished building target: $@'
	@echo ' '
	
# Build util binaries.
util: $(UTIL_BINS)

# Build utilities.
./$(BUILD_DIR)/bin/lm_%: ./$(BUILD_DIR)/utils/c/lm_%.o $(UTIL_OBJS) $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	mkdir -p $(@D)
	$(LD) $(LDFLAGS) $< $(UTIL_OBJS) $(OBJS) $(LM_LIB_DIRS) $(LM_LIBS) -o "$@" 
	@echo 'Finished building target: $@'
	@echo ' '

# Build protcol buffers source files.
protobuf: $(PROTOBUF_GEN_SRC)


# Build timing binaries.
time: $(TIMING_BINS)

# Build timing.
./$(BUILD_DIR)/timing/c/%: ./$(BUILD_DIR)/timing/c/%.o $(TIMING_OBJS) $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	$(LD) $(LDFLAGS) $< $(TIMING_OBJS) $(LIB_DIRS) $(LIBS) -o "$@"
	mkdir -p ./$(BUILD_DIR)/bin/timing
	cp "$@" ./$(BUILD_DIR)/bin/timing
	@echo 'Finished building target: $@'
	@echo ' '

./$(BUILD_DIR)/timing/cuda/%: ./$(BUILD_DIR)/timing/cuda/%.o $(TIMING_OBJS) $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	$(LD) $(LDFLAGS) $< $(TIMING_OBJS) $(LIB_DIRS) $(LIBS) -o "$@" 
	mkdir -p ./$(BUILD_DIR)/bin/timing
	cp "$@" ./$(BUILD_DIR)/bin/timing
	@echo 'Finished building target: $@'
	@echo ' '


# Build testing.
check: test
	@cd $(BUILD_DIR)/bin/testing; for t in $(TESTING_BINS); do echo "----------------------------------------\nUnit test: `basename $$t`"; ./`basename $$t`; echo "----------------------------------------"; done;

test: $(TESTING_BINS)
./$(BUILD_DIR)/testing/c/%: ./$(BUILD_DIR)/testing/c/%.o $(TESTING_OBJS) $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	$(LD) $(LDFLAGS) $< $(TESTING_OBJS) $(OBJS) $(TESTING_LIB_DIRS) $(LIB_DIRS) $(TESTING_LIBS) $(LIBS) -o "$@"
	mkdir -p ./$(BUILD_DIR)/bin/testing
	cp "$@" ./$(BUILD_DIR)/bin/testing
	@echo 'Finished building target: $@'
	@echo ' '

./$(BUILD_DIR)/testing/cuda/%: ./$(BUILD_DIR)/testing/cuda/%.o $(TESTING_OBJS) $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Linker'
	$(LD) $(LDFLAGS) $< $(TESTING_OBJS) $(OBJS) $(TESTING_LIB_DIRS) $(LIB_DIRS) $(TESTING_LIBS) $(LIBS) -o "$@" 
	mkdir -p ./$(BUILD_DIR)/bin/testing
	cp "$@" ./$(BUILD_DIR)/bin/testing
	@echo 'Finished building target: $@'
	@echo ' '


# Build the vmd molfile plugin.
vmdplugin: $(VMD_PLUGIN)
	@echo $(VMD_PLUGIN)

$(VMD_PLUGIN): $(VMD_PLUGIN_OBJS) $(VMD_PLUGIN_DEP_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: $(SHLD)'
	mkdir -p $(@D)
	$(SHLD) $(SHLDFLAGS) $(VMD_PLUGIN_OBJS) $(VMD_PLUGIN_DEP_OBJS) $(HDF5_LIB_DIR) $(PROTOBUF_LIB_DIR) $(HDF5_LIB) $(PROTOBUF_LIB) $(SHLDOPTO) "$@"
	@echo 'Finished building target: $@'

# Build the python shared lib
ifeq ($(USE_PROF),1)
pylm: 
	@echo 'disable profiling (USE_PROF) to build _lm.so'
else
pylm: ./$(BUILD_DIR)/lib/_lm.so
	@echo 'Python interface (lm.py, _lm.so) bulit in ./$(BUILD_DIR)/lib'
endif

# Clean up any generated files.
clean:
	@echo 'Cleaning'
	$(RM) $(BUILD_DIR)

install: $(INSTALL_TARGETS)

# Install the programs.
installlm:
	@echo 'Installing LM binaries'
	install -d $(INSTALL_PREFIX)/bin
	install -m 0755 $(MAIN) $(UTIL_BINS) $(INSTALL_PREFIX)/bin
	@echo 'Installing LM libraries'
	install -d $(INSTALL_PREFIX)/lib/lm
	install -m 0644 $(INSTALL_LIBS) $(INSTALL_PREFIX)/lib/lm
	install -d $(INSTALL_PREFIX)/lib/python/pySTDLM
	install -d $(INSTALL_PREFIX)/lib/python/pyLM
	install -m 0644 ./src/python/pyLM/* $(INSTALL_PREFIX)/lib/python/pyLM
	install -m 0644 ./src/python/pySTDLM/* $(INSTALL_PREFIX)/lib/python/pySTDLM
ifeq ($(USE_MPI),1)
	install -m 0755 $(BUILD_DIR)/bin/lm-mpmpd $(INSTALL_PREFIX)/bin
endif


installvmd:
	@echo 'Installing VMD plugin'
	install -d $(VMD_INSTALL_DIR)
	install -m 0644 $(VMD_PLUGIN) $(VMD_INSTALL_DIR)
	
# Build a distribution package.
dist:
	@echo 'Packaging distribution'
	rm -rf tmp
	mkdir -p tmp/lm-2.0
	cp Makefile tmp/lm-2.0
	cp -r docs tmp/lm-2.0
	cp -r profile tmp/lm-2.0
	cp -r src tmp/lm-2.0
	cp -r testing tmp/lm-2.0
	cp -r timing tmp/lm-2.0
	cp -r utils tmp/lm-2.0
	cp -r vmd tmp/lm-2.0
	tar zcvf lm-2.0_src.tgz  --exclude="*/.*" -C tmp lm-2.0
	rm -rf tmp
	
# Including any generated dependencies.
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif


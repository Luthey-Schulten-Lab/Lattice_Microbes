
# Add inputs and outputs from these tool invocations to the build variables.
VMD_PLUGIN := \
$(BUILD_DIR)/lib/lmplugin.so

VMD_PLUGIN_OBJS := \
./$(BUILD_DIR)/vmd/c/lmplugin.o

VMD_PLUGIN_DEP_OBJS := \
./$(BUILD_DIR)/src/c/lm/Print.o \
./$(BUILD_DIR)/src/protobuf/lm/io/DiffusionModel.pb.o \
./$(BUILD_DIR)/src/protobuf/lm/io/SpatialModel.pb.o \
./$(BUILD_DIR)/src/protobuf/lm/io/ReactionModel.pb.o \
./$(BUILD_DIR)/src/c/lm/io/hdf5/SimulationFile.o \
./$(BUILD_DIR)/src/c/lm/rdme/ByteLattice.o \
./$(BUILD_DIR)/src/c/lm/rdme/Lattice.o

CPP_DEPS += \
./$(BUILD_DIR)/vmd/c/lmplugin.d

# Each subdirectory must supply rules for building sources it contributes
./$(BUILD_DIR)/vmd/c/%.o: ./vmd/c/%.cpp 
	@echo 'Building file: $<'
	@echo 'Invoking: $(CXX)'
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) $(CXXDEPENDFLAGS) -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

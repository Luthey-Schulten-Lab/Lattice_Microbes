
# Add inputs and outputs from these tool invocations to the build variables 
TESTING_BINS += \
./$(BUILD_DIR)/testing/c/boost_test \
./$(BUILD_DIR)/testing/c/lm/cme/cme_solver \
./$(BUILD_DIR)/testing/c/lm/cme/gillespie_d_solver \
./$(BUILD_DIR)/testing/c/lm/cme/fluctuating_nr_solver \
./$(BUILD_DIR)/testing/c/lm/io/hdf5/simulation_file \
./$(BUILD_DIR)/testing/c/lm/main/data_output_queue \
./$(BUILD_DIR)/testing/c/lm/main/resource_allocator \
./$(BUILD_DIR)/testing/c/lm/rdme/byte_lattice \
./$(BUILD_DIR)/testing/c/lm/rdme/lattice \
./$(BUILD_DIR)/testing/c/lm/reaction/reaction_queue \
./$(BUILD_DIR)/testing/c/lm/rng/xorshift

#./$(BUILD_DIR)/testing/c/lm/main/signal_handler \

TESTING_OBJS := \
./$(BUILD_DIR)/testing/c/TestingHelper.o

CPP_DEPS += \
./$(BUILD_DIR)/testing/c/boost_test.d \
./$(BUILD_DIR)/testing/c/TestingHelper.d \
./$(BUILD_DIR)/testing/c/lm/cme/gillespie_d_solver.d \
./$(BUILD_DIR)/testing/c/lm/io/hdf5/simulation_file.d \
./$(BUILD_DIR)/testing/c/lm/main/data_output_queue.d \
./$(BUILD_DIR)/testing/c/lm/main/resource_allocator.d \
./$(BUILD_DIR)/testing/c/lm/main/signal_handler.d \
./$(BUILD_DIR)/testing/c/lm/rdme/byte_lattice.d \
./$(BUILD_DIR)/testing/c/lm/rdme/lattice.d \
./$(BUILD_DIR)/testing/c/lm/reaction/reaction_queue.d \
./$(BUILD_DIR)/testing/c/lm/rng/xorshift.d

# Each subdirectory must supply rules for building sources it contributes
./$(BUILD_DIR)/testing/c/%.o: ./testing/c/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: $(CXX)'
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(TESTING_CXXFLAGS) $(TESTING_INCLUDE_DIRS) $(INCLUDE_DIRS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

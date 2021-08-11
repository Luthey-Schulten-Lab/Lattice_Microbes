
# Add inputs and outputs from these tool invocations to the build variables 
TIMING_BINS +=

TIMING_OBJS +=

CPP_DEPS +=
 
# Each subdirectory must supply rules for building sources it contributes
./$(BUILD_DIR)/timing/c/%.o: ./timing/c/%.cpp 
	@echo 'Building file: $<'
	@echo 'Invoking: $(CXX)'
	mkdir -p ./$(BUILD_DIR)/timing/c
	$(CXX) $(PROFNOFILE_CXXFLAGS) $(TIMING_INCLUDE_DIRS) $(INCLUDE_DIRS) -DPROF_OUT_FILE=timing.`basename $@ .o`.prof -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

./$(BUILD_DIR)/timing/c/%.o: ./timing/c/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: $(CC)'
	mkdir -p ./$(BUILD_DIR)/timing/c
	$(CC) $(PROFNOFILE_CCFLAGS) $(TIMING_INCLUDE_DIRS) $(INCLUDE_DIRS) -DPROF_OUT_FILE=timing.`basename $@ .o`.prof -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

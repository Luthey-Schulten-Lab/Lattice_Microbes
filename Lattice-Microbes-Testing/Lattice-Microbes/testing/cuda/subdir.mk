
# Add inputs and outputs from these tool invocations to the build variables 
TESTING_BINS += \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window_periodic \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window_periodic_unaligned \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_y_window \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_y_window_periodic \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_z_window \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_z_window_periodic \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/packed_site_skip \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_2d_dev/copy_xy_window \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_2d_dev/copy_xy_window_periodic \
./$(BUILD_DIR)/testing/cuda/lm/rng/xorwow

#./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_3d_dev/copy_xyz_window

TESTING_OBJS += \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window_periodic.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window_periodic_unaligned.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_y_window.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_y_window_periodic.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_z_window.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_z_window_periodic.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/packed_site_skip.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_2d_dev/copy_xy_window.cu.o \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_2d_dev/copy_xy_window_periodic.cu.o

#./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_3d_dev/copy_xyz_window.cu.o

CPP_DEPS += \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window_periodic.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window_periodic.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window_periodic_unaligned.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_x_window_periodic_unaligned.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_y_window.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_y_window.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_y_window_periodic.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_y_window_periodic.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_z_window.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_z_window.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_z_window_periodic.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/copy_z_window_periodic.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/packed_site_skip.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_1d_dev/packed_site_skip.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_2d_dev/copy_xy_window.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_2d_dev/copy_xy_window_periodic.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_2d_dev/copy_xy_window.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_2d_dev/copy_xy_window_periodic.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rng/xorwow.d

#./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_3d_dev/copy_xyz_window.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_3d_dev/copy_xyz_window.cu.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_3dw_dev/copy_xyz_window.d \
./$(BUILD_DIR)/testing/cuda/lm/rdme/lattice_sim_3dw_dev/copy_xyz_window.cu.d
 
# Each subdirectory must supply rules for building sources it contributes
./$(BUILD_DIR)/testing/cuda/%.o: ./testing/cuda/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: $(CXX)'
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(TESTING_CXXFLAGS) $(TESTING_INCLUDE_DIRS) $(INCLUDE_DIRS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

./$(BUILD_DIR)/testing/cuda/%.cu.o: ./testing/cuda/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: $(NVCC)'
	mkdir -p $(@D)
	$(CUDA_NVCC) $(CUDA_FLAGS) $(TESTING_INCLUDE_DIRS) $(INCLUDE_DIRS) -M "$<" -odir $(@D) -o "$(@:%.o=%.d)" 
	$(CUDA_NVCC) $(CUDA_FLAGS) $(TESTING_INCLUDE_DIRS) $(INCLUDE_DIRS) -c "$<" -o "$@" 
ifeq ($(GENERATE_CUDA_PTX_CODE),1)
	$(CUDA_NVCC) $(CUDA_FLAGS) $(TESTING_INCLUDE_DIRS) $(INCLUDE_DIRS) -ptx --opencc-options=-LIST:source=on  "$<" -o "$(@:%.o=%.ptx)" 
endif
ifeq ($(GENERATE_CUDA_BIN_CODE),1)
	$(CUDA_NVCC) $(CUDA_FLAGS) $(TESTING_INCLUDE_DIRS) $(INCLUDE_DIRS) -cubin "$<" -o "$(@:%.o=%.cubin)" 
endif
ifeq ($(GENERATE_CUDA_ASM_CODE),1)
	/usr/local/decuda/decuda.py -o "$(@:%.o=%.asm)" "$(@:%.o=%.cubin)"
endif
	@echo 'Finished building: $<'
	@echo ' '

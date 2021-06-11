
# Add inputs and outputs from these tool invocations to the build variables 
TIMING_BINS += \
./$(BUILD_DIR)/timing/cuda/lm/rdme/rng_generate

#./$(BUILD_DIR)/timing/cuda/bit_packed_diffusion_1d \
#./$(BUILD_DIR)/timing/cuda/byte_diffusion_1d \
#./$(BUILD_DIR)/timing/cuda/byte_diffusion_2d \
#./$(BUILD_DIR)/timing/cuda/byte_diffusion_3dw \
#./$(BUILD_DIR)/timing/cuda/kernel_launch \

TIMING_OBJS +=

CPP_DEPS += \
./$(BUILD_DIR)/timing/cuda/lm/rdme/rng_generate.d

#./$(BUILD_DIR)/timing/cuda/bit_packed_diffusion_1d.d \
#./$(BUILD_DIR)/timing/cuda/byte_diffusion_1d.d \
#./$(BUILD_DIR)/timing/cuda/byte_diffusion_2d.d \
#./$(BUILD_DIR)/timing/cuda/byte_diffusion_3dw.d \
#./$(BUILD_DIR)/timing/cuda/kernel_launch.d \
 
# Each subdirectory must supply rules for building sources it contributes
./$(BUILD_DIR)/timing/cuda/%.o: ./timing/cuda/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: $(NVCC)'
	mkdir -p $(@D)
	$(CUDA_NVCC) $(TIMING_CUDA_FLAGS) -DPROF_OUT_FILE=timing.`basename $@ .o`.prof $(TIMING_INCLUDE_DIRS) $(INCLUDE_DIRS) -M "$<" -odir $(@D) -o "$(@:%.o=%.d)" 
	$(CUDA_NVCC) $(TIMING_CUDA_FLAGS) -DPROF_OUT_FILE=timing.`basename $@ .o`.prof $(TIMING_INCLUDE_DIRS) $(INCLUDE_DIRS) -c "$<" -o "$@" 
ifeq ($(GENERATE_CUDA_PTX_CODE),1)
	$(CUDA_NVCC) $(TIMING_CUDA_FLAGS) -DPROF_OUT_FILE=timing.`basename $@ .o`.prof $(TIMING_INCLUDE_DIRS) $(INCLUDE_DIRS) -ptx --opencc-options=-LIST:source=on  "$<" -o "$(@:%.o=%.ptx)" 
endif
ifeq ($(GENERATE_CUDA_BIN_CODE),1)
	$(CUDA_NVCC) $(TIMING_CUDA_FLAGS) -DPROF_OUT_FILE=timing.`basename $@ .o`.prof $(TIMING_INCLUDE_DIRS) $(INCLUDE_DIRS) -cubin "$<" -o "$(@:%.o=%.cubin)" 
endif
ifeq ($(GENERATE_CUDA_ASM_CODE),1)
	/usr/local/decuda/decuda.py -o "$(@:%.o=%.asm)" "$(@:%.o=%.cubin)"
endif
	@echo 'Finished building: $<'
	@echo ' '

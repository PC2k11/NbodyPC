################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/BarnesHut.c \
../src/exception.c \
../src/genlib.c \
../src/random.c \
../src/stack.c 

OBJS += \
./src/BarnesHut.o \
./src/exception.o \
./src/genlib.o \
./src/random.o \
./src/stack.o 

C_DEPS += \
./src/BarnesHut.d \
./src/exception.d \
./src/genlib.d \
./src/random.d \
./src/stack.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



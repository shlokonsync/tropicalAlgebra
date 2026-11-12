# Host Compiler & Flags
CC ?= gcc
CFLAGS ?= -Wall -Wextra -Wpedantic -Wconversion -Werror -O2 -std=c99

# Cortex-M4 Cross-Compilation Toolchain
ARM_CC ?= clang
ARM_AR ?= llvm-ar
ARM_NM ?= llvm-nm
ARM_SIZE ?= llvm-size

# ARM Cortex-M4 Target Flags (STM32F407: ARMv7E-M, Cortex-M4 with FPU)
# -ffreestanding -fno-builtin guarantee zero libc dependencies
ARM_CFLAGS ?= --target=arm-none-eabi -mcpu=cortex-m4 -mthumb -mfloat-abi=hard -mfpu=fpv4-sp-d16 \
              -Wall -Wextra -Wpedantic -Wconversion -Werror -O2 -std=c99 \
              -ffreestanding -fno-builtin

HOST_SRCS = tropical.c karp.c test_week1.c
HOST_OBJS = $(HOST_SRCS:.c=.o)
HOST_TARGET = test_week1

ARM_SRCS = tropical.c karp.c
ARM_OBJS = tropical_arm.o karp_arm.o
ARM_ASM  = tropical.s karp.s
ARM_LIB  = libtropical.a

all: $(HOST_TARGET)

# Host Test Executable
$(HOST_TARGET): $(HOST_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c tropical.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(HOST_TARGET)
	./$(HOST_TARGET)

# Cortex-M4 Cross-Compilation
arm: $(ARM_LIB)

$(ARM_LIB): $(ARM_OBJS)
	$(ARM_AR) rcs $@ $^

%_arm.o: %.c tropical.h
	$(ARM_CC) $(ARM_CFLAGS) -c $< -o $@

asm: $(ARM_ASM)

%.s: %.c tropical.h
	$(ARM_CC) $(ARM_CFLAGS) -S $< -o $@

# Verify Cortex-M4 Artifacts (No libc dependencies, Section Footprints)
check-arm: arm asm
	@echo "==============================================================="
	@echo "  ARM Cortex-M4 Symbol Table Verification (Undefined Symbols)  "
	@echo "==============================================================="
	$(ARM_NM) -u $(ARM_LIB)
	@echo "==============================================================="
	@echo "  ARM Cortex-M4 Memory Footprint (Static Buffers in .bss)      "
	@echo "==============================================================="
	$(ARM_SIZE) $(ARM_OBJS) $(ARM_LIB)

clean:
	@python -c "import os, glob; [os.remove(f) for pat in ('*.o', '*.s', '*.a', '*.exe', '$(HOST_TARGET)') for f in glob.glob(pat)]"

.PHONY: all test arm asm check-arm clean

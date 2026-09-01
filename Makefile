CC ?= cc
CFLAGS ?= -O2 -std=c99 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Isrc
LDFLAGS ?=
LDLIBS ?=
LDLIBS += -lm

SRC_DIR := src
BUILD_DIR := build
SIMD_BUILD_DIR := build/simd
SIMD_CFLAGS ?= -O3 -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512

PROGRAM := $(BUILD_DIR)/subseqhash3
SOURCES := main.c subseqhash3.c table_file.c ssh3_simulate.c \
	ssh3_tables.c ssh3_seeds.c ssh3_collisions.c ssh3_test.c
OBJECTS := $(SOURCES:%.c=$(BUILD_DIR)/%.o)
DEPS := $(OBJECTS:.o=.d)

.PHONY: all simd test clean clean-simd

all: $(PROGRAM)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(PROGRAM): $(OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test: $(PROGRAM)
	$(PROGRAM) test

simd: clean-simd
	$(MAKE) BUILD_DIR=$(SIMD_BUILD_DIR) \
		CFLAGS="$(CFLAGS) $(SIMD_CFLAGS)" all

clean-simd:
	rm -f $(SIMD_BUILD_DIR)/*.o $(SIMD_BUILD_DIR)/*.d $(SIMD_BUILD_DIR)/subseqhash3
	rmdir $(SIMD_BUILD_DIR) 2>/dev/null || true

clean: clean-simd
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.d $(PROGRAM)
	rmdir $(BUILD_DIR) 2>/dev/null || true

-include $(DEPS)

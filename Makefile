CC ?= cc
CFLAGS ?= -O2 -std=c99 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Isrc
LDFLAGS ?=
LDLIBS ?=

SRC_DIR := src
BUILD_DIR := build
SIMD_BUILD_DIR := build/simd
SIMD_CFLAGS ?= -O3 -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512

I8_TEST := $(BUILD_DIR)/test_subseqhash_edits_i8
I8_BENCHMARK := $(BUILD_DIR)/benchmark_subseqhash_edits_i8
I8_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_i8_collisions
I8_SEED_GENERATOR := $(BUILD_DIR)/generate_subseqhash_i8_seeds
EDIT_PAIR_SIMULATOR := $(BUILD_DIR)/simulate_edit_pairs
I8_OBJECTS := $(BUILD_DIR)/subseqhash_edits.o \
	$(BUILD_DIR)/subseqhash_edits_i8.o

.PHONY: all simd test benchmark benchmark-i8-collisions \
	generate-subseqhash-i8-seeds simulate-edit-pairs clean clean-simd

all: $(I8_TEST) $(I8_BENCHMARK) $(I8_COLLISION_BENCHMARK) \
	$(I8_SEED_GENERATOR) $(EDIT_PAIR_SIMULATOR)

simd: clean-simd
	$(MAKE) BUILD_DIR=$(SIMD_BUILD_DIR) \
		CFLAGS="$(CFLAGS) $(SIMD_CFLAGS)" all

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/subseqhash_edits.o: $(SRC_DIR)/subseqhash_edits.c \
	$(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/subseqhash_edits_i8.o: $(SRC_DIR)/subseqhash_edits_i8.c \
	$(SRC_DIR)/subseqhash_edits_i8.h $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash_edits_i8.o: \
	$(SRC_DIR)/test_subseqhash_edits_i8.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(I8_TEST): $(I8_OBJECTS) $(BUILD_DIR)/test_subseqhash_edits_i8.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test: $(I8_TEST)
	$(I8_TEST)

$(BUILD_DIR)/benchmark_subseqhash_edits_i8.o: \
	$(SRC_DIR)/benchmark_subseqhash_edits_i8.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(I8_BENCHMARK): $(I8_OBJECTS) $(BUILD_DIR)/benchmark_subseqhash_edits_i8.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

benchmark: $(I8_BENCHMARK)
	$(I8_BENCHMARK)

$(BUILD_DIR)/benchmark_i8_collisions.o: \
	$(SRC_DIR)/benchmark_i8_collisions.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(I8_COLLISION_BENCHMARK): $(I8_OBJECTS) $(BUILD_DIR)/benchmark_i8_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

benchmark-i8-collisions: $(I8_COLLISION_BENCHMARK)

$(BUILD_DIR)/generate_subseqhash_i8_seeds.o: \
	$(SRC_DIR)/generate_subseqhash_i8_seeds.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(I8_SEED_GENERATOR): $(I8_OBJECTS) \
	$(BUILD_DIR)/generate_subseqhash_i8_seeds.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

generate-subseqhash-i8-seeds: $(I8_SEED_GENERATOR)

$(BUILD_DIR)/simulate_edit_pairs.o: $(SRC_DIR)/simulate_edit_pairs.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(EDIT_PAIR_SIMULATOR): $(BUILD_DIR)/simulate_edit_pairs.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

simulate-edit-pairs: $(EDIT_PAIR_SIMULATOR)

clean-simd:
	rm -f $(SIMD_BUILD_DIR)/*.o \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_i8 \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_i8 \
		$(SIMD_BUILD_DIR)/benchmark_i8_collisions \
		$(SIMD_BUILD_DIR)/generate_subseqhash_i8_seeds \
		$(SIMD_BUILD_DIR)/simulate_edit_pairs
	rmdir $(SIMD_BUILD_DIR) 2>/dev/null || true

clean: clean-simd
	rm -f $(BUILD_DIR)/*.o $(I8_TEST) $(I8_BENCHMARK) \
		$(I8_COLLISION_BENCHMARK) $(I8_SEED_GENERATOR) \
		$(EDIT_PAIR_SIMULATOR)
	rmdir $(BUILD_DIR) 2>/dev/null || true

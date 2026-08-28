CC ?= cc
CFLAGS ?= -O2 -std=c99 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Isrc
LDFLAGS ?=
LDLIBS ?=

SRC_DIR := src
BUILD_DIR := build
SIMD_BUILD_DIR := build/simd
SIMD_CFLAGS ?= -O3 -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512
TARGET := $(BUILD_DIR)/test_subsequence_scores
BENCHMARK := $(BUILD_DIR)/benchmark_subsequence_scores
WITNESS_GENERATOR := $(BUILD_DIR)/generate_witnesses
COLLISION_ESTIMATOR := $(BUILD_DIR)/estimate_hash_collisions
SUBSEQHASH_TEST := $(BUILD_DIR)/test_subseqhash
SUBSEQHASH_EDITS_TEST := $(BUILD_DIR)/test_subseqhash_edits
EDIT_PAIR_SIMULATOR := $(BUILD_DIR)/simulate_edit_pairs
EDIT_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_edit_collisions
SUBSEQHASH_ADJACENT_TEST := $(BUILD_DIR)/test_subseqhash_adjacent
ADJACENT_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_adjacent_collisions
SUBSEQHASH_EDITS_I16_TEST := $(BUILD_DIR)/test_subseqhash_edits_i16
SUBSEQHASH_EDITS_I8_TEST := $(BUILD_DIR)/test_subseqhash_edits_i8
SUBSEQHASH2_EDITS_I8_TEST := $(BUILD_DIR)/test_subseqhash2_edits_i8
SUBSEQHASH2_SEED_GENERATOR := $(BUILD_DIR)/generate_subseqhash2_seeds
SUBSEQHASH_I8_SEED_GENERATOR := $(BUILD_DIR)/generate_subseqhash_i8_seeds
SUBSEQHASH_EDITS_A10_TEST := $(BUILD_DIR)/test_subseqhash_edits_a10
SUBSEQHASH_EDITS_I16_BENCHMARK := $(BUILD_DIR)/benchmark_subseqhash_edits_i16
SUBSEQHASH_EDITS_I8_BENCHMARK := $(BUILD_DIR)/benchmark_subseqhash_edits_i8
SUBSEQHASH_EDITS_A10_BENCHMARK := $(BUILD_DIR)/benchmark_subseqhash_edits_a10
I16_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_i16_collisions
I8_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_i8_collisions
I8_TWOPASS_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_i8_twopass_collisions
A10_TWOPASS_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_a10_twopass_collisions
MAX_I8_TEST := $(BUILD_DIR)/test_subseqhash_edits_max_i8
MAX_I8_BENCHMARK := $(BUILD_DIR)/benchmark_subseqhash_edits_max_i8
MAX_I8_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_max_i8_collisions
SUBSEQHASH_BENCHMARK := $(BUILD_DIR)/benchmark_subseqhash
LIB_OBJECT := $(BUILD_DIR)/subsequence_scores.o

.PHONY: all simd test benchmark benchmark-subseqhash benchmark-subseqhash-edits-i16 benchmark-edit-collisions benchmark-i16-collisions benchmark-i8-collisions benchmark-adjacent-collisions witnesses collision-estimator simulate-edit-pairs generate-subseqhash2-seeds generate-subseqhash-i8-seeds sanitize clean clean-simd

all: $(TARGET) \
	$(BENCHMARK) \
	$(WITNESS_GENERATOR) \
	$(COLLISION_ESTIMATOR) \
	$(SUBSEQHASH_TEST) \
	$(SUBSEQHASH_EDITS_TEST) \
	$(EDIT_PAIR_SIMULATOR) \
	$(EDIT_COLLISION_BENCHMARK) \
	$(SUBSEQHASH_ADJACENT_TEST) \
	$(ADJACENT_COLLISION_BENCHMARK) \
	$(SUBSEQHASH_EDITS_I16_TEST) \
	$(SUBSEQHASH_EDITS_I8_TEST) \
	$(SUBSEQHASH2_EDITS_I8_TEST) \
	$(SUBSEQHASH2_SEED_GENERATOR) \
	$(SUBSEQHASH_I8_SEED_GENERATOR) \
	$(SUBSEQHASH_EDITS_A10_TEST) \
	$(SUBSEQHASH_EDITS_I16_BENCHMARK) \
	$(SUBSEQHASH_EDITS_I8_BENCHMARK) \
	$(SUBSEQHASH_EDITS_A10_BENCHMARK) \
	$(I16_COLLISION_BENCHMARK) \
	$(I8_COLLISION_BENCHMARK) \
	$(I8_TWOPASS_COLLISION_BENCHMARK) \
	$(A10_TWOPASS_COLLISION_BENCHMARK) \
	$(MAX_I8_TEST) \
	$(MAX_I8_BENCHMARK) \
	$(MAX_I8_COLLISION_BENCHMARK) \
	$(SUBSEQHASH_BENCHMARK)

# Build every program that exercises the SIMD subseqhash-edits implementation.
# A separate object directory prevents non-SIMD objects from a prior ordinary
# build from being silently reused.
simd: clean-simd
	$(MAKE) BUILD_DIR=$(SIMD_BUILD_DIR) \
		CFLAGS="$(CFLAGS) $(SIMD_CFLAGS)" \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_i16 \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_i8 \
		$(SIMD_BUILD_DIR)/test_subseqhash2_edits_i8 \
		$(SIMD_BUILD_DIR)/generate_subseqhash2_seeds \
		$(SIMD_BUILD_DIR)/generate_subseqhash_i8_seeds \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_i16 \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_i8 \
		$(SIMD_BUILD_DIR)/benchmark_i16_collisions \
		$(SIMD_BUILD_DIR)/benchmark_i8_collisions \
		$(SIMD_BUILD_DIR)/benchmark_i8_twopass_collisions \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_a10 \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_a10 \
		$(SIMD_BUILD_DIR)/benchmark_a10_twopass_collisions \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_max_i8 \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_max_i8 \
		$(SIMD_BUILD_DIR)/benchmark_max_i8_collisions

$(BUILD_DIR):
	mkdir -p $@

$(LIB_OBJECT): $(SRC_DIR)/subsequence_scores.c $(SRC_DIR)/subsequence_scores.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subsequence_scores.o: $(SRC_DIR)/test_subsequence_scores.c $(SRC_DIR)/subsequence_scores.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TARGET): $(LIB_OBJECT) $(BUILD_DIR)/test_subsequence_scores.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test: $(TARGET) $(SUBSEQHASH_TEST) $(SUBSEQHASH_EDITS_TEST) $(SUBSEQHASH_ADJACENT_TEST) $(SUBSEQHASH_EDITS_I16_TEST) $(SUBSEQHASH_EDITS_I8_TEST) $(SUBSEQHASH2_EDITS_I8_TEST) $(SUBSEQHASH_EDITS_A10_TEST) $(MAX_I8_TEST)
	./$(TARGET)
	./$(SUBSEQHASH_TEST)
	./$(SUBSEQHASH_EDITS_TEST)
	./$(SUBSEQHASH_ADJACENT_TEST)
	./$(SUBSEQHASH_EDITS_I16_TEST)
	./$(SUBSEQHASH_EDITS_I8_TEST)
	./$(SUBSEQHASH2_EDITS_I8_TEST)
	./$(SUBSEQHASH_EDITS_A10_TEST)
	./$(MAX_I8_TEST)

$(BUILD_DIR)/benchmark_subsequence_scores.o: $(SRC_DIR)/benchmark_subsequence_scores.c $(SRC_DIR)/subsequence_scores.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BENCHMARK): $(LIB_OBJECT) $(BUILD_DIR)/benchmark_subsequence_scores.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

benchmark: $(BENCHMARK)
	./$(BENCHMARK)

$(BUILD_DIR)/generate_witnesses.o: $(SRC_DIR)/generate_witnesses.c $(SRC_DIR)/subsequence_scores.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(WITNESS_GENERATOR): $(LIB_OBJECT) $(BUILD_DIR)/generate_witnesses.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

witnesses: $(WITNESS_GENERATOR)
	./$(WITNESS_GENERATOR)

$(BUILD_DIR)/hash_collision_estimator.o: $(SRC_DIR)/hash_collision_estimator.c $(SRC_DIR)/hash_collision_estimator.h $(SRC_DIR)/subseqhash.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/estimate_hash_collisions.o: $(SRC_DIR)/estimate_hash_collisions.c $(SRC_DIR)/hash_collision_estimator.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(COLLISION_ESTIMATOR): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/hash_collision_estimator.o $(BUILD_DIR)/estimate_hash_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

collision-estimator: $(COLLISION_ESTIMATOR)

$(BUILD_DIR)/subseqhash.o: $(SRC_DIR)/subseqhash.c $(SRC_DIR)/subseqhash.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash.o: $(SRC_DIR)/test_subseqhash.c $(SRC_DIR)/subseqhash.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_TEST): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/test_subseqhash.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/subseqhash_edits.o: $(SRC_DIR)/subseqhash_edits.c $(SRC_DIR)/subseqhash_edits.h $(SRC_DIR)/subseqhash.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash_edits.o: $(SRC_DIR)/test_subseqhash_edits.c $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_EDITS_TEST): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/test_subseqhash_edits.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/subseqhash_edits_i16.o: $(SRC_DIR)/subseqhash_edits_i16.c $(SRC_DIR)/subseqhash_edits_i16.h $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash_edits_i16.o: $(SRC_DIR)/test_subseqhash_edits_i16.c $(SRC_DIR)/subseqhash_edits_i16.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_EDITS_I16_TEST): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i16.o $(BUILD_DIR)/test_subseqhash_edits_i16.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/benchmark_subseqhash_edits_i16.o: $(SRC_DIR)/benchmark_subseqhash_edits_i16.c $(SRC_DIR)/subseqhash_edits_i16.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_EDITS_I16_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i16.o $(BUILD_DIR)/benchmark_subseqhash_edits_i16.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

benchmark-subseqhash-edits-i16: $(SUBSEQHASH_EDITS_I16_BENCHMARK)
	./$(SUBSEQHASH_EDITS_I16_BENCHMARK)

$(BUILD_DIR)/subseqhash_edits_a10.o: $(SRC_DIR)/subseqhash_edits_a10.c $(SRC_DIR)/subseqhash_edits_a10.h $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash_edits_a10.o: $(SRC_DIR)/test_subseqhash_edits_a10.c $(SRC_DIR)/subseqhash_edits_a10.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_EDITS_A10_TEST): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_a10.o $(BUILD_DIR)/test_subseqhash_edits_a10.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/benchmark_subseqhash_edits_a10.o: $(SRC_DIR)/benchmark_subseqhash_edits_a10.c $(SRC_DIR)/subseqhash_edits_a10.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_EDITS_A10_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_a10.o $(BUILD_DIR)/benchmark_subseqhash_edits_a10.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/subseqhash_edits_i8.o: $(SRC_DIR)/subseqhash_edits_i8.c $(SRC_DIR)/subseqhash_edits_i8.h $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash_edits_i8.o: $(SRC_DIR)/test_subseqhash_edits_i8.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_EDITS_I8_TEST): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i8.o $(BUILD_DIR)/test_subseqhash_edits_i8.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/subseqhash2_edits_i8.o: $(SRC_DIR)/subseqhash2_edits_i8.c $(SRC_DIR)/subseqhash2_edits_i8.h $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash2_edits_i8.o: $(SRC_DIR)/test_subseqhash2_edits_i8.c $(SRC_DIR)/subseqhash2_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH2_EDITS_I8_TEST): $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash2_edits_i8.o $(BUILD_DIR)/test_subseqhash2_edits_i8.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/generate_subseqhash2_seeds.o: $(SRC_DIR)/generate_subseqhash2_seeds.c $(SRC_DIR)/subseqhash2_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH2_SEED_GENERATOR): $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash2_edits_i8.o $(BUILD_DIR)/generate_subseqhash2_seeds.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

generate-subseqhash2-seeds: $(SUBSEQHASH2_SEED_GENERATOR)

$(BUILD_DIR)/generate_subseqhash_i8_seeds.o: $(SRC_DIR)/generate_subseqhash_i8_seeds.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_I8_SEED_GENERATOR): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i8.o $(BUILD_DIR)/generate_subseqhash_i8_seeds.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

generate-subseqhash-i8-seeds: $(SUBSEQHASH_I8_SEED_GENERATOR)

$(BUILD_DIR)/benchmark_subseqhash_edits_i8.o: $(SRC_DIR)/benchmark_subseqhash_edits_i8.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_EDITS_I8_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i8.o $(BUILD_DIR)/benchmark_subseqhash_edits_i8.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/benchmark_i16_collisions.o: $(SRC_DIR)/benchmark_i16_collisions.c $(SRC_DIR)/subseqhash_edits_i16.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(I16_COLLISION_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i16.o $(BUILD_DIR)/benchmark_i16_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

benchmark-i16-collisions: $(I16_COLLISION_BENCHMARK)

$(BUILD_DIR)/benchmark_i8_collisions.o: $(SRC_DIR)/benchmark_i8_collisions.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(I8_COLLISION_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i8.o $(BUILD_DIR)/benchmark_i8_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

benchmark-i8-collisions: $(I8_COLLISION_BENCHMARK)

$(BUILD_DIR)/subseqhash_adjacent.o: $(SRC_DIR)/subseqhash_adjacent.c $(SRC_DIR)/subseqhash_adjacent.h $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash_adjacent.o: $(SRC_DIR)/test_subseqhash_adjacent.c $(SRC_DIR)/subseqhash_adjacent.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_ADJACENT_TEST): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_adjacent.o $(BUILD_DIR)/test_subseqhash_adjacent.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/simulate_edit_pairs.o: $(SRC_DIR)/simulate_edit_pairs.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(EDIT_PAIR_SIMULATOR): $(BUILD_DIR)/simulate_edit_pairs.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

simulate-edit-pairs: $(EDIT_PAIR_SIMULATOR)

$(BUILD_DIR)/benchmark_edit_collisions.o: $(SRC_DIR)/benchmark_edit_collisions.c $(SRC_DIR)/subseqhash.h $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(EDIT_COLLISION_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/benchmark_edit_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

benchmark-edit-collisions: $(EDIT_COLLISION_BENCHMARK)

$(BUILD_DIR)/benchmark_adjacent_collisions.o: $(SRC_DIR)/benchmark_adjacent_collisions.c $(SRC_DIR)/subseqhash_adjacent.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(ADJACENT_COLLISION_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_adjacent.o $(BUILD_DIR)/benchmark_adjacent_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

benchmark-adjacent-collisions: $(ADJACENT_COLLISION_BENCHMARK)

$(BUILD_DIR)/benchmark_subseqhash.o: $(SRC_DIR)/benchmark_subseqhash.c $(SRC_DIR)/subseqhash.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/benchmark_subseqhash.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

benchmark-subseqhash: $(SUBSEQHASH_BENCHMARK)
	./$(SUBSEQHASH_BENCHMARK)

$(BUILD_DIR)/benchmark_a10_twopass_collisions.o: $(SRC_DIR)/benchmark_a10_twopass_collisions.c $(SRC_DIR)/subseqhash_edits_a10.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(A10_TWOPASS_COLLISION_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_a10.o $(BUILD_DIR)/benchmark_a10_twopass_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/benchmark_i8_twopass_collisions.o: $(SRC_DIR)/benchmark_i8_twopass_collisions.c $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(I8_TWOPASS_COLLISION_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i8.o $(BUILD_DIR)/benchmark_i8_twopass_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/subseqhash_edits_max_i8.o: $(SRC_DIR)/subseqhash_edits_max_i8.c $(SRC_DIR)/subseqhash_edits_max_i8.h $(SRC_DIR)/subseqhash_edits.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subseqhash_edits_max_i8.o: $(SRC_DIR)/test_subseqhash_edits_max_i8.c $(SRC_DIR)/subseqhash_edits_max_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(MAX_I8_TEST): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_max_i8.o $(BUILD_DIR)/test_subseqhash_edits_max_i8.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/benchmark_subseqhash_edits_max_i8.o: $(SRC_DIR)/benchmark_subseqhash_edits_max_i8.c $(SRC_DIR)/subseqhash_edits_max_i8.h $(SRC_DIR)/subseqhash_edits_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(MAX_I8_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_i8.o $(BUILD_DIR)/subseqhash_edits_max_i8.o $(BUILD_DIR)/benchmark_subseqhash_edits_max_i8.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

$(BUILD_DIR)/benchmark_max_i8_collisions.o: $(SRC_DIR)/benchmark_max_i8_collisions.c $(SRC_DIR)/subseqhash_edits_max_i8.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(MAX_I8_COLLISION_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/subseqhash_edits.o $(BUILD_DIR)/subseqhash_edits_max_i8.o $(BUILD_DIR)/benchmark_max_i8_collisions.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

sanitize: | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) -O1 -g -std=c99 -Wall -Wextra -Wpedantic \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		$(SRC_DIR)/subsequence_scores.c $(SRC_DIR)/test_subsequence_scores.c \
		-o $(BUILD_DIR)/test_subsequence_scores_sanitized
	./$(BUILD_DIR)/test_subsequence_scores_sanitized

clean-simd:
	rm -f $(SIMD_BUILD_DIR)/*.o \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_i16 \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_i8 \
		$(SIMD_BUILD_DIR)/test_subseqhash2_edits_i8 \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_i16 \
		$(SIMD_BUILD_DIR)/generate_subseqhash2_seeds \
		$(SIMD_BUILD_DIR)/generate_subseqhash_i8_seeds \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_i8 \
		$(SIMD_BUILD_DIR)/benchmark_i16_collisions \
		$(SIMD_BUILD_DIR)/benchmark_i8_collisions \
		$(SIMD_BUILD_DIR)/benchmark_i8_twopass_collisions \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_a10 \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_a10 \
		$(SIMD_BUILD_DIR)/benchmark_a10_twopass_collisions \
		$(SIMD_BUILD_DIR)/test_subseqhash_edits_max_i8 \
		$(SIMD_BUILD_DIR)/benchmark_subseqhash_edits_max_i8 \
		$(SIMD_BUILD_DIR)/benchmark_max_i8_collisions
	rmdir $(SIMD_BUILD_DIR) 2>/dev/null || true

clean: clean-simd
	rm -f $(BUILD_DIR)/*.o $(TARGET) $(BENCHMARK) $(WITNESS_GENERATOR) \
		$(COLLISION_ESTIMATOR) \
		$(SUBSEQHASH_TEST) \
		$(SUBSEQHASH_EDITS_TEST) \
		$(EDIT_PAIR_SIMULATOR) \
		$(EDIT_COLLISION_BENCHMARK) \
		$(SUBSEQHASH_ADJACENT_TEST) \
		$(ADJACENT_COLLISION_BENCHMARK) \
		$(SUBSEQHASH_EDITS_I16_TEST) \
		$(SUBSEQHASH_EDITS_I8_TEST) \
		$(SUBSEQHASH2_EDITS_I8_TEST) \
		$(SUBSEQHASH2_SEED_GENERATOR) \
		$(SUBSEQHASH_I8_SEED_GENERATOR) \
		$(SUBSEQHASH_EDITS_I16_BENCHMARK) \
		$(SUBSEQHASH_EDITS_I8_BENCHMARK) \
		$(I16_COLLISION_BENCHMARK) \
		$(I8_COLLISION_BENCHMARK) \
		$(I8_TWOPASS_COLLISION_BENCHMARK) \
		$(SUBSEQHASH_EDITS_A10_TEST) \
		$(SUBSEQHASH_EDITS_A10_BENCHMARK) \
		$(A10_TWOPASS_COLLISION_BENCHMARK) \
		$(MAX_I8_TEST) \
		$(MAX_I8_BENCHMARK) \
		$(MAX_I8_COLLISION_BENCHMARK) \
		$(SUBSEQHASH_BENCHMARK) \
		$(BUILD_DIR)/test_subsequence_scores_sanitized
	rmdir $(BUILD_DIR) 2>/dev/null || true

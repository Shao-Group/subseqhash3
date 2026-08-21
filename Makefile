CC ?= cc
CFLAGS ?= -O2 -std=c99 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Isrc
LDFLAGS ?=
LDLIBS ?=

SRC_DIR := src
BUILD_DIR := build
TARGET := $(BUILD_DIR)/test_subsequence_scores
BENCHMARK := $(BUILD_DIR)/benchmark_subsequence_scores
WITNESS_GENERATOR := $(BUILD_DIR)/generate_witnesses
COLLISION_ESTIMATOR := $(BUILD_DIR)/estimate_hash_collisions
SUBSEQHASH_TEST := $(BUILD_DIR)/test_subseqhash
SUBSEQHASH_EDITS_TEST := $(BUILD_DIR)/test_subseqhash_edits
EDIT_PAIR_SIMULATOR := $(BUILD_DIR)/simulate_edit_pairs
EDIT_COLLISION_BENCHMARK := $(BUILD_DIR)/benchmark_edit_collisions
SUBSEQHASH_ADJACENT_TEST := $(BUILD_DIR)/test_subseqhash_adjacent
SUBSEQHASH_BENCHMARK := $(BUILD_DIR)/benchmark_subseqhash
LIB_OBJECT := $(BUILD_DIR)/subsequence_scores.o

.PHONY: all test benchmark benchmark-subseqhash benchmark-edit-collisions witnesses collision-estimator simulate-edit-pairs sanitize clean

all: $(TARGET) \
	$(BENCHMARK) \
	$(WITNESS_GENERATOR) \
	$(COLLISION_ESTIMATOR) \
	$(SUBSEQHASH_TEST) \
	$(SUBSEQHASH_EDITS_TEST) \
	$(EDIT_PAIR_SIMULATOR) \
	$(EDIT_COLLISION_BENCHMARK) \
	$(SUBSEQHASH_ADJACENT_TEST) \
	$(SUBSEQHASH_BENCHMARK)

$(BUILD_DIR):
	mkdir -p $@

$(LIB_OBJECT): $(SRC_DIR)/subsequence_scores.c $(SRC_DIR)/subsequence_scores.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subsequence_scores.o: $(SRC_DIR)/test_subsequence_scores.c $(SRC_DIR)/subsequence_scores.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TARGET): $(LIB_OBJECT) $(BUILD_DIR)/test_subsequence_scores.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test: $(TARGET) $(SUBSEQHASH_TEST) $(SUBSEQHASH_EDITS_TEST) $(SUBSEQHASH_ADJACENT_TEST)
	./$(TARGET)
	./$(SUBSEQHASH_TEST)
	./$(SUBSEQHASH_EDITS_TEST)
	./$(SUBSEQHASH_ADJACENT_TEST)

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

$(BUILD_DIR)/benchmark_subseqhash.o: $(SRC_DIR)/benchmark_subseqhash.c $(SRC_DIR)/subseqhash.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SUBSEQHASH_BENCHMARK): $(BUILD_DIR)/subseqhash.o $(BUILD_DIR)/benchmark_subseqhash.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -lm -o $@

benchmark-subseqhash: $(SUBSEQHASH_BENCHMARK)
	./$(SUBSEQHASH_BENCHMARK)

sanitize: | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) -O1 -g -std=c99 -Wall -Wextra -Wpedantic \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		$(SRC_DIR)/subsequence_scores.c $(SRC_DIR)/test_subsequence_scores.c \
		-o $(BUILD_DIR)/test_subsequence_scores_sanitized
	./$(BUILD_DIR)/test_subsequence_scores_sanitized

clean:
	rm -f $(BUILD_DIR)/*.o $(TARGET) $(BENCHMARK) $(WITNESS_GENERATOR) \
		$(COLLISION_ESTIMATOR) \
		$(SUBSEQHASH_TEST) \
		$(SUBSEQHASH_EDITS_TEST) \
		$(EDIT_PAIR_SIMULATOR) \
		$(EDIT_COLLISION_BENCHMARK) \
		$(SUBSEQHASH_ADJACENT_TEST) \
		$(SUBSEQHASH_BENCHMARK) \
		$(BUILD_DIR)/test_subsequence_scores_sanitized
	rmdir $(BUILD_DIR) 2>/dev/null || true

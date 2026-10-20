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
	ssh3_tables.c ssh3_seeds.c ssh3_collisions.c ssh3_test.c ssh3_pivots.c \
	pivot_table_file.c ssh3_repeats.c ssh3_deletions.c ssh3_deletions_cli.c
OBJECTS := $(SOURCES:%.c=$(BUILD_DIR)/%.o)
DEPS := $(OBJECTS:.o=.d)

.PHONY: all simd test bench-one check-one bench-pivots check-pivots bench-pivots-speed check-pivots-speed clean clean-simd
.PHONY: bench-pivots-single check-pivots-single check-repeats

all: $(PROGRAM)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(PROGRAM): $(OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test: $(PROGRAM) $(BUILD_DIR)/check-repeats $(BUILD_DIR)/check-deletions
	$(PROGRAM) test
	$(BUILD_DIR)/check-deletions 2>$(BUILD_DIR)/deletions-test.log || { cat $(BUILD_DIR)/deletions-test.log >&2; exit 1; }
	$(BUILD_DIR)/check-repeats 2>$(BUILD_DIR)/repeats-test.log || { cat $(BUILD_DIR)/repeats-test.log >&2; exit 1; }

$(BUILD_DIR)/check-repeats: tests/repeats.c $(BUILD_DIR)/subseqhash3.o $(BUILD_DIR)/ssh3_pivots.o $(BUILD_DIR)/pivot_table_file.o $(BUILD_DIR)/table_file.o $(BUILD_DIR)/ssh3_tables.o $(BUILD_DIR)/ssh3_repeats.o $(BUILD_DIR)/ssh3_seeds.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-repeats: $(BUILD_DIR)/check-repeats
	$(BUILD_DIR)/check-repeats 2>$(BUILD_DIR)/repeats-test.log || { cat $(BUILD_DIR)/repeats-test.log >&2; exit 1; }

$(BUILD_DIR)/bench-one: bench/one_substitution.c $(BUILD_DIR)/subseqhash3.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

bench-one: $(BUILD_DIR)/bench-one
	$(BUILD_DIR)/bench-one

check-one: $(BUILD_DIR)/bench-one
	$(BUILD_DIR)/bench-one --check-only

$(BUILD_DIR)/bench-pivots: bench/pivots.c $(BUILD_DIR)/ssh3_pivots.o $(BUILD_DIR)/subseqhash3.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

bench-pivots: $(BUILD_DIR)/bench-pivots
	$(BUILD_DIR)/bench-pivots

check-pivots: $(BUILD_DIR)/bench-pivots
	$(BUILD_DIR)/bench-pivots --check-only

# Immutable pre-optimization pivot reference, compiled with the SAME flags.
# Requires the repository history containing ca35567; no source checkout edits.
$(BUILD_DIR)/pivots-baseline.o: src/ssh3_pivots.h src/subseqhash3.h | $(BUILD_DIR)
	git show ca35567:src/ssh3_pivots.c | $(CC) $(CPPFLAGS) $(CFLAGS) \
	  -Dsubseqhash3_pivot_tables_init=baseline_pivot_tables_init \
	  -Dsubseqhash3_pivot_tables_destroy=baseline_pivot_tables_destroy \
	  -Dsubseqhash3_pivot_workspace_create=baseline_pivot_workspace_create \
	  -Dsubseqhash3_pivot_workspace_destroy=baseline_pivot_workspace_destroy \
	  -Dsubseqhash3_compute_pivot_seeds=baseline_compute_pivot_seeds \
	  -Dsubseqhash3_compute_pivot_repeats=baseline_compute_pivot_repeats \
	  -Dsubseqhash3_compute_pivot_repeats_with_workspace=baseline_compute_pivot_repeats_with_workspace \
	  -x c -c - -o $@

$(BUILD_DIR)/bench-pivots-speed: bench/pivots_speed.c $(BUILD_DIR)/ssh3_pivots.o $(BUILD_DIR)/pivots-baseline.o $(BUILD_DIR)/subseqhash3.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

bench-pivots-speed: $(BUILD_DIR)/bench-pivots-speed
	$(BUILD_DIR)/bench-pivots-speed

check-pivots-speed: $(BUILD_DIR)/bench-pivots-speed
	$(BUILD_DIR)/bench-pivots-speed --check-only

# Restored single-threaded baseline immediately before the second round.
$(BUILD_DIR)/pivots-round1.o: bench/reference/pivots_round1.c src/ssh3_pivots.h src/subseqhash3.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	  -Dsubseqhash3_pivot_tables_init=baseline_pivot_tables_init \
	  -Dsubseqhash3_pivot_tables_destroy=baseline_pivot_tables_destroy \
	  -Dsubseqhash3_pivot_workspace_create=baseline_pivot_workspace_create \
	  -Dsubseqhash3_pivot_workspace_destroy=baseline_pivot_workspace_destroy \
	  -Dsubseqhash3_compute_pivot_seeds=baseline_compute_pivot_seeds \
	  -Dsubseqhash3_compute_pivot_repeats=baseline_compute_pivot_repeats \
	  -Dsubseqhash3_compute_pivot_repeats_with_workspace=baseline_compute_pivot_repeats_with_workspace \
	  -c $< -o $@

$(BUILD_DIR)/bench-pivots-single: bench/pivots_speed.c $(BUILD_DIR)/ssh3_pivots.o $(BUILD_DIR)/pivots-round1.o $(BUILD_DIR)/subseqhash3.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -DPIVOT_BASELINE_LABEL='"restored-single-thread-round1"' $(LDFLAGS) $^ $(LDLIBS) -o $@

bench-pivots-single: $(BUILD_DIR)/bench-pivots-single
	$(BUILD_DIR)/bench-pivots-single

check-pivots-single: $(BUILD_DIR)/bench-pivots-single
	$(BUILD_DIR)/bench-pivots-single --check-only

simd: clean-simd
	$(MAKE) BUILD_DIR=$(SIMD_BUILD_DIR) \
		CFLAGS="$(CFLAGS) $(SIMD_CFLAGS)" all

clean-simd:
	rm -f $(SIMD_BUILD_DIR)/*.o $(SIMD_BUILD_DIR)/*.d $(SIMD_BUILD_DIR)/subseqhash3 $(SIMD_BUILD_DIR)/bench-one $(SIMD_BUILD_DIR)/bench-pivots $(SIMD_BUILD_DIR)/bench-pivots-speed $(SIMD_BUILD_DIR)/bench-pivots-single $(SIMD_BUILD_DIR)/check-repeats $(SIMD_BUILD_DIR)/repeats-test.log $(SIMD_BUILD_DIR)/check-deletions $(SIMD_BUILD_DIR)/bench-deletions $(SIMD_BUILD_DIR)/deletions-test.log
	rmdir $(SIMD_BUILD_DIR) 2>/dev/null || true

clean: clean-simd
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.d $(PROGRAM) $(BUILD_DIR)/bench-one $(BUILD_DIR)/bench-pivots $(BUILD_DIR)/bench-pivots-speed $(BUILD_DIR)/bench-pivots-single $(BUILD_DIR)/check-repeats $(BUILD_DIR)/repeats-test.log $(BUILD_DIR)/check-deletions $(BUILD_DIR)/bench-deletions $(BUILD_DIR)/deletions-test.log
	rmdir $(BUILD_DIR) 2>/dev/null || true

-include $(DEPS)

.PHONY: check-deletions
$(BUILD_DIR)/check-deletions: tests/deletions.c $(BUILD_DIR)/ssh3_deletions.o $(BUILD_DIR)/ssh3_deletions_cli.o $(BUILD_DIR)/ssh3_pivots.o $(BUILD_DIR)/subseqhash3.o $(BUILD_DIR)/pivot_table_file.o $(BUILD_DIR)/ssh3_tables.o $(BUILD_DIR)/table_file.o $(BUILD_DIR)/ssh3_repeats.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-deletions: $(BUILD_DIR)/check-deletions
	$(BUILD_DIR)/check-deletions 2>$(BUILD_DIR)/deletions-test.log || { cat $(BUILD_DIR)/deletions-test.log >&2; exit 1; }

# Optional comparison against the user's unmodified SubseqHash2 checkout.
SUBSEQHASH2_DIR ?= subseqhash2
CXX ?= c++
DELETION_CXXFLAGS ?= $(filter-out -std=c99,$(CFLAGS)) -std=c++11
$(BUILD_DIR)/bench-deletions: bench/deletions.cpp $(BUILD_DIR)/ssh3_deletions.o $(BUILD_DIR)/ssh3_pivots.o $(BUILD_DIR)/subseqhash3.o $(SUBSEQHASH2_DIR)/src/subseqhash2seeding.cpp $(SUBSEQHASH2_DIR)/src/seeding.cpp
	$(CXX) $(CPPFLAGS) -I$(SUBSEQHASH2_DIR)/src $(DELETION_CXXFLAGS) -include cstdint -include cstdio -include string $(LDFLAGS) $^ $(LDLIBS) -o $@

.PHONY: bench-deletions check-deletions-reference
bench-deletions: $(BUILD_DIR)/bench-deletions
	$(BUILD_DIR)/bench-deletions

check-deletions-reference: $(BUILD_DIR)/bench-deletions
	$(BUILD_DIR)/bench-deletions --check-only

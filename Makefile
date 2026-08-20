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
LIB_OBJECT := $(BUILD_DIR)/subsequence_scores.o

.PHONY: all test benchmark witnesses sanitize clean

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $@

$(LIB_OBJECT): $(SRC_DIR)/subsequence_scores.c $(SRC_DIR)/subsequence_scores.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_subsequence_scores.o: $(SRC_DIR)/test_subsequence_scores.c $(SRC_DIR)/subsequence_scores.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TARGET): $(LIB_OBJECT) $(BUILD_DIR)/test_subsequence_scores.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test: $(TARGET)
	./$(TARGET)

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

sanitize: | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) -O1 -g -std=c99 -Wall -Wextra -Wpedantic \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		$(SRC_DIR)/subsequence_scores.c $(SRC_DIR)/test_subsequence_scores.c \
		-o $(BUILD_DIR)/test_subsequence_scores_sanitized
	./$(BUILD_DIR)/test_subsequence_scores_sanitized

clean:
	rm -f $(BUILD_DIR)/*.o $(TARGET) $(BENCHMARK) $(WITNESS_GENERATOR) \
		$(BUILD_DIR)/test_subsequence_scores_sanitized
	rmdir $(BUILD_DIR) 2>/dev/null || true

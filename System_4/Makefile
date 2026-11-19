# System 4 -- Fixed-Point & Bit-Accurate C Model
CC      ?= gcc
CSTD    ?= -std=c99
WARN    := -Wall -Wextra -Wpedantic
OPT     ?= -O2
INC     := -Iinclude -Iinclude/core -Iinclude/model -Iinclude/analysis -Ithird_party/icd
CFLAGS  := $(CSTD) $(WARN) $(OPT) $(INC)
LDLIBS  := -lm
BUILD   := build
OBJ     := $(BUILD)/obj
BIN     := $(BUILD)/bin

SRC     := $(wildcard src/core/*.c) $(wildcard src/model/*.c) $(wildcard src/analysis/*.c)
LIB_OBJ := $(patsubst src/%.c,$(OBJ)/%.o,$(SRC))
HDRS    := $(wildcard include/*.h include/*/*.h)

UNIT    := test_fixed_math test_divider test_cordic test_adapter test_quant_stress test_models
UNIT_BIN:= $(addprefix $(BIN)/,$(UNIT))
PROOF   := $(BIN)/system4_proof

.PHONY: all lib tests test proof analysis clean
all: lib tests $(PROOF) $(BIN)/run_analysis

lib: $(LIB_OBJ)

$(OBJ)/%.o: src/%.c $(HDRS)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_%: tests/unit/test_%.c $(LIB_OBJ) tests/unit/test_util.h
	@mkdir -p $(BIN)
	$(CC) $(CFLAGS) -Itests/unit $< $(LIB_OBJ) -o $@ $(LDLIBS)

$(PROOF): tests/proof/system4_proof.c $(LIB_OBJ) $(HDRS)
	@mkdir -p $(BIN)
	$(CC) $(CFLAGS) $< $(LIB_OBJ) -o $@ $(LDLIBS)

$(BIN)/run_analysis: tools/run_analysis.c $(LIB_OBJ) $(HDRS)
	@mkdir -p $(BIN)
	$(CC) $(CFLAGS) $< $(LIB_OBJ) -o $@ $(LDLIBS)

tests: $(UNIT_BIN)

test: all
	./scripts/run_all_tests.sh

proof: $(PROOF)
	$(PROOF) results

analysis: $(BIN)/run_analysis
	$(BIN)/run_analysis results 60
	python3 scripts/plot_results.py results

clean:
	rm -rf $(BUILD)

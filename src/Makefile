# A052075 — Primes p such that nextprime(p) is substring of p^3
# Makefile v2.0 — Post multi-LLM review
#
# Targets:
#   make           — optimized build (auto-detect CPU)
#   make zen4      — optimized for Ryzen 9 7940HS / Zen 4 [F1, F2]
#   make debug     — sanitizers for debugging
#   make verify    — build + oracle verification (15/15)
#   make test      — completeness test (first 7 terms)
#   make test8     — extended test (first 8 terms)
#   make clean     — remove binaries and temp files

CC = gcc
CFLAGS_BASE = -Wall -Wextra -Wpedantic -std=c11
CFLAGS_OPT = $(CFLAGS_BASE) -O2 -march=native
CFLAGS_ZEN4 = $(CFLAGS_BASE) -O3 -march=znver4 -mtune=znver4 -flto  # [F1, F2]
CFLAGS_OMP = -fopenmp
CFLAGS_DBG = $(CFLAGS_BASE) -O0 -g -fsanitize=address,undefined
LIBS = -lgmp -lm

SRC = a052075_v4.c
BIN = a052075
BIN_DBG = a052075_debug

.PHONY: all zen4 debug verify test test8 clean

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS_OPT) $(CFLAGS_OMP) -o $@ $< $(LIBS)
	@echo "Built: $@ (optimized, -march=native)"

zen4: $(SRC)
	$(CC) $(CFLAGS_ZEN4) $(CFLAGS_OMP) -o $(BIN) $< $(LIBS)
	@echo "Built: $(BIN) (Zen 4 optimized, -march=znver4 -flto)"

debug: $(SRC)
	$(CC) $(CFLAGS_DBG) -o $(BIN_DBG) $< $(LIBS)
	@echo "Built: $(BIN_DBG) (debug + sanitizers)"

verify: $(BIN)
	./$(BIN) --verify

test: $(BIN)
	@rm -f a052075_checkpoint.txt
	./$(BIN) 2 600000
	@echo ""
	@echo "Expected: exactly 7 hits (a(1) through a(7))"

test8: $(BIN)
	@rm -f a052075_checkpoint.txt
	./$(BIN) 2 60000000
	@echo ""
	@echo "Expected: exactly 8 hits (a(1) through a(8))"

clean:
	rm -f $(BIN) $(BIN_DBG) a052075_checkpoint.txt a052075_log_*.csv

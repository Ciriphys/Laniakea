CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native
LIBS = -lm

LANIAKEA_DIR = Build/Laniakea
SOLITON_DIR  = Build/SolitonSim
BIN_DIR      = Bin

LANIAKEA_SRC_DIR = $(LANIAKEA_DIR)/Source
LANIAKEA_INC_DIR = $(LANIAKEA_DIR)/Include

SOLITON_SRC_DIR = $(SOLITON_DIR)/Source
SOLITON_INC_DIR = $(SOLITON_DIR)/Include

LANIAKEA_EXE = $(BIN_DIR)/Laniakea
SOLITON_EXE  = $(BIN_DIR)/SolitonSim

LANIAKEA_SRCS = \
	$(LANIAKEA_SRC_DIR)/uldm.cpp \
	$(LANIAKEA_SRC_DIR)/vec3.cpp \
	$(LANIAKEA_SRC_DIR)/solver.cpp \
	$(LANIAKEA_SRC_DIR)/parser.cpp \
	$(LANIAKEA_SRC_DIR)/ffmpeg.cpp

LANIAKEA_CXXFLAGS = $(CXXFLAGS) -fopenmp -DFFMPEG -I$(LANIAKEA_INC_DIR)
LANIAKEA_LDFLAGS  = -L/usr/local/lib
LANIAKEA_LIBS     = -lfftw3 -lfftw3_omp -fopenmp $(LIBS)

SOLITON_SRCS = \
	$(SOLITON_SRC_DIR)/soliton.cpp \
	$(SOLITON_SRC_DIR)/vec4.cpp \
	$(SOLITON_SRC_DIR)/integrator.cpp

SOLITON_CXXFLAGS = $(CXXFLAGS) -DTOL=1e-12 -I$(SOLITON_INC_DIR)
SOLITON_LIBS     = $(LIBS)

.PHONY: all clean Laniakea SolitonSim run-utils

all: Laniakea SolitonSim run-utils

Laniakea: $(LANIAKEA_EXE)

$(LANIAKEA_EXE): $(LANIAKEA_SRCS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(LANIAKEA_CXXFLAGS) $^ -o $@ \
	$(LANIAKEA_LDFLAGS) $(LANIAKEA_LIBS)

SolitonSim: $(SOLITON_EXE)

$(SOLITON_EXE): $(SOLITON_SRCS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(SOLITON_CXXFLAGS) $^ -o $@ \
	$(SOLITON_LIBS)

run-utils:
	printf '#!/bin/sh\nexec ./Bin/Laniakea "$$@"\n' > run-lnk
	printf '#!/bin/sh\nexec ./Bin/SolitonSim "$$@"\n' > run-sltn
	chmod +x run-lnk run-sltn

clean:
	rm -rf $(BIN_DIR) run-lnk run-sltn

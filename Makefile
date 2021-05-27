ROOT_DIR=$(shell git rev-parse --show-toplevel)
BIN_DIR=$(ROOT_DIR)/bin
SRC_DIR=$(ROOT_DIR)/src

INCLUDE_DIR=$(ROOT_DIR)/include
BATCH_DYNAMIC_DIR=$(INCLUDE_DIR)/batch_dynamic_connectivity
PARALLEL_EULER_TOUR_TREE_DIR=$(INCLUDE_DIR)/parallel_euler_tour_tree
SEQUENCE_DIR=$(INCLUDE_DIR)/sequence
UTIL_DIR=$(INCLUDE_DIR)/utilities
BENCHMARK_DIR=$(INCLUDE_DIR)/benchmark
PARLAY_DIR=$(ROOT_DIR)/external/parlay
GTEST_DIR=/usr/src/gtest

INC=$(INCLUDE_DIR) $(BATCH_DYNAMIC_DIR) $(PARALLEL_EULER_TOUR_TREE_DIR) $(SEQUENCE_DIR) $(UTIL_DIR) $(BENCHMARK_DIR) $(PARLAY_DIR) $(GTEST_DIR)
INC_PARAMS=$(INC:%=-I%)

CXX=/home/sualehasif/code/opencilk/opencilk-project/build/bin/clang++
CXXFLAGS=-std=c++17 -Wall -mcx16 $(INC_PARAMS) -I$(SRC_DIR) -MMD -MP
LDFLAGS=-fopencilk $(INC_PARAMS)
PARALLEL_FLAGS=-fopencilk -DPARLAY_CILK

TARGET=test_connectivity
ALL=$(TARGET)

OBJS=$(SRC_DIR)/$(TARGET).o \
	$(INCLUDE_DIR)/parallel_euler_tour_tree/src/edge_map.o \
	$(INCLUDE_DIR)/parallel_euler_tour_tree/src/euler_tour_sequence.o \
	$(INCLUDE_DIR)/parallel_euler_tour_tree/src/euler_tour_tree.o \
	$(INCLUDE_DIR)/sequence/parallel_skip_list/src/skip_list_base.o \

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(PARALLEL_FLAGS) -c -o $@ $<

-include $(TARGET).d

$(BIN_DIR)/$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

all: $(ALL)

.PHONY: clean

RM=rm

clean:
	$(RM) \
		$(OBJS) \
		$(patsubst %.o,%.d,$(OBJS)) \
        	$(BIN_DIR)/$(TARGET) \

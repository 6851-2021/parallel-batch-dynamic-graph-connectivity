ROOT_DIR=$(shell git rev-parse --show-toplevel)
BIN_DIR=$(ROOT_DIR)/bin
SRC_DIR=$(ROOT_DIR)/src
INCLUDE_DIR=$(ROOT_DIR)/include

CXX=/home/sualehasif/code/opencilk/opencilk-project/build/bin/clang++
CXXFLAGS=-std=c++14 -Wall -mcx16 -I$(INCLUDE_DIR) -I$(SRC_DIR) -MMD -MP
LDFLAGS=-fopencilk
PARALLEL_FLAGS=-fopencilk

TARGET=test_connectivity

OBJS=$(TARGET).o \
	$(INCLUDE_DIR)/parallel_euler_tour_tree/src/edge_map.o \
	$(INCLUDE_DIR)/parallel_euler_tour_tree/src/euler_tour_tree.o \
	$(INCLUDE_DIR)/sequence/parallel_skip_list/src/skip_list_base.o \

$(BIN_DIR)/$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(PARALLEL_FLAGS) -c -o $@ $<

-include $(TARGET).d

.PHONY: clean

RM=rm

clean:
	$(RM) \
		$(OBJS) \
		$(patsubst %.o,%.d,$(OBJS)) \
        	$(BIN_DIR)/$(TARGET) \

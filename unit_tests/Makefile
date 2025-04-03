BIN := CUModPTests
CXX := g++

SRCS := \
	test_two_convolution_poly_mul.cpp

OBJS := $(SRCS:.cpp=.o)

CXXFLAGS := \
	-Wall \
	-Werror \
	-std=c++20 \

INCLUDE_DIRS := \
	../include \

CPPFLAGS := \
	-DLINUXINTEL64 \
	-DPRINT_DEBUG
CPPFLAGS += $(addprefix -I, $(INCLUDE_DIRS))

LIB_DIRS := \
	../src

LDFLAGS += $(addprefix -L, $(LIB_DIRS))

LIB_NAMES := gtest gtest_main cumodp gmp gmpxx cudart
LDFLAGS += $(addprefix -l, $(LIB_NAMES))

debug: CXXFLAGS += -g
debug: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(BIN)

.PHONY: all clean
BIN := TwoConvPolyMulBench
CXX := g++

SRCS := \
	bench_two_convolution_poly_mul.cpp

OBJS := $(SRCS:.cpp=.o)

CXXFLAGS := \
	-Wall \
	-Werror \
	-std=c++14\
	-O2 \
	-g

INCLUDE_DIRS := \
	../include \

CPPFLAGS := \
	-DLINUXINTEL64 \

CPPFLAGS += $(addprefix -I, $(INCLUDE_DIRS))

LIB_DIRS := \
	../src

LDFLAGS += $(addprefix -L, $(LIB_DIRS))

LIB_NAMES := gmp gmpxx cudart cumodp 
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

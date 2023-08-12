
# ---------------------------------------------------------------------
#   Makefile for SLEMM 
#   
#   Supported platforms: Linux
# ---------------------------------------------------------------------

EFILE = slemm

# Intel oneAPI DPC++/C++ Compiler
CXX = icpx

# Compiler flags
CXXFLAGS = -Wall -O3 -qmkl -qopenmp -std=c++11
# Includes
# Change the EIGEN path before compiling.
INCLUDES = -I/home/jjiang26/software/eigen-3.3.9/

SRC_DIR := ./src
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
HDRS := $(wildcard $(SRC_DIR)/*.h)

OBJS = $(SRCS:%.cpp=%.o)

all : $(EFILE) 

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(EFILE) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EFILE) $(OBJS)

.PHONY: clean
clean: 
	rm -r $(EFILE) $(OBJS)


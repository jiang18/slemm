
# ---------------------------------------------------------------------
#   Makefile for SLEMM 
#   
#   Supported platforms: Linux
# ---------------------------------------------------------------------

# Directory of the target
OUTPUT = slemm

# Compiler
CXX = icpc

# Compiler flags
CXXFLAGS = -Wall -O3 -mkl -qopenmp -std=c++11 

# Includes
# Change the EIGEN path before compiling.
INCLUDES = -I/usr/include/eigen-3.3.9/

HDRS = StrFunc.h chisq.h gamma_inc.h reml.h snphwe.h slemm.h wls.h
SRCS = StrFunc.cpp chisq.cpp gamma_inc.cpp input.cpp option.cpp output.cpp pred.cpp reml.cpp snphwe.cpp weight.cpp wls.cpp
 
OBJS = $(SRCS:.cpp=.o)

all : $(OUTPUT) 

%.o: %.cpp ${HDRS}
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OUTPUT) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJS)

clean: 
	rm -f *.o


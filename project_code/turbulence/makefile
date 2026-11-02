# get all submodule objects and flags
SRCS = $(wildcard src/*.cpp)
OBJS = $(patsubst src/%.cpp, bin/%.o, $(SRCS))
FLAGS = -Wall -Wextra -pedantic -I../pv_coupling/src/ -I../face_addressed_mesh/src/

# compile all submodule objects
all: $(OBJS)

bin/tensor.o: src/tensor.cpp src/tensor.hpp \
								../face_addressed_mesh/src/point.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/turbulence2D%.o: src/turbulence2D%.cpp \
						src/turbulence2D.hpp src/tensor.hpp \
						../pv_coupling/src/pvCoupling2D.hpp
	g++ -g -c $< -o $@ $(FLAGS)
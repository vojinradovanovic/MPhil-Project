# get all submodule objects and flags
SRCS = $(wildcard src/*.cpp)
OBJS = $(patsubst src/%.cpp, bin/%.o, $(SRCS))
FLAGS = -Wall -Wextra -pedantic -I../face_addressed_mesh/src/ -I../pv_coupling/src
FACE_MESH_OBJS = $(wildcard ../face_addressed_mesh/bin/*.o)
PV_COUPLING_OBJS = $(wildcard ../pv_coupling/bin/*.o)

# compile all submodule objects
all: $(OBJS)

# rules for all files
bin/%.o: src/%.cpp ../pv_coupling/src/pvCoupling2D.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/%.exe: bin/%.o $(FACE_MESH_OBJS) $(PV_COUPLING_OBJS)
	g++ $^ -o $@
# define variables
FLAGS = -Wall -Wextra -pedantic -I./face_addressed_mesh/src/ -I./pv_coupling/src -I./turbulence/src
FACE_MESH_OBJS = $(wildcard ./face_addressed_mesh/bin/*.o)
PV_COUPLING_OBJS = $(wildcard ./pv_coupling/bin/*.o)
TURBULENCE_OBJS = $(wildcard ./turbulence/bin/*.o)

bin/main.o: main.cpp turbulence/src/turbulence2D.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/main.exe: bin/main.o $(FACE_MESH_OBJS) $(PV_COUPLING_OBJS) $(TURBULENCE_OBJS)
	g++ $^ -o $@
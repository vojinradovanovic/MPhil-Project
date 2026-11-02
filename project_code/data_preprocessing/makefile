# define variables
FLAGS = -Wall -Wextra -pedantic -I../face_addressed_mesh/src/ \
								-I../pv_coupling/src \
								-I../turbulence/src
FACE_MESH_OBJS = $(wildcard ../face_addressed_mesh/bin/*.o)
PV_COUPLING_OBJS = $(wildcard ../pv_coupling/bin/*.o)
TURBULENCE_OBJS = $(wildcard ../turbulence/bin/*.o)

bin/helpers.o: src/helpers.cpp src/helpers.hpp ../turbulence/src/turbulence2D.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/write_output_data.o: src/write_output_data.cpp src/helpers.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/write_output_data.exe: bin/write_output_data.o bin/helpers.o $(FACE_MESH_OBJS) $(PV_COUPLING_OBJS) $(TURBULENCE_OBJS)
	g++ $^ -o $@

bin/compare_data.o: src/compare_data.cpp src/helpers.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/compare_data.exe: bin/compare_data.o bin/helpers.o $(FACE_MESH_OBJS) $(PV_COUPLING_OBJS) $(TURBULENCE_OBJS)
	g++ $^ -o $@
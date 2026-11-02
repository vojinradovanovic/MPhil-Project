# get all submodule objects and flags
SRCS = $(wildcard src/*.cpp)
OBJS = $(patsubst src/%.cpp, bin/%.o, $(SRCS))
FLAGS = -Wall -Wextra -pedantic -I../face_addressed_mesh/src/

# compile all submodule objects
all: $(OBJS)

# rules for each .o file
bin/compressed_row_matrix.o: src/compressed_row_matrix.cpp \
							src/compressed_row_matrix.hpp src/physics_constants.hpp ../face_addressed_mesh/src/constants.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/linear_system.o: src/linear_system.cpp \
						src/linear_system.hpp src/compressed_row_matrix.hpp src/physics_constants.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/pvCoupling2D_%.o: src/pvCoupling2D_%.cpp \
						src/pvCoupling2D.hpp src/compressed_row_matrix.hpp ../face_addressed_mesh/src/face_mesh.hpp \
						src/linear_system.hpp src/physics_constants.hpp
	g++ -g -c $< -o $@ $(FLAGS)
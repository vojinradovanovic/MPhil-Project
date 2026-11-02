# get all submodule objects
SRCS = $(wildcard src/*.cpp)
OBJS = $(patsubst src/%.cpp, bin/%.o, $(SRCS))

# compile all submodule objects
all: $(OBJS)

# rules for each .o file
bin/point.o: src/point.cpp src/point.hpp src/constants.hpp
	g++ -g -c $< -o $@ -Wall -Wextra -pedantic

bin/face.o: src/face.cpp src/face.hpp src/point.hpp src/constants.hpp
	g++ -g -c $< -o $@ -Wall -Wextra -pedantic

bin/cell.o: src/cell.cpp src/cell.hpp src/point.hpp src/constants.hpp
	g++ -g -c $< -o $@ -Wall -Wextra -pedantic

bin/bdr_patch.o: src/bdr_patch.cpp src/bdr_patch.hpp
	g++ -g -c $< -o $@ -Wall -Wextra -pedantic

bin/face_mesh.o: src/face_mesh.cpp src/face_mesh.hpp src/bdr_patch.hpp src/cell.hpp src/face.hpp src/point.hpp src/constants.hpp
	g++ -g -c $< -o $@ -Wall -Wextra -pedantic
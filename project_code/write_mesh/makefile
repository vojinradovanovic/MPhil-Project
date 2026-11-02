# define variables
FLAGS = -Wall -Wextra -pedantic

bin/write.o: src/write.cpp src/write.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/main.o: src/main.cpp src/write.hpp
	g++ -g -c $< -o $@ $(FLAGS)

bin/main.exe: bin/main.o bin/write.o
	g++ $^ -o $@
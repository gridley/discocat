INC_FLAGS = -lstdc++fs -I/usr/include/eigen3
FLAGS_DBG = -g -O0 -fstack-protector -Wall $(INC_FLAGS)
FLAGS_OPT = -O4 -ffast-math $(INC_FLAGS)
FLAGS = $(FLAGS_OPT)
main: discocat.C lsquadrature
	g++ discocat.C lsquadrature.o $(FLAGS)
lsquadrature: lsquadrature.C lsquadrature.h
	g++ -c lsquadrature.C $(FLAGS)

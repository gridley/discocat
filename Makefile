FLAGS_DBG = -g -O0 -lstdc++fs -fstack-protector
FLAGS_OPT = -O4 -ffast-math -lstdc++fs
FLAGS = $(FLAGS_OPT)
main: discocat.C lsquadrature
	g++ discocat.C lsquadrature.o $(FLAGS)
lsquadrature: lsquadrature.C lsquadrature.h
	g++ -c lsquadrature.C $(FLAGS)

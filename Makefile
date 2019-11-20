main: arrt.C
	# g++ discocat.C -g -O0 -lgd -lstdc++fs -fstack-protector
	g++ discocat.C -O4 -ffast-math -lgd -lstdc++fs

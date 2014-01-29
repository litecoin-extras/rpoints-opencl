all:
	g++ -O3 -march=native -std=c++0x -DUSE_GMP -D_LINUX main.cpp targets.cpp -lgmp -lclsocket -lpthread -o rpoints

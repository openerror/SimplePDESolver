ALL: MFSolve

# ARMA_LIB="/usr/local/Cellar/armadillo/8.400.0/lib"
# ARMA_INC="/usr/local/Cellar/armadillo/8.400.0/include"

MFSolve: ./reactions.cpp ./main.cpp
	g++ -O3 -std=c++11 reactions.cpp main.cpp -larmadillo -o MFSolve

clean:
	rm MFSolve

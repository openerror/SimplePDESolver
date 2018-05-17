ALL: MFSolve

ARMA_LIB="/Users/everestl/CustomG++Libs/lib"
ARMA_INC="/Users/everestl/CustomG++Libs/include"

MFSolve: ./reactions.cpp ./main.cpp
	g++-7 -O3 -fopenmp -L $(ARMA_LIB) -I $(ARMA_INC)\
	 -std=c++11 reactions.cpp main.cpp -larmadillo -o MFSolve

clean:
	rm MFSolve
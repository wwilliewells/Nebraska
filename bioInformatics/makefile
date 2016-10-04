all: neighborJoin

INCLUDE = -I/include

neighborJoin: joinNeighbor.o matrix.o
	g++ -o neighborJoin joinNeighbor.o matrix.o 

joinNeighbor.o: joinNeighbor.cpp
	g++ -c joinNeighbor.cpp $(INCLUDE)

matrix.o: matrix.cpp
	g++ -c matrix.cpp -I $(INCLUDE)

clean:
	rm *.o neighborJoin

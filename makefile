all: GWOC

GWOC: RngStream.o GWO-C.o
	g++ RngStream.o GWO-C.o -o GWOC

RngStream.o: RngStream.c
	g++ -c RngStream.c

GWO-C.o: GWO-C.cpp
	g++ -c GWO-C.cpp

clean:
	rm -rf *o GWOC
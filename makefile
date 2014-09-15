all: GWOC

GWOC: RngStream.o GWO-C.o
	g++ RngStream.o GWO-C.o -o GWOC

RngStream.o: RngStream.c
	g++ -c RngStream.c

GWO-C.o: GWO-C.c
	g++ -c GWO-C.c

clean:
	rm -rf *o GWOC
CC = g++    
CFLAGS = -O3 -std=c++11

all: difsimplificator

difsimplificator: simplificator.cpp io.o EasyBMP.o
	$(CC) $(CFLAGS) simplificator.cpp EasyBMP.o io.o  -o difsimplificator

io.o: include/io.cpp
	$(CC) $(CFLAGS) -c include/io.cpp


EasyBMP.o: externals/EasyBMP.cpp
	$(CC) $(CFLAGS) -c externals/EasyBMP.cpp

clean:
	rm -f *.o
	rm -f median
	rm -f customImage
	rm -f *.bmp

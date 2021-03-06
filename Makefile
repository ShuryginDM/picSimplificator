CC = g++    
CFLAGS = -O3 -std=c++11

all: difsimplificator

difsimplificator: main.cpp simplificator.o 
	$(CC) $(CFLAGS) main.cpp simplificator.o EasyBMP.o io.o -o difsimplificator

simplificator.o: src/simplificator.cpp io.o EasyBMP.o
	$(CC) $(CFLAGS) src/simplificator.cpp EasyBMP.o io.o  -c 

io.o: include/io.cpp
	$(CC) $(CFLAGS) -c include/io.cpp


EasyBMP.o: externals/EasyBMP.cpp
	$(CC) $(CFLAGS) -c externals/EasyBMP.cpp

clean:
	rm -f *.o
	rm -f median
	rm -f customImage
	rm -f *.bmp

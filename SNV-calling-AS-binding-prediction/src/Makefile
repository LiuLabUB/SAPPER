CC = g++ 
CFLAGS = -g -Wall -O3
OBJS = main.o io.o analysis.o swalign.o math.o

SNVAS : $(OBJS)
	$(CC) $(CFLAGS) -lz -o SNVAS $(OBJS)

math.o : math.cpp math.hpp
	$(CC) -c $(CFLAGS) math.cpp

swalign.o : swalign.cpp swalign.hpp
	$(CC) -c $(CFLAGS) swalign.cpp

analysis.o : analysis.cpp analysis.hpp
	$(CC) -c $(CFLAGS) analysis.cpp

io.o : io.cpp io.hpp swalign.hpp analysis.hpp
	$(CC) -c $(CFLAGS) io.cpp

main.o : main.cpp io.hpp
	$(CC) -c $(CFLAGS) main.cpp

clean :
	-rm -f SNVAS $(OBJS)

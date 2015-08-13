CC = g++ -c -O3
OBJS = main.o io.o analysis.o swalign.o math.o


SNVAS : $(OBJS)
	g++ -O3 -lz -o SNVAS $(OBJS)

math.o : math.cpp math.hpp
	$(CC) math.cpp
	
swalign.o : swalign.cpp swalign.hpp
	$(CC) swalign.cpp

analysis.o : analysis.cpp analysis.hpp
	$(CC) analysis.cpp

io.o : io.cpp io.hpp swalign.hpp analysis.hpp
	$(CC) io.cpp

main.o : main.cpp io.hpp
	$(CC) main.cpp

clean :
	-rm -f SNVAS $(OBJS)


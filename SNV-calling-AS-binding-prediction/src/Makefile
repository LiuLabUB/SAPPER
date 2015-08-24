MAKE = make
MKDIR = mkdir
CC = g++ 
CFLAGS = -g -Wall -O3
OBJS = main.o io.o analysis.o swalign.o math.o

.PHONY:all

all: bin SNVAS SNVAS_filter SNVAS_cutoff_statistics
	cd fermi && ${MAKE}
	mv fermi/SNVAS_fermi bin/ #we renamed the name `SNVAS_fermi` because we modifed the codes.
	mv SNVAS SNVAS_filter bin/
	mv SNVAS_cutoff_statistics bin/

bin:
	${MKDIR} bin	

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

SNVAS_filter : SNVAS_filter.cpp
	$(CC) $(CFLAGS) -lz -o SNVAS_filter SNVAS_filter.cpp
	
SNVAS_cutoff_statistics : SNVAS_cutoff_statistics.cpp
	$(CC) $(CFLAGS) -lz -o SNVAS_cutoff_statistics SNVAS_cutoff_statistics.cpp

clean :
	cd fermi && ${MAKE} clean
	-rm -f $(OBJS)
	-rm -f bin/SNVAS_fermi bin/SNVAS bin/SNVAS_filter bin/SNVAS_cutoff_statistics
	-rmdir bin


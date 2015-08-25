MYNAME = SNVAS
MAKE = make
MKDIR = mkdir
CC = g++ 
CFLAGS = -g -Wall -O3
OBJS = main.o io.o analysis.o swalign.o math.o

.PHONY:all

all: bin ${MYNAME} ${MYNAME}_filter ${MYNAME}_stat
	cd fermi && ${MAKE}
	mv fermi/${MYNAME}_fermi bin/ #we renamed the name `${MYNAME}_fermi` because we modifed the codes.
	mv ${MYNAME} ${MYNAME}_filter bin/
	mv ${MYNAME}_stat bin/

bin:
	${MKDIR} bin	

${MYNAME} : $(OBJS)
	$(CC) $(CFLAGS) -lz -o ${MYNAME} $(OBJS)

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

${MYNAME}_filter : ${MYNAME}_filter.cpp
	$(CC) $(CFLAGS) -lz -o ${MYNAME}_filter SNVAS_filter.cpp

${MYNAME}_stat : SNVAS_cutoff_statistics.cpp
	$(CC) $(CFLAGS) -lz -o ${MYNAME}_stat SNVAS_cutoff_statistics.cpp

clean :
	cd fermi && ${MAKE} clean
	-rm -f $(OBJS)
	-rm -f bin/${MYNAME}_fermi bin/${MYNAME} bin/${MYNAME}_filter bin/${MYNAME}_cutoff_statistics
	-rmdir bin


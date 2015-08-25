MYNAME = SNVAS
MAKE = make
MKDIR = mkdir
CC = g++ 
CFLAGS = -g -Wall -O3
OBJS = io.o analysis.o swalign.o math.o call.o filter.o stat.o SNVAS_common.o

.PHONY:all

all: bin ${MYNAME} 
	cd fermi && ${MAKE}
	mv fermi/${MYNAME}_fermi bin/ #we renamed the name `${MYNAME}_fermi` because we modifed the codes.
	mv ${MYNAME} bin/

bin:
	${MKDIR} bin	

${MYNAME} : $(OBJS)
	$(CC) $(CFLAGS) -lz -o ${MYNAME} $(OBJS) main.cpp

math.o : math.cpp math.hpp
	$(CC) -c $(CFLAGS) math.cpp

swalign.o : swalign.cpp swalign.hpp
	$(CC) -c $(CFLAGS) swalign.cpp

analysis.o : analysis.cpp analysis.hpp
	$(CC) -c $(CFLAGS) analysis.cpp

io.o : io.cpp io.hpp swalign.hpp analysis.hpp
	$(CC) -c $(CFLAGS) io.cpp

call.o : call.cpp io.hpp
	$(CC) -c $(CFLAGS) call.cpp

filter.o : filter.cpp SNVAS_common.hpp
	$(CC) -c $(CFLAGS) filter.cpp

stat.o : stat.cpp SNVAS_common.hpp
	$(CC) -c $(CFLAGS) stat.cpp

SNVAS_common.o : SNVAS_common.cpp SNVAS_common.hpp
	$(CC) -c $(CFLAGS) SNVAS_common.cpp

clean :
	cd fermi && ${MAKE} clean
	-rm -f $(OBJS)
	-rm -f bin/${MYNAME}_fermi bin/${MYNAME} bin/* 
	-rmdir bin


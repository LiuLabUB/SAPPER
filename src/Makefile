APPNAME = SAPPER
MAKE = make
MKDIR = mkdir
CPPC = g++ 
CC = gcc
CPPFLAGS = -g -Wall -O3 -ffast-math -std=c++11 -Wsign-compare
CPPOBJS = io.o analysis.o swalign.o math.o call.o filter.o stat.o SAPPER_common.o
LIBS=	-lpthread -lm -lz

.PHONY:all

all: ${APPNAME}

${APPNAME}:$(CPPOBJS)
	$(MAKE) -C fermi/
	$(CPPC) $(CPPFLAGS) $(LIBS) -o $(APPNAME) $(CPPOBJS) main.cpp fermi/fermi.a

math.o : math.cpp math.hpp
	$(CPPC) -c $(CPPFLAGS) math.cpp

swalign.o : swalign.cpp swalign.hpp
	$(CPPC) -c $(CPPFLAGS) swalign.cpp

analysis.o : analysis.cpp analysis.hpp
	$(CPPC) -c $(CPPFLAGS) analysis.cpp

io.o : io.cpp io.hpp swalign.hpp analysis.hpp
	$(CPPC) -c $(CPPFLAGS) io.cpp

call.o : call.cpp io.hpp
	$(CPPC) -c $(CPPFLAGS) call.cpp

filter.o : filter.cpp SAPPER_common.hpp
	$(CPPC) -c $(CPPFLAGS) filter.cpp

stat.o : stat.cpp SAPPER_common.hpp
	$(CPPC) -c $(CPPFLAGS) stat.cpp

SAPPER_common.o : SAPPER_common.cpp SAPPER_common.hpp
	$(CPPC) -c $(CPPFLAGS) SAPPER_common.cpp

clean :
	$(MAKE) -C fermi/ clean
	-rm -f $(CPPOBJS)
	-rm -f $(APPNAME)


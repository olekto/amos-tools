LIBS = -lCommon -lAMOS -lz -lc -lm
CXX=g++
CC=g++
CXXFLAGS=-Wall -Wno-missing-declarations -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
PROG=extractScaffold
INCDIR=-I/usit/titan/u1/olekto/bin/amos_git/include/AMOS
LIBDIR=-L/usit/titan/u1/olekto/bin/amos_git/lib/AMOS
all : $(PROG)

OBJS = extractScaffold.o

#DEPS = Makefile

extractScaffold.o : extractScaffold.cc 
	$(CC) $(CXXFLAGS) $(INCDIR) -c extractScaffold.cc
#	$(CC) $(CXXFLAGS) $(DEPS) $(INCDIR) -c extractScaffold.cc

$(PROG) : $(OBJS)
	$(CC) $(CXXFLAGS)  -o $(PROG) $(OBJS) $(LIBS) $(INCDIR) $(LIBDIR)

clean :
	rm -f *.o *~ $(PROG)


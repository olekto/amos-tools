LIBS = -lCommon -lAMOS -lz -lc -lm
CXX=g++
CC=g++
CXXFLAGS=-Wall -Wno-missing-declarations -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
PROG=extractScaffold bank2sam
#You will have to change these paths to something that suits your system
INCDIR=-I/projects/454data/bin/amos_git/include/AMOS
LIBDIR=-L/projects/454data/bin/amos_git/lib/AMOS
all : $(PROG)

OBJS = extractScaffold.o bank2sam.o

bank2sam.o : bank2sam.cc 
	$(CC) $(CXXFLAGS) $(INCDIR) -c bank2sam.cc

extractScaffold.o : extractScaffold.cc 
	$(CC) $(CXXFLAGS) $(INCDIR) -c extractScaffold.cc

extractScaffold :
	$(CC) $(CXXFLAGS)  -o extractScaffold extractScaffold.o  $(LIBS) $(INCDIR) $(LIBDIR)

bank2sam :
	$(CC) $(CXXFLAGS)  -o bank2sam bank2sam.o  $(LIBS) $(INCDIR) $(LIBDIR)

$(PROG) : $(OBJS) 

clean :
	rm -f *.o *~ $(PROG)


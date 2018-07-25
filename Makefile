CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

CPPFLAGS := $(shell root-config --cflags) $(STDINCDIR)
LDFLAGS := -Xlinker -rpath . $(shell root-config --glibs) $(STDLIBDIR)

#CPPFLAGS += -g
#CPPFLAGS += -O1

#TARGET = transparency

#SRC = transparency.C

#OBJ = $(SRC:.C=.o)

#all : plots

#plots : plots.o
#	$(LD) $(CPPFLAGS) -o plots plots.o $(LDFLAGS)

#transparency : transparency.o
#	$(LD) $(CPPFLAGS) -o transparency transparency.o $(LDFLAGS)


#%.o : %.C
#	$(CXX) $(CPPFLAGS) -o $@ -c $<

#clean :
#	rm -f *.o $(TARGET) *~

#CXX = $(shell root-config --cxx)
#LD = $(shell root-config --ld)

#CPPFLAGS := $(shell root-config --cflags) $(STDINCDIR)
#LDFLAGS := -Xlinker -rpath . $(shell root-config --glibs) $(STDLIBDIR)

#CPPFLAGS += -g
#CPPFLAGS += -O1

#TARGET = transparency

#SRC = transparency.C

#OBJ = $(SRC:.C=.o)

#all : plots

#plots : plots.o
#	$(LD) $(CPPFLAGS) -o plots plots.o $(LDFLAGS)

#transparency : transparency.o
#	$(LD) $(CPPFLAGS) -o transparency transparency.o $(LDFLAGS)


#%.o : %.C
#	$(CXX) $(CPPFLAGS) -o $@ -c $<

#clean :
#	rm -f *.o $(TARGET) *~
#pentaquarkMatrix.cpp Isobar.cpp main.cpp pentaquarkMatrix.h Isobar.h

main : Isobar.cpp Isobar.h main.cpp PentaquarkMatrix.cpp PentaquarkMatrix.h SpecialFunctions.h
	g++ -o main Isobar.cpp main.cpp  PentaquarkMatrix.cpp  -I/usr/include/gsl/  -lgsl -lgslcblas -lm

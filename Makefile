

CXX_PROG=zos_getentropy

# *******************************************************
# ************************ USS **************************
# *******************************************************

CXXFLAGS+=-O0 -std=gnu++20 -fzos-le-char-mode=ascii $(INCLUDE)
CC=ibm-clang64
CXX=ibm-clang++64
LD=ibm-clang++64

INCLUDE+=-I. ../include -I/usr/local/include -I/usr/include

CXXSOURCES = zos_getentropy.cpp 
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)  


$(CXX_PROG): $(CXXOBJECTS)
	$(CXX) -o $@ $(CXXOBJECTS) 

all:: clean $(CXX_PROG) 
	
clean::
	$(RM) $(CXXOBJECTS) $(CXX_PROG) 

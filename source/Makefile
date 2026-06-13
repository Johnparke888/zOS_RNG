CXX_PROG=zos_getentropy

# *******************************************************
# ************************ USS **************************
# *******************************************************

RM=rm -f

CC=ibm-clang64
CXX=ibm-clang++64
LD=ibm-clang++64

INCLUDE+=-I. -I../include -I/usr/local/include -I/usr/include

CXXFLAGS+=-O0 -std=gnu++23 -fzos-le-char-mode=ascii -mnocsect -m64 \
-fno-short-enums -faligned-allocation $(INCLUDE)

CFLAGS+=-O0 -std=c17 -fzos-le-char-mode=ascii -mnocsect -m64 $(INCLUDE)

CXXSOURCES = zos_entropy.cpp zos_statistical.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)

CSOURCES = zos_generators.c
COBJECTS = $(CSOURCES:.c=.o)

all:: clean $(CXX_PROG)

$(CXX_PROG): $(CXXOBJECTS) $(COBJECTS)
	$(LD) -o $@ $(CXXOBJECTS) $(COBJECTS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean::
	$(RM) $(CXXOBJECTS) $(COBJECTS) $(CXX_PROG)

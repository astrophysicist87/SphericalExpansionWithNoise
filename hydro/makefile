SHELL=/bin/sh

SRCS= \
main.cpp \
gauss_quadrature.cpp \
lib.cpp \
legendre.cpp

HDRS= \
gauss_quadrature.h \
defs1.h \
defs2.h \
HBT_and_emission_density.h \
Phi_functions_and_derivatives.h \
Psi_functions_and_derivatives.h \
lib.h \
legendre.h \
Stopwatch.h

MAKEFILE=makefile

COMMAND=run.e

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
#CFLAGS=  -pg -g -O3
CFLAGS=  -g -O3
WARNFLAGS= -ansi -pedantic -Wall -W \
   -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
   -Wwrite-strings -fshort-enums -fno-common 
LDFLAGS= -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include

$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(CFLAGS) $(LIBS)

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 

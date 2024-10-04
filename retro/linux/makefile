#
# cwt -- tiny complex (analitic) audio wave conwerter
# posix / gcc / FFT(W) only
# (c) 2014 Rat and Catcher Tech.
#

PROGRAM     = cwt

# installation prefix of FFTW
FFTW_PREFIX = $(HOME)/opt

CC          = gcc
LD          = gcc
CCOPT       = -Wall -Winline -O3 -malign-double
CCINC       = -I$(FFTW_PREFIX)/include
CCLIB       = -L$(FFTW_PREFIX)/lib64

LIBRARIES   = -l:libfftw3_threads.a -l:libfftw3.a -lpthread -lm

# own includes
INCL        = cwt.h cwave.h crc32.h
# own objects
OBJS        = cwt.o cw_fft.o helpers.o crc32.o
# 

# other stuff
.SUFFIXES:
.SUFFIXES: .c .h .o

default:
	@echo
	@echo "Usage: make all | $(PROGRAM) | install"
	@echo

# the targets

all: $(PROGRAM)

install: $(PROGRAM)
	rm -f $(HOME)/bin/$(PROGRAM)
	mv $(PROGRAM) $(HOME)/bin
	rm -f *.o

$(PROGRAM): $(OBJS)
	$(LD) $(OBJS) -o $(PROGRAM) $(CCLIB) $(LIBRARIES)

crc32.o: crc32.c crc32.h
	$(CC) -c crc32.c $(CCOPT) $(CCINC)

helpers.o: helpers.c $(INCL)
	$(CC) -c helpers.c $(CCOPT) $(CCINC)

cw_fft.o: cw_fft.c $(INCL)
	$(CC) -c cw_fft.c $(CCOPT) $(CCINC)

cwt.o: cwt.c $(INCL)
	$(CC) -c cwt.c $(CCOPT) $(CCINC)

clean:
	-rm -f *.o
	-rm -f $(PROGRAM)

.PHONY: default all

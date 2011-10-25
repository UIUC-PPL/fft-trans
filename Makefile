include ../config.mk

CC = mpicxx
OPTS = -O3
FFTW3DIR = /expand/home/arya/fftw-3.3

INC = -I$(FFTW3DIR)/include
LIBDIR = -L$(FFTW3DIR)/lib
LIBS = -lfftw3 -lm

REFINC = -I$(FFTW3DIR)/include
REFLIBDIR = -L$(FFTW3DIR)/lib
REFLIBS = -lfftw3_mpi -lfftw3 -lm

all: fft1d fft_ref

test: all
	./test.sh 4 64

fft1d: fft1d.o
	$(CHARMC) $(OPTS) -language charm++ -o fft1d fft1d.o $(LIBDIR) $(LIBS)

fft1d.o: fft1d.cc fft1d.decl.h fft1d.def.h
	$(CHARMC) $(OPTS) -c fft1d.cc $(INC)

fft1d.decl.h fft1d.def.h: fft1d.ci
	$(CHARMC) $(OPTS) fft1d.ci

fft_ref: fft_ref.o
	${CC} $(OPTS) fft_ref.o -o fft_ref $(REFLIBDIR) $(REFLIBS)

fft_ref.o: fft_ref.cc
	${CC} $(OPTS) -c fft_ref.cc $(REFINC)

clean:
	rm -f *dump* *.decl.h *.def.h conv-host *.o *~

clobber: clean
	rm -f fft1d fft_ref charmrun
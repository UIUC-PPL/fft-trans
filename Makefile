OPTS	= -O3

CHARMC	= $(HOME)/charm/bin/charmc $(OPTS)
CC      = mpicxx $(OPTS)

INC =
LIBDIR =
LIBS = -lfftw3 -lm

REFINC = -I/expand/home/arya/fftw-3.3/include
REFLIBDIR = -L/expand/home/arya/fftw-3.3/lib
REFLIBS = -lfftw3_mpi -lfftw3 -lm

all: fft1d fft_ref

fft1d: fft1d.o
	$(CHARMC) -language charm++ -o fft1d fft1d.o $(LIBDIR) $(LIBS)

fft1d.o: fft1d.cc fft1d.decl.h fft1d.def.h
	$(CHARMC) -c fft1d.cc $(INC)

fft1d.decl.h fft1d.def.h: fft1d.ci
	$(CHARMC) fft1d.ci

fft_ref: fft_ref.o
	${CC} fft_ref.o -o fft_ref $(REFLIBDIR) $(REFLIBS)

fft_ref.o: fft_ref.cc
	${CC} -c fft_ref.cc $(REFINC)

clean:
	rm -f *.decl.h *.def.h conv-host *.o fft1d fft_ref charmrun *~

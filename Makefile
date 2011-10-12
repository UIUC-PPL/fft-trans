OPTS	= -O3
CHARMC	= $(HOME)/charm/bin/charmc $(OPTS)
CC=mpicxx

OBJS = fft1d.o

all: fft1d fft_ref

fft1d: $(OBJS)
	$(CHARMC) -language charm++ -o fft1d $(OBJS) -lfftw3 -lm

projections: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections -lz -o fft1d.prj $(OBJS)

summary: $(OBJS)
	$(CHARMC) -language charm++ -tracemode summary -lz -o fft1d.sum $(OBJS)

fft1d.decl.h: fft1d.ci
	$(CHARMC)  fft1d.ci

fft_ref: fft_ref.o
	${CC} fft_ref.o -o fft_ref -L/expand/home/arya/fftw-3.3/lib -lfftw3_mpi -lfftw3 -lm

fft_ref.o: fft_ref.cpp
	${CC} -c fft_ref.cpp -I/expand/home/arya/fftw-3.3/include

clean:
	rm -f *.decl.h *.def.h conv-host *.o fft1d fft1d.prj charmrun fft_ref *~

fft1d.o: fft1d.C fft1d.decl.h
	$(CHARMC) -c fft1d.C

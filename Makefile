OPTS	= -O3
CHARMC	= $(HOME)/charm-production/bin/charmc $(OPTS)
CC=mpicxx

FFTWPATH = /soft/apps/fftw-3.1.2-double
INC = -I$(FFTWPATH)/include
LIBS = -L$(FFTWPATH)/lib -lfftw3 -lm
CHARMLIBS = -module NDMeshStreamer -module completion

OBJS = fft1d.o

all: fft1d projections fft_bench

fft_bench: fft_bench.o
	${CC} fft_bench.o -o fft_bench $(LIBS)

fft_bench.o: fft_bench.cpp
	${CC} -c fft_bench.cpp $(INC)

fft1d: $(OBJS)
	$(CHARMC) -language charm++ -o fft1d $(OBJS) $(LIBS) $(CHARMLIBS)

projections: fft1d.prj fft1d.sum
fft1d.prj: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections $(LIBS) $(CHARMLIBS) -lz -o fft1d.prj $(OBJS)

fft1d.sum: $(OBJS)
	$(CHARMC) -language charm++ -tracemode summary $(LIBS) $(CHARMLIBS) -o fft1d.sum $(OBJS)

fft1d.decl.h: fft1d.ci
	$(CHARMC)  fft1d.ci

cleanproj:
	rm -f *.log *.sts *.projrc

clean:
	rm -f *.decl.h *.def.h conv-host *.o fft1d fft1d.prj fft1d.sum fft_bench charmrun *~

fft1d.o: fft1d.cc fft1d.decl.h
	$(CHARMC) -c fft1d.cc $(INC)

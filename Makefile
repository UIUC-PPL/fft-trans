OPTS	= -O3
CHARMC	?= $(HOME)/charm-debug/bin/charmc $(OPTS)
CC=mpicxx
LIBS = -lfftw3 -lm
CHARMLIBS = -module NDMeshStreamer -module completion
FFTW3 ?= $(HOME)/fftw-3.3

OBJS = fft1d.o

all: fft1d

fft1d: $(OBJS)
	$(CHARMC) -language charm++ -o fft1d $(OBJS) $(LIBS) $(CHARMLIBS)

projections: fft1d.prj fft1d.sum
fft1d.prj: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections $(LIBS) $(CHARMLIBS) -lz -o fft1d.prj $(OBJS)

fft1d.sum: $(OBJS)
	$(CHARMC) -language charm++ -tracemode summary $(LIBS) $(CHARMLIBS) -o fft1d.sum $(OBJS)

fft1d.decl.h: fft1d.ci
	$(CHARMC)  fft1d.ci

fft.decl.h: fft.ci
	$(CHARMC) fft.ci

cleanproj:
	rm -f *.log *.sts *.projrc

clean:
	rm -f *.decl.h *.def.h conv-host *.o fft1d fft1d.prj fft1d.sum charmrun *~

fft1d.o: fft1d.cc fft1d.decl.h fft.decl.h
	$(CHARMC) -c fft1d.cc

OPTS	= -O3
CHARMC	?= $(HOME)/charm-debug/bin/charmc $(OPTS)
CC=mpicxx
LIBS = -lfftw3 -lm
CHARMLIBS = -module NDMeshStreamer -module completion
FFTW3 ?= $(HOME)/fftw-3.3

OBJS = main.o

all: main

main: $(OBJS)
	$(CHARMC) -language charm++ -o main $(OBJS) $(LIBS) $(CHARMLIBS)

projections: main.prj main.sum
main.prj: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections $(LIBS) $(CHARMLIBS) -lz -o main.prj $(OBJS)

main.sum: $(OBJS)
	$(CHARMC) -language charm++ -tracemode summary $(LIBS) $(CHARMLIBS) -o main.sum $(OBJS)

main.decl.h: main.ci
	$(CHARMC)  main.ci

fft.decl.h: fft.ci
	$(CHARMC) fft.ci

cleanproj:
	rm -f *.log *.sts *.projrc

clean:
	rm -f *.decl.h *.def.h conv-host *.o main main.prj main.sum charmrun *~

main.o: main.cc main.decl.h fft.decl.h
	$(CHARMC) -c main.cc

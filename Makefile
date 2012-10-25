OPTS	= -O3
CHARMHOME ?= $(HOME)/charms/charm/mpi-linux-x86_64
CHARMC  = $(CHARMHOME)/bin/charmc $(OPTS)
FFTW3   ?= $(HOME)/fftw-3.3
INC     = -I$(FFTW3)/include
LIBS    = -L$(FFTW3)/lib -lfftw3 -lm
CHARMLIBS = -module NDMeshStreamer -module completion

OBJS = main.o

all: main

main: $(OBJS)
	$(CHARMC) -language charm++ -o main $(OBJS) $(LIBS) $(CHARMLIBS)

projections: main.prj main.sum
main.prj: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections $(LIBS) $(CHARMLIBS) -lz -o main.prj $(OBJS)

main.sum: $(OBJS)
	$(CHARMC) -language charm++ -tracemode summary $(LIBS) $(CHARMLIBS) -o main.sum $(OBJS)

main.o: main.cc main.decl.h fft.decl.h fftData.decl.h
	$(CHARMC) -c main.cc $(INC)

main.decl.h fftData.decl.h: main.ci
	$(CHARMC)  main.ci

fft.decl.h: fft.ci
	$(CHARMC) fft.ci

test: main
	./test.sh 2 512

cleanproj:
	rm -f *.log *.sts *.projrc

clean:
	rm -f *.decl.h *.def.h conv-host *.o main main.prj main.sum charmrun *~

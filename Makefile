OPTS	= -O3
CHARMC	?= $(HOME)/charm/bin/charmc $(OPTS)
LIBS    = -lfftw3 -lm
CHARMLIBS = -module NDMeshStreamer -module completion
FFTW3   ?= $(HOME)/fftw-3.3

all: main

main: main.o
	$(CHARMC) -language charm++ -o main main.o $(LIBS) $(CHARMLIBS)

main.o: main.cc main.decl.h fft.decl.h fftData.decl.h
	$(CHARMC) -c main.cc

main.decl.h fftData.decl.h: main.ci
	$(CHARMC)  main.ci

fft.decl.h: fft.ci
	$(CHARMC) fft.ci

test: main
	./test.sh 2 512

cleanproj:
	rm -f *.log *.sts *.projrc

clean:
	rm -f *.decl.h *.def.h conv-host *.o main charmrun *~

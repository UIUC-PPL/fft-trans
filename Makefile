OPTS	= -O3
CHARMHOME ?= $(HOME)/charm
CHARMC  = $(CHARMHOME)/bin/charmc $(OPTS)

ESSLPATH = /soft/apps/ESSL-4.4
WRAPPATH = $(HOME)/esslfftw
INC = -I$(WRAPPATH)/include -I$(ESSLPATH)/include
BGP_LIBS = -L$(ESSLPATH)/lib \
           -L$(WRAPPATH)/lib \
           -L/bgsys/ibm_compilers/sles10/prod/opt/ibmcmp/xlf/bg/11.1/bglib/ \
           -L/soft/apps/ibmcmp/xlsmp/bg/1.7/bglib \
           -L/soft/apps/ibmcmp/xlf/bg/11.1/bglib
LIBS = $(BGP_LIBS) -lfftw3_esslbg -lesslbg -lmass -lxlfmath -lxlf90_r -lxlsmp -lxlomp_ser -lpthread
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

cleanproj:
	rm -f *.log *.sts *.projrc

clean:
	rm -f *.decl.h *.def.h conv-host *.o main main.prj main.sum charmrun *~

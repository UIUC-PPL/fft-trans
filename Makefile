OPTS	= -O3
CHARMC	= $(HOME)/charm-production/bin/charmc $(OPTS)
CC=mpixlcxx $(OPTS)

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

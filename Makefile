OPTS	= -O3
CHARMC	= $(HOME)/charm/bin/charmc $(OPTS)
CC=mpixlcxx $(OPTS)

#FFTWPATH = /soft/apps/fftw-3.1.2-double
ESSLPATH = /soft/apps/ESSL-4.4
WRAPPATH = $(HOME)/esslfftw
#INC = -I$(FFTWPATH)/include
#LIBS = -L$(FFTWPATH)/lib -lfftw3 -lm -lz
INC = -I$(WRAPPATH)/include -I$(ESSLPATH)/include
BGP_LIBS = -L$(ESSLPATH)/lib \
           -L$(WRAPPATH)/lib \
           -L/bgsys/ibm_compilers/sles10/prod/opt/ibmcmp/xlf/bg/11.1/bglib/ \
           -L/soft/apps/ibmcmp/xlsmp/bg/1.7/bglib \
           -L/soft/apps/ibmcmp/xlf/bg/11.1/bglib
LIBS = $(BGP_LIBS) -lfftw3_esslbg -lesslbg -lmass -lxlfmath -lxlf90_r -lxlsmp -lxlomp_ser -lpthread

OBJS = fft1d.o

all: fft1d fft_bench #fft_ref fft1d.prj

fft_bench: fft_bench.o
	${CC} fft_bench.o -o fft_bench $(LIBS)

fft_bench.o: fft_bench.cpp
	${CC} -c fft_bench.cpp $(INC)

fft1d: $(OBJS)
	$(CHARMC) -language charm++ -o fft1d $(OBJS) $(LIBS) -module MeshStreamer

projections: fft1d.prj
fft1d.prj: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections $(LIBS) -lz -o fft1d.prj $(OBJS)

summary: $(OBJS)
	$(CHARMC) -language charm++ -tracemode summary $(LIBS) -o fft1d.sum $(OBJS)

fft1d.decl.h: fft1d.ci
	$(CHARMC)  fft1d.ci

fft_ref: fft_ref.o
	${CC} fft_ref.o -o fft_ref -L/expand/home/arya/fftw-3.3/lib -lfftw3_mpi -lfftw3 -lm

fft_ref.o: fft_ref.cpp
	${CC} -c fft_ref.cpp -I/expand/home/arya/fftw-3.3/include

cleanproj:
	rm -f *.log *.sts *.projrc

clean:
	rm -f *.decl.h *.def.h conv-host *.o fft1d fft1d.prj fft_bench charmrun fft_ref *~

fft1d.o: fft1d.C fft1d.decl.h
	$(CHARMC) -c fft1d.C $(INC)

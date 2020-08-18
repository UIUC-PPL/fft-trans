OPTS = -O3

CHARM_DIR = $(HOME)/work/charm
CHARMC = $(CHARM_DIR)/bin/charmc $(OPTS)
CHARM_INC = -I$(CHARM_DIR)/include

NVCC_FLAGS = -c -std=c++11 -use_fast_math -lineinfo $(OPTS)
CC = mpicxx
MPI_LIBS = -lfftw3 -lm
CHARM_LIBS = -lfftw3 -lm
#CHARM_LIBS = -lcufft

OBJS = fft1d.o fft1dcu.o

all: fft1d fft_ref fft_bench

fft_bench: fft_bench.o
	${CC} fft_bench.o -o fft_bench $(MPI_LIBS)

fft_bench.o: fft_bench.cpp
	${CC} -c fft_bench.cpp $(INC)

fft1d: $(OBJS)
	$(CHARMC) -language charm++ -o fft1d $(OBJS) $(CHARM_LIBS)

projections: fft1d.prj
fft1d.prj: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections $(CHARM_LIBS) -lz -o fft1d.prj $(OBJS)

summary: $(OBJS)
	$(CHARMC) -language charm++ -tracemode summary $(CHARM_LIBS) -o fft1d.sum $(OBJS)

fft1d.ci: fft1d.h

fft1d.decl.h: fft1d.ci
	$(CHARMC) fft1d.ci

fft1d.def.h: fft1d.ci

fft_ref: fft_ref.o
	${CC} fft_ref.o -o fft_ref -L/expand/home/arya/fftw-3.3/lib -lfftw3_mpi -lfftw3 -lm

fft_ref.o: fft_ref.cpp
	${CC} -c fft_ref.cpp -I/expand/home/arya/fftw-3.3/include

cleanproj:
	rm -f *.log *.sts *.projrc

clean:
	rm -f *.decl.h *.def.h conv-host *.o fft1d fft1d.prj fft_bench charmrun fft_ref *~

fft1d.o: fft1d.C fft1d.decl.h fft1d.def.h fft1d.h
	$(CHARMC) -c fft1d.C

fft1dcu.o: fft1d.cu fft1d.h
	nvcc -o $@ $(NVCC_FLAGS) $(CHARM_INC) $<

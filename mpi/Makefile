CHARMHOME?=$(HOME)/charm/mpi-linux-x86_64/
CHARMC=$(CHARMHOME)/bin/charmc $(OPTS)
OPTS = -O3

CXX=mpicc $(OPTS)
FFTW3DIR?=$(HOME)/fftw-3.3

INC = -I$(FFTW3DIR)/include
LIBDIR = -L$(FFTW3DIR)/lib 
LIBS = -lfftw3 -lm -module NDMeshStreamer -module completion

all: fft_charm_mpi

fft_charm_mpi: fft_charm_mpi.cc libmodulefft1d_mpi_wrapper.a
	$(CXX) -c fft_charm_mpi.cc -o fft_charm_mpi.o -I$(CHARMHOME)/include $(INC)
	$(CHARMC) -mpi -o fft_charm_mpi fft_charm_mpi.o -L./ -module fft1d_mpi_wrapper \
	$(LIBDIR) $(LIBS) 

libmodulefft1d_mpi_wrapper.a: fft1d_mpi_wrapper.cc fft1d_mpi_wrapper.decl.h fft1d_mpi_wrapper.def.h  fft.decl.h
	$(CHARMC) $(OPTS) -o libmodulefft1d_mpi_wrapper.a fft1d_mpi_wrapper.cc $(INC)

fft1d_mpi_wrapper.decl.h fft1d_mpi_wrapper.def.h: fft1d_mpi_wrapper.ci
	$(CHARMC) $(OPTS) fft1d_mpi_wrapper.ci

fft.decl.h: fft.ci
	$(CHARMC) fft.ci

test: fft_charm_mpi
	mpirun -np 8 ./fft_charm_mpi 1024 1

clean:
	rm -f *dump* *.decl.h *.def.h conv-host *.o *~

clobber: clean
	rm -f fft_charm_mpi charmrun libmodulefft1d_mpi_wrapper.a *.o charmrun  libmodulefft1d_mpi_wrapper.a

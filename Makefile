OPTS	= -O3
CHARMC	= $(HOME)/charm/bin/charmc $(OPTS)

OBJS = fft1d.o

all: fft1d

fft1d: $(OBJS)
	$(CHARMC) -language charm++ -o fft1d $(OBJS) -lfftw3 -lm

projections: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections -lz -o fft1d.prj $(OBJS)

summary: $(OBJS)
	$(CHARMC) -language charm++ -tracemode summary -lz -o fft1d.sum $(OBJS)

fft1d.decl.h: fft1d.ci
	$(CHARMC)  fft1d.ci

clean:
	rm -f *.decl.h *.def.h conv-host *.o fft1d fft1d.prj charmrun *~

fft1d.o: fft1d.C fft1d.decl.h
	$(CHARMC) -c fft1d.C

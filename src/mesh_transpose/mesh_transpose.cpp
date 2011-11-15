#include <fftw3.h>
#include <MeshStreamer.h>
#include <TopoManager.h>
#include <gluon/gluon.h>

#include <transpose.h>

#define N2 100
#define NCHARE 2
#define BUFSIZE N2*N2/NCHARE/NCHARE

struct fftBuf
{
	int iter;
	int source;
	fftw_complex data[BUFSIZE];
};

#include "mesh_transpose.decl.h"

struct fftMsg:
			public CMessage_fftMsg
{
	int source;
	fftw_complex *data;
};

struct mesh_transpose:
//		public CBase_mesh_transpose,
		public MeshStreamerClient<fftBuf>
{
	mesh_transpose_SDAG_CODE
	
	//Externally configured (see bellow)
	uint32_t numChares;
	void size ( uint32_t size ) 
	{
		numChares = size;
	}
	
	int iteration, count;
	uint64_t n;
	fftw_plan p1;
	fftBuf *msg;
	fftw_complex *in, *out;
	int thisIndex;
	MeshStreamerMessage<fftBuf> *msg1;

	mesh_transpose()
	{
		__sdag_init();

		thisIndex = CkMyPe();

		msg = new fftBuf;
	}
	
	void init ( uint64_t N )
	{
		this->N = N;
		n = N*N/numChares;
	}

	void sendTranspose(fftw_complex *src_buf) {
		// All-to-all transpose by constructing and sending
		// point-to-point messages to each chare in the array.
		for (int i = thisIndex; i < thisIndex+numChares; i++) {
			//  Stagger communication order to avoid hotspots and the
			//  associated contention.
			int k = i % numChares;
			for (int j = 0, l = 0; j < N/numChares; j++)
				memcpy(msg->data[(l++)*N/numChares], src_buf[k*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);

			// Tag each message with the iteration in which it was
			// generated, to prevent mis-matched messages from chares that
			// got all of their input quickly and moved to the next step.
			// Runtime system takes ownership of messages once they're sent
			msg->iter = iteration;
			msg->source = thisIndex;
			((MeshStreamer<fftBuf> *)CkLocalBranch(aggregator))->insertData(*msg, k);
		}

		contribute(CkCallback(CkIndex_Main::startFlush(), mainProxy));
	}

	void receiveCombinedData(MeshStreamerMessage<fftBuf> *msg) {
		msg1 = msg;
		for (int i = 0; i < msg->numDataItems; i++) {
			fftBuf *m1 = &((msg->data)[i]);
			fftMsg *m2 = new fftMsg;
			m2->source = m1->source;
			m2->data = m1->data;
			CkSetRefNum(m2, m1->iter);
			processData(m2);
		}
	}

	void applyTranspose(fftMsg *m)
	{
		int k = m->source;
		for (int j = 0, l = 0; j < N/numChares; j++) {
			for (int i = 0; i < N/numChares; i++) {
				out[k*N/numChares+(i*N+j)][0] = m->data[l][0];
				out[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
			}
		}
		delete m;
	}

	void calcResidual()
	{
		double infNorm = 0.0;

		srand48(thisIndex);
		for (int i = 0; i < n; i++) {
			out[i][0] = out[i][0]/(N*N) - drand48();
			out[i][1] = out[i][1]/(N*N) - drand48();

			double mag = sqrt(pow(out[i][0], 2) + pow(out[i][1], 2));
			if (mag > infNorm) infNorm = mag;
		}

		double r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));

		CkCallback cb(CkReductionTarget(Main, printResidual), mainProxy);
		contribute(sizeof(double), &r, CkReduction::max_double, cb);
	}
};



#include "mesh_transpose.def.h"

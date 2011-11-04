#include <vector>

#include <gluon/gluon.h>

#include <transpose.h>

#include <p2p_transpose.decl.h>

struct fftMsg : public CMessage_fftMsg {
	int source;
	fftw_complex *data;
};


struct p2p_transpose : public CBase_p2p_transpose
{
	//TODO: cfg
	transpose_callback* from_fft;
	uint32_t numChares;
	uint64_t N, n;

	void init(uint64_t N) {
		this->N = N;
		n = N*N/numChares;
		msgs = new fftMsg*[numChares];
		for (int i = 0; i < numChares; i++) {
			msgs[i] = new (n/numChares) fftMsg;
			msgs[i]->source = thisIndex;
		}
	}

	fftMsg **msgs;
	int count;
	fftw_complex* out_buf;

	p2p_transpose() {
		__sdag_init();
	}
	p2p_transpose(CkMigrateMessage *) { }

	p2p_transpose_SDAG_CODE

	void sendTranspose(int iteration, fftw_complex *src_buf, fftw_complex* out_buf, transpose_callback& callback) {
		this->from_fft = &callback;
		this->out_buf = out_buf;

		// All-to-all transpose by constructing and sending
		// point-to-point messages to each chare in the array.
		for (int i = thisIndex; i < thisIndex+numChares; i++) {
			//  Stagger communication order to avoid hotspots and the
			//  associated contention.
			int k = i % numChares;
			for (int j = 0, l = 0; j < N/numChares; j++)
				memcpy(msgs[k]->data[(l++)*N/numChares], src_buf[k*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);

			// Tag each message with the iteration in which it was
			// generated, to prevent mis-matched messages from chares that
			// got all of their input quickly and moved to the next step.
			CkSetRefNum(msgs[k], iteration);
			thisProxy[k].getTranspose(msgs[k]);
			// Runtime system takes ownership of messages once they're sent
			msgs[k] = NULL;
		}

		thisProxy[thisIndex].receive(iteration);
	}

	void applyTranspose(fftMsg *m, fftw_complex *out) {
		int k = m->source;
		for (int j = 0, l = 0; j < N/numChares; j++)
			for (int i = 0; i < N/numChares; i++) {
				out[k*N/numChares+(i*N+j)][0] = m->data[l][0];
				out[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
			}

		// Save just-received messages to reuse for later sends, to
		// avoid reallocation
		delete msgs[k];
		msgs[k] = m;
		msgs[k]->source = thisIndex;
	}
};

#include "p2p_transpose.def.h"

GCMP_A(p2p_transpose)
	G_PROPERTY(uint32_t, numChares);
	G_CHARM_APROVIDE(transpose, to_fft);
GEND

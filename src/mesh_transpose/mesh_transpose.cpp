#include <MeshStreamer.h>
#include <gluon/gluon.h>

#include <transpose_to_fft.h>

#include "mesh_transpose.decl.h"

struct mesh_transpose:
		public CBase_mesh_transpose,
		virtual public transpose_to_fft
{
	fftw_complex *buf;
	int count;
	
	mesh_transpose() {
		__sdag_init();
		buf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n/numChares+1));
	}
	
	mesh_transpose_SDAG_CODE
	
	void sendTranspose(int iteration, fftw_complex *src_buf, fftw_complex* out, fft_to_transpose& callback) {
		// All-to-all transpose by constructing and sending
		// point-to-point messages to each chare in the array.
		for (int i = thisIndex; i < thisIndex+numChares; i++) {
			//  Stagger communication order to avoid hotspots and the
			//  associated contention.
			int k = i % numChares;
			for (int j = 0, l = 0; j < N/numChares; j++)
				memcpy(buf[(l++)*N/numChares+1], src_buf[k*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);
			
			buf[0][0] = iteration;
			buf[0][1] = thisIndex;
			
			aggregator.ckLocalBranch()->insertData((double *)buf, k);
		}
		
		if (thisIndex == 0)
			CkStartQD(CkCallback(CkIndex_MeshStreamer<double>::flushDirect(), aggregator));
	}
	
	void applyTranspose(MeshStreamerMessage<double> *m, fftw_complex *out) {
		int k = ( (fftw_complex*) m->data)[0][1];
		for (int j = 0, l = 1; j < N/numChares; j++)
			for (int i = 0; i < N/numChares; i++) {
				out[k*N/numChares+(i*N+j)][0] = ( (fftw_complex*) m->data)[l][0];
				out[k*N/numChares+(i*N+j)][1] = ( (fftw_complex*) m->data)[l++][1];
			}
			
			delete m;
	}
};

#include "mesh_transpose.def.h"

/* This file contains a  f1qclean function of FFTQuiver plugin for vapoursynth
// Row by row the image is transormed into frequency domain, frequency filtered and 
transformed back into row. In addition to a large number of Butterworth
filters, filter can be custom designed.

  This plugin needs any one of libfftw3f-3.dll,
  FFTW3 dll, fftw.dll to reside in path (may be windows\system32 folder)
  
Author V.C.Mohan. 
jun 2015, 14 sep 2020, 26 May 2021

  Copyright (C) <2014 - 2021>  <V.C.Mohan>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    see <http://www.gnu.org/licenses/>.
*/	
//---------------------------------------------------------------------------
/*
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"

#include  <math.h>
#define _USE_MATH_DEFINES
#include <vector>
#include "windows.h"
#include "fftwlite.h"
#include "FQDomainHelper.h"
*/
//-----------------------------------------------------------------------------


typedef struct
{
		VSNodeRef *node;
		const VSVideoInfo *vi;	
		
		int span;		// number of filters or for custo pairs specified
		int from;
		int upto;
		int option;
		int limit; // original value reduced to this %age  
		int frequency[10];
		int nfrequencies;
		float* ampSquareBuf;	//  freq response of filter 
		
		fftwf_plan 	pf, pin;		// fftwf creates a plan of process. pointer to it

		int		wbest;	// wbest dimension for speed
		
		int freqWidth;

#include "fftLateBindingClassParams.cpp"
		float* inBuf;
		fftwf_complex* outBuf;
		float**sortBuf;
	
} F1QClean;	
	
class LesserThan {
public:

	template <typename finc>
	bool operator()(const finc* f1, const finc* f2)	
	{
		return *f1 < *f2;
	}
};

//std::sort(sortBuf, sortBuf + span, LesserThan());
template <typename finc>
void f1qCleanProcessFull(F1QClean* d, const finc * sp, finc *dp, const int pitch,
			const int wd , const int ht,  finc min, finc max );
void getAmpSqValues(float *ampSquareBuf, fftwf_complex * outBuf, int freqWidth);

void cleanOutBuf(fftwf_complex* outBuf, float* ampSquareBuf, 
				float **sortBuf, int span,int from,int upto, int freqWidth);
void limitMaxAmplitudeInSpan(fftwf_complex* outBuf, int frequency, int span, int limit);

void scaleValues(fftwf_complex* outBuf, int freqWidth, float scale);


// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC f1qcleanInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, 
								VSCore *core, const VSAPI *vsapi) 
{
    F1QClean *d = (F1QClean *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);

	int * facbuf = (int *) vs_aligned_malloc <int>(sizeof( int) *64, 32);	//maximum 64 factors, in this buf values filled are pairs of first is factor, second is dividend to be factored. At 
								// a value of 1 no more factors	
	//	wbest dimensions for speed. make sure starting with even number for width
	int wdEven = ((d->vi->width + 3) >> 2) << 2;
	d->wbest = getBestDim(wdEven, facbuf);
	
	vs_aligned_free(facbuf);

	d->freqWidth = d->wbest / 2 + 1;
	d->span = ((d->span * d->freqWidth) / NYQUIST) | 1;

	if (d->option == 2)
	{
		d->from = (d->from * d->freqWidth) / NYQUIST;
		d->upto = (d->upto * d->freqWidth) / NYQUIST;
	}
	else
	{
		for (int i = 0; i < d->nfrequencies; i ++)

			d->frequency[i] = (d->frequency[i] * d->freqWidth) / NYQUIST;

	}
#include "ConstructorCodeForLateBindingfft.cpp"

	if (!ok)
	{
		vsapi->setError(out, "F1QClean or F1QLimit: could not load any of the dll or get required fnctions");
		if (d->hinstLib != NULL)
			FreeLibrary(d->hinstLib);
		vsapi->freeNode(d->node);
		return;
	}	
	 // create fft plans. Requires buffers 
	d->inBuf =  (float *)d->fftwf_malloc (sizeof(float) * d->wbest );

	d->outBuf = (fftwf_complex*) d->fftwf_malloc (sizeof(fftwf_complex) * d->freqWidth);

	if (d->option == 2) // F2QClean auto sets it to 2
	{
		d->ampSquareBuf = (float*)d->fftwf_malloc(sizeof(float) * d->freqWidth);	

		d->sortBuf = (float**)vs_aligned_malloc(sizeof(float*) * d->span, 32);
	}
			// get fft sine cosine config buffers allocated by plans	
	d->pf = d-> fftwf_plan_dft_r2c_1d( d->wbest, d->inBuf, d->outBuf, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	d->pin = d->fftwf_plan_dft_c2r_1d( d->wbest, d->outBuf, d->inBuf, FFTW_MEASURE | FFTW_DESTROY_INPUT);
	
}

//---------------------------------------------------------------------------------------------------------------------------------


template <typename finc>

void f1qCleanProcessFull(F1QClean* d, const finc * sp, finc * dp, const int pitch,
					 const int wd ,const int ht, finc min, finc max)
{
	float scale = 1.0f / d->wbest;

	for(int h = 0; h < ht ; h++)
	{		
		getRowInput(d->inBuf, sp, d->wbest, wd );

		d->fftwf_execute(d->pf);

		scaleValues(d->outBuf, d->freqWidth, scale);

		if (d->option == 2)		// F2QClean
		{
			//search  spectrum between from and upto
			getAmpSqValues(d->ampSquareBuf, d->outBuf, d->freqWidth);

			cleanOutBuf(d->outBuf, d->ampSquareBuf, d->sortBuf, d->span, d->from, d->upto, d->freqWidth);
		}

		else if (d->option == 1)		// F2QLimit
		{
			for (int i = 0; i < d->nfrequencies; i++)
			{
				limitMaxAmplitudeInSpan(d->outBuf, d->frequency[i], d->span, d->limit);
			}

		}

		d->fftwf_execute(d->pin);

		getRowOutput(d->inBuf, dp, wd, min, max);

		sp += pitch;
		dp += pitch;

	}	
}
//....................................................................
void getAmpSqValues(float *ampSquareBuf, fftwf_complex* outBuf, int freqWidth)
{
	for (int i = 0; i < freqWidth; i++)
	{
		ampSquareBuf[i] = getAmpSquareOfComplex(outBuf + i);
	}
}
void scaleValues( fftwf_complex* outBuf, int freqWidth, float scale)
{
	for (int i = 0; i < freqWidth; i++)
	{
		outBuf[i][0] *= scale;
		outBuf[i][1] *= scale;
	}
}

void cleanOutBuf(fftwf_complex* outBuf, float* ampSquareBuf,
	float** sortBuf, int span,int from,int upto, int freqWidth)
{
	
	int center = span / 2;

	for (int w = from - center; w < upto - center; w++)
	{
		for (int s = 0; s < span; s++)
		{
			sortBuf[s] = ampSquareBuf + w + s;
		}
		std::sort(sortBuf, sortBuf + span, LesserThan());
		int median = (int)(sortBuf[center] - ampSquareBuf);
		if (ampSquareBuf[w + center] > *(sortBuf [center]))
		{
			outBuf[w + center][0] = outBuf[median][0];
			outBuf[w + center][1] = outBuf[median][1];
		}
	}
	
}

void limitMaxAmplitudeInSpan(fftwf_complex * outBuf, int frequency, int span, int limit)
{
	float max = 0.0f;
	float local;
	float lim = limit / 100.0f;
	int point = 0;

	for (int i = frequency - span; i < frequency + span; i++)
	{
		local = getAmpSquareOfComplex(outBuf + i);
		if (local > max)
		{
			max = local;
			point = i;
		}
	}
	// point and neighbours are zeroed
	for (int i = point - 1; i <= point + 1; i++)
	{
		outBuf[i][0] *= lim;
		outBuf[i][1] *= lim;
	}
}


//---------------------------------------------------------------------------------------------------------
// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC f1qcleanGetFrame(int n, int activationReason, void **instanceData, void **frameData,
						VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
{
    F1QClean *d = (F1QClean *) * instanceData;

    if (activationReason == arInitial)
	{
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } 
	else if (activationReason == arAllFramesReady)
	{
		
		const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);

		const VSFormat* fi = d->vi->format;
		// process R,G,B or  Y component only will be processed
		int nplanes = fi->colorFamily == cmRGB ? 3 : fi->numPlanes;
		
		VSFrameRef* dst = vsapi->copyFrame(src, core);

		int nbytes = fi->bytesPerSample;
		int nbits = fi->bitsPerSample;

		for (int p = 0; p < nplanes; p++)
		{
			int ht = vsapi->getFrameHeight(src, p);
			int wd = vsapi->getFrameWidth(src, p);

			const uint8_t* srcp = vsapi->getReadPtr(src, p);
			int src_stride = vsapi->getStride(src, p);
			uint8_t* dstp = vsapi->getWritePtr(dst, p);
			int pitch = src_stride / nbytes;

			if (fi->sampleType == stInteger && nbits == 8)
			{
				uint8_t max = (1 << nbits) - 1, min = 0;

				f1qCleanProcessFull(d, srcp, dstp, pitch,
					wd, ht, min, max);
			}

			else if (fi->sampleType == stInteger && nbits > 8)
			{
				uint16_t* dp = (uint16_t*)dstp;

				const uint16_t* sp = (const uint16_t*)srcp;

				uint16_t max = (1 << nbits) - 1, min = 0;

				f1qCleanProcessFull(d, sp, dp, pitch,
					wd, ht, min, max);
			}

			else // float
			{
				float* dp = (float*)dstp;

				const float* sp = (const float*)srcp;

				float max = 1.0f, min = 0.0f;

				f1qCleanProcessFull(d, sp, dp, pitch,
					wd, ht, min, max);
			}
		}
			// Release the source frame
		vsapi->freeFrame(src);				
		return dst;
    }

    return 0;
}

// Free all allocated data on filter destruction
static void VS_CC f1qcleanFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
    F1QClean *d = (F1QClean *)instanceData;
    vsapi->freeNode(d->node);
	if (d->option == 2)
	{
		d->fftwf_free(d->ampSquareBuf);
		vs_aligned_free(d->sortBuf);
	}
	d->fftwf_free(d->inBuf);
	d->fftwf_free(d->outBuf);
	d->fftwf_destroy_plan(d->pf);
	d->fftwf_destroy_plan(d->pin);
	if (d->hinstLib != NULL)
		FreeLibrary(d->hinstLib);
	
    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC f1qcleanCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
    F1QClean d;
    F1QClean *data;
    int err;
	//int temp;
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    // In this first version we only want to handle 8bit integer formats. Note that
    // vi->format can be 0 if the input clip can change format midstream.
    if (!isConstantFormat(d.vi) || d.vi->width == 0 || d.vi->height == 0 
		|| (d.vi->format->colorFamily != cmYUV 	&& d.vi->format->colorFamily != cmGray
			 &&  d.vi->format->colorFamily != cmRGB) )
	{
        vsapi->setError(out, "F1QClean: only RGB, Yuv or Gray color constant formats and const frame dimensions input supported");
        vsapi->freeNode(d.node);
        return;
    }
	
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "F1QClean: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string, the only
    // reason this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
	d.option = 2;
	d.limit = 0;
	d.span = 3;
	d.from = int64ToIntS(vsapi->propGetInt(in, "span", 0, &err));
	if (err)
		d.span = 5;
	else if (d.span < 3 || d.span > 63 || (d.span & 1 ) == 0 )
	{
		vsapi->setError(out, "F1QClean: span must be be odd number between 3 and 63 ");
		vsapi->freeNode(d.node);
		return;
	}
	
	d.from = int64ToIntS(vsapi->propGetInt(in, "fromf", 0, &err));
	if (err)
		d.from = 30;
	
	else if (d.from < 10 + d.span / 2 || d.from > NYQUIST / 2 - 11 - d.span / 2)
	{
		vsapi->setError(out, "F1QClean: fromf must be be between 10 + half of span and less than  245 - half of span ");
		vsapi->freeNode(d.node);
		return;
	}

	d.upto = int64ToIntS(vsapi->propGetInt(in, "upto", 0, &err));
	if (err)
		d.upto = NYQUIST -  10 - d.span;
	
	else if (d.upto < d.from + d.span  || d.upto > NYQUIST - 10 - d.span)
	{
		vsapi->setError(out, "F1QClean: upto can be between fromf + span to  502 - span  ");
		vsapi->freeNode(d.node);
		return;
	}
	

    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = (F1QClean *) malloc(sizeof(d));
    *data = d;

	vsapi->createFilter(in, out, "F1QClean", f1qcleanInit, f1qcleanGetFrame, f1qcleanFree, fmParallelRequests, 0, data, core);

}

// registerFunc("F1QClean", "clip:clip;span:int:opt;fromf:int:opt;upto:int:opt;", f1qcleanCreate, 0, plugin);

// This function is responsible for validating arguments and creating a new filter
static void VS_CC f1qlimitCreate(const VSMap* in, VSMap* out, void* userData, VSCore* core, const VSAPI* vsapi)
{
	F1QClean d;
	F1QClean* data;
	int err;
	//int temp;
	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);
	// vi->format can be 0 if the input clip can change format midstream.
	if (!isConstantFormat(d.vi) || d.vi->width == 0 || d.vi->height == 0
		|| (d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray
			&& d.vi->format->colorFamily != cmRGB))
	{
		vsapi->setError(out, "F1QLimit: only RGB, Yuv or Gray color constant formats and const frame dimensions input supported");
		vsapi->freeNode(d.node);
		return;
	}

	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "F1QLimit: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}

	// If a property read fails for some reason (index out of bounds/wrong type)
	// then err will have flags set to indicate why and 0 will be returned. This
	// can be very useful to know when having optional arguments. Since we have
	// strict checking because of what we wrote in the argument string, the only
	// reason this could fail is when the value wasn't set by the user.
	// And when it's not set we want it to default to enabled.
	d.option = 1;	
	d.span = int64ToIntS(vsapi->propGetInt(in, "span", 0, &err));
	if (err)
		d.span =  15;
	else if (d.span < 3 || d.span > NYQUIST / 8 || (d.span & 1) == 0)
		{
			vsapi->setError(out, "F1QLimit: span must  be odd number  3 to 63");
			vsapi->freeNode(d.node);
			return;
		}
	d.span |= 1;	// make it odd number

	d.limit = int64ToIntS(vsapi->propGetInt(in, "limit", 0, &err));
	if (err)
		d.limit = 50;
	else if (d.limit < 0 || d.limit > 99)
	{
		vsapi->setError(out, "F1QLimit: limit percentage value can be 0 to 99");
		vsapi->freeNode(d.node);
		return;
	}
	
	d.nfrequencies = vsapi->propNumElements(in, "freqs");
	if (d.nfrequencies == 0 || d.nfrequencies > 10)
	{
		vsapi->setError(out, "F1QLimit: for option 1, at least one and not more than 10 freqs must be specified in the array");
		vsapi->freeNode(d.node);
		return;
	}

	for (int i = 0; i < d.nfrequencies; i++)
	{
		d.frequency[i] = int64ToIntS(vsapi->propGetInt(in, "freqs", i, 0));
		if (d.frequency[i] < 10 + d.span / 2 || d.frequency[i] > NYQUIST - 10 - d.span / 2)
		{
			vsapi->setError(out, "F1QLimit:  freqs must be between 10 + half of span and 502 - half of span");
			vsapi->freeNode(d.node);
			return;
		}

	}
	

	// I usually keep the filter data struct on the stack and don't allocate it
	// until all the input validation is done.
	data = (F1QClean*)malloc(sizeof(d));
	*data = d;

	vsapi->createFilter(in, out, "F1QLimit", f1qcleanInit, f1qcleanGetFrame, f1qcleanFree, fmParallelRequests, 0, data, core);

}

// registerFunc("F1QLimit", "clip:clip;span:int:opt;freqs:int[];", f1qlimitCreate, 0, plugin);
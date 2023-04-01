/******************************************************************************
F2QBokeh filter plugin for vapoursynth( 32 & 64 bit) by V.C.Mohan
This filter operates in freq domain (2d), detects parts in focus.
Parts in focus are inserted into a blurred clip to mimic Bokeh effect.

This plugin needs any one of libfftw3f-3.dll 32bit and 64bit of FFTW.org to reside in path
(may be windows\system32 folder, or wow)

Author V.C.Mohan.
20 Dec 2020. 25 May 2021

Copyright (C) <2020-2021>  <V.C.Mohan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is at
see <http://www.gnu.org/licenses/>.
For details of how to contact author see <http://www.avisynth.nl/users/vcmohan/vcmohan.html>
********************************************************************************/
/******************************************************************************
f2qBokeh filter plugin for vapoursynth by V.C.Mohan
This filter operates in freq domain (2d) and f2qbokehs
linear (motion) or circular (focus) styles within a window

Author V.C.Mohan.
 20 Nov 2020

  Copyright (C) < 2008- 2020>  <V.C.Mohan>

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, version 3 of the License.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	A copy of the GNU General Public License is at
	see <http://www.gnu.org/licenses/>.
	For details of how to contact author see <http://www.avisynth.org/vcmohan>
********************************************************************************/
/*
#include <Windows.h>
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
#include <math.h>
//#include <complex>
#include "fftwlite.h"

#include "factorize.cpp"
*/
#include "FQDomainHelper.h"


//----------------------------------------------------------------------------------------

typedef struct
{
	VSNodeRef* node;
	VSNodeRef* nodeB;
	const VSVideoInfo* vi;
	float thresh;	// threshold above which is blur
	int grid;
	int rgb[3];
	int yuv[3];

	int block;		//  block X block is dim of fft
	bool proc[3];
	int bestR;		// freq width for real data will be  bwd / 2 + 1.and used for filter design 
	int f2size;		// complex freq size
	int* gridLUT;	// lookup table for grid from input frame
	int* circleLUT;	// look up table for copying circle from input to output frame
	int noffsets, count;// number of points in lookup tables

	fftwf_plan  pf, pinv;	// forward and inverse fft
#include "fftLateBindingClassParams.cpp"
	float* inBuf;
	fftwf_complex* outBuf;

}F2QBokehData;

//------------------------------------------------------------------------------
void AutoCorrelate(fftwf_complex* Afreq, int block, bool center);
float NormalizeSpectrum(fftwf_complex* buf, int size, bool center);

template <typename finc>
bool isBlockSharp(F2QBokehData * d, float* buf, fftwf_complex* freq, 
	const finc* ap, int apitch);

template <typename finc>
float makeZeroMean(float* in, const finc* sp,int pitch, int* offsets,int grid, int block);

bool checkBlur(float* in, int block, float thresh);

//-----------------------------------------------------------------------------------------------

bool  checkBlur(float* in, int block, float thresh)
{
	int centerOffset = (block / 2) * block + block / 2;
	float maxGo = *(in + centerOffset) * thresh;
	return (maxGo > in[centerOffset - block - 1] || maxGo > in[centerOffset - 1]
		|| maxGo > in[centerOffset + block - 1] || maxGo > in[centerOffset - block]);
}

//--------------------------------------------------------------------------------------------------------------
template <typename finc>
float  makeZeroMean(float* in, const finc* sp,int pitch, int* offsets, int grid, int block)
{
	float mean = getMeanValue(sp, offsets, grid * grid);
	
	float* inBuf = in;

	for (int h = 0; h < grid; h++)
	{
		for (int w = 0; w < grid; w++)
		{
			inBuf[w] = (float)sp[w] - mean;
		}
		for (int w = grid; w < block; w++)
		{
			// fill zeroes
			inBuf[w] = 0.0f;
		}

		inBuf += block;
		sp +=  pitch;
	}
	inBuf = in + grid * block;

	for (int h = grid; h < block; h++)
	{
		// fill zeroes
		for (int w = 0; w < block; w++)
		{
			inBuf[w] = 0.0f;
		}
	}
	return mean;
}

//-----------------------------------------------------------------------------------------------------------

void  AutoCorrelate(fftwf_complex* Afreq, int fsize, bool center)
{
	
	float mult = 1.0 / (fsize);

	// complex multiply with conjugate and scale down to compensate fft upscaling
	
	float real;

	if (center)
	{
		int i = 1;

		for (int w = 0; w < fsize; w++)
		{			
			real = (Afreq[w][0] * Afreq[w][0] + Afreq[w][1] * Afreq[w][1]);

			Afreq[w][0] = real * mult * i;
			Afreq[w][1] = 0.0f;

			i = -i;	//			
		}
	}
	else // no centering
	{		

		for (int n = 0; n < fsize; n++)
		{

			real = (Afreq[n][0] * Afreq[n][0] + Afreq[n][1] * Afreq[n][1]);

			Afreq[n][0] = real * mult;
			Afreq[n][1] = 0.0;
		}
	}
}
//---------------------------------------------------------------------------
float NormalizeSpectrum(fftwf_complex* buf, int fsize, bool center)
{
	// normalizing auto correlation spectrum
	float maximum = buf[0][0];

	// get  max values
	for (int n = 0; n < fsize; n++)
	{
		if (maximum < buf[n][0])
			maximum = buf[n][0];

	}
	// normalize
	if (maximum > 0.0001f)
	{
		float mult = 1.0 / (maximum);

		if (center)
		{
			int i = 1;
			for (int n = 0; n < fsize; n++)
			{
				buf[n][0] *= mult * i;
				i = -i;
			}
		}
		else
		{
			for (int n = 0; n < fsize; n++)
			{
				buf[n][0] *= mult;
			}
		}

	}
	return maximum;
}

//--------------------------------------------------------------------------------------------
template <typename finc>
bool isBlockSharp(F2QBokehData * d,float* buf, fftwf_complex* Afreq,
	const finc* ap, int apitch)
{

	float mean = makeZeroMean(buf, ap,apitch, d->gridLUT, d->grid, d->block);
	// 2d forward transform into Afreq complex buffer
	//d->fftwf_execute_dft_r2c(d->pf, buf, Afreq);
	d->fftwf_execute(d->pf);
	// auto correlation in freq domain is multiplication. Also scale down as scaleup happens in FFT
	AutoCorrelate(Afreq, d->f2size, false);
	float max = NormalizeSpectrum(Afreq, d->f2size, true);
	// avoid zero or same color in block
	if (max > 0.0001f)
	{
		// inverse fft
		//d->fftwf_execute_dft_c2r(d->pinv, Afreq, buf);
		d->fftwf_execute(d->pinv);
		// check if in focus
		return  checkBlur(buf, d->block, d->thresh);
	}
	return false;
}


/*************************************************/
static void VS_CC f2qbokehInit(VSMap* in, VSMap* out, void** instanceData, VSNode* node, VSCore* core, const VSAPI* vsapi)
{
	F2QBokehData* d = (F2QBokehData*)*instanceData;
	vsapi->setVideoInfo(d->vi, 1, node);
	
	d->block =  ( ( d->grid + 7) >> 3) << 3;
		// for RGB or Y planes
	d->bestR = d->block / 2 + 1;
	d->f2size = d->block * d->bestR;
	int bsize = d->block * d->block;
#include "ConstructorCodeForLateBindingfft.cpp"	
	// buffers 
	d->inBuf = (float*)d->fftwf_malloc(sizeof(float) * bsize);

	d->outBuf = (fftwf_complex*)d->fftwf_malloc(sizeof(fftwf_complex) * d->f2size);

	// We require one forward for padded size   and two inverse transforms( one of best size and other for padded . As our dimensions are good (multiple of 2x3x5 measure is used.
	d->pf = d->fftwf_plan_dft_r2c_2d(d->block, d->block, d->inBuf, d->outBuf, FFTW_MEASURE);
	// inverse so complex to real used
	d->pinv = d->fftwf_plan_dft_c2r_2d(d->block, d->block, d->outBuf, d->inBuf, FFTW_MEASURE);

	//d->fftwf_free(inBuf);
	//d->fftwf_free(outBuf);

	const VSFrameRef* srcA = vsapi->getFrame(0, d->node, NULL, 0);
	const VSFormat* fi = d->vi->format;
	
	int nBytes = fi->bytesPerSample;
	int pitch = vsapi->getStride(srcA, 0) / nBytes;
	int gPoints = d->grid * d->grid;

	d->gridLUT = (int*)vs_aligned_malloc <int>(sizeof(int) * gPoints * 5, 32);
	d->circleLUT = d->gridLUT + gPoints;

	d->noffsets = makeRectGridLUT(d->gridLUT, pitch, d->grid);
	// will have grid as radius. So a larger area than grid 8 grid will be treated as sharp
	int cOffset = (d->grid / 2) * pitch + d->grid / 2;
	d->count = makeCircularLUT(d->circleLUT, pitch, d->grid , cOffset);
	vsapi->freeFrame(srcA);

	if (d->noffsets != gPoints || d->count >=  4 * gPoints)
	{
		vs_aligned_free(d->gridLUT);
		vsapi->setError(out, "bokeh: noffsets or count are in error");
		vsapi->freeNode(d->node);
		vsapi->freeNode(d->nodeB);
		return;

	}
	
}
//-----------------------------------------------------------------------------------------------

// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef* VS_CC f2qbokehGetFrame(int n, int activationReason, void** instanceData,
	void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	F2QBokehData* d = (F2QBokehData*)*instanceData;

	if (activationReason == arInitial)
	{
		// Request the source frame on the first call
		vsapi->requestFrameFilter(n, d->node, frameCtx);
		vsapi->requestFrameFilter(n, d->nodeB, frameCtx);
	}
	else if (activationReason == arAllFramesReady)
	{

		const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);
		// The reason we query this on a per frame basis is because we want our filter
		// to accept clips with varying dimensions. If we reject such content using d->vi
		// would be better.
		const VSFormat* fi = d->vi->format;
		
		// When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
		// supply the "dominant" source frame to copy properties from. Frame props
		// are an essential part of the filter chain and you should NEVER break it.
		const VSFrameRef* blr = vsapi->getFrameFilter(n, d->nodeB, frameCtx);
		VSFrameRef* dst = vsapi->copyFrame(blr, core);

		float* inBuf = d->inBuf; // (float*)d->fftwf_malloc(sizeof(float) * d->block * d->block);

		fftwf_complex* outBuf = d->outBuf; // (fftwf_complex*)d->fftwf_malloc(sizeof(fftwf_complex) * d->f2size);

		//int* offsets = (int*)vs_aligned_malloc<int>(sizeof(int) * d->block * d->block, 32);		

		
		int nbytes = fi->bytesPerSample;
		const uint8_t* srcp[] = { NULL, NULL, NULL, NULL };
		uint8_t* dstp[] = { NULL, NULL, NULL, NULL };
		int stride = vsapi->getStride(src, 0);;
		int pitch = vsapi->getStride(src, 0) / nbytes;
		int np = fi->numPlanes > 3 ? 3 : fi->numPlanes;

		bool proc[] = { true, true, true };

		for (int p = 0; p < 3; p++)
		{
			if (fi->colorFamily == cmRGB)
			{

				if (d->rgb[p] == 0) proc[2 - p] = false;

			}
			else if (fi->colorFamily == cmYUV)
			{
				if (d->yuv[p] == 0) proc[p] = false;
			}
		}

		for (int plane = 0; plane < np; plane++)
		{
			srcp[plane] = vsapi->getReadPtr(src, plane);
			dstp[plane] = vsapi->getWritePtr(dst, plane);
		}
		int ptrInc = (d->grid / 2) * stride;
		
		int pht = vsapi->getFrameHeight(src, 0);
		int pwd = vsapi->getFrameWidth(src, 0);				

		int count = d->count, noffsets = d->noffsets;

		// initial position offset from zero 
		int iOffset = (d->grid ) * stride;

		for (int p = 0; p < np; p++)
		{
			srcp[p] += iOffset;
			dstp[p] += iOffset;
		}

		for (int h = d->grid ; h <= pht - d->grid - d->grid / 2 - 1; h += d->grid / 2)
		{
			for (int w = d->grid ; w <= pwd - d->grid - d->grid / 2  - 1; w += d->grid / 2)
			{
				bool sharp =  false;
				// process image	
				if (fi->sampleType == stInteger)
				{
					if (nbytes == 1)
					{
						for (int p = 0; p < np; p++)
						{
							if( proc[p])

								sharp = isBlockSharp(d, inBuf, outBuf,
									srcp[p] + w, pitch);
							if (sharp)
								break;
						}

						if (sharp)
						{
							// copy circular area
							for (int p = 0; p < np; p++)
							{
								for (int i = 0; i < count; i++)

									*(dstp[p] + w + d->circleLUT[i])
									= *(srcp[p] + w + d->circleLUT[i]);
							}
						}
					}

					else if (nbytes == 2)
					{
						for (int p = 0; p < np; p++)
						{
							if( proc[p])

								sharp = isBlockSharp(d, inBuf, outBuf,
									(uint16_t*)(srcp[p]) + w, pitch);
							if (sharp)
								break;
						}

						if (sharp)
						{
							// copy circular area
							for (int p = 0; p < np; p++)
							{
								for (int i = 0; i < count; i++)

									*(((uint16_t*)(dstp[p])) + w + d->circleLUT[i])
									= *(((uint16_t*)(srcp[p])) + w + d->circleLUT[i]);
							}
						}
					}
				}
				

				else if (nbytes == 4)
				{
					for (int p = 0; p < np; p++)
					{
						if( proc[p])

							sharp = isBlockSharp(d, inBuf, outBuf,
								(float*)(srcp[p]) + w, pitch);
						if (sharp)
							break;
					}

					if (sharp)
					{
						// copy circular area
						for (int p = 0; p < np; p++)
						{
							for (int i = 0; i < count; i++)

								*(((float*)(dstp[p])) + w + d->circleLUT[i])
								= *(((float*)(srcp[p])) + w + d->circleLUT[i]);
						}
					}
				}

			}
			for (int p = 0; p < np; p++)
			{
				srcp[p] += ptrInc;
				dstp[p] += ptrInc;
			}
		}

		// Release the source frames
		vsapi->freeFrame(src);
		vsapi->freeFrame(blr);

		// release fft buffers
		//d->fftwf_free(inBuf);
		//d->fftwf_free(outBuf);
		
		// A reference is consumed when it is returned, so saving the dst reference somewhere
		// and reusing it is not allowed.
		return dst;
	}

	return 0;

}
//----------------------------------------------------------------
// Free all allocated data on filter destruction
static void VS_CC f2qbokehFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	F2QBokehData* d = (F2QBokehData*)instanceData;

	
	d->fftwf_destroy_plan(d->pf);
	d->fftwf_destroy_plan(d->pinv);	
	d->fftwf_free(d->inBuf);
	d->fftwf_free(d->outBuf);

	if (d->hinstLib != NULL)
		FreeLibrary(d->hinstLib);
	vsapi->freeNode(d->node);
	vsapi->freeNode(d->nodeB);

	vs_aligned_free(d->gridLUT);
	//vs_aligned_free(d->circleLUT); part of gridLUT creation
	free(d);
}


/***************************************************************/


// This function is responsible for validating arguments and creating a new filter
static void VS_CC f2qbokehCreate(const VSMap* in, VSMap* out, void* userData, VSCore* core, const VSAPI* vsapi)
{
	F2QBokehData d;
	F2QBokehData* data;
	int err;
	
	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);
	// vi->format can be 0 if the input clip can change format midstream.
	if (!isConstantFormat(d.vi) || d.vi->format->colorFamily == cmCompat
		|| d.vi->format->colorFamily == pfRGBH || d.vi->format->colorFamily == pfYUV444PH
		|| d.vi->format->colorFamily == pfGrayH)
	{
		vsapi->setError(out, "f2qBokeh: clip must have constant dimensions and in YUV or RGB or Grey format. Half float formats not allowed  ");
		vsapi->freeNode(d.node);
		return;
	}
	if ( d.vi->format-> colorFamily == cmYUV && d.vi->format->subSamplingH != 0 
		|| d.vi->format->subSamplingW != 0)
	{
		vsapi->setError(out, "f2qBokeh: YUV format clip must be type YUV444 i.e. no subsampling.");
		vsapi->freeNode(d.node);
		return;
	}

	d.nodeB = vsapi->propGetNode(in, "clipb", 0, 0);
	const VSVideoInfo* bvi = vsapi->getVideoInfo(d.nodeB);
	if (!isSameFormat(d.vi, bvi))
	{
		vsapi->setError(out, "f2qBokeh: blur clip must have identical format with the input clip  ");
		vsapi->freeNode(d.node);
		vsapi->freeNode(d.nodeB);
		return;

	}

	// If a property read fails for some reason (index out of bounds/wrong type)
	// then err will have flags set to indicate why and 0 will be returned. This
	// can be very useful to know when having optional arguments. Since we have
	// strict checking because of what we wrote in the argument string, the only
	// reason this could fail is when the value wasn't set by the user.
	// And when it's not set we want it to default to enabled.

	d.grid = vsapi->propGetInt(in, "grid", 0, &err);
	if (err)
	{
		d.grid = 16;
	}
	else if (d.grid < 3 || d.grid > 64)
	{
		vsapi->setError(out, "f2qBokeh:grid must be between 3 and 64");
		vsapi->freeNode(d.node);
		vsapi->freeNode(d.nodeB);
		return;
	}
	d.thresh = vsapi->propGetFloat(in, "thresh", 0, &err);
	if (err)
	{
		d.thresh = 0.45f;
	}
	else if (d.thresh < 0.0f || d.thresh > 1.0f)
	{
		vsapi->setError(out, "f2qBokeh:value of thresh must be between 0 qnd 1.0");
		vsapi->freeNode(d.node);
		vsapi->freeNode(d.nodeB);
		return;
	}


	if (d.vi->format->colorFamily == cmRGB)
	{
		int count = vsapi->propNumElements(in, "rgb");
		if (count > 3)
		{
			vsapi->setError(out, "f2qBokeh: rgb array cannot have more than 3 entries.");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.nodeB);
			return;
		}
		else if (count < 1)
		{
			d.rgb[0] = 1;
			d.rgb[1] = 1;
			d.rgb[2] = 1;
		}

		for (int p = 0; p < 3; p++)
		{
			d.rgb[p] = vsapi->propGetInt(in, "rgb", p, &err);

			 if (err)
			{
				d.rgb[p] = d.rgb[p - 1];
			}

			else if (d.rgb[p] < 0 || d.rgb[p] > 1)
			{
				vsapi->setError(out, "f2qBokeh: rgb array can have values of 0 or 1 only.");
				vsapi->freeNode(d.node);
				vsapi->freeNode(d.nodeB);
				return;
			}

		}
		if (d.rgb[0] == 0 && d.rgb[1] == 0 && d.rgb[2] == 0)
		{
			vsapi->setError(out, "f2qBokeh: rgb array all values should not be 0");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.nodeB);
			return;
		}
	}
	else if (d.vi->format->colorFamily == cmYUV)
	{
		int count = vsapi->propNumElements(in, "yuv");
		if (count > 3)
		{
			vsapi->setError(out, "f2qBokeh: yuv array cannot have more than 3 entries.");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.nodeB);
			return;
		}
		else if (count < 1)
		{
			d.yuv[0] = 1;
			d.yuv[1] = 1;
			d.yuv[2] = 1;
		}
		for (int p = 0; p < 3; p++)
		{
			d.yuv[p] = vsapi->propGetInt(in, "yuv", p, &err);

			if (err)
			{
				d.yuv[p] = d.yuv[p - 1];
			}

			else if (d.yuv[p] < 0 || d.yuv[p] > 1)
			{
				vsapi->setError(out, "f2qBokeh: yuv array can have values of 0 or 1 only.");
				vsapi->freeNode(d.node);
				vsapi->freeNode(d.nodeB);
				return;
			}

		}
		if (d.yuv[0] == 0 && d.yuv[1] == 0 && d.yuv[2] == 0)
		{
			vsapi->setError(out, "f2qBokeh: yuv array all values should not be 0");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.nodeB);
			return;
		}
	}
	

	// I usually keep the filter data struct on the stack and don't allocate it
	// until all the input validation is done.
	data = (F2QBokehData*)malloc(sizeof(d));
	*data = d;

	// Creates a new filter and returns a reference to it. Always pass on the in and out
	// arguments or unexpected things may happen. The name should be something that's
	// easy to connect to the filter, like its function name.
	// The three function pointers handle initialization, frame processing and filter destruction.
	// The filtermode is very important to get right as it controls how threading of the filter
	// is handled. In general you should only use fmParallel whenever possible. This is if you
	// need to modify no shared data at all when the filter is running.
	// For more complicated filters, fmParallelRequests is usually easier to achieve as it can
	// be prefetched in parallel but the actual processing is serialized.
	// The others can be considered special cases where fmSerial is useful to source filters and
	// fmUnordered is useful when a filter's state may change even when deciding which frames to
	// prefetch (such as a cache filter).
	// If your filter is really fast (such as a filter that only resorts frames) you should set the
	// nfNoCache flag to make the caching work smoother.

	vsapi->createFilter(in, out, "f2qBokeh", f2qbokehInit, f2qbokehGetFrame, f2qbokehFree, fmParallel, 0, data, core);

}
/*
// The following function is the function that actually registers the filter in vapoursynth
// It is called automatically, when the plugin is loaded to see which functions this filter contains.


	registerFunc("f2qBokeh", "clip:clip;clipb:clip;grid:int:opt;thresh:float:opt;rgb:int[]:opt;yuv:int[]:opt;", f2qbokehCreate, 0, plugin);
*/

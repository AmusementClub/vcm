// This file contains a  f2quiver filter which transforms frame into 2D freq domain through fft
// frequency filters included are low cut, high cut, band pass, band stop and noth. 
// 2D domain symmetry of circular, vertical, horizntal or point can be opted,
// The frequency at which these are positioned are to be specified. 
// Depending on the strength of filter Gauss or butterworth filters are designed. 
// upto 16 such filters can be specified simultaneously. All will be cascaded
// In the test mode the Frequency spectrum is displayed in left part of frame. 
// On the right Designed filters are displayed.
// the basics of the filter api.
// This file may make more sense when
// read from the bottom and up.
/*
This plugin needs any one of libfftw3f-3.dll, FFTW3 dll, fftw.dll to reside in path
(may be windows\system32 folder)

Author V.C.Mohan.
Jun 2015, 18 May 2021

Copyright (C) <2006, 2021>  <V.C.Mohan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is at
see <http://www.gnu.org/licenses/>.


********************************************************************************/
/*
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
#include  <math.h>
#include <vector>
#include "windows.h"
#include "fftwlite.h"
#include "Factorize.cpp"
#include "F2QFilters.h"
#include "F2QuiverSpectralDisplay.h"
*/
typedef struct {
		VSNodeRef *node;
		const VSVideoInfo *vi;
		
		bool  	test;		// is this a test?
		float 	gamma;		// for spectrum scaling
		bool 	morph;		// is homomorphic process required
		int		Fspec[60];
		int 	npoints;		// number of filters or for custo pairs specified		
		float	*logLUT;		// log look up table for morph
		bool	ham;		// flag for hamming the filter 
		int		frad;		// radius of filter in image domain

		// late binding of fft dll
#include "fftLateBindingClassParams.cpp"
		float* inBuf;
		fftwf_complex* outBuf;
		float  * FreqFilter;	//  freq response of filter and powerspectrum buffer pointers
		fftwf_plan 	pf, pinv;		// fftwf creates a plan of process. pointer to it	
		int		hbest, wbest, frqwidth;	// best dimension for speed




	
} F2QuiverData;
//------------------------------------------------------------------------------
template <typename finc>
void displaySpectraAndFilters2D(F2QuiverData* d, float* inBuf, fftwf_complex* outBuf, const finc* sp, finc * dp,
								int pitch, int ht, int wd, finc Grey, finc max);
template <typename finc>
void getFilteredOutput2D(F2QuiverData* d, float* inBuf, fftwf_complex* outBuf,
	const finc* sp, finc* dp, const int pitch, const int ht, const int wd, finc min, finc max);

void designFilter2D(int* FSpec, int npts, float* filter, int ht, int wd);

//====================================================================================

// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC f2quiverInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    F2QuiverData *d = (F2QuiverData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);
	

	int * facbuf =  vs_aligned_malloc<int> (sizeof(int) * 64, 32);	//maximum 64 factors, first is factor, second is dividend to be factored. At 
								// a value of 1 no more factors
	
			//	best dimensions for speed
	// make sure we have even numbers as starting values of width and height
	int wdEven = ((d->vi->width + 3) >> 2) << 2;
	int htEven = ((d->vi->height + 3) >> 2) << 2;
			// As we are not filtering we use nearest to frame frame
	d->wbest = getBestDim(wdEven + ADDSAFE, facbuf);
	d->hbest = getBestDim(htEven + ADDSAFE, facbuf);

	vs_aligned_free(facbuf);

	if( d->test)
	{
		d->wbest = getBestDim(d->vi->width , facbuf);
		d->hbest = getBestDim(d->vi->height ,facbuf);
	}
	d->frqwidth = (d->wbest / 2) + 1;
	int f2qsize = d->hbest * d->frqwidth;

#include "ConstructorCodeForLateBindingfft.cpp"


	if (!ok)
	{
		vsapi->setError(out, "vcm.f2quiver: could not load any of the fft dll or get required fnctions");
		if (d->hinstLib != NULL)
			FreeLibrary(d->hinstLib);
		vsapi->freeNode(d->node);
		return;
	}
	// buffers 
	d->inBuf = (float*)d->fftwf_malloc(sizeof(float) * d->wbest * d->hbest);
	
	d->outBuf = (fftwf_complex*)d->fftwf_malloc (sizeof(fftwf_complex) * f2qsize);//  is only a safeguard not really reqd
	

	if(d->inBuf == NULL || d->outBuf == NULL )
	{
		vsapi->setError(out, "F2Quiver: unexpectedly buffers not allocated error");
		vsapi->freeNode(d->node);
		if (d->hinstLib != NULL)
			FreeLibrary(d->hinstLib);
		free(d);
		return;
	}
	

	d->pf = d->fftwf_plan_dft_r2c_2d(d->hbest, d->wbest, d->inBuf, d->outBuf, FFTW_MEASURE );
	d->pinv = d->fftwf_plan_dft_c2r_2d(d->hbest, d->wbest, d->outBuf, d->inBuf, FFTW_MEASURE );

	if(  d->pf == NULL || d->pinv == NULL)
	{
		vsapi->setError(out, "F2Quiver: unexpected  fft plans  error");
		vsapi->freeNode(d->node);
		d->fftwf_free(d->inBuf);
		d->fftwf_free(d->outBuf);
		if (d->hinstLib != NULL)
			FreeLibrary(d->hinstLib);
		free(d);
		return;
	}
	int nbits = d->vi->format->bitsPerSample;

	d->logLUT = NULL;
	if(d->morph && nbits <= 12 )
	{
		// only up to 12 bits input this buffer is used
		int nval = 1 << nbits;

		d->logLUT = vs_aligned_malloc <float>( sizeof(float) * nval, 32);

		for (int i = 0; i < nval; i ++)
		{
			d->logLUT[i] = log((float)i + 1.0f);	// 1 added to prevent log 0. In output reduce by 2.0
		}
	}
	
	d->FreqFilter = vs_aligned_malloc <float>(sizeof(float) * f2qsize, 32);

		// initialize 
	for (int i = 0; i < f2qsize; i ++)
	{
		// initialize with value of 1. All cascading filters multiply the value
		d->FreqFilter[i] = 1.0f;
	}
	
	designFilter2D(d->Fspec, d->npoints, d->FreqFilter, d-> hbest, d->frqwidth );

		// normalize filter
	float * filt = d->FreqFilter;
	float fmax = filt[0];

	for(int i = 0; i < f2qsize; i ++)
	{			
		if (fmax < filt[i])
			fmax = filt[i];	
	}

	float fscale = d->test ? 1.0f / fmax : 1.0f / (fmax * d->wbest * d->hbest); 
	
		// normalize for display
		
	for ( int i = 0; i < f2qsize; i ++)
	{
		
		filt[i] *= fscale;
	}
	if (d->ham)
	{
#include "hammingCodeInsert.cpp"
	}
	/*if (! d->test)
	{
		// the filter was designed with origin at center of frame. We need to reposition it
		// transform to space domain , remove centering 
		// we can limit length of filter, at cost of effectiveness. 
		// Ensure no wrap around convolution takes place
		// transfer filter to complex buffer. Our filter is zero phase, so imaginary part = 0;

		for ( int i = 0; i < f2qsize; i ++)
		{
			d->outBuf[i][0] = d->FreqFilter[i];
			d->outBuf[i][1] = 0.0;
		}

		// now inverse fft
		
		d->fftwf_execute(d->pinv);

		// remove centering sign
		float * imageSpace = d->inBuf;
		int start = 1;

		for ( int h = 0; h < d->hbest; h ++)
		{
			int wstart = start;
			for (int w = 0; w < d->wbest; w += 2)
			{
					
				imageSpace[w] *= wstart;  // getSign(h, w);
				imageSpace[w + 1] *= wstart;
				wstart = -wstart;
			}
			imageSpace += d->wbest;
			start = -start;
		}

		// zero out central part leaving a margin of frad;
		// = inBuf + d->filterRadius * d->wbest;

		for ( int h = d->filterRadius; h < d->hbest - d->filterRadius; h ++)
		{
			for ( int w = d->filterRadius; w < d->wbest - d->filterRadius; w ++)
			{
				imageSpace[w] = 0.0;
			}
			imageSpace += d->wbest;
		}

		// now transform (back) to frequency domain		
		d->fftwf_execute(d->pf);

		fmax = 0.0;
		for (int i = 0; i < f2qsize; i++)
		{
			if (fmax < d->outBuf[i][0])
				fmax = d->outBuf[i][0];
				
		}

		// transfer real part to freq buffer. We also remove scaling of hbest * wbest due to inverse and forward fft
		fscale = 1.0f / (fmax * d->hbest * d->wbest);

		for( int i = 0; i < f2qsize; i ++)
		{
			d->FreqFilter[i] = d->outBuf[i][0] * fscale;
		}
			
	}*/	// not test
}

//--------------------------------------------------------------------------
// butterworth 2D filters
void designFilter2D( int * spec, int npts,float * freqFilt, int ht, int wd)
{
	
	for ( int i = 0; i < npts; i += 5)
	{	
		int symmetry = spec[i]; // 1 circular, 2 horizontal, 3 vertical 4  Fan  through axis 5 point
		//	int type = spec[i + 1];	// 1 locut, 2 hicut, 3 band pass 4 band stop 5 Notch
		switch (symmetry)
		{
			case 1 : ButterworthCircular2D(freqFilt, spec + i, ht, wd);
					break;
			case 2: ButterworthHorizontal2D(freqFilt, spec + i, ht, wd);
					break;
			case 3: ButterworthVertical2D(freqFilt, spec + i, ht, wd);
					break;
			case 4: FanFilter2D(freqFilt, spec + i, ht, wd);
					break;
			case 5: pointNotchFilter2D(freqFilt, spec + i, ht, wd);
		}					
	}
}
//-------------------------------------------------------------------------------------------------
template <typename finc>
void displaySpectraAndFilters2D(F2QuiverData * d, float * inBuf, fftwf_complex * outBuf, const finc* sp, finc* dp,
	int pitch, int ht, int wd, finc Grey, finc max)
{
	displayFilter2D(dp + wd / 2, pitch, d->FreqFilter, d->frqwidth, wd / 2, ht, max, d->wbest * d->hbest);

	if (d->morph)

		getHMRealInput2D(inBuf, sp, pitch, ht, wd, d->hbest, d->wbest, true, d->logLUT);
	else
		getRealInput2D(inBuf, sp, pitch, ht, wd, d->hbest, d->wbest, true);

	d->fftwf_execute(d->pf);

	// gets a normalized powerspectrum
	float* powerSpectrum = inBuf;	// inBuf is available
	getPowerSpectrum2D(powerSpectrum, outBuf, d->hbest * d->frqwidth);

	displayPowerSpectrum2D(dp, pitch, powerSpectrum, wd / 2, ht, d->frqwidth, max, d->gamma);

	// draw scales on powerspectrum
	drawHorizontalRuler2D(dp, pitch, wd, d->frqwidth, max, NYQUIST);
	drawVerticalRuler2D(dp, pitch, ht, d->hbest, max, NYQUIST);

	// draw scales on filter
	drawHorizontalRuler2D(dp + wd / 2, pitch, wd, d->frqwidth, Grey, NYQUIST);
	drawVerticalRuler2D(dp + wd / 2, pitch, ht, d->hbest, Grey, NYQUIST);
}
//-----------------------------------------------------------------------------------------------
template <typename finc>
void getFilteredOutput2D(F2QuiverData* d, float* inBuf, fftwf_complex* outBuf,
	const finc* sp, finc* dp, const int pitch, const int ht, const int wd, finc min, finc max)
{
	if (d->morph)
		getHMRealInput2D(inBuf, sp, pitch, ht, wd, d->hbest, d->wbest, false, d->logLUT);
	else
		getRealInput2D(inBuf, sp, pitch, ht, wd, d->hbest, d->wbest, false);

	d->fftwf_execute(d->pf);


	ApplyFilter2D(outBuf, d->FreqFilter, d->hbest, d->frqwidth);

	d->fftwf_execute(d->pinv);

	if (d->morph)
		getHMRealOutput2D(inBuf, dp, pitch, ht, wd, d->hbest, d->wbest, min, max);
	else
		getRealOutput2D(inBuf, dp, pitch, ht, wd, d->hbest, d->wbest, min, max);

}
//--------------------------------------------------------------------------------------------------------------------------
static const VSFrameRef* VS_CC f2qtestGetFrame(int n, int activationReason, void** instanceData,
	void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	F2QuiverData* d = (F2QuiverData*)*instanceData;

	if (activationReason == arInitial) {
		// Request the source frame on the first call
		vsapi->requestFrameFilter(n, d->node, frameCtx);
	}
	else if (activationReason == arAllFramesReady)
	{
		const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);
		// The reason we query this on a per frame basis is because we want our filter
		// to accept clips with varying dimensions. If we reject such content using d->vi
		// would be better.
		const VSFormat* fi = d->vi->format;
		int ht = vsapi->getFrameHeight(src, 0);
		int wd = vsapi->getFrameWidth(src, 0);
		VSFrameRef* dst = vsapi->copyFrame(src, core);
		// float input
		int iwidth = d->wbest;
		int owidth = 2 + d->wbest / 2;

		//float* inBuf = (float*)d->fftwf_malloc(sizeof(float) * iwidth * d->hbest);
		//fftwf_complex* outBuf = (fftwf_complex*)d->fftwf_malloc(sizeof(fftwf_complex) * owidth * d->hbest);
		// process Green of RGB or Y of YUV or Grey
		int plane = fi->colorFamily == cmRGB ? 1 : 0;

		const uint8_t* srcp = vsapi->getReadPtr(src, plane);
		int src_stride = vsapi->getStride(src, plane);
		uint8_t* dstp = vsapi->getWritePtr(dst, plane);
		int pitch = src_stride / fi->bytesPerSample;
		int nbits = fi->sampleType == stInteger ? fi->bitsPerSample : 0;
		int nbytes = fi->bytesPerSample;
		int npl = fi->numPlanes > 3 ? 3 : fi->numPlanes;

		float* powerSpectrum = d->inBuf; // in will be available  // (float *) malloc (sizeof(float) * d->hbest * (d->wbest / 2 + 1) );

		if (fi->sampleType == stInteger && nbits == 8)
		{
			uint8_t max = (1 << nbits) - 1;
			uint8_t Grey = 1 << (nbits - 1);
			uint8_t min = 0;
			// display filter on right
			
			displaySpectraAndFilters2D(d, d->inBuf, d->outBuf, srcp, dstp,
				pitch, ht, wd, Grey, max);
			
			if (fi->colorFamily != cmRGB)
			{
				for (int plane = 1; plane < npl; plane++)
				{

					uint8_t* dstp = vsapi->getWritePtr(dst, plane);
					int dst_stride = vsapi->getStride(dst, plane);
					int ht = vsapi->getFrameHeight(dst, plane);
					int wd = vsapi->getFrameWidth(dst, plane);
					int pitch = dst_stride / nbytes;

					// black background 

					fillPlaneWithVal(dstp, pitch, wd, ht, Grey);
				}
			}
		}

		else if (fi->sampleType == stInteger && nbits > 8)
		{
			const uint16_t* sp = (const uint16_t*)srcp;
			uint16_t* dp = (uint16_t*)dstp;
			uint16_t max = (1 << nbits) - 1;
			uint16_t Grey = 1 << (nbits - 1);
			uint16_t min = 0;
			// display filter on right

			displaySpectraAndFilters2D(d, d->inBuf, d->outBuf, sp, dp,
				pitch, ht, wd, Grey, max);
			if (fi->colorFamily != cmRGB)
			{
				for (int plane = 1; plane < npl; plane++)
				{

					uint16_t* dp = (uint16_t*)vsapi->getWritePtr(dst, plane);
					int dst_stride = vsapi->getStride(dst, plane);
					int ht = vsapi->getFrameHeight(dst, plane);
					int wd = vsapi->getFrameWidth(dst, plane);
					int pitch = dst_stride / nbytes;

					// black background 

					fillPlaneWithVal(dp, pitch, wd, ht, Grey);
				}
			}

		}

		else if (fi->sampleType == stFloat)
		{
			const float* sp = (const float*)srcp;
			float* dp = (float*)dstp;
			float max = 1.0;
			float Grey = 0.0;
			float min = 0;
			// display filter on right

			displaySpectraAndFilters2D(d, d->inBuf, d->outBuf, sp, dp,
				pitch, ht, wd, Grey, max);
			if (fi->colorFamily != cmRGB)
			{
				for (int plane = 1; plane < npl; plane++)
				{

					float* dp = (float*)vsapi->getWritePtr(dst, plane);
					int dst_stride = vsapi->getStride(dst, plane);
					int ht = vsapi->getFrameHeight(dst, plane);
					int wd = vsapi->getFrameWidth(dst, plane);
					int pitch = dst_stride / nbytes;

					// black background 

					fillPlaneWithVal(dp, pitch, wd, ht, Grey);
				}
			}
		}


		if (fi->colorFamily == cmRGB)
		{
			// copy Green on to Blu and Red planes
			vs_bitblt(vsapi->getWritePtr(dst, 0), vsapi->getStride(dst, 0),
				vsapi->getWritePtr(dst, 1), vsapi->getStride(dst, 1),
				wd * nbytes, ht);
			vs_bitblt(vsapi->getWritePtr(dst, 2), vsapi->getStride(dst, 2),
				vsapi->getWritePtr(dst, 1), vsapi->getStride(dst, 1),
				wd * nbytes, ht);
		}
		
		vsapi->freeFrame(src);
		return (dst);
	}
	return 0;
}



//---------------------------------------------------------------------------------------------------------------------------
static const VSFrameRef *VS_CC f2quiverGetFrame(int n, int activationReason, void **instanceData, 
						void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    F2QuiverData *d = (F2QuiverData *) * instanceData;

    if (activationReason == arInitial) {
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } 
	else if (activationReason == arAllFramesReady) 
	{
        const VSFrameRef *src = vsapi->getFrameFilter(n, d->node, frameCtx);
        // The reason we query this on a per frame basis is because we want our filter
        // to accept clips with varying dimensions. If we reject such content using d->vi
        // would be better.
        const VSFormat *fi = d->vi->format;
        int ht = vsapi->getFrameHeight(src, 0);
        int wd = vsapi->getFrameWidth(src, 0);
        VSFrameRef *dst = vsapi->copyFrame(src, core);
		// float input
		int iwidth = d->wbest;
		int owidth = 2 + d->wbest / 2;		

		//inBuf = (float*)d->fftwf_malloc (sizeof(float) * iwidth * d->hbest);
		//outBuf = (fftwf_complex*)d->fftwf_malloc(sizeof(fftwf_complex)* owidth * d->hbest);
		// process all of RGB and Y of YUV or Gray
        int np = fi->colorFamily == cmRGB ? 3 : 1;
		int nbits = fi->sampleType == stInteger ? fi->bitsPerSample : 0;
		int nbytes = fi->bytesPerSample;

		for (int p = 0; p < np; p++)
		{
			const uint8_t* srcp = vsapi->getReadPtr(src, p);
			int src_stride = vsapi->getStride(src, p);
			uint8_t* dstp = vsapi->getWritePtr(dst, p);
			int pitch = src_stride / fi->bytesPerSample;

			if (fi->sampleType == stInteger && nbits == 8)
			{
				uint8_t max = (1 << nbits) - 1, min = 0;

				getFilteredOutput2D(d, d->inBuf, d->outBuf, srcp, dstp, pitch, ht, wd, min, max);

			}

			else if (fi->sampleType == stInteger && nbytes == 2)
			{
				uint16_t max = (1 << nbits) - 1, min = 0;
				uint16_t* dp = (uint16_t*)dstp;
				const uint16_t* sp = (const uint16_t*)srcp;

				getFilteredOutput2D(d, d->inBuf, d->outBuf, sp, dp, pitch, ht, wd, min, max);
			}

			else if (fi->sampleType == stFloat)
			{
				float max = 1.0, min = 0.0f;
				float* dp = (float*)dstp;
				const float* sp = (const float*)srcp;

				getFilteredOutput2D(d, d->inBuf, d->outBuf, sp, dp, pitch, ht, wd, min, max);
			}
		}
        

        // Release the source frame
        vsapi->freeFrame(src);

        // A reference is consumed when it is returned, so saving the dst reference somewhere
        // and reusing it is not allowed.
        return dst;
    }

    return 0;
}
//-----------------------------------------------------------------------------------------------------
// Free all allocated data on filter destruction
static void VS_CC f2quiverFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    F2QuiverData *d = (F2QuiverData *)instanceData;
    vsapi->freeNode(d->node);
	if (d->FreqFilter != NULL)
		vs_aligned_free(d->FreqFilter);
	if( d->logLUT != NULL)
		vs_aligned_free(d->logLUT);
	d->fftwf_destroy_plan (d->pf);
	d->fftwf_destroy_plan ( d->pinv);
	d->fftwf_free(d->inBuf);
	d->fftwf_free(d->outBuf);

	if (d->hinstLib != NULL)
		FreeLibrary(d->hinstLib);

    free(d);
}

//--------------------------------------------------------------------------------------------------------------
// This function is responsible for validating arguments and creating a new filter
static void VS_CC f2quiverCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
    F2QuiverData d;
    F2QuiverData *data;
    int err;
	int temp;
	
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

   
	if (!isConstantFormat(d.vi) && (d.vi->format->colorFamily != cmRGB 
						&& d.vi->format->colorFamily != cmYUV 	
						&& d.vi->format->colorFamily != cmGray))
	{
        vsapi->setError(out, "F2Quiver: only constant format RGB YUV or Gray  input supported");
        vsapi->freeNode(d.node);
        return;
    }
	
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "F2Quiver: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	temp = !!vsapi->propGetInt(in, "ham", 0, &err);
	if (err)
		d.ham = false;
	else
		d.ham = temp == 0 ? false : true;
	if (d.ham)
	{
		int maxfrad = (d.vi->height > d.vi->width ? d.vi->height : d.vi->width);
		temp = vsapi->propGetInt(in, "frad", 0, &err);
		if (err)
			temp = 32;
		else
		{
			if (temp < 10 || temp > 50)
			{
				vsapi->setError(out, "F2Quiver: frad being %age of min dimension be 10 to 50 only");
				vsapi->freeNode(d.node);
				return;
			}
		}
		d.frad = (maxfrad * temp) / 100;
	}

	temp = !!vsapi->propGetInt(in, "test", 0, &err);
	if(err)
		d.test = false;
	else
		d.test = temp == 0 ? false: true;
	

	d.npoints = vsapi->propNumElements(in, "fspec");

	if(d.npoints < 5 || (d.npoints % 5 ) != 0 || d.npoints > 60)
	{
		
		vsapi->setError(out, "F2Quiver: fspec at least one and upto 12 filter specifications be given. Each filter is specified as a set of 5 integer values.");
		vsapi->freeNode(d.node);
		return;
	}

	for ( int i = 0; i < d.npoints; i ++)
		d.Fspec[i] = vsapi->propGetInt(in, "fspec", i, 0);
	
	for ( int i = 0; i < d.npoints; i += 5)
	{
		int nfilt = i / 5 + 1;
		int sym = d.Fspec[i];
		int type = d.Fspec[i + 1];
		int frq1 = d.Fspec[i + 2];
		int frq2 = d.Fspec[i + 3];
		int fan1 = d.Fspec[i + 2];
		int fan2 = d.Fspec[i + 3];
		int deg = d.Fspec[i + 4];

		if (sym < 1 || sym > 5)
		{
		
			
			vsapi->setError(out, "F2Quiver:the first number specifying symmetry of each filter must  be from 1 to 5");
			vsapi->freeNode(d.node);
			return;
		}
		if (type < 1 || type > 5)
		{			
			vsapi->setError(out, "F2Quiver:the second number specifying type of each filter must be from 1 to 5");
			vsapi->freeNode(d.node);
			return;
		}

		if (sym < 4 && (type < 3 || type == 5))	// locut, high cut or notch
		{
			if (frq1 < 0 || frq1 >= NYQUIST / 2)
			{
				
				vsapi->setError(out, "F2Quiver:the third number  specifying freq1 of each filter must be between 0 and %d");
				vsapi->freeNode(d.node);
				return;
			}
		}
		else if (sym < 4 && type < 5)	// band pass, band reject
		{
			if (frq1 < 0 || frq2 > NYQUIST / 2 || frq2 <= frq1)
			{
				
				vsapi->setError(out, "F2Quiver:for  symetry 4 and type 5 filter 3rd value freq1 must be positive and less than the 4th value freq2.");
				vsapi->freeNode(d.node);
				return;
			}
		}
		else if (sym == 4)	// fan filter
		{
			if (type == 5 || type < 3)	// only one angle used
			{
				if (fan1 < 1 || fan1 > 179 || (fan1 > 89 && fan1 < 91))
				{
					vsapi->setError(out, "F2Quiver:for the fan filter the third number specifying angle1 must be between 1 and 89 or 91 to 179.");
					vsapi->freeNode(d.node);
					return;
				}
			}
			else
			{
				// band pass or reject
				if (fan1 < 1 || fan1 > 178 || (fan1 > 88 && fan1 < 91))
				{
					vsapi->setError(out, "F2Quiver:for the fan filter the angle1 must be between 1 and 88 or 91 and 178");
					vsapi->freeNode(d.node);	
					return;
				}
				if (fan2 <= fan1 || (fan1 < 89 && fan2 > 89) || fan2 > 179)
				{
					vsapi->setError(out, "F2Quiver:for the fan filter the third number specifying angle1 must be between 1 and 88 or 91 to 178.");
					vsapi->freeNode(d.node);
					return;
				}
			}
		}

		else if (sym == 5)
		{
			if (frq1 < -NYQUIST / 2 || frq1 > NYQUIST / 2 || frq2 < -NYQUIST / 2 || frq2 > NYQUIST / 2)
			{
				vsapi->setError(out, "F2Quiver:for the filterof symetry 5, the xfreq  and yfreq must be within + or - 250");
				vsapi->freeNode(d.node);
				return;
			}
		}

		if (deg < 1 || deg > 24)
		{			
			vsapi->setError(out, "F2Quiver: the fifth value of each filter  must be from 1 to 24 only.");
			vsapi->freeNode(d.node);
			return;
		}
	}

	temp = vsapi->propGetInt(in, "morph", 0, &err);

	if(err)
	{ 
		d.morph = false;
	}
	else
	{
		if ( temp < 0 || temp > 1)
		{
			vsapi->setError(out, "F2Quiver: hm must be  either 0 or 1");
			vsapi->freeNode(d.node);
			return;
		}

		d.morph = temp == 0 ? false: true;
	}

	if( d.test)
	{
		d.gamma = vsapi->propGetFloat(in, "gamma", 0, &err);
		if(err)
			d.gamma = 0.05;
		else
		{
			if ( d.gamma < 1e-5 || d.gamma > 5.0)
			{
				vsapi->setError(out, "F2Quiver: gamma must be non zero +ve and less than 5.0");
				vsapi->freeNode(d.node);
				return;
			}

		}
	}

	
	
	
    data = (F2QuiverData *) malloc(sizeof(d));
    *data = d;
	if (d.test)
		vsapi->createFilter(in, out, "F2Quiver", f2quiverInit, f2qtestGetFrame, f2quiverFree, fmParallelRequests, 0, data, core);
	else
		vsapi->createFilter(in, out, "F2Quiver", f2quiverInit, f2quiverGetFrame, f2quiverFree, fmParallelRequests, 0, data, core);
}



// This file contains a  F2QLimit filter which transforms frame into 2D freq domain through fft
// and limits values after searching for local maximum
/*
This plugin needs any one of libfftw3f-3.dll, FFTW3 dll, fftw.dll to reside in path
(may be windows\system32 folder)

Author V.C.Mohan.
Jun 2015, 22 May 2021

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
		
		int	fspec[60];
		int npoints;		// number of filters or for custo pairs specified
		int	grid;		// square from -grid to grid horizontal and vertical search 
		int	inner;		// - inner to inner limited
		int warn;		// warning level 0 no warning, 1 warns if close to origin, 2 warns if close to axis
		int xgrid, ygrid;
		int xinner, yinner;
						// late binding of fft dll
#include "fftLateBindingClassParams.cpp"
		float* inBuf;
		fftwf_complex* outBuf;	
		fftwf_plan 	pf, pinv;		// fftwf creates a plan of process. pointer to it	
		int	 hbest, wbest;	// best dimension for speed
		int frqwidth;		// width in freq domain;
} F2QLimitData;
//------------------------------------------------------------------------------
template <typename finc>
void F2QLimitProcess(F2QLimitData *d,finc* wp, const int wpitch,
	const int ht, const int wd, 
	const finc min, const finc max);
int getOffsetOfMaxInSearchArea(fftwf_complex* center, int pitch,
								int x, int y, int xgrid, int ygrid);
void limitValuesInner(fftwf_complex* center,  int pitch, int xinner,int yinner, int limit);

void applyLimits(F2QLimitData* d, int i);

//====================================================================================

// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC f2qlimitInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    F2QLimitData *d = (F2QLimitData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);	

	int * facbuf =  vs_aligned_malloc<int> (sizeof(int) * 64, 32);	//maximum 64 factors, first is factor, second is dividend to be factored. At 
								
				
	// make sure we have even numbers as starting values of width and height
	int wdEven = ((d->vi->width + 3) >> 2) << 2;
	int htEven = ((d->vi->height + 3) >> 2) << 2;
	//	best dimensions for speed// a value of 1 no more factors		
	d->wbest = getBestDim(wdEven + ADDSAFE, facbuf);
	d->hbest = getBestDim(htEven + ADDSAFE, facbuf);

	vs_aligned_free(facbuf);
	
	d->frqwidth = (d->wbest / 2) + 1;
	int f2qsize = d->hbest * d->frqwidth;
	d->xgrid = (d->grid * d->frqwidth) / NYQUIST;
	d->ygrid = (d->grid * d->hbest) / NYQUIST;
	d->xinner = (d->inner * d->frqwidth) / NYQUIST;
	d->yinner = (d->inner * d->hbest) / NYQUIST;
	// set up points to limit
	for (int i = 0; i < d->npoints; i += 3)
	{
		d->fspec[i] = (d->fspec[i] * d->frqwidth) / NYQUIST;	// freq along width
		d->fspec[i + 1] = (d->fspec[i + 1] * d->hbest) / NYQUIST; // freq along height
		
		if (d->warn > 0)
		{
			bool xtrue = (d->fspec[i] + d->xgrid >= 0 && d->fspec[i] - d->xgrid <= 0);
			bool ytrue = (d->fspec[i + 1] + d->ygrid >= 0 && d->fspec[i + 1] - d->ygrid <= 0);
			if (d->warn == 1)
			{
				if (xtrue && ytrue)
				{
					vsapi->setError(out, "F2QLimit: search area includes origin");
					vsapi->freeNode(d->node);					
					free(d);
					return;
				}
			}
			else if (xtrue || ytrue)
			{
				vsapi->setError(out, "F2QLimit: search area includes axis");
				vsapi->freeNode(d->node);
				free(d);
				return;
			}

		}
	}

#include "ConstructorCodeForLateBindingfft.cpp"

	if (!ok)
	{
		vsapi->setError(out, "F2QLimit: could not load any of the fft dll or get required fnctions");
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
		vsapi->setError(out, "F2QLimit: unexpectedly buffers not allocated error");
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
		vsapi->setError(out, "F2QLimit: unexpected  fft plans  error");
		vsapi->freeNode(d->node);
		d->fftwf_free(d->inBuf);
		d->fftwf_free(d->outBuf);
		if (d->hinstLib != NULL)
			FreeLibrary(d->hinstLib);
		free(d);
		return;
	}
	

	
}
//-----------...........................................

int getOffsetOfMaxInSearchArea(fftwf_complex* center, int pitch, 
								int x, int y, int xgrid, int ygrid)
{
	int searchOffset = y * pitch + x;
	//int hmax, wmax;
	int offsetMax;
	float maximum = 0.0f, local = 0.0f;
	for (int h = - ygrid; h <=  ygrid; h++)
	{
		for (int w = - xgrid; w <=  xgrid; w++)
		{
			local = getAmpSquareOfComplex(center + searchOffset + h * pitch + w);
			if (local > maximum)
			{
				maximum = local;
				
				offsetMax = searchOffset + h * pitch + w;
			}
		}
	}
	return offsetMax;
}

void limitValuesInner(fftwf_complex* maxC, int pitch, int xinner, int yinner, int limit)
{
	for (int h = -yinner; h <= yinner; h++)
	{
		for (int w = -xinner; w <= xinner; w++)
		{
			maxC[h * pitch + w][0] *= limit / 100.0f;
			maxC[h * pitch + w][1] *= limit / 100.0f;			
		}
	}
}
template <typename finc>
void F2QLimitProcess(F2QLimitData * d, finc* wp, const int wpitch,
	const int ht, const int wd, const finc min, const finc max)	
{
	// get input in float form and padded with required zeroes
	getRealInput2D(d->inBuf, wp, wpitch,
		ht, wd, d->hbest, d->wbest, true);
	// forward fft transform
	d->fftwf_execute(d->pf);
	// search grid around points, apply limits to inner square
	for (int i = 0; i < d->npoints; i += 3)
		applyLimits(d, i);
	// inverse fft
	d->fftwf_execute(d->pinv);

	removeInputCentering(d->inBuf, d->wbest, d->hbest ); // also scales down as fft twice will scale up by wbest * hbest
	
	getRealOutput2D(d->inBuf, wp, wpitch, 
		ht, wd, d->hbest, d->wbest, min, max);
	
}
//................................................................
void applyLimits(F2QLimitData* d, int nf)
{
	float max = 0.0f;
	
	int centerOffset = (d->hbest / 2) * d->frqwidth + d->frqwidth / 2;

	int offMax = getOffsetOfMaxInSearchArea(d->outBuf + centerOffset,
		d->frqwidth, d->fspec[nf], d->fspec[nf + 1], d->xgrid, d->ygrid);	

	limitValuesInner(d->outBuf + centerOffset + offMax, 
						d->frqwidth, d->xinner, d->yinner, d->fspec[nf + 2]);
		// diagonally opposite
	limitValuesInner(d->outBuf + centerOffset - offMax,
		d->frqwidth, d->xinner,d->yinner, d->fspec[nf + 2]);
	
}
//......................................................................

//---------------------------------------------------------------------------------------------------------------------------
static const VSFrameRef *VS_CC f2qlimitGetFrame(int n, int activationReason, void **instanceData, 
						void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    F2QLimitData *d = (F2QLimitData *) * instanceData;

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
		
		// process all of RGB and Y of YUV or Gray
        int np = fi->colorFamily == cmRGB ? 3 : 1;
		int nbits = fi->sampleType == stInteger ? fi->bitsPerSample : 0;
		int nbytes = fi->bytesPerSample;

		for (int p = 0; p < np; p++)
		{
			const uint8_t* srcp = vsapi->getReadPtr(src, p);
			int dst_stride = vsapi->getStride(dst, p);
			uint8_t* dstp = vsapi->getWritePtr(dst, p);
			int pitch = dst_stride / fi->bytesPerSample;

			if (fi->sampleType == stInteger && nbits == 8)
			{
				uint8_t max = (1 << nbits) - 1, min = 0;				

				F2QLimitProcess(d, dstp, pitch, ht, wd,  min, max);
			}

			else if (fi->sampleType == stInteger && nbytes == 2)
			{
				uint16_t max = (1 << nbits) - 1, min = 0;
				uint16_t* dp = (uint16_t*)dstp;
				F2QLimitProcess(d, dp, pitch, ht, wd,  min, max);
			}

			else if (fi->sampleType == stFloat)
			{
				float max = 1.0, min = 0.0f;
				float* dp = (float*)dstp;
				F2QLimitProcess(d, dp, pitch, ht, wd,  min, max);
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
static void VS_CC f2qlimitFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    F2QLimitData *d = (F2QLimitData *)instanceData;
    vsapi->freeNode(d->node);
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
static void VS_CC f2qlimitCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
    F2QLimitData d;
    F2QLimitData *data;
    int err;
	int temp;
	
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

   
	if (!isConstantFormat(d.vi) && (d.vi->format->colorFamily != cmRGB 
						&& d.vi->format->colorFamily != cmYUV 	
						&& d.vi->format->colorFamily != cmGray))
	{
        vsapi->setError(out, "F2QLimit: only constant format RGB YUV or Gray  input supported");
        vsapi->freeNode(d.node);
        return;
    }
	
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "F2QLimit: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	
		
	d.grid = int64ToIntS(vsapi->propGetInt(in, "grid", 0, &err));
	if (err)
		d.grid = 10;
	else if (d.grid < 1 || d.grid > 50)
	{
		vsapi->setError(out, "F2QLimit: grid specifies search area and be 1 to 50 only");
		vsapi->freeNode(d.node);
		return;
	}
	d.inner = int64ToIntS(vsapi->propGetInt(in, "inner", 0, &err));
	if (err)
		d.inner = d.grid / 10;
	else if (d.inner < 0 || d.inner > d.grid )
	{
		vsapi->setError(out, "F2QLimit: inner area of limiting can be 0 to value of grid only");
		vsapi->freeNode(d.node);
		return;
	}
	d.warn = int64ToIntS(vsapi->propGetInt(in, "warn", 0, &err));
	if (err)
		d.warn = 1;
	else if (d.warn < 0 || d.warn > 2)
	{
		vsapi->setError(out, "F2QLimit: warn level can be 0 or 1 or 2 only");
		vsapi->freeNode(d.node);
		return;
	}

	d.npoints = vsapi->propNumElements(in, "fspec");

	if(d.npoints < 3 || (d.npoints % 3 ) != 0 || d.npoints > 60)
	{		
		vsapi->setError(out, "F2QLimit: fspec at least one and upto 20 filter specifications be given. Each filter is specified as a set of 5 integer values.");
		vsapi->freeNode(d.node);
		return;
	}

	for ( int i = 0; i < d.npoints; i ++)
		d.fspec[i] = int64ToIntS(vsapi->propGetInt(in, "fspec", i, 0));
	
	for (int i = 0; i < d.npoints; i += 3)
	{		
		//int wfreq = d.fspec[i];
		//int hfreq = d.fspec[i + 1];
		//int limit = d.fspec[i + 2];

		if (d.fspec[i] < 0 || d.fspec[i] >= NYQUIST / 2 )
		{
			vsapi->setError(out, "F2QLimit:the first number horizontal freq of each filter must  be from 0 to nyquist / 2 here  250");
			vsapi->freeNode(d.node);
			return;
		}
		if (d.fspec[i + 1] < - NYQUIST / 2 || d.fspec[i + 1] >= NYQUIST /2)
		{
			vsapi->setError(out, "F2QLimit:the second number vertical freq of each filter must be -nyqiust/2 to nyquist / 2 here -250 to 250");
			vsapi->freeNode(d.node);
			return;
		}

		if (d.fspec[i + 2] < 0 || d.fspec[i + 2] > 99)
		{
			vsapi->setError(out, "F2QLimit:the third number  specifying freq1 of each filter must be between 0 and %d");
			vsapi->freeNode(d.node);
			return;

		}
	}		
	
    data = (F2QLimitData *) malloc(sizeof(d));
    *data = d;
	
	vsapi->createFilter(in, out, "F2QLimit", f2qlimitInit, f2qlimitGetFrame, f2qlimitFree, fmParallelRequests, 0, data, core);
}



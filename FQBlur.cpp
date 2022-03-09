/******************************************************************************
FQBlur filter plugin for vapoursynth by V.C.Mohan
This filter operates in freq domain (2d) and blurs 
linear (motion) or circular (focus) styles within a window

Author V.C.Mohan. 
12 june 2015, 26 May 2021
Copyright (C) < 2008- 2021>  <V.C.Mohan>

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
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------





//----------------------------------------------------------------------------------------

typedef struct 
{
    VSNodeRef *node;
    const VSVideoInfo *vi;				
		bool line;			// true for line or false for circular
		int xcoord;			// right end x coordinate or blur radius in pixels
		int ycoord;			// right end y coordinate for line being true		
	
		int wbest;				// nearest higher width for best speed of fft
		int hbest;				// nearest higher height for best speed of fft
		int wbestUV;			// for YUV with subsamplingW
		int hbestUV;			// for YUV with subsamplingH
		float * FreqFilter;	// This is the place designed filter resides
		float* FreqFilterUV;
		fftwf_plan  pf, pinv;	// forward and inverse fft plans
		fftwf_plan  pfUV, pinvUV;	// for subsampled planes
		int bestR, bestRUV;				// width for real data will be half of bwd + 1.and used for filter design 

#include "fftLateBindingClassParams.cpp"
		float* inBuf;
		fftwf_complex* outBuf;
		
}F2QBlurData;	
/*************************************************/
template <typename finc>
void blurPlane2D(F2QBlurData* d, float* inBuf, fftwf_complex* outBuf, float* filter,
	fftwf_plan pf, fftwf_plan pinv, const finc* sp, finc* dp,
	int pitch, int height, int width, int bestY, int bestX, finc min, finc max);
void positionBlurFilter(fftwf_complex* fout, float* Filter, int bestx, int besty);
//---------------------------------------------------------------------------------
static void VS_CC f2qblurInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    F2QBlurData *d = (F2QBlurData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);

    const VSFormat *fi = d->vi->format;
       
			// frame dimensions
	int fwd = ( (d->vi->width + 3) >> 2) << 2;
    int fht = ( (d->vi->height + 3) >> 2) << 2;
	int fwdUV = fwd, fhtUV = fht;
				// filter dimensions	

	int *factorsbuf = vs_aligned_malloc <int>(sizeof( int) *64, 32);	//maximum 64 factors, first is factor, second is dividend to be factored. At 

	d->wbest = getBestDim(fwd + ADDSAFE, factorsbuf);

	d->hbest = getBestDim(fht + ADDSAFE, factorsbuf);
	int subH = 0;
	int subW = 0;

	if (fi->colorFamily == cmYUV)
	{
		subW = fi->subSamplingW;
		subH = fi->subSamplingH;

		fwdUV = fwd >> subW;
		fhtUV = fht >> subH;

		d->wbestUV = getBestDim(fwdUV + ADDSAFE, factorsbuf);
		d->hbestUV = getBestDim(fhtUV + ADDSAFE, factorsbuf);
	}
	else
	{
		d->wbestUV = d->wbest;
		d->hbestUV = d->hbest;
	}

	vs_aligned_free(factorsbuf);

	d->bestR = d->wbest/2 + 1;
	d->bestRUV = d->wbestUV / 2 + 1;

#include "ConstructorCodeForLateBindingfft.cpp"

	if (!ok)
	{
		vsapi->setError(out, "vcm.fqblur: could not load any of the dll or get required fnctions");
		if (d->hinstLib != NULL)
			FreeLibrary(d->hinstLib);
		vsapi->freeNode(d->node);		
		return;
	}

	int nbits=fi->bitsPerSample;	
	
	int isize = d->hbest * d->wbest;
	int fqsize = d->hbest * d->bestR;
	int isizeUV = d->hbestUV * d->wbestUV;
	int fqsizeUV = d->hbestUV * d->bestRUV;
	// buffers 
	d->inBuf = (float*)d->fftwf_malloc(sizeof(float) * isize);

	d->outBuf = (fftwf_complex*)d->fftwf_malloc (sizeof(fftwf_complex) * fqsize);// +1 is only a safeguard not really reqd

		
			// creates forward and inverse fft plans

	d->FreqFilter = (float*)d->fftwf_malloc(sizeof(float) * fqsize);

			//  forward for padded size complex to complex  
	d->pf = d->fftwf_plan_dft_r2c_2d(d->hbest, d->wbest, d->inBuf, d->outBuf,  FFTW_MEASURE);
			// inverse 

	d->pinv = d->fftwf_plan_dft_c2r_2d(d->hbest, d->wbest, d->outBuf, d->inBuf,  FFTW_MEASURE);


	if (subH != 0 || subW != 0)
	{
		d->FreqFilterUV = (float*)d->fftwf_malloc(sizeof(float) * fqsizeUV);
		d->pfUV = d->fftwf_plan_dft_r2c_2d(d->hbestUV, d->wbestUV, d->inBuf, d->outBuf, FFTW_MEASURE);
		d->pinvUV = d->fftwf_plan_dft_c2r_2d(d->hbestUV, d->wbestUV, d->outBuf, d->inBuf, FFTW_MEASURE);
	}
	else
	{
		d->FreqFilterUV = d->FreqFilter;
		d-> pfUV = d->pf;
		d->pinvUV = d->pinv;
	}

	int count = DrawPSF(d->inBuf, d->line, d-> xcoord, d-> ycoord, d->wbest, d->hbest, 0.0);	// we can add spike to mellow inversion
					// forward transform
	d->fftwf_execute(d->pf);
			// make it zero phase
	for( int i = 0; i < fqsize; i ++)
			
		d->outBuf[i][0] = sqrt (getAmpSquareOfComplex(d->outBuf + i));

	positionBlurFilter(d->outBuf, d->FreqFilter, d->wbest, d->hbest);

	if (fi->subSamplingH != 0 || fi->subSamplingW != 0)
	{
		if (d->line)
			DrawPSF(d->inBuf, d->line, d->xcoord >> subW, d->ycoord >> subH, d->wbestUV, d->hbestUV, 0.0);
		else
			DrawCircularPSFUV(d->inBuf, d->xcoord, d->wbestUV, d->hbestUV, subW, subH);
					// forward transform
		d->fftwf_execute(d->pfUV);
		// make it zero phase
		for (int i = 0; i < fqsizeUV; i++)

			d->outBuf[i][0] = sqrt(getAmpSquareOfComplex(d->outBuf + i));

		positionBlurFilter(d->outBuf, d->FreqFilterUV, d->wbestUV, d->hbestUV);

	}
 	
	//d->fftwf_free(inBuf);
	//d->fftwf_free(outBuf);

}
//--------------------------------------------------------------------------------------------
void positionBlurFilter(fftwf_complex* fout, float* Filter, int bestx, int besty)
{
	float scale = 1.0f / (bestx * besty);

	// transfer to filter buffer including scaling and removing filtering

	for (int h = 0; h < besty; h++)
	{
		for (int w = 0; w < bestx / 2 + 1; w++)
		{
			Filter[w] = fout[w][0] * scale;
		}
		Filter += bestx / 2 + 1;
		fout += bestx / 2 + 1;
	}

}
//--------------------------------------------------------------------------------------
template <typename finc>
void blurPlane2D(F2QBlurData * d, float * inBuf, fftwf_complex * outBuf, float * filter,
	fftwf_plan pf, fftwf_plan pinv, const finc* sp, finc* dp,
		 int pitch, int height, int width, int bestY, int bestX,  finc min, finc max )
{
	getRealInput2D(inBuf, sp, pitch, height, width, bestY, bestX, false);

	d->fftwf_execute_dft_r2c(pf, inBuf, outBuf);

	ApplyFilter2D(outBuf, filter, bestY, (bestX / 2 ) + 1);

	d->fftwf_execute_dft_c2r(pinv, outBuf, inBuf);

	getRealOutput2D(inBuf, dp, pitch,
		height, width, bestY, bestX, min, max);
}

//-----------------------------------------------------------------------------------------------

// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC f2qblurGetFrame(int n, int activationReason, void **instanceData, 
		void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    F2QBlurData *d = (F2QBlurData *) * instanceData;

    if (activationReason == arInitial) 
	{
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
        int height = vsapi->getFrameHeight(src, 0);
        int width = vsapi->getFrameWidth(src, 0);
        // When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
        // supply the "dominant" source frame to copy properties from. Frame props
        // are an essential part of the filter chain and you should NEVER break it.
        VSFrameRef *dst = vsapi->copyFrame(src, core);
		
		//float* inBuf = d->inBuf; // vs_aligned_malloc <float>(sizeof(float) * d->wbest * d->hbest, 32);

		//fftwf_complex* outBuf = d->outBuf;// vs_aligned_malloc<fftwf_complex>(sizeof(fftwf_complex) * (d->wbest / 2 + 1) * d->hbest, 64);

		// It's processing loop time!
        // Loop over all the planes
        int nplanes = fi->numPlanes > 3 ? 3 : fi->numPlanes;

        for (int plane = 0; plane < nplanes; plane ++) 
		{
			
            const uint8_t *sp = vsapi->getReadPtr(src, plane);
            int src_stride = vsapi->getStride(src, plane);
            uint8_t *dp = vsapi->getWritePtr(dst, plane);
            int dst_stride = vsapi->getStride(dst, plane); // note that if a frame has the same dimensions and format, the stride is guaranteed to be the same. int dst_stride = src_stride would be fine too in this filter.
            // Since planes may be subsampled you have to query the height of them individually
            int height = vsapi->getFrameHeight(src, plane);            
            int width = vsapi->getFrameWidth(src, plane);
			int nbytes = fi->bytesPerSample;
			int spitch = src_stride / nbytes;
			int dpitch = dst_stride / nbytes;
			
		// process image	
			if(fi->sampleType == stInteger)
			{
				if(fi->bitsPerSample == 8)
				{
					uint8_t min = 0;
					uint8_t max = 255;

					if (plane == 0)

						blurPlane2D(d, d->inBuf, d->outBuf, d->FreqFilter,
							d->pf, d->pinv, sp, dp,
							spitch, height, width, d->hbest, d->wbest, min, max);
					else
						blurPlane2D(d, d->inBuf, d->outBuf, d->FreqFilterUV,
							d->pfUV, d->pinvUV, sp, dp,
							spitch, height, width, d->hbestUV, d->wbestUV, min, max);
				}
				else
				{	// 16 bit data
					int nbits = fi->bitsPerSample;
					uint16_t min = 0;
					uint16_t max = (1 << nbits) - 1;

					const uint16_t *srcp = (uint16_t *) sp;
					uint16_t * dstp = (uint16_t *) dp;

					if (plane == 0)

						blurPlane2D(d, d->inBuf, d->outBuf, d->FreqFilter,
							d->pf, d->pinv, srcp, dstp,
							spitch, height, width, d->hbest, d->wbest, min, max);
					else

						blurPlane2D(d, d->inBuf, d->outBuf, d->FreqFilterUV,
							d->pfUV, d->pinvUV, srcp, dstp,
							spitch, height, width, d->hbestUV, d->wbestUV, min, max);
				}
			}
			else
			{
					// float data
					float min = 0.0f;
					float max = 1.0f;
					if (fi->colorFamily == cmYUV && plane != 0)
					{
						min = -0.5f;
						max = 0.5f;
					}
					const float *srcp = (float *) sp;
					float * dstp = (float *) dp;

					if (plane == 0)

						blurPlane2D(d, d->inBuf, d->outBuf, d->FreqFilter,
							d->pf, d->pinv, srcp, dstp,
							spitch, height, width, d->hbest, d->wbest, min, max);
					else

						blurPlane2D(d, d->inBuf, d->outBuf, d->FreqFilterUV,
							d->pfUV, d->pinvUV, srcp, dstp,
							spitch, height, width, d->hbestUV, d->wbestUV, min, max);
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
//----------------------------------------------------------------
// Free all allocated data on filter destruction
static void VS_CC f2qblurFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    F2QBlurData *d = (F2QBlurData *)instanceData;
	if (d->FreqFilterUV != d->FreqFilter)
		d->fftwf_free(d->FreqFilterUV);
	d->fftwf_free(d->FreqFilter);
	d->fftwf_free(d->inBuf);
	d->fftwf_free(d->outBuf);
	if (d->pf != d->pfUV)
		d->fftwf_destroy_plan(d->pfUV);
	d->fftwf_destroy_plan(d->pf);
	if (d->pinv != d->pinvUV)
		d->fftwf_destroy_plan(d->pinvUV);
	d->fftwf_destroy_plan(d->pinv);

	FreeLibrary(d->hinstLib);
    vsapi->freeNode(d->node);

    free(d);
}


/***************************************************************/


// This function is responsible for validating arguments and creating a new filter
static void VS_CC f2qblurCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) 
{
    F2QBlurData d;
    F2QBlurData *data;
    int err;
	int temp;
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    // In this first version we only want to handle 8bit integer formats. Note that
    // vi->format can be 0 if the input clip can change format midstream.
    if (!isConstantFormat(d.vi)  )
	{
        vsapi->setError(out, "F2QBlur: clip must have constant dimensions and in YUV or RGB or Grey format  ");
        vsapi->freeNode(d.node);
        return;
    }
	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "F2QBlur: RGB, YUV and Gray color formats only for input allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "F2QBlur: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string, the only
    // reason this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
    temp = !!int64ToIntS(vsapi->propGetInt(in, "line", 0, &err));
    if (err)
	{
        d.line = true;
	}
	else
    // Let's pretend the only allowed values are 1 or 0...
		if (temp < 0 || temp > 1)
	{
		vsapi->setError(out, "F2QBlur: line must be 0 (for circular blur) 1(for linear blur) ");
		vsapi->freeNode(d.node);
		return;
	}
	else
	{
		if( temp == 1)
			d.line = true;
		else
			d.line = false;
	}
	

	d.xcoord = int64ToIntS(vsapi->propGetInt(in, "x", 0, &err));
    if (err)
	{
        d.xcoord = 2;
	}
	else
    //  the only allowed values are 
		if ( (d.line && d.xcoord < 0) || d.xcoord > d.vi->width / 8 || ( !d.line && d.xcoord < 1))
	{
		vsapi->setError(out, "F2QBlur: x coordinate can have a value from 0 for line and 1 for circular blur to 1/8th frame width only ");
		vsapi->freeNode(d.node);
		return;
	}

	d.ycoord = int64ToIntS(vsapi->propGetInt(in, "y", 0, &err));
    if (err)
	{
        d.ycoord = 2;
	}
	else
    //  the only allowed values are 
		if ( d.ycoord <  - d.vi->height / 8 ||  d.ycoord  > d.vi->height / 8)
	{
		vsapi->setError(out, "F2QBlur: y coordinate can have a value between plus and minus 1/8th frame height only ");
		vsapi->freeNode(d.node);
		return;
	}
	if ( d.xcoord == 0 && d.ycoord == 0)
	{
		vsapi->setError(out, "F2QBlur: both x and y coordinate must not be zeroes ");
		vsapi->freeNode(d.node);
		return;
	}

    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = (F2QBlurData*)malloc(sizeof(d));
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

    vsapi->createFilter(in, out, "F2QBlur", f2qblurInit, f2qblurGetFrame, f2qblurFree, fmParallelRequests, 0, data, core);

}
/*
// The following function is the function that actually registers the filter in AviSynth
// It is called automatically, when the plugin is loaded to see which functions this filter contains.


    registerFunc("fqBlur", "clip:clip;line:int:opt;x:int:opt;y:int:opt;", Create_FQRestore, 0);
*/			

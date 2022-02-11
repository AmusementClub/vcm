/******************************************************************************
FQSharp filter plugin for vapoursynth by V.C.Mohan
This filter operates in freq domain (2d) and improves image having blurred image
 due to camera /motion.
 
   From the given x and y coordinates and type of blur an inverse filter using
   theoritical PSF is designed and applied.

  White noise, radius of inverse operator in spatial domain are to control 
  ringing. Scale is applied to bring image to acceptable levels.  

Author V.C.Mohan. 
11 jun 2015, 22 May 2021
  Copyright (C) <2008-2021>  <V.C.Mohan>

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
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
#include <math.h>
//#include <complex>
#include "fftwlite.h"

#include "factorize.cpp"
*/

typedef struct 
{
    VSNodeRef *node;
    const VSVideoInfo *vi;
	bool line;	//True: blur is linear . false : circular
    int xcoord, ycoord; // if linear blur line end coordinates (symmetrical about origin). If circular radius
	float wn;	// white noise to stabilize inversion.
	int frad;	// radius of designed filter if ham is true
	float scale; // scaling to compensate losses due to process. trial and error
	bool plane[3];
	bool ham;	// flag for hamming filter used along frad
	
	int wbest;				// nearest higher width for best speed of fft
	int hbest;				// nearest higher height for best speed of fft
	int fsize;				// input data size
	int fqsize;				// size after fft in freq domain 
	int frqwidth;
	float * FreqFilter;		// This is the place designed filter resides	
						
	fftwf_plan  pf, pinv;			// forward and inverse fft plans for full frames
#include "fftLateBindingClassParams.cpp"
	float* inBuf;
	fftwf_complex* outBuf;

} F2QSharpData;
//--------------------------------------------------------------------------------------
template <typename finc>
void sharpenPlane(F2QSharpData* d, const finc* sp, finc* dp, int pitch,
					int height, int width, finc min, finc max);
//----------------------------------------------------------------------------


// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC f2qsharpInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    F2QSharpData *d = (F2QSharpData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);

    const VSFormat *fi = d->vi->format;
	if (fi->colorFamily == cmYUV && (fi->subSamplingH != 0 || fi->subSamplingW != 0))
	{
		d->plane[1] = false;
		d->plane[2] = false;
	}
       
			// frame dimensions
	int fwd = d->vi->width;
    int fht = d->vi->height;

	if( d->ham)
			// filter radius
		d->frad = ((fwd > fht ? fht : fwd) * d->frad)/ 200;


	int *factorsbuf = new int[64];	//maximum 64 factors, first is factor, second is dividend to be factored. At 
	// make even number
	d->wbest = getBestDim( (((fwd + 3 ) >> 2 ) << 2) + ADDSAFE, factorsbuf);

	d->hbest = getBestDim((((fht + 3) >> 2) << 2) + ADDSAFE, factorsbuf);

	delete []factorsbuf;

	d->frqwidth = d->wbest/2 + 1;

	int nbits = fi->bitsPerSample;	
	
	d->fsize = d->hbest * d->wbest;
	int fqsize = d->hbest * d->frqwidth;

#include "ConstructorCodeForLateBindingfft.cpp"

	if (!ok)
	{
		vsapi->setError(out, "FQSharp: could not load any of the fft dll or get required fnctions");
		if (d->hinstLib != NULL)
			FreeLibrary(d->hinstLib);
		vsapi->freeNode(d->node);
		free(d);
		return;
	}
	 
	d->inBuf = (float*)d->fftwf_malloc(sizeof(float) * d->fsize);

	d->outBuf = (fftwf_complex*)d->fftwf_malloc (sizeof(fftwf_complex) * d->fqsize);// +1 is only a safeguard not really reqd
	
				// creates forward and inverse fft plans
	d->FreqFilter = (float*)d->fftwf_malloc(sizeof(float) * d->fqsize);
			//  forward for padded size complex to complex  
	d->pf = d->fftwf_plan_dft_r2c_2d(d->hbest, d->wbest, d->inBuf, d->outBuf,  FFTW_MEASURE);
			// inverse 

	d->pinv = d->fftwf_plan_dft_c2r_2d(d->hbest, d->wbest, d->outBuf, d->inBuf,  FFTW_MEASURE);

			// draw PSF for the specified   blur value. 

	DrawPSF(d->inBuf, d->line, d-> xcoord, d-> ycoord, d->wbest, d->hbest, d->wn);	// we can add spike wn to mellow inversion


				// forward transform
	d->fftwf_execute(d->pf);

			// make it zero phase

	for( int i = 0; i < d->fqsize; i ++)
	{			
		d->outBuf[i][0] = sqrt (getAmpSquareOfComplex(d->outBuf + i));
		d->outBuf[i][1] = 0.0;
	}		
		//	DesignInverse

	DesignInverse(d->outBuf, d->FreqFilter, d->wn, d->frqwidth, d->hbest,  d->scale);	

	if (d->ham)
	{
#include "hammingCodeInsert.cpp"
	}
/*	
		
			// hamming window as per filter radius in spatial domain
	F2QhammingWindowing(d->inBuf, d->wbest, d->wbest, d->hbest, d->frad);

			// forward transform into freq domain
	d->fftwf_execute(d->pf);
			// zero phase and also scale as it has gone twice transform
	float scaler = 1.0 / ( d->hbest * d->wbest);

	for(int i = 0; i < d->frqwidth * d->hbest; i ++)
		d->FreqFilter[i] = scaler * sqrt (getAmpSquareOfComplex( d->outBuf + i) ); 
*/		
				// completed initializing
}
//------------------------------------------------------------------------------------------------
template <typename finc>
void sharpenPlane(F2QSharpData* d, const finc * sp, finc * dp, int pitch,
					 int height, int width, finc min, finc max)					
{
	getRealInput2D(d->inBuf, sp, pitch, height, width, d->hbest, d->wbest, false);

	d->fftwf_execute(d->pf);

	ApplyFilter2D(d->outBuf, d->FreqFilter, d->hbest, d->frqwidth);

	d->fftwf_execute(d->pinv);

	getRealOutput2D(d->inBuf, dp, pitch, height, width,
		 d->hbest, d->wbest, min, max);
}

//-----------------------------------------------------------------------------------------------

// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC f2qsharpGetFrame(int n, int activationReason, void **instanceData, 
		void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    F2QSharpData *d = (F2QSharpData *) * instanceData;

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

        VSFrameRef *dst = vsapi->copyFrame( src, core);
		
        int nplanes = fi->numPlanes > 3 ? 3 : fi->numPlanes;

		for (int plane = 0; plane < nplanes; plane++)
		{
			if ( ! d->plane[plane]) continue;
			
            const uint8_t *sp = vsapi->getReadPtr(src, plane);
            uint8_t *dp = vsapi->getWritePtr(dst, plane);
            int stride = vsapi->getStride(dst, plane); // note that if a frame has the same dimensions and format, the stride is guaranteed to be the same. int dst_stride = src_stride would be fine too in this filter.
            // Since planes may be subsampled you have to query the height of them individually
            int height = vsapi->getFrameHeight(src, plane);            
            int width = vsapi->getFrameWidth(src, plane);
			int kb = fi->bytesPerSample;
			int pitch = stride / kb;
			
		// process image	
			if(fi->sampleType == stInteger)
			{
				if(fi->bitsPerSample == 8)
				{
					uint8_t min = 0;
					uint8_t max = 255;

					sharpenPlane(d, sp, dp, pitch, height, width, min, max);
				}
				else
				{	// 16 bit data
					int nb = fi->bitsPerSample;
					uint16_t min = 0;
					uint16_t max = (1 << nb) - 1;
					const uint16_t *srcp = (uint16_t *) sp;
					uint16_t * dstp = (uint16_t *) dp;
					
					sharpenPlane(d, srcp, dstp, pitch, height, width, min, max);
				}
			}
			else
			{
				// float data
				float min = 0.0f;
				float max = 1.0f;
				if (fi->colorFamily == cmYUV && plane > 0)
				{
					min = -0.5f;
					max = 0.5f;
				}

				const float* srcp = (float*)sp;
				float* dstp = (float*)dp;

				sharpenPlane(d, srcp, dstp, pitch, height, width, min, max);
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
static void VS_CC f2qsharpFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    F2QSharpData *d = (F2QSharpData *)instanceData;
	
	d->fftwf_destroy_plan(d->pf);
	d->fftwf_destroy_plan(d->pinv);
    vsapi->freeNode(d->node);
	// release buffers
	d->fftwf_free(d->inBuf);
	d->fftwf_free(d->outBuf);
	d->fftwf_free(d->FreqFilter);	
	FreeLibrary(d->hinstLib);
    free(d);
}
//......................................................................................

// This function is responsible for validating arguments and creating a new filter
static void VS_CC f2qsharpCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    F2QSharpData d;
    F2QSharpData *data;
    int err;
	int temp;
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    // In this first version we only want to handle 8bit integer formats. Note that
    // vi->format can be 0 if the input clip can change format midstream.
	if (!isConstantFormat(d.vi) && d.vi->format->colorFamily != cmRGB
		&& d.vi->format->colorFamily != cmYUV
		&& d.vi->format->colorFamily != cmGray)
	{
        vsapi->setError(out, "F2QSharp: Input clip must have constant dimensions and in YUV or RGB or Grey format");
        vsapi->freeNode(d.node);
        return;
    }

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string, the only
    // reason this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
    temp = !!vsapi->propGetInt(in, "line", 0, &err);
    if (err)
	{
        d.line = false;
	}
	else
    // Let's pretend the only allowed values are 1 or 0...
		if (temp < 0 || temp > 1)
	{
		vsapi->setError(out, "F2QSharp: line must be 0 (for circular blur) 1(for linear blur) ");
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
	d.wn = (float)vsapi->propGetFloat(in, "wn", 0, &err);
    if (err)
	{
        d.wn = 0.05f;
	}
	else if (d.wn < 0.0001f || d.wn > 0.99f)
	{
		vsapi->setError(out, "F2QSharp: white noise wn value can only be between 0.0001 and 0.99  ");
		vsapi->freeNode(d.node);
		return;
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
		vsapi->setError(out, "F2QSharp: x coordinate can have a value from 0 for line and 1 for circular blur to 1/8th frame width only ");
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
		vsapi->setError(out, "F2QSharp: y coordinate can have a value between plus and minus 1/8th frame height only ");
		vsapi->freeNode(d.node);
		return;
	}
	if ( d.xcoord == 0 && d.ycoord == 0)
	{
		vsapi->setError(out, "F2QSharp: both x and y coordinate must not be zeroes ");
		vsapi->freeNode(d.node);
		return;
	}
	
	temp = !!int64ToIntS(vsapi->propGetInt(in, "ham", 0, &err));
	d.ham = err ||temp == 0 ? false : true;

	if (d.ham)
	{
		d.frad = int64ToIntS(vsapi->propGetInt(in, "frad", 0, &err));
		if (err)
		{
			d.frad = 30;
		}
		else if (d.frad < 10 || d.frad > 50)
		{
			vsapi->setError(out, "F2QSharp: filter radius %age of smaller dimension of frame, can have a value of 10 to 50 only ");
			vsapi->freeNode(d.node);
			return;
		}
	}

	d.scale = (float)vsapi->propGetFloat(in, "scale", 0, &err);
    if (err)
	{
        d.scale = 0.45f;
	}
	else if (d.scale < 0.000001f || d.scale > 1000000.0f)
	{
		vsapi->setError(out, "F2QSharp: scale value must be between 0.000001 and 1000000");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->colorFamily == cmRGB)
	{
		int ntemp = vsapi->propNumElements(in, "rgb");
		if (ntemp == 0)
		{
			d.plane[0] = true;
			d.plane[1] = true;
			d.plane[2] = true;
		}
		else if ( ntemp > 3)
		{ 
			vsapi->setError(out, "F2QSharp: rgb array can not have more than 3 values");
			vsapi->freeNode(d.node);
			return;
		}
		else
		{
			temp = !!int64ToIntS(vsapi->propGetInt(in, "rgb", 0, &err));
			d.plane[0] = temp == 0 ? false : true;
		}

		for (int i = 1; i < 3; i++)
		{
			temp = !!int64ToIntS(vsapi->propGetInt(in, "rgb", i, &err));
			if (err)
				d.plane[i] = d.plane[i - 1];
			else
				d.plane[i] = temp == 0 ? false : true;

		}
		// as rgb planes are in the order of bgr internally correct here
		bool rgb = d.plane[0];
		d.plane[0] = d.plane[2];
		d.plane[2] = rgb;
	}

	else if (d.vi->format->colorFamily == cmYUV)
	{
		int ntemp = vsapi->propNumElements(in, "yuv");
		if (ntemp == 0)
		{
			d.plane[0] = true;
			d.plane[1] = false;
			d.plane[2] = false;
		}
		else if (ntemp > 3)
		{
			vsapi->setError(out, "F2QSharp: yuv array can not have more than 3 values");
			vsapi->freeNode(d.node);
			return;
		}
		else
		{
			temp = !!int64ToIntS(vsapi->propGetInt(in, "yuv", 0, &err));
			d.plane[0] = temp == 0 ? false : true;
		}

		for (int i = 1; i < 3; i++)
		{
			temp = !!int64ToIntS(vsapi->propGetInt(in, "yuv", i, &err));
			if (err)
				d.plane[i] = d.plane[i - 1];
			else
				d.plane[i] = temp == 0 ? false : true;

		}
	}

    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = (F2QSharpData*)malloc(sizeof(d));
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

	vsapi->createFilter(in, out, "F2QSharp", f2qsharpInit, f2qsharpGetFrame, f2qsharpFree, fmParallelRequests, 0, data, core);

}
/*
// The following function is the function that actually registers the filter in AviSynth
// It is called automatically, when the plugin is loaded to see which functions this filter contains.


    registerFunc("Sharp", "clip:clip;line:int:opt;wn:float:opt;x:int:opt;y:int:opt;ham:int:opt;frad:int:opt;scale:float:opt;rgb:int[]:opt;yuv:int[]:opt;", Create_FQSharp, 0);
					
    


*/
/******************************************************************************
(Adaptive)median filter plugin for Vapoursynth by V.C.Mohan

  Using the maximum, minimum and median values in a grid using the rules:
  1. if median is more than minimum and less than maximium
			then if the grid center value is >than min and < than max,the value is not changed
			else  it is replaced by median.  in each grid of 3 x 3 pixels
  2 else if grid size is less than maxgrid, it is increased to next step and process repeated.
		else value is unchanged.
 No input clip format restriction
   Good for additive noise  gaussian or salt and pepper noise elimination. Preserves edges,
   removes impulsive noises at the same time smoothing image.
   Rules implemented for processing. 
starting value of grid is 3 i.e. 3x3 pixels. 
Compute max, min and median of the grid 
1. if median is more than minimum and less than maximium:- 
 if the grid center value is greater than minimum and less than maximum,the value is not changed 
 Otherwise center value of grid is replaced by median. 
2 else:- 
 if grid size is less than maxgrid, it is increased to next step and process repeated. 
 Otherwise center value is unchanged.

As per theory 
Quote "Ordinary median filters perform well as long as the spatial density of impulsive noise is small. Adaptive Median filter can handle impulsive noise having larger probablity. An additional benefit is this seeks to preserve detail while smoothing nonimpulse noise something that the traditional median filter does not do.

The algorithm described has three purposes. 
1.To remove salt and pepper (impulse) noise. 
2.To smooth other noise which may not be impulsive 
3.To reduce distortionsuch as excessive thinning or thickening of object boundaries. 
" unquote


 
 Author V.C.Mohan

  Aug 22, 2020, 8 nov 2021
*************************************************************************************************/ 

//////////////////////////////////////////
// This file contains  adaptivemedian

/*
#include <stdlib.h>
#include <algorithm>
#include "VapourSynth.h"
#include "VSHelper.h"
*/
//--------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
template <typename finc>
void getMinAndMax(finc* minmaxmed, const finc* tbuf, const int length);
template <typename finc>
void AdMed(const finc* fp, finc* wp, const int pitch,
	const int bwd, const int bht, const int mingrid, const int maxgrid, int * LUT,
	const finc min, const finc max);

int createOffsetLUT(int* LUT, int pitch, int kb, int maxGrid);
template <typename finc>
void fillUpBuffer(const finc* pindex, int sIndex, int eIndex, const int* LUT, finc* tbuf);

//............................................................................................

int createOffsetLUT(int* LUT, int pitch, int kb, int maxGrid)
{
	// creates offsets for  grid. By negating second half is obtained
	int index = 0;
	for (int g = 1; g <= maxGrid; g += 2)
	{
		if (g == 1)
		{
			// center point
			LUT[index] = 0;
			index++;
		}
		else
		{
			
			for (int i = -g / 2; i < g / 2; i++)
			{
				LUT[index] = (pitch * (g / 2) + kb * i);
				index++;
				LUT[index] = -LUT[index - 1];
				index++;
			
				LUT[index] = -kb * (g / 2) + pitch * i;
				index++;
				LUT[index] = -LUT[index - 1];
				index++;
			}
		}
	}
	return index;

}
template <typename finc>
void getMinAndMax(finc* minmaxmed, const finc * tbuf, const int length)
{
	finc minimum = minmaxmed[0];
	finc maximum = minmaxmed[1];


	for (int i = 0; i < length; i++)
	{
		if (minimum > tbuf[i])
			minimum = tbuf[i];
		if (maximum < tbuf[i])
			maximum = tbuf[i];
	}
	minmaxmed[0] = minimum;
	minmaxmed[1] = maximum;
}


template <typename finc>
void fillUpBuffer(const finc *pindex, int sIndex, int eIndex, const int * LUT, finc * tbuf)
{
	for (int i = sIndex; i < eIndex; i++)
	{
		tbuf[i] = *(pindex + LUT[i]);
	}

}


template<typename finc>
void AdMed ( const finc *fp, finc  *wp, const int pitch, 
			const int bwd, const int bht, const int mingrid,const int maxgrid, int * LUT,
	const finc min, const finc max)
{
	finc tolr = (finc)(max * 0.005f);

	finc* tbuf = (finc*)vs_aligned_malloc( sizeof(finc) * maxgrid * maxgrid , 32);

	int gsize = mingrid;// will be increased if required as adaptation
	

	for(int h = mingrid/2 ;	h < bht - mingrid/2; h ++)
	{
		for(int w = mingrid/2; w < bwd - mingrid/2; w ++)
		{	
			finc minmaxmed[] = { max, min, min }; // values will be changed in processing
			int maxgsize = maxgrid;	// will be modified to possible value					

			// get the maximum size to which grid can increase at this point
			// At frame margins it may not increase to max size
			
			while( h < maxgsize / 2 || bht - h - 1 < maxgsize / 2 
					|| w < maxgsize / 2 || bwd - w - 1 < maxgsize / 2)
				
				maxgsize -= 2;				
			
			gsize = mingrid;
			
			while (gsize <= maxgsize)
			{
				int start = gsize == 3? 0: (gsize - 2) * (gsize - 2);
				
				int full = gsize * gsize;
				
				fillUpBuffer(fp + h * pitch + w, start, full, LUT, tbuf);
				// get minimum, maximum
				getMinAndMax<finc>(minmaxmed, tbuf + start, full);

				if (minmaxmed[1] - minmaxmed[0] <= tolr)
				{					
					break;	// uniform color
				}

				std::nth_element(tbuf, tbuf + full / 2, tbuf + full - 1);
				minmaxmed[2] = tbuf[full / 2];
				
				// check if median is within max and min value or at the extreme				
				finc pval = fp[h * pitch + w];  // value at center of grid

				 if (minmaxmed[ 2] > minmaxmed[0] && minmaxmed[ 2] < minmaxmed[ 1])
				 {
					 if (pval > minmaxmed[0] && pval < minmaxmed[1])
					 {
						
						 break;		//  looks good. retain original value
					 }
				
					 // looks like salt or pepper or noise. So replace with median					
					wp[h * pitch + w] = (finc)minmaxmed[2];
					break;
					
				 }
				 
				gsize += 2;
			}
		
		}	// for int w =

	}	//for h=
	vs_aligned_free(tbuf);
}

//---------------------------------------------------------------------------------------------

typedef struct {
				VSNodeRef *node;
				const VSVideoInfo *vi;
				int maxgrid;
				int yy; 
				int uu; 
				int vv;	
				int* LUT;
				int* uvLUT;
				bool uvbuf;
				} AdaptiveMedianData;

// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC adaptivemedianInit
				(VSMap *in, VSMap *out, void **instanceData, 
				VSNode *node, VSCore *core, const VSAPI *vsapi) 
{
    AdaptiveMedianData *d = (AdaptiveMedianData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);
	d->LUT = (int*)vs_aligned_malloc<int>(sizeof(int) * d->maxgrid * d->maxgrid, 32);
	const VSFormat* fi = d->vi->format;
	int nbytes = fi->bytesPerSample;	
	
	const VSFrameRef* src = vsapi->getFrame(0, d->node, NULL, 0);
	
	const int pitch = vsapi->getStride(src, 0) / nbytes;
	
		// create LUT
	int count = createOffsetLUT(d->LUT, pitch, 1, d->maxgrid);

	d->uvbuf = false;

	if (d->uu != 0 || d->vv != 0)
	{
		if ( fi->subSamplingH == 0 || fi->subSamplingW == 0)
		
			d->uvLUT = d->LUT;
		else
		{
			const int uvpitch = vsapi->getStride(src, 1) / nbytes;
			d->uvLUT = (int*)vs_aligned_malloc<int>(sizeof(int) * d->maxgrid * d->maxgrid, 32);
			d->uvbuf = true;

			int count = createOffsetLUT(d->uvLUT, uvpitch, 1, d->maxgrid);
		}

	}
	vsapi->freeFrame(src);	
	
}

// This is the main function that gets called when a frame should be produced. It will in most cases get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all input frames you need. Always do it i ascending order to play nice with the
// upstream filters.
// Once all frames are ready the the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC adaptivemedianGetFrame
				(int n, int activationReason, void **instanceData, 
				void **frameData, VSFrameContext *frameCtx, 
				VSCore *core, const VSAPI *vsapi)
{
    AdaptiveMedianData *d = (AdaptiveMedianData *) * instanceData;

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
		int opt = fi ->colorFamily == cmRGB? 7 :  (((d->vv << 2) |  d->uu) << 1) | d->yy; 
		
        // When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
        // supply the "dominant" source frame to copy properties from. Frame props
        // are an essential part of the filter chain and you should NEVER break it.
        VSFrameRef *dst = vsapi->copyFrame( src, core);

        // It's processing loop time!
        // Loop over all the planes
        
        for (int plane = 0; plane < fi->numPlanes; plane++)
		{

			const uint8_t *srcp = vsapi->getReadPtr(src, plane);
            int spitch = vsapi->getStride(src, plane);
            uint8_t *dstp = vsapi->getWritePtr(dst, plane);
			int bht = vsapi->getFrameHeight(src, plane);            
            int bwd = vsapi->getFrameWidth(src, plane);
			int nbytes = fi->bytesPerSample;
			int pitch = spitch / nbytes;
			int nbits = fi->bitsPerSample;
			
			if (fi->sampleType == stInteger)
			{
				if(nbytes == 1)
				{
					// 8 bit integer samples. used to get only min and max in grid					
					uint8_t min = 0, max = (uint8_t)255;
					
					int mingrid = 3;

					if(((opt >> plane) && 1) == 1  )
					{
						if ( plane == 0 || fi->colorFamily == cmRGB)
							AdMed ( srcp, dstp, pitch, bwd, bht, mingrid, d->maxgrid, d->LUT, min, max);
						else
							AdMed(srcp, dstp, pitch, bwd, bht, mingrid, d->maxgrid, d->uvLUT, min, max);

					}
				}

			
				else if (nbytes == 2)	// 9 to 16 bit integer samples
				{

					const uint16_t *sp = (uint16_t *) srcp;
					// only max and min possible. Actual does not affect
					uint16_t min = 0, max = (uint16_t)(( 1 << nbits) - 1);
					
					uint16_t *dp = (uint16_t*) dstp;
					
					int mingrid = 3;
					if (((opt >> plane) && 1) == 1)
					{
						if (plane == 0 || fi->colorFamily == cmRGB)
							AdMed(sp, dp, pitch, bwd, bht, mingrid, d->maxgrid, d->LUT, min, max);
						else
							AdMed(sp, dp, pitch, bwd, bht, mingrid, d->maxgrid, d->uvLUT, min, max);
					}

				}	//16 bit integer

			}	// integer formats

			else	// float
			{
				
				const float *sp = (float *) srcp;
           
				float *dp = (float*) dstp;
				float min = -0.5f, max = 1.0f;				
				int mingrid = 3;
				if (((opt >> plane) && 1) == 1)
				{
					if (plane == 0 || fi->colorFamily == cmRGB)
						AdMed(sp, dp, pitch, bwd, bht, mingrid, d->maxgrid, d->LUT, min, max);
					else
						AdMed(sp, dp, pitch, bwd, bht, mingrid, d->maxgrid, d->uvLUT, min, max);
				}

			}	// float

		}	// for plane


        // Release the source frame
        vsapi->freeFrame(src);

        // A reference is consumed when it is returned so saving the dst ref somewhere
        // and reusing it is not allowed.
		
        return dst;
    }

    return 0;
}

// Free all allocated data on filter destruction
static void VS_CC adaptivemedianFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    AdaptiveMedianData *d = (AdaptiveMedianData *)instanceData;
    vsapi->freeNode(d->node);
	vs_aligned_free(d->LUT);
	if (d->uvbuf)
		vs_aligned_free(d->uvLUT);

    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC adaptivemedianCreate(const VSMap *in, 
					VSMap *out, void *userData, VSCore *core,
					const VSAPI *vsapi) 
{
    AdaptiveMedianData d;
    AdaptiveMedianData *data;
    int temp;
    int err;

    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "Median: RGB, YUV and Gray color formats only for input allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "Median: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string the only reason
    // this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
    
	d.maxgrid = int64ToIntS(vsapi->propGetInt(in, "maxgrid", 0, &err));

    if (err)
        d.maxgrid = 5;

    // the only allowed values are 3 to 9 odd numbers...
    if (d.maxgrid < 3 || d.maxgrid > 11 || (d.maxgrid & 1) == 0) 
	{
        vsapi->setError(out, "Median: maxgrid value can be odd number 3 to 11 only");
        vsapi->freeNode(d.node);
        return;
    }

	
	temp = 0;
	temp = vsapi->propNumElements(in, "plane");
	if (temp <= 0)
	{
		d.yy = 1;
		if (d.vi->format->colorFamily == cmYUV)
		{
			d.uu = 0;
			d.vv = 0;
		}
		else
		{
			d.uu = 1;
			d.vv = 1;
		}
	}
	else if (temp < 3)
	{
		vsapi->setError(out, "Median: values of each of 3 planes as one or zero must be specified");
		vsapi->freeNode(d.node);
		return;
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			if (int64ToIntS(vsapi->propGetInt(in, "plane", i, 0)) < 0
				|| int64ToIntS(vsapi->propGetInt(in, "plane", i, 0)) > 1)
			{
				vsapi->setError(out, "Median:  value of each of 3 planes as one or zero must be specified");
				vsapi->freeNode(d.node);
				return;
			}
		}


		d.yy = int64ToIntS(vsapi->propGetInt(in, "plane", 0, 0));
		d.uu = int64ToIntS(vsapi->propGetInt(in, "plane", 1, 0));
		d.vv = int64ToIntS(vsapi->propGetInt(in, "plane", 2, 0));
	}


	if (d.yy == 0 && d.uu == 0 && d.vv == 0)
	{
		vsapi->setError(out, "Median: At least one of plane must be set to 1 ");
		vsapi->freeNode(d.node);
		return;
	}
	
	
    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
	
    data =   (AdaptiveMedianData * ) malloc(sizeof(d));
    *data = d;

    // Create a new filter and returns a reference to it. Always pass on the in and out
    // arguments or unexpected things may happen. The name should be something that's
    // easy to connect to the filter, like its function name.
    // The three function pointers handle initialization, frame processing and filter destruction.
    // The filtermode is very important to get right as it controls how threading of the filter
    // is handled. In general you should only use fmParallel whenever possible. This is if you
    // need to modify no shared data at all when the filter is running.
    // For more complicated filters fmParallelRequests is usually easier to achieve as an
    // be prefetched in parallel but the actual processing is serialized.
    // The others can be considered special cases where fmSerial is useful to source filters and
    // fmUnordered is useful when a filter's state may change even when deciding which frames to
    // prefetch (such as a cache filter).
    // If you filter is really fast (such as a filter that only resorts frames) you should set the
    // nfNoCache flag to make the caching work smoother.
    vsapi->createFilter(in, out, "Median", adaptivemedianInit, 
								adaptivemedianGetFrame, adaptivemedianFree, 
								fmParallel, 0, data, core);
    return;
}

//////////////////////////////////////////
// Init

// This is the entry point that is called when a plugin is loaded. You are only supposed
// to call the two provided functions here.
// configFunc sets the id, namespace, and long name of the plugin (the last 3 arguments
// never need to be changed for a normal plugin).
//
// id: Needs to be a "reverse" url and unique among all plugins.
//   It is inspired by how android packages identify themselves.
//   If you don't own a domain then make one up that's related
//   to the plugin name.
//
// namespace: Should only use [a-z_] and not be too long.
//
// full name: Any name that describes the plugin nicely.
//
// registerFunc is called once for each function you want to register. Function names
// should be PascalCase. The argument string has this format:
// name:type; or name:type:flag1:flag2....;
// All argument name should be lowercase and only use [a-z_].
// The valid types are int,float,data,clip,frame,func. [] can be appended to allow arrays
// of type to be passed (numbers:int[])
// The available flags are opt, to make an argument optional, empty, which controls whether
// or not empty arrays are accepted and link which will not be explained here.
/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("vc.mohan.dns", "adapt", "VapourSynth AdaptiveMedian", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("AdaptiveMedian", "clip:clip;maxgrid:int:opt;yy:int:opt;uu:int:opt;vv:int:opt;", adaptivemedianCreate, 0, plugin);
}
*/
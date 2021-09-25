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

  Aug 22, 2020
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
void AdMed ( const finc * , finc  *, const int , const	int , const int , const  int , const int, finc tol );

template<typename finc>
void AdMed ( const finc *fp, finc  *wp, const int pitch, 
			const int bwd, const int bht, const int mingrid,const int maxgrid)
{

	finc lbuf[100];

	int gsize = mingrid;// will be increased if required as adaptation

	for(int h = mingrid/2 ;	h < bht - mingrid/2; h ++)
	{

		for(int w = mingrid/2; w < bwd - mingrid/2; w ++)
		{						
			int maxgsize = maxgrid;	// will be modified to possible value					

			// get the maximum size to which grid can increase at this point
			
			while( h < maxgsize / 2 || bht - h - 1 < maxgsize / 2 
					|| w < maxgsize / 2 || bwd - w - 1 < maxgsize / 2)
				
				maxgsize --;				
			
			gsize = mingrid;

			while (gsize <= maxgsize	)
			{
				// fill the buffer with grid values

				int i = 0;

				int hhindex = ( h - gsize / 2) * pitch;

				for(int hh =  h - gsize / 2; hh <= h + gsize / 2 ; hh ++)
				{
					//int hhindex = hh * fpitch

					for(int ww = w - gsize / 2; ww <= w + gsize / 2; ww ++)
					{

						lbuf[i ] = fp[hhindex + ww];
						i ++;
					}

					hhindex += pitch;

				}

				// now sort in ascending order			
				std:: sort( lbuf, lbuf + i );
				// check if median is within max and min value or at the extreme
				finc pval = fp[h * pitch + w];  // value at center of grid
			
			//	 finc lmin = lbuf[0];
			//	 finc lmax = lbuf[i-1];
			//	 finc lmed = lbuf[i / 2];
				 if (lbuf[0] == lbuf[i - 1])
					 break;	// uniform color

				 if (lbuf[i / 2] > lbuf[0] && lbuf[i / 2] < lbuf[i - 1])
				 {
					 if (pval > lbuf[0] && pval < lbuf[i - 1])
						 break;		// retain original value
				
					 // looks like salt or pepper or noise. So replace with med
					 wp[h * pitch + w] = lbuf[i / 2];
					 break;

				 }
				 
				gsize += 2;
			}

		
		}	// for int w =

	}	//for h=
}

//---------------------------------------------------------------------------------------------

typedef struct {
				VSNodeRef *node;
				const VSVideoInfo *vi;
				int maxgrid;
				int yy; 
				int uu; 
				int vv;				
				
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
		char opt = fi ->colorFamily == cmRGB? 0x111 :  (((d->vv << 1) |  d->uu) << 1) | d->yy; 
		
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
			int samplesize = fi->bytesPerSample;
			int pitch = spitch / samplesize;
			
			if (fi->sampleType == stInteger)
			{
				if(fi->bitsPerSample == 8)
				{
					// 8 bit integer samples					
					
					int mingrid = 3;

					if(((opt >> plane) & 1) == 1  )
					{
						AdMed ( srcp, dstp, pitch, bwd, bht, mingrid, d->maxgrid);						
					}
				}

			
				else	// 10 or 16 bit integer samples
				{

					const uint16_t *sp = (uint16_t *) srcp;
					
					uint16_t *dp = (uint16_t*) dstp;
					
					int mingrid = 3;

					if(((opt >> plane) & 1) == 1  )
					{
						AdMed ( sp, dp, pitch,  bwd, bht,  mingrid, d->maxgrid);
					}

				}	//16 bit integer

			}	// integer formats

			else	// float
			{
				
				const float *sp = (float *) srcp;
           
				float *dp = (float*) dstp;
							
				int mingrid = 3;

				if(((opt >> plane) & 1) == 1  )
				{
					AdMed ( sp, dp, pitch,  bwd, bht,  mingrid, d->maxgrid);					
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

	if (d.vi->format->colorFamily == cmCompat)
	{
		vsapi->setError(out, "median: Compat format is not supported");
		vsapi->freeNode(d.node);
		return;
	}

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string the only reason
    // this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
    
	d.maxgrid = vsapi->propGetInt(in, "maxgrid", 0, &err);

    if (err)
        d.maxgrid = 5;

    // the only allowed values are 3 to 9 odd numbers...
    if (d.maxgrid < 3 || d.maxgrid > 9 || (d.maxgrid & 1) == 0) 
	{
        vsapi->setError(out, "median: maxgrid value can be  3 or 5 or 7 or 9 only");
        vsapi->freeNode(d.node);
        return;
    }
	temp = 0;
	temp = vsapi->propNumElements(in, "plane");
	if( temp <= 0)
	{
		d.yy = 1, d.uu = 0, d.vv = 0;
	}
	else if(temp < 3)
	{
		vsapi->setError(out, "median: values of each of 3 planes as one or zero must be specified");
        vsapi->freeNode(d.node);
        return;
    }
	else
	{
		for ( int i = 0; i < 3; i ++)
		{
			if( vsapi->propGetInt(in, "plane", i, 0) < 0 || vsapi->propGetInt(in, "plane", i, 0) > 1)
			{
				vsapi->setError(out, "median:  value of each of 3 planes as one or zero must be specified");
				vsapi->freeNode(d.node);
				return;
			}
		}
	

		d.yy = vsapi->propGetInt(in, "plane", 0, 0);
		d.uu = vsapi->propGetInt(in, "plane", 1, 0);
		d.vv = vsapi->propGetInt(in, "plane", 2, 0);
	}
   
    
    if (d.yy == 0 && d.uu == 0 && d.vv == 0)  
	{
        vsapi->setError(out, "median: At least one of plane must be set to 1 ");
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
    vsapi->createFilter(in, out, "median", adaptivemedianInit, 
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
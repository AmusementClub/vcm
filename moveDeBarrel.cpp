/******************************************************************************
DeBarrel filter plugin for vapoursynth by V.C.Mohan
DeBarrel corrects input image from barrel or pin cushion type distortion
given parametrs a, b, c in mode 1 and c in mode 2. In test mode draws barrel or pincushion type 
lines on input frame corresponding to the constants . This feature enables to determine values for a,b and or c.

distortion where  x axis distortion  is independant of y axis distortion also can be corrected.
(included at a request from an user)

vhr parameter corrects for unequal distortions along height and width of frame 
 seen on some anamorphic projections.
 
 Author V.C.Mohan

  sep 2014, modified on 20 Aug 2020 20 Oct 2021

Copyright (C) <2020 - 2021>  <V.C.Mohan>

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
#include "math.h"
#include "interpolationMethods.h"
*/
typedef struct {

		VSNodeRef *node;
		const VSVideoInfo *vi;
				
		float abc[3];
		int q;
		bool pin;	// whether pincushion type  or  barrel type
		bool test;	// whether a test run		
		int method;		// formula type 1. uses a,b and c 2. Uses c ( and yc) only			
		int dots;	// in test
		int ddensity;
		float dim;
		
		float * iCoeff;	// interpolation coefficients	
		int quantile;		// accuracy of interpolation
		int span;			// span of lanczos 6 (X6)
		int* xyAndQ;
		int nEntries;
		
		uint8_t col[12];
	
}DeBarrelData;	

	
/*--------------------------------------------------
 * The following is the implementation 
 * of the defined functions.
 --------------------------------------------------*/
//Here is the init code used			
			
static void VS_CC debarrelInit(VSMap *in, VSMap *out, void **instanceData,
	VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    DeBarrelData *d = (DeBarrelData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);

	// input frame dimensions
	const int width = d->vi->width;
	const int height = d->vi->height;
	const VSFormat* fi = d->vi->format;
	int nbytes = fi->bytesPerSample;
	int nbits = fi->bitsPerSample;
	d->quantile = 64;
	d->nEntries = d->test ? 2 : d->q == 1 ? 3 : 4;

	d->xyAndQ = (int*)vs_aligned_malloc<int>(sizeof(int) * (width / 2) * (height / 2) * d->nEntries, 32);

	int* xyQ = d->xyAndQ;
	float xy[2];
	int x, y, qx, qy;

	d->iCoeff = NULL;

	if (!d->test)
		d->iCoeff = setInterpolationScheme(d->q, d->quantile, &d->span);	
	
	
	if (!d->pin)
	{
		for (int i = 0; i < 3; i++)
			d->abc[i] = -d->abc[i];			
	}
	float rNorm = (float)sqrt((double)(width * width + height * height) / 4);
	// create a look up table
	for ( int h = 0; h < height / 2; h ++)
	{
		int hoff = h * width / 2 *  d->nEntries;

		for( int w = 0; w < width / 2; w ++)
		{
			int woff = w * d->nEntries;
			getSourceCoord(xy, (float) w, (float)h, d->method, rNorm, d->abc);

			x = (int)floor(xy[0]);
			y = (int)floor(xy[1]);

			
			if (x >= width / 2 || y >= height / 2 || x < 0 || y < 0)
			{
				xyQ[hoff + woff] = -1;
			}
			else
			{
				xyQ[hoff + woff] = x;
				xyQ[hoff + woff + 1] = y;

				if (!d->test)
				{
					// calculate nearest quantile of the fraction
					qx = (int)((xy[0] - x) * d->quantile);
					qy = (int)((xy[1] - y) * d->quantile);				
				
					if (d->q > 1)
					{
						xyQ[hoff + woff + 2] = qx;
						xyQ[hoff + woff + 3] = qy;
					}
					else
						// index value
						xyQ[hoff + woff + 2] = bestOfNineIndex(qx, qy, d->quantile);
				}
			}
		}
	}

	// color to blacken out of area points
	uint8_t bgr[] = { 0,0,0 }, yuv[] = { 16,128,128 };

	if (d->test)
	{
		// will have white dots
		d->ddensity = (5 - d->dots) * 16;
		bgr[0] = 255;
		bgr[1] = 255;
		bgr[2] = 255;
	}

	convertBGRforInputFormat(d->col, bgr, fi);

}

	// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC debarrelGetFrame(int n, int activationReason, void **instanceData, 
	void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
{
    DeBarrelData *d = (DeBarrelData *) * instanceData;

    if (activationReason == arInitial)
	{
		
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    }
	else if (activationReason == arAllFramesReady) 
	{
        const VSFrameRef *src = vsapi->getFrameFilter(n, d->node, frameCtx);
		VSFrameRef *dst;
        // The reason we query this on a per frame basis is because we want our filter
        // to accept clips with varying dimensions. If we reject such content using d->vi
        // would be better.
        const VSFormat *fi = d->vi->format;
		int height = vsapi->getFrameHeight(src, 0);
        int width = vsapi->getFrameWidth(src, 0);		

		// When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
        // supply the "dominant" source frame to copy properties from. Frame props
        // are an essential part of the filter chain and you should NEVER break it.       

		int nbits = (fi->sampleType == stInteger) ? fi->bitsPerSample : 0;
		int	nbytes = fi->bytesPerSample;

		if (d->test)
		{
			// copy input frame on to dst
			dst = vsapi->copyFrame(src, core);
			// we will draw a distortion grid on to the frame

			for (int plane = 0; plane < fi->numPlanes; plane++)
			{
				const uint8_t* sp = vsapi->getReadPtr(src, plane);
				uint8_t* dp = vsapi->getWritePtr(dst, plane);
				int spitch = vsapi->getStride(src, plane) / nbytes; // this is in pixels. Required for uint16 and float pointers
				// note that if a frame has the same dimensions and format, the stride is guaranteed to be the same. int dst_stride = src_stride would be fine too in this filter.
				int dpitch = spitch;

				int iCenter = height / 2 * dpitch + width / 2;	// offset to center	
				if (fi->colorFamily == cmRGB)
				{
					if (nbytes == 1)
						dimplaneRGB(dp, sp, spitch, width, height, d->dim);
					else if (nbytes == 2)
						dimplaneRGB((uint16_t*)dp, (uint16_t*)sp, spitch, width, height, d->dim);
					else if (nbytes == 4)
						dimplaneRGB((float*)dp, (float*)sp, spitch, width, height, d->dim);
				}

				else if (plane == 0 && fi->colorFamily == cmYUV)
				{
					if (nbytes == 1)
					{
						uint8_t limit = (uint8_t)16;
						dimplaneYUV(dp, dp, dpitch, width, height, d->dim, limit);
					}
					else if (nbytes == 2)
					{
						uint16_t limit = (uint16_t)(16 << (nbits - 8));
						dimplaneYUV((uint16_t*)dp, (uint16_t*)dp, dpitch, width, height, d->dim, limit);
					}
					else if (nbytes == 4)
						dimplaneYUV((float*)dp, (float*)dp, dpitch, width, height, d->dim, 0.0f);
				}
				int p = plane;
				// we will mark every  pixelperdot
				for (int h = d->ddensity / 2; h < height / 2; h += d->ddensity)
				{
					int hoff = h * width / 2 * d->nEntries;

					for (int w = d->ddensity / 2; w < width / 2; w += d->ddensity)
					{
						int woff = w * d->nEntries;
						// for each point get the coordinates of distorted image
						int x = d->xyAndQ[hoff + woff];
						if (x < 0)
							continue;
						int y = d->xyAndQ[hoff + woff + 1];

						// paint dots 4 fold symmetry.
						if (nbytes == 1)
							paint4FoldSym(dp + iCenter, dpitch, 1, x, y, d->col[p]);
						else if (nbytes == 2)
							paint4FoldSym((uint16_t*)dp + iCenter, dpitch, 1, x, y, *((uint16_t*)d->col + p));
						else if (nbytes == 4)
							paint4FoldSym((float*)dp + iCenter, dpitch, 1, x, y, *((float*)d->col + p));
					}
				}
			}
		}						

		else	// not test. normal processing
		{
			dst = vsapi->newVideoFrame(fi, width, height, src, core);

			for (int plane = 0; plane < fi->numPlanes; plane++)
			{
				const uint8_t *sp = vsapi->getReadPtr(src, plane);
				int spitch = vsapi->getStride(src, plane) / nbytes;
				uint8_t *dp = vsapi->getWritePtr(dst, plane);
				int dpitch = vsapi->getStride(dst, plane) / nbytes;

				int iCenter = (height / 2) * spitch + width / 2;	// offset from left top to center of frame
				int oCenter = iCenter;
				uint8_t min8 = 0, max8 = (uint8_t)255;
				uint16_t min16 = (uint16_t)(fi->colorFamily == cmYUV ? 16 << (nbits - 8) : 0);
				uint16_t max16 = (uint16_t)((fi->colorFamily == cmYUV ? 235 : 255 << (nbits - 8)) << (nbits - 8));
				float minf = 0, maxf = 1.0f;

				if (plane > 0 && fi->colorFamily == cmYUV)
				{
					minf = -0.5f;
					maxf = 0.5f;
				}
				int x, y, qx, qy, span2 = d->span / 2;
				int hoff, woff;
				int index;
				int p = plane;

				for (int h = 0; h < height / 2; h++)
				{
					hoff = h * d->nEntries * width / 2;

					for (int w = 0; w < width / 2; w++)
					{
						woff = w * d->nEntries;
						// for each point get the coordinates of distorted image

						x = d->xyAndQ[hoff + woff];

						if (x < 0)
						{
							
							// points are outside  src or fish. So make the point black
							if (nbytes == 1)
								paint4FoldSym(dp + oCenter, dpitch, 1, w, h, d->col[p]);
							else if (nbytes == 2)
								paint4FoldSym((uint16_t*)dp + oCenter, dpitch, 1, w, h, *((uint16_t*)d->col + p));
							else if (nbytes == 4)
								paint4FoldSym((float*)dp + oCenter, dpitch, 1, w, h, *((float*)d->col + p));
						}

						else
						{
							
							y = d->xyAndQ[hoff + woff + 1];

							if (d->q > 1)
							{
								qx = d->xyAndQ[hoff + woff + 2];
								qy = d->xyAndQ[hoff + woff + 3];
							}
							else
								index = d->xyAndQ[hoff + woff + 2];

							if (x >= width / 2 - span2 - 1 || y >= height / 2 - span2 - 1)
							{
								
								//  interpolation does not have sufficient points
								// points are within src frame
								if (nbytes == 1)
								{
									// near point
									copy4FoldSym(dp + oCenter, dpitch, sp + iCenter, spitch, 1, w, h, x, y);
								}
								else if (nbytes == 2)
								{
									copy4FoldSym((uint16_t*)dp + oCenter, dpitch, (uint16_t*)sp + iCenter, spitch, 1, w, h, x, y);
								}
								else if (nbytes == 4)
								{

									// near point
									copy4FoldSym((float*)dp + oCenter, dpitch, (float*)sp + iCenter, spitch, 1, w, h, x, y);
								}
							}

							else
							{
								
								// sufficient points for interpolation are available
								if (nbytes == 1)
								{
									if (d->q == 1)
									{
										//continue;
										interpolate9pt4FoldSym(dp + oCenter, dpitch, sp + iCenter, spitch,
											1, w, h, x, y, index);
									}
									else
									{
										
										//bilinear 2x2 or cubic 4x4 or lanczos 6x6
										interpolate4FoldSym(dp + oCenter, dpitch, sp + iCenter, spitch,
											1, w, h, x, y, qx, qy, d->span, d->iCoeff, min8, max8);
									}

								}

								else if (nbytes == 2)
								{
									if (d->q == 1)
										interpolate9pt4FoldSym((uint16_t*)dp + oCenter, dpitch, (uint16_t*)sp + iCenter, spitch,
											1, w, h, x, y, index);
									else

										interpolate4FoldSym((uint16_t*)dp + oCenter, dpitch, (uint16_t*)sp + iCenter, spitch,
											1, w, h, x, y, qx, qy, d->span, d->iCoeff, min16, max16);

								}

								else if (nbytes == 4)
								{
									if (d->q == 1)
										interpolate9pt4FoldSym((float*)dp + oCenter, dpitch, (float*)sp + iCenter, spitch,
											1, w, h, x, y, index);
									else

										interpolate4FoldSym((float*)dp + oCenter, dpitch, (float*)sp + iCenter, spitch,
											1, w, h, x, y, qx, qy, d->span, d->iCoeff, minf, maxf);

								}
							}
						}
					}
				}

			}	// for plane
		}	// not test
		
			// activation reason arAllFramesReady processing completed
		vsapi->freeFrame (src);
		return dst;		
	}	// activation reason
	return 0;
}
//---------------------------------------------------------------------
static void VS_CC debarrelFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	DeBarrelData* d = (DeBarrelData*)instanceData;
	vsapi->freeNode(d->node);

	vs_aligned_free(d->xyAndQ);
	if (!d->iCoeff == NULL)
		vs_aligned_free(d->iCoeff);

	free(d);
}


/***************************************************************/
// This is the function that created the filter, when the filter has been called.
// This can be used for simple parameter checking, so it is possible to create different filters,

static void VS_CC debarrelCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
    DeBarrelData d;
    DeBarrelData *data;
    int err;
	int temp;	// used to convert int to bool
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    // In this first version we only want to handle 8bit integer formats. Note that
    // vi->format can be 0 if the input clip can change format midstream.
    if (!isConstantFormat(d.vi) || d.vi->width == 0 || d.vi->height == 0 || d.vi->format->subSamplingH != 0 || d.vi->format->subSamplingW != 0)
	{
        vsapi->setError(out, "DeBarrel: only RGB or those YUV formats that have no subsampling are supported. Frame dimensions should remain constant");
        vsapi->freeNode(d.node);
        return;
    }

	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "DeBarrel: RGB, YUV and Gray color formats only for input allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "DeBarrel: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}

	d.method = int64ToIntS(vsapi->propGetInt(in, "method", 0, &err));
	if (err)
		d.method = 2;
	if (d.method < 1 || d.method > 2)
	{
		vsapi->setError(out, "DeBarrel: method can be 1 or 2 ");
		vsapi->freeNode(d.node);
		return;
	}

	if (d.method == 1)
	{
		for (int i = 0; i < 3; i++)
		{
			d.abc[i] = (float)vsapi->propGetFloat(in, "abc", i, 0);
			if (d.abc[i] < 0.0 || d.abc[i] > 0.5)
			{
				vsapi->setError(out, "DeBarrel: abc[] values can be zero to less than 0.5 only");
				vsapi->freeNode(d.node);
				return;
			}
		}

		if (d.abc[0] + d.abc[1] + d.abc[2] > 1.0)
		{
			vsapi->setError(out, "DeBarrel: sum of all three abc array values must be less than 1.0 ");
			vsapi->freeNode(d.node);
			return;
		}
	}

	else
	{
		for (int i = 0; i < 3; i++)
		{
			d.abc[i] = (float)vsapi->propGetFloat(in, "abc", i, 0);
		}

		if (d.abc[2] < 0.0 || d.abc[2] >= 0.5f)
		{
			vsapi->setError(out, "DeBarrel: third value of abc[]  can be zero to less than 0.5");
			vsapi->freeNode(d.node);
			return;
		}
	}

	temp = !!int64ToIntS(vsapi->propGetInt(in, "pin", 0, &err));

	if (err)
		temp = 0;

	if (temp == 0)
		d.pin = false;
	else
		d.pin = true;
	d.q = int64ToIntS(vsapi->propGetInt(in, "q", 0, &err));
	if (err)
		d.q = 1;
	else if(d.q < 0 || d.q > 4)
	{
		vsapi->setError(out, "DeBarrel: q  can be 1 to 4 only");
		vsapi->freeNode(d.node);
		return;
	}
	
	// d.test
	temp = !!int64ToIntS(vsapi->propGetInt(in, "test", 0, &err));

	if (err)
		temp = 0;

	if (temp == 0)
		d.test = false;
	else
		d.test = true;

	if (d.test)
	{
		d.dots = int64ToIntS(vsapi->propGetInt(in, "dots", 0, &err));

		if (err)
			d.dots = 2;
		if (d.dots < 0 || d.dots > 4)			
		{
			vsapi->setError(out, "DeBarrel: dots can be 1 to 4 only.");
			vsapi->freeNode(d.node);
			return;
		}

		d.dim = (float)( 1.0 - vsapi->propGetFloat(in, "dots", 0, &err));
		if (err)
			d.dim = 0.75;
		else if (d.dim < 0 || d.dim > 1.0f)
		{
			vsapi->setError(out, "DeBarrel: dim can be 0 to 1.0 only.");
			vsapi->freeNode(d.node);
			return;
		}

	}	// if test


    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string, the only
    // reason this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
 
	
	data = (DeBarrelData *) malloc(sizeof(d));
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
    vsapi->createFilter(in, out, "DeBarrel", debarrelInit, debarrelGetFrame, debarrelFree, fmParallel, 0, data, core);
}
//////////////////////////////////////////
// Init

// This is the entry point that is called when a plugin is loaded. You are only supposed
// to call the two provided functions here.
// configFunc sets the id, namespace, and long name of the plugin (the last 3 arguments
// never need to be changed for a normal plugin).
//
// id: Needs to be a "reverse" url and unique among all plugins.
// It is inspired by how android packages identify themselves.
// If you don't own a domain then make one up that's related
// to the plugin name.
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
// or not empty arrays are accepted
/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin) {
    configFunc("in.vcmohan.repo", "geom", "VapourSynth DeBarrel plugin", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("DeBarrel", "clip:clip;abc:float[];method:int:opt;pin:int:opt;q:int:opt;"
	"test:int:opt;dots:data:opt; dim:float:opt;", debarrelCreate, 0, plugin);
}

*/


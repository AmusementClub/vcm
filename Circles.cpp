/******************************************************************************
Circles filter plugin for Vapoursynth+ 32/64 bit version
Circles draws concentric circles with given origin coordinates and marks in bold specified diameter. 
 Thread agnostic (operates under multi thread mode)
 Author V.C.Mohan

 10 Sep 2021
Copyright (C) <2021>  <V.C.Mohan>

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

//#include "windows.h"
//#include <stdint.h>const finc* sp,
// #define _USE_MATH_DEFINES
//#include "math.h"
//#include "InterpolationPack.h"
//
//#include "VapourSynth.h"
//

//-------------------------------------------------------------------------
typedef struct {
	VSNodeRef* node;
	const VSVideoInfo* vi;

	// circles origin coordinates
	int origin_y;
	int origin_x;	
	int dots;	// dot density 1 for 16, 2 for 12, 3 for 8, 4 for 4.	
	float dim;	// imaged dimmed to this level
	int fdia;	// circle with this diameter will be high lighted
	int cint;	// interval between circles
	int ddensity;	// pixel interval of dots of
	unsigned char col[16];		// color components for fill
	unsigned char bgr[3];
} CirclesData;




/*--------------------------------------------------
 * The following is the implementation
 * of the defined functions.
 --------------------------------------------------*/
 //Here is the acutal constructor code used

static void VS_CC circlesInit(VSMap* in, VSMap* out, void** instanceData, VSNode* node, VSCore* core, const VSAPI* vsapi)
{
	CirclesData* d = (CirclesData*)*instanceData;
	vsapi->setVideoInfo(d->vi, 1, node);

	int ht = d->vi->height;
	int wd = d->vi->width;
	const VSFormat* fi = d->vi->format;
	int nbytes = fi->bytesPerSample;
	int nbits = fi->bitsPerSample;
	
	d->ddensity = (5 - d->dots) * 16;
	//uint8_t bgr[] = { 0,0,0 };
	uint8_t  yuv[3];
	BGR8_YUV(yuv, d->bgr, 8);
	uint8_t* bgr = d->bgr;

	for (int k = 0; k < 3; k++)
	{
		if (nbytes == 1)
		{
			if (fi->colorFamily == cmRGB)
				d->col[k] = bgr[k];
			else
				d->col[k] = yuv[k];
		}
		else if (nbytes == 2)
		{
			if (fi->colorFamily == cmRGB)
				*((uint16_t*)d->col + k) = (uint16_t)(bgr[k] << (nbits - 8));
			else
				*((uint16_t*)d->col + k) = (uint16_t)(yuv[k] << (nbits - 8));
		}
		else // float
		{
			if (fi->colorFamily == cmRGB)
				*((float*)d->col + k) = (float)(bgr[k] / 255.0f);
			else
			{
				if (k == 0)
					*((float*)d->col + k) = (float)((yuv[k] - 16) / 220.0f);
				else
					*((float*)d->col + k) = (float)((yuv[k] - 128) / 220.0f);
			}
		}

	}
}

//------------------------------------------------------------------------------------------------

static const VSFrameRef* VS_CC circlesGetFrame(int n, int activationReason, void** instanceData,
	void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	CirclesData* d = (CirclesData*)*instanceData;

	if (activationReason == arInitial)
	{
		vsapi->requestFrameFilter(n, d->node, frameCtx);
	}
	else if (activationReason == arAllFramesReady)
	{
		const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);
		VSFrameRef* dst;
		const VSFormat* fi = d->vi->format;
		int height = vsapi->getFrameHeight(src, 0);
		int width = vsapi->getFrameWidth(src, 0);
		int nbits = fi->bitsPerSample;
		int nbytes = fi->bytesPerSample;
		//will not process A plane
		int np = fi->numPlanes > 3 ? 3 : fi->numPlanes;
		
		int kb = 1;

			// get src on which dots will be overlain			
		dst = vsapi->copyFrame(src, core);
		
		int frsq = d->fdia * d->fdia / 4;

		for (int p = 0; p < np; p++)
		{

			const uint8_t* sp = vsapi->getReadPtr(src, p);
			uint8_t* dp = vsapi->getWritePtr(dst, p);
			int pitch = vsapi->getStride(src, p) / nbytes;
			
			if ( fi->colorFamily == cmRGB)
			{
					if (nbytes == 1)
						dimplaneRGB(dp,sp, pitch, width, height, d->dim);
					else if (nbytes == 2)
						dimplaneRGB((uint16_t*)dp, (uint16_t*)sp, pitch, width, height, d->dim);
					else if (nbytes == 4)
						dimplaneRGB((float*)dp, (float *)sp, pitch, width, height, d->dim);
			}
			else if (p == 0 && fi->colorFamily == cmYUV)
			{
				if (nbytes == 1)
				{
					uint8_t limit = 16;
					dimplaneYUV(dp, sp, pitch, width, height, d->dim, limit);
				}
				else if (nbytes == 2)
				{
					uint16_t limit = (uint16_t)(16 << (nbits - 8));
					dimplaneYUV((uint16_t*)dp, (uint16_t*)sp, pitch, width, height, d->dim, limit);
				}
				else if (nbytes == 4)
				{
					float limit = 0.0f;
					dimplaneYUV((float*)dp, (float*)sp, pitch, width, height, d->dim, limit);
				}
			}
			int fdmax = VSMAX(VSMAX(abs(d->origin_y - height), abs(d->origin_y)), VSMAX(abs(d->origin_x - width), abs(d->origin_x)))  ;
			// draw concentric circles
			for (int fd = 0; fd <= fdmax; fd += d->cint)
			{
				int inc = fd == 0 ? 1 : 6 - d->dots;
				int frad = fd == 0 ? d->fdia / 2 : fd;

				for (int alfa = 1; alfa < 90; alfa += inc)
				{
					double cosalfa = cos(alfa * M_PI / 180);
					double sinalfa = sin(alfa * M_PI / 180);

					int x = (int)(frad * cosalfa);
					int y = (int)(frad * sinalfa);

					for (int i = -1; i < 2; i += 2)
					{
						for (int j = -1; j < 2; j += 2)
						{
							int w = d->origin_x + x * j;
							int h = d->origin_y + y * i;
							if (w >= 0 && w < width && h >= 0 && h < height)
							{
								if (nbytes == 1)
									*(dp + h * pitch + w) = d->col[p];

								else if (nbytes == 2)
									*((uint16_t*)dp + h * pitch + w) = *((uint16_t*)d->col + p);
								else if (nbytes == 4)
									*((float*)dp + h * pitch + w) = *((float*)d->col + p);
							}
						}
					}
				}
			}			
				
		}
		vsapi->freeFrame(src);
		return dst;
	}
	return 0;
}

/***************************************************************/
static void VS_CC circlesFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	CirclesData* d = (CirclesData*)instanceData;
	vsapi->freeNode(d->node);
	
	free(d);
}

static void VS_CC circlesCreate(const VSMap* in, VSMap* out, void* userData,
	VSCore* core, const VSAPI* vsapi)
{
	CirclesData d;
	CirclesData* data;
	int err;
	int temp;

	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);

	// In this first version we only want to handle 8bit integer formats. Note that
	// vi->format can be 0 if the input clip can change format midstream.
	if (!isConstantFormat(d.vi) || d.vi->width == 0 || d.vi->height == 0
		|| (d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray
			&& d.vi->format->colorFamily != cmRGB))
	{
		vsapi->setError(out, "Circles: only RGB, Yuv or Gray color constant formats and const frame dimensions input supported");
		vsapi->freeNode(d.node);
		return;
	}
	

	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "Circles: half float input not allowed.");
		vsapi->freeNode(d.node);
		return;
	}

	// If a property read fails for some reason (index out of bounds/wrong type)
	// then err will have flags set to indicate why and 0 will be returned. This
	// can be very useful to know when having optional arguments. Since we have
	// strict checking because of what we wrote in the argument string, the only
	// reason this could fail is when the value wasn't set by the user.
	// And when it's not set we want it to default to enabled.	
	
	d.origin_x = int64ToIntS(vsapi->propGetInt(in, "xo", 0, &err));
	if (err)
		d.origin_x = d.vi->width / 2;

	d.origin_y = int64ToIntS(vsapi->propGetInt(in, "yo", 0, &err));
	if (err)
		d.origin_y = d.vi->height / 2;

	d.fdia = int64ToIntS(vsapi->propGetInt(in, "frad", 0, &err)) * 2;
	if (err)
		d.fdia = (d.vi->height > d.vi->width ? d.vi->width : d.vi->height) * 2;

	else if (d.fdia < 128 )
	{
		vsapi->setError(out, "Circles: frad must be at least 64  ");
		vsapi->freeNode(d.node);
		return;
	}
	if ( (d.origin_x + d.fdia < 0 || d.origin_x - d.fdia > d.vi->width) 
		|| (d.origin_y + d.fdia < 0 || d.origin_y - d.fdia > d.vi->height))
	{
		vsapi->setError(out, "Circles: fdia and origin takes the fisheye image outside frame  ");
			vsapi->freeNode(d.node);
			return;
	}
	d.cint = int64ToIntS(vsapi->propGetInt(in, "cint", 0, &err));
	if (err)
		d.cint = 50;

	else if (d.cint < 10 || d.cint > VSMAX(d.vi->width, d.vi->height) / 2 )
	{
		vsapi->setError(out, "Circles: cint must be at least 10 and not more than half of max dimension of frame");
		vsapi->freeNode(d.node);
		return;
	}
	
	d.dots = int64ToIntS(vsapi->propGetInt(in, "dots", 0, &err));
	if (err)
		d.dots = 2;
	else if (d.dots < 1 || d.dots > 4)
	{
		vsapi->setError(out, "Circles: dots must be 1 to 4 only ");
		vsapi->freeNode(d.node);
		return;
	}

	d.dim = (float)(1.0 - vsapi->propGetFloat(in, "dim", 0, &err));
	if (err)
		d.dim = 0.25f;
	if (d.dim < 0.0f || d.dim > 1.0f)
	{
		vsapi->setError(out, "Circles: dim must be from 0 to 1.0 only ");
		vsapi->freeNode(d.node);
		return;
	}
	for (int i = 0; i < 3; i++)
	{
		temp = int64ToIntS(vsapi->propGetInt(in, "rgb", 0, &err));
		if (err)
			d.bgr[2 - i] = (uint8_t)255;
		else if (temp < 0 || temp > 255)
		{
			vsapi->setError(out, "Circles: rgb array values must be from 0 to 255 only ");
			vsapi->freeNode(d.node);
			return;
		}
		else d.bgr[2 - i] = (uint8_t)temp;
	}
	
	// I usually keep the filter data struct on the stack and don't allocate it
// until all the input validation is done.
	data = (CirclesData*)malloc(sizeof(d));
	*data = d;

	vsapi->createFilter(in, out, "Circles", circlesInit, circlesGetFrame, circlesFree, fmParallel, 0, data, core);

}

// registerFunc("Circles", "clip:clip;xo:int:opt;yo:int:opt;fdia:int:opt;cint:int:opt;dots:int:opt;rgb:int[]:opt;dim:float:opt;", circlesCreate, 0, plugin);



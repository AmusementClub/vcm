/******************************************************************************
Fisheye filter plugin for Vapoursynth+ 32/64 bit version
Fisheye corrects input image recorded as  circular (fishEye) projection
to a rectangular image
Fish eye correction as in the paper "A Flexible Architecture for Fisheye
Correction in Automotive Rear-View Cameras" of ALTERA Manipal DOT NET of 2008
Material from Wikipedia used for barrel and pincushion corrections
 Thread agnostic (operates under multi thread mode)
 Author V.C.Mohan

 15 oct 2021
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
	const VSVideoInfo* ivi;
	VSVideoInfo vi;
	
	// starting or constant for x and y or for x alone coefficients
	int origin_y;
	int origin_x;
	bool test;	// whether a test run
	
	int dots;	// dot density 1 for 16, 2 for 12, 3 for 8, 4 for 4.
	int method;
	float dim;
	double fov;
	double rix;
	int frad, fdia;
	int oRadius;
	int q;		// type of interpolation. 0 near pt, 1 ninept 2x2, 2 bilenear 2x2, 3 cubic 3x3, 4 lanczos 6x6
	bool sqr;	// is squaring circle required?
	float* iCoeff;
	int ddensity;	// pixel interval of dots of 
	int quantile;	// Accuracy of fraction to which values are interpolated
	int span;		// interpolation function 1d span or taps
	unsigned char col[16];		// color components for fill
	int* xyAndQ;
	double rNorm;
	int nEntries;
} FisheyeData;


//void getSquircleUV(float* cUV, int sqW, int sqH, int rad);
//#include "CorrectLD.cpp"
#include "FisheyePart.cpp"

/*--------------------------------------------------
 * The following is the implementation
 * of the defined functions.
 --------------------------------------------------*/
 //Here is the acutal constructor code used

static void VS_CC fisheyeInit(VSMap* in, VSMap* out, void** instanceData, VSNode* node, VSCore* core, const VSAPI* vsapi)
{
	FisheyeData* d = (FisheyeData*)*instanceData;	
	d->frad = d->fdia / 2;
	d->oRadius = d->frad;

	double focal = getFocalLength(d->frad, d->method, d->fov);

	d->oRadius = getOutputRadius(d->frad, focal, d->rix);
	
	if (!d->test)
	{
		d->vi.format = d->ivi->format;
		d->vi.numFrames = d->ivi->numFrames;
		d->vi.fpsDen = d->ivi->fpsDen;
		d->vi.fpsNum = d->ivi->fpsNum;
		d->vi.flags = d->ivi->flags;
		d->vi.height = d->oRadius * 2;
		d->vi.width = d->oRadius * 2;

		vsapi->setVideoInfo(&d->vi, 1, node);
	}
	else
	{
		// in test  frame dimensions remain unaltered
		vsapi->setVideoInfo(d->ivi, 1, node);
	}
	
	
	// input frame dimensions
	int swidth = d->ivi->width;
	int sheight = d->ivi->height;
	int frsq = d->frad * d->frad;
	// output
	int ht = d->oRadius * 2;
	int wd = d->oRadius * 2;
	const VSFormat* fi = d->vi.format;
	int nbytes = fi->bytesPerSample;
	int nbits = fi->bitsPerSample;
	d->quantile = 64;
	int nEntries = d->test ? 2 :d->q == 1? 3: 4;	
		
	d->xyAndQ = (int*)vs_aligned_malloc<int>(sizeof(int) * d->oRadius * d->oRadius * nEntries, 32);

	int* xyQ = d->xyAndQ;
	float xy[2];
	int x, y, qx, qy;	

	d->iCoeff = NULL;

	if ( ! d->test)
		d->iCoeff = setInterpolationScheme(d->q, d->quantile, &d->span);

	d->rNorm = (double)d->frad;
	
	//terms Barrel and Pincushion distortions are usually associated to longer focal lengths 
	// and fish eye to wide angle(shorter focal length) lenses
	
	float cUV[2]; // corresponding circle coordinates
	for (int h = 0; h < d->oRadius; h++)
	{

		for (int w = 0; w < d->oRadius; w++)
		{
			if (d->sqr )
			{
				getSquircleUV(cUV, w, h, d->oRadius);
			}
			else
			{
				cUV[0] = (float)w;
				cUV[1] = (float)h;
			}


			getSourceXY(xy, cUV[0], cUV[1], d->method, focal, d->rNorm, d->rix);

			x = (int)floor(xy[0]);
			y = (int)floor(xy[1]);

			int off = nEntries * (h * d->oRadius + w);

			if (x >= swidth / 2 || y >= sheight / 2 || x < 0 || y < 0
				|| x * x + y * y > frsq)
			{
				xyQ[off] = - 1;
			}
			else
			{
				// calculate nearest quantile of the fraction
				qx = (int)((xy[0] - x) * d->quantile);
				qy = (int)((xy[1] - y) * d->quantile);
				xyQ[off] = x;
				xyQ[off + 1] = y;

				if (!d->test)
				{
					if (d->q > 1)
					{
						xyQ[off + 2] = qx;
						xyQ[off + 3] = qy;
					}
					else
						// index value
						xyQ[off + 2] = bestOfNineIndex(qx, qy, d->quantile);
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
//------------------------------------------------------------------------------------------------

static const VSFrameRef* VS_CC fisheyeGetFrame(int n, int activationReason, void** instanceData,
	void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	FisheyeData* d = (FisheyeData*)*instanceData;

	if (activationReason == arInitial)
	{
		vsapi->requestFrameFilter(n, d->node, frameCtx);
	}
	else if (activationReason == arAllFramesReady)
	{
		const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);
		VSFrameRef* dst;
		const VSFormat* fi = d->ivi->format;
		int sheight = vsapi->getFrameHeight(src, 0);
		int swidth = vsapi->getFrameWidth(src, 0);
		int nbits = fi->bitsPerSample;
		int nbytes = fi->bytesPerSample;
		//will not process A plane
		int np = fi->numPlanes > 3 ? 3 : fi->numPlanes;
		int dwidth = d->vi.width;
		int dheight = d->vi.height;
		
		int kb = 1;

		if (d->test)
			// get src on which dots will be overlain			
			dst = vsapi->copyFrame(src, core);
		else
			dst = vsapi->newVideoFrame(fi, dwidth, dheight, src, core);

		dwidth = vsapi->getFrameWidth(dst, 0);
		dheight = vsapi->getFrameHeight(dst, 0);

		int frsq = d->frad * d->frad;

		for (int p = 0; p < np; p++)
		{

			const uint8_t* sp = vsapi->getReadPtr(src, p);
			uint8_t* dp = vsapi->getWritePtr(dst, p);
			int spitch = vsapi->getStride(src, p) / nbytes;
			int dpitch = vsapi->getStride(dst, p) / nbytes;
			// number of entries per row in the xyAndQ buffer
			int nEntries = d->test ? 2 : d->q == 1 ? 3 : 4;

			if (d->test)
			{
				
				if (fi->colorFamily == cmRGB)
				{
					if (nbytes == 1)
						dimplaneRGB(dp, sp, spitch, swidth, sheight, d->dim);
					else if (nbytes == 2)
						dimplaneRGB((uint16_t*)dp, (uint16_t*)sp, spitch, swidth, sheight, d->dim);
					else if (nbytes == 4)
						dimplaneRGB((float*)dp, (float*)sp, spitch, swidth, sheight, d->dim);
				}

				else if ( p == 0 && fi->colorFamily == cmYUV)
				{
					if (nbytes == 1)
					{
						uint8_t limit = (uint8_t)16;
						dimplaneYUV(dp, dp, dpitch, dwidth, dheight, d->dim, limit);
					}
					else if (nbytes == 2)
					{
						uint16_t limit = (uint16_t)(16 << (nbits - 8));
						dimplaneYUV((uint16_t*)dp, (uint16_t*)dp, dpitch, dwidth, dheight, d->dim, limit);
					}
					else if (nbytes == 4)
						dimplaneYUV((float*)dp, (float*)dp, dpitch, dwidth, dheight, d->dim, 0.0f);
				}

				
				int frsq = d->frad * d->frad;
				int iCenter = d->origin_y * spitch + d->origin_x;
				
					// we will put dots
				for (int h = d->ddensity / 2; h < d->oRadius; h += d->ddensity)
				{
					int hoff = nEntries * h * d->oRadius;

					for (int w = d->ddensity / 2; w < d->oRadius; w += d->ddensity)
					{
						int woff = nEntries * w;

						int x = d->xyAndQ[hoff + woff];
						int y = d->xyAndQ[hoff + woff + 1];
						// ensure points are within frame
						if (x >= d->frad || y >= d->frad || x < 0 || y < 0)
							continue;
						// ensure points are within src / fish eye						
						if ( x * x + y * y <= frsq)
						{
							if (nbytes == 1)
								paint4FoldSym(dp + iCenter, dpitch,1, x, y, d->col[p]);
							else if (nbytes == 2)
								paint4FoldSym( (uint16_t*)dp + iCenter, dpitch,1, x, y, *((uint16_t*)d->col + p) );
							else if (nbytes == 4)
								paint4FoldSym((float*)dp + iCenter, dpitch,1, x, y, *((float*)d->col + p) );
						}
						
					}

				}
				
			}	// if test

			else	// not test. normal processing
			{	
				int iCenter = d->origin_y * spitch + d->origin_x;
				int oCenter = d->oRadius * dpitch + d->oRadius;
				
				uint8_t min8 = 0, max8 = (uint8_t)255;
				uint16_t min16 = (uint16_t)(fi->colorFamily == cmYUV ? 16 << (nbits - 8) : 0);
				uint16_t max16 = (uint16_t)((fi->colorFamily == cmYUV ? 235 : 255 << (nbits - 8)) << (nbits - 8));
				float minf = 0, maxf = 1.0f;

				if (p > 0 && fi->colorFamily == cmYUV)
				{
					minf = -0.5f;
					maxf = 0.5f;
				}
				int x, y, qx, qy,  span2 = d->span / 2;
				int offh, offw;
				int index;
				int frsq = d->frad * d->frad;

				for (int h = 0; h < d->oRadius - 1; h++)
				{
					offh = h * d->oRadius * nEntries;

					for (int w = 0; w < d->oRadius - 1; w++)
					{				
						offw = nEntries * w;

						x = d->xyAndQ[offh + offw];
						y = d->xyAndQ[offh + offw + 1];

						if (d->q > 1)
						{
							qx = d->xyAndQ[offh + offw + 2];
							qy = d->xyAndQ[offh + offw + 3];
						}
						else
							index = d->xyAndQ[offh + offw + 2];
						// was checked in init
						//if (x >= swidth / 2   || y >= sheight / 2  || x < 0 || y < 0
						//	|| x * x + y * y > frsq)
						if ( x < 0)
						{
							// points are outside  src or fish. So make the point black
							if (nbytes == 1)
								paint4FoldSym(dp + oCenter, dpitch,1, w, h, d->col[p]);
							else if (nbytes == 2)
								paint4FoldSym((uint16_t*)dp + oCenter, dpitch, 1, w, h, *((uint16_t*)d->col + p));
							else if (nbytes == 4)
								paint4FoldSym((float*)dp + oCenter, dpitch, 1, w, h, *((float*)d->col + p));
						}

						else if (x >= swidth / 2 - span2 - 1 || y >= sheight / 2 - span2 -1 )
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
									interpolate9pt4FoldSym(dp + oCenter, dpitch, sp + iCenter, spitch,
										1, w, h, x, y,	index);
								else
									//bilinear 2x2 or cubic 4x4 or lanczos 6x6
									interpolate4FoldSym(dp + oCenter, dpitch, sp + iCenter, spitch,
										 1, w, h, x, y, qx, qy,  d->span, d->iCoeff, min8, max8);
																	
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
									1,w, h, x, y, qx, qy, d->span, d->iCoeff, minf, maxf);
								
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
static void VS_CC fisheyeFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	FisheyeData* d = (FisheyeData*)instanceData;
	vsapi->freeNode(d->node);
	
		vs_aligned_free(d->xyAndQ);
		if (!d->iCoeff == NULL)
			vs_aligned_free(d->iCoeff);
	
	free(d);
}

static void VS_CC fisheyeCreate(const VSMap* in, VSMap* out, void* userData,
	VSCore* core, const VSAPI* vsapi)
{
	FisheyeData d;
	FisheyeData* data;
	int err;
	int temp;

	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.ivi = vsapi->getVideoInfo(d.node);

	// In this first version we only want to handle 8bit integer formats. Note that
	// vi->format can be 0 if the input clip can change format midstream.
	if (!isConstantFormat(d.ivi) || d.ivi->width == 0 || d.ivi->height == 0
		|| (d.ivi->format->colorFamily != cmYUV && d.ivi->format->colorFamily != cmGray
			&& d.ivi->format->colorFamily != cmRGB))
	{
		vsapi->setError(out, "Fisheye: only RGB, Yuv or Gray color constant formats and const frame dimensions input supported");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.ivi->format->colorFamily == cmYUV && (d.ivi->format->subSamplingH != 0 || d.ivi->format->subSamplingW != 0))
	{
		vsapi->setError(out, "Fisheye: for YUV input only YUV444 allowed");
		vsapi->freeNode(d.node);
		return;
	}

	if (d.ivi->format->sampleType == stFloat && d.ivi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "Fisheye: half float input not allowed.");
		vsapi->freeNode(d.node);
		return;
	}

	// If a property read fails for some reason (index out of bounds/wrong type)
	// then err will have flags set to indicate why and 0 will be returned. This
	// can be very useful to know when having optional arguments. Since we have
	// strict checking because of what we wrote in the argument string, the only
	// reason this could fail is when the value wasn't set by the user.
	// And when it's not set we want it to default to enabled.

	d.method = int64ToIntS(vsapi->propGetInt(in, "method", 0, &err));
	if (err)
		d.method = 3;
	if (d.method < 1 || d.method > 5)
	{
		vsapi->setError(out, "Fisheye: method must be between 1 and 5 ");
		vsapi->freeNode(d.node);
		return;
	}
	
	d.origin_x = int64ToIntS(vsapi->propGetInt(in, "xo", 0, &err));
	if (err)
		d.origin_x = d.ivi->width / 2;

	d.origin_y = int64ToIntS(vsapi->propGetInt(in, "yo", 0, &err));
	if (err)
		d.origin_y = d.ivi->height / 2;

	int radius = VSMAX( VSMAX(abs(d.origin_x), abs(d.ivi->width - d.origin_x) ),
		VSMAX(abs(d.origin_y), abs(d.ivi->height - d.origin_y) ) );
	


	d.frad = int64ToIntS(vsapi->propGetInt(in, "frad", 0, &err));	
	if (err)
		d.frad = radius; // (d.ivi->width > d.ivi->height ? d.ivi->height : d.ivi->width) / 2;

	else if (d.frad < 64)
	{
		vsapi->setError(out, "Fisheye: frad must be at least 64 ");
		vsapi->freeNode(d.node);
		return;
	}

	d.fov = (double)vsapi->propGetFloat(in, "fov", 0, &err);

	if (err)
		d.fov = 120.0;
	else if (d.fov < 40 || d.fov > 170)
	{
		vsapi->setError(out, "Fisheye: fov can be 40 to 170 only ");
		vsapi->freeNode(d.node);
		return;
	}

	bool insideFrame = false, partiallyInside = false;
	if (d.origin_x - d.frad >= 0 && d.origin_x + d.frad < d.ivi->width && d.origin_y - d.frad >= 0 && d.origin_y + d.frad < d.ivi->height)
	
		insideFrame = true;
	else if (((d.origin_x < 0 && d.origin_x + d.frad > 0) || (d.origin_x > d.ivi->width && d.origin_x - d.frad < d.ivi->width)
		|| (d.origin_x >= 0 && d.origin_x < d.ivi->width ))
		&& ((d.origin_y < 0 && d.origin_y + d.frad > 0) || (d.origin_y > d.ivi->height && d.origin_y - d.frad < d.ivi->height)
			|| (d.origin_y >= 0 && d.origin_y < d.ivi->height )))
	
		partiallyInside = true;
	else
	{
		vsapi->setError(out, "Fisheye: origin and frad must ensure at least part of fisheye image is inside frame ");
		vsapi->freeNode(d.node);
		return;
	}
	d.fdia = d.frad * 2;

	
	temp = !!int64ToIntS(vsapi->propGetInt(in, "sqr", 0, &err));
	if (err)
		d.sqr = true;
	else
		d.sqr = temp == 0 ? false : true;
	

	d.rix = (double)(vsapi->propGetFloat(in, "rix", 0, &err));
	if (err)
		d.rix = 1.15;
	if (d.rix < 1.0 || d.rix > 1.5)
	{
		vsapi->setError(out, "Fisheye: rix must be 1.0 to 1.5 ");
		vsapi->freeNode(d.node);
		return;
	}


	temp = !!int64ToIntS(vsapi->propGetInt(in, "test", 0, &err));
	if (err)
		d.test = false;
	else
		d.test = temp == 0 ? false : true;
	if (d.test)
	{
		
		d.dots = int64ToIntS(vsapi->propGetInt(in, "dots", 0, &err));
		if (err)
			d.dots = 2;
		else if (d.dots < 1 || d.dots > 4)
		{
			vsapi->setError(out, "Fisheye: dots must be 1 to 4 only ");
			vsapi->freeNode(d.node);
			return;
		}

		d.dim = (float)(1.0 - vsapi->propGetFloat(in, "dim", 0, &err));
		if (err)
			d.dim = 0.75f;
		if (d.dim < 0.0f || d.dim > 1.0f)
		{
			vsapi->setError(out, "Fisheye: dim must be from 0 to 1.0 only ");
			vsapi->freeNode(d.node);
			return;
		}
		
	}
	else
	{
		
		d.q = int64ToIntS(vsapi->propGetInt(in, "q", 0, &err));
		if (err)
			d.q = 1;
		else if (d.q < 1 || d.q > 4)
		{
			vsapi->setError(out, "Fisheye: q must be 1 to 4 only ");
			vsapi->freeNode(d.node);
			return;
		}
	}
	

	// I usually keep the filter data struct on the stack and don't allocate it
// until all the input validation is done.
	data = (FisheyeData*)malloc(sizeof(d));
	*data = d;

	if (insideFrame)
		vsapi->createFilter(in, out, "Fisheye", fisheyeInit, fisheyeGetFrame, fisheyeFree, fmParallel, 0, data, core);
	else
		vsapi->createFilter(in, out, "FisheyePart", fisheyepartInit, fisheyepartGetFrame, fisheyepartFree, fmParallel, 0, data, core);
}

// registerFunc("Fisheye", "clip:clip;method:int:opt;xo:int:opt;yo:int:opt;frad:int:opt;sqr:int:opt;rix:float:opt;fov:float:opt;test:int:opt;q:int:opt;dots:int:opt;", fisheyeCreate, 0, plugin);



/******************************************************************************
Fisheye filter plugin for Vapoursynth+ 32/64 bit version
Fisheye corrects input image recorded as  circular (fishEye) projection
to a rectangular image
Fish eye correction as in the paper "A Flexible Architecture for Fisheye
Correction in Automotive Rear-View Cameras" of ALTERA Manipal DOT NET of 2008
Material from Wikipedia used for barrel and pincushion corrections
 Thread agnostic (operates under multi thread mode)
 Author V.C.Mohan

 14 Oct 2021
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


/*--------------------------------------------------
 * The following is the implementation
 * of the defined functions.
 --------------------------------------------------*/
 //Here is the acutal constructor code used

static void VS_CC fisheyepartInit(VSMap* in, VSMap* out, void** instanceData, VSNode* node, VSCore* core, const VSAPI* vsapi)
{
	FisheyeData* d = (FisheyeData*)*instanceData;	
	vsapi->setVideoInfo(d->ivi, 1, node);
	double focal = getFocalLength(d->frad, d->method, d->fov);

	d->oRadius = getOutputRadius(d->frad, focal, d->rix);
	// input frame dimensions
	int swidth = d->ivi->width;
	int sheight = d->ivi->height;
	int frsq = d->frad * d->frad;
	// output
	
	const VSFormat* fi = d->ivi->format;
	int nbytes = fi->bytesPerSample;
	int nbits = fi->bitsPerSample;
	d->quantile = 64;
	d->nEntries = d->test ? 2 :d->q == 1? 3: 4;	
		
	d->xyAndQ = (int*)vs_aligned_malloc<int>(sizeof(int) * swidth * sheight * d->nEntries, 32);
	d->rNorm = 1.0; // value not used in this part
	int* xyQ = d->xyAndQ;
	float xy[2];
	int x, y, qx, qy;	

	d->iCoeff = NULL;

	if ( ! d->test)
		d->iCoeff = setInterpolationScheme(d->q, d->quantile, &d->span);

	float rNorm = (float)d->frad;
	
	//terms Barrel and Pincushion distortions are usually associated to longer focal lengths 
	// and fish eye to wide angle(shorter focal length) lenses
	
	//float cUV[2]; // corresponding circle coordinates
	for (int h = 0; h < sheight; h++)
	{
		float hh = (float)(h - d->origin_y);

		for (int w = 0; w < swidth; w++)
		{
			float ww = (float)(w - d->origin_x);

			getSourceXY(xy, ww, hh, d->method, focal, rNorm, d->rix);

			x = (int)(floor(xy[0]));
			y = (int)(floor(xy[1]));

			// calculate nearest quantile of the fraction
			qx = (int)((xy[0] - x) * d->quantile);
			qy = (int)((xy[1] - y) * d->quantile);

			x += d->origin_x;
			y += d->origin_y;

			int off = d->nEntries * (h * swidth + w);


			if (x >= swidth  || y >= sheight  || x < 0 || y < 0
				|| (x - d->origin_x) * (x - d->origin_x) + (y - d->origin_y) * (y - d->origin_y) > frsq)
			{
				xyQ[off] = - 1;
			}
			else
			{
				
				xyQ[off] = x;
				xyQ[off + 1] = y;

				if (!d->test)
				{
					if (d->q > 1)
					{
						xyQ[off + 2] = qx;
						xyQ[off + 3] = qy;
					}
					else if (d->q == 1)
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

//...............................................................
static const VSFrameRef* VS_CC fisheyepartGetFrame(int n, int activationReason, void** instanceData,
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
				
		int kb = 1;

		if (d->test)
			// get src on which dots will be overlain			
			dst = vsapi->copyFrame(src, core);
		else
			dst = vsapi->newVideoFrame(fi, swidth, sheight, src, core);

		int dwidth = vsapi->getFrameWidth(dst, 0);
		int dheight = vsapi->getFrameHeight(dst, 0);

		int frsq = d->frad * d->frad;

		for (int p = 0; p < np; p++)
		{

			const uint8_t* sp = vsapi->getReadPtr(src, p);
			uint8_t* dp = vsapi->getWritePtr(dst, p);
			int spitch = vsapi->getStride(src, p) / nbytes;
			int dpitch = vsapi->getStride(dst, p) / nbytes;
			// number of entries per row in the xyAndQ buffer
			int nEntries = d->nEntries;

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
				
					// we will put dots
				for (int h = d->ddensity / 2; h < sheight; h += d->ddensity)
				{
					int hoff = nEntries * h * swidth;

					for (int w = d->ddensity / 2; w < swidth; w += d->ddensity)
					{
						int woff = nEntries * w;

						int x = d->xyAndQ[hoff + woff];
						int y = d->xyAndQ[hoff + woff + 1];

						// ensure points are within src / fish eye						
						if (x >= 0)
						{
							
							if (nbytes == 1)
								*(dp + y * dpitch + x) = d->col[p];
							else if (nbytes == 2)
								*( (uint16_t*)dp + y * dpitch + x) = *((uint16_t*)d->col + p);
							else if (nbytes == 4)
								*((float*)dp + y * dpitch + x) = *((float*)d->col + p);
						}						
					}
				}				
			}	// if test

			else	// not test. normal processing
			{	
				
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

				for (int h = 0; h < dheight - 1; h++)
				{
					offh = h * dwidth * nEntries;

					for (int w = 0; w < dwidth - 1; w++)
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
								*(dp + h * dpitch + w) = d->col[p];
							else if (nbytes == 2)
								*((uint16_t*)dp + h * dpitch + w) = *((uint16_t*)d->col + p);
							else if (nbytes == 4)
								*((float*)dp + h * dpitch + w) = *((float*)d->col + p);
						}

						else if (x >= swidth  - span2 - 1 || y >= sheight  - span2 -1
							|| x < span2 + 1 || y < span2 + 1)
						{
							//  interpolation does not have sufficient points
						// points are within src frame
							if (nbytes == 1)
							{
								// near point
								*(dp + h * dpitch + w) = *( sp + y * spitch + w);
							}
							else if (nbytes == 2)
							{
								*( (uint16_t*)dp + h * dpitch + w) = *((uint16_t*)sp + y * spitch + w);
							}
							else if (nbytes == 4)
							{
								
								// near point
								*((float*)dp + h * dpitch + w) = *((float*)sp + y * spitch + w);
							}
						}
						
						else if (x < swidth - span2 - 1 && y < sheight - span2 - 1 && x > span2 + 1 && y > span2 + 1)
						{
							
							// sufficient points for interpolation are available
							if (nbytes == 1)
							{
								if (d->q == 1)
									*(dp + h * dpitch + w) = bestOfNine(sp, spitch,1, x, y, index);
								else
									//bilinear 2x2 or cubic 4x4 or lanczos 6x6
									*(dp + h * dpitch + w) = clamp(LaQuantile(sp + y * spitch + x, spitch, 
										  qx, qy, d->span, d->iCoeff), min8, max8);
																	
							}

							else if (nbytes == 2)
							{
								if (d->q == 1)
									*((uint16_t*)dp + h * dpitch + w) = bestOfNine((uint16_t*)sp, spitch,1, x, y, index);
								else
								
									*((uint16_t*)dp + h * dpitch + w) = clamp(LaQuantile((uint16_t*)sp + y * spitch + x, spitch,
										qx, qy, d->span, d->iCoeff), min16, max16);
								
							}

							else if (nbytes == 4)
							{
								if (d->q == 1)
									*((float*)dp + h * dpitch + w) = bestOfNine((float*)sp, spitch,1, x, y, index);
								else
								
									*((float*)dp + h * dpitch + w) = clamp(LaQuantile((float*)sp + y * spitch + x, spitch,
										qx, qy, d->span, d->iCoeff), minf, maxf);
								
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
static void VS_CC fisheyepartFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	FisheyeData* d = (FisheyeData*)instanceData;
	vsapi->freeNode(d->node);
	
		vs_aligned_free(d->xyAndQ);
		if (!d->iCoeff == NULL)
			vs_aligned_free(d->iCoeff);
	
	free(d);
}



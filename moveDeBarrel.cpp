/******************************************************************************
DeBarrel filter plugin for vapoursynth by V.C.Mohan
DeBarrel corrects input image from barrel or pin cushion type distortion
given parametrs a, b, c. In test mode draws barrel or pincushion type 
lines on input frame with a, b, c . This feature enables to determine values for a,b and c.

distortion where  x axis distortion  is independant of y axis distortion also can be corrected.
(included at a request from an user)

vhr parameter corrects for unequal distortions along height and width of frame 
 seen on some anamorphic projections.
 
 Author V.C.Mohan

  sep 2014, modified on 20 Aug 2020

Copyright (C) <2020>  <V.C.Mohan>

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
		bool pin;	// whether pincushion type  or  barrel type
		bool test;	// whether a test run
		float vhRatio;	// ratio of h/v distortion
		bool ind;		// independant x and y distortions
			 
		float yabc[3]; // coefficients for y // starting  for  y
		bool ypin;	// in case of ind y has pincushion distortion otherwise barrel type
		float eabc[3];
		float eyabc[3]; // for test process end values if ind
		uint8_t col[3];		// color of dots in test
		int pixelsPerDot;	// in test
		
		float * xbuf, *ybuf; // LUT for correctios
		float * lanczos6;	// interpolation coefficients
	
		int quantile;		// accuracy of interpolation
		int span;			// span of lanczos 6 (X6)
	
}DeBarrelData;


	

	void makeLUT (float * buf, float a, float b, float c, float d, int x, int y,
									float vhrsqr, int rr,bool xtab );
	template<typename finc>
	void symetricMark4(finc * wp, int hOffset, int wOffset, finc color);
	template <typename finc>
	void SymmetricInterpolatedValues(finc * wp, int hOffset, int wOffset,
		const finc * sp, int hNearOffset, int wNearOffset, int spitch,
		int span, int qx, int qy, int q, float * lanc6, finc min, finc max);
		
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
	
	if (!d->test)
	{
		d->quantile = 32;	// up to 1/32 pixel position interpolation 
		d->span = 6;		// 6 X 6 point interpolation

		d->lanczos6 = vs_aligned_malloc <float>(d->span * (d->quantile + 1) * sizeof(float), 32);	// 1/quantile accuracy 6 coeff per quantile
		// create lanczos coefficients
		LanczosCoeff(d->lanczos6, d->span, d->quantile);

		for (int i = 0; i < 3; i++)
		{
			if (!d->pin)
				d->abc[i] = -d->abc[i];
			if (!d->ypin)
				d->yabc[i] = -d->yabc[i];
		}


		const int width = d->vi->width;
		const int height = d->vi->height;
		// if horizontal distortion less than vertical then vhRatio is more than 1.0
		// float vhRatio = 1.0;

		int hadjusted = d->ind ? height : int(height * d->vhRatio) & 0xfffffffe;

		int RR = hadjusted > width ? width / 2 : hadjusted / 2;	// normalised radius will be 1.0

		d->xbuf = vs_aligned_malloc <float>(width / 2 * height / 2 * sizeof(float), 32);

		d->ybuf = vs_aligned_malloc <float>(width / 2 * height / 2 * sizeof(float), 32);

		float ad = 1.0f - (d->abc[0] + d->abc[1] + d->abc[2]);

		float ydd = 1.0f - (d->yabc[0] + d->yabc[1] + d->yabc[2]);

		float vhRsq = d->vhRatio * d->vhRatio;

		// create a look up table

		makeLUT(d->xbuf, d->abc[0], d->abc[1], d->abc[2], ad, width / 2, height / 2,
			vhRsq, RR, true);

		if (!d->ind)
		{

			makeLUT(d->ybuf, d->abc[0], d->abc[1], d->abc[2], ad, width / 2, height / 2,
				vhRsq, RR, false);
		}
		else
		{

			makeLUT(d->ybuf, d->yabc[0], d->yabc[1], d->yabc[2], ydd, width / 2, height / 2,
				vhRsq, RR, false);
		}
	}

}
//--------------------------------------------------------------------------------------------

	static void VS_CC debarrelFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
	{
		DeBarrelData *d = (DeBarrelData *)instanceData;
		vsapi->freeNode(d->node);
		if (!d->test)
		{
			vs_aligned_free(d->lanczos6);
			vs_aligned_free(d->xbuf);
			vs_aligned_free(d->ybuf);
		}

		free(d);
	}
//---------------------------------------------------------------------
// creates a look up table with new values for x and y 
void makeLUT (float * buf, float a, float b, float c, float d,
					   int x, int y, float vhrsqr, int rr,bool xtab)
{	

	for ( int h = 0; h < y; h++)
	{

		for( int w = 0; w < x; w ++)
		{

			float rnormal = sqrt( h * h * vhrsqr + w * w ) / rr ;

			float rnew =  (((a * rnormal  

								+ b) * rnormal  

								+ c )* rnormal 

								+ d );
			
			if (xtab)
				
				buf [h * x + w] =  rnew * w;

			else 
				
				buf [h * x + w] =  rnew * h;
		}
	}
}

//------------------------------------------------------------------------------------
template<typename finc>
void symetricMark4(finc * wp, int hOffset, int wOffset, finc color)
{
	*(wp + hOffset + wOffset) = color;
	*(wp + hOffset - wOffset) = color;
	*(wp - hOffset + wOffset) = color;
	*(wp - hOffset - wOffset) = color;
}
 
//------------------------------------------------------------------------------------------------

template <typename finc>
void SymmetricInterpolatedValues(finc * wp, int hOffset, int wOffset,
	const finc * sp, int hNearOffset, int wNearOffset, int spitch,
	int span, int qx, int qy, int q, float * lanc6, finc min, finc max)
{
	*(wp + hOffset + wOffset) =
		clamp(LaQuantile(sp + hNearOffset + wNearOffset, spitch,
		span, qx, qy, lanc6), min, max);

	*(wp + hOffset - wOffset) =
		clamp(LaQuantile(sp + hNearOffset - (wNearOffset + 1), spitch,
		span, q - qx, qy, lanc6), min, max);

	*(wp - hOffset + wOffset) =
		clamp(LaQuantile(sp - (hNearOffset + spitch) + wNearOffset, spitch,
		span, qx, q - qy, lanc6), min, max);

	*(wp - hOffset - wOffset) =
		clamp(LaQuantile(sp - (hNearOffset + spitch) - (wNearOffset + 1), spitch,
		span, q - qx, q - qy, lanc6), min, max);
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

		float dd = 1.0f - ( d->abc[0] + d->abc[1] + d->abc[2]);
		float ydd = 1.0f - (d->yabc[0] + d->yabc[1] + d->yabc[2]);

		// When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
        // supply the "dominant" source frame to copy properties from. Frame props
        // are an essential part of the filter chain and you should NEVER break it.
       

		int nbits = (fi->sampleType == stInteger) ? fi->bitsPerSample : 0;
		int	nbytes = fi->bytesPerSample;

		// if horizontal distortion less than vertical then vhRatio is more than 1.0

		int hadjusted = d->ind ? height : int(height * d->vhRatio) & 0xfffffffe;

		int RR = hadjusted > width ? width / 2 : hadjusted / 2;	// normalised radius will be 1.0		

		float vhRsq = d->vhRatio * d->vhRatio;

		if (d->test)
		{

			uint8_t col[] = { d->col[0], d->col[1], d->col[2] };
			float fcol[] = { d->col[0] / 256.0f, d->col[1] / 256.0f, d->col[2] / 256.0f };

			if (fi->colorFamily == cmYUV)
			{
				BGRtoYUV(d->col, col);
				fcol[1] -= 0.5f;
				fcol[2] -= 0.5f;
			}
			// copy input frame on to dst
			dst = vsapi->copyFrame(src, core);			
			// we will draw a distortion grid on to the frame
			// determine  coefficients for current frame
			float nabc[3], nyabc[3];
			for (int i = 0; i < 3; i++)
			{
				nabc[i] = (d->abc[i] + n *(d->eabc[i] - d->abc[i]) / d->vi->numFrames);
				if (d->ind)
					nyabc[i] = (d->yabc[i] + n *(d->eyabc[i] - d->yabc[i]) / d->vi->numFrames);
				else
					nyabc[i] = nabc[i];
			}
			
			float ndd = 1.0f - (nabc[0] + nabc[1] + nabc[2]);
			float nydd = 1.0f - (nyabc[0] + nyabc[1] + nyabc[2]);

			for (int plane = 0; plane < fi->numPlanes; plane++)
			{
				uint8_t *dstp = vsapi->getWritePtr(dst, plane);
				int dpitch = vsapi->getStride(dst, plane) / nbytes; // this is in pixels. Required for uint16 and float pointers
				// note that if a frame has the same dimensions and format, the stride is guaranteed to be the same. int dst_stride = src_stride would be fine too in this filter.

				int center = height / 2 * dpitch + width / 2;	// offset to center		

				// we will mark every  pixelperdot
				for (int h = 0; h < height / 2; h += d->pixelsPerDot)
				{

					for (int w = 0; w < width / 2; w += d->pixelsPerDot)
					{
						// for each point get the coordinates of distorted image

						float rnormal = sqrt(h * h * vhRsq + w * w) / RR;

						float rnew = (((nabc[0] * rnormal	// ar^3

							+ nabc[1]) * rnormal		// br^2  

							+ nabc[2])* rnormal		// cr

							+ ndd);

						float wnew = rnew * w;
						// nearest point
						int wnearest = wnew + 0.5f;		// nearest int we will mark

						if (wnearest >= width / 2)

							continue;

						if (d->ind)
						{

							rnew = (((nyabc[0] * rnormal

								+ nyabc[1]) * rnormal

								+ nyabc[2])* rnormal

								+ nydd);
						}

						float hnew = rnew * h;

						int hnearest = hnew + 0.5f;

						if (hnearest >= height / 2)

							continue;

						if (fi->sampleType == stInteger)
						{

							if (nbytes == 1)
							{
								// on the frame we will mark the nearest points to indicate distortion grid
								symetricMark4(dstp + center, hnearest * dpitch, wnearest, (unsigned char)col[plane]);

							}

							else
							{
								uint16_t * dp = (uint16_t *)dstp;
								// on the frame we will mark the nearest points to indicate distortion grid
								symetricMark4(dp + center, hnearest * dpitch, wnearest, (uint16_t)(col[plane] << (nbits - 8)));

							}

						}

						else // floating point
						{

							float * dp = (float *)dstp;
							// on the frame we will mark the nearest points to indicate distortion grid
							symetricMark4(dp + center, hnearest * dpitch, wnearest, fcol[plane]);

						}
					} // for w
				} // for h

			}	// for plane
		}	// if test

		else	// not test. normal processing
		{
			uint8_t col[] = { 0, 0, 0 };
			float fcol[] = { 0, 0, 0 };
			if (fi->colorFamily == cmYUV)
			{
				col[1] = 127, col[2] = 127;
			}

			dst = vsapi->newVideoFrame(fi, width, height, src, core);

			

			for (int plane = 0; plane < fi->numPlanes; plane++)
			{
				const uint8_t *srcp = vsapi->getReadPtr(src, plane);
				int spitch = vsapi->getStride(src, plane) / nbytes;
				uint8_t *dstp = vsapi->getWritePtr(dst, plane);
				int dpitch = vsapi->getStride(dst, plane) / nbytes;

				int center = (height / 2) * spitch + width / 2;	// offset from left top to center of frame

				for (int h = 0; h < height / 2; h++)
				{
					for (int w = 0; w < width / 2; w++)
					{
						// for each point get the coordinates of distorted image

						float wnew = d->xbuf[h * width / 2 + w];

						int wnearest = wnew;

						float hnew = d->ybuf[h * width / 2 + w];

						int hnearest = hnew;

						if (wnearest >= width / 2 - 3 || hnearest >= height / 2 - 3)
						{
							// we have no points for interpolation. Make this part black
							//May be could do a nearest point to get a few more points. Later to do
							if (nbytes == 1)
							{
								// on the frame we will mark the nearest points to indicate distortion grid
								symetricMark4(dstp + center, h * dpitch, w, col[plane]);


							}

							else  if (nbytes == 2)
							{
								symetricMark4((uint16_t *)dstp + center, h * dpitch, w, (uint16_t)(col[plane] << (nbits - 8)));


							}

							else //	// float
							{
								symetricMark4((float*)dstp + center, h * dpitch, w, fcol[plane]);

							}

							continue;
						}

						int quantilex = d->quantile * (wnew - wnearest);

						int quantiley = d->quantile * (hnew - hnearest);

						if (nbytes == 1)
						{
							uint8_t min = 0, max = 255;
							// get the interpolated value , clamped to the permitted range of values
							SymmetricInterpolatedValues(dstp + center, h* dpitch, w,
								srcp + center, hnearest * spitch, wnearest, spitch,
								d->span, quantilex, quantiley, d->quantile, d->lanczos6, min, max);

						}

						else if (nbytes == 2)
						{
							uint16_t min = 0, max = (1 << nbits) - 1;
							uint16_t * dp = (uint16_t *)dstp;
							uint16_t * sp = (uint16_t *)srcp;
							// get the interpolated value , clamped to the permitted range of values
							SymmetricInterpolatedValues(dp + center, h* dpitch, w,
								sp + center, hnearest * spitch, wnearest, spitch,
								d->span, quantilex, quantiley, d->quantile, d->lanczos6, min, max);

						}

						else	//	// float
						{
							float max = plane == 0 ? 1.0f : fi->colorFamily == cmYUV ? 0.5f : 1.0f;
							float min = plane == 0 ? 0 : fi->colorFamily == cmYUV ? -0.5f : 0.0;
							float * dp = (float *)dstp;
							float * sp = (float *)srcp;
							// get the interpolated value , clamped to the permitted range of values
							SymmetricInterpolatedValues(dp + center, h* dpitch, w,
								sp + center, hnearest * spitch, wnearest, spitch,
								d->span, quantilex, quantiley, d->quantile, d->lanczos6, min, max);

						}

					}	// for (int w...)

				}// for (int h...)

			}	// for plane
		}	// not test
		
			// activation reason arAllFramesReady processing completed
		vsapi->freeFrame (src);
		return dst;		
	}	// activation reason
	return 0;
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
	for (int i = 0; i < 3; i++)
	{
		d.abc[i] = vsapi->propGetFloat(in, "abc", i, 0);
		if (d.abc[i] < 0.0 || d.abc[i] > 0.5)
		{
			vsapi->setError(out, "deBarrel: abc[] values must all be less than 0.5");
			vsapi->freeNode(d.node);
			return;
		}
	}

	if (d.abc[0] + d.abc[1] + d.abc[2] > 1.0)
	{
		vsapi->setError(out, "deBarrel: sum of all three abc values must be less than 1.0 ");
		vsapi->freeNode(d.node);
		return;
	}

	temp = !!vsapi->propGetInt(in, "pin", 0, &err);

	if (err)
		temp = 0;

	if (temp < 0 || temp > 1)
	{
		vsapi->setError(out, "deBarrel: pin  value if 0  indicate barrel type and 1  pincushion type distortion  ");
		vsapi->freeNode(d.node);
		return;
	}
	if (temp == 0)
		d.pin = false;
	else
		d.pin = true;

	temp = !!vsapi->propGetInt(in, "yind", 0, &err);

	if (err)
		temp = 0;

	if (temp < 0 || temp > 1)
	{
		vsapi->setError(out, "deBarrel: yind can have a value of 0 to indicate identical  or 1 for independant x and y  distortion  ");
		vsapi->freeNode(d.node);
		return;
	}

	if (temp == 0)
		d.ind = false;
	else
		d.ind = true;

	if (d.ind)
	{
		//d.ypin 
		temp = !!vsapi->propGetInt(in, "ypin", 0, &err);

		if (err)
			temp = 0;

		if (temp < 0 || temp > 1)
		{
			vsapi->setError(out, "deBarrel: ypin can have a value of 0 to indicate barrel type and 1 for pincushion type distortion  ");
			vsapi->freeNode(d.node);
			return;
		}

		if (temp == 0)
			d.ypin = false;
		else
			d.ypin = true;

		for (int i = 0; i < 3; i++)
		{
			d.yabc[i] = vsapi->propGetFloat(in, "yabc", i, &err);

			if (err)
				d.yabc[i] = d.abc[i];
			else if (d.yabc[i] < 0.0 || d.yabc[i] > 0.5)
			{
				vsapi->setError(out, "deBarrel: all three values of yabc must  be between 0.0 and 0.5 ");
				vsapi->freeNode(d.node);
				return;
			}
		}
		if (d.yabc[0] + d.yabc[1] + d.yabc[2] > 1.0)
		{
			vsapi->setError(out, "deBarrel: sum of all three values of yabc must  be less than 1.0 ");
			vsapi->freeNode(d.node);
			return;
		}

	}	// if ind

	d.vhRatio = vsapi->propGetFloat(in, "vhr", 0, &err);

	if (err)
		d.vhRatio = 1.0f;

	if (d.vhRatio > 10 || d.vhRatio < 0.1)
	{
		vsapi->setError(out, "deBarrel: vhr can be between 0.1 and 10.0 only ");
		vsapi->freeNode(d.node);
		return;
	}

	// d.test
	temp = !!vsapi->propGetInt(in, "test", 0, &err);

	if (err)
		temp = 0;

	if (temp < 0 || temp > 1)
	{
		vsapi->setError(out, "deBarrel: test can have a value of 0 to indicate regular process or 1 for distortion grid display ");
		vsapi->freeNode(d.node);
		return;
	}
	if (temp == 0)
		d.test = false;
	else
		d.test = true;

	if (d.test)
	{
		/* commented out to enable use of Functional Tools also
		if (d.vi->numFrames < 10 || d.vi->numFrames > 100)
		{
			vsapi->setError(out, "deBarrel: in test mode clip length must be between 10 and 100");
			vsapi->freeNode(d.node);
			return;
		}
		*/

		for (int i = 0; i < 3; i++)
		{
			d.eabc[i] = vsapi->propGetFloat(in, "eabc", i, &err);

			if (err)
				d.eabc[i] = d.abc[i];
			else if (d.eabc[i] < 0.0 || d.eabc[i] > 0.5)
			{
				vsapi->setError(out, "deBarrel: all three values of eabc must  be between 0.0 and 0.5 ");
				vsapi->freeNode(d.node);
				return;
			}
		}
		if (d.eabc[0] + d.eabc[1] + d.eabc[2] > 1.0)
		{
			vsapi->setError(out, "deBarrel: sum of all three values of eabc must  be less than 1.0 ");
			vsapi->freeNode(d.node);
			return;
		}
		
		if (d.ind)
		{
			for (int i = 0; i < 3; i++)
			{
				d.eyabc[i] = vsapi->propGetFloat(in, "eyabc", i, &err);

				if (err)
					d.eyabc[i] = d.yabc[i];
				else if (d.eyabc[i] < 0.0 || d.eyabc[i] > 0.5)
				{
					vsapi->setError(out, "deBarrel: all three values of eyabc must  be between 0.0 and 0.5 ");
					vsapi->freeNode(d.node);
					return;
				}
			}
			if (d.eyabc[0] + d.eyabc[1] + d.eyabc[2] > 1.0)
			{
				vsapi->setError(out, "deBarrel: sum of all three values of eyabc must  be less than 1.0 ");
				vsapi->freeNode(d.node);
				return;
			}
			/* commented out to enable use through FunctionalTools

			if (d.abc[0] == d.eabc[0] && d.abc[1] == d.eabc[1] && d.abc[2] == d.eabc[2]
				&& d.yabc[0] == d.eyabc[0] && d.yabc[1] == d.eyabc[1] && d.yabc[2] == d.eyabc[2])
			{
				{
					vsapi->setError(out, "deBarrel: All input a,b,c parameters are same as end parameters. At least one should be different ");
					vsapi->freeNode(d.node);
					return;
				}
			}
			*/
		}	// if yind

		uint8_t col[] = { 0,0,0 };
		temp = vsapi->propGetInt(in, "rgb", 0, &err);
		if (err)
			d.col[0] = 0;
		else if (temp < 0 || temp > 255)
		{
			vsapi->setError(out, "deBarrel: rgb color components must be between 0 and 255 ");
			vsapi->freeNode(d.node);
			return;
		}
		else
			col[0] = temp;

		for (int i = 1; i < 3; i++)
		{
			temp = vsapi->propGetInt(in, "rgb", i, &err);
			if (err)
				col[i] = col[i - 1];
			else if (temp < 0 || temp > 255)
			{
				vsapi->setError(out, "deBarrel: rgb color components must be between 0 and 255 ");
				vsapi->freeNode(d.node);
				return;
			}
			else
				d.col[i] = temp;

		}
		// internally col is stored as bgr. So correct it
		d.col[0] = col[2];
		d.col[1] = col[1];
		d.col[2] = col[0];
		

		const char* pd = vsapi->propGetData(in, "dots", 0, &err);

		if (err)
			d.pixelsPerDot = 8;
		else if (strcmp(pd, "many") == 0)
			d.pixelsPerDot = 4;
		else if (strcmp(pd, "mod") == 0)
			d.pixelsPerDot = 8;
		else if (strcmp(pd, "few") == 0)
			d.pixelsPerDot = 12;
		else
		{
			vsapi->setError(out, "deBarrel: dots can be \"many\" or \"mod\" or \"few\" only.");
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
    vsapi->createFilter(in, out, "deBarrel", debarrelInit, debarrelGetFrame, debarrelFree, fmParallel, 0, data, core);
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
    registerFunc("DeBarrel", "clip:clip;abc:float[];vhr:float:opt;pin:int:opt;"
	"yind:int:opt;ypin:int:opt;yabc:float[]:opt;test:int:opt;"
	"eabc:float[]:opt;eyabc:float[]:opt;"
	"rgb:int[]:opt;dots:data:opt", debarrelCreate, 0, plugin);
}

*/


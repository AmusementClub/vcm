/******************************************************************************
Correlation filter plugin for vapoursynth( 32 & 64 bit) by V.C.Mohan
This filter operates in freq domain (2d), takes 2 clips A and B as input
and  displays  cross correlation of A and B clips.

This plugin needs any one of libfftw3f-3.dll 32bit and 64bit of FFTW.org to reside in path
(may be windows\system32 folder, or wow)

Author V.C.Mohan.
28 July 2020, 22 May 2021

Copyright (C) <2020- 2021>  <V.C.Mohan>

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

/*
#include "windows.h"
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
#include <complex>
#include  <math.h>
//#include "fftw3.h"	// fftwlite used instead as I can understand it
#include "fftwlite.h"
#include "Factorize.cpp"
#define NYQUIST 512		// arbitrary definition of Nyquist spatial frequency
*/

#include <stdio.h>
typedef struct
{
	VSNodeRef *node[2];
	VSVideoInfo vi;
	const VSVideoInfo *avi;
	// float gamma;			// correlation scaling. Not used		
	bool center;			// display centered on frame? else left top corner reference
	bool txt;				// is a correlation file with x and y to be created?
	const char * filename;	// file having x and y  in text format. name of text output file with full path
	FILE *ofile;	// 
	int cx, cy;				// max search window for correlation width max +- x, height max +- y
	int sf, ef, every;
	int wbest;				// nearest higher width for best speed of fft
	int hbest;				// nearest higher height for best speed of fft
	int freqWidth;			// freq width for real data will be  wbest / 2 + 1.and used for filter design 
	int f2size;
	//int * fbuf;			//  x and y coord of max corr for each frame are stored here.
	
	fftwf_plan  pf, pinv;	// forward and inverse fft plansnal float[2] format.
#include "fftLateBindingClassParams.cpp"
	float* inBuf;
	fftwf_complex* outBuf;
	fftwf_complex* Bfreq;
}F2QCorrData;
//------------------------------------------------------------------
	
	void xCorrelate(fftwf_complex* Afreq, fftwf_complex* Bfreq, int frqwidth, int hbest, bool center);
	int xNormGamma(float * buf, int size);

	template <typename finc>
	void xTransferToDst(finc * dp, int dpitch, 
		int dwd, int dht, float * buf, int bestwd, int bestht, finc max);
	template <typename finc>
	int xFullProcess(F2QCorrData *d, 
		finc * dp, const finc * ap, const finc * bp, int abpitch, 
		int pwd, int pht, int dpitch, finc max);
	template <typename finc>
	void drawOrigin(float* inBuf, int wbest, int hbest, finc val);


/*************************************************
* The following is the implementation
* of the defined functions.
*************************************************/
//Here is the acutal constructor code used

static void VS_CC f2qcorrInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
	F2QCorrData *d = (F2QCorrData *)* instanceData;
	//vsapi->setVideoInfo(d->vi, 1, node);	
	const int fht = ((d->avi->height + 3) >> 2) << 2;
	const int fwd = ((d->avi->width + 3) >> 2) << 2;
	int* facbuf = (int*)vs_aligned_malloc(sizeof(int) * 64, 32);
	
	d->wbest = getBestDim(fwd, facbuf);
	d->hbest = getBestDim(fht, facbuf);
	vs_aligned_free(facbuf);
	d->vi.height = d->hbest;
	d->vi.width = d->wbest;
	d->vi.format = d->avi->format;
	d->vi.numFrames = d->avi->numFrames;
	d->vi.fpsDen = d->avi->fpsDen;
	d->vi.fpsNum = d->avi->fpsNum;
	d->vi.flags = d->avi->flags;
	vsapi->setVideoInfo(&d->vi, 1, node);
	d->freqWidth = d->wbest / 2 + 1;	// for real data this is the output width in the first transform

	d->f2size = d->hbest * d->freqWidth;
 	if (d->txt)
	{		
		errno_t err = fopen_s(&d->ofile, d->filename, "w");
		if (err != 0)
		{
			vsapi->setError(out, "FQCorr:init: could not open output text file");
			vsapi->freeNode(d->node[0]);
			vsapi->freeNode(d->node[1]);
			free(d);
			return;
		}
		// heading and explanation of file
		fprintf(d->ofile, "Correlation Shifts determined between sf = %d and ef = %d at intervals of %d \n",
			d->sf, d->ef, d->every);
		fprintf(d->ofile, " Frame number fn, x shift, y shift in search area and fx and fy shifts in frame\n");
		fprintf(d->ofile, " fn\tx\ty\tfx\tfy\n");
	
	}
#include "ConstructorCodeForLateBindingfft.cpp"

	if (!ok)
	{
		vsapi->setError(out, "FQCorr: could not load fft dll");
		vsapi->freeNode(d->node[0]);
		vsapi->freeNode(d->node[1]);
		if (d->txt && d->ofile != NULL)
			fclose(d->ofile);
		free (d);
		return;
	}
	// buffers 
	d->inBuf = (float*)d->fftwf_malloc(sizeof(float) * d->wbest * d->hbest);
	d->outBuf = (fftwf_complex*)d->fftwf_malloc(sizeof(fftwf_complex) * d->f2size);
	d->Bfreq = (fftwf_complex*)d->fftwf_malloc(sizeof(fftwf_complex) * d->f2size);

	
	// We require forward  and inverse transform plans.
	d->pf = d->fftwf_plan_dft_r2c_2d(d->hbest, d->wbest, d->inBuf, d->outBuf, FFTW_ESTIMATE);
	// inverse so complex to real used
	d->pinv = d->fftwf_plan_dft_c2r_2d(d->hbest, d->wbest, d->outBuf, d->inBuf, FFTW_ESTIMATE);	
	
}

//-----------------------------------------------------------------------------------------------------------

void xCorrelate(fftwf_complex* Afreq, fftwf_complex* Bfreq, float scale, int frqwidth, int hbest, bool center)
{
	// float scale = 1.0 / (wbest * hbest);
	// complex multiply with conjugate and scale down to compensate fft upscaling
	int fsize = hbest * frqwidth;
	float real, imag;	
	int wi = 1;
	if (center)
	{
		for (int n = 0; n < fsize; n++)
		{
			real = (Afreq[n][0] * Bfreq[n][0] + Afreq[n][1] * Bfreq[n][1]) * scale;
			imag = (Afreq[n][0] * Bfreq[n][1] - Afreq[n][1] * Bfreq[n][0]) * scale;
			Afreq[n][0] = real * wi;
			Afreq[n][1] = imag * wi;
			wi = -wi;
		}		
	}
	else // no centering
	{		

		for (int n = 0; n < fsize; n++)
		{
			real = (Afreq[n][0] * Bfreq[n][0] + Afreq[n][1] * Bfreq[n][1]) * scale;
			imag = (Afreq[n][0] * Bfreq[n][1] - Afreq[n][1] * Bfreq[n][0]) * scale;
			Afreq[n][0] = real ;
			Afreq[n][1] = imag ;
			
		}
	}
}
//---------------------------------------------------------------------------
int xNormGamma(float * buf,  int size)
{
	float maximum = buf[0];
	float minimum = buf[0];
	int xyPlus = 0;
	// get min and max values
	for (int n = 0; n < size; n++)
	{
		if (maximum < buf[n])
		{
			maximum = buf[n];
			xyPlus = n;
		}
		else if (minimum > buf[n])
			minimum = buf[n];
	}
	// normalize. Avoid divide by zero.
	if (maximum > minimum)
	{
		float mult = 1.0f / (maximum - minimum);	// mult is faster than div

		for (int n = 0; n < size; n++)
		{
			// for some reason this gamma if more than 1 does not work 
			// I expect normalized values less than 1 should become much less in value for gamma more than 1, but remain 1.0. Why? 
		//	buf[n] = pow((buf[n] - minimum) * mult, gamma);
			buf[n] = (buf[n] - minimum) * mult;
		}
	}

	else
	{
		for (int n = 0; n < size; n++)
			buf[n] = 0;

	}

	return xyPlus;
}
//---------------------------------------------------------------------------------------------
template<typename finc>
void xTransferToDst(finc * dp, int dpitch, 
	int dwd, int dht, float * buf, int bestwd, int bestht, finc max)
{
	for (int h = 0; h < dht; h++)
	{		
		for (int w = 0; w < dwd; w++)
		{
			dp[w] = (finc)(buf[w] * max);	// converts normalized data to actual values
		}
		dp += dpitch;
		buf += bestwd;
	}
}
template <typename finc>
void drawOrigin(float* inBuf, int wbest, int hbest, finc val)
{
	int origin = hbest / 2 * wbest + wbest / 2;

	int length = 20;
	for (int i = -length; i < length; i++)
		inBuf[origin + i * wbest] = val;
	for (int i = -length; i < length; i++)
		inBuf[origin + i ] = val;

}
//--------------------------------------------------------------------------------------------
template <typename finc>
int xFullProcess(F2QCorrData* d,  // forward fft of the two input frames
	finc * dp, 	// output frame
	const finc * ap,  const finc * bp, int abpitch, // input A and B frames, pitch of frames
	int pwd, int pht, int dpitch, finc max) //  frame dimensions and max pixel value
		// dimensions computed as best used for fft

{
	// transfer A frame data to fft input float buffer
	getRealInput2D(d->inBuf, ap, abpitch, pht, pwd, d->hbest, d->wbest, false);
	
	// 2d forward transform into Afreq complex buffer
	d->fftwf_execute_dft_r2c(d->pf, d->inBuf, d->outBuf);
	// transfer B Frame data to fft input float buffer
	getRealInput2D(d->inBuf, bp, abpitch, pht, pwd, d->hbest, d->wbest, false);
	
	// 2d forward transform into Bfreq complex buffer
 	d->fftwf_execute_dft_r2c(d->pf, d->inBuf, d->Bfreq);
	// as during fft scaleup occurs, need to scale down
	float scale = 1.0f / (d->wbest * d->hbest);	
	// cross correlation in freq domain is multiplication. Result in Afreq 
	xCorrelate(d->outBuf, d->Bfreq, scale, d->hbest, d->freqWidth, true);
	// inverse fft
	d->fftwf_execute(d->pinv);
	// normalize and apply gamma scaling. gamma part deleted as not working per expectation.
	int xyPlus = xNormGamma(d->inBuf, d->wbest * d->hbest);
	// transfer data to output frame
	//xTransferToDst(dp, abpitch, pwd, pht, d->inBuf, d->wbest, d->hbest, max);
	xTransferToDst(dp, dpitch, d->wbest, d->hbest, d->inBuf, d->wbest, d->hbest, max);
	
	return xyPlus;
}

//---------------------------------------------------------------------------------------------
static const VSFrameRef *VS_CC f2qcorrGetFrame(int n, int activationReason, void **instanceData,
	void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
{
	F2QCorrData *d = (F2QCorrData *)* instanceData;

	if (activationReason == arInitial) {
		// Request the source frames on the first call
		vsapi->requestFrameFilter(n, d->node[0], frameCtx);
		vsapi->requestFrameFilter(n, d->node[1], frameCtx);
	}
	else if (activationReason == arAllFramesReady)
	{
		const VSFrameRef *srcA = vsapi->getFrameFilter(n, d->node[0], frameCtx);
		const VSFrameRef *srcB = vsapi->getFrameFilter(n, d->node[1], frameCtx);
		// The reason we query this on a per frame basis is because we want our filter
		// to accept clips with varying dimensions. If we reject such content using d->vi
		// would be better.
		
		int ht = vsapi->getFrameHeight(srcA, 0);
		int wd = vsapi->getFrameWidth(srcA, 0);
		int owd = d->vi.width;
		int oht = d->vi.height;
		//VSFrameRef *dst;
		const VSFormat *fi = d->vi.format;
		VSFrameRef *dst = vsapi->newVideoFrame(fi, owd, oht, srcA, core);
				// use green for RGB and Y for YUV
		int plane = fi->colorFamily == cmRGB ? 1 : 0;	

		const uint8_t *srcpA = vsapi->getReadPtr(srcA, plane);
		const uint8_t *srcpB = vsapi->getReadPtr(srcB, plane);
		int src_stride = vsapi->getStride(srcA, plane);
		int dst_stride = vsapi->getStride(dst, plane);
		uint8_t *dstp = vsapi->getWritePtr(dst, plane);
		int pitch = src_stride / fi->bytesPerSample;
		int dpitch = dst_stride / fi->bytesPerSample;
		int nbits = fi->sampleType == stInteger ? fi->bitsPerSample : 0;
		int nbytes = fi->bytesPerSample;
		int pht = vsapi->getFrameHeight(srcA, plane);
		int pwd = vsapi->getFrameWidth(srcA, plane);

		int xyPlusW;

		if (nbytes == 1)
		{
			// 8 bit data
			uint8_t min = 0, max = 255;

			xyPlusW = xFullProcess(d, dstp,
								srcpA, srcpB,  pitch, 
								pwd, pht, dpitch, max);
		}

		else if (nbytes == 2)
		{
			// 16 bit data
			uint16_t min = 0, max = (1 << nbits) - 1;

			xyPlusW = xFullProcess(d,  (uint16_t*)dstp,
								(uint16_t *)srcpA, (uint16_t *)srcpB, pitch,
								pwd, pht, dpitch, max);
		}
		else
		{
			// float data
			float min = 0, max = 1.0f;

			xyPlusW = xFullProcess(d, (float*)dstp,
				(float*)srcpA, (float*)srcpB, pitch,
				pwd, pht, dpitch, max);
		}

		// Release the source frame
		vsapi->freeFrame(srcA);
		vsapi->freeFrame(srcB);
		
		if (fi->numPlanes > 1)
		{
			for (plane = 0; plane < fi->numPlanes; plane++)
			{
				if ((fi->colorFamily == cmRGB && plane == 1) || (fi->colorFamily == cmYUV && plane == 0))
					continue;
				uint8_t *dstp = vsapi->getWritePtr(dst, plane);
				int dpitch = vsapi->getStride(dst, plane) / fi->bytesPerSample;
				//	int spitch = src_stride / fi->bytesPerSample;
				int nbits = fi->sampleType == stInteger ? fi->bitsPerSample : 0;
				int nbytes = fi->bytesPerSample;
				
				int dht = vsapi->getFrameHeight(dst, plane);
				int dwd = vsapi->getFrameWidth(dst, plane);

				// zero out chroma enable correlation display
				if (nbytes == 1)
				{

					uint8_t grey = fi->colorFamily == cmRGB ? 0 : 128;
					fillPlaneWithVal(dstp, dpitch, dwd, dht, grey);
					//xFillPlaneWithVal(dstp, dpitch, dwd, dht, grey);
				}
				else if (nbytes == 2)
				{
					uint16_t grey = fi->colorFamily == cmRGB ? 0 : 1 << (nbits - 1);

					fillPlaneWithVal((uint16_t*)dstp, dpitch, dwd, dht, grey);
				}
				else
				{ // float data
					float grey = 0.0;

					fillPlaneWithVal((float*)dstp, dpitch, dwd, dht, grey);
				}
			}	/// for plane = 1;....
		}	// if num planes > 1
		if (fi->colorFamily == cmRGB)
		{
			// copy Green on to Blu and Red planes
			vs_bitblt(vsapi->getWritePtr(dst, 0), vsapi->getStride(dst, 0),
				vsapi->getWritePtr(dst, 1), vsapi->getStride(dst, 1),
				wd * nbytes, ht);
			vs_bitblt(vsapi->getWritePtr(dst, 2), vsapi->getStride(dst, 2),
				vsapi->getWritePtr(dst, 1), vsapi->getStride(dst, 1),
				wd * nbytes, ht);
		}

		if (d->txt && ((n - d->sf) % d->every) == 0 && n <= d->ef )
		{
			// find x and y coordinates of max correlation
			float * topLeft = d->inBuf + (d->hbest / 2 - d->cy) * d->wbest + d->wbest / 2 - d->cx;// left top of rectangle
			float max = 0;
			int xcoord = 0, ycoord = 0;
			// search in a rectangle of -cx, -cy to cx, cy
			for (int h = 0; h <= 2 * d->cy; h++)
			{
				for (int w = 0; w <= 2 * d->cx; w++)
				{
					if (topLeft[w] > max)
					{
						max = topLeft[w];
						xcoord = w;
						ycoord = h;
					}
				}
				topLeft += d->wbest;
			}
			int x = xcoord - d->cx, y = ycoord - d->cy;
			int fw = (xyPlusW % d->wbest) - d->wbest /2;
			int fh = (xyPlusW / d->wbest) - d->hbest / 2;
			//ofile << "Correlation  Frame number, xcoord, ycoord, frame xcoord, frame ycoord" << nl;
			int row = (n - d->sf) / d->every;
			// after every 10 th row a double line interval
			if (row % 10 == 0)
				fprintf(d->ofile, "\n");

			fprintf(d->ofile, " n %d\t x %d\t y %d\t fx %d\tfh %d\n",
				n, x, y, fw, fh);
			
		}

		if (n >= d->ef && d->ofile != NULL)
			fclose(d->ofile);


		return dst;
	}

	return 0;
}

// Free all allocated data on filter destruction
static void VS_CC f2qcorrFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
	F2QCorrData *d = (F2QCorrData *)instanceData;
	vsapi->freeNode(d->node[0]);
	vsapi->freeNode(d->node[1]);
	d-> fftwf_free(d->inBuf);
	d->fftwf_free(d->outBuf);
	d->fftwf_free(d->Bfreq);
	d->fftwf_destroy_plan(d->pf);
	d->fftwf_destroy_plan(d->pinv);
	if (d->txt && d->ofile != NULL)
		fclose(d->ofile);
	free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC f2qcorrCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
	F2QCorrData d;
	F2QCorrData *data;
	int err;
	int temp;
	// Get a clip reference from the input arguments. This must be freed later.
	d.node[0] = vsapi->propGetNode(in, "clip", 0, 0);
	d.avi = vsapi->getVideoInfo(d.node[0]);
	
	// get second clip
	d.node[1] = vsapi->propGetNode(in, "bclip", 0, 0);
	const VSVideoInfo *bvi = vsapi->getVideoInfo(d.node[1]);

	if (!isSameFormat(d.avi, bvi) || d.avi->numFrames != bvi->numFrames  )
	{
		vsapi->setError(out, "F2QCorr: both clips must be of same format, length and frame dimensions ");
		vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
		return;
	}
	if (d.avi->format->colorFamily == cmCompat)
	{
		vsapi->setError(out, "F2QCorr: compat format is not accepted. Only Planar format clips can be input ");
		vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
		return;
	}
	// In this first version we only want to handle 8bit integer formats. Note that
	// vi->format can be 0 if the input clip can change format midstream.
	if (!isConstantFormat(d.avi) || d.avi->width == 0 || d.avi->height == 0
		|| d.avi->width != bvi->width || d.avi->height != bvi->height)
	{
		vsapi->setError(out, "F2QCorr: only constant format and const frame dimensions input supported");
		vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
		return;
	}
	// set width and height of output frame
	//d.vi.height = d.hbest;
	//d.vi.width = d.wbest;

	// If a property read fails for some reason (index out of bounds/wrong type)
	// then err will have flags set to indicate why and 0 will be returned. This
	// can be very useful to know when having optional arguments. Since we have
	// strict checking because of what we wrote in the argument string, the only
	// reason this could fail is when the value wasn't set by the user.
	// And when it's not set we want it to default to enabled.
	temp = !!int64ToIntS(vsapi->propGetInt(in, "txt", 0, &err));
	if (err || temp == 0)
		d.txt = false;
	else
		d.txt = true;
	if (d.txt)
	{		
		temp = int64ToIntS(vsapi->propGetInt(in, "cx", 0, &err));
		if (err)
			d.cx = 20;
		else if (abs(temp) < 2 && abs(temp) > d.avi->width / 8)
		{
			vsapi->setError(out, "F2QCorr: absolute values of cx must be between 2 and 1/8 frame wwidth");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
		else
			d.cx = abs(temp);
		temp = int64ToIntS(vsapi->propGetInt(in, "cy", 0, &err));
		if (err)
			d.cy = d.cx <= d.avi->height / 8 ? d.cx : d.avi->height / 8;
		else if (abs(temp) < 2 && abs(temp) > d.avi->height / 8)
		{
			vsapi->setError(out, "F2QCorr: absolute values of cy must be between 2 and 1/8 frame height");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
		else
			d.cy = abs(temp);

		d.sf = int64ToIntS(vsapi->propGetInt(in, "sf", 0, &err));
		if (err)
			d.sf = 0;
		else if (d.sf < 0 || d.sf >= d.avi->numFrames - 1)
		{
			vsapi->setError(out, "F2QCorr: sf must be within clip");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
		d.ef = int64ToIntS(vsapi->propGetInt(in, "ef", 0, &err));
		if (err)
			d.ef = d.avi->numFrames - 1;
		else if (d.ef < d.sf || d.ef >= d.avi->numFrames)
		{
			vsapi->setError(out, "F2QCorr: ef must not be less than sf and must be within clip");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}

		temp = d.ef - d.sf <= 1000 ? 1 : (d.ef - d.sf) / 1000 + 1;
		d.every = int64ToIntS(vsapi->propGetInt(in, "every", 0, &err));
		if (err)
			d.every = temp;
		else if (d.every < temp || d.every >= d.ef - d.sf)
		{
			vsapi->setError(out, "F2QCorr: every should not result in either zero or over 1000 records");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
	
		const char * fn = vsapi->propGetData(in, "filename", 0, &err);
		if (err)
		{
			vsapi->setError(out, "F2QCorr: valid File name with full path must be specified");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
		
		d.filename = fn;
		
	}
		
	// I usually keep the filter data struct on the stack and don't allocate it
	// until all the input validation is done.
	data = (F2QCorrData *)malloc(sizeof(d));
	*data = d;

	vsapi->createFilter(in, out, "F2QCorr", f2qcorrInit, f2qcorrGetFrame, f2qcorrFree, fmParallelRequests, 0, data, core);
}

//registerFunc("fqCorr", "clip:clip;bclip:clip;cx:int:opt;cy:int:opt;
//		txt:int:opt;filename:data:opt;sf:int:opt;ef:int:opt;every:int:opt;", corrCreate, 0, plugin);



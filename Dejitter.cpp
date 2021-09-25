/******************************************************************************
DeJitter filter function of vcmisc plugin for vapoursynth.
by V.C.Mohan

Rows in a frame are displaced horizontally randomly by some video tape vhs players
during transcription due to hardware malfunction(probably some noise is treated as sync signal).
This plugin tries to correct it, by reading the values and guessing the shift.
Thread safe
Author V.C.Mohan
Created on 2 Aug 2020

Copyright (C) < 2020>  <V.C.Mohan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is at
http://www.gnu.org/licenses/.


Author : V.C.Mohan


*************************************************************************************************/
//------------------------------------------------------------------------------  
/*
#include <stdlib.h>
#include "windows.h"
#include <stdint.h>
#include "VapourSynth.h"
#include "VsHeper.h"
#include "colorconverter.h"
*/

typedef struct {
	VSNodeRef *node;
	const VSVideoInfo *vi;
	int jmax;	// max jitter to correct
	float thresh;	// threshold value for row start detection
	int wsyn;	// width of sync signal to be skipped
}DejitterData;

// support functions 
	template <typename finc>
	int	getShift(const finc * rp, int maxw, finc th);
	template <typename finc>
	void shiftRow( finc * rp, int shift, int w);
	template <typename finc>
	void fillColor(finc * rp, int shift, int w, finc color);
	template <typename finc>
	void extendColor(finc * rp, int shift, int w);	

	//-------------------------------------------------------------------------------
	template <typename finc>
	int	getShift(const finc * rp, int maxw, finc th)
	{

		for (int i = 0; i < maxw; i++)

			if (rp[i] >= th)

				return i;

		return maxw;
	}
	//-------------------------------------------------------------------------------
	template <typename finc>

	void shiftRow(finc * rp, int shift, int w)
	{
		for (int i = 0; i < (w - shift); i++)
		{
			rp[i] = rp[(i + shift)];
		}

	}
	//---------------------------------------------------------------------------------
	template <typename finc>

	void fillColor(finc * rp, int shift, int w, finc color)
	{
		for (int i = w - shift; i < w; i++)

			rp[i] = color;
	}
	//--------------------------------------------------------------------------------
	template <typename finc>

	void extendColor(finc * rp, int shift, int w)
	{
		for (int i = w - shift; i < w; i++)

			rp[i] = rp[(w - shift)];
	}


	
/*************************************************
 * The following is the implementation 
 * of the defined functions.
 *************************************************/
	static void VS_CC dejitterInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
	{
		DejitterData *d = (DejitterData *)* instanceData;
		vsapi->setVideoInfo(d->vi, 1, node);
		
	}
/***************************************************************/	

	/***************************************************************/
	// This is the main function that gets called when a frame should be produced. It will, in most cases, get
	// called several times to produce one frame. This state is being kept track of by the value of
	// activationReason. The first call to produce a certain frame n is always arInitial. In this state
	// you should request all the input frames you need. Always do it in ascending order to play nice with the
	// upstream filters.
	// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
	// do the actual processing.
	static const VSFrameRef *VS_CC dejitterGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
	{
		DejitterData *d = (DejitterData *)* instanceData;

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
			VSFrameRef *dst = vsapi->copyFrame(src, core);
			int nplanes = fi->numPlanes;
			int nbytes = fi->bytesPerSample;
			int nbits = fi->bitsPerSample;
			
			unsigned char * dp[] = { NULL, NULL, NULL, NULL };
			int pitch[4];

			int pht[4], pwd[4];
			for (int p = 0; p < nplanes; p++)
			{
				pitch[p] = vsapi->getStride(dst, p) / nbytes;
				dp[p] = vsapi->getWritePtr(dst, p);				

			}

			vsapi->freeFrame(src);

			for (int h = 0; h < height; h++)
			{
				int shift = 0;

				int np = (fi->colorFamily == cmRGB) ? 3 : 1;	//check on green or luma
				for (int p = 0; p < np; p++)
				{
					int shiftp;
					if (nbytes == 1)
					{

						shiftp = getShift(dp[p], d->jmax, (uint8_t)((1 << nbits) * d->thresh));

					}

					else if (nbytes == 2)
					{

						shiftp = getShift((uint16_t*)dp[p], d->jmax, (uint16_t)((1 << nbits) * d->thresh));

					}
					else // float
					{

						shiftp = getShift((float*)dp[p], d->jmax, d->thresh);

					}
					if (shiftp > shift && shiftp < d->jmax)
						shift = shiftp;
				}
				
				// check if the shift has been detected
				if (shift >0 && shift < d->jmax )
				{
					shift += d->wsyn;	// as sync signal might have been detected go past by its width

					for (int p = 0; p < nplanes; p++)
					{
						if (nbytes == 1)
						{
							uint8_t col = (fi->colorFamily == cmRGB) ? 0 : p > 0 ? 127 : 0;

							shiftRow(dp[p], shift, width);

							fillColor(dp[p], shift, width, col);
						}

						else if (nbytes == 2)
						{
							uint16_t col = (fi->colorFamily == cmRGB) ? 0 : p > 0 ? 127 << (nbits - 8) : 0;
							// if 16 bit samples						
							shiftRow((uint16_t *)dp[p], shift, width);

							fillColor((uint16_t *)dp[p], shift, width, col);

						}
						else	//nbytes == 4 float
						{
							float col = 0.0; // (fi->colorFamily == cmRGB) ? 0 : p > 0 ? 0 : 0;
							// if 32 bit samples						
							shiftRow((float *)dp[p], shift, width);

							fillColor((float *)dp[p], shift, width, col);

						}
					}	// for p =0;
				} 	// if shift > 0
				for (int p = 0; p < nplanes; p++)
				{
					
					dp[p] += pitch[p] * nbytes;
				}
			}	// for h =0;

			
			return dst;
		}
		return 0;
	}

	//----------------------------------------------------------------------
	// Free all allocated data on filter destruction
	static void VS_CC dejitterFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
	{
		DejitterData *d = (DejitterData *)instanceData;
		vsapi->freeNode(d->node);
		free(d);
		
	}
	//------------------------------------------------------------------

	// This function is responsible for validating arguments and creating a new filter
	static void VS_CC dejitterCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
	{
		DejitterData d;
		//   JitterData *data;
		int err;
		int temp;

		// Get a clip reference from the input arguments. This must be freed later.
		d.node = vsapi->propGetNode(in, "clip", 0, 0);
		d.vi = vsapi->getVideoInfo(d.node);

		// In this integer and float. Note that
		// vi->format can be 0 if the input clip can change format midstream.
		if (!isConstantFormat(d.vi)) // || d.vi->format->sampleType != stInteger || d.vi->format->bitsPerSample != 8)
		{
			vsapi->setError(out, "deJitter: format of clip must be constant");
			vsapi->freeNode(d.node);
			return;
		}

		if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV  && d.vi->format->colorFamily != cmGray)
		{
			vsapi->setError(out, "deJitter: only RGB or YUV or Gray color formats allowed");
			vsapi->freeNode(d.node);
			return;
		}

		if (d.vi->format->subSamplingH != 0 || d.vi->format->subSamplingW != 0)
		{
			vsapi->setError(out, "deJitter: color planes should have no subsampling");
			vsapi->freeNode(d.node);
			return;
		}

		temp = vsapi->propGetInt(in, "jmax", 0, &err);

		if (err)
			d.jmax = 40;
		else if (temp < 1 || temp > d.vi->width / 4)
		{
			vsapi->setError(out, "deJitter: jmax value must be between 1 and quarter of frame width");
			vsapi->freeNode(d.node);
			return;
		}
		d.jmax = temp;

		temp = vsapi->propGetInt(in, "wsyn", 0, &err);

		if (err)
			d.wsyn = 20;
		else if (temp < 0 || temp > 40)
		{
			vsapi->setError(out, "deJitter: wsyn width of sync signal be between 0 and 40");
			vsapi->freeNode(d.node);
			return;
		}
		else
			d.wsyn = temp;
		d.thresh = vsapi->propGetFloat(in, "thresh", 0, &err);

		if (err)
			d.thresh = 0.08f;
		else if (d.thresh < 0.01f || d.thresh > 0.5f)
		{
			vsapi->setError(out, "deJitter: thresh can be between 0.01 and 0.5 only");
			vsapi->freeNode(d.node);
			return;
		}
		
		DejitterData * data = (DejitterData *)malloc(sizeof(d));

		*data = d;


		/***************************************************************/
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
		vsapi->createFilter(in, out, "deJitter", dejitterInit, dejitterGetFrame, dejitterFree, fmParallel, 0, data, core);
	}

	//////////////////////////////////////////
	// Initial entry

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
	VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
	{
	configFunc("in.vcmohan.Jitter", "vcm", "VapourSynth Jitter plugin", VAPOURSYNTH_API_VERSION, 1, plugin);
	registerFunc("dejitter", "clip:clip;jmax:int:opt;wsyn:int:opt;thresh:float:opt;", dejitterCreate, 0, plugin);
	}
	*/

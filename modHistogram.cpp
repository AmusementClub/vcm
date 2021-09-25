/*********************************************************************************************
HistogramAdjust filter plugin for vapoursynth by V.C.Mohan

	Full frame  or windowed equalization or matching Histograms.

	For matching either a table or a frame of a clip can be used.

	YUV  and Gray scale formats only are accepted as input.

	Thread safe.

	22 jun 2015 20 Aug 2020
	 Copyright (C) <2005-2020>  <V.C.Mohan>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    http://www.gnu.org/licenses/.
	
 Author V.C.Mohan

 
*************************************************************************************************/ 
/*
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
#include "HistogramAdjustHelper.cpp"
*/
typedef struct {
				VSNodeRef *node[2];
				const VSVideoInfo *vi[2];
				int type;			// 1. equalize,2. match clip, 3 match %age table 4. match cummulative table
				int table[40];				// histogram values from this table  to be used for matching
				int nentries;
				int mf;				// frame number of example clip to be used for matching Histogram
				bool window;		// is windowed equalization required? true yes
		
				int limit;			// in windowed processing max %age increase allowed
		
				float *histbuf, *matchbuf;			// buffers to get frame, local and matching clip  histograms;
				float *lbuf, *cbuf;				// buffer for lookup table to match histograms

				

}HistogramAdjustData;

//-------------------------------------------------------------------------------------------------

static void VS_CC histogramadjustInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    HistogramAdjustData *d = (HistogramAdjustData *) * instanceData;
    vsapi->setVideoInfo(d->vi[0], 1, node);

	int nb = d->vi[0]->format->bitsPerSample;
	int maxvalue = nb <= 12 ? 1 << nb : 4096;	// corresponds to 12 bit depth

	d->matchbuf = NULL;

	if ( d->type != 1)
		// not equalization
		d->matchbuf = vs_aligned_malloc <float>(sizeof( float) * maxvalue , 32);	 
		
	if (d->type == 2)	// matching with given frame of a clip
	{		
					// so get the frame of that clip		
		const VSFrameRef *matchf = vsapi->getFrame(d->mf, d->node[1], NULL, NULL);
			
        // The reason we query this on a per frame basis is because we want our filter
        // to accept clips with varying dimensions. If we reject such content using d->vi
        // would be better.
        const VSFormat *matchfi = d->vi[1]->format;
        int mht = vsapi->getFrameHeight(matchf, 0);
        int mwd = vsapi->getFrameWidth(matchf, 0);
		const uint8_t *mfp = vsapi->getReadPtr(matchf, 0);
        int m_stride = vsapi->getStride(matchf, 0);
		int pitch = m_stride / matchfi->bytesPerSample;
            
		if (matchfi->sampleType == stInteger)
		{
			int mb = matchfi->bitsPerSample; 			

			if ( mb == 8)
			{
					// get its histogram
				getHistFromWindow8(mfp, pitch, mb,  mwd, mht, d->matchbuf);
			}

			else 
			{
				const uint16_t * mp = (const uint16_t *) mfp;				

				getHistFromWindow16( mp, pitch , mb, mwd, mht, d->matchbuf);
			}			
			
		}

		else		// floating format. limit bins to 4096
		{
			const float * mp = (const float *) mfp;
			
			getHistFromWindowf( mp, pitch , 12, mwd, mht, d->matchbuf);			
		}

		sigmaHist(d->matchbuf, d->matchbuf,maxvalue);

		vsapi->freeFrame( matchf);			

		vsapi->freeNode (d->node[1]);
	
	}
	else if (d->type == 3)		 
	{
		// construct histogram from table of luma and population %ages
		getHistFromTable(d->table, d->nentries, maxvalue, d->matchbuf);
		// progrssive sum. 
		sigmaHist(d->matchbuf, d->matchbuf, maxvalue);
	}
	else if (d->type == 4)
	{
		// construct histogram from table of luma and cummulative population %ages
		int count = getHistCummTable(d->table, d->nentries, maxvalue, d->matchbuf);

		if (count != maxvalue)
		{
			vsapi->setError(out, "in correct count ");
			free(d->matchbuf);
			vsapi->freeNode(d->node[0]);
			return;
		}
		
	}
}

//-------------------------------------------------------------------------------------

static void VS_CC histogramadjustFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
    HistogramAdjustData *d = (HistogramAdjustData *)instanceData;
    vsapi->freeNode(d->node[0]);
	if (d->matchbuf != NULL)
		vs_aligned_free (d->matchbuf);

    free(d);
}
//-------------------------------------------------------------------------------------

 // This is the main function that gets called when a frame should be produced-> It will in most cases get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all input frames you need-> Always do it i ascending order to play nice with the
// upstream filters.
// Once all frames are ready the the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC histogramadjustGetFrame
				(int n, int activationReason, void **instanceData, 
				void **frameData, VSFrameContext *frameCtx, 
				VSCore *core, const VSAPI *vsapi)
{
    HistogramAdjustData *d = (HistogramAdjustData *) * instanceData;

    if (activationReason == arInitial) 
	{
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node[0], frameCtx);
    }
	else if (activationReason == arAllFramesReady) 
	{
        const VSFrameRef *src = vsapi->getFrameFilter(n, d->node[0], frameCtx);
       
        const VSFormat *fi = d->vi[0]->format;

		int nbits = fi->bitsPerSample;
		
        VSFrameRef *dst = vsapi->newVideoFrame(fi, d->vi[0]->width, d->vi[0]->height, src, core);
		int fwidth;	// number of samples to get. Will be multiplied number of bytes per sample by NewVideoFrame
		if(fi->sampleType== stInteger)
		{
			fwidth = nbits > 12 ? 4096 : 1 << nbits;
			
		}
		else
			fwidth  = 4096;	// floating point
		
		float *histbuf = vs_aligned_malloc<float>(sizeof(float) * fwidth, 32);
		
		for(int plane = 0; plane < fi->numPlanes; plane ++)
		{
			const uint8_t *srcp = vsapi->getReadPtr(src, plane);
			int src_stride = vsapi->getStride(src, plane);
			uint8_t *dstp = vsapi->getWritePtr(dst, plane);
			int bht = vsapi->getFrameHeight(src, plane);            
			int bwd = vsapi->getFrameWidth(src, plane);
			int samplesize = fi->bytesPerSample; 
			int pitch = src_stride / samplesize;
			if (plane > 0  && fi->colorFamily == cmYUV)
			{
				vs_bitblt(dstp, src_stride, srcp, src_stride, bwd * samplesize, bht);

				continue;
			}
	
			bool windowed = false;
				// check whether frame dimensions are adequate for windowed processing
			if (d->window == 1 && bht >= 200 && bwd >= 200	)
			{			
				windowed = true;
			}
		
				// processing frame
			if( ! windowed)
			{
					// full frame processing
				if (fi->sampleType == stInteger)
				{
					int nb = fi->bitsPerSample;

					int maxval = fwidth;

					if(nb == 8)
					{
						const uint8_t * sp = (const uint8_t *)srcp;
						uint8_t * dp = (uint8_t *)dstp;				

						getHistFromWindow8( sp, pitch, nb,  bwd , bht , histbuf);

						sigmaHist(histbuf, histbuf, maxval);

						if (d->type == 1)
						{

							fillAdjustedValues8(sp , dp , pitch, bwd, bht, nb, histbuf,  d->limit );
						}
						else
							fillMatchedValues8(sp , dp , pitch , bwd, bht, nb, histbuf, d->matchbuf,  d->limit );

					}

					else
					{
						const uint16_t * sp =  (const uint16_t *)srcp;
						uint16_t * dp = (uint16_t *)dstp;					
					
						getHistFromWindow16(  sp, pitch , nb,  bwd , bht , histbuf);
					
						sigmaHist(histbuf, histbuf, maxval);

						if( d-> type == 1)
						{

							fillAdjustedValues16(sp , dp , pitch , bwd, bht, nb, histbuf,  d->limit );
						}
						else

							fillMatchedValues16(sp , dp , pitch , bwd, bht, nb, histbuf,  d->matchbuf,  d->limit );

					}
				}

				else	// floating point
				{
					int nb = 12;	// limit to 4096 values
				
					const float * sp = (const float *)srcp;
					float * dp = (float *) dstp;				

					getHistFromWindowf( sp, pitch , nb,  bwd , bht , histbuf);

					sigmaHist(histbuf, histbuf, 4096);

					if (d->type == 1)
					{

						fillAdjustedValuesf(sp , dp , pitch , bwd, bht, nb, histbuf,  d->limit );
					}

					else

						fillMatchedValuesf(sp , dp , pitch , bwd, bht, nb, histbuf, d->matchbuf,  d->limit );
				}
			}

			else //if (windowed)
			{
				int nb = fi->bitsPerSample;
					// window dimensions
				int wside =  nb > 10? 128 : 64;	// 128 x 128 or  64 x 64 ensure sufficient statistics
				int innerw = wside / 4;	// 32x32 or 16x16
				int iwinstart = (wside - innerw) /2 ;	// 48 or  24
				int iwinend = wside - innerw / 2; // 112 or 56		
				
				int woverlap = innerw; // window overlap dimensions

				int offset, wstart, wd, ht;

				for( int w = 0; w < bwd - wside; w += woverlap)
				{
					for ( int h = 0; h < bht - wside; h += woverlap)
					{			
				
						if ( h == 0)
						{	// top row processing 
							if ( w <= bwd - wside * 2)
							{
								// left top corner or less than right top corner
								wstart = w == 0? 0 : w + iwinstart;
								wd = w == 0 ? iwinend : innerw;
								ht = iwinend;
								offset = wstart;
							
							}

							else if( w > bwd - wside * 2 )
							{
								// top right corner
								wstart = w + iwinstart;
								wd = bwd - w - iwinstart;
								ht = iwinend;
								offset = wstart;
						
							}
						}	// h == 0?

						else if ( h > bht - 2 * wside)
						{// bottom row processing
							if ( w <= bwd - wside * 2)
							{
								// left bottom corner  up to right bottom corner
								wstart = w == 0? 0 : w + iwinstart;
								wd = w == 0 ? iwinend : innerw;
								ht = bht - h - iwinstart;
								offset = ( h + iwinstart ) * pitch + wstart;
							
							}

							else if ( w > bwd - 2 * wside)
							{
								// right bottom corner
								wstart =  w + iwinstart;
								wd =  bwd - w - iwinstart;
								ht =  bht - h - iwinstart;
								offset = ( h + iwinstart ) * pitch + wstart;

							
							}
						}
						//  process left margin excluding corners
						else if ( w == 0)
						{
							wstart = 0;
							wd = iwinend;
							ht = innerw;
							offset = (h + iwinstart ) * pitch;
						
						
						}
							// right margin
						else if ( w > bwd - 2 * wside)
						{
							wstart = w + iwinstart;
							wd = bwd - w - innerw;
							ht = innerw;
							offset = ( h + iwinstart ) * pitch + wstart;
						}

						else		//well inside frame
						{
							wstart = w + iwinstart;
							wd = innerw;
							ht = innerw;
							offset = ( h + iwinstart ) * pitch + wstart;
						}

						if (fi->sampleType == stInteger)
						{
							int nb = fi->bitsPerSample;				

							if(fi->bitsPerSample == 8)						
							{
								const uint8_t *sp = (const uint8_t *) srcp; 
								uint8_t *dp = (uint8_t *) dstp;

								getHistFromWindow8( sp + h * pitch + w , pitch, nb, wside, wside, histbuf);

								sigmaHist( histbuf, histbuf, 1 << nb );

								fillAdjustedValues8(sp + offset, dp + offset, pitch, wd, ht, nb, histbuf,  d->limit );			

							}

							else
							{
								// 10 or 12 or 16 bit depth samples
								const uint16_t * sp = (const uint16_t *) srcp;
								uint16_t * dp = (uint16_t *) dstp;
								int nmax = nb <= 12 ? 1 << nb : 4096;
								getHistFromWindow16( sp + h * pitch + w , pitch, nb, wside, wside, histbuf);

								sigmaHist( histbuf, histbuf, nmax );
							
								fillAdjustedValues16(sp + offset, dp + offset, pitch, wd, ht, nb, histbuf,  d->limit );

							}
						}

						else
						{// floating point

							const float * sp = (const float *) srcp;
							float * dp = (float *) dstp;
							int nb = 12;

							getHistFromWindowf( sp + h * pitch + w , pitch, nb, wside, wside, histbuf);

							sigmaHist( histbuf, histbuf, 1 << nb );

							fillAdjustedValuesf(sp + offset, dp + offset, pitch, wd, ht, nb, histbuf,  d->limit );
						}

					}	//for  int h

				}	// for int w
			}	// window

		}
	
		vsapi->freeFrame (src);	
		vs_aligned_free(histbuf);
		return dst;	
	}	// get frame

	return 0;
}


// Calls the constructor with the arguments provied->
// This function is responsible for validating arguments and creating a new filter
static void VS_CC histogramadjustCreate(const VSMap *in, 
					VSMap *out, void *userData, VSCore *core,
					const VSAPI *vsapi) 
{

    HistogramAdjustData d;
    HistogramAdjustData *data;
//  VSNodeRef *cref;
    int err;
	int temp;


    // Get a clip reference from the input arguments. This must be freed later.
    d.node[0] = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi[0] = vsapi->getVideoInfo(d.node[0]);
	d.node[1] = NULL;
    // In this first version we only want to handle 8 bit integer formats. Note that
    // vi->format can be 0 if the input clip can change format midstream.
    //if (!isConstantFormat(d->vi) || d->vi->format->sampleType != stInteger || d->vi->format->bitsPerSample != 8) {
    //    vsapi->setError(out, "HistogramAdjust: only constant format 8 bit integer input supported");
    //    vsapi->freeNode(d->node);
    //    return;

	if (!isConstantFormat(d.vi[0]) )
	{
		vsapi->setError(out,"hist accepts const  format clip only");
		vsapi->freeNode(d.node[0]);
		return;
    }

	if(  d.vi[0]->format->colorFamily != cmYUV && d.vi[0]->format->colorFamily != cmGray)
	{
		vsapi->setError(out,"hist accepts YUV and Gray formats only");
		vsapi->freeNode(d.node[0]);
		return;
    }

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned-> This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string the only reason
    // this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled->
    d.type = vsapi->propGetInt(in, "type",0, &err);
	if(err)
		d.type = 1;
	if(d.type < 1 || d.type > 3)
	{
		// 1 =  equalisation, 2 = match using table values 3 = match using a specified frame of matching clip
        vsapi->setError(out, "hist: type can have value of 1 for equalization, 2 match with frame of mclip, 3.match table of luma and population %ages or 4 match with table of luma and cummulative population %ages.y");
        vsapi->freeNode(d.node[0]);
        return;
    }

	if (d.type == 2)
	{
		

		int m = vsapi->propNumElements(in, "table");

		if ( m < 4 || (m & 1) != 0 || m > 40 )
		{
			vsapi->setError(out,"hist atleast 2 and not more than 20 pairs   of table values for matching to be specified");
			vsapi->freeNode(d.node[0]);
			
			return;
		}

		for(int i = 0; i < m; i ++)
			d.table[i] = vsapi->propGetInt(in, "table", i, 0);

		for( int i = 0, j = -1, k = 0; i < m; i += 2)
		{
			// check pairs of values
			if (d.type == 3)
			{
				if (d.table[i] <= j || d.table[i] > 100 || d.table[i + 1] < 0 || d.table[i + 1] > 100)
				{
					vsapi->setError(out, "hist first member values of  pairs must be in ascending order and not more than 100. the second value of pair must be between 0 and 100");
					vsapi->freeNode(d.node[0]);
					return;
				}
			}

			else if (d.type == 4)
			{
				
				if (d.table[i] <= j || d.table[i] > 100 || d.table[i + 1] < k || d.table[i + 1] > 100)
				{
					vsapi->setError(out, "hist  luma and cummulative populatio values of  pairs. Both must be in ascending order 0 to 100.");
					vsapi->freeNode(d.node[0]);
					return;
				}
				
			}

			j = d.table[i];
			k = d.table[i + 1];
		}

		d.nentries = m;

	}	

	else if(d.type == 3)
	{		
		 d.node[1] = vsapi->propGetNode(in, "clipm", 0, &err);
		
		 if (err)
		 {
			d.node[1] = vsapi->cloneNodeRef(d.node[0]);	// default use first clip only
			
		 }

		 d.vi[1] =  vsapi->getVideoInfo(d.node[1]);

		 if (!isConstantFormat(d.vi[1]) )
		{
			vsapi->setError(out,"hist accepts for matching const YUV  format clipm only");
			
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}	
		
		if(d.vi[1]->format->colorFamily != cmYUV)
		{
			vsapi->setError(out,"hist accepts for matching clipm YUV  formats only. ");			
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
		
		d.mf = vsapi->propGetInt(in, "fm", 0, &err);
		
		if (err)
			d.mf = 0;

		if( d.mf < 0 || d.mf >= d.vi[1]->numFrames )
		{
			vsapi->setError(out,"hist clip for matching have fewer frames than Frame number specified");
			
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
	}
	temp = !!vsapi->propGetInt(in, "window", 0, &err);

	if (err)

		d.window = false;	// default use full frame

	else if(temp < 0 || temp > 1)
	{
		vsapi->setError(out,"hist window can have a value of 0 or 1 only ");
		
		vsapi->freeNode(d.node[0]);
		if (d.node[1] != NULL)
			vsapi->freeNode(d.node[1]);
		
		return;
	}
	if (temp == 0)
		d.window = false;
	else
		d.window = true;

		
	d.limit = vsapi->propGetInt(in, "limit", 0, &err);

	if (err)
	
		d.limit  = 0;	// default 100% adjust fully

	
	if (d.limit < 0 || d.limit > 99 ) 
	{
        vsapi->setError(out, "hist: limit a %age value  can be an integer from 0 to 99 only");
        
		vsapi->freeNode(d.node[0]);
		if (d.node[1] != NULL)
			vsapi->freeNode(d.node[1]);
        return;
    }

    // I usually keep the filter data struct on the stack and don't allocate it

    // until all the input validation is done.

  data =  (HistogramAdjustData *) malloc(sizeof(d));
  *data = d;

  

    // Create a new filter and returns a reference to it. Always pass on the in and out
    // arguments or unexpected things may happen. The name should be something that's
    // easy to connect to the filter, like its function name.
    // The three function pointers handle initialization, frame processing and filter destruction.
    // The filtermode is very important to get right as it controls how threading of the filter
    // is handled-> In general you should only use fmParallel whenever possible. This is if you
    // need to modify no shared data at all when the filter is running.
    // For more complicated filters fmParallelRequests is usually easier to achieve as an
    // be prefetched in parallel but the actual processing is serialized->
    // The others can be considered special cases where fmSerial is useful to source filters and
    // fmUnordered is useful when a filter's state may change even when deciding which frames to
    // prefetch (such as a cache filter).
    // If you filter is really fast (such as a filter that only resorts frames) you should set the
    // nfNoCache flag to make the caching work smoother.
    vsapi->createFilter(in, out, "HistogramAdjust", histogramadjustInit, 
								histogramadjustGetFrame, histogramadjustFree, 
								fmParallel, 0, data, core);
   // return;


}

// This is the entry point that is called when a plugin is loaded-> You are only supposed
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
    configFunc("vc.mohan.histo", "cmv", "VapourSynth Histogram Adjuster", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("hist", "clip:clip;clipm:clip:opt;type:int:opt;table:int[]:opt;mf:int:opt;window:int:opt;limit:int:opt", histogramadjustCreate, 0, plugin);
}
*/

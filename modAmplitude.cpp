/******************************************************************************
Amp(litude) filter plugin for vapoursynth by V.C.Mohan
Amplitude uses watershed algorithm after vincent- soille for segmentation of image.
Using basin numbers tagged by watershed as guide, the values are smoothed. Watershed is used for sharpening
No restriction on input format
Author V.C.Mohan
Oct, 2014, 20 Aug 2020

  Copyright (C) <2014 - 2020>  <V.C.Mohan>

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
#include <algorithm>
#include <queue>
#include "math.h"
#include <vector>
#include <stack>		
#include <functional>
#include "WSSegment.cpp"
#include "VapourSynth.h"
#include "VSHelper.h"
*/
//---------------------------------------------------------------------------------
typedef struct 
{
    VSNodeRef *node[2];	
    const VSVideoInfo *vi[2];

	bool connect4; // if false connect 8 used
	int sh[3];
    int sm[3];
	bool sclip; // use separate clip for segment

} AmplitudeData;
//---------------------------------------------------------------------------------------
// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC amplitudeInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) 
{
    AmplitudeData *d = (AmplitudeData *) * instanceData;
    vsapi->setVideoInfo(d->vi[0], 1, node);
	
}
//--------------------------------------------------------------------------


template <typename finc>
void SegmentSmooth(finc * dstr, int dspitch, 
				   finc * buff, int * tagg, int pitch,				   
				    int wdth, int hgt, int smf);

template <typename finc>
void SegmentSharp( finc * dstr, int dspitch, const finc * frp, int frpitch, 
								 char * wsh, int wspitch, int wdt, int hgt,int sh, finc min, finc max);

template < typename finc>
void fillBuffer(finc * bufp, int bpitch, const finc * srp, int srpitch,int wd, int ht);

template <typename finc>
void preSegmentProcess(finc * bufr, const int bpitch,
							  const finc * srp, const int srpitch,
							  const int wd, const int ht,
							  const finc **pix);

template <typename finc>
void postSegmentProcess(finc * wkp, const int wkpitch,
									  const finc * frp, const int frpitch,									  
									  const int wd, const int ht,
									  finc * bufr,int * tagg, char * wsh, const int bpitch,
									  int smooth, int sharp, finc min, finc max);
//-----------------------------------------------------------------------------

class LessThan{
	public :

	template <typename finc>
	bool operator()( const finc *f1,const finc *f2)
	//bool operator()( const unsigned char *f1,const unsigned char *f2)
				{ return *f1 < *f2; }
};
/*
template <typename finc>
bool LessThan(const finc * f1, const finc * f2)
{
	return *f1 < *f2;
}
*/
//-----------------------------------------------------------------------------
template <typename finc>
void SegmentSmooth(finc * dstr, int dspitch, 
				   finc * buff, int * tagg, int pitch,				   
				    int wdth, int hgt, int smf)
{
		// smooths values that are within same basins whose numbers >> sqz

	for(int h = 0; h < hgt; h++)
	{

		for(int w = 0; w < wdth ; w++)
		{
			int count = smf;

			float sum = (float)(buff[ w ] * smf);

			int basin = tagg[ w] ;

			for ( int hh = -1, hhpitch = -pitch; hh < 2; hh ++, hhpitch += pitch)
			{

				for ( char ww = - 1; ww < 2; ww++)
				{

					if( h + hh >= 0 && h + hh < hgt && w + ww >= 0 && w + ww < wdth)
					{
						if( tagg [  hhpitch + w + ww]  == basin)
						{
							if( hh == 0 || ww == 0)	// near 4 points or itself
							{
								sum += (buff[ hhpitch + (w + ww)]) * 2;

								count += 2;
							}

							else
							{// corner points so farther
								sum += buff[hhpitch + (w + ww)];

								count ++;
							}
								
						}
					}
				}
			}
					
			dstr [w ] = (finc)(sum / count);
		}

		dstr += dspitch;
		tagg += pitch;
		buff += pitch;
	}

}

//-------------------------------------------------------------------------------------
template <typename finc>
void SegmentSharp( finc * dstr, int dspitch, const finc * frp, int frpitch, 
								 char * wsh, int wspitch, int wdt, int hgt,int sh, finc min, finc max)
{
	for(int h = 0;	h < hgt ; h++)
	{
	
		for(int w = 0; w < wdt  ; w++)
		{	
			if(wsh[ w ] == 0)		// yes this is a edge
			{
				float val = (float)((frp[w] * ( 100 + sh))/100.0) ;

				dstr [w] = (finc)(val  > max ? max : val < min ? min : val);
			}
		}

		dstr += dspitch;
		frp += frpitch;
		wsh += wspitch;
	}
	
}

//....................................................................................

template < typename finc>
void fillBuffer(finc * bufp, int bpitch, const finc * srp, int srpitch,int wd, int ht)
{
	
	for (int h = 0;	h < ht; h++)
	{

		for ( int w = 0; w < wd; w++)
		{
						
			bufp[ w ] = srp[w ];
		}

		bufp += bpitch;
		srp += srpitch;
	}
}
//------------------------------------------------------------------------------
template <typename finc>
void preSegmentProcess(finc * bufr, const int bpitch,
							  const finc * srp, const int srpitch,
							  const int wd, const int ht,
							  const finc **pix)
{

		// fill buf with values of  frame either from clip or sclip

	fillBuffer(bufr, bpitch, srp, srpitch, wd, ht);
			
	int psize = wd * ht;

	for(int h = 0; h < psize; h ++)
						
		// sort buffer fill with pointers to data
		pix[h] = bufr + h ;

		// default sort is in ascending order.  this is done in place
		// As we are sorting pointers, our LessThan function is supplied
			
	std::sort(pix, pix + psize , LessThan());
}

//...............................................................................

//............................................................................................
template <typename finc>
void postSegmentProcess(finc * wkp, const int wkpitch,
									  const finc * frp, const int frpitch,									  
									  const int wd, const int ht,
									  finc * bufr,int * tagg, char * wsh, const int bpitch,
									  int smooth, int sharp, finc min, finc max)
{
		// now we can Segment smooth

	if ( smooth > 0)				

		SegmentSmooth(wkp, wkpitch,  bufr, tagg, bpitch, wd, ht, 10 - smooth);

				// sharpen if required
	if ( sharp != 0)
	{
					// As we have not filled output frame so far and sharp fills only a few values
		

		SegmentSharp( wkp, wkpitch, frp , frpitch, 
						wsh, bpitch, wd, ht, sharp , min, max);

	}


}
//...............................................................................
// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC amplitudeGetFrame(int n, int activationReason, void **instanceData,
							void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    AmplitudeData *d = (AmplitudeData *) * instanceData;

    if (activationReason == arInitial)
	{
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node[0], frameCtx);

		if( d->sclip)

			vsapi->requestFrameFilter(n, d->node[1], frameCtx);
    }
	else if (activationReason == arAllFramesReady)
	{
        const VSFrameRef *src = vsapi->getFrameFilter(n, d->node[0], frameCtx);

		const VSFormat *fi = d->vi[0]->format;
        int height = vsapi->getFrameHeight(src, 0);
        int width = vsapi->getFrameWidth(src, 0);
		
		const VSFrameRef * s2frame = d->sclip? vsapi->getFrameFilter(n, d->node[1], frameCtx) : src;
        
        VSFrameRef *dst = vsapi->newVideoFrame(fi, width, height, src, core);
        // The reason we query this on a per frame basis is because we want our filter
        // to accept clips with varying dimensions. If we reject such content using d->vi
        // would be better.
       

		int psize = height * width ;
		// get a work frame to serve as buffers
		VSFrameRef *work = vsapi->newVideoFrame(fi, width * ( sizeof(int) + sizeof(int) + 1 + sizeof(void*) + fi->bytesPerSample), height, src, core);
		uint8_t* wp = vsapi->getWritePtr(work, 0);
		// for large buffers it seems asking a new video frame is faster
		//	vs_aligned_malloc <uint8_t>(psize * (sizeof(int) + sizeof(int) + 8 + sizeof(void*) + fi->bytesPerSample), 32);
			
		// specify pointers to individual buffers
		int * dist = (int *)wp;			// - wp
		//	int * dist = ( int *) malloc (sizeof(int) * psize);
		int * tag = dist + psize;// = wp + sizeof int * psize
	//	int * tag = ( int *) malloc (sizeof(int) * psize);
		char * ws = (char *)(tag + psize);		// wp + (sizeof int + sizeof int) * psize
	//	char * ws = ( char *) malloc(psize);
		const uint8_t ** pixelsort8 = (const uint8_t **)(ws + psize);	//  wp + (sizeof int + sizeof int) * psize) + psize
	//	const uint8_t ** pixelsort8 = ( const uint8_t **) malloc(psize * sizeof(void *));
		uint8_t * buf8 = (uint8_t *)(pixelsort8 + psize);	//  wp + (sizeof int + sizeof int) * psize) + psize + (sixe of void *) * psize
															// wp + psize( sizeof int + sizeofint + 1 + sizeof void*)
	//	uint8_t * buf8 = ( uint8_t *) malloc(psize * fi->bytesPerSample);
		

        // It's processing loop time!
        // Loop over all the planes
     
        for (int plane = 0; plane < fi->numPlanes; plane++) 
		{
            const uint8_t *srcp = vsapi->getReadPtr(src, plane);
            int src_stride = vsapi->getStride(src, plane);
            uint8_t *dstp = vsapi->getWritePtr(dst, plane);
            int dst_stride = vsapi->getStride(dst, plane); // note that if a frame has the same dimensions and format, the stride is guaranteed to be the same. int dst_stride = src_stride would be fine too in this filter.
            // Since planes may be subsampled you have to query the height of them individually
            int ht = vsapi->getFrameHeight(src, plane);            
            int wd = vsapi->getFrameWidth(src, plane);
            int kb = fi->bytesPerSample;
			int pitch = src_stride / kb;
			const uint8_t * s2ptr = d->sclip ? vsapi->getReadPtr(s2frame, plane) : srcp;
			
			if ( d->sm[plane] == 0 && d->sh[plane] != 0)
			{
				vs_bitblt(dstp, dst_stride, srcp, src_stride, wd * fi->bytesPerSample, ht);
			}

			
			if( d->sm[plane] != 0 || d->sh[plane] != 0)
			{
				if ( fi->sampleType == stInteger)
				{
					int nb = fi->bitsPerSample;

					if ( nb == 8)
					{
						const uint8_t * sp = srcp;
						uint8_t * dp = dstp;
						// char * ws = ws8;
						const uint8_t ** pixelsort = pixelsort8;
						uint8_t * buf = buf8;
						const uint8_t * s2p = s2ptr;
						uint8_t min = 0, max = (1 << nb) - 1;
					
						// uses appropriate frame as d->sclip points to that clip frame
						preSegmentProcess(buf, wd, s2p , pitch,
									 wd, ht,  pixelsort);

						WSSegment( pixelsort, buf, tag, dist, ws, ht, wd, d->connect4);

						if ( d->sclip)

							fillBuffer(buf,wd, sp, pitch,wd, ht);

						postSegmentProcess(dp , pitch, sp, pitch, wd, ht,
									  buf, tag, ws, wd, d->sm[plane],d-> sh[plane], min, max ); 
			
					}

					else // if nb > 8
					{
	
						const uint16_t * sp = ( uint16_t *)srcp;
						uint16_t * dp = ( uint16_t *)dstp;
						
						const uint16_t ** pixelsort = (const uint16_t **) pixelsort8;
						uint16_t * buf = ( uint16_t *) buf8;
						const uint16_t * s2p = ( uint16_t *)s2ptr;
						uint16_t min = 0, max = (1 << nb) - 1;
						
						// uses appropriate frame as d->sclip points to that clip frame
						preSegmentProcess(buf, wd, s2p , pitch,	wd, ht,  pixelsort);

						WSSegment( pixelsort, buf, tag, dist, ws, ht, wd, d->connect4);

						if ( d->sclip)

							fillBuffer(buf,wd, sp, pitch,wd, ht);

						postSegmentProcess(dp , pitch, sp, pitch, wd, ht,
									  buf, tag, ws, wd, d->sm[plane],d-> sh[plane], min, max ); 
			
					}

				}

				else // float
				{
					const float * sp = ( float *)srcp;
						float * dp = (float *)dstp;
					//	float * ws = ( float *)ws8;
						const float ** pixelsort = (const float **) pixelsort8;
						float * buf = ( float *) buf8; 
						const float * s2p = ( float *)s2ptr;
						float min = plane == 0 || fi->colorFamily == cmRGB ? 0.0f : -0.5f;
						float max = plane == 0 || fi->colorFamily == cmRGB ? 1.0f :  0.5f;
					// uses appropriate frame as d->sclip points to that clip frame
						preSegmentProcess(buf, wd, s2p , pitch, wd, ht,  pixelsort);

						WSSegment( pixelsort, buf, tag, dist, ws, ht, wd, d->connect4);

						if ( d->sclip)

							fillBuffer(buf,wd, sp, pitch,wd, ht);

						postSegmentProcess(dp , pitch, sp, pitch, wd, ht,
									  buf, tag, ws, wd, d->sm[plane],d-> sh[plane],min, max ); 
			
					}
			}

			else if( d->sm[plane] == 0 && d->sh[plane] == 0) // this plane not opted for
			{
				vs_bitblt(dstp, dst_stride, srcp, src_stride, wd * fi->bytesPerSample, ht);
			}
		}
	
	//	free(ws);
	//	free(pixelsort8);
	//	free (buf8);
	//	free (dist);
	//	free(tag);
		vsapi->freeFrame(work);
	//	vs_aligned_free(wp);
		vsapi->freeFrame(src);
		if (d->sclip)
			vsapi->freeFrame(s2frame);

		return dst;
	}

	return 0;
}



// Free all allocated data on filter destruction
static void VS_CC amplitudeFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
    AmplitudeData *d = (AmplitudeData *)instanceData;
    vsapi->freeNode(d->node[0]);
	if(d->sclip)
		vsapi->freeNode(d->node[1]);
    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC amplitudeCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) 
{
    AmplitudeData d;
    AmplitudeData *data;
    int err;
	int temp;
    // Get a clip reference from the input arguments. This must be freed later.
    d.node[0] = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi[0] = vsapi->getVideoInfo(d.node[0]);
	if (d.vi[0]->format->colorFamily != cmRGB && d.vi[0]->format->colorFamily != cmYUV && d.vi[0]->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "Amp: RGB, YUV and Gray color formats only for input allowed ");
		vsapi->freeNode(d.node[0]);
		return;
	}
	if (d.vi[0]->format->sampleType == stFloat && d.vi[0]->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "Amp: Half float formats not allowed ");
		vsapi->freeNode(d.node[0]);
		return;
	}

	temp = !!int64ToIntS(vsapi->propGetInt(in, "connect4", 0, &err));

	if(err)
		d.connect4 = true;
	else 
	{
		if ( temp < 0 || temp > 1)
		{
			vsapi->setError(out, "Amp: connect4 must  be 0  to use 8 connect or 1 for  4 connect");
			vsapi->freeNode(d.node[0]);			
			return;
		}

		d.connect4 = temp == 0 ? false : true;
	}
	temp = vsapi->propNumElements(in,"sh");
	if (temp == -1 || temp > 3)
	{
			vsapi->setError(out, "Amp: sh array must specify not more than 3 and at least first of 3 values corresponding to 3 planes");
			vsapi->freeNode(d.node[0]);			
			return;
	}
	else
	{
		
		for (int i = 0; i < temp; i++)
		{
			d.sh[i] = int64ToIntS(vsapi->propGetInt(out, "sh", i, 0));
			
		}
		for (int i = temp; i < 3; i++)
		{
			d.sh[i] = d.sh[i - 1];
		}

		for (int i = 0; i < 3; i++)
		{
			if(d.sh[i] < -5 || d.sh[i] > 5)
			{
				vsapi->setError(out, "Amp: sh values must be between - 5 and 5. If 0 no sharpening will be done");
				vsapi->freeNode(d.node[0]);
				return;
			}

		}
	}

	temp = vsapi->propNumElements(in,"sm");
	if (temp == -1 || temp > 3)
	{
		vsapi->setError(out, "Amp: sm array must specify not more than 3 and at least first of 3 values corresponding to 3 planes");
		vsapi->freeNode(d.node[0]);
		return;
	}
	else
	{
		for (int i = 0; i < temp; i++)
		{
			d.sm[i] = int64ToIntS(vsapi->propGetInt(out, "sm", i, 0));

		}
		for (int i = temp; i < 3; i++)
		{
			d.sm[i] = d.sm[i - 1];
		}

		for (int i = 0; i < 3; i++)
		{
			if (d.sm[i] < 0 || d.sm[i] > 10)
			{
				vsapi->setError(out, "Amp: sm values must be between 0 and 5. If 0 no smoothening will be done");
				vsapi->freeNode(d.node[0]);
				return;
			}

		}
	}

	temp = 0;
	for (int i = 0; i < 3; i++)
	{
		temp += abs(d.sh[i]) + d.sm[0];
	}

	if( temp == 0)
	{
		vsapi->setError(out, "Amp: all sh and sm values are zero so no processing is opted");
		vsapi->freeNode(d.node[0]);			
		return;
	}


	temp = !!int64ToIntS(vsapi->propGetInt(in, "useclip",0,&err));
	if(err)

		d.sclip = false;
	else
	{		

		d.sclip = temp == 0? false : true;
	}

    
	if (d.sclip)
	{
		d.node[1] = vsapi->propGetNode(in, "sclip", 0, &err);

		if(err)
		{
			vsapi->setError(out, "Amp: sclip must be specified for useclip option");
			vsapi->freeNode(d.node[0]);
			return;
		}

		 d.vi[1] = vsapi->getVideoInfo(d.node[1]);

		if (!isConstantFormat(d.vi[0]) || !isSameFormat(d.vi[0], d.vi[1]) || d.vi[0]->numFrames != d.vi[1]->numFrames )
		{
			vsapi->setError(out, "Amp: for use clip option both clips must have constant and identical formats and same number of frames");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode( d.node[1]);
			return;
		}
    }

	

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string, the only
    // reason this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
   

    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = (AmplitudeData *) malloc (sizeof(d));
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
    vsapi->createFilter(in, out, "Amp", amplitudeInit, amplitudeGetFrame, amplitudeFree, fmParallel, 0, data, core);
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
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("in.vcmohan.amp", "vcmod", "VapourSynth Amplitude modification by watershed segmentation", VAPOURSYNTH_API_VERSION, 1, plugin); 
	registerFunc("Amplitude", "clip:clip;useclip:int:opt;sclip:clip:opt;connect4:int:opt;sh:int[];sm:int[];", amplitudeCreate, 0, plugin);
}
*/
/*---------------------------------------------------------------------------------------
vcm is a vapoursynth plugin and has several functions that modify, move pixel values. Also
It has frequency domain functions as well as some utility and test functions
 included functions are:
 a).functions that modify pixel values to reduce noise or introduce blur:-
 amp, variance, hist, gBlur, mBlur, median, saltPepper, veed, fan, neural
 b). functions that move pixels:-
 deBarrel, rotate, reform, deJitter, correctLD
 c) functions of miscellaneous nature;-
 jitter, pattern, grid, bokeh, colorbox
 d) Functions that use Freq Domain :- f1quiver, f1qtest, f2quiver, f2qtest, 
	f2qblur, f2qcorrelation, f2qsharp, f2qbokeh
 

Author V.C.Mohan.
Date 24 Dec 2020, 31 Dec 2020, 31 May 2021, 5 Sep 2021
copyright 2015- 2021

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    see <http://www.gnu.org/licenses/>.
	
-----------------------------------------------------------------------------*/
#include <stdlib.h>
#ifdef _WIN32
#include <Windows.h>
#else
#include <dlfcn.h>
#include <string>
#endif

#include <stdlib.h>
#include <algorithm>
#include <queue>
#include <fstream>
#define _USE_MATH_DEFINES
#include "math.h"

#include <vector>
#include <stack>		
#include <functional>

#include "WSSegment.cpp"
#include "VapourSynth.h"
#include "VSHelper.h"
#include "modHistogramHelper.cpp"
#include <ctime>

#include "UnitSq2Quad_matrix.cpp"
#include "interpolationMethods.h"
#include "colorconverter.h"
#include "statsAndOffsets.h"
# include "ConvertBGRforInput.h"
#include "Squircles.h"
#include "FisheyeMethods.h"
#include "FourFoldSymmetricMarking.h"

#include "fftwlite.h"
//#include "fftw3.h"
#include "fillPlaneWithVal.h"
#define NYQUIST 512
#define ADDSAFE 64
#include "F2QFilters.h"
#include "F2QuiverSpectralDisplay.h"
#include "FQDomainHelper.h"
#include "statsAndOffsets.h"
#include "spotlightDim.h"

#include "colorBox.cpp"
#include "Grid.cpp"
#include "Pattern.cpp"
#include "Jitter.cpp"
#include "Dejitter.cpp"
#include "Circles.cpp"

#include "F1Quiver.cpp"
#include "F2Quiver.cpp"
#include "FQBlur.cpp"
#include "FQSharp.cpp"
#include "FQCorrelation.cpp"
#include "F2QBokeh.cpp"
#include "F2QLimit.cpp"
#include "F1QClean.cpp"

#include "moveRotate.cpp"
#include "moveDeBarrel.cpp"
#include "moveReform.cpp"
#include "Fisheye.cpp"

#include "Mean.cpp"
#include "modFan.cpp"
#include "modAmplitude.cpp"
#include "modHistogram.cpp"
#include "modMedian.cc"

#include "modGBlur.cpp"
#include "modMBlur.cpp"
#include "modVariance.cpp"
#include "modSaltPepper.cpp"
#include "modVeed.cpp"
#include "modNeural.cpp"
#include "modBokeh.cpp"
#include "StepFilter.cpp"


VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("in.vcmohan.cm", "vcm", "VapourSynth Plugin by vcmohan ", VAPOURSYNTH_API_VERSION, 1, plugin);	

	registerFunc("Amp", "clip:clip;useclip:int:opt;sclip:clip:opt;connect4:int:opt;sh:int[];sm:int[];", amplitudeCreate, 0, plugin);
	registerFunc("Fan", "clip:clip;span:int:opt;edge:int:opt;plus:int:opt;minus:int:opt;uv:int:opt;", fanCreate, 0, plugin); 
	registerFunc("Hist", "clip:clip;clipm:clip:opt;type:int:opt;table:int[]:opt;mf:int:opt;window:int:opt;limit:int:opt", histogramadjustCreate, 0, plugin);
	registerFunc("Median", "clip:clip;maxgrid:int:opt;plane:int[]:opt;", adaptivemedianCreate, 0, plugin);
	registerFunc("GBlur", "clip:clip;ksize:int:opt;sd:float:opt;", gblurCreate, 0, plugin);
	registerFunc("MBlur", "clip:clip;type:int:opt;x:int:opt;y:int:opt;", mblurCreate, 0, plugin);
	registerFunc("Neural", "clip:clip;txt:data:opt;fname:data:opt;tclip:clip:opt;xpts:int:opt;ypts:int:opt;tlx:int:opt;tty:int:opt;trx:int:opt;tby:int:opt;iter:int:opt;"
							"bestof:int:opt;wset:int:opt;rgb:int:opt;", neuralCreate, 0, plugin);
	registerFunc("Variance", "clip:clip;lx:int;wd:int;ty:int;ht:int;fn:int:opt;uv:int:opt;xgrid:int:opt;ygrid:int:opt;", varianceCreate, 0, plugin);
	registerFunc("SaltPepper", "clip:clip;planes:int[]:opt;tol:int:opt;avg:int:opt", saltpepperCreate, 0, plugin);
	registerFunc("Veed", "clip:clip;str:int:opt;rad:int:opt;planes:int[]:opt;plimit:int[]:opt;mlimit:int[]:opt;", veedCreate, 0, plugin);
	registerFunc("Mean", "clip:clip;grid:int:opt;tol:float:opt;", meanCreate, 0, plugin);

	registerFunc("F1Quiver", "clip:clip;filter:int[];morph:int:opt;custom:int:opt;test:int:opt;strow:int:opt;nrows:int:opt;gamma:float:opt;", f1quiverCreate, 0, plugin);
	registerFunc("F1QClean", "clip:clip;span:int:opt;fromf:int:opt;upto:int:opt;", f1qcleanCreate, 0, plugin);
	registerFunc("F1QLimit", "clip:clip;span:int:opt;limit:int:opt;freqs:int[];", f1qlimitCreate, 0, plugin);
	registerFunc("F2Quiver", "clip:clip;frad:int:opt;ham:int:opt;test:int:opt;morph:int:opt;gamma:float:opt;fspec:int[];", f2quiverCreate, 0, plugin);
	registerFunc("F2QLimit", "clip:clip;grid:int:opt;inner:int:opt;warn:int:opt;fspec:int[];", f2qlimitCreate, 0, plugin);
	registerFunc("F2QBlur", "clip:clip;line:int:opt;x:int:opt;y:int:opt;", f2qblurCreate, 0, plugin);
	registerFunc("F2QSharp", "clip:clip;line:int:opt;wn:float:opt;x:int:opt;y:int:opt;frad:int:opt;ham:int:opt;;scale:float:opt;rgb:int[]:opt;yuv:int[]:opt", f2qsharpCreate, 0, plugin);
	registerFunc("F2QCorr", "clip:clip;bclip:clip;cx:int:opt;cy:int:opt;txt:int:opt;filename:data:opt;sf:int:opt;ef:int:opt;every:int:opt;", f2qcorrCreate, 0, plugin);
	registerFunc("F2QBokeh", "clip:clip;clipb:clip;grid:int:opt;thresh:float:opt;rgb:int[]:opt;yuv:int[]:opt;", f2qbokehCreate, 0, plugin);

	registerFunc("Rotate", "clip:clip;bkg:clip;angle:float;dinc:float:opt;lx:int:opt;wd:int:opt;ty:int:opt;ht:int:opt;"
							"axx:int:opt;axy:int:opt;intq:int:opt;", rotateCreate, 0, plugin);
	registerFunc("DeBarrel", "clip:clip;abc:float[];method:int:opt;pin:int:opt;q:int:opt;test:int:opt;dots:data:opt;dim:float:opt;", debarrelCreate, 0, plugin);
	
	registerFunc("Reform", "clip:clip;bkg:clip;intq:int:opt;norm:int:opt;rect:float[]:opt;quad:float[];q2r:int:opt;", reformCreate, 0, plugin);
	
	registerFunc("Fisheye", "clip:clip;method:int:opt;xo:int:opt;yo:int:opt;frad:int:opt;sqr:int:opt;"
		"rix:float:opt;fov:float:opt;test:int:opt;dim:float:opt;q:int:opt;dots:int:opt;", fisheyeCreate, 0, plugin);
	

	registerFunc("ColorBox", "format:int:opt;luma:int:opt;nbw:int:opt;nbh:int:opt", colorBoxCreate, 0, plugin);
	registerFunc("Grid", "clip:clip;lineint:int:opt;bold:int:opt;vbold:int:opt;color:int[]:opt;bcolor:int[]:opt;vbcolor:int[]:opt;style:int:opt;", gridCreate, 0, plugin);
	registerFunc("Pattern", "clip:clip;type:int:opt;orient:int:opt;spk:int:opt;spike:float:opt;wl:int:opt;x:int:opt;y:int:opt;rad:int:opt;stat:int:opt;overlay:float:opt;bgr:int[]:opt;", patternCreate, 0, plugin);
	registerFunc("Jitter", "clip:clip;type:int:opt;jmax:int:opt;dense:data:opt;wl:int:opt;stat:int:opt;speed:data:opt;", jitterCreate, 0, plugin);
	registerFunc("DeJitter", "clip:clip;jmax:int:opt;wsyn:int:opt;thresh:float:opt;", dejitterCreate, 0, plugin);
	registerFunc("Bokeh", "clip:clip;clipb:clip;grid:int:opt;thresh:float:opt;rgb:int[]:opt;yuv:int[]:opt;", bokehCreate, 0, plugin);
	registerFunc("StepFilter", "clip:clip;add:int:opt;boost:float:opt;"
		"segmenthor:int:opt;segmentvert:int:opt;limit:int:opt;", stepfilterCreate, 0, plugin);
	registerFunc("Circles", "clip:clip;xo:int:opt;yo:int:opt;frad:int:opt;cint:int:opt;dots:int:opt;rgb:int[]:opt;dim:float:opt;", circlesCreate, 0, plugin);

}
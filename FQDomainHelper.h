#pragma once
#ifndef FREQDOMAIN_HELPER_FUNCTIONS
#define FREQDOMAIN_HELPER_FUNCTIONS

void ApplyFilter(fftwf_complex* fout, float* Filter, int wd, int ht);

/*template <typename finc>
void getRealInput(float* data, const finc* fptr, int pitch,
	int wd, int ht, int wpad, int hpad, bool centered);


template <typename finc>
void getRealOutput(float* data, finc* fptr, int pitch,
	int wd, int ht, int wpad, bool cent, finc min, finc max);*/

	//	int getSign ( int h, int w);

	int getSign(int i);

	int DrawPSF(float* psf, bool linear, int xval, int yval, int bestx, int besty,float spike = 0.0);

	int DrawCircularPSFUV(float* psf, int radius, int bestxUV, int bestyUV, int subX = 0, int subY = 0);

	void DesignInverse(fftwf_complex* fout, float* Filter, float wn,
		int bestx, int besty, float scale);

	void GetFactors(int n, int* facbuf);
	int getBestDim(int dimension, int* facbuf, int nfact = 64);
	// F1Quiver

	void f1BuildFilterCascade(float * FreqFilter, int * filterSpec, int nfft, int npoints);

	void f1BbuildCustomFilter(float * FrqFilt, int * specs, int nfft, int nval);
	template <typename finc>
	void getRowInput(float *data,const finc * rowptr, int nft, int wd);
	template <typename finc>
	void getRowOutput(float *data, finc * rowptr, int wd,finc min, finc max);
	void F1ApplyFilter(fftwf_complex *freqbuf, float * filtbuf, int nfreq);
	template <typename finc>
	void getRowMorphInput(float *data,const finc * rowptr, int nft, int wd,
		bool center = false, int start = 1, float* logLUT = NULL);
	template <typename finc>
	void getRowMorphOutput(float *data, finc * rowptr, int wd,finc min, finc max);
	template <typename finc>
	void f1DisplayHorizontalScale(int nyq, int best, int panelh, int wd, int pitch, finc* dp, finc max);
	template <typename finc>
	finc fclamp(float val, finc min, finc max);

/// F2Quiver
	template <typename finc>
	void getRealInput2D(float* dp, const finc* fptr, int pitch, int ht,
		int wd, int hbest, int wbest, bool centered);
	template <typename finc>
	void getHMRealInput2D(float* data, const finc* fptr, int pitch, int ht, int wd,
							int hbest, int wbest, bool centered, float * logLUT);
	//void getHMRealInput2D8bit(float* dp, const uint8_t* fptr, int pitch, int ht,
	//	int wd, int hbest, int wbest, bool centered, float* logLUT);
	template <typename finc>
	void getRealOutput2D(float* data, finc* fptr, int pitch,
		int ht, int wd, int hbest, int wbest,finc min, finc max);
	template <typename finc>
	void getHMRealOutput2D(float* data, finc* fptr, int pitch, int ht, int wd, int hbest, int wbest,finc min, finc max);	
	
	void ApplyFilter2D(fftwf_complex* out, float* frqFilter, int hbest, int frqwidth);
	
	int getSign(int h, int w);

	float getAmpSquareOfComplex(fftwf_complex* point);
/*
	// FQCorr
	template <typename finc>
	void xFillPlaneWithVal(finc* wp, const int wpitch, const int wd, const int ht, finc val);
	template <typename finc>
	// convert input unsigned char data to  float type
	void xGetRealInput(float* data, const finc* fptr,
		int pitch, int wd, int ht,
		int wpad, int hpad, bool centered);

	void xCorrelate(fftwf_complex* Afreq, fftwf_complex* Bfreq, int fsize);
	int xNormGamma(float* buf, int size);
*/
	template <typename finc>
	void xTransferToDst(finc* dp, int dpitch,
		int dwd, int dht, float* buf, int bestwd, int bestht, finc max);
		// FQSharp and F2Quiver
	void F2QhammingWindowing(float* cosBell, int pitch, int width, int height, int rfilt);
		//F2QLimit
	void removeInputCentering(float* inBuf, int wbest, int hbest);
//--------------------------------------------------------------------------------------------
void GetFactors(int n, int* facbuf)
{
	// finds factors and fills facbuf with factor, remainder. maximum allowed or 32 pairs
	// factors are 4,2,3,5,....  remaining primes. 
	int p = 4;
	double floor_sqrt;
	floor_sqrt = floor(sqrt((double)n));

	/*factor out powers of 4, powers of 2, then any remaining primes */
	do
	{
		while (n % p)
		{
			switch (p)
			{
			case 4: p = 2; break;
			case 2: p = 3; break;
			default: p += 2; break;
			}
			if (p > floor_sqrt)
				p = n;          /* no more factors, skip to end */
		}
		n /= p;
		*facbuf++ = p;
		*facbuf++ = n;
	} while (n > 1);
}

int getBestDim(int dimension, int* facbuf, int nfact)
{
	// returns nearest larger value having factors limited  2, 3, and 5
	int n = dimension;
	int largest = 7;
	size_t i;

	while (largest > 5)
	{
		GetFactors(n, facbuf);
		// check value of largest factor
		for (i = 0; i < nfact; i += 2)
			if (facbuf[i + 1] == 1)	// 
			{
				largest = facbuf[i];
				break;
			}

		if (largest > 5)
			n += 4;
	}
	return (n);
}

template <typename finc>
finc fclamp(float val, finc min, finc max)
{
	return (finc)(val < min ? min : val > max ? max : val);
}

float getAmpSquareOfComplex(fftwf_complex* point)
{
	return  (point[0][0] * point[0][0] + point[0][1] * point[0][1]);

}
//-----------------------------------------------------------------------------
// one dimension functions
void f1BuildCustomFilter(float* FrqFilt, int* specs, int nfft, int nval)
{
	// extend towards zero
	if (specs[0] > 0)
	{
		int f = (specs[0] * nfft / 2) / NYQUIST;
		for (int j = 0; j < f; j++)	// extend first value towards zero
		{

			FrqFilt[j] = 0.01f * specs[1];	// it is a %age value. hence 0.01
		}
	}
	// extend till end
	if (specs[nval - 2] < NYQUIST)
	{
		int f = (specs[nval - 2] * nfft / 2) / NYQUIST;
		for (int j = f; j < nfft / 2; j++)	// extend last value to end

			FrqFilt[j] = 0.01f * specs[nval - 1];
	}

	for (int i = 0; i < nval - 2; i += 2)
	{
		int f1 = (specs[i] * nfft / 2) / NYQUIST;
		int f2 = (specs[i + 2] * nfft / 2) / NYQUIST;
		//  linear interpolation
		for (int j = f1; j < f2; j++)
		{
			FrqFilt[j] = 0.01f * (specs[i + 1] + ((j - f1) * (specs[i + 3] - specs[i + 1])
				/ (f2 - f1)));
		}

	}

	FrqFilt[nfft / 2] = FrqFilt[nfft / 2 - 1];	// nyquist freq bin

		//  apply 5 point averaging 	
	for (int i = 2; i < nfft / 2 - 2; i++)
	{
		float frq = 0.0;
		for (int k = -2; k < 3; k++)
		{
			frq += FrqFilt[i + k];
		}

		FrqFilt[i] = frq / 5;
	}
}
//---------------------------------------------------------------------------------

void f1BuildFilterCascade(float* FreqFilter, int* filterSpec, int nfft, int npoints)
{
	// build cascaded butterworth filter.
	// array first val:- 0: reduction 1 highcut, 2 lowcut, 3 band pass, 4  band stop
	// array second val :- frequency
	// array 3rd val :- band width (applicable for band pass only)
	// array 4th val :- degree of sharpness. max 12. Note degree = 1 is Gaussian
	for (int i = 0; i < npoints; i += 4)
	{
		int type = filterSpec[i];

		float freq = filterSpec[i + 1] * nfft / (2 * NYQUIST);

		float bandwidth = (freq * filterSpec[i + 2]) / 100.0f;		

		int degree = 2 * filterSpec[i + 3];

		if (type == 0)
		{			
			float freq2 = filterSpec[i + 2] * nfft / (2 * NYQUIST);

			for (int i = freq; i <= freq2; i++)
			{
				FreqFilter[i] *= 1.0 / (1 + degree);
			}
		}

		else if (type == 1)
		{
			// high frequencies are filtered out

			for (int j = 1; j < nfft / 2; j++)
			{
				FreqFilter[j] *= 1.0f / (1.0f + pow(j / freq, degree));
			}
		}

		else if (type == 2)
		{
			// low frequencies filtered off
			for (int j = 1; j < nfft / 2; j++)
			{
				FreqFilter[j] *= 1.0f / (1.0f + pow(freq / j, degree));
			}

			FreqFilter[0] = 0.0;
		}

		else if (type == 3)
		{
			// band pass filter
			for (int j = 1; j < nfft / 2; j++)
			{
				FreqFilter[j] *= 1.0f / ((1.0f + pow(j / (freq + bandwidth), degree)) * (1.0f + pow((freq - bandwidth) / j, degree)));


			}
		}

		else if (type == 4)
		{
			// notch/ band reject filter
			for (int j = 1; j < nfft / 2; j++)
			{
				FreqFilter[j] *= float(1.0 - 1.0 / ((1.0 + pow(j / (freq + bandwidth), degree)) * (1.0 + pow((freq - bandwidth) / j, degree))));

			}
		}

		FreqFilter[nfft / 2] = FreqFilter[nfft / 2 - 1];	// nyquist freq
	}
}
//------------------------------------------------------------------------
void F1ApplyFilter(fftwf_complex* freqbuf, float* filter, int nfreq)
{
	for (int i = 0; i < nfreq; i++)
	{
		freqbuf[i][0] *= filter[i];

		freqbuf[i][1] *= filter[i];
	}
}
//-----------------------------------------------------------------------------
template <typename finc>

void getRowInput(float* data, const finc* rowptr, int nft, int wd)
{
	// for apparent vertical noise, row data is needed 

	for (int i = 0; i < wd; i++)

		data[i] = rowptr[i];

	// zero padding 
	for (int i = wd; i < nft; i++)

		data[i] = 0.0;
}
//------------------------------------------------------------------------------
template <typename finc>
void getRowMorphInput(float* data, const finc* rowptr, int nft, int wd, bool center, int start, float * logLUT)
{
	//  float * logLUT will be default null for float and more than 12 bit value input
	// start value default 1 or minus 1 for centering transformed spectrum
	// centered default value false. whether centering required
	if (center)
	{
		if (logLUT == NULL)
		{
			for (int i = 0; i < wd; i++)
			{

				data[i] = start * log((float)rowptr[i]);
				start = -start;
			}
		}
		else
		{
			for (int i = 0; i < wd; i++)
			{
				data[i] = start * logLUT[(int)rowptr[i]];
				start = -start;
			}
		}
	}

	else // if (!center)
	{
		if (logLUT == NULL)
		{
			for (int i = 0; i < wd; i++)
			{
				data[i] =  log((float)rowptr[i]);				
			}
		}
		else
		{
			for (int i = 0; i < wd; i++)
			{
				data[i] = logLUT[(int)rowptr[i]];				
			}
		}
	}


	for (int i = wd; i < nft; i++)

		data[i] = 0.0;
}
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
template <typename finc>
void getRowOutput(float* data, finc* rowptr, int wd,finc min, finc max)
{
	
	for (int i = 0; i < wd; i++)
	{
		float val = data[i];	// scale down

		rowptr[i] = fclamp(val, min, max);


	}

}
//---------------------------------------------------------------------------
template <typename finc>
void getRowMorphOutput(float* data, finc* rowptr, int wd,finc min, finc max)
{
	
	for (int i = 0; i < wd; i++)
	{
		float expo = exp(data[i]);

		rowptr[i] = fclamp(expo, min, max);
	}
}
//-------------------------------------------------------------------------------------------------------------------------
template <typename finc>
void f1DisplayHorizontalScale(int nyq, int best, int panelh, int wd, int pitch, finc* dp, finc max)
{
	for (int i = 0; i < nyq; i++)
	{
		int ws = (i * best) / (2 * nyq);
		// display horizontal scale for freq
		if ((i % 100) == 0)
		{
			for (int h = 0; h < 10; h++)
			{

				dp[(h + panelh + 5) * pitch + ws] = max;
			}
		}

		else if ((i % 50) == 0)
		{
			if (wd >= nyq / 20)	// to ensure readability
			{
				for (int h = 0; h < 6; h++)
				{

					dp[(h + panelh + 6) * pitch + ws] = (4 * max) / 5;
				}
			}
		}
		else if ((i % 10) == 0)
		{
			if (wd >= nyq / 2)
			{
				for (int h = 0; h < 3; h++)
				{

					dp[(h + panelh + 7) * pitch + ws] = (3 * max) / 4;
				}
			}
		}

	}

}
//-------------------------------------------------------------------------------------------
// Two demensional fft functions
//---------------------------------------------------------------------------------------------
			//Draws the PSF linear or circular at center of the psf buffer 
		// draws psf in float data
int  DrawPSF(float* psf, bool linear, int xval, int yval, int bestx, int besty, float spike)
{

	int count;

	// zero psf area
	for (int h = 0; h < bestx * besty; h++)

		psf[h] = 0.0;

	if (linear)
	{
		int length = abs(yval) > xval ? abs(yval) : xval;

		count = 2 * length + 1;

		// draw the blur line at the center of frames best sizes

		if (abs(yval) > xval)
		{
			if (yval < 0)
			{
				yval = -yval;

				xval = -xval;
			}


			for (int h = -abs(yval); h <= abs(yval); h++)
			{

				int w = (h * xval) / yval;		// get nearest integer

				int fraction = abs(h * xval) % abs(yval);		// get fractional part

					// distribute amp in ratio of fraction

				psf[(besty / 2 + h) * bestx + (bestx / 2 + w)]
					= (1.0f * (abs(yval) - fraction)) / (count * abs(yval));

				if (h * xval > 0)

					psf[(besty / 2 + h) * bestx + (bestx / 2 + w + 1)]
					= (1.0f * fraction) / (count * abs(yval));

				else	// if( h  < 0)

					psf[(besty / 2 + h) * bestx + (bestx / 2 + w - 1)] = (1.0f * fraction) / (count * abs(yval));



			}
		}


		else	// xval is greater than yval

		{
			//xval is always a positive


			for (int w = -xval; w <= xval; w++)
			{

				int h = (w * yval) / xval;

				int fraction = abs((w * yval)) % xval;		// get fractional part

				psf[(besty / 2 + h) * bestx + (bestx / 2 + w)]
					= (1.0f * (xval - fraction)) / (count * xval);

				if (w * yval > 0)

					psf[(besty / 2 + h + 1) * bestx + (bestx / 2 + w)]
					= (1.0f * fraction) / (count * xval);

				else	// if( w * yval < 0)

					psf[(besty / 2 + h - 1) * bestx + (bestx / 2 + w)]
					= (1.0f * fraction) / (count * xval);

			}

		}

	}

	else		// circular

	{

		// draw the blur circle at the center of frame

		count = 0;

		for (int h = -xval; h <= xval; h++)

			for (int w = -xval; w <= xval; w++)

				if (h * h + w * w <= xval * xval)

					count++;

		// incase of deblur there is a spike added to center value
		// actually all values except center are reduced.

		for (int h = -xval; h <= xval; h++)

			for (int w = -xval; w <= xval; w++)

				if (h * h + w * w <= xval * xval)

					psf[(besty / 2 + h) * bestx + (bestx / 2 + w)] = 1.0f / count;

		// center value will have spike to provide white noise and reduce instability
		if (spike > 0.001)

			psf[(besty / 2) * bestx + (bestx / 2)] = (1.0f + spike) / count;
	}

	return count;
}
//----------------------------------------------------------------------------------------------------------
int DrawCircularPSFUV(float* psf, int rad, int bestx, int besty, int subX, int subY)
{
	//zero out buffer
	for (int h = 0; h < (bestx) * (besty); h++)

		psf[h] = 0.0;

	int count = 0;
	int andX = (1 << subX) - 1;
	int andY = (1 << subY) - 1;
	int xval = rad >> subX;
	int yval = rad >> subY;
	int rsq = rad * rad;

	for (int h = -yval; h <= yval; h++)
	{
		int hsq = (h << subY) * (h << subY);

		for (int w = -xval; w <= xval; w++)
		{
			int wsq = (w << subX) * (w << subX);

			if (hsq + wsq <= rsq)
			{
				psf[(besty / 2 + h) * bestx + (bestx / 2 + w)] = 1.0f;

				count++;
			}
		}
	}

	for (int i = 0; i < bestx * besty; i++)

		psf[i] /= count;

	return count;
}

//------------------------------------------------------------------------------------------------------------
void DesignInverse(fftwf_complex* fout, float* Filter,
	float wn, int bestx, int besty,
	float scale)
{
	// the forward transform of PSF is in fout. Only real  positive values
				// get max value
	float mval = fout[0][0];


	for (int h = 0; h < bestx * besty; h++)


		// get max value. Possibly the zeroth val is max but lets find
		if ((fout[h][0]) > mval)

			mval = fout[h][0];

	mval *= wn;		// this is the min value we will accept for inversion

					//  apply freq mask, and wn of inversion
	// PSF during forward transform is scaled up by sqrt(bestx * besty)
	// while inverting it it is 1/ this value. so no sscaling for transform is needed.
	float scaler = scale * (1.0 + wn) / (bestx * besty) / (1.0f - wn); // This trial/error derived approximation for scaler
									// as the higher wn is frequencies are inverted less
									// and we lose amplitude
					// we add white noise and  scaler . 
	for (int h = 0; h < bestx * besty; h++)


		Filter[h] = scaler / (fout[h][0] + mval);



}

//------------------------------------------------------------------------------------------------------------
// F2Quiver uses these
//---------------------------------------------------------------------------------------------
	// use this when data type is  complex 
int getSign(int i)
{
	return (i & 1) == 0 ? 1 : -1;

}

int getSign(int h, int w)
{

	return ((h + (w >> 1)) & 1) == 0 ? 1 : -1;
}


//---------------------------------------------------------------------------------------------
template <typename finc>
void getRealInput2D(float* in, const finc* ptr, int pitch, int ht,
	int wd, int hbest, int wbest, bool centered)
{
	// convert frame y values to float and keep in data buffer
	float* data = in;
	if (centered)
	{
		// values multiplied by -1^(x+y) i.e sign to get the spectogram centered in frame
		for (size_t h = 0; h < ht; h++)
		{
			for (size_t w = 0; w < wd; w++)
			{

				data[w] = getSign(h, w) * ptr[w];
			}
			data += wbest;
			ptr += pitch;
		}
	}


	else		// not centered. So as it is 	
	{
		// values multiplied by -1^(x+y) i.e sign to get the spectogram centered in frame
		for (size_t h = 0; h < ht; h++)
		{
			for (size_t w = 0; w < wd; w++)
			{

				data[w] = ptr[w];
			}
			data += wbest;
			ptr += pitch;
		}
	}

	// fill with zeroes rest of buffer
	data = in + ht * wbest;
	for (size_t h = ht; h < hbest; h++)
	{
		for (size_t w = 0; w < wbest; w++)
		{

			data[w] = 0.0;
		}

		data += wbest;
	}

	//	right margin

	for (size_t w = wd; w < wbest; w++)
	{
		data = in + w;

		for (size_t h = 0; h < hbest; h++)
		{

			data[0] = 0.0;

			data += wbest;

		}
	}

}
//----------------------------------------------------------------------------
template <typename finc>
void getHMRealInput2D(float* in, const finc* ptr, int pitch, int ht,
	int wd, int hbest, int wbest, bool centered, float * logLUT)
{
	int start = 1;
	// convert frame y values to float and keep in data buffer
	float* data = in;
	if (centered)
	{
		if (logLUT == NULL)
		{
			// values multiplied by -1^(x+y) i.e sign to get the spectogram centered in frame
			for (size_t h = 0; h < ht; h++)
			{
				for (size_t w = 0; w < wd; w++)
				{

					data[w] = start * log(2.0 + ptr[w]);
					start = -start;
				}
				data += wbest;
				ptr += pitch;
			}
		}

		else //  logLUT is available
		{
			// values multiplied by -1^(x+y) i.e sign to get the spectogram centered in frame
			for (size_t h = 0; h < ht; h++)
			{
				for (size_t w = 0; w < wd; w++)
				{

					data[w] = start * logLUT[(int) ptr[w]];
					start = -start;
				}
				data += wbest;
				ptr += pitch;
			}
		}
	}


	else		// not centered. So as it is 	
	{
		// values multiplied by -1^(x+y) i.e sign to get the spectogram centered in frame
		for (size_t h = 0; h < ht; h++)
		{
			for (size_t w = 0; w < wd; w++)
			{
				if ( logLUT == NULL)
					data[w] = log(2.0 + ptr[w]);
				else
					data[w] = logLUT[ (int)ptr[w]];
			}
			data += wbest;
			ptr += pitch;
		}
	}

	// fill with zeroes rest of buffer
	for (size_t h = ht; h < hbest; h++)
	{
		for (size_t w = 0; w < wbest; w++)
		{

			data[w] = 0.0;
		}

		data += wbest;
	}

	//	right margin

	for (size_t w = wd; w < wbest; w++)
	{
		data = in + w;

		for (size_t h = 0; h < hbest; h++)
		{

			data[0] = 0.0;

			data += wbest;

		}
	}

}

//===============================================================================
template <typename finc>
void getRealOutput2D(float* in, finc* ptr, int pitch, int ht,
	int wd, int hbest, int wbest,finc min,  finc max)
{
	for (int h = 0; h < ht; h++)
	{
		for (int w = 0; w < wd; w++)
		{
			ptr[w] = fclamp(in[w], min, max);
		}
		ptr += pitch;
		in += wbest;
	}
}

//===============================================================================
template <typename finc>
void getHMRealOutput2D(float* in, finc* dp, int pitch, int ht,
	int wd, int hbest, int wbest,finc min, finc max)
{
	for (int h = 0; h < ht; h++)
	{
		for (int w = 0; w < wd; w++)
		{
			float val = exp(in[w]) - 2;

			dp[w] = fclamp(val, min, max);
		}
		dp += pitch;
		in += wbest;
	}
}
//----------------------------------------------------------------------------------------

void ApplyFilter2D(fftwf_complex* out, float* filter, int hbest, int frqwd)
{
	// applies the designed filter.in freq domain just multiplication
	// of freq response of designed filter with freq transform of input
	int nval = hbest * frqwd;

	for (int i = 0; i < nval; i++)
	{
		out[i][0] *= filter[i];
		out[i][1] *= filter[i];
	}

}

//-------------------------------------------------------------------------------------------------------


void ApplyFilter(fftwf_complex* fout, float* Filter, int wd, int ht)
{
	// applies the designed filter.in freq domain. Scalar multiply
	// of freq response of designed filter with freq transform of input

//float scale = 1.0 / (hbest * wbest );
	int nval = ht * wd;
	for (int h = 0; h < nval; h++)

	{
		fout[h][0] *= Filter[h];
		fout[h][1] *= Filter[h];

	}

}
//-------------------------------------------------------------------------------------------------------------------------
void F2QhammingWindowing(float* cosBell, int pitch, int width, int height, int rfilt)
{

	// design a cosine bell hamming windowing function
	int rfiltsq = rfilt * rfilt;
	//	float rfilt = sqrt(rfiltsq );	

	for (int h = 0; h < height / 2; h++)
	{
		for (int w = 0; w < width / 2; w++)
		{
			if (h * h + w * w <= rfiltsq)
			{
				float radial = sqrt((float)h * h + w * w);
				float bell = 0.46 + 0.54 * cos((M_PI * radial) / rfilt);

				for (int hh = -1; hh <= 1; hh += 2)
				{
					for (int ww = -1; ww <= 1; ww += 2)
					{

						cosBell[(height / 2 + h * hh) * pitch + width / 2 + w * ww] *= bell;
					}
				}
			}

			else
			{
				for (int hh = -1; hh <= 1; hh += 2)
				{
					for (int ww = -1; ww <= 1; ww += 2)
					{

						cosBell[(height / 2 + h * hh) * pitch + width / 2 + w * ww] *= 0.0;
					}
				}
			}
		}
	}
}
//............................................................
void removeInputCentering(float* inBuf, int wbest, int hbest)
{	
	float scale = 1.001 / (hbest * wbest);
	for (int h = 0; h < hbest; h++)
	{
		for (int w = 0; w < wbest; w++)
		{
			inBuf[w] *= scale * getSign(h, w);
		}
		inBuf += wbest;
	}
}
//......................................................................
#endif


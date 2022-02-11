#ifndef INTERPOLATION_METHODS_V_C_MOHAN
#define INTERPOLATION_METHODS_V_C_MOHAN
//-------------------------------------------------------------------
void LanczosCoeff(float* cbuf, const int span, const int quantp);

float sinc(float f);

void CubicIntCoeff(float* cbuf, const int quantiles);

void LinearIntCoeff(float* cbuf, const int quantiles);
float * setInterpolationScheme(const int q, const int quantile, int *dspan);

//-------------------------------------------------------------------
template <typename finc>
float LaQuantile(const finc* point, const int spitch, const int kb,
	const int span, const int qx, const int qy, const float* lbuf);

template <typename finc>
bool needNotInterpolate(const finc* sp, int const pitch, const int kb);

template <typename finc>
float alongLineInterpolate(const finc* point, const int step,
	const int span, const int quant, const float* iBuf);

int bestOfNineIndex(float fx, float fy);
int bestOfNineIndex(int qx, int qy, int quantile);
template <typename finc>
finc bestOfNine(const finc* fp, const int fpitch, const int kb,
	const int x, const int y, const int index);

//restricts values to min and max
//---------------------------------------------------------------------------------
template <typename  finc>
finc clamp(const float val, const finc min, const finc max);
//--------------------------------------------------------------------------------
#ifndef CLAMP_V_C_MOHAN
#define CLAMP_V_C_MOHAN
template <typename  finc>
finc clamp(const float val, const finc min, const finc max)
{

	return (finc)(val < min ? min : val > max ? max : val);
}
#endif
//----------------------------------------------------------------------------------
float fclamp(const float val, const float min, const float max)
{
	return val < min ? min : val > max ? max : val;
}

//-------------------------------------------------------------------

template <typename finc>
// this is general. can be used for any span 2 4 6 and bytes per pixel planar formats

float LaQuantile(const finc* point, const int spitch, const int kb,
	const int span, const int qx, const int qy, const float* lbuf)
{
	// nb = bit depth of sample. for float 0
	// span = 6 for 6 x 6 point interpolation
	// point is input nearest (floor) pixel
	// qx, qy are quantiles in x and y 
	// lbuf has precomputed coefficients for quantiles
	if (span == 0)
	{
		// near point
		return (float)(*point);
	}
	float xy[6];	//
	point += (-span / 2 + 1) * spitch;

	const float* lbufr = lbuf + span * qx;

	for (int h = 0; h < span; h++)
	{
		xy[h] = 0;

		for (int w = 0; w < span; w++)
		{
			xy[h] += point[(w - span / 2 + 1) * kb] * lbufr[w];
		}

		point += spitch;
	}

	float sum = 0;

	lbufr = lbuf + span * qy;

	for (int h = 0; h < span; h++)
	{
		sum += xy[h] * lbufr[h];
	}

	return sum;

}
//-------------------------------------------------------------------------------------------------	

template <typename finc>
bool needNotInterpolate(const finc* sp, const int pitch, const int kb)
{
	return (*sp == *(sp + kb) && *sp == *(sp + pitch) && *sp == *(sp + pitch + kb));
}
//-------------------------------------------------

void LinearIntCoeff(float* cbuf, const int quantiles)
{
	float q = 0, qinc = 1.0f / quantiles;

	for (int i = 0; i < 2 * quantiles + 2; i += 2)
	{
		cbuf[i] = 1.0f - q;
		cbuf[i + 1] = q;

		q += qinc;
	}
}
//------------------------------------------------
template <typename finc>
float alongLineInterpolate(const finc* point, const int step,
	const int span, const int quant, const float* iBuf)
{
	float sum = 0.0;

	for (int i = 0; i < span; i++)
	{
		sum += iBuf[span * quant + i] * point[(i + 1 - span / 2) * step];
	}
	return sum;
}
//...................................................................

/* For cubic interpolation generates coefficients for each quantile
 between 0.0 and 3.0
 cbuf should be 4 * quantiles per unit interval  * 3 plus 4
 given a quantile from a0,, multiply values at a0, a1, a2, a3 to get interpolated value
*/
/*
void CubicIntCoeff( int * cbuf, int quantiles,  int prec)
{
	float inc = 1.0 / quantiles;

	float startx = 0;

	for( int i = 0, i <= 3 * quantiles; i ++)
	{
		float x = startx;
		float xsq = x * x;
		float xcube = x * xsq;

		cbuf[4 * i    ] = ( 1.0 - 29 * x /18 + 4 * xsq / 6 - xcube / 18) * prec;	// a0 multiplier

		cbuf[4 * i + 1] = ( 42 * x / 18 - 9 * xsq / 6 + xcube / 6) * prec;			// a1 multiplier

		cbuf[4 * i + 2] = ( - 15* x / 18 + xsq - xcube / 6) * prec;					// a2 multiplier

		cbuf[4 * i + 3] = ( x / 9 - xsq / 6 + xcube / 18) * prec;					// a3 multiplier

		startx += inc;
	}
}

*/

// Alternate formulation
/*
 given p0, p1, p2 p3
 at   x-1, x0, x1, x2
 for any given x ( between x0 and x1)
 f( p0,p1,p2,p3,x) = ( -1/2 p0 + 3/2 p1 - 3/2 p2 + 1/2 p3) * x*x*x
					+ ( p0 -5/2 p1 + 2 p2 - 1/2 p3) * x *x
					+(- 1/2 p0 + 1/2 p2) * x
					+ p1
Or as coefficients at discrete quantiles of x
	( -x*x*x + 2*x*x - x) /2  for p0
	(3*x*x* - 5*x*x + 2) /2 for p1
	(-3*x*x*x + 4*x*x + x) /2 for p2
	( x*X*X - x*x) / 2 for p3

cbuf should be 4 * quantiles  plus 4
 given a quant from x0,, multiply values at
 cbuf[4 * quant} with p0, p1, p2 and p3
 at x-1, x0, x1, x2 add and divide by prec to get interpolated value
 at that quant
 prec (between 1 << (1 to  23) enables integer arithmetic avoiding float
 final result is left shifted same amount
*/

void CubicIntCoeff(float* cbuf, const int quantiles)
{


	float inc = 1.0f / quantiles;

	float x = 0;

	for (int i = 0; i <= 4 * quantiles; i += 4)
	{

		float xsq = x * x;

		float xcube = x * xsq;

		cbuf[i] = ((-x + 2 * xsq - xcube));	// p0 multiplier

		cbuf[i + 1] = ((2 - 5 * xsq + 3 * xcube));	// p1 multiplier

		cbuf[i + 2] = ((x + 4 * xsq - 3 * xcube));// p2 multiplier

		cbuf[i + 3] = ((-xsq + xcube));			// p3  multiplier

		float sum = cbuf[i] + cbuf[i + 1] + cbuf[i + 2] + cbuf[i + 3];
		//	normalize

		cbuf[i] = (cbuf[i]) / sum;

		cbuf[i + 1] = (cbuf[i + 1]) / sum;

		cbuf[i + 2] = (cbuf[i + 2]) / sum;

		cbuf[i + 3] = (cbuf[i + 3]) / sum;


		x += inc;
	}
}

//-------------------------------------------------------------------

void LanczosCoeff(float* cbuf, const int span, const int quantp)
{
	/*
	Fills buffer with coefficients
	cbuf: in which coefficients to be placed. Must be (quantp + 1) span
	span:  lanczos  number of values used for one interpolation 4 or 6
	quantp: quantile precision interval for interpolation

	*/
	for (int s = 0; s < span; s++)
	{
		// zero out first and last quantile ( zero value quantile)

		cbuf[s] = 0;
		cbuf[quantp * span + s] = 0;
	}

	cbuf[span / 2 - 1] = 1.0;	// so that the nearest pixel value is used 

	cbuf[quantp * span + span / 2] = 1.0f;

	float fraction = 1.0f / quantp;

	float frac = fraction;

	for (int s = span; s < span * quantp; s += span)
	{
		float csum = 0.0f;		// sum of coefficients should equal 1.0 

		for (int i = 0; i < span; i++)
		{
			// calculte lanczos coefficients
			cbuf[s + i] = sinc(span / 2 - 1 - i + frac) * sinc((span / 2 - 1 - i + frac) / (span / 2));

			csum += cbuf[s + i];
		}
		for (int i = 0; i < span; i++)
		{
			cbuf[s + i] *= (1.0f / csum);
		}
		frac += fraction;
	}
}

//-----------------------------------------------------------------------------------

float sinc(float f)
{
	if (f < 0.0)

		f = -f;

	f *= (float)M_PI;	// pi

	return (f != 0) ? (float)sin(f) / f : 1.0f; ;
}
//......................................................................
int bestOfNineIndex(float fx, float fy)
{
	int index = 0;
	if (fx < 0.25f)
		index = 0;
	else if (fx < 0.75f)
		index = 1;
	else
		index = 2;

	if (fy < 0.25f)
		index += 0;
	else if (fy < 0.75f)
		index += 3;
	else
		index += 6;

	return index;
}


int bestOfNineIndex(int qx, int qy, int quantile)
{
	int index = 0;
	if (qx < quantile / 4)
		index = 0;
	else if (qx < (3 * quantile) / 4)
		index = 1;
	else
		index = 2;

	if (qy < quantile / 4)
		index += 0;
	else if (qy < (3 * quantile) / 4)
		index += 3;
	else
		index += 6;

	return index;
}
template <typename finc>
finc bestOfNine(const finc* fp, const int fpitch, const int kb,
	const int x, const int y, const int index)
{

	finc xnear = *(fp + y * fpitch + x * kb);
	finc ynear = xnear;
	finc xnext = *(fp + y * fpitch + (x + 1) * kb);
	finc ynext = *(fp + (y + 1) * fpitch + x * kb);
	finc xynext = *(fp + (y + 1) * fpitch + (x + 1) * kb);

	switch (index)
	{
	case 0: return xnear;
	case 1: return (xnear + xnext) / 2;
	case 2: return xnext;
	case 3: return ynear;
	case 4: return (xnear + xnext + ynear + ynext) / 4;
	case 5: return (xnext + xynext) / 2;
	case 6: return ynext;
	case 7: return (ynext + xynext) / 2;
	case 8: return xynext;
	}
	return 0;
}

//----------------------------------------------------------------------------------------------	
float * setInterpolationScheme(const int q, const int quantile, int *dspan)
{
	// presets.		
	float* iCoeff = NULL;
	int span;
	if (q == 1)
	{
		span = 2;	// ninept 
		iCoeff = NULL;
	}
	else if (q == 2)
	{
		span = 2;
		iCoeff = new float[span * (quantile + 1)];
		LinearIntCoeff(iCoeff, quantile);
	}
	else if (q == 3)
	{
		span = 4; // cubic
		iCoeff = new float[span * (quantile + 1)];
		CubicIntCoeff(iCoeff, quantile);
	}
	else if (q == 4)
	{
		span = 6; // lanczos is 6 x 6
		iCoeff = new float[span * (quantile + 1)];	// 1/64 accuracy 6 coeff per quantile
		// create lanczos coefficients for every quantile
		LanczosCoeff(iCoeff, span, quantile);
	}
	*dspan = span;
	return iCoeff;
}
#endif // !INTERPOLATION_STUFF_V_C_MOHAN


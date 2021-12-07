#ifndef FOURFOLD_SYMMETRIC_MARKING_H_V_C_MOHAN
#define FOURFOLD_SYMMETRIC_MARKING_H_V_C_MOHAN
// using four fold symmetry output is marked
// requires interpolationmethods.h for some functions
// ensure no access violation occurs
// dp and sp points to the centers, w and h are x and y distances from center
// span 2 is bilinear, 4 is bicubic and 6 is lanczos 6x6 interpolation
// quantile is quantized accuracy (1/q) of interpolation. LaQuantile is an interpolater.
//iCoeff has (quantiles + 1)* span precomputed coefficients.
template <typename finc>
void paint4FoldSym(finc* dp, const int dpitch,const int kb, const int dw, const int dh, const finc col);
template <typename finc>
void paint4FoldSymIfChecks(finc* dp, const int dpitch, const int kb,
	const int width, const int height, const int centerx, const int centery,
	const int dw, const int dh, const int sw, const int sh, const finc col);
template <typename finc>
void copy4FoldSym(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int dw, const int dh, const int sw, const int sh);
template <typename finc>
void copy4FoldSymIfChecks(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int width, const int height, const int centerx, const int centery,
	const int dw, const int dh, const int sw, const int sh);
template <typename finc>
void interpolate4FoldSym(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int dw, const int dh, const int sw, const int sh,
	const int qx, const int qy, const int span,
	const float* iCoeff, finc min, finc max);
template <typename finc>
void interpolate4FoldSymIfChecks(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int width, const int height, const int centerx, const int centery,
	const int dw, const int dh, const int sw, const int sh,
	const int qx, const int qy, const int span,
	const float* iCoeff, finc min, finc max);
template <typename finc>
void interpolate9pt4FoldSym(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int dw, const int dh, const int sw, const int sh,
	const int index);
template <typename finc>
void interpolate9pt4FoldSymIfChecks(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int width, const int height, const int centerx, const int centery,
	const int dw, const int dh, const int sw, const int sh,
	const int index);
//.......................................................................................
template <typename finc>
void paint4FoldSym(finc* dp, const int dpitch, const int kb, const int dw, const int dh, const finc col)
{
	// dp points to center of circle / square
	*(dp + dh * dpitch + kb * dw) = col;
	*(dp + dh * dpitch - kb * dw) = col;
	*(dp - dh * dpitch + kb * dw) = col;
	*(dp - dh * dpitch - kb * dw) = col;
}

template <typename finc>
void copy4FoldSym(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int dw, const int dh, const int sw, const int sh)
{
	// dp and sp points to centers of output and input circles
	*(dp + dh * dpitch + kb * dw) = *(sp + sh * spitch + kb * sw);
	*(dp + dh * dpitch - kb * dw) = *(sp + sh * spitch - kb * sw);
	*(dp - dh * dpitch + kb * dw) = *(sp - sh * spitch + kb * sw);
	*(dp - dh * dpitch - kb * dw) = *(sp - sh * spitch - kb * sw);

}
template <typename finc>
void interpolate4FoldSym(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int dw, const int dh, const int sw, const int sh,
	const int qx, const int qy,  const int span,
	const float* iCoeff, finc min, finc max)
{
	// dp and sp at centers of output and input circles;
	*(dp + dh * dpitch + kb * dw) =
		clamp(LaQuantile(sp + sh * spitch + kb * sw, spitch,kb,
			span, qx, qy, iCoeff), min, max);
	*(dp + dh * dpitch - kb * dw) =
		clamp(LaQuantile(sp + sh * spitch  - kb * sw, spitch, -kb,
			span, qx, qy, iCoeff), min, max);
	*(dp - dh * dpitch + kb * dw) =
		clamp(LaQuantile(sp - sh * spitch + kb * sw, -spitch, kb,
			span, qx, qy, iCoeff), min, max);
	*(dp - dh * dpitch - kb * dw) =
		clamp(LaQuantile(sp - sh * spitch - kb * sw, -spitch, -kb,
			span, qx, qy, iCoeff), min, max);

}

template <typename finc>
void interpolate9pt4FoldSym(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int dw, const int dh, const int sw, const int sh,
	const int index)
{
	// both h and w positive
	*(dp + dh * dpitch + kb * dw) = bestOfNine(sp, spitch, kb, sw, sh, index);
	// both h & w negative
	*(dp - dh * dpitch - kb * dw) = bestOfNine(sp, -spitch, -kb, sw, sh, index);
	//  h negative w positive
	*(dp - dh * dpitch + kb * dw) = bestOfNine(sp, -spitch, kb,	sw, sh, index);
	// h  positive w negative
	*(dp + dh * dpitch - kb * dw) = bestOfNine(sp, spitch, -kb,	sw, sh, index);

}

template <typename finc>
void paint4FoldSymIfChecks(finc* dp, const int dpitch, const int kb, 
	const int width, const int height, const int centerx, const int centery,
	const int dw, const int dh, const int sw, const int sh, const finc col)
{
	// dp points to center of circle / square
	if (sh + centery < height)
	{
		if (sw + centerx < width)
			*(dp + dh * dpitch + kb * dw) = col;

		if (centerx - sw >= 0)
			*(dp + dh * dpitch - kb * dw) = col;
	}
	if (centery - dh >= 0)
	{
		if (sw + centerx < width)
			*(dp - dh * dpitch + kb * dw) = col;

		if (centerx - sw >= 0)
			*(dp - dh * dpitch - kb * dw) = col;
	}
}

template <typename finc>
void copy4FoldSymIfChecks(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int width, const int height, const int centerx, const int centery,
	const int dw, const int dh, const int sw, const int sh)
{
	// dp and sp points to centers of output and input circles
	
	if (sh + centery < height)
	{
		if (sw + centerx < width)
			*(dp + dh * dpitch + kb * dw) = *(sp + sh * spitch + kb * sw);

		if (centerx - sw >= 0)
			*(dp + dh * dpitch - kb * dw) = *(sp + sh * spitch - kb * sw);
	}
	if (centery - sh >= 0)
	{
		if (sw + centerx < width)
			*(dp - dh * dpitch + kb * dw) = *(sp - sh * spitch + kb * sw);

		if (centerx - sw >= 0)
			*(dp - dh * dpitch - kb * dw) = *(sp - sh * spitch - kb * sw);
	}

}
template <typename finc>
void interpolate4FoldSymIfChecks(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int width, const int height, const int centerx, const int centery,
	const int dw, const int dh, const int sw, const int sh,
	const int qx, const int qy, const int span,
	const float* iCoeff, finc min, finc max)
{
	// dp and sp at centers of output and input circles;
	if (sh + centery < height)
	{
		if (sw + centerx < width)
			*(dp + dh * dpitch + kb * dw) =
			clamp(LaQuantile(sp + sh * spitch + kb * sw, spitch, kb,
				span, qx, qy, iCoeff), min, max);

		if (centerx - sw >= 0)
			*(dp + dh * dpitch - kb * dw) =
			clamp(LaQuantile(sp + sh * spitch - kb * sw, spitch, -kb,
				span, qx, qy, iCoeff), min, max);
	}
	if (centery - sh >= 0)
	{
		if (sw + centerx < width)
			*(dp - dh * dpitch + kb * dw) =
			clamp(LaQuantile(sp - sh * spitch + kb * sw, -spitch, kb,
				span, qx, qy, iCoeff), min, max);

		if (centerx - sw >= 0)
			*(dp - dh * dpitch - kb * dw) =
			clamp(LaQuantile(sp - sh * spitch - kb * sw, -spitch, -kb,
				span, qx, qy, iCoeff), min, max);
	}
	
}

template <typename finc>
void interpolate9pt4FoldSymIfChecks(finc* dp, const int dpitch, const finc* sp, const int spitch, const int kb,
	const int width, const int height, const int centerx, const int centery,
	const int dw, const int dh, const int sw, const int sh,
	const int index)
{
	if (sh + centery < height)
	{
		if (sw + centerx < width)
			// both h and w positive
			*(dp + dh * dpitch + kb * dw) = bestOfNine(sp, spitch, kb,
				sw, sh, index);
		if (centerx - sw >= 0)
			// h  positive w negative
			*(dp + dh * dpitch - kb * dw) = bestOfNine(sp, spitch, -kb,
				sw, sh, index);
	}

	if (centery - sh >= 0)
	{
		if (sw + centerx < width)
			//  h negative w positive
			*(dp - dh * dpitch + kb * dw) = bestOfNine(sp, -spitch, kb,
				sw, sh, index);
		if (centerx - sw >= 0)
			// both h & w negative
			*(dp - dh * dpitch - kb * dw) = bestOfNine(sp, -spitch, -kb,
				sw, sh, index);
	}
}
#endif

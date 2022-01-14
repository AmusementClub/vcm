#pragma once
#ifndef SPOTLIGHT_DIM_A_COLOR_PLANE_V_C_MOHAN
#define SPOTLIGHT_DIM_A_COLOR_PLANE_V_C_MOHAN

// dims a color plane 

template <typename finc>
void dimplaneRGB(finc* dp, const finc* sp, int pitch,
						int wd, int ht, float dim);

template <typename finc>
void dimplaneYUV(finc* dp, const finc* sp, int pitch,
	int wd, int ht, float dim, finc limit);

template <typename finc>
void YUVspotLight(finc* dp, const finc* sp, int pitch,
	int x, int y, int r, int wd, int ht, int subW, int subH,
	finc gray, finc color);

template <typename finc>
void RGBspotLight(finc* dp, const finc* sp, int pitch,
	int x, int y, int r, int wd, int ht, finc color);


//---------------------------------------------------------------
template <typename finc>
void dimplaneRGB(finc* dp, const finc* sp, int pitch,
	int wd, int ht, float dim)
{
	for (int h = 0; h < ht; h++)
	{
		for (int w = 0; w < wd; w++)
		{
			dp[w] = (finc)(sp[w] * dim);
		}
		dp += pitch;
		sp += pitch;
	}
}

template <typename finc>
void dimplaneYUV(finc* dp, const finc* sp, int pitch,
	int wd, int ht, float dim, finc limit)
{
	for (int h = 0; h < ht; h++)
	{
		for (int w = 0; w < wd; w++)
		{
			dp[w] = (finc)((sp[w] - limit) * dim) + limit;
		}
		dp += pitch;
		sp += pitch;
	}
}
//------------------------------------------------------------
template <typename finc>
void RGBspotLight(finc* dp, const finc* sp, int pitch,
	int x, int y, int r, int wd, int ht, finc color)
{
	int sx = VSMIN(VSMAX(x - r, 0), wd - 1);
	int ex = VSMIN(VSMAX(x + r, 0), wd - 1);

	int sy = VSMIN(VSMAX(y - r, 0), ht - 1);
	int ey = VSMIN(VSMAX(y + r, 0), ht - 1);
	int rsq = r * r;

	for (int h = sy; h < ey; h++)
	{
		int hsq = (h - y) * (h - y);

		for (int w = sx; w < ex; w++)
		{
			if (hsq + (w - x) * (w - x) <= rsq)
			{				
				finc col = *(sp + (h * pitch + w));
				if ( color < col)
				{

					*(dp + h * pitch + w) = color;
				}
				else
				{
					*(dp + h * pitch + w) = col;
				}
			}
		}
	}

}
//--------------------------------------------------------------------------------------------------
template <typename finc>
void YUVspotLight(finc* dp, const finc* sp, int pitch,
	int x, int y, int r, int wd, int ht, int subW, int subH,
	finc gray, finc color)
{	
	// x, y, r , wd, ht are values for y plane. 
	// gray will be 0 for y and some value for u, v
	int sx = (VSMIN(VSMAX(x - r, 0), wd - 1)) >> subW;
	int ex = (VSMIN(VSMAX(x + r, 0), wd - 1)) >> subW;

	int sy = (VSMIN(VSMAX(y - r, 0), ht - 1)) >> subH;
	int ey = (VSMIN(VSMAX(y + r, 0), ht - 1)) >> subH;
	// initialize pointers

	int rsq = r * r;
	int rsubW = r >> subW, rsubH = r >> subH;

	//finc yuv[] = { yuvCol[3 * colFlag], yuvCol[3 * colFlag + 1], yuvCol[3 * colFlag + 2] };

	for (int h = sy; h < ey; h++)
	{
		int hsq = ( ((h << subH) - y) * ((h << subH) - y) );

		for (int w = sx; w < ex; w++)
		{
			int wsq = ( ((w << subW) - x) * ((w << subW) - x) );

			if (hsq + wsq <= rsq)
			{				
				finc col = *(sp + h * pitch + w);
				*(dp + h * pitch + w)
					= col >= gray && color > gray ||  
					 col <= gray && color < gray ? (col + color)/2
					:(col);				
			}
		}
	}

}

#endif


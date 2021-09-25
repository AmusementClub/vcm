// convert rgb integer to components red, green, blue and to Y, U, V
#ifndef COLOR_COMPONENTS_CONVERSION_H_V_C_MOHAN
#define COLOR_COMPONENTS_CONVERSION_H_V_C_MOHAN
//-------------------------------------------------------
void colorToBGRtoYUV(int color, unsigned char * bgr, unsigned char * yuv);
void colorToBGR(int color, unsigned char * bgr);
void BGRtoYUV(const unsigned char * bgr, unsigned char * yuv);
void YUVtoBGR(const unsigned char * yuv, unsigned char * bgr);
int clampcolor(int i);
//-------------------------------------------------------

void colorToBGRtoYUV(int color, unsigned char * bgr, unsigned char * yuv)
{
	// separate components
	colorToBGR(color, bgr);

	BGRtoYUV(bgr, yuv);
}
//------------------------------------------------------------------------
void colorToBGR(int color, unsigned char * bgr)
{
	// separate components
	//blue 
	bgr[0] = color & 255;
	//green
	bgr[1] = (color & 0xff00) >> 8;
	//red
	bgr[2] = (color & 0xff0000) >> 16;

}
void BGRtoYUV(const unsigned char * bgr, unsigned char * yuv)
{
	// convert bgr values  to Y u v
	yuv[0] = (unsigned char)(0.299 * bgr[2] + 0.587 * bgr[1] + 0.114 * bgr[0]);
	yuv[1] = (unsigned char)(- 0.169 * bgr[2] - 0.332 * bgr[1] + 0.5 * bgr[0] + 128);
	yuv[2] = (unsigned char)(0.5 * bgr[2] - bgr[1] * 0.419 - bgr[0] * 0.0813 + 128);
}
//------------------------------------------------------------------------------

void YUVtoBGR(const unsigned char * yuv, unsigned char * bgr)
{
	//bgr[2] = yuv[0] + 1.4075*(yuv[2] - 128);
	//bgr[1] = yuv[0] - 0.3455*(yuv[1] - 128) - 0.7169*(yuv[2] - 128);
	//bgr[0] = yuv[0] + 1.779*(yuv[1] - 128);

	//bgr[2] = yuv[0] + 1.14075*(yuv[2] - 128);
	//bgr[1] = yuv[0] - 0.3955*(yuv[1] - 128) - 0.581*(yuv[2] - 128);
	//bgr[0] = yuv[0] + 2.033*(yuv[1] - 128);

	int c = yuv[0] - 16, d = yuv[1] - 128, e = yuv[2] - 128;;

	bgr[2] = clampcolor((298 * c + 409 * e + 128) >> 8);
	bgr[1] = clampcolor((298 * c - 100 * d - 208 * e + 128) >> 8);
	bgr[0] = clampcolor((298 * c + 516 * d + 128) >> 8);

}

int clampcolor(int i)
{
	return (i < 0 ? 0 : i > 255 ? 255 : i);
}
//------------------------------------------------------------
#endif
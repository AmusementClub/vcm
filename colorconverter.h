// convert rgb integer to components red, green, blue and to Y, U, V
#ifndef COLOR_COMPONENTS_CONVERSION_H_V_C_MOHAN
#define COLOR_COMPONENTS_CONVERSION_H_V_C_MOHAN
//-------------------------------------------------------
void colorToBGRtoYUV(int color, unsigned char * bgr, unsigned char * yuv);
void colorToBGR(int color, unsigned char * bgr);
void BGRtoYUV(const unsigned char * bgr, unsigned char * yuv);
void YUVtoBGR(const unsigned char * yuv, unsigned char * bgr);



int clampRGBcolor(int i);
void YUVfromBGR(unsigned char *YUV, const unsigned char* BGR);
void BGRfromYUV(unsigned char* BGR, unsigned char* YUV);
void BGR8YUV(uint8_t* yuv, const uint8_t* bgr);
void BGR16YUV(uint16_t* yuv, const uint16_t* bgr, int nbits);
void BGR32YUV(float* yuv, const float* bgr);
void YUV8BGR(uint8_t* bgr, const uint8_t* yuv);
void YUV16BGR(uint16_t* bgr, const uint16_t* yuv, int nbits);
void YUV32BGR(float* bgr, const float* yuv);
int clampYUVcolor(int i);
// Always conversion to or from 8 bit bgr color 
void YUV_BGR8(uint8_t* bgr, uint8_t* yuv, int nbits);
void YUV_BGR8(uint8_t* bgr, uint16_t* yuv, int nbits);
void YUV_BGR8(uint8_t* bgr, float* yuv, int nbits);

void BGR8_YUV(uint8_t* yuv, uint8_t* bgr, int nbits);
void BGR8_YUV(uint16_t* yuv, uint8_t* bgr, int nbits);
void BGR8_YUV(float* yuv, uint8_t* bgr, int nbits);
template <typename finc>
finc Clamp(float val, finc min, finc max);

float clamp_color(float val, float min, float max)
{
	return val < min ? min : val > max ? max : val;
}

template <typename finc>
finc Clamp(float val, finc min, finc max)
{
	return (finc)(val < min ? min : val > max ? max : val);
}
//-------------------------------------------------------
void RGBandYUVfromColor(int col, unsigned char* BGR, unsigned char* YUV)
{
	colorToBGR(col, BGR);
	YUVfromBGR(YUV, BGR);
}

//-------------------------------------------------------------------------
void colorToBGRtoYUV(int color, unsigned char * bgr, unsigned char * yuv)
{
	// separate components
	colorToBGR(color, bgr);

	BGRtoYUV(bgr, yuv);
}
//------------------------------------------------------------------------
void colorToBGR(int color, unsigned char * bgr)
{
	// color specified as RGB
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
	yuv[0] = (unsigned char)clampYUVcolor((int)(0.299 * bgr[2] + 0.587 * bgr[1] + 0.114 * bgr[0]));
	yuv[1] = (unsigned char )clampYUVcolor((int)(- 0.169 * bgr[2] - 0.332 * bgr[1] + 0.5 * bgr[0] + 128));
	yuv[2] = (unsigned char)clampYUVcolor((int)(0.5 * bgr[2] - bgr[1] * 0.419 - bgr[0] * 0.0813 + 128));
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

bgr[2] = (unsigned char)clampYUVcolor((298 * c + 409 * e + 128) >> 8);
bgr[1] = (unsigned char)clampYUVcolor((298 * c - 100 * d - 208 * e + 128) >> 8);
bgr[0] = (unsigned char)clampYUVcolor((298 * c + 516 * d + 128) >> 8);

}

int clampYUVcolor(int i)
{
	return (i < 16 ? 16 : i > 235 ? 235 : i);
}

int clampRGBcolor(int i)
{
	return (i < 0 ? 0 : i > 255 ? 255 : i);
}
int clamp(int i)
{
	return (i < 0 ? 0 : i > 255 ? 255 : i);
}
//-----------------------------------------------------------------------------------------------
// changed from GITHUB
void YUVfromBGR(unsigned char* YUV, const unsigned char* BGR)
{
	YUV[0] = (unsigned char)clampYUVcolor((int)(0.257 * BGR[2] + 0.504 * BGR[1] + 0.098 * BGR[0] + 16));
	YUV[1] = (unsigned char)clampYUVcolor((int)(-0.148 * BGR[2] - 0.291 * BGR[1] + 0.439 * BGR[0] + 128));
	YUV[2] = (unsigned char)clampYUVcolor((int)(0.439 * BGR[2] - 0.368 * BGR[1] - 0.071 * BGR[0] + 128));
}
void BGRfromYUV(unsigned char* BGR, unsigned char* YUV)
{
	double Y = YUV[0] - 16;
	double U = YUV[1] - 128;
	double V = YUV[2] - 128;
	BGR[2] = (unsigned char)(1.164 * Y + 1.596 * V);
	BGR[1] = (unsigned char)(1.164 * Y - 0.392 * U - 0.813 * V);
	BGR[0] = (unsigned char)(1.164 * Y + 2.017 * U);
}
//------------------------------------------------------------
/*
 * RGB -> YUV conversion macros
 */
#define RGB2Y(r, g, b) (uint8_t)(((66 * (r) + 129 * (g) +  25 * (b) + 128) / 256) +  16)
#define RGB2U(r, g, b) (uint8_t)(((-38 * (r) - 74 * (g) + 112 * (b) + 128) / 256) + 128)
#define RGB2V(r, g, b) (uint8_t)(((112 * (r) - 94 * (g) -  18 * (b) + 128) / 256) + 128)
/*
YUV[0] = (unsigned char)clampYUVcolor((0.257 * BGR[2] + 0.504 * BGR[1] + 0.098 * BGR[0] + 16));
YUV[1] = (unsigned char)clampYUVcolor((-0.148 * BGR[2] - 0.291 * BGR[1] + 0.439 * BGR[0] + 128));
YUV[2] = (unsigned char)clampYUVcolor((0.439 * BGR[2] - 0.368 * BGR[1] - 0.071 * BGR[0] + 128));
*/
#define ADD16Y(n) ( 1 << (n - 4))
#define ADD16UV(n) (1 << (n  - 1))

#define BGR16Y(yuv,bgr, n) yuv[0] = (uint16_t)(0.257 * bgr[2] + 0.504 * bgr[1] + 0.098 * bgr[0] + ADD16Y(n))
#define BGR16U(yuv,bgr, n) yuv[1] = (uint16_t)(-0.148 * bgr[2] - 0.291 * bgr[1] + 0.439 * bgr[0] + ADD16UV(n))
#define BGR16V(yuv,bgr, n) yuv[2] = (uint16_t)(0.439 * bgr[2] - 0.368 * bgr[1] - 0.071 * bgr[0] + ADD16UV(n))

#define BGR32Y(yuv, bgr) yuv[0] = (float)(0.257 * bgr[2] + 0.504 * bgr[1] + 0.098 * bgr[0] )
#define BGR32U(yuv, bgr) yuv[1] = (float)(-0.148 * bgr[2] - 0.291 * bgr[1] + 0.439 * bgr[0] )
#define BGR32V(yuv, bgr) yuv[2] = (float)(0.439 * bgr[2] - 0.368 * bgr[1] - 0.071 * bgr[0] )
 /*
  *   macros that take the original Y, U, and V values
  */
#define YUV2R(y, u, v) clamp((298 * ((y)-16) + 409 * ((v)-128) + 128) / 256)
#define YUV2G(y, u, v) clamp((298 * ((y)-16) - 100 * ((u)-128) - 208 * ((v)-128) + 128) / 256)
#define YUV2B(y, u, v) clamp((298 * ((y)-16) + 516 * ((u)-128) + 128) / 256)
  //----------------------------------------------------------------------------------------
  /*
   * YUV -> RGB conversion macros
   */
   /* "Optimized" macros that take specialy prepared Y, U, and V values:
	*  C = Y - 16
	*  D = U - 128
	*  E = V - 128
	*/
#define YUV2RO(C, D, E) clamp((298 * (C) + 409 * (E) + 128) >> 8)
#define YUV2GO(C, D, E) clamp((298 * (C) - 100 * (D) - 208 * (E) + 128) >> 8)
#define YUV2BO(C, D, E) clamp((298 * (C) + 516 * (D) + 128) >> 8)
	/*
	 *   macros that take RGB and give individual components
	 */
#define RGB2R(RGB) RGB & 0xFF
#define RGB2G(RGB) (RGB >> 8) & 0xFF
#define RGB2B(RGB) (RGB >> 16) & 0xFF
#define RGB2BGR(RGB) (RGB2B << 16) | (RGB2G << 8) | (RGB2R)

	 //-------------------------------------------------------------
void RGB8ToYUV(uint8_t r, uint8_t g, uint8_t b, uint8_t* y, uint8_t* u, uint8_t* v)
{
	*y = RGB2Y((int)r, (int)g, (int)b);
	*u = RGB2U((int)r, (int)g, (int)b);
	*v = RGB2V((int)r, (int)g, (int)b);
}

void BGR2YUV(uint8_t* bgr, uint8_t* yuv)
{
	RGB8ToYUV(bgr[2], bgr[1], bgr[0], yuv, yuv + 1, yuv + 2);
}
void YUV2BGR(uint8_t* yuv, uint8_t* bgr)
{
	bgr[0] = YUV2B(yuv[0], yuv[1], yuv[2]);
	bgr[1] = YUV2G(yuv[0], yuv[1], yuv[2]);
	bgr[2] = YUV2R(yuv[0], yuv[1], yuv[2]);
}

void BGR8YUV(uint8_t* yuv, uint8_t* bgr)
{
	*yuv	=	 RGB2Y((int)(bgr[2]), (int)(bgr[1]), (int)(bgr[0]) );
	*(yuv + 1) = RGB2U((int)(bgr[2]), (int)(bgr[1]), (int)(bgr[0]) );
	*(yuv + 2) = RGB2V((int)(bgr[2]), (int)(bgr[1]), (int)(bgr[0]) );

}

void BGR16YUV(uint16_t* yuv, const uint16_t* bgr, int nbits)
{
	BGR16Y(yuv,bgr, nbits);
	BGR16U(yuv,bgr, nbits);
	BGR16V(yuv,bgr, nbits);

}

void BGR32YUV(float* yuv, const float* bgr)
{
	BGR32Y(yuv, bgr);
	BGR32U(yuv, bgr);
	BGR32V(yuv, bgr);

}
void YUV8BGR(uint8_t* bgr, const uint8_t* yuv)
{
	/*float y = yuv[0] - 16;
	//float u = yuv[1] - 128;
	//float v = yuv[2] - 128;
	//bgr[2] = (unsigned char)(1.164 * y + 1.596 * v);
	//bgr[1] = (unsigned char)(1.164 * y - 0.392 * u - 0.813 * v);
	//bgr[0] = (unsigned char)(1.164 * y + 2.017 * u);*/
	uint8_t y = yuv[0] - 16;
	uint8_t u = yuv[1] - 128;
	uint8_t v = yuv[2] - 128;
	bgr[0] = YUV2BO(y, u, v);
	bgr[1] = YUV2GO(y, u, v);
	bgr[2] = YUV2RO(y, u, v);
}

void YUV16BGR(uint16_t* bgr, const uint16_t* yuv, int nbits)
{
	float y = (float)(yuv[0] - (16 << (nbits - 8)) );
	float u = (float)(yuv[1] - (128 << (nbits - 8)));
	float v = (float)(yuv[2] - (128 << (nbits - 8)));
	bgr[2] = (uint16_t)(1.164 * y + 1.596 * v);
	bgr[1] = (uint16_t)(1.164 * y - 0.392 * u - 0.813 * v);
	bgr[0] = (uint16_t)(1.164 * y + 2.017 * u);
}


void YUV32BGR (float * bgr, const float* yuv)
{
	float y = yuv[0];
	float u = yuv[1];
	float v = yuv[2];
	bgr[2] = (float)(1.164 * y + 1.596 * v);
	bgr[1] = (float)(1.164 * y - 0.392 * u - 0.813 * v);
	bgr[0] = (float)(1.164 * y + 2.017 * u);
}


void YUV_BGR8(uint8_t* bgr, uint8_t* yuv, int nbits)
{
	bgr[2] = YUV2R(yuv[0], yuv[1], yuv[2]);
	bgr[1] = YUV2G(yuv[0], yuv[1], yuv[2]);
	bgr[0] = YUV2B(yuv[0], yuv[1], yuv[2]);
}

void YUV_BGR8(uint8_t* bgr, uint16_t* yuv, int nbits)
{
	bgr[2] = YUV2R(yuv[0] >> (nbits - 8), yuv[1] >> (nbits - 8), yuv[2] >> (nbits - 8));
	bgr[1] = YUV2G(yuv[0] >> (nbits - 8), yuv[1] >> (nbits - 8), yuv[2] >> (nbits - 8));
	bgr[0] = YUV2B(yuv[0] >> (nbits - 8), yuv[1] >> (nbits - 8), yuv[2] >> (nbits - 8));
}

void YUV_BGR8(uint8_t* bgr, float* yuv, int nbits)
{
	uint8_t y = (uint8_t)(yuv[0] * 255);
	uint8_t u = (uint8_t)(yuv[1] * 255 + 127);
	uint8_t v = (uint8_t)(yuv[2] * 255 + 127);

	bgr[2] = (uint8_t)YUV2R(y, u, v);
	bgr[1] = (uint8_t)YUV2G(y, u, v);
	bgr[0] = (uint8_t)YUV2B(y, u, v);
}
void BGR8_YUV(uint8_t* yuv, uint8_t* bgr, int nbits)
{
	yuv[0] = RGB2Y(bgr[2], bgr[1], bgr[0]);
	yuv[1] = RGB2U(bgr[2], bgr[1], bgr[0]);
	yuv[2] = RGB2V(bgr[2], bgr[1], bgr[0]);
}

void BGR8_YUV(uint16_t* yuv, uint8_t* bgr, int nbits)
{
	yuv[0] = (RGB2Y(bgr[2], bgr[1], bgr[0])) << (nbits - 8);
	yuv[1] = (RGB2U(bgr[2], bgr[1], bgr[0])) << (nbits - 8);
	yuv[2] = (RGB2V(bgr[2], bgr[1], bgr[0])) << (nbits - 8);
}

void BGR8_YUV(float* yuv, uint8_t* bgr, int nbits)
{
	yuv[0] = (RGB2Y(bgr[2], bgr[1], bgr[0])) / 255.0f;
	yuv[1] = (RGB2U(bgr[2], bgr[1], bgr[0])) / 255.0f - 0.5f;
	yuv[2] = (RGB2V(bgr[2], bgr[1], bgr[0])) / 255.0f - 0.5f;
}
#endif

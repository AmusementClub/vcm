#ifndef CONVERT_BGR_FOR_INPUT_FORMAT_V_C_MOHAN
#define CONVERT_BGR_FOR_INPUT_FORMAT_V_C_MOHAN
// input 8 bit bgr array. color array must be 12 bytes
void convertBGRforInputFormat(uint8_t* color, const uint8_t* bgr, const VSFormat* fi);

void convertBGRforInputFormat(uint8_t* color, const uint8_t* bgr, const VSFormat* fi)
{

	uint8_t yuv[3];
	
	BGRtoYUV( bgr, yuv);

	int nbytes = fi->bytesPerSample;
	int nbits = fi->bitsPerSample;

	for (int k = 0; k < 3; k++)
	{
		if (nbytes == 1)
		{
			if (fi->colorFamily == cmRGB)
				color[k] = bgr[k];
			else
				color[k] = yuv[k];
		}
		else if (nbytes == 2)
		{
			if (fi->colorFamily == cmRGB)
				*((uint16_t*)color + k) = (uint16_t)(bgr[k] << (nbits - 8));
			else
				*((uint16_t*)color + k) = (uint16_t)(yuv[k] << (nbits - 8));
		}
		else // float
		{
			if (fi->colorFamily == cmRGB)
				*((float*)color + k) = (float)(bgr[k] / 255.0f);
			else
			{
				if (k == 0)
					*((float*)color + k) = (float)((yuv[k] - 16) / 220.0f);
				else
					*((float*)color + k) = (float)((yuv[k] - 128) / 220.0f);
			}
		}

	}
}


#endif
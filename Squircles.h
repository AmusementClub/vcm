#ifndef SQUIRCLES_H_V_C_MOHAN
#define SQUIRCLES_H_V_C_MOHAN
// Assume #define _USE_MATH_DEFINES
// All coordinates are for unit circle ( rad = 1) and square (-1,1, -1, 1)
void SquareToCircularDisc(float* uv, const float x, const float y);
void CircularDiscToSquare(float* xy, const float cU, const float cV);
// normalize,convert and get back in input units
void getSquircleUV(float* cUV, const int w, const int h, const int rad);
void getxySquircle(float* xy, const float cU, const float cV, const int rad);
// 
// mapping a circular disc to a square region
// input: (u,v) coordinates in the circle
// output: xy[2] coordinates in the square
void CircularDiscToSquare( float * xy, const float u, const float v)
{
    double u2 = (double)(u * u);
    double v2 = (double)(v * v);
    double twosqrt2 = 2.0 * M_SQRT2;
    double subtermx = 2.0 + u2 - v2;
    double subtermy = 2.0 - u2 + v2;
    double termx1 = subtermx + u * twosqrt2;
    double termx2 = subtermx - u * twosqrt2;
    double termy1 = subtermy + v * twosqrt2;
    double termy2 = subtermy - v * twosqrt2;
    xy[0] = (float)(0.5 * sqrt(termx1) - 0.5 * sqrt(termx2));
    xy[1] = (float)(0.5 * sqrt(termy1) - 0.5 * sqrt(termy2));

}


// 
// mapping a square region to a circular disc
// input: (x,y) coordinates in the square
// output: uv[2] coordinates in the circle
void SquareToCircularDisc(float * uv, const float x, const float y)
{
    uv[0] = (float)(x * sqrt(1.0 - (y * y) / 2.0));
    uv[1] = (float)(y * sqrt(1.0 - (x * x) / 2.0));
}

void getSquircleUV(float* cUV, const int w, const int h, const int rad)
{
    // normalize
    float ww = (float)w / (float)rad;
    float hh = (float)h / (float)rad;
    SquareToCircularDisc(cUV, ww, hh);
    // get back in original units
    cUV[0] *= rad;
    cUV[1] *= rad;
}

void getxySquircle(float* xy, const float cU, const float cV, const int rad)
{
    // normalize
    float cUU = cU / rad;
    float cVV = cV / rad;
    CircularDiscToSquare(xy, cUU, cVV);
    
    // get back in original units
    xy[0] *= rad;
    xy[1] *= rad;
}
#endif


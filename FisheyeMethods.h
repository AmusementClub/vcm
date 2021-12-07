#pragma once
#ifndef FISHEYE_METHODS_H_V_C_MOHAN
#define FISHEYE_METHODS_H_V_C_MOHAN
/* this file has versions for calculating sourc coordinates
for corresponding output coordinates for correcting or
producing fisheye images
methods are :-

1 Fish eye ortho correction
2 Fisheye linear correction
3.Fisheye equisolid correction
4.Fisheye panoramic correction
5 Fisheye radial correction
6 Produce Fisheye ortho
7 Produce Fisheye linear
8 Produce Fisheye equisolid
9 Produce Fisheye panoramic
10 Produce Fisheye radial
 barrel and Pincushion 2 methods:-
 Barrel and pin cushion  type distortion correction with a,b,c,d old method
Barrel and PinCushion type distortion correction with new c */

//  rNorm radius for normalization used by methods 
void getSourceCoord(float* xy, const float ww, const float hh, const int method,
	const float rNormal, const float* abc);
// Results in xy array of 2. ww and xx are coordinates
// focal is focal length and rix is akin to refractive index
// used by fisheye methods
void getSourceXY(float* xy, const float ww, const float hh, const int method,
	const double focal, const double rNorm, const double rix);
// method 1,6 ortho, 2,7 linear, 3,8 equisolid,4,9 Panoramic, 5,10 radial
double getFocalLength(const int frad, const int method, const double fov);
int getOutputRadius(const int frad, const double focal,  const double rix);

int getOutputRadius(const int frad, const double focal,  const double rix)
{
	
	double alfa = atan(frad / focal);
	double theta = asin(sin(alfa) * rix);
	return  (int)(focal * tan(theta));

}
double getFocalLength(const int frad, const int method, const double fov)
{
	// type 1,6 ortho, 2,7 linear, 3,8 equisolid,4,9 Panoramic, 5,10 radial
	double focal;
	double fovrad = fov * M_PI / 360.0;
	switch (method)
	{
	case 1:
	{ }
	case 6 :
	{
		focal = (frad / sin(fovrad));
		break;
	}
	case 2 :
	{}
	case 7:
	{
		focal = (frad / (fovrad)  );
		break;
	}
	case 3:
	{}
	case 8 :
	{
		focal = (frad / (2 * sin(fovrad / 2) ) );
		break;
	}
	case 4:
	{}
	case 9: 
	{
		focal = (frad / (2 * tan(fovrad / 2)));
		break;
	}
	case 5 :
	{}
	case 10:
	{
		focal = (frad / tan(fovrad));
		break;
	}
	}
	return focal;
}

void getSourceXY(float* xy, const float ww, const float hh, const int method,
	const double focal, const double rNorm, const double rix)
{
	
	double h = (double)hh;
	double w = (double)ww + 0.5;// to ensure no 0 / 0  tan does not become infinity

	switch (method)
	{
	case  11:
	{

		// barrel
		double rsq = hh * hh + ww * ww;
		double denom = 1.0 - focal * rsq / rNorm;
		xy[0] = (float)(ww / denom);
		xy[1] = (float)(hh / denom);
		break;
	}

	case 12:
	{

		// pin cushion
		double rsq = hh * hh + ww * ww;
		double denom = 1.0 + focal * rsq / rNorm;
		xy[0] = (float)(ww / denom);
		xy[1] = (float)(hh / denom);
		break;
	}
	case 1:
	{
		//ortho f sin@
		double rpersp = sqrt(1.0 * (h * h + w * w));
		double theta = atan(rpersp / focal);
		double alfa = asin(sin(theta) / rix);

		double rfish = focal * sin(alfa);

		xy[0] = (float)(rfish * w / rpersp);
		xy[1] = (float)(rfish * h / rpersp);
		break;
	}
	case 2:
	{
		//linear f@
		double rpersp = sqrt(1.0 * (h * h + w * w));
		double theta = atan(rpersp / focal);
		double alfa = asin(sin(theta) / rix);

		double rfish = focal * alfa;
		xy[0] = (float)(rfish * w / rpersp);
		xy[1] = (float)(rfish * h / rpersp);
		break;
	}
	case 3:
	{
		//equisolid 2 f sin (@ / 2)
		double rpersp = sqrt(1.0 * (h * h + w * w));
		double theta = atan(rpersp / focal);
		double alfa = asin(sin(theta) / rix);

		double rfish = 2 * focal * sin(alfa / 2);
		xy[0] = (float)(rfish * w / rpersp);
		xy[1] = (float)(rfish * h / rpersp);
		break;
	}
	case 4:
	{
		//panoramic 2 f tan(@ / 2)
		double rpersp = sqrt(1.0 * (h * h + w * w));
		double theta = atan(rpersp / focal);
		double alfa = asin(sin(theta) / rix);

		double rfish = 2 * focal * tan(alfa / 2);
		xy[0] = (float)(rfish * w / rpersp);
		xy[1] = (float)(rfish * h / rpersp);
		break;
	}
	case 5:
	{
		//from  fish to radial
		//tan theta get source
		double rpersp = sqrt(1.0 * (h * h + w * w));
		double theta = atan(rpersp / focal);
		//double alfa = atan(tanalfa);
		double sinalfa = sin(theta) / rix;
		double alfa = asin(sinalfa);
		double rfish = focal * tan(alfa);
		xy[0] = (float)(rfish * w / rpersp);
		xy[1] = (float)(rfish * h / rpersp);
		break;
	}
	case 6:
	{
		//for  ortho fish f sin@ get source ortho coord
		double rfish = sqrt(1.0 * (h * h + w * w));
		double alfa = atan(rfish / focal);
		double sintheta = sin(alfa) * rix;	// same as rfish * rix / focal		
		double theta = asin(sintheta);
		double rortho = focal * tan(theta);

		xy[0] = (float)(rortho * w / rfish);
		xy[1] = (float)(rortho * h / rfish);
		break;
	}
	case 7:
	{
		//for linear fish f @ get source 
		double rfish = sqrt(1.0 * (h * h + w * w));
		double alfa = (rfish / focal);
		double sintheta = sin(alfa) * rix;
		double theta = asin(sintheta);
		double rlinear = focal * tan(theta);

		xy[0] = (float)(rlinear * w / rfish);
		xy[1] = (float)(rlinear * h / rfish);
		break;
	}
	case 8:
	{
		//for equi solid fish 2 f sin (@/2) get source
		double rfish = sqrt(1.0 * (h * h + w * w));
		double sinalfa2 = rfish / (2 * focal);
		double alfa2 = asin(sinalfa2);
		double alfa = 2 * alfa2;
		double sintheta = sin(alfa) * rix;
		double theta = asin(sintheta);
		double rsolid = focal * tan(theta);

		xy[0] = (float)(rsolid * w / rfish);
		xy[1] = (float)(rsolid * h / rfish);
		break;
	}

	case 9:
	{
		//from panoramic to fish
		//for panoramic fish 2 f tan (@/2) get source
		double rfish = sqrt(1.0 * (h * h + w * w));
		double tanalfa2 = rfish / (2 * focal);
		double alfa2 = atan(tanalfa2);
		double alfa = 2 * alfa2;
		double sintheta = sin(alfa) * rix;
		double theta = asin(sintheta);
		double rpan = focal * tan(theta);


		xy[0] = (float)(rpan * w / rfish);
		xy[1] = (float)(rpan * h / rfish);
		break;
	}
	case 10:
	{
		//from Radial to fish
		//for panoramic fish 2 f tan (@/2) get source
		double rfish = sqrt(1.0 * (h * h + w * w));
		double tanalfa = rfish / (focal);
		double alfa = atan(tanalfa);
		double sintheta = sin(alfa) * rix;
		//double sinalfa = rix * asin(sintheta);
		double theta = asin(sintheta);
		double rRadial = focal * tan(theta);

		xy[0] = (float)(rRadial * w / rfish);
		xy[1] = (float)(rRadial * h / rfish);
		break;
	}


	default:
		break;
	}
}

void getSourceCoord(float* xy, const float ww, const float hh, const int method,
	const float rNormal, const float *abc)
{
	if (method == 1)
	{	
		//Brown Conrady model
		float rNorm =  (float)sqrt((double)(hh * hh  + ww * ww)) / rNormal;
		float rnew = (((abc[0] * rNorm

			+ abc[1]) * rNorm

			+ abc[2]) * rNorm

			+ (1.0f - abc[0] - abc[1] - abc[2]) );	// d

			xy[0] = rnew * ww;
			xy[1] = rnew * hh;
	}
	else if (method == 2)
	{
		// division model
		float radsq = (float)(hh * hh + ww * ww)/ (rNormal * rNormal);
		float denom = ( 1.0f 
			- (((abc[0] * radsq
			+ abc[1]) * radsq
			+ abc[2]) * radsq) );
		xy[0] = ww / denom;
		xy[1] = hh / denom;
	}
}
#endif

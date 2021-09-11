
//-------------------------------------------------------------------------------
#ifndef BUTTERWORTHFILTERSFORF2QUIVER_V_C_MOHAN
#define BUTTERWORTHFILTERSFORF2QUIVER_V_C_MOHAN
// Butterworth filter functions required by F2Quiver of FFTQuiver plugin
//--------------------------------------------------------------------------------
// declarations
	void ButterworthCircular2D(float *filt, int * fspec, int ht,int wd);
	void ButterworthHorizontal2D(float *filt, int * fspec, int ht,int wd);
	void ButterworthVertical2D(float *filt, int * fspec, int ht,int wd);
	void pointNotchFilter2D( float * freqFilter, int* spec, int ht, int wd);
	void FanFilter2D( float * freqFilter, int* spec, int ht, int wd);
	//------------------------------------------------------------------------------

	// Definitions
//-----------------------------------------------------------------------------------

//====================================================================================================================================
void FanFilter2D( float * freqFilter, int* spec, int ht, int wd)
{
	float pi = M_PI;
	int type = spec[1];

	float theta1 = spec[2] * pi / 180;	// fan angle start
	float theta2 = spec[3] * pi / 180;	// fan angle end

	int degree = 2 * spec[4];

	if(  spec[2] <= 88  )	// values range from 0 to 180
	{
		
		for(int w = 0; w < wd; w ++)
		{
			float hfreq = (w - wd / 2) * tan (theta2) ; // sin(theta2) / cos(theta2);
			float lfreq = (w - wd / 2) * tan (theta1); //sin(theta1) / cos(theta1);
			float * flt = freqFilter;
			float dsqrh = hfreq * hfreq, dsqrl = lfreq * lfreq;
			
			for( int h = 0; h < ht; h ++) 
			{
				float hsqr = (h - ht / 2) * (h - ht / 2);
				float mfactor = 1.0;
				if( (h - ht/2) * (w - wd /2) > 0) 
				{				

					if ( dsqrh > 0)
						mfactor *=   1.0 / ( 1.0 + pow ( hsqr /dsqrh, degree) );
					else  					mfactor = 0;
					if (hsqr > 0)
						mfactor *= 1.0 / ( 1.0 + pow( dsqrl / hsqr, degree) ) ;
				else 
					mfactor = 0;

				flt[w] *= 1.0 - mfactor;
				}																								

				flt += wd;
			} 
		}
	}

		else if(  spec[2] > 90  )	// values range from 0 to 180
		{
		
			for(int w = 0; w < wd; w ++)
			{
				float lfreq = (w - wd / 2) * tan (pi - theta2) ; // sin(theta2) / cos(theta2);
				float hfreq = (w - wd / 2) * tan (pi - theta1); //sin(theta1) / cos(theta1);
				//float * flt = freqFilter;
				float dsqrh = hfreq * hfreq, dsqrl = lfreq * lfreq;
				float * flt = freqFilter + (ht - 1 ) * wd;

				for( int h = ht - 1; h >= 0; h --) 
				{
					float hsqr = (h - ht / 2) * (h - ht / 2);
					float mfactor = 1.0;
					if( (h - ht/2) * (w - wd /2) < 0) 
					{				

					if ( dsqrh > 0)
						mfactor *=   1.0 / ( 1.0 + pow ( hsqr /dsqrh, degree) );
					else
						mfactor = 0;
					if (hsqr > 0)
						mfactor *= 1.0 / ( 1.0 + pow( dsqrl / hsqr, degree) ) ;
					else 
						mfactor = 0;

					flt[w] *= 1.0 - mfactor;
				}																								

				flt -= wd;
			} 
		}
	}	
}
//-------------------------------------------------------------------------------------------------------------------------------------------------
void pointNotchFilter2D( float * freqFilter, int* spec, int ht, int wd)
{
	// point centered notch filter
	float xf = spec[2] * wd / (NYQUIST);

	float yf = spec[3] * ht / (NYQUIST );

	int degree = 2 * spec[4];

	// calculate notch radial sqr
	float d0sq = sqrt (xf*xf + yf*yf);
	d0sq = (ht * ht + wd * wd)/ 4;
	for(int h = 0; h < ht; h ++)
	{
		for(int w = 0; w < wd; w ++)
		{
			float d1 = ( (h - ht / 2 - yf) * (h - ht / 2 - yf) + (w - wd / 2 - xf) * (w - wd / 2 - xf) );
			float d2 = ( (h - ht / 2 + yf) * (h - ht / 2 + yf) + (w - wd / 2 + xf) * (w - wd / 2 + xf) );
		
			if(d1 * d2 > 0.0 )
		
				freqFilter[ w] *= 1.0 / (( 1.0 + pow(1 + 1/d1, degree)) * (1.0 +  pow(1 + 1.0/d2, degree)) );
			else
				freqFilter[ w] = 0.0;
				
		}

		freqFilter += wd;
	}
}
			
//--------------------------------------------------------------------------------

void ButterworthCircular2D(float *filt, int * fspec, int ht,int wd)
{
	int type = fspec[1];	// 1 locut, 2 hicut, 3 band pass 4 band stop 5 Notch
	int vf1 = (fspec[2] * wd ) /(NYQUIST);
	int hf1 = (fspec[2] * ht ) / (NYQUIST);
	int vf2 = (fspec[3] * wd ) /( NYQUIST);
	int hf2 = (fspec[3] * ht ) /(NYQUIST);
	int degree = fspec[4] == 0 ? 1 : 2 * fspec[4];
	float dsqr = (vf1 * vf1); 
	float dsqrband =(vf2 * vf2) ; 
	float ratiow2h = (1.0 * wd * wd) /(ht * ht);
	switch (type)
	{
		case 1 :	// locut or high freq pass
		{
			

			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr = ( h  - ht / 2 ) * ( h  - ht / 2 ) * ratiow2h;				

				for ( int w = 0; w < wd ; w ++)
				{
					float d2sqr = ( w - wd / 2 ) * ( w - wd / 2 ); 

					if ( d1sqr + d2sqr > 0 )

						filt[w] *= 1.0 / (1.0 + pow( dsqr / ( d1sqr + d2sqr) , degree) );
					else
						filt[w] = 0;
				}

				filt += wd;
			}

			break;
		}

		case 2 : // highcut or low freq pass 1.0 - locut filter
		{

			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr = ( h - ht / 2 ) * ( h - ht / 2 ) * ratiow2h;				

				for ( int w = 0; w < wd ; w ++)
				{
					float d2sqr = ( w - wd / 2 ) * ( w - wd / 2 ); 

					if ( dsqr > 0 )

						filt[w] *=  1.0 / (1.0 + pow( ( d1sqr + d2sqr)/ dsqr , degree) );
					else
						filt[w] = 0;
				}

				filt += wd;
			}

			break;
		}

		case 3 :	// band pass (locut * hicut)
		{
			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr = ( h - ht / 2 ) * ( h - ht / 2 ) * ratiow2h;				
			//	float d1sqrband = ( h - ht/ 2 ) * ( h - ht / 2 );
				float mfactor = 1.0;

				for ( int w = 0; w < wd ; w ++)
				{
					float d2sqr = ( w - wd / 2 ) * ( w - wd / 2 ); 
			//		float d2sqrband = ( w - wd / 2 ) * ( w - wd / 2 );
					if ( d1sqr + d2sqr > 0 )

						mfactor = 1.0 / (1.0 + pow( dsqr / ( d1sqr + d2sqr) , degree) );
					else
						mfactor = 0;

					if ( dsqr > 0 )

						mfactor *=  1.0 / ( 1.0 + pow (  ( d1sqr + d2sqr) / dsqrband , degree) );

					filt[w] *= mfactor;
				}

				filt += wd;
			}

			break;
		}

		case 4 :	// band reject (locut * hicut)
		{
			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr = ( h - ht / 2 ) * ( h - ht / 2 ) * ratiow2h;				
			//	float d1sqrband = ( h - ht/ 2 ) * ( h - ht / 2 );
				float mfactor = 1.0;

				for ( int w = 0; w < wd ; w ++)
				{
					float d2sqr = ( w - wd / 2 ) * ( w - wd / 2 ); 
				//	float d2sqrband = ( w - wd / 2 ) * ( w - wd / 2 );
					if ( d1sqr + d2sqr > 0 )

						mfactor = 1.0 / (1.0 + pow( dsqr / ( d1sqr + d2sqr) , degree) );
					else
						mfactor = 0;

					if ( dsqrband > 0 )

						mfactor *= 1.0 / ( 1.0 + pow (  ( d1sqr + d2sqr) / dsqrband , degree) );

					filt[w] *= 1.0 - mfactor;
				}

				filt += wd;
			}

			break;
		}

		case 5 :	// notch reject (locut * hicut)
		{
			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr = ( h - ht / 2 ) * ( h - ht / 2 ) * ratiow2h;				
			//	float d1sqrband = ( h - ht/ 2 ) * ( h - ht / 2 );
				float mfactor = 1.0;

				for ( int w = 0; w < wd ; w ++)
				{
					float d2sqr = ( w - wd / 2 ) * ( w - wd / 2 ); 
				//	float d2sqrband = ( w - wd / 2 ) * ( w - wd / 2 );
					if ( d1sqr + d2sqr > 0 )
					{

						mfactor = 1.0 / (1.0 + pow( dsqr / ( d1sqr + d2sqr) , degree) );
					
						mfactor *= 1.0 - 1.0 / ( 1.0 + pow ( dsqr / ( d1sqr + d2sqr) , degree) );
					}
					else
						mfactor = 0;

					filt[w] *= 1.0 - mfactor;
				}

				filt += wd;
			}

		}
		
		default : {}


	}
}

//--------------------------------------------------------------------------------

void ButterworthHorizontal2D(float *filt, int * fspec, int ht,int wd)
{
	int type = fspec[1];	// 1 locut, 2 hicut, 3 band pass 4 band stop 5 Notch
	int vf1 = (fspec[2] * ht ) /(NYQUIST);
//	int hf1 = (fspec[2] * ht ) / NYQUIST;
	int vf2 = (fspec[3] * ht ) /( NYQUIST);
//	int hf2 = (fspec[3] * ht ) / NYQUIST;
	int degree = fspec[4] == 0 ? 1 : 2 * fspec[4];
	float dsqr = (vf1 * vf1 );
	float dsqrband =(vf2 * vf2  );
	
	switch (type)
	{
		case 1 :	// locut or high freq pass
		{
			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr = ( h - ht / 2 ) * ( h - ht / 2) ;				
				float mfactor = 1.0;
				
					if ( d1sqr > 0 )

						mfactor *= 1.0 / (1.0 + pow( dsqr / ( d1sqr ) , degree) );
					else
						mfactor = 0;

				for ( int w = 0; w < wd ; w ++)
				{
					filt[w] *= mfactor;
				}

				filt += wd;
			}

			break;
		}

		case 2 : // highcut or low freq pass 1.0 - locut filter
		{

			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr =  ( h - ht / 2 )*( h - ht / 2)  ;				
				float mfactor = 1.0;
				
				if ( d1sqr > 0 )

					mfactor =   1.0 / (1.0 + pow( d1sqr / ( dsqr) , degree) );

				else

					mfactor = 1.0;
					
				for ( int w = 0; w < wd ; w ++)
				{
					filt[w] *= mfactor;
				}

				filt += wd;
			}

			break;
		}

		case 3 :	// band pass (locut * hicut)
		{
			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr =  ( h - ht / 2 )*( h - ht / 2)  ;				
				//float d1sqrband = ( h - ht/ 2 ) * ( h - ht / 2 );
				float mfactor = 1.0;

				
					if ( d1sqr  > 0 )					

						mfactor = 1.0 / (1.0 + pow( dsqr / ( d1sqr ) , degree) );
					else

						mfactor = 0.0;
					
					mfactor *= 1.0 / ( 1.0 + pow ( d1sqr / ( dsqrband ) , degree) );
					
				for ( int w = 0; w < wd ; w ++)
				{

					filt[w] *= mfactor;
				}

				filt += wd;
			}

			break;
		}

		case 4 :	// band reject (locut * hicut)
		{
			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr =  ( h - ht / 2 )*( h - ht / 2)  ;				
			//	float d1sqrband = ( h - ht/ 2 ) * ( h - ht / 2 );
				float mfactor = 1.0;

				
					if ( d1sqr  > 0 )

						mfactor = 1.0 / (1.0 + pow( dsqr / ( d1sqr ) , degree) );
					else
						mfactor = 0.0;

					if ( dsqrband  > 0 )

						mfactor *=  1.0 / ( 1.0 + pow ( d1sqr / ( dsqrband ) , degree) );
					else
						mfactor = 0.0;

				for ( int w = 0; w < wd ; w ++)
				{

					filt[w] *= 1.0 - mfactor;
				}

				filt += wd;
			}

			break;
		}

		case 5 :	// notch reject (locut * hicut)
		{
			for (int h = 0; h < ht ; h ++)
			{
				float d1sqr =  ( h - ht / 2 )*( h - ht / 2)  ;															
			//	float d1sqrband = ( h - ht/ 2 ) * ( h - ht / 2 );
				float mfactor = 1.0;

				
					
					if ( d1sqr  > 0 )

						mfactor = 1.0 / (1.0 + pow( dsqr / ( d1sqr ) , degree) );
					else
						mfactor = 0.0;

					if ( dsqr > 0 )

						mfactor *= 1.0 / ( 1.0 + pow ( d1sqr / ( dsqr ) , degree) );

					else
						mfactor = 0.0;

				for ( int w = 0; w < wd ; w ++)
				{

					filt[w] *= 1.0 - mfactor;
				}

				filt += wd;
			}

		}
		
		default : {}


	}
}

//--------------------------------------------------------------------------------

void ButterworthVertical2D(float *filt, int * fspec, int ht,int wd)
{
	int type = fspec[1];	// 1 locut, 2 hicut, 3 band pass 4 band stop 5 Notch
	//int vf1 = (fspec[2] * ht ) /(2 * NYQUIST);
	int hf1 = (fspec[2] * wd ) / NYQUIST; 
	//int vf2 = (fspec[3] * ht ) /(2 * NYQUIST);
	int hf2 = (fspec[3] * wd ) / NYQUIST; 
	int degree = fspec[4] == 0 ? 1 : 2 * fspec[4];
	float dsqr = ( hf1 * hf1 );
	float dsqrband =( hf2 * hf2 );

	switch (type)
	{
		case 1 :	// locut or high freq pass
		{			

			for (int w = 0; w < wd ; w ++)
			{
				float * flt = filt + w;

				float d2sqr = ( w - wd / 2 ) * ( w - wd / 2 ) ;

				float mfactor;

				
					if (  d2sqr > 0 )

						mfactor =  1.0 / (1.0 + pow( dsqr / (  d2sqr) , degree) );
					else
						mfactor = 0;

				for ( int h = 0; h < ht ; h ++)
				{
					flt[0] *= mfactor;

					flt += wd;
				}
			}

			break;
		}

		case 2 : // highcut or low freq pass 1.0 - locut filter
		{

			for (int w = 0; w < wd ; w ++)
			{
				float * flt = filt + w;

				float d2sqr =  ( w - wd / 2 ) * ( w - wd / 2 ) ; 				

				float mfactor;

					if (  d2sqr > 0 ) 					

						mfactor =   1.0 / (1.0 + pow( d2sqr / ( dsqr) , degree) );
					else
						mfactor = 1.0;

				for ( int h = 0; h < ht ; h ++)
				{
					flt[0] *= mfactor;

					flt += wd;
				}
			}

			break;
		}

		case 3 :	// band pass (locut * hicut)
		{
			for (int w = 0; w < wd ; w ++)
			{
				float * flt = filt + w;

				float d2sqr = ( w - wd / 2 ) * ( w - wd / 2 ) ; 				
			//	float d2sqrband = abs( w - wd / 2 );

				float mfactor = 1.0;				
					
					if ( d2sqr > 0 )

						mfactor = 1.0 /  (1.0 + pow( dsqr / ( d2sqr) , degree) );
							   
					else
						mfactor = 0.0;

					if ( dsqrband > 0)

						mfactor *= 1.0 / ( 1.0 + pow ( d2sqr / ( dsqrband ) , degree) )  ;
				//	else
				//		mfactor = 1.00;
					
				for ( int h = 0; h < ht ; h ++)
				{

					flt[0] *= mfactor;				

					flt += wd;
				}
			}

			break;
		}

		case 4 :	// band reject (locut * hicut)
		{
			for (int w = 0; w < wd ; w ++)
			{
				float * flt = filt + w;

				float d2sqr = ( w - wd / 2 )* ( w - wd / 2 ) ; 				
			//	float d2sqrband = abs( w - wd / 2 );

				float mfactor = 1.0;										
				if (d2sqr > 0)

					mfactor =  1.0 /  (1.0 + pow( dsqr / ( d2sqr) , degree) );	
				else
					mfactor = 0;

					if(dsqrband > 0)

						mfactor *= 1.0 / ( 1.0 + pow ( d2sqr / ( dsqrband) , degree) ) ;

					else
						mfactor = 0;
					
				for ( int h = 0; h < ht ; h ++)
				{

					flt[0] *= 1 -  mfactor;				

					flt += wd;
				}
			}

			break;
		}

		case 5 :	// notch reject (locut * hicut)
		{
			
			float mfactor = 1.0;

			for (int w = 0; w < wd ; w ++)
			{
				float * flt = filt + w;

				float d2sqr = ( w - wd / 2 ) * ( w - wd / 2 ) ; 			

				float mfactor = 1.0;				
					
					if ( d2sqr > 0 )					

						mfactor = 1.0 /  (1.0 + pow( dsqr / ( d2sqr) , degree) );
							 
					else
						mfactor = 0.0;

					if ( dsqr > 0)

						mfactor *=  1.0 / ( 1.0 + pow ( d2sqr / ( dsqr ) , degree) );
					else
						mfactor = 0;
					
				for ( int h = 0; h < ht ; h ++)
				{

					flt[0] *= 1.0 - mfactor;				

					flt += wd;
				}
			}

			break;
		}
		
		default : {}


	}
}



#endif
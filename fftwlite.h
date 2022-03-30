// Lite version of fftw header on base of fftw3.h
// some needed fftwf typedefs added  for delayed loading
// (by Fizick) modified by vcmohan
// 
#ifndef __FFTWLITE_H__
#define __FFTWLITE_H__

#ifndef STATIC_FFTW3
typedef float fftwf_complex[2];
typedef struct fftwf_plan_s  *fftwf_plan;
#else
#include <fftw3.h>
#endif

typedef fftwf_complex* (*fftwf_malloc_proc)(size_t n); 

typedef void (*fftwf_free_proc) (void *ptr);
		// 1d transforms
typedef fftwf_plan (*fftwf_plan_dft_r2c_1d_proc)( int wbest, float *in, fftwf_complex *out, int flags);
typedef fftwf_plan (*fftwf_plan_dft_c2r_1d_proc)( int wbest, fftwf_complex *out, float *in, int flags);
		// 2d transforms
typedef fftwf_plan (*fftwf_plan_dft_r2c_2d_proc) (int winy, int winx, float *in, fftwf_complex *out, int flags);
typedef fftwf_plan (*fftwf_plan_dft_c2r_2d_proc) (int winy, int winx, fftwf_complex *out, float *in, int flags);

typedef fftwf_plan (*fftwf_plan_dft_2d_proc)(int winy, int winx, fftwf_complex *in, fftwf_complex *out, int inv, int flags);

//typedef fftwf_plan (*fftwf_plan_many_dft_r2c_proc) (int rank, const int *n,	int howmany,  float *in, const int *inembed, int istride, int idist, fftwf_complex *out, const int *onembed, int ostride, int odist, unsigned flags);
//typedef fftwf_plan (*fftwf_plan_many_dft_c2r_proc) (int rank, const int *n,	int howmany,  fftwf_complex *out, const int *inembed, int istride, int idist, float *in, const int *onembed, int ostride, int odist, unsigned flags);

typedef void (*fftwf_destroy_plan_proc) (fftwf_plan);

typedef void (*fftwf_execute_dft_r2c_proc) (fftwf_plan, float *realdata, fftwf_complex *fftsrc);
typedef void (*fftwf_execute_dft_c2r_proc) (fftwf_plan, fftwf_complex *fftsrc, float *realdata);
typedef void (*fftwf_execute_proc) (fftwf_plan);

typedef int  (*fftwf_init_threads_proc) (void);
typedef void (*fftwf_plan_with_nthreads_proc) (int nthreads);
typedef void (*fftwf_cleanup_threads_proc)(void);
typedef void (*fftwf_cleanup_proc)(void);

#ifndef STATIC_FFTW3
#define FFTW_MEASURE (0U)
#define FFTW_ESTIMATE (1U << 6)
#define FFTW_DESTROY_INPUT (1U << 0)
#endif
#endif
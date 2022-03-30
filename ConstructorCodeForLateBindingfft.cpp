
#ifndef STATIC_FFTW3
// late binding

bool ok = true;
HINSTANCE hinstLib = LoadLibrary(TEXT("libfftw3f-3.dll"));

if (hinstLib == NULL)
hinstLib = LoadLibrary(TEXT("fftw3.dll"));

if (hinstLib == NULL)
hinstLib = LoadLibrary(TEXT("fftw.dll"));

if (hinstLib != NULL)
{
	d->hinstLib = hinstLib;
	d->fftwf_free = (fftwf_free_proc)GetProcAddress(d->hinstLib, "fftwf_free");
	d->fftwf_malloc = (fftwf_malloc_proc)GetProcAddress(d->hinstLib, "fftwf_malloc");
	// for 1 d fft
	d->fftwf_plan_dft_r2c_1d = (fftwf_plan_dft_r2c_1d_proc)GetProcAddress(d->hinstLib, "fftwf_plan_dft_r2c_1d");
	d->fftwf_plan_dft_c2r_1d = (fftwf_plan_dft_c2r_1d_proc)GetProcAddress(d->hinstLib, "fftwf_plan_dft_c2r_1d");
	// for 2d fft
	d->fftwf_plan_dft_r2c_2d = (fftwf_plan_dft_r2c_2d_proc)GetProcAddress(d->hinstLib, "fftwf_plan_dft_r2c_2d");
	d->fftwf_plan_dft_c2r_2d = (fftwf_plan_dft_c2r_2d_proc)GetProcAddress(d->hinstLib, "fftwf_plan_dft_c2r_2d");

	

	d->fftwf_destroy_plan = (fftwf_destroy_plan_proc)GetProcAddress(d->hinstLib, "fftwf_destroy_plan");
	d->fftwf_execute = (fftwf_execute_proc)GetProcAddress(d->hinstLib, "fftwf_execute");

	d->fftwf_execute_dft_r2c = (fftwf_execute_dft_r2c_proc)GetProcAddress(d->hinstLib, "fftwf_execute_dft_r2c");
	d->fftwf_execute_dft_c2r = (fftwf_execute_dft_c2r_proc)GetProcAddress(d->hinstLib, "fftwf_execute_dft_c2r");
	// required for internal threading
//	fftwf_init_threads = (fftwf_init_threads_proc)GetProcAddress(hinstLib, "fftwf_init_threads");
//	fftwf_plan_with_nthreads = (fftwf_plan_with_nthreads_proc)GetProcAddress(hinstLib, "fftwf_plan_with_nthreads");
//	fftwf_cleanup_threads = (fftwf_cleanup_threads_proc)GetProcAddress(hinstLib, "fftwf_cleanup_threads");
//	fftwf_cleanup = (fftwf_cleanup_proc)GetProcAddress(hinstLib, "fftwf_cleanup");

}
if (d->hinstLib == NULL || d->fftwf_free == NULL || d->fftwf_malloc == NULL
	|| d->fftwf_plan_dft_r2c_1d == NULL || d->fftwf_plan_dft_c2r_1d == NULL
	|| d->fftwf_plan_dft_r2c_2d == NULL || d->fftwf_plan_dft_c2r_2d == NULL
	|| d->fftwf_destroy_plan == NULL || d->fftwf_execute == NULL
	|| d->fftwf_execute_dft_r2c == NULL || d->fftwf_execute_dft_c2r == NULL
	// required for internal threading
//	|| d->fftwf_init_threads == NULL || d->fftwf_plan_with_nthreads == NULL
//	|| d->fftwf_cleanup_threads == NULL || d->fftwf_cleanup == NULL


)
ok = false;
	
// 
/*
code for fft dll internal threading
numCPU = 1;
int nThreadMult = 1;	// one thread per cpu

SYSTEM_INFO sysinfo;
GetSystemInfo(&sysinfo);

numCPU = sysinfo.dwNumberOfProcessors;
// check number of processors

d->fftwf_init_threads();

d->fftwf_plan_with_nthreads(nThreadMult * numCPU);
*/
#else // STATIC_FFTW3
	d->hinstLib = NULL;

	d->fftwf_free = (fftwf_free_proc)fftwf_free;
	d->fftwf_malloc = (fftwf_malloc_proc)fftwf_malloc;
	// for 1d fft
	d->fftwf_plan_dft_r2c_1d = (fftwf_plan_dft_r2c_1d_proc)fftwf_plan_dft_r2c_1d;
	d->fftwf_plan_dft_c2r_1d = (fftwf_plan_dft_c2r_1d_proc)fftwf_plan_dft_c2r_1d;
	// for 2d fft
	d->fftwf_plan_dft_r2c_2d = (fftwf_plan_dft_r2c_2d_proc)fftwf_plan_dft_r2c_2d;
	d->fftwf_plan_dft_c2r_2d = (fftwf_plan_dft_c2r_2d_proc)fftwf_plan_dft_c2r_2d;

	d->fftwf_destroy_plan = (fftwf_destroy_plan_proc)fftwf_destroy_plan;
	d->fftwf_execute = (fftwf_execute_proc)fftwf_execute;

	d->fftwf_execute_dft_r2c = (fftwf_execute_dft_r2c_proc)fftwf_execute_dft_r2c;
	d->fftwf_execute_dft_c2r = (fftwf_execute_dft_c2r_proc)fftwf_execute_dft_c2r;

	bool ok = true;
#endif

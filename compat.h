// Unix compatibility layer.
#include <dlfcn.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

typedef void *HINSTANCE;
#define TEXT(x) x

#ifdef __APPLE__
#define SO_SUFFIX ".dylib"
#else
#define SO_SUFFIX ".so"
#endif

static inline void *LoadLibrary(const char *name) {
	if (strcmp(name, "libfftw3f-3.dll") == 0)
		name = "libfftw3f" SO_SUFFIX;
	else if (strcmp(name, "fftw.dll") == 0)
		name = "libfftw3" SO_SUFFIX;
	else
		assert(0 &&  "unknown dll");
	return dlopen(name, RTLD_LAZY | RTLD_LOCAL);
}

static inline void *GetProcAddress(void *handle, const char *name) {
	return dlsym(handle, name);
}
static inline void FreeLibrary(void *handle) {
	dlclose(handle);
}

static inline void *_mm_malloc(int size, int align) {
	void *ptr = NULL;
	if (posix_memalign(&ptr, align, size))
		return NULL;
	return ptr;
}

static inline void _mm_free(void *p) {
	free(p);
}

int fopen_s(FILE **fp, const char *filename, const char *mode) {
	*fp = fopen(filename, mode);
	return *fp == NULL;
}

int strcpy_s(char *dest, int sz, const char *src) {
	assert(strlen(src) < sz);
	strcpy(dest, src);
	return 0;
}

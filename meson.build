project('vcmod', 'cpp',
  default_options: ['buildtype=release', 'b_ndebug=if-release', 'c_std=c99', 'cpp_std=c++11'],
  meson_version: '>=0.48.0'
)

sources = [
  'vcm.cpp',
]

add_project_arguments('-DSTATIC_FFTW3', language: 'cpp')

if build_machine.system() == 'windows'
  fftw3f_dep = dependency('fftwf')
  deps = [fftw3f_dep]
  install_dir = 'installed'  # dummy
else
  vapoursynth_dep = dependency('vapoursynth').partial_dependency(compile_args: true, includes: true)
  fftw3f_dep = dependency('fftw3f')
  deps = [vapoursynth_dep, fftw3f_dep]
  install_dir = join_paths(vapoursynth_dep.get_pkgconfig_variable('libdir'), 'vapoursynth')
  add_project_arguments('-include', 'compat.h', language: 'cpp')
endif

shared_module('vcmod', sources,
  dependencies: deps,
  #install: true,
  #install_dir: install_dir,
  gnu_symbol_visibility: 'hidden'
)

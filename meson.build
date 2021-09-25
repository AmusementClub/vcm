project('vcmod', 'cpp',
  default_options: ['buildtype=release', 'b_ndebug=if-release', 'c_std=c99', 'cpp_std=c++11'],
  meson_version: '>=0.48.0'
)

add_project_arguments('-include', 'compat.h', language: 'cpp')

sources = [
  'vcm.cpp',
]

#vapoursynth_dep = dependency('vapoursynth').partial_dependency(compile_args: true, includes: true)

shared_module('vcmod', sources,
  #dependencies: vapoursynth_dep,
  #install: true,
  #install_dir: join_paths(vapoursynth_dep.get_pkgconfig_variable('libdir'), 'vapoursynth'),
  gnu_symbol_visibility: 'hidden'
)

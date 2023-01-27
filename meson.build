project(
	'crp',
	'cpp',
	version: '0.1',
	default_options: [
		'cpp_std=c++2a',
		'werror=false',
	],
    meson_version : '>= 0.52.0',
)

# for simd instructions
add_project_arguments('-march=native', language: 'cpp')

lib_includes = include_directories('.')

if get_option('enable_kahip')
  add_project_arguments('-DENABLE_KAHIP=1', language: 'cpp')

  kahip_inc = include_directories('./lib/kahip/KaHIP/installed/include')
  kahip_lib = meson.get_compiler('cpp').find_library('kahip', dirs: meson.source_root() + '/lib/kahip/KaHIP/installed/lib')
  kahip = declare_dependency(include_directories: kahip_inc, dependencies: kahip_lib)
else
  kahip = dependency()
  add_project_arguments('-DENABLE_KAHIP=0', language: 'cpp')
endif

subdir('src')
subdir('tests')

executable('sp', 'src/main_sp.cpp', dependencies: crp_dep)

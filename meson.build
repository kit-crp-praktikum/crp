project(
	'crp',
	'cpp',
	version: '0.1',
	default_options: [
		'cpp_std=c++2a',
		'werror=false',
	],
)

lib_includes = include_directories('.')

subdir('src')
subdir('tests')

executable('sp', 'src/main_sp.cpp')
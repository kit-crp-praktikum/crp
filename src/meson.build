# We build the project in two components:
# - CRP library: includes everything about the algorithm except a main() function and CLI interface
# - CRP binary: main function and CLI interface, TBD later.
project_inc_dirs = include_directories('.')

# for simd instructions
add_project_arguments('-march=native', language: 'cpp')

if get_option('enable_print_edges')
  message('enable print edges')
  add_project_arguments('-DPRINT_EDGES=1', language: 'cpp')
else
  message('disable print edges')
endif

# List of source files for the crp library
lib_source_files = [
    'graph.cpp',
    'partitioner/bfs-partitioner.cpp',
    'partitioner/kahip-wrapper.cpp',
    'partitioner/geo-data.cpp',
    'crp/crp-main.cpp',
    'crp/crp-prepare.cpp',
    'crp/crp-query.cpp',
    'crp/crp-customization.cpp',
    'path-unpacker.cpp',
]

crp_lib_includes = [project_inc_dirs, lib_includes]

# Declare omp
omp = dependency('openmp')

# Declare that we want to create a library named 'crp'
crp_library = library('crp', lib_source_files,
    include_directories: crp_lib_includes, dependencies: [omp, kahip])

# A dependency object for our library, makes it easy to use for tests and executables
crp_dep = declare_dependency(link_with: crp_library,
    include_directories: crp_lib_includes)

# Main executable
executable('crp', 'main.cpp', dependencies: [crp_dep, omp])
executable('crp-nora', 'main-partitioner.cpp', dependencies: crp_dep)
executable('drank-gen', 'dijkstra-rank-generator.cpp', dependencies: crp_dep)
executable('partition-extractor', 'partition-extractor.cpp', dependencies: crp_dep)
executable('generate-partition-from-punch', 'punch-partition.cpp', dependencies: crp_dep)

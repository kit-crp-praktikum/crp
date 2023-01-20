# We build the project in two components:
# - CRP library: includes everything about the algorithm except a main() function and CLI interface
# - CRP binary: main function and CLI interface, TBD later.
project_inc_dirs = include_directories('.')


# List of source files for the crp library
lib_source_files = [
    'graph.cpp',
    'partitioner/bfs-partitioner.cpp',
    'partitioner/geo-data.cpp',
    'crp/crp-main.cpp',
    'crp/crp-prepare.cpp',
    'crp/crp-query.cpp',
    'crp/crp-customization.cpp',
    'path-unpacker.cpp',
]

crp_lib_includes = [project_inc_dirs, lib_includes]

# Declare that we want to create a library named 'crp'
crp_library = library('crp', lib_source_files,
    include_directories: crp_lib_includes)

# A dependency object for our library, makes it easy to use for tests and executables
crp_dep = declare_dependency(link_with: crp_library,
    include_directories: crp_lib_includes)

# Main executable
executable('crp', 'main.cpp', dependencies: crp_dep)
executable('crp-nora', 'main-partitioner.cpp', dependencies: crp_dep)
executable('drank-gen', 'dijkstra-rank-generator.cpp', dependencies: crp_dep)

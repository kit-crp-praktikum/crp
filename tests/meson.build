doctest = dependency('doctest')

whole_algorithm = executable('whole-algorithm-test',
    ['whole-algorithm-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that algorithms have the same output', whole_algorithm)

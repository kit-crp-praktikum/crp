doctest = dependency('doctest')

whole_algorithm = executable('whole-algorithm-test',
    ['whole-algorithm-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that algorithms have the same output', whole_algorithm)

overlay_structure = executable('os-test',
    ['os-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that the overlay structure works', overlay_structure)


customization = executable('cust-test',
    ['cust-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that the customization works', customization)

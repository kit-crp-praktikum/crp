doctest = dependency('doctest')

basic_algos = executable('basic-algorithms-test',
    ['basic-algorithms-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that basic algorithms work correctly', basic_algos)

whole_algorithm = executable('whole-algorithm-test',
    ['whole-algorithm-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that CRP works correctly', whole_algorithm)

overlay_structure = executable('os-test',
    ['os-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that the overlay structure works', overlay_structure)

customization = executable('cust-test',
    ['cust-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that the customization works', customization)

path_unpack = executable('path-unpacking',
    ['path-unpacking-test.cpp'], dependencies: [crp_dep, doctest])
test('Test that basic path unpacking works', path_unpack)

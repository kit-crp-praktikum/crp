#pragma once
#include <cstdlib>
#include <filesystem>
#include <iostream>

inline int parse_integer_or_bail(char *param)
{
    char *end;
    int result = strtol(param, &end, 10);
    if (*end != '\0')
    {
        std::cout << "Invalid number given: \"" << param << "\"" << std::endl;
        std::exit(-1);
    }

    return result;
}

inline double parse_double_or_bail(char *param)
{
    char *end;
    double result = strtod(param, &end);
    if (*end != '\0')
    {
        std::cout << "Invalid double given: \"" << param << "\"" << std::endl;
        std::exit(-1);
    }

    return result;
}

inline void check_file_exists_or_bail(std::string file)
{
    if (!std::filesystem::exists(file))
    {
        std::cout << "Missing graph data file " << file << std::endl;
        std::exit(-1);
    }
}

inline void check_input_directory_valid(std::string data_dir, std::string weight, bool requires_geo_data)
{
    auto root = std::filesystem::path{data_dir};
    check_file_exists_or_bail(root / "first_out");
    check_file_exists_or_bail(root / "head");
    check_file_exists_or_bail(root / weight);
    if (requires_geo_data)
    {
        check_file_exists_or_bail(root / "latitude");
        check_file_exists_or_bail(root / "longitude");
    }
}

inline int find_argument_index(int argc, char **argv, char short_arg, std::string long_arg)
{
    std::string sh = std::string("-") + short_arg;
    std::string ln = "--" + long_arg;

    for (int i = 1; i < argc; i++)
    {
        if (argv[i] == sh || argv[i] == ln)
        {
            return i;
        }
    }

    return -1;
}

// Find argument which MUST be provided
inline int find_required_argument(int argc, char **argv, char short_arg, std::string long_arg, bool requires_value)
{
    auto pos = find_argument_index(argc, argv, short_arg, long_arg);

    if (pos == -1)
    {
        std::cout << "Need to specify " << long_arg << std::endl;
        std::exit(-1);
    }

    if (requires_value && pos == argc - 1)
    {
        std::cout << "Missing value for argument " << long_arg << std::endl;
        std::exit(-1);
    }

    return pos;
}

#pragma once

#include "src/data-types.h"
#include <string>
#include <vector>

namespace partitioner
{
class GeoData
{
  public:
    GeoData();

    // Load a lat/lon from files in the format of the praktikum.
    GeoData(std::string latitude, std::string longitude);

    // Load a lat/lon from a vector, useful for simpler tests.
    GeoData(std::vector<int> latitude, std::vector<int> longitude);

  public:
    std::vector<int> latitude;
    std::vector<int> longitude;
};
} // namespace partitioner
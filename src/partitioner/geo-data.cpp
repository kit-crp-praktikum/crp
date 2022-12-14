#include "geo-data.h"
#include "src/data-types.h"
#include "lib/vector_io.h"

partitioner::GeoData::GeoData() {}

partitioner::GeoData::GeoData(std::string latitude, std::string longitude)
{
    this->latitude = load_vector<int>(latitude);
    this->longitude = load_vector<int>(longitude);
}

partitioner::GeoData::GeoData(
    std::vector<int> lat, std::vector<int> lon) : latitude(lat), longitude(lon) {}
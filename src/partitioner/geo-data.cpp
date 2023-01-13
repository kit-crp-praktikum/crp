#include "geo-data.h"
#include "lib/vector_io.h"
#include "src/data-types.h"

partitioner::GeoData::GeoData()
{
}

partitioner::GeoData::GeoData(std::string latitude, std::string longitude)
{
    this->latitude = load_vector<float>(latitude);
    this->longitude = load_vector<float>(longitude);
}

partitioner::GeoData::GeoData(std::vector<float> lat, std::vector<float> lon) : latitude(lat), longitude(lon)
{
}

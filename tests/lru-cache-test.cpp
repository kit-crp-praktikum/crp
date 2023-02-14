#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "crp/lru.hpp"

TEST_CASE("Test LRU")
{
    crp::LRUCache<int, int> cache(3);

    REQUIRE(!cache.get_value(1));
    cache.push_value(1, 1);
    cache.push_value(2, 4);
    cache.push_value(3, 9);

    REQUIRE(*cache.get_value(3) == 9);
    REQUIRE(*cache.get_value(1) == 1);
    REQUIRE(*cache.get_value(2) == 4);

    cache.push_value(4, 16);

    REQUIRE(*cache.get_value(1) == 1);
    REQUIRE(*cache.get_value(2) == 4);
    REQUIRE(!cache.get_value(3)); // 3 was LRU element, so it should have been removed
    REQUIRE(*cache.get_value(4) == 16);
}

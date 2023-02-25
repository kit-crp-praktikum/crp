#pragma once

#include <cassert>
#include <iostream>
#include <list>
#include <map>
#include <optional>

namespace crp
{

/**
 * Implements a cache which maps keys to values using the LRU strategy.
 * This means that the cache is initialized with a maximum size N.
 * When the N+1'th element is added, the element which was least recently used is removed.
 *
 * Implementation strategy: keep all values in a list sorted by LRU order, and also a map from keys to list
 * iterators.
 */
template <class Key, class Value> class LRUCache
{
  public:
    LRUCache(std::size_t size_limit) : limit(size_limit)
    {
        assert(size_limit > 0);
    }

    ~LRUCache()
    {
        if (total_queries > 0)
        {
            std::cerr << "lru_cache_entries=" << indices.size() << std::endl;
            std::cerr << "lru_cache_hits=" << hits << std::endl;
            std::cerr << "lru_cache_queries=" << total_queries << std::endl;
        }
    }

    Value *get_value(const Key &key)
    {
        ++total_queries;
        auto it = indices.find(key);
        if (it != indices.end())
        {
            ++hits;
            auto list_it = it->second;
            // Push to the front
            if (list_it != values.begin())
            {
                values.splice(values.begin(), values, list_it, std::next(list_it));
            }

            return &list_it->first;
        }

        return nullptr;
    }

    void push_value(const Key &key, Value &&value)
    {
        if (values.size() == limit)
        {
            remove(std::prev(values.end()));
        }

        values.insert(values.begin(), {value, key});
        indices[key] = values.begin();
    }

  private:
    using ListType = std::list<std::pair<Value, Key>>;
    using ListIterator = typename ListType::iterator;
    size_t hits = 0;
    size_t total_queries = 0;

    void remove(ListIterator it)
    {
        indices.erase(it->second);
        values.erase(it);
    }

    ListType values;
    std::map<Key, ListIterator> indices;
    std::size_t limit;
};

} // namespace crp

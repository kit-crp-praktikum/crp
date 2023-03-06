#pragma once

#include <cassert>
#include <iostream>
#include <list>
#include <map>
#include <optional>
#include <vector>

namespace crp
{
namespace detail
{
template <class T> size_t estimate_memory(const T &t)
{
    return sizeof(std::decay_t<T>);
}

template <class T> size_t estimate_memory(const std::vector<T> &t)
{
    return t.size() * sizeof(std::decay_t<T>);
}
} // namespace detail

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
            std::cerr << "lru_cache_limit=" << limit << std::endl;
            std::cerr << "lru_cache_entries=" << indices.size() << std::endl;
            std::cerr << "lru_cache_hits=" << hits << std::endl;
            std::cerr << "lru_cache_queries=" << total_queries << std::endl;

            size_t total_size = 0;
            for (const auto &it : values)
            {
                total_size += detail::estimate_memory(it.first);
                total_size += detail::estimate_memory(it.second);
            }

            std::cerr << "lru_cache_memory=" << total_size << std::endl;
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
        if (auto existing = get_value(key))
        {
            *existing = value;
            return;
        }

        assert(indices.count(key) == 0);
        if (values.size() == limit)
        {
            remove(std::prev(values.end()));
        }

        values.insert(values.begin(), {value, key});
        indices[key] = values.begin();
    }

    void reset_stats()
    {
        hits = 0;
        total_queries = 0;
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

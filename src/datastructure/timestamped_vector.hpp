#pragma once

#include <cstdint>
#include <limits>
#include <vector>

/**
 * A vector which keeps timestamps for each entry.
 * Entries which are not from the current timestamp will return a default value.
 * Reset in O(1), by incrementing the current timestamp.
 *
 * Usage:
 * int n = 100, def = 1;
 * TimestampedVector<int> v(n, def);
 * assert(v[0] == def);
 * v[0] = 2;
 * assert(v[0] == 2);
 * v.reset();
 * assert(v[0] == def);
 */
template <class T> class TimestampedVector
{
    using Time = uint32_t;

  public:
    TimestampedVector(std::size_t size, T default_value)
        : data(size, default_value), timestamps(size, 0), default_value(default_value), current_time(1)
    {
    } // initially everything is empty

    std::size_t size() const
    {
        return data.size();
    }

    // maybe a bit faster than [], because no branching is involved
    void set(std::size_t index, T value)
    {
        data[index] = value;
        timestamps[index] = current_time;
    }

    // increment current time. in case of overflow, reset everything.
    void reset()
    {
        if (current_time == std::numeric_limits<Time>::max())
        {
            current_time = 0;
            std::fill(data.begin(), data.end(), default_value);
        }
        else
        {
            current_time++;
        }
    }

    typename std::vector<T>::reference operator[](std::size_t index)
    {
        if (timestamps[index] != current_time)
        {
            data[index] = default_value;
            timestamps[index] = current_time;
        }
        return data[index];
    }

  private:
    std::vector<T> data;
    std::vector<Time> timestamps;
    T default_value;
    Time current_time;
};

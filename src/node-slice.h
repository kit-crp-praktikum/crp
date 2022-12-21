#pragma once

#include "data-types.h"
#include <iterator>
#include <utility>

namespace detail
{
template <class Graph> class NodeSlice
{
  public:
    NodeSlice(const Graph *g_, NodeId u_) : g(g_), u(u_)
    {
    }

    class NodeSliceIterator
    {
        const Graph *g;
        int32_t pos;

      public:
        using iterator_category = std::output_iterator_tag;

        NodeSliceIterator(const Graph *g_, int p) : g(g_), pos(p)
        {
        }
        NodeSliceIterator operator++(int) /* postfix */
        {
            return NodeSliceIterator(g, pos++);
        }

        NodeSliceIterator &operator++() /* prefix */
        {
            ++pos;
            return *this;
        }

        std::pair<NodeId, Distance> operator*() const
        {
            return {g->head[pos], g->weights[pos]};
        }

        bool operator==(const NodeSliceIterator &rhs) const
        {
            return pos == rhs.pos;
        }

        bool operator!=(const NodeSliceIterator &rhs) const
        {
            return pos != rhs.pos;
        }
    };

    using iterator = NodeSliceIterator;
    using const_iterator = NodeSliceIterator;

    iterator begin()
    {
        return NodeSliceIterator(g, g->first_out[u]);
    }
    const_iterator begin() const
    {
        return NodeSliceIterator(g, g->first_out[u]);
    }
    iterator end()
    {
        return NodeSliceIterator(g, g->first_out[u + 1]);
    }
    const_iterator end() const
    {
        return NodeSliceIterator(g, g->first_out[u + 1]);
    }

  private:
    const Graph *g;
    NodeId u;
};
} // namespace detail

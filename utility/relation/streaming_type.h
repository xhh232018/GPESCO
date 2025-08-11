#ifndef RAPIDMATCH_STREAMING_TYPE_H
#define RAPIDMATCH_STREAMING_TYPE_H

#include <cstdint>
#include <vector>

#include "../sparsepp/spp.h"
#include "../han/utils/util.hpp"
/**
 * begin vertex, end vertex
 */

/**
 * Edges and the corresponding view.
 */
typedef std::pair<std::vector<Edge>, uint32_t> MappedViews;

namespace std
{
    template<>
    struct hash<Edge>{
        std::size_t operator()(Edge const &l) const{
            std::size_t seed = 0;
            spp::hash_combine(seed, l.first);
            spp::hash_combine(seed, l.second);
            return seed;
        }
    };

    template<>
    struct hash<pEdge>{
        std::size_t operator()(pEdge const &l) const{
            std::size_t seed = 0;
            spp::hash_combine(seed, l.first);
            spp::hash_combine(seed, l.second);
            return seed;
        }
    };

}

#endif //RAPIDMATCH_STREAMING_TYPE_H

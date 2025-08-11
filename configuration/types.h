#ifndef SUBGRAPHMATCHING_TYPES_H
#define SUBGRAPHMATCHING_TYPES_H

#include <cstdint>
#include <stdlib.h>
#include <queue>

#include "../utility/sparsepp/spp.h"

typedef unsigned int ui;

typedef uint32_t VertexID;
typedef ui LabelID;

typedef std::pair<uint32_t, uint32_t> pEdge;

namespace std
{
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

enum MatchingIndexType {
    VertexCentric = 0,
    EdgeCentric = 1
};

enum RelationStructure{
    EncodedTrieRelation,
    HashRelation,
    TrieRelation
};

struct edge {
    uint32_t vertices_[2];
};

class TreeNode {
public:
    VertexID id_;
    VertexID parent_;
    ui level_;
    ui under_level_count_;
    ui children_count_;
    ui bn_count_;
    ui fn_count_;
    VertexID* under_level_;
    VertexID* children_;
    VertexID* bn_;
    VertexID* fn_;
    size_t estimated_embeddings_num_;
public:
    TreeNode() {
        id_ = 0;
        under_level_ = NULL;
        bn_ = NULL;
        fn_ = NULL;
        children_ = NULL;
        parent_ = 0;
        level_ = 0;
        under_level_count_ = 0;
        children_count_ = 0;
        bn_count_ = 0;
        fn_count_ = 0;
        estimated_embeddings_num_ = 0;
    }

    ~TreeNode() {
        delete[] under_level_;
        delete[] bn_;
        delete[] fn_;
        delete[] children_;
    }

    void initialize(const ui size) {
        under_level_ = new VertexID[size];
        bn_ = new VertexID[size];
        fn_ = new VertexID[size];
        children_ = new VertexID[size];
    }
};

class Edges {
public:
    ui* offset_;
    ui* edge_;
    ui vertex_count_;
    ui edge_count_;
    ui max_degree_;
public:
    Edges() {
        offset_ = NULL;
        edge_ = NULL;
        vertex_count_ = 0;
        edge_count_ = 0;
        max_degree_ = 0;
    }

    ~Edges() {
        delete[] offset_;
        delete[] edge_;
    }
};

struct result_struct {
    uint64_t inc_cnt = 0;
    uint32_t head;
    uint32_t tail;

    result_struct(uint64_t ic, uint32_t hd, uint32_t tl) {
        inc_cnt = ic;
        head = hd;
        tail = tl;
    }
};

class res_comp {
public:
    bool operator()(const result_struct &r1, const result_struct &r2) {
        return r1.inc_cnt > r2.inc_cnt;
    }
};

typedef std::priority_queue<result_struct, std::vector<result_struct>, res_comp> Res_PQ;

#endif //SUBGRAPHMATCHING_TYPES_H

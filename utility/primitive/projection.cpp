#include "projection.h"
#include <cstring>
#include <algorithm>
void projection::execute(edge_relation *relation, uint32_t kp, uint32_t* &res, uint32_t &res_cnt) {
    edge* edges = relation->edges_;
    uint32_t edge_size = relation->size_;

    memset(index_, 0, sizeof(bool) * index_size_);
    uint32_t cnt = 0;
    for (uint32_t i = 0; i < edge_size; ++i) {
        uint32_t u = edges[i].vertices_[kp];

        if (!index_[u]) {
            index_[u] = true;
            buffer_[cnt++] = u;
        }
    }

    std::sort(buffer_, buffer_ + cnt);

    res = new uint32_t[cnt];
    memcpy(res, buffer_, sizeof(uint32_t) * cnt);
    res_cnt = cnt;
}

void projection::store_res(uint32_t* &res, uint32_t &res_cnt){

    res = new uint32_t[init_size_];
    memcpy(res, buffer_, sizeof(uint32_t) * init_size_);
    res_cnt = init_size_;

}

void projection::init(edge_relation *relation, uint32_t kp){
    edge* edges = relation->edges_;
    uint32_t edge_size = relation->size_;

    memset(index_, 0, sizeof(bool) * index_size_);
    uint32_t cnt = 0;
    for (uint32_t i = 0; i < edge_size; ++i) {
        uint32_t u = edges[i].vertices_[kp];

        if (!index_[u]) {
            index_[u] = true;
            buffer_[cnt++] = u;
        }
    }

    std::sort(buffer_, buffer_ + cnt);

    this->init_size_ = cnt;
}

void projection::update(edge_relation *relation, uint32_t kp){
    edge* edges = relation->edges_;
    uint32_t edge_size = relation->size_;

    memset(index_curr_, 0, sizeof(bool) * index_size_);
    uint32_t cnt = 0, valid_cnt = 0;

    for (uint32_t i = 0; i < edge_size; ++i) {
        uint32_t u = edges[i].vertices_[kp];

        if (!index_curr_[u]) {
            index_curr_[u] = true;
        }
    }

    //update state
    for (uint32_t i = 0; i < init_size_; ++i) {
        uint32_t candidate = buffer_[i];
        if(index_[candidate] && index_curr_[candidate]){
            buffer_[valid_cnt++] = candidate;
        }
    }

    init_size_ = valid_cnt;
}
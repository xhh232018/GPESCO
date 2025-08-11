#ifndef SUBGRAPHMATCHING_PROJECTION_H
#define SUBGRAPHMATCHING_PROJECTION_H


#include <cstdint>
#include <relation/edge_relation.h>
/// Will sort the result after projection.
class projection {
private:
    bool* index_;
    bool* index_curr_;
    uint32_t* buffer_;
    uint32_t index_size_;
    uint32_t init_size_;

public:
    projection(uint32_t index_size) {
        index_size_ = index_size;
        buffer_ = new uint32_t[index_size_];
        index_ = new bool[index_size_];
        index_curr_ = new bool[index_size_];
    }
    ~projection() {
        delete[] buffer_;
        delete[] index_;
        delete[] index_curr_;
    }

    void execute(edge_relation *relation, uint32_t kp, uint32_t* &res, uint32_t& res_cnt);

    void init(edge_relation *relation, uint32_t kp);

    void update(edge_relation *relation, uint32_t kp);

    void store_res(uint32_t* &res, uint32_t& res_cnt);

};


#endif //SUBGRAPHMATCHING_PROJECTION_H

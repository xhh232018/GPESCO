#ifndef SUBGRAPHMATCHING_ENCODERSP_H
#define SUBGRAPHMATCHING_ENCODERSP_H

#include <relation/catalog.h>

class encoderSP {
public:
    double encoding_time_;
    double build_relation_time_;
    double projection_time_;
    double convert_to_sparse_bp_time_;
    uint32_t core_size;

private:
    uint32_t max_vertex_id_;
    uint32_t* temp_buffer1_;
    uint32_t* temp_buffer2_;
    //uint32_t* candidates_buffer;
    const Graph* query_graph_;
    uint32_t* join_plan_;
    bool** local_zero_deg;

private:
    void convert_to_trie_relation(catalog *storage);
    void convert_to_trie_relation(catalog *storage, uint32_t u, uint32_t v);

    void candidate_space_shrinking(catalog *storage);
    void space_init_proj(catalog *storage);

    void convert_to_hash_relation(catalog *storage);
    void convert_to_hash_relation(catalog *storage, uint32_t u, uint32_t v);
    

    void convert_to_encoded_relation(catalog *storage);
    void convert_to_encoded_relation(catalog *storage, uint32_t u, uint32_t v);
    void convert_to_encoded_relation(catalog *storage, uint32_t u, uint32_t v, uint32_t index);

    void convert_trie_relation_to_sparse_bitmap(catalog *storage);
    void convert_encoded_relation_to_sparse_bitmap(catalog *storage);

    void convert_hash_relation_to_sparse_bitmap(catalog *storage);

public:
    encoderSP(const Graph* query_graph, uint32_t max_vertex_id) :
            encoding_time_(0), build_relation_time_(0), projection_time_(0), convert_to_sparse_bp_time_(0),
            max_vertex_id_(max_vertex_id) {

        query_graph_ = query_graph;
        temp_buffer1_ = new uint32_t[max_vertex_id_];
        temp_buffer2_ = new uint32_t[max_vertex_id_];
        //candidates_buffer = new uint32_t[max_vertex_id_];
        core_size = query_graph->get2CoreSize();
        local_zero_deg = new bool*[core_size];
        for(int i = 0; i < core_size;++i){
            local_zero_deg[i] = new bool[max_vertex_id_];
            memset(local_zero_deg[i], 0, sizeof(bool) * max_vertex_id_);
        }
        memset(temp_buffer1_, 0, max_vertex_id_ * sizeof(uint32_t));
    }

    ~encoderSP() {
        delete[] temp_buffer1_;
        delete[] temp_buffer2_;
        for(int i = 0; i < core_size;++i){
            delete[] local_zero_deg[i];
        }
        delete[] local_zero_deg;
        //delete[] candidates_buffer;
    }

    void refine_candidate_set(catalog *storage);

    void execute_iso(catalog *storage, RelationStructure relation_type, bool enable_sparsebp, uint32_t *join_plan);

    void execute(catalog *storage, RelationStructure relation_type, bool enable_sparsebp, uint32_t *join_plan);

    void execute_with_shrinking(catalog *storage, RelationStructure relation_type, bool enable_sparsebp, uint32_t *join_plan, int k);

    void execute(catalog *storage, RelationStructure relation_type, bool enable_sparsebp, uint32_t *join_plan, uint32_t u, uint32_t v);

    void execute_nc(catalog *storage, RelationStructure relation_type, bool enable_sparsebp, uint32_t *join_plan);
    
    void print_metrics();
};


#endif //SUBGRAPHMATCHING_ENCODER_H

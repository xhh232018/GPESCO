
#ifndef SUBGRAPHMATCHING_LEAPFROGTRIEJOIN_H
#define SUBGRAPHMATCHING_LEAPFROGTRIEJOIN_H

#include "../relation/catalog.h"
#include "../relation/flat_relation.h"
#include <bitset>
class leapfrogtriejoin {

#ifdef ENABLE_OUTPUT
    uint32_t* current_position_;
#endif

public:
    uint64_t* set_intersection_count_;
    uint64_t* set_intersection_cost_;
    uint64_t* intermediate_result_count_;
    uint64_t fail_si_count_;
    uint64_t iso_conflict_count_;
    uint64_t invalid_core_pr_count_;
    uint64_t invalid_leaf_pr_count_;
    uint64_t* enter_result_count_;
    uint64_t* enter_ptb_count_;
    spp::sparse_hash_map<pEdge,uint64_t> res_tmp;

    double prob = 1;
private:
    flat_relation* input_;
    flat_relation* output_;
    uint32_t* vertex_ordering_;
    uint32_t num_vertex_;
    uint32_t num_core_vertex_;
    uint32_t lc_exist_count;
    catalog* catalog_;
    uint64_t output_count_limit_;
    bool materialize_;
    bool begin_node_first;
    uint64_t count_;

    uint32_t* idx_;
    uint32_t* num_local_candidates_;
    uint32_t* num_ptb_candidates_;
    uint32_t* idx_embedding_;
    uint32_t* embedding_;
    uint32_t* si_buffer_;
    uint32_t* si_buffer_precomputed;
    uint32_t** local_candidates_;
    uint32_t** ptb_candidates_;
    uint32_t* last_level_candidates_;

    uint32_t* bns_offset_;
    uint32_t* bns_;
    uint32_t* bns_depth_;
    uint32_t* fns_offset_;
    uint32_t* fns_;
    uint32_t* fns_depth_;

    std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> end_node_bn;

#ifdef INTERSECTION_CACHE
    enum INTERSECTION_CACHE_TYPE {
        none = 0, full = 1, partial = 2
    };
    uint32_t* cached_bns_offset_;
    uint32_t* cached_bns_;
    uint32_t* cached_bns_depth_;
    uint32_t* no_cached_bns_offset_;
    uint32_t* no_cached_bns_;
    uint32_t* no_cached_bns_depth_;
    uint32_t* mapping_to_cached_bns_;
    uint32_t** cached_candidates_;
    uint32_t* cached_candidates_count_;
    INTERSECTION_CACHE_TYPE* intersection_cache_type_;
#endif

#ifdef FAILING_SET_PRUNING
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors_;
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set_;
    std::unordered_map<VertexID, VertexID> reverse_embedding_;
#endif

#ifdef SPARSE_BITMAP
    BSRSet* cached_bsr_sets_;
    BSRSet buffer_bsr_set1_;
    BSRSet buffer_bsr_set2_;
#endif

#ifndef HOMOMORPHISM
    bool* visited_;
#endif

private:
#if RELATION_STRUCTURE == 0
    void set_idx(uint32_t depth);
#endif

    uint32_t compute_local_candidates(uint32_t depth);

    uint32_t pre_compute_ptb_candidates(uint32_t depth);

    uint32_t pre_compute_ptb_candidates_nc(uint32_t depth);

    uint32_t compute_local_candidates_core_ptb(uint32_t depth);

    uint32_t compute_local_candidates_noncore_ptb(uint32_t bn, uint32_t u, uint32_t key,uint32_t depth);

    uint32_t compute_local_candidates_core_ptb_nc(uint32_t depth);

    uint32_t compute_local_candidates_noncore_ptb_nc(uint32_t depth);

    void initialize();

    void initialize_bn_fn();

    void probe();

    bool exist_check(uint32_t u, uint32_t *candidates, uint32_t count);

    bool enumerate_non_core_results_lae(uint32_t start_depth, uint32_t max_depth);

    bool enumerate_non_core_results_iee(uint32_t start_depth, uint32_t max_depth);

    bool enumerate_non_core_results_je(uint32_t start_depth, uint32_t max_depth);

    bool enumerate_non_core_results_drq_single(uint32_t start_depth, uint32_t max_depth);

    bool enumerate_non_core_results_probe_mat(uint32_t start_depth, uint32_t max_depth);

    bool enumerate_non_core_results_pe(uint32_t start_depth, uint32_t max_depth);

    bool enumerate_non_core_results_rwc(uint32_t start_depth, uint32_t max_depth);

    bool rwd_non_core(uint32_t start_depth, uint32_t max_depth);

    bool rwch_core(uint32_t start_depth, uint32_t max_depth);

    bool rwch_non_core(uint32_t start_depth, uint32_t max_depth);

    bool rwch_non_core_single(uint32_t start_depth, uint32_t max_depth);

    void rwch_non_core_sampling(uint32_t start_depth, uint32_t max_depth);

    bool rwcw_non_core_epl(uint32_t start_depth, uint32_t max_depth);

    bool rwcw_non_core_final(uint32_t start_depth, uint32_t max_depth);

#ifdef INTERSECTION_CACHE
    void initialize_intersection_cache();

    bool is_intersection_cache_valid(uint32_t depth);
#endif

#ifdef FAILING_SET_PRUNING
    void initialize_failing_set_pruning();

    void compute_ancestors();
#endif

    void clear();

    uint32_t* get_backward_neighbors(uint32_t depth) {
        return bns_ + bns_offset_[depth];
    }

    uint32_t* get_backward_neighbors_depth(uint32_t depth) {
        return bns_depth_ + bns_offset_[depth];
    }

    uint32_t get_backward_neighbors_count(uint32_t depth) {
        return bns_offset_[depth + 1] - bns_offset_[depth];
    }

    uint32_t* get_forward_neighbors(uint32_t depth) {
        return fns_ + fns_offset_[depth];
    }

    uint32_t* get_forward_neighbors_depth(uint32_t depth) {
        return fns_depth_ + fns_offset_[depth];
    }

    uint32_t get_forward_neighbors_count(uint32_t depth) {
        return fns_offset_[depth + 1] - fns_offset_[depth];
    }

public:
    leapfrogtriejoin(flat_relation* input, flat_relation* output, uint32_t* vertex_ordering,
                     uint32_t vertex_count, uint32_t core_vertex_count, catalog* catalog,
                     uint64_t output_count_limit = std::numeric_limits<uint64_t>::max()) :
                             input_(input), output_(output), vertex_ordering_(vertex_ordering),
                             num_vertex_(vertex_count), num_core_vertex_(core_vertex_count), catalog_(catalog),
                             output_count_limit_(output_count_limit), materialize_(false), count_(0) {
        initialize();
        initialize_bn_fn();

#ifdef INTERSECTION_CACHE
        initialize_intersection_cache();
#endif

#ifdef FAILING_SET_PRUNING
        initialize_failing_set_pruning();
#endif
    }

    ~leapfrogtriejoin() {
        clear();
    }

    uint64_t execute();

    uint64_t execute_je();

    uint64_t execute_lae();

    uint64_t execute_iee();

    uint64_t execute_pe();

    uint64_t execute_mat_probe();

    uint64_t execute_rwc();

    uint64_t execute_rwch();

    uint64_t execute_rwd();

    uint64_t execute_rwdw_epl();

    uint64_t execute_rwdw_final();

};


#endif //SUBGRAPHMATCHING_LEAPFROGTRIEJOIN_H

#ifndef SUBGRAPHMATCHING_DRQENUM_H
#define SUBGRAPHMATCHING_DRQENUM_H

#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <computesetintersection.h>

#include "matchingcommand.h"
#include "graph/graph.h"
#include "relation/catalog.h"
#include "execution_tree_generator.h"
#include "pretty_print.h"
#include "global_variables.h"
#include "preprocessor.h"
#include "encoder.h"
#include "query_plan_generator.h"


class drqEnum
{
public:
       /* data */
    std::string tk;
    std::string output_file;
    bool enable_sparsebp;
    bool enable_elimination;
    RelationStructure relation_type;
    uint32_t topk;
    uint64_t time_limit;
    uint64_t output_limit;
    Res_PQ rpq;

public:
   
    drqEnum(std::string tk_, std::string of, bool es, bool ee, RelationStructure rt, uint32_t k, uint64_t tl, uint64_t ol):
    tk(tk_), output_file(of), enable_sparsebp(es), enable_elimination(ee), relation_type(rt), topk(k), time_limit(tl), output_limit(ol) {}

    void getVEList(Graph* query_graph, std::vector<uint32_t> &component, 
                std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> &vertex_list,
                std::vector<std::pair<uint32_t, uint32_t>> &edge_list,
                std::unordered_map<uint32_t,uint32_t> &rvs_map);

    void adjust_order(Graph* large_comp, Graph* small_comp, std::vector<uint32_t> &big_order,  std::vector<uint32_t> &small_order,
                     uint32_t big_depth, uint32_t small_depth, std::vector<uint32_t> &new_order, uint32_t &total_core,
                      std::unordered_map<uint32_t,uint32_t> &rvs_map_big,
                      std::unordered_map<uint32_t,uint32_t> &rvs_map_small);
    
    void bfs(std::vector<std::vector<uint32_t>> &graph, uint32_t S, std::vector<uint32_t>& par, std::vector<uint32_t>& dist);

    void candidates_pair(catalog* storage,std::vector<uint32_t> &new_order);

    void printShortDist(std::vector<std::vector<uint32_t>> &forward_nbr_vec, uint32_t source, uint32_t dest, uint32_t graph_size,std::vector<uint32_t>& path);

    void adj_exe(Graph* comp_graph,std::vector<uint32_t> &old_order,uint32_t ext_depth, std::vector<uint32_t> &new_order,uint32_t &core_num_aft);

    void common_label_depth(Graph* component_graph, catalog* catalog, std::vector<uint32_t> &order);

    void common_label_detection(Graph* query_graph, catalog* storage,std::vector<uint32_t> &comp_small, std::vector<uint32_t> &comp_big);

    void cld_pe(Graph* component_graph, catalog* catalog, std::vector<uint32_t> &order);

    void execute_pe(std::vector<uint32_t> &component_small,std::vector<uint32_t> &component_large,Graph* query_graph, Graph* data_graph, uint32_t u, uint32_t v);

    void execute_iee(std::vector<uint32_t> &component_small,std::vector<uint32_t> &component_large,Graph* query_graph, Graph* data_graph, uint32_t u, uint32_t v);

    void execute_within_time_limit(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit);

    void execute_within_time_limit_iee(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit);

    void execute_within_time_limit_probe_mat(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit);

    void execute_within_time_limit_pe(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit);

    void count_single(Graph* data_graph, uint32_t u_label, catalog* storage);

    void count_product(Graph* data_graph, catalog* storage_small,catalog* storage_big);

    void bound_pruning(Graph* data_graph,std::vector<std::pair<uint32_t,uint32_t>> &candidate_space,
                        catalog* storage_small, catalog* storage_big);

};

#endif
#ifndef SUBGRAPHMATCHING_CRQSAMPLING_H
#define SUBGRAPHMATCHING_CRQSAMPLING_H

#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <computesetintersection.h>
#include <random>
#include <cmath>

#include "matchingcommand.h"
#include "graph/graph.h"
#include "relation/catalog.h"
#include "execution_tree_generator.h"
#include "pretty_print.h"
#include "global_variables.h"
#include "preprocessor.h"
#include "encoderSP.h"
#include "query_plan_generator.h"
#include "utility/leapfrogtriejoin/leapfrogtriejoin.h"  

class crqSampling
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
    crqSampling(std::string tk_, std::string of, bool es, bool ee, RelationStructure rt, uint32_t k, uint64_t tl, uint64_t ol):
    tk(tk_), output_file(of), enable_sparsebp(es), enable_elimination(ee), relation_type(rt), topk(k), time_limit(tl), output_limit(ol) {}
    
    void RWC(Graph* query_graph, Graph* data_graph, uint32_t u, uint32_t v);

    void execute_within_time_limit_rwc(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit);

    void execute_within_time_limit_rwch(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit);

    void exact_result(catalog* storage);

    void exact_result_rwch(catalog* storage);

    void RWCH(Graph* query_graph, Graph* data_graph, uint32_t u, uint32_t v);
   
};

#endif
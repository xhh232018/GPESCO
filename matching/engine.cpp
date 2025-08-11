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
#include "crqEnum.h"
#include "drqEnum.h"
#include "crqSampling.h"
#include "drqSampling.h"

bool checkComponent(Graph* query_graph,std::vector<std::vector<uint32_t>> &connected_components){
    uint32_t n = query_graph->getVerticesCount();
    std::vector<bool> visited(n, false);

    std::queue<uint32_t> q;
    bool flag = false;
    for (uint32_t u = 0; u < n; ++u) {
        if (!visited[u]) {
            std::vector<uint32_t> component;
            component.push_back(u);
            q.push(u);
            visited[u] = true;

            while (!q.empty()) {
                uint32_t uu = q.front();
                q.pop();
                uint32_t uu_nbrs_count;
                auto uu_nbrs = query_graph->getVertexNeighbors(uu, uu_nbrs_count);

                for (uint32_t i = 0; i < uu_nbrs_count; ++i) {
                    uint32_t uuu = uu_nbrs[i];
                    if (!visited[uuu]) {
                        q.push(uuu);
                        component.push_back(uuu);
                        visited[uuu] = true;
                    }
                }
            }

            connected_components.emplace_back(component);
        }
    }
    if(connected_components.size() > 1) flag = true;
    return flag;
}

void loadQueryGraphMeta(const std::string& file_path, ui &num_ptbs, std::vector<std::pair<VertexID,VertexID>> &edgeMap){
    std::ifstream infile(file_path);
    char type;
    ui tmp;
    while (infile >> type) {
        if(type == 't'){
            infile >> tmp >> num_ptbs;
        }
        else if (type == 'e') { // Read edge.
            VertexID begin;
            VertexID end;
            infile >> begin >> end;
            edgeMap.push_back(make_pair(begin,end));
        }
    }
    infile.close();
}

int main(int argc, char** argv) {
    MatchingCommand command(argc, argv);
    std::string input_query_folder = command.getQueryGraphFolderPath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string input_time_limit = command.getTimeLimit();
    std::string input_crq = command.getCRQPath();
    std::string input_drq = command.getDRQPath();
    std::string input_qlist = command.getQList();
    std::string input_enable_preprocessor = command.getPreprocessor();
    std::string input_output_path = command.getOutputPath();
    std::string input_topK = command.getTopK();
    /**
     * Output the command line information.
     */
    std::cout << "Command Line:" << std::endl;
    
    std::cout << "\tData Graph: " << input_data_graph_file << std::endl;
    std::cout << "\tQuery Graph Folder: " << input_query_folder << std::endl;
    std::cout << "\tQuery List File: " << input_qlist << std::endl;
    std::cout << "\tCRQ Algorithm: " << input_crq << std::endl;
    std::cout << "\tDRQ Algorithm: " << input_drq << std::endl;
    std::cout << "\tTime Limit (seconds): " << input_time_limit << std::endl;
    std::cout << "\tEnable Preprocessor: " << input_enable_preprocessor << std::endl;
    std::cout << "\tOutput Path: " << input_output_path << std::endl;
    std::cout << "\ttopK: " << input_topK  << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;

    /**
     * Output the configuration.
     */
    std::cout << "Configuration:" << std::endl;
#if RELATION_STRUCTURE == 0
    RelationStructure relation_type = RelationStructure::EncodedTrieRelation;
    std::cout << "\tRelation Structure: Encoded trie relation" << std::endl;
#elif RELATION_STRUCTURE == 1
    RelationStructure relation_type = RelationStructure::HashRelation;
    std::cout << "\tRelation Structure: Hash relation" << std::endl;
#elif RELATION_STRUCTURE == 2
    RelationStructure relation_type = RelationStructure::TrieRelation;
    std::cout << "\tRelation Structure: Trie relation" << std::endl;
#endif

#ifdef HOMOMORPHISM
    std::cout << "\tEmbedding Structure: Subgraph Homomorphism " << std::endl;
#else
    std::cout << "\tEmbedding Structure: Subgraph Isomorphism " << std::endl;
#endif

#ifdef COUNT_RESULTS
    std::cout << "\tCount Results: True" << std::endl;
#else
    std::cout << "\tCount Results: False" << std::endl;
#endif

#if HYBRID == 0
    std::cout << "\tHybrid Set Intersection: Enable" << std::endl;
#else
    std::cout << "\tHybrid Set Intersection: Disable" << std::endl;
#endif

#if SI == 0
    std::cout << "\tSIMD Set Intersection: AVX2" << std::endl;
#elif SI == 1
    std::cout << "\tSIMD Set Intersection: AVX512" << std::endl;
#elif SI == 2
    std::cout << "\tSIMD Set Intersection: NONE" << std::endl;
#endif

    std::cout << "\tHash Table Ratio of Pair-wise Join: " << HASH_TABLE_RATIO << std::endl;

#ifdef INTERSECTION_CACHE
    std::cout << "\tIntersection Cache: Enabled" << std::endl;
#else
    std::cout << "\tIntersection Cache: Disabled" << std::endl;
#endif

#ifdef FAILING_SET_PRUNING
    std::cout << "\tFailing Set Pruning: Enabled" << std::endl;
#else
    std::cout << "\tFailing Set Pruning: Disabled" << std::endl;
#endif

#ifdef OUTPUT_OPTIMIZATION
    std::cout << "\tOutput Optimization: Enabled" << std::endl;
#else
    std::cout << "\tOutput Optimization: Disabled" << std::endl;
#endif

#ifdef SPARSE_BITMAP
    bool enable_sparsebp = true;
    std::cout << "\tSparse Bitmap: Enabled" << std::endl;
#else
    bool enable_sparsebp = false;
    std::cout << "\tSparse Bitmap: Disabled" << std::endl;
#endif

    bool enable_preprocessor = true;
    if (input_enable_preprocessor == "false") {
        enable_preprocessor = false;
    }
    else if (input_enable_preprocessor == "true") {
        enable_preprocessor = true;
    }
    uint64_t time_limit = 0;
    sscanf(input_time_limit.c_str(), "%zu", &time_limit);

    uint64_t output_limit = std::numeric_limits<uint64_t>::max();

    uint32_t topk = 1;
    sscanf(input_topK.c_str(), "%zu", &topk);

    std::cout << "--------------------------------------------------------------------" << std::endl;
    /**
     * Load input graphs.
     */
    std::cout << "Load graphs..." << std::endl;
    std::ifstream ifs(input_qlist);
    std::string line;
    auto data_graph = new Graph(true);
    data_graph->loadGraphFromFile(input_data_graph_file);


    while (getline(ifs, line)) {
        std::istringstream iss(line);
        std::string tk;
        iss >> tk;
        std::string query_file = input_query_folder + tk;

        ui num_ptbs; 
        std::vector<std::pair<uint32_t,uint32_t>>pEdgeMap;
        loadQueryGraphMeta(query_file,num_ptbs,pEdgeMap);

        auto query_graph_ori = new Graph(true);
        query_graph_ori->loadGraphFromFile(query_file);
        query_graph_ori->buildCoreTable();

        for(ui i = 0; i < num_ptbs;++i){

            VertexID begin_ = pEdgeMap[i].first;
            VertexID end_ = pEdgeMap[i].second;

            auto query_graph = new Graph(true);
            query_graph->constuctPerturbedQueryG(query_file,i,begin_,end_);
            query_graph->buildCoreTable();

            std::vector<std::vector<uint32_t>> component_vec;

            if(checkComponent(query_graph,component_vec)) {//drq
                if(input_drq == "ProbeEnum"){
                    drqEnum drq(tk, input_output_path, enable_sparsebp, enable_preprocessor, relation_type, topk, time_limit, output_limit);
                    if(component_vec.front().size() < component_vec.back().size()) {
                        drq.execute_pe(component_vec.front(), component_vec.back(), query_graph, data_graph, begin_, end_);
                    }
                    else{
                        drq.execute_pe(component_vec.back(), component_vec.front(), query_graph, data_graph, begin_, end_);
                    }
                }
                else if(input_drq == "IEEnum"){
                    drqEnum drq(tk, input_output_path, enable_sparsebp, enable_preprocessor, relation_type, topk, time_limit, output_limit);
                    if(component_vec.front().size() < component_vec.back().size()) {
                        drq.execute_iee(component_vec.front(), component_vec.back(), query_graph, data_graph, begin_, end_);
                    }
                    else{
                        drq.execute_iee(component_vec.back(), component_vec.front(), query_graph, data_graph, begin_, end_);
                    }
                }
                else if(input_drq == "RWD"){
                    drqSampling drq(tk, input_output_path, enable_sparsebp, enable_preprocessor, relation_type, topk, time_limit, output_limit);
                    drq.RWD(query_graph_ori, query_graph, data_graph, begin_, end_);
                }
                else if(input_drq == "RWDW"){
                    drqSampling drq(tk, input_output_path, enable_sparsebp, enable_preprocessor, relation_type, topk, time_limit, output_limit);
                    drq.RWDW(component_vec.front(), component_vec.back(), query_graph, data_graph, begin_, end_);
                }
                else{
                    std::cerr << "Invalid DRQ algorithm." << std::endl;
                    exit(1);
                }
            }
            else{//crq
                if(input_crq == "JointEnum"){
                    crqEnum crq(tk, input_output_path, enable_sparsebp, enable_preprocessor, relation_type, topk, time_limit, output_limit);
                    crq.execute_je(query_graph, data_graph, begin_, end_);
                }
                else if(input_crq == "LookAheadEnum"){
                    crqEnum crq(tk, input_output_path, enable_sparsebp, enable_preprocessor, relation_type, topk, time_limit, output_limit);
                    crq.execute_lae(query_graph, data_graph, begin_, end_);
                }
                else if(input_drq == "RWC"){
                    crqSampling crq(tk, input_output_path, enable_sparsebp, enable_preprocessor, relation_type, topk, time_limit, output_limit);
                    crq.RWC(query_graph, data_graph, begin_, end_);
                }
                else if(input_drq == "RWCH"){
                    crqSampling crq(tk, input_output_path, enable_sparsebp, enable_preprocessor, relation_type, topk, time_limit, output_limit);
                    crq.RWCH(query_graph, data_graph, begin_, end_);
                }
                else{
                    std::cerr << "Invalid CRQ algorithm." << std::endl;
                    exit(1);
                }
            }

            delete query_graph;
        }

        delete query_graph_ori;
    }

    delete data_graph;
    std::cout << "End." << std::endl;
    return 0;
}

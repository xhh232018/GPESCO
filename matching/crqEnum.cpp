#include "crqEnum.h"

void crqEnum::adjust_order(std::vector<uint32_t> &order, Graph* query_graph, uint32_t &begin_depth){
    if(begin_depth > 0) {
        uint32_t bn_cnt = 0;
        uint32_t v = order[begin_depth];
        uint32_t pos = 0;
        for(uint32_t i = 0; i < begin_depth;++i){
            uint32_t u = order[i];
            if(query_graph->checkEdgeExistence(u,v)){
                bn_cnt++;
                pos = i + 1;
            }
            if(bn_cnt >= 2) break;
        }

        if(bn_cnt > 1){
            uint32_t tmp = order[begin_depth];
            for(int i = query_graph->getVerticesCount() - 1; i > pos;i--){
                if(i <= begin_depth){
                    order[i] = order[i - 1];
                }
            }
            order[pos] = tmp;
            begin_depth = pos;
        }
    }
}

void crqEnum::execute_within_time_limit_lae(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute_lae(*catalog, output_limit);;
    });

    std::future_status status;
    do {
        status = future.wait_for(std::chrono::seconds(time_limit));
        if (status == std::future_status::deferred) {
            std::cout << "Deferred\n";
            exit(-1);
        } else if (status == std::future_status::timeout) {
            g_exit = true;
        }
    } while (status != std::future_status::ready);    
}

void crqEnum::execute_lae(Graph *query_graph, Graph *data_graph, uint32_t begin, uint32_t end){
    ofstream myout;
    myout.open(output_file,std::ofstream::app);

    // Execute preprocessor.
    auto pp = new preprocessor();
    auto storage = new catalog(query_graph, data_graph);
    storage->alg = "LookAheadEnum";

    pp->execute(query_graph, data_graph, storage, enable_elimination, begin, end);
    pp->print_metrics();

    delete pp;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Generate query plan..." << std::endl; 

    std::vector<std::vector<uint32_t>> spectrum;
    query_plan_generator::generate_query_plan_with_nd(query_graph, storage, spectrum);
    query_plan_generator::print_metrics(); 

    //begin_end depth
    for(int i = 0; i < spectrum.back().size();++i){
        if(spectrum.back()[i] == begin) storage->begin_depth = i;
        if(spectrum.back()[i] == end) storage->end_depth = i;
    }

    if(storage->begin_depth > storage->end_depth) std::swap(storage->begin_depth,storage->end_depth);

    //adjust the spectrum
    adjust_order(spectrum.back(),query_graph,storage->begin_depth);

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Encode..." << std::endl;
    auto en = new encoder(query_graph, data_graph->getVerticesCount());
    en->execute(storage, relation_type, enable_sparsebp, spectrum.back().data());
    en->print_metrics();
    delete en;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Enumerate..." << std::endl;

    uint64_t embedding_count = 0;
    double enumeration_time_in_ns = 0;

    auto tree = execution_tree_generator::generate_single_node_execution_tree(spectrum.back());

    execute_within_time_limit_lae(tree, storage, output_limit, time_limit);

    embedding_count = tree->get_output_count();
    enumeration_time_in_ns = tree->get_execution_time();

    printf("Enumerate time (seconds): %.6lf\n",NANOSECTOSEC(enumeration_time_in_ns));
    printf("#Embeddings: %zu\n", embedding_count);

    delete tree;

    /**
     * Counting
     */

    for(auto &kv_pair:storage->res_map){
        uint32_t small_node = kv_pair.first.first;
        uint32_t big_node = kv_pair.first.second;
        uint32_t tmp_cnt = kv_pair.second;

        if(rpq.size() < topk){
            rpq.push(result_struct(tmp_cnt,small_node,big_node));
        }
        else{
            if(tmp_cnt > rpq.top().inc_cnt){
                rpq.pop();
                rpq.push(result_struct(tmp_cnt,small_node,big_node));
            }
        }
    }

    while (!rpq.empty()){
        result_struct rs = rpq.top();

        uint32_t max_head = rs.head;
        uint32_t max_tail = rs.tail;
        uint64_t max = rs.inc_cnt;

        rpq.pop();

        myout << tk << ","  << spectrum.back()[storage->begin_depth] << "," << spectrum.back()[storage->end_depth] 
                << "," << max_head << "," << max_tail << "," << max << "," << NANOSECTOSEC(enumeration_time_in_ns) << endl;
    }

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;
    delete storage;

    myout.close();
}

void crqEnum::execute_within_time_limit_je(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute_je(*catalog, output_limit);;
    });

    std::future_status status;
    do {
        status = future.wait_for(std::chrono::seconds(time_limit));
        if (status == std::future_status::deferred) {
            std::cout << "Deferred\n";
            exit(-1);
        } else if (status == std::future_status::timeout) {
            g_exit = true;
        }
    } while (status != std::future_status::ready);    
}

void crqEnum::execute_je(Graph *query_graph, Graph *data_graph, uint32_t begin, uint32_t end) {
    ofstream myout;
    myout.open(output_file,std::ofstream::app);

    // Execute preprocessor.
    auto pp = new preprocessor();
    auto storage = new catalog(query_graph, data_graph);
    storage->alg = "JointEnum";

    pp->execute(query_graph, data_graph, storage, enable_elimination, begin, end);
    pp->print_metrics();

    delete pp;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Generate query plan..." << std::endl; 

    std::vector<std::vector<uint32_t>> spectrum;
    query_plan_generator::generate_query_plan_with_nd(query_graph, storage, spectrum);
    query_plan_generator::print_metrics(); 

    //begin_end depth
    for(int i = 0; i < spectrum.back().size();++i){
        if(spectrum.back()[i] == begin) storage->begin_depth = i;
        if(spectrum.back()[i] == end) storage->end_depth = i;
    }

    if(storage->begin_depth > storage->end_depth) std::swap(storage->begin_depth,storage->end_depth);

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Encode..." << std::endl;
    auto en = new encoder(query_graph, data_graph->getVerticesCount());
    en->execute(storage, relation_type, enable_sparsebp, spectrum.back().data());
    en->print_metrics();
    delete en;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Enumerate..." << std::endl;

    uint64_t embedding_count = 0;
    double enumeration_time_in_ns = 0;

    auto tree = execution_tree_generator::generate_single_node_execution_tree(spectrum.back());

    execute_within_time_limit_je(tree, storage, output_limit, time_limit);

    embedding_count = tree->get_output_count();
    enumeration_time_in_ns = tree->get_execution_time();

    printf("Enumerate time (seconds): %.6lf\n",NANOSECTOSEC(enumeration_time_in_ns));
    printf("#Embeddings: %zu\n", embedding_count);

    delete tree;

    /**
     * Counting
     */

    for(auto &kv_pair:storage->res_map){
        uint32_t small_node = kv_pair.first.first;
        uint32_t big_node = kv_pair.first.second;
        uint32_t tmp_cnt = kv_pair.second;

        if(rpq.size() < topk){
            rpq.push(result_struct(tmp_cnt,small_node,big_node));
        }
        else{
            if(tmp_cnt > rpq.top().inc_cnt){
                rpq.pop();
                rpq.push(result_struct(tmp_cnt,small_node,big_node));
            }
        }
    }

    while (!rpq.empty()){
        result_struct rs = rpq.top();

        uint32_t max_head = rs.head;
        uint32_t max_tail = rs.tail;
        uint64_t max = rs.inc_cnt;

        rpq.pop();

        myout << tk << ","  << spectrum.back()[storage->begin_depth] << "," << spectrum.back()[storage->end_depth] 
                << "," << max_head << "," << max_tail << "," << max << "," << NANOSECTOSEC(enumeration_time_in_ns) << endl;
    }

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;
    delete storage;
    
    myout.close();
}
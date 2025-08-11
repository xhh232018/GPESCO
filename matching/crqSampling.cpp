#include "crqSampling.h"

void crqSampling::execute_within_time_limit_rwc(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async,[tree, catalog, output_limit](){
        return tree->execute_rwc(*catalog, output_limit);
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

void crqSampling::execute_within_time_limit_rwch(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async,[tree, catalog, output_limit](){
        return tree->execute_rwch(*catalog, output_limit);
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

void crqSampling::exact_result(catalog* storage){ 
    for(auto &kv_pair:storage->res_map){
        uint32_t head = kv_pair.first.first;
        uint32_t tail = kv_pair.first.second;
        Edge max_e = {head,tail};
        uint64_t cnt = kv_pair.second / storage->total_sp_num;
        uint64_t suc_cnt = 0;
        auto iter = storage->success_map.find(max_e);
        if(iter == storage->success_map.end()){
            continue;
        }
        else{
            suc_cnt = iter->second;
        }
        if(rpq.size() < topk && suc_cnt >= 50){
            rpq.push(result_struct(cnt,head,tail));
        }
        else{
            if(cnt > rpq.top().inc_cnt && suc_cnt >= 50){
                rpq.pop();
                rpq.push(result_struct(cnt,head,tail));
            }
        }
    }
}

void crqSampling::exact_result_rwch(catalog* storage){

    bool last_step = (storage->max_depth + 1 == storage->query_graph_->getVerticesCount()) || (storage->max_depth == storage->query_graph_->getVerticesCount());

    if(last_step){
        for(auto &kv_pair:storage->res_map){
            uint32_t small_node = kv_pair.first.first;
            uint32_t big_node = kv_pair.first.second;
            uint64_t tmp_cnt = kv_pair.second;

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
    }
    else{
        for(auto &kv_pair:storage->res_map){
            uint32_t small_node = kv_pair.first.first;
            uint32_t big_node = kv_pair.first.second;
            uint64_t tmp_cnt = kv_pair.second;

            uint64_t success_cnt = 0;
            auto iter = storage->success_map.find(kv_pair.first);
            if(iter != storage->success_map.end()) success_cnt = iter->second;

            if(rpq.size() < topk && success_cnt > 30){
                rpq.push(result_struct(tmp_cnt,small_node,big_node));
            }
            else{
                if(tmp_cnt > rpq.top().inc_cnt && success_cnt > 30){
                    rpq.pop();
                    rpq.push(result_struct(tmp_cnt,small_node,big_node));
                }
            }
        }
    }
}

void crqSampling::RWC(Graph* query_graph, Graph* data_graph, uint32_t begin_, uint32_t end_){
    ofstream myout;
    myout.open(output_file,std::ofstream::app);

    // Execute preprocessor.
    auto pp = new preprocessor();
    auto storage = new catalog(query_graph, data_graph);
    storage->alg = "RWC";
    pp->execute(query_graph, data_graph, storage, enable_elimination, begin_, end_);
    pp->print_metrics();
    delete pp;


    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Generate query plan..." << std::endl;
    
    std::vector<std::vector<uint32_t>> spectrum;

    query_plan_generator::generate_query_plan_with_nd(query_graph, storage, spectrum);
    query_plan_generator::print_metrics();

    //begin_end depth
    for(int i = 0; i < spectrum.back().size();++i){
        if(spectrum.back()[i] == begin_) storage->begin_depth = i;
        if(spectrum.back()[i] == end_) storage->end_depth = i;
    }

    if(storage->begin_depth > storage->end_depth) std::swap(storage->begin_depth,storage->end_depth);

    std::cout << "--------------------------------------------------------------------" << std::endl;

    std::cout << "Encode..." << std::endl;
    auto en = new encoderSP(query_graph, data_graph->getVerticesCount());
    en->execute_with_shrinking(storage, relation_type, enable_sparsebp, spectrum.back().data(),3);
    en->print_metrics();
    delete en;

    std::cout << "--------------------------------------------------------------------" << std::endl;

    auto tree = execution_tree_generator::generate_single_node_execution_tree(spectrum.back());
    execute_within_time_limit_rwc(tree, storage, output_limit, time_limit);
    delete tree;

    exact_result(storage);
    //output enum
        
    while (!rpq.empty()){
        result_struct rs = rpq.top();

        uint32_t max_head = rs.head;
        uint32_t max_tail = rs.tail;

        rpq.pop();
        myout << tk << ","  << begin_ << "," << end_ << "," << max_head << "," << max_tail << "," << rs.inc_cnt << endl;

    }

    myout.close();
    delete storage;
}

void crqSampling::RWCH(Graph* query_graph, Graph* data_graph, uint32_t begin_, uint32_t end_){
    ofstream myout;
    myout.open(output_file,std::ofstream::app);

    // Execute preprocessor.
    auto pp = new preprocessor();
    auto storage = new catalog(query_graph, data_graph);
    storage->alg = "RWCH";
    pp->execute(query_graph, data_graph, storage, enable_elimination, begin_, end_);
    pp->print_metrics();
    delete pp;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Generate query plan..." << std::endl;
    
    std::vector<std::vector<uint32_t>> spectrum;

    query_plan_generator::generate_query_plan_with_nd(query_graph, storage, spectrum, begin_, end_);
    query_plan_generator::print_metrics();

    //begin_end depth
    for(int i = 0; i < spectrum.back().size();++i){
        if(spectrum.back()[i] == begin_) storage->begin_depth = i;
        if(spectrum.back()[i] == end_) storage->end_depth = i;
    }

    if(storage->begin_depth > storage->end_depth) std::swap(storage->begin_depth,storage->end_depth);

    std::cout << "--------------------------------------------------------------------" << std::endl;

    std::cout << "Encode..." << std::endl;
    auto en = new encoderSP(query_graph, data_graph->getVerticesCount());
    en->execute_with_shrinking(storage, relation_type, enable_sparsebp, spectrum.back().data(),3);
    en->print_metrics();
    delete en;

    storage->sampling_num = 10000;

    auto tree = execution_tree_generator::generate_single_node_execution_tree(spectrum.back());
    execute_within_time_limit_rwch(tree, storage, output_limit, time_limit);
    delete tree;

    exact_result_rwch(storage);
    //output enum
        
    while (!rpq.empty()){
        result_struct rs = rpq.top();

        uint32_t max_head = rs.head;
        uint32_t max_tail = rs.tail;

        rpq.pop();
        myout << tk << ","  << begin_ << "," << end_ << "," << max_head << "," << max_tail << "," << rs.inc_cnt << endl;

    }
    
    delete storage;
    myout.close();
}
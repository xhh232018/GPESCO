#include "drqSampling.h"

void drqSampling::adjust_order(Graph* large_comp, Graph* small_comp, std::vector<uint32_t> &big_order,  std::vector<uint32_t> &small_order,
                     uint32_t big_depth, uint32_t small_depth, std::vector<uint32_t> &new_order, uint32_t &total_core,
                      std::unordered_map<uint32_t,uint32_t> &rvs_map_big,
                      std::unordered_map<uint32_t,uint32_t> &rvs_map_small){
    
    std::vector<uint32_t> new_order_big, new_order_small;
    uint32_t core_big_aft, core_small_aft;
    adj_exe(large_comp,big_order,big_depth,new_order_big,core_big_aft);
    adj_exe(small_comp,small_order,small_depth,new_order_small,core_small_aft);

    //combine - [u,v, core,non-core]
    total_core = core_big_aft + core_small_aft;

    //init
    new_order.push_back(rvs_map_small[new_order_small.front()]);
    new_order.push_back(rvs_map_big[new_order_big.front()]);
    
    for(int i = 1; i < core_big_aft;++i){
        new_order.push_back(rvs_map_big[new_order_big[i]]);
    }

    for(int i = 1; i < core_small_aft; ++i){
        new_order.push_back(rvs_map_small[new_order_small[i]]);
    }

    for(int i = core_big_aft; i < new_order_big.size();++i){
        new_order.push_back(rvs_map_big[new_order_big[i]]);
    }

    for(int i = core_small_aft; i < new_order_small.size();++i){
        new_order.push_back(rvs_map_small[new_order_small[i]]);
    }
    
}


void drqSampling::count_single(Graph* data_graph, uint32_t u_label, catalog* storage){
    uint32_t can_cnt;
    const ui* candidate_endnode = data_graph->getVerticesByLabel(u_label,can_cnt);

    for(auto &kv_pair:storage->res_map_iso){
        uint32_t head = kv_pair.first;
        uint64_t cnt = kv_pair.second / storage->freq_map_iso[head];
        uint64_t start_cnt = kv_pair.second;

        uint32_t succ_cnt = 0;
        auto iter = storage->success_map_iso.find(head);
        if(iter != storage->success_map_iso.end()){
            succ_cnt = iter->second;
        }

        for(int i = 0; i < can_cnt; ++i){
            if(rpq.size()== topk && cnt <= rpq.top().inc_cnt) break;
            uint32_t tail = candidate_endnode[i];
            if(!data_graph->checkEdgeExistence(head,tail)){
                if(rpq.size() < topk && succ_cnt >= 50){
                    rpq.push(result_struct(cnt,head,tail));
                }
                else{
                    if(cnt > rpq.top().inc_cnt && succ_cnt >= 50){
                        rpq.pop();
                
                        //get end node
                        uint32_t tail;
                        for(int i = 0; i < can_cnt; ++i){
                            uint32_t end_candidate = candidate_endnode[i];
                            if(!data_graph->checkEdgeExistence(head,end_candidate)){

                                tail = end_candidate;
                                break;
                            }
                        }
                        rpq.push(result_struct(cnt,head,tail));
                    }
                }
            }
        }
    }

}

void drqSampling::getVEList(Graph* query_graph, std::vector<uint32_t> &component, 
                std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> &vertex_list,
                std::vector<std::pair<uint32_t, uint32_t>> &edge_list,
                std::unordered_map<uint32_t,uint32_t> &rvs_map){

    for(int i = 0; i < component.size(); ++i){
        uint32_t u = component[i];
        uint32_t label = query_graph->getVertexLabel(u);
        uint32_t deg = 0;
        for(int j = 0; j < component.size(); ++j){
            if(j == i) continue;
            if(query_graph->checkEdgeExistence(u,component[j])){
                deg++;
                if(i < j) edge_list.push_back({i,j});
            } 
        }
        vertex_list.push_back({i,label,deg});
        rvs_map[i] = u;
    }
}

void drqSampling::bfs(std::vector<std::vector<uint32_t>> &graph, uint32_t S, std::vector<uint32_t>& par, std::vector<uint32_t>& dist){
    std::queue<uint32_t> q;
    dist[S] = 0;
    q.push(S);

    while (!q.empty()) {
        // Pop the node at the front of the queue
        int node = q.front();
        q.pop();

        // Explore all the neighbours of the current node
        for (auto neighbour : graph[node]) {
            // Check if the neighbouring node is not visited
            if (dist[neighbour] == 1e9) {
                // Mark the current node as the parent of
                // the neighbouring node
                par[neighbour] = node;
                // Mark the distance of the neighbouring
                // node as distance of the current node + 1
                dist[neighbour] = dist[node] + 1;
                // Insert the neighbouring node to the queue
                q.push(neighbour);
            }
        }
    }
}

void drqSampling::printShortDist(std::vector<std::vector<uint32_t>> &forward_nbr_vec, uint32_t source, uint32_t dest, uint32_t graph_size,std::vector<uint32_t>& path){
    std::vector<uint32_t> par(graph_size, -1);

    std::vector<uint32_t> dist(graph_size, 1e9);

    bfs(forward_nbr_vec, source, par, dist);

    if (dist[dest] == 1e9) {
        cout << "Source and Destination are not connected";
        return;
    }

    uint32_t currentNode = dest;
    path.push_back(dest);
    while (par[currentNode] != -1) {
        path.push_back(par[currentNode]);
        currentNode = par[currentNode];
    }
}

void drqSampling::adj_exe(Graph* comp_graph,std::vector<uint32_t> &old_order,uint32_t ext_depth, std::vector<uint32_t> &new_order,uint32_t &core_num_aft){
    std::vector<std::vector<uint32_t>> forward_nbr_vec(old_order.size());

    for(int i = 0; i < old_order.size();++i){
        uint32_t u = old_order[i];
        for(int j = i + 1; j < old_order.size();++j){
            uint32_t v = old_order[j];
            if(comp_graph->checkEdgeExistence(u,v)){
                forward_nbr_vec[u].emplace_back(v);
            }
        }
    }

    uint32_t source = old_order[0];
    uint32_t destination = old_order[ext_depth];
    std::vector<uint32_t> path;

    printShortDist(forward_nbr_vec, source, destination, old_order.size(),path);

    //std::cout<< "Adjust Begins!"<< endl;
    uint32_t core_num_old = comp_graph->get2CoreSize();
    uint32_t core_num_current = core_num_old;

    for(int i = 0; i < path.size();++i){
        uint32_t u = path[i];
        auto idx = std::find(old_order.begin(),old_order.end(),u) - old_order.begin();
        if(idx >= core_num_old) core_num_current++;
        new_order.emplace_back(u);
    }
    for(auto &node:old_order){
        auto iter = std::find(new_order.begin(),new_order.end(),node);
        if(iter == new_order.end()) new_order.emplace_back(node);
    }
    if(core_num_old == 0) core_num_current = 1;

    core_num_aft = core_num_current;
}

void drqSampling::adjust_single(Graph* large_comp,std::vector<uint32_t> &big_order, uint32_t big_depth, std::vector<uint32_t> &new_order,uint32_t &total_core){

    std::vector<uint32_t> new_order_big;
    uint32_t core_big_aft;
    adj_exe(large_comp,big_order,big_depth,new_order_big,core_big_aft);

    total_core = core_big_aft;

    for(int k = 0; k < new_order_big.size();++k){
        new_order.emplace_back(new_order_big[k]);
    }
}

void drqSampling::common_label_detection(Graph* query_graph, catalog* storage,std::vector<uint32_t> &comp_small, std::vector<uint32_t> &comp_big){
    for(int k = 0; k < comp_small.size(); ++k){
        uint32_t label = query_graph->getVertexLabel(comp_small[k]);
        bool same = false;
        for(int q = 0; q < comp_big.size(); ++q){
            uint32_t label_ref = query_graph->getVertexLabel(comp_big[q]);
            if(label == label_ref){
                same = true;
                break;
            }
        }
        if(same){
            auto iter = std::find(storage->common_label_vec.begin(),storage->common_label_vec.end(),label);
            if(iter == storage->common_label_vec.end()) storage->common_label_vec.emplace_back(label);
        }
    }
}

void drqSampling::common_label_depth(Graph* component_graph, catalog* catalog, std::vector<uint32_t> &order){

    catalog->cl_depth_vec.resize(catalog->common_label_vec.size());
    for(int i = 0; i < catalog->common_label_vec.size();++i){
        uint32_t ref_label = catalog->common_label_vec[i];
        for(int k = 0; k < order.size();++k){
            uint32_t u = order[k];
            uint32_t label_u = component_graph->getVertexLabel(u);
            if(label_u == ref_label) catalog->cl_depth_vec[i].emplace_back(k);
        }
    }
}

void drqSampling::adjust_order(Graph* query_graph,std::vector<uint32_t> &order){
    uint32_t last_vertex = order.back();
    if(query_graph->getVertexDegree(last_vertex) > 0){
        uint32_t zero_position;
        for(int i = 0; i < order.size(); ++i){
            if(query_graph->getVertexDegree(order[i]) == 0){
                zero_position = i;
                break;
            }
        }
        uint32_t zero_vtx = order[zero_position];
        for(int k = zero_position; k < query_graph -> getVerticesCount() - 1; ++k){
            order[k] = order[k + 1];
        }
        order[query_graph -> getVerticesCount() - 1] = zero_vtx;
    }
}

void drqSampling::execute_within_time_limit_rwd(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute_rwd(*catalog, output_limit);
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

void drqSampling::execute_within_time_limit_rwdw_epl(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute_rwdw_epl(*catalog, output_limit);
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

void drqSampling::execute_within_time_limit_rwdw_final(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute_rwdw_final(*catalog, output_limit);
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


void drqSampling::RWD(Graph* query_graph_ori, Graph* query_graph, Graph* data_graph, uint32_t begin_, uint32_t end_){
    ofstream myout;
    myout.open(output_file,std::ofstream::app);
    uint32_t core_num = query_graph_ori->get2CoreSize();

    auto pp_ori = new preprocessor();
    auto storage_ori = new catalog(query_graph_ori, data_graph);

    pp_ori->execute(query_graph_ori, data_graph, storage_ori, true);
    pp_ori->print_metrics();
    delete pp_ori;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Generate original query plan..." << std::endl;

    std::vector<std::vector<uint32_t>> spectrum_ori;
    query_plan_generator::generate_query_plan_with_nd(query_graph_ori, storage_ori, spectrum_ori);
    query_plan_generator::print_metrics();   

    delete storage_ori;

    auto pp = new preprocessor();
    catalog* storage = new catalog(query_graph, data_graph);

    pp->execute(query_graph, data_graph, storage, enable_elimination, begin_, end_); 
    pp->print_metrics();
    delete pp;

    uint32_t deg_begin = query_graph -> getVertexDegree(begin_);
    uint32_t deg_end = query_graph -> getVertexDegree(end_);

    if(deg_begin == 0 || deg_end == 0) adjust_order(query_graph,spectrum_ori.back());

    //begin_end depth
    for(int i = 0; i < spectrum_ori.back().size();++i){
        if(spectrum_ori.back()[i] == begin_) storage->begin_depth = i;
        if(spectrum_ori.back()[i] == end_) storage->end_depth = i;
    }

    if(storage->begin_depth > storage->end_depth) std::swap(storage->begin_depth,storage->end_depth);
    storage->combine_core = core_num; 

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Encode..." << std::endl;
    auto en = new encoderSP(query_graph, data_graph->getVerticesCount());
    en->execute_nc(storage, relation_type, enable_sparsebp, spectrum_ori.back().data());
    en->print_metrics();
    delete en;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    auto tree = execution_tree_generator::generate_single_node_execution_tree(spectrum_ori.back());
    execute_within_time_limit_rwd(tree, storage, output_limit, time_limit);
    delete tree;

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
        if(rpq.size() < topk){
            rpq.push(result_struct(cnt,head,tail));
        }
        else{
            if(cnt > rpq.top().inc_cnt){
                rpq.pop();
                rpq.push(result_struct(cnt,head,tail));
            }
        }
    }

    while (!rpq.empty()){
        result_struct rs = rpq.top();

        uint32_t max_head = rs.head;
        uint32_t max_tail = rs.tail;

        rpq.pop();
        myout << tk << ","  << begin_ << "," << end_ << "," << max_head << "," << max_tail << "," << rs.inc_cnt << endl;

    }

    myout.close();
}

void drqSampling::RWDW(std::vector<uint32_t> &component_small,std::vector<uint32_t> &component_large,Graph* query_graph, Graph* data_graph, uint32_t begin_, uint32_t end_){
    ofstream myout;
    myout.open(output_file,std::ofstream::app);
    
    if(component_small.size() == 1){
        uint32_t total_limit = 200;
        uint32_t node = component_small.front();
        uint32_t iso_label = query_graph->getVertexLabel(node);
        Res_PQ can_pq;

        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list;
        std::unordered_map<uint32_t, uint32_t> rvs_map;

        uint64_t epl_limit = time_limit / 2;

        getVEList(query_graph,component_large,vertex_list, edge_list, rvs_map);
        auto component_graph = new Graph(true);
        component_graph->loadGraphFromMemory(vertex_list,edge_list);
        component_graph->buildCoreTable();

        auto storage = new catalog(component_graph, data_graph);
        storage->alg = "rwdw";
        common_label_detection(query_graph,storage,component_small,component_large);

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Preprocess..." << std::endl;
        // Execute preprocessor.
        auto pp = new preprocessor();
        pp->execute(component_graph, data_graph, storage, enable_elimination);
        pp->print_metrics();

        delete pp;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Generate query plan..." << std::endl;

        std::vector<std::vector<uint32_t>> spectrum;
        query_plan_generator::generate_query_plan_with_nd(component_graph, storage, spectrum);
        query_plan_generator::print_metrics();

        for(int i = 0; i < spectrum.back().size();++i){
            uint32_t u_trun = spectrum.back()[i];
            uint32_t u_ori = rvs_map[u_trun];
            if(u_ori == begin_ || u_ori == end_){
                storage->end_depth = i;
                break;
            }
        }

        common_label_depth(component_graph,storage,spectrum.back());

        std::vector<uint32_t> new_order_small;
        uint32_t total_core = 0;
        adjust_single(component_graph,spectrum.back(),storage->end_depth,new_order_small,total_core);
        storage->combine_core = total_core;
        storage->anc = true;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Encode..." << std::endl;

        auto en = new encoderSP(component_graph, data_graph->getVerticesCount());
        en->execute(storage, relation_type, enable_sparsebp, new_order_small.data());
        en->print_metrics();
        delete en;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        auto tree_sp = execution_tree_generator::generate_single_node_execution_tree(new_order_small);
        execute_within_time_limit_rwdw_epl(tree_sp,storage,output_limit,epl_limit);
        delete tree_sp;

        //create distribution
        auto storage_total = new catalog(component_graph, data_graph);
        
        for(auto &kv_pair: storage->res_node){
            uint64_t est = kv_pair.second / storage->total_sp_num;
            uint64_t success = storage->res_success[kv_pair.first];
            if(est > 0 && success > 50){
                if(can_pq.size() < total_limit){
                    can_pq.push(result_struct(est,kv_pair.first,0));
                }
                else{
                    if(est > can_pq.top().inc_cnt){
                        can_pq.pop();
                        can_pq.push(result_struct(est,kv_pair.first,0));
                    }
                }
            }
        }
    
        if(can_pq.empty()){
            std::cout << "No candidate nodes!" << std::endl;
            return;
        }

        while (!can_pq.empty()){
            result_struct rs = can_pq.top();
            can_pq.pop();
            storage_total->iso_vec.emplace_back(rs.head);
            storage_total->weight_vec_big.emplace_back(rs.inc_cnt);
        }

        storage_total->iso = true;

        std::vector<uint32_t> combine_order;
        total_core = 0;
        adjust_single(component_graph,spectrum.back(),storage->end_depth,combine_order,total_core);
        auto pp_total = new preprocessor();
        pp_total->execute(component_graph, data_graph, storage_total, false);
        storage_total->combine_core = total_core;
        delete pp_total;

        auto en_total = new encoderSP(component_graph, data_graph->getVerticesCount());
        en_total->execute_iso(storage_total, relation_type, enable_sparsebp, combine_order.data());
        en_total->print_metrics();
        delete en_total;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Sampling Combination ISO..." << std::endl;
        //Sampling Again!
        auto tree_comb = execution_tree_generator::generate_single_node_execution_tree(combine_order);
        execute_within_time_limit_rwdw_final(tree_comb,storage_total,output_limit,epl_limit);
        delete tree_comb;

        count_single(data_graph,iso_label,storage_total);

        while (!rpq.empty()){
                
            result_struct rs = rpq.top();

            uint32_t max_head = rs.head;
            uint32_t max_tail = rs.tail;
            uint64_t max = rs.inc_cnt;

            rpq.pop();

            myout << tk << ","  << begin_ << "," << end_<< "," << max_head << "," << max_tail << "," << max  << endl;

        }

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Release Memory" << std::endl;
        delete storage;
        delete component_graph;
        delete storage_total;
    }
    else{
        double k_s = 0, k_l = 0;
        k_s = (double)component_small.size() / query_graph -> getVerticesCount();
        k_l = (double)component_small.size() / query_graph -> getVerticesCount();
        uint64_t epl_time = time_limit / 2;

        //construct small component
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list_small;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list_small;
        std::unordered_map<uint32_t, uint32_t> rvs_map_small;

        getVEList(query_graph,component_small,vertex_list_small, edge_list_small, rvs_map_small);

        auto component_small_graph = new Graph(true);
        component_small_graph->loadGraphFromMemory(vertex_list_small,edge_list_small);
        component_small_graph->buildCoreTable();
        auto storage_small = new catalog(component_small_graph, data_graph);

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Preprocess..." << std::endl;
        // Execute preprocessor.
        auto pp_small = new preprocessor();
        pp_small->execute(component_small_graph, data_graph, storage_small, enable_elimination);
        pp_small->print_metrics();
        delete pp_small; 

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Generate query plan..." << std::endl;

        std::vector<std::vector<uint32_t>> spectrum_small;
        query_plan_generator::generate_query_plan_with_nd(component_small_graph, storage_small, spectrum_small);
        query_plan_generator::print_metrics();

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Encode..." << std::endl;

        for(int i = 0; i < spectrum_small.back().size();++i){
            uint32_t u_trun = spectrum_small.back()[i];
            uint32_t u_ori = rvs_map_small[u_trun];
            if(u_ori == begin_ || u_ori == end_){
                storage_small->end_depth = i;
                break;
            }
        }

        std::vector<uint32_t> new_order_small, new_order_big;
        uint32_t total_core = 0;
        adjust_single(component_small_graph,spectrum_small.back(),storage_small->end_depth,new_order_small,total_core);
        storage_small->combine_core = total_core;
        storage_small->anc = true;

        auto en_small = new encoderSP(component_small_graph, data_graph->getVerticesCount());
        en_small->execute(storage_small, relation_type, enable_sparsebp, new_order_small.data());
        en_small->print_metrics();
        delete en_small;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Sampling Small Part..." << std::endl;

        auto tree = execution_tree_generator::generate_single_node_execution_tree(new_order_small);
        execute_within_time_limit_rwdw_epl(tree, storage_small, output_limit, k_s * epl_time);
        delete tree;

        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list_big;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list_big;
        std::unordered_map<uint32_t, uint32_t> rvs_map_big;

        getVEList(query_graph,component_large,vertex_list_big, edge_list_big, rvs_map_big);

        auto component_big_graph = new Graph(true);
        component_big_graph->loadGraphFromMemory(vertex_list_big,edge_list_big);
        component_big_graph->buildCoreTable();

        auto storage_big = new catalog(component_big_graph, data_graph);

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Preprocess..." << std::endl;
        // Execute preprocessor.
        auto pp_big = new preprocessor();
        pp_big->execute(component_big_graph, data_graph, storage_big, enable_elimination);
        pp_big->print_metrics();
        delete pp_big;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Generate query plan..." << std::endl;

        std::vector<std::vector<uint32_t>> spectrum_big;
        query_plan_generator::generate_query_plan_with_nd(component_big_graph, storage_big, spectrum_big);
        query_plan_generator::print_metrics();

        for(int i = 0; i < spectrum_big.back().size();++i){
            uint32_t u_trun = spectrum_big.back()[i];
            uint32_t u_ori = rvs_map_big[u_trun];
            if(u_ori == begin_ || u_ori == end_){
                storage_big->end_depth = i;
                break;
            }
        }

        total_core = 0;
        adjust_single(component_big_graph,spectrum_big.back(),storage_big->end_depth,new_order_big,total_core);
        storage_big->combine_core = total_core;
        storage_big->anc = true;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Encode..." << std::endl;

        auto en_big = new encoderSP(component_big_graph, data_graph->getVerticesCount());
        en_big->execute(storage_big, relation_type, enable_sparsebp, new_order_big.data());
        en_big->print_metrics();
        delete en_big;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Sampling large Part..." << std::endl;

        auto tree_big = execution_tree_generator::generate_single_node_execution_tree(new_order_big);
        execute_within_time_limit_rwdw_epl(tree_big, storage_big, output_limit, k_l * epl_time);
        delete tree_big;

        auto storage_total = new catalog(query_graph, data_graph);
        uint32_t comp_limit = 200,total_limit = 1000;
        Res_PQ ptb_candidates_small, ptb_candidates_big;

        for(auto &kv_pair: storage_small->res_node){
            uint32_t small_node = kv_pair.first;
            uint64_t small_cnt = kv_pair.second / storage_small->total_sp_num;
            uint64_t success = storage_small->res_success[kv_pair.first];
            if(small_cnt > 0){
                if(ptb_candidates_small.size() < comp_limit && success > 50){
                    ptb_candidates_small.push(result_struct(small_cnt,small_node,0));
                }
                else{
                    if(small_cnt > ptb_candidates_small.top().inc_cnt && success > 50){
                        ptb_candidates_small.pop();
                        ptb_candidates_small.push(result_struct(small_cnt,small_node,0));
                    }
                }
            }
        } 

        for(auto &kv_pair: storage_big->res_node){
            uint32_t big_node = kv_pair.first;
            uint64_t big_cnt = kv_pair.second / storage_big->total_sp_num;
            uint64_t success = storage_big->res_success[kv_pair.first];
            if(big_cnt > 0){
                if(ptb_candidates_big.size() < comp_limit && success > 50){
                    ptb_candidates_big.push(result_struct(big_cnt,big_node,0));
                }
                else{
                    if(big_cnt > ptb_candidates_big.top().inc_cnt && success > 50){
                        ptb_candidates_big.pop();
                        ptb_candidates_big.push(result_struct(big_cnt,big_node,0));
                    }
                }
            }
        }

        if(!ptb_candidates_big.empty() && !ptb_candidates_small.empty()){
            Res_PQ can_pq;

            //candidate generation
            std::vector<uint32_t> candidate_small,candidate_big;
            std::vector<uint64_t> weight_small,weight_big;
            while (!ptb_candidates_small.empty()){
                result_struct rs = ptb_candidates_small.top();
                ptb_candidates_small.pop();
                candidate_small.push_back(rs.head);
                weight_small.push_back(rs.inc_cnt);
            }

            while (!ptb_candidates_big.empty()){
                result_struct rs = ptb_candidates_big.top();
                ptb_candidates_big.pop();
                candidate_big.push_back(rs.head);
                weight_big.push_back(rs.inc_cnt);
            }

            for(int j = 0;j < candidate_small.size();++j){
                uint32_t small_node = candidate_small[j];
                uint64_t small_cnt = weight_small[j];
                for(int k = 0; k < candidate_big.size();++k){
                    uint32_t big_node = candidate_big[k];
                    if(data_graph->checkEdgeExistence(small_node,big_node)) continue;
                    uint64_t big_cnt = weight_big[k];
                    uint64_t cnt = big_cnt * small_cnt;
                    if(can_pq.size() < total_limit){
                        can_pq.push(result_struct(cnt,small_node,big_node));
                    }
                    else{
                        if(cnt > can_pq.top().inc_cnt){
                            can_pq.pop();
                            can_pq.push(result_struct(cnt,small_node,big_node));
                        }
                    }
                }
            }

            while(!can_pq.empty()){
                result_struct rs = can_pq.top();
                uint32_t small_node = rs.head;
                uint32_t big_node = rs.tail;
                uint64_t cnt = rs.inc_cnt;
                can_pq.pop();
                storage_total->node_small.emplace_back(small_node);
                storage_total->node_big.emplace_back(big_node);
                storage_total->weight_vec_big.emplace_back(cnt);
            }

            std::vector<uint32_t> combine_order;
            uint32_t t_core = 0;
            adjust_order(component_big_graph,component_small_graph,spectrum_big.back(),spectrum_small.back(),
                            storage_big->end_depth,storage_small->end_depth,combine_order, t_core, rvs_map_big,rvs_map_small);

            auto pp_total = new preprocessor();
            pp_total->execute(query_graph, data_graph, storage_total, false, begin_, end_);
            storage_total->combine_core = t_core;
            delete pp_total;

            auto en_total = new encoderSP(query_graph, data_graph->getVerticesCount());
            en_total->execute_iso(storage_total, relation_type, enable_sparsebp, combine_order.data());
            en_total->print_metrics();
            delete en_total;

            std::cout << "--------------------------------------------------------------------" << std::endl;
            std::cout << "Sampling Combine..." << std::endl;
            //Sampling Over Distribution
            auto tree_comb = execution_tree_generator::generate_single_node_execution_tree(combine_order);
            execute_within_time_limit_rwdw_final(tree_comb,storage_total,output_limit,epl_time);
            delete tree_comb;

            /**
             * Counting
             */
            std::cout << "--------------------------------------------------------------------" << std::endl;
            std::cout << "Counting..." << std::endl;

            for(auto &kv_pair:storage_total->res_map){
                uint32_t head = kv_pair.first.first;
                uint32_t tail = kv_pair.first.second;
                Edge max_e = {head,tail};
                uint64_t cnt = kv_pair.second / storage_total->freq_map[max_e];
                uint64_t suc_cnt = 0;
                auto iter = storage_total->success_map.find(max_e);
                if(iter == storage_total->success_map.end()){
                    continue;
                }
                else{
                    suc_cnt = iter->second;
                }
                if(rpq.size() < topk && suc_cnt > 100){
                    rpq.push(result_struct(cnt,head,tail));
                }
                else{
                    if(cnt > rpq.top().inc_cnt && suc_cnt > 100){
                        rpq.pop();
                        rpq.push(result_struct(cnt,head,tail));
                    }
                }
            }                                                
        }

        delete storage_small;
        delete component_small_graph;
        delete storage_big;
        delete component_big_graph;
        delete query_graph;
        delete storage_total;

        while (!rpq.empty()){
            result_struct rs = rpq.top();

            uint32_t max_head = rs.head;
            uint32_t max_tail = rs.tail;

            rpq.pop();
            myout << tk << ","  << begin_ << "," << end_ << "," << max_head << "," << max_tail << "," << rs.inc_cnt << endl;

        }

    }

    myout.close();
}
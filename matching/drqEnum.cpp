#include "drqEnum.h"

void drqEnum::candidates_pair(catalog* storage,std::vector<uint32_t> &new_order){
    uint32_t u_small = new_order.front();
    uint32_t u_big = new_order[1];
    std::vector<uint32_t> small_can,big_can;
    for(auto &pair:storage->candidate_pairs){
        small_can.emplace_back(pair.first);
        big_can.emplace_back(pair.second);
    }
    //delete replcates
    sort(small_can.begin(),small_can.end());
    small_can.erase(unique(small_can.begin(), small_can.end() ), small_can.end());
    sort(big_can.begin(),big_can.end());
    big_can.erase(unique(big_can.begin(), big_can.end() ), big_can.end());

    storage->candidate_sets_[u_small] = new uint32_t[small_can.size()];
    memcpy(storage->candidate_sets_[u_small], small_can.data(), sizeof(uint32_t) * small_can.size());
    storage->num_candidates_[u_small] = small_can.size();
    
    storage->candidate_sets_[u_big] = new uint32_t[big_can.size()];
    memcpy(storage->candidate_sets_[u_big], big_can.data(), sizeof(uint32_t) * big_can.size());
    storage->num_candidates_[u_big] = big_can.size();
}

void drqEnum::bfs(std::vector<std::vector<uint32_t>> &graph, uint32_t S, std::vector<uint32_t>& par, std::vector<uint32_t>& dist){
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

void drqEnum::printShortDist(std::vector<std::vector<uint32_t>> &forward_nbr_vec, uint32_t source, uint32_t dest, uint32_t graph_size,std::vector<uint32_t>& path){
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

    // // printing path from source to destination
    // for (int i = path.size() - 1; i >= 0; i--)
    //     cout << path[i] << " ";
}

void drqEnum::adj_exe(Graph* comp_graph,std::vector<uint32_t> &old_order,uint32_t ext_depth, std::vector<uint32_t> &new_order,uint32_t &core_num_aft){
    
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

void drqEnum::adjust_order(Graph* large_comp, Graph* small_comp, std::vector<uint32_t> &big_order,  std::vector<uint32_t> &small_order,
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

void drqEnum::bound_pruning(Graph* data_graph,std::vector<std::pair<uint32_t,uint32_t>> &candidate_space,
                            catalog* storage_small, catalog* storage_big){
    uint64_t lower_global = 0;
    std::vector<uint64_t> upper;
    std::vector<std::pair<uint32_t,uint32_t>> candidates;

    std::priority_queue<uint64_t, std::vector<uint64_t>, std::greater<uint64_t>> lower_bound_pq;

    for(auto &kv_pair_small:storage_small -> res_map_single){
        uint64_t small_count = kv_pair_small.second;
        uint32_t small_node = kv_pair_small.first;

        for(auto &kv_pair_big:storage_big -> res_map_single){
            uint64_t big_count = kv_pair_big.second;
            uint64_t upper_current = small_count * big_count;

            if(lower_bound_pq.size() < topK){
                lower_global = 0;
            }
            else{
                lower_global = lower_bound_pq.top();
            }

            if(upper_current <= lower_global) continue;
            uint32_t big_node = kv_pair_big.first;
            if(data_graph->checkEdgeExistence(small_node,big_node) || small_node == big_node) continue;
            candidates.push_back({small_node,big_node});
            upper.emplace_back(upper_current);

            uint64_t lower_current = upper_current;
            uint64_t same_cnt = 0;
            //exclusion get lower bound of current combination
            for(int i = 0; i < storage_small->common_label_vec.size();++i){
                for(auto &kv_pair_cl_small:storage_small->common_label_map[small_node][i]){
                    uint32_t cl_node_small = kv_pair_cl_small.first;
                    uint32_t cl_node_small_cnt = kv_pair_cl_small.second;

                    auto iter = storage_big->common_label_map[big_node][i].find(cl_node_small);
                    if(iter != storage_big->common_label_map[big_node][i].end()) same_cnt += cl_node_small_cnt * iter->second;
                }
            }
            if(same_cnt >= lower_current) lower_current = 0;
            else lower_current -= same_cnt;

            if (lower_bound_pq.size() < topK) {
                lower_bound_pq.push(lower_current);
            } else {
                if (lower_current > lower_bound_pq.top()) {
                    lower_bound_pq.pop();
                    lower_bound_pq.push(lower_current);
                }
            }
        }

    }

    lower_global = lower_bound_pq.top();

    for(int i = 0; i < candidates.size();++i){
        if (upper[i] >= lower_global) {
            candidate_space.emplace_back(candidates[i]);
        }
    }

    while (!lower_bound_pq.empty())
    {
        lower_bound_pq.pop();
    }

}

void drqEnum::count_product(Graph* data_graph, catalog* storage_small,catalog* storage_big){
    for(auto &kv_pair_small:storage_small->res_map_single){
        uint32_t small_node = kv_pair_small.first;
        uint32_t small_cnt = kv_pair_small.second;
        for(auto &kv_pair_big:storage_big->res_map_single){
            uint32_t big_node = kv_pair_big.first;
            uint64_t tmp_cnt = small_cnt * kv_pair_big.second;
            
            if(!data_graph->checkEdgeExistence(small_node,big_node)){
                if(rpq.size() < topK){
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
    }
}

void drqEnum::count_single(Graph* data_graph, uint32_t u_label, catalog* storage){
    uint32_t can_cnt;
    const ui* candidate_endnode = data_graph->getVerticesByLabel(u_label,can_cnt);
    for(auto &kv_pair:storage->res_map_single){
        uint32_t head = kv_pair.first;
        uint64_t start_cnt = kv_pair.second;

        for(int i = 0; i < can_cnt; ++i){
            if(rpq.size()== topk && start_cnt <= rpq.top().inc_cnt) break;
            uint32_t end_candidate = candidate_endnode[i];
            if(!data_graph->checkEdgeExistence(head,end_candidate)){
                if(!storage->common_label_vec.empty()){
                    auto iter = storage->common_label_map[head].front().find(end_candidate);
                    uint32_t final_cnt = start_cnt;
                    if(iter == storage->common_label_map[head].front().end()) final_cnt = start_cnt;
                    else final_cnt = start_cnt - iter->second;

                    if(rpq.size() < topk){
                        rpq.push(result_struct(final_cnt,head,end_candidate));
                    }
                    else if(final_cnt > rpq.top().inc_cnt){
                        rpq.pop();
                        rpq.push(result_struct(final_cnt,head,end_candidate));
                    }
                }
                else{
                    if(rpq.size() < topk){
                        rpq.push(result_struct(kv_pair.second,head,end_candidate));
                    }
                    else if(kv_pair.second > rpq.top().inc_cnt){
                        rpq.pop();
                        rpq.push(result_struct(kv_pair.second,head,end_candidate));
                    }
                }
            }
        }
    }
}

void drqEnum::execute_within_time_limit(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute(*catalog, output_limit);;
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

void drqEnum::getVEList(Graph* query_graph, std::vector<uint32_t> &component, 
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

void drqEnum::common_label_depth(Graph* component_graph, catalog* catalog, std::vector<uint32_t> &order){

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

void drqEnum::cld_pe(Graph* component_graph, catalog* catalog, std::vector<uint32_t> &order){
    for(int i = 0; i < catalog->common_label_vec.size();++i){
        uint32_t ref_label = catalog->common_label_vec[i];
        if(i == catalog->end_depth) continue;
        for(int k = 0; k < order.size();++k){
            uint32_t u = order[k];
            uint32_t label_u = component_graph->getVertexLabel(u);
            if(label_u == ref_label) catalog->cld_pe.emplace_back(k);
        }
    }
    catalog->cld_pe.emplace_back(catalog->end_depth);
    catalog->probe_depth = catalog->cld_pe.size();
}

void drqEnum::common_label_detection(Graph* query_graph, catalog* storage,std::vector<uint32_t> &comp_small, std::vector<uint32_t> &comp_big){
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

void drqEnum::execute_within_time_limit_iee(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute_iee(*catalog, output_limit);;
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

void drqEnum::execute_within_time_limit_pe(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute_pe(*catalog, output_limit);;
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

void drqEnum::execute_within_time_limit_probe_mat(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit){
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute_probe_mat(*catalog, output_limit);;
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

void drqEnum::execute_iee(std::vector<uint32_t> &component_small,std::vector<uint32_t> &component_large,Graph* query_graph, Graph* data_graph, uint32_t begin, uint32_t end){
    
    ofstream myout;
    myout.open(output_file,std::ofstream::app);

    if(component_small.size() == 1){

        uint32_t node = component_small.front();
        uint32_t iso_label = query_graph->getVertexLabel(node);

        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list;
        std::unordered_map<uint32_t, uint32_t> rvs_map;

        getVEList(query_graph,component_large,vertex_list, edge_list, rvs_map);

        auto component_graph = new Graph(true);
        component_graph->loadGraphFromMemory(vertex_list,edge_list);
        component_graph->buildCoreTable();

        auto storage = new catalog(component_graph, data_graph);
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

        common_label_depth(component_graph,storage,spectrum.back());


        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Encode..." << std::endl;

        auto en = new encoder(component_graph, data_graph->getVerticesCount());
        en->execute(storage, relation_type, enable_sparsebp, spectrum.back().data());
        en->print_metrics();
        delete en;

        for(int i = 0; i < spectrum.back().size();++i){
            uint32_t u_trun = spectrum.back()[i];
            
            uint32_t u_ori = rvs_map[u_trun];

            if(u_ori == begin || u_ori == end){
                storage->end_depth = i;
                break;
            }
        }

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Enumerate..." << std::endl;

        auto tree = execution_tree_generator::generate_single_node_execution_tree(spectrum.back());

        execute_within_time_limit(tree, storage, output_limit, time_limit);

        double enumeration_time_in_ns = tree->get_execution_time();

        delete tree;

        count_single(data_graph,iso_label,storage);

        delete component_graph;
        delete storage;

        while (!rpq.empty()){
            result_struct rs = rpq.top();

            uint32_t max_head = rs.head;
            uint32_t max_tail = rs.tail;
            uint64_t max = rs.inc_cnt;

            rpq.pop();

            myout << tk << ","  << begin<< "," << end << "," << max_head << "," << max_tail << "," << max << "," <<NANOSECTOSEC(enumeration_time_in_ns) << endl;
        }
    }
    else{
        //construct small component
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list_small;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list_small;
        std::unordered_map<uint32_t, uint32_t> rvs_map_small;

        getVEList(query_graph,component_small,vertex_list_small, edge_list_small, rvs_map_small);

        auto component_small_graph = new Graph(true);
        component_small_graph->loadGraphFromMemory(vertex_list_small,edge_list_small);
        component_small_graph->buildCoreTable();

        auto storage_small = new catalog(component_small_graph, data_graph);
        storage_small->alg = "drq";
        common_label_detection(query_graph,storage_small,component_small,component_large);

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

        common_label_depth(component_small_graph,storage_small,spectrum_small.back());

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Encode..." << std::endl;

        auto en_small = new encoder(component_small_graph, data_graph->getVerticesCount());
        en_small->execute(storage_small, relation_type, enable_sparsebp, spectrum_small.back().data());
        en_small->print_metrics();
        delete en_small;

        for(int i = 0; i < spectrum_small.back().size();++i){
            uint32_t u_trun = spectrum_small.back()[i];
            uint32_t u_ori = rvs_map_small[u_trun];
            if(u_ori == begin || u_ori == end){
                storage_small->end_depth = i;
                break;
            }
        }

        auto tree_small = execution_tree_generator::generate_single_node_execution_tree(spectrum_small.back());

        execute_within_time_limit(tree_small, storage_small, output_limit, time_limit);

        double enumeration_time_in_ns_small = tree_small->get_execution_time();

        delete tree_small;

        time_limit -= NANOSECTOSEC(enumeration_time_in_ns_small);

        //construct big component enumerate with probing
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list_big;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list_big;
        std::unordered_map<uint32_t, uint32_t> rvs_map_big;

        getVEList(query_graph,component_large,vertex_list_big, edge_list_big, rvs_map_big);

        auto component_big_graph = new Graph(true);
        component_big_graph->loadGraphFromMemory(vertex_list_big,edge_list_big);
        component_big_graph->buildCoreTable();

        auto storage_big = new catalog(component_big_graph, data_graph);
        storage_big->alg = "drq";
        common_label_detection(query_graph,storage_big,component_small,component_large);

        std::cout << "--------------------------------------------------------------------" << std::endl;

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

        common_label_depth(component_big_graph,storage_big,spectrum_big.back());

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Encode..." << std::endl;

        auto en_big = new encoder(component_big_graph, data_graph->getVerticesCount());
        en_big->execute(storage_big, relation_type, enable_sparsebp, spectrum_big.back().data());
        en_big->print_metrics();        
        delete en_big;

        storage_big->common_label_map = storage_small->common_label_map;
        storage_big->res_map_single = storage_small->res_map_single;

        for(int i = 0; i < spectrum_big.back().size();++i){
            uint32_t u_trun = spectrum_big.back()[i];
            uint32_t u_ori = rvs_map_big[u_trun];
            if(u_ori == begin || u_ori == end){
                storage_big->end_depth = i;
                break;
            }
        }

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Enumerate..." << std::endl;

        auto tree_big = execution_tree_generator::generate_single_node_execution_tree(spectrum_big.back());

        execute_within_time_limit(tree_big, storage_big, output_limit, time_limit);

        double enumeration_time_in_ns_big = tree_big->get_execution_time();

        delete tree_big;

        time_limit -= NANOSECTOSEC(enumeration_time_in_ns_big);

        uint64_t enum_time = NANOSECTOSEC(enumeration_time_in_ns_big)+ NANOSECTOSEC(enumeration_time_in_ns_small);

        if(storage_big->common_label_vec.empty()) {//Cartesian Product
            count_product(data_graph,storage_small,storage_big);
            //print topk
            while (!rpq.empty()){
                result_struct rs = rpq.top();

                uint32_t max_head = rs.head;
                uint32_t max_tail = rs.tail;
                uint64_t max = rs.inc_cnt;

                rpq.pop();

                myout << tk << ","  << begin << "," << end << "," << max_head << "," << max_tail << "," << max << "," << enum_time << endl;
            }

        }
        else{//ie bound pruning
            auto storage_total = new catalog(query_graph, data_graph);
            storage_total->alg = "drq";
            auto start_bound = std::chrono::high_resolution_clock::now();
            bound_pruning(data_graph,storage_total->candidate_pairs,storage_small,storage_big);
            auto end_bound = std::chrono::high_resolution_clock::now();
            double bound_time = NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end_bound - start_bound).count());

            time_limit -= bound_time;
            enum_time += bound_time;

            std::vector<uint32_t> combine_order;
            uint32_t total_core;
            adjust_order(component_big_graph,component_small_graph,spectrum_big.back(),spectrum_small.back(),
                            storage_big->end_depth,storage_small->end_depth,combine_order, total_core, rvs_map_big,rvs_map_small);
            
            storage_total->combine_core = total_core;

            auto pp_total = new preprocessor();
            pp_total->execute(query_graph, data_graph, storage_total, false);
            delete pp_total;

            candidates_pair(storage_total,combine_order);

            auto en_iee = new encoder(query_graph, data_graph->getVerticesCount());
            en_iee->execute_iee(storage_total, relation_type, enable_sparsebp, combine_order.data());
            en_iee->print_metrics();
            delete en_iee;

            auto tree_iee = execution_tree_generator::generate_single_node_execution_tree(combine_order);
            execute_within_time_limit_iee(tree_iee,storage_total,output_limit,time_limit);
            double iee_time = NANOSECTOSEC(tree_iee->get_execution_time());
            delete tree_iee;

            enum_time += iee_time;
            //counting topk
            for(auto &kv_pair:storage_total->res_map){
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

            //print topk
            while (!rpq.empty()){
                result_struct rs = rpq.top();

                uint32_t max_head = rs.head;
                uint32_t max_tail = rs.tail;
                uint64_t max = rs.inc_cnt;

                rpq.pop();

                myout << tk << ","  << begin << "," << end<< "," << max_head << "," << max_tail << "," << max << ","<< enum_time << endl;
            }
            delete storage_total;       
        }

        delete component_small_graph;
        delete component_big_graph;
        delete storage_big;
        delete storage_small;
    }

    myout.close();
}

void drqEnum::execute_pe(std::vector<uint32_t> &component_small,std::vector<uint32_t> &component_large,Graph* query_graph, Graph* data_graph, uint32_t begin, uint32_t end){
    
    ofstream myout;
    myout.open(output_file,std::ofstream::app);

    //small = 1
    if(component_small.size() == 1){

        uint32_t node = component_small.front();
        uint32_t iso_label = query_graph->getVertexLabel(node);

        uint32_t can_cnt;
        const ui* candidate_endnode = data_graph->getVerticesByLabel(iso_label,can_cnt);
        
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list;
        std::unordered_map<uint32_t, uint32_t> rvs_map;

        getVEList(query_graph,component_large,vertex_list, edge_list, rvs_map);

        auto component_graph = new Graph(true);
        component_graph->loadGraphFromMemory(vertex_list,edge_list);
        component_graph->buildCoreTable();

        auto storage = new catalog(component_graph, data_graph);
        storage->alg = "drq";

        for(uint32_t k = 0; k < can_cnt;++k){
            storage->probe_instances.emplace_back(candidate_endnode[k]);
        }
        storage->probe_depth = 1;
        storage->probe_cnt = can_cnt;

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

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Encode..." << std::endl;

        auto en = new encoder(component_graph, data_graph->getVerticesCount());
        en->execute(storage, relation_type, enable_sparsebp, spectrum.back().data());
        en->print_metrics();
        delete en;

        for(int i = 0; i < spectrum.back().size();++i){
            uint32_t u_trun = spectrum.back()[i];
            
            uint32_t u_ori = rvs_map[u_trun];

            if(u_ori == begin || u_ori == end){
                storage->end_depth = i;
                break;
            }
        }

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Enumerate..." << std::endl;

        auto tree = execution_tree_generator::generate_single_node_execution_tree(spectrum.back());

        execute_within_time_limit_pe(tree, storage, output_limit, time_limit);
        double enum_time = NANOSECTOSEC(tree->get_execution_time());
        delete tree;

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
        //print topk
        while (!rpq.empty()){
            result_struct rs = rpq.top();

            uint32_t max_head = rs.head;
            uint32_t max_tail = rs.tail;
            uint64_t max = rs.inc_cnt;

            rpq.pop();

            myout << tk << ","  << begin<< "," << end << "," << max_head << "," << max_tail << "," << max << "," << enum_time << endl;
        
        }
        delete storage;
        delete component_graph;

    }
    else{
        //construct small component
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list_small;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list_small;
        std::unordered_map<uint32_t, uint32_t> rvs_map_small;

        getVEList(query_graph,component_small,vertex_list_small, edge_list_small, rvs_map_small);

        auto component_small_graph = new Graph(true);
        component_small_graph->loadGraphFromMemory(vertex_list_small,edge_list_small);
        component_small_graph->buildCoreTable();

        auto storage_small = new catalog(component_small_graph, data_graph);
        storage_small->alg = "drq";
        cld_pe(component_small_graph,storage_small,component_small);

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

        auto en_small = new encoder(component_small_graph, data_graph->getVerticesCount());
        en_small->execute(storage_small, relation_type, enable_sparsebp, spectrum_small.back().data());
        en_small->print_metrics();
        delete en_small;

        auto tree_mat = execution_tree_generator::generate_single_node_execution_tree(spectrum_small.back());

        execute_within_time_limit_probe_mat(tree_mat, storage_small, output_limit, time_limit);

        double enumeration_time_in_ns_small = tree_mat->get_execution_time();

        delete tree_mat;

        time_limit -= NANOSECTOSEC(enumeration_time_in_ns_small);

        //construct big component enumerate with probing
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> vertex_list_big;
        std::vector<std::pair<uint32_t, uint32_t>> edge_list_big;
        std::unordered_map<uint32_t, uint32_t> rvs_map_big;

        getVEList(query_graph,component_large,vertex_list_big, edge_list_big, rvs_map_big);

        auto component_big_graph = new Graph(true);
        component_big_graph->loadGraphFromMemory(vertex_list_big,edge_list_big);
        component_big_graph->buildCoreTable();

        auto storage_big = new catalog(component_big_graph, data_graph);
        storage_big->alg = "drq";

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
            if(u_ori == begin || u_ori == end){
                storage_big->end_depth = i;
                break;
            }
        }

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Encode..." << std::endl;

        auto en_big = new encoder(component_big_graph, data_graph->getVerticesCount());
        en_big->execute(storage_big, relation_type, enable_sparsebp, spectrum_big.back().data());
        en_big->print_metrics();        
        delete en_big;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Probe..." << std::endl;
        storage_big->probe_depth = storage_small->probe_depth;
        storage_big->probe_instances.swap(storage_small->probe_instances);
        storage_big->probe_cnt = storage_small->probe_cnt;

        auto tree_pe = execution_tree_generator::generate_single_node_execution_tree(spectrum_big.back());

        execute_within_time_limit_pe(tree_pe, storage_big, output_limit, time_limit);

        double enumeration_time_in_ns_big = tree_pe->get_execution_time();

        delete tree_pe;

        double enum_time = NANOSECTOSEC(enumeration_time_in_ns_big) + NANOSECTOSEC(enumeration_time_in_ns_small);

        for(auto &kv_pair:storage_big->res_map){
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
        //print topk
        while (!rpq.empty()){
            result_struct rs = rpq.top();

            uint32_t max_head = rs.head;
            uint32_t max_tail = rs.tail;
            uint64_t max = rs.inc_cnt;

            rpq.pop();
            
            myout << tk << ","  << begin<< "," << end << "," << max_head << "," << max_tail << "," << max << "," << enum_time << endl;
        
        }

        delete storage_big;
        delete storage_small;
        delete component_big_graph;
        delete component_small_graph;
    }

    myout.close();
}
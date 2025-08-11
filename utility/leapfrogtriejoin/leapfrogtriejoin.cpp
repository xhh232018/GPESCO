#include "global_variables.h"
#include "leapfrogtriejoin.h"
#include "../computesetintersection.h"
#include <random>
#include <queue>

bool leapfrogtriejoin::rwcw_non_core_final(uint32_t start_depth, uint32_t max_depth){
    uint32_t ext_depth = start_depth;

    std::random_device rd;
    std:mt19937 gen_nc(rd());

    std::uniform_int_distribution<> d(0,num_local_candidates_[start_depth] - 1);
    int random = d(gen_nc);

    uint32_t v = si_buffer_[random];
    if(visited_[v]){
        return false;
    }
    visited_[v] = true;
    embedding_[start_depth] = v;
    ext_depth++;

    prob = prob / (double)num_local_candidates_[start_depth];

    uint32_t cur_depth = start_depth + 1;
    uint32_t u = vertex_ordering_[cur_depth];
    uint32_t bn = get_backward_neighbors(cur_depth)[0];
    uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
    uint32_t key = embedding_[bns_depth];
    local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);

    bool valid_nc = true;

    if (cur_depth == max_depth - 1) {
        uint32_t temp_count = num_local_candidates_[cur_depth];
        uint32_t* temp_buffer = local_candidates_[cur_depth];
        if(temp_count > 0) {

            prob = prob / (double)temp_count;
            
            std::uniform_int_distribution<> d(0,temp_count - 1);
            int random = d(gen_nc);

            uint32_t temp_v = temp_buffer[random];
            if(!visited_[temp_v]){
                embedding_[cur_depth] = temp_v;
                if(catalog_->iso){
                    auto iter = catalog_->res_map_iso.find(embedding_[0]);
                    if(iter == catalog_->res_map_iso.end()) catalog_->res_map_iso[embedding_[0]] = 1/prob;
                    else iter->second += 1/prob;

                    auto iter_succ = catalog_->success_map_iso.find(embedding_[0]);
                    if(iter_succ == catalog_->success_map_iso.end()) catalog_->success_map_iso[embedding_[0]] = 1;
                    else iter_succ->second += 1;
                }
                else{
                    Edge ptb = {embedding_[0],embedding_[1]};
                    auto iter = catalog_->res_map.find(ptb);
                    if(iter == catalog_->res_map.end()) catalog_->res_map[ptb] = 1/prob;
                    else iter->second += 1/prob;

                    //#success
                    auto iter_succ = catalog_->success_map.find(ptb);
                    if(iter_succ == catalog_->success_map.end()) catalog_->success_map[ptb] = 1;
                    else iter_succ->second += 1;
                }
            }
        }
    }
    else{
        while (valid_nc)
        {
           while (cur_depth <= max_depth - 2){
                if(num_local_candidates_[cur_depth] == 0) {
                    valid_nc = false;
                    break;
                }
                prob = prob / (double)num_local_candidates_[cur_depth];

                std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                int random = d(gen_nc);

                uint32_t v = local_candidates_[cur_depth][random];
                
                if (visited_[v]) {
                    valid_nc = false;
                    break;
                }

                visited_[v] = true;
                embedding_[cur_depth] = v;
                ext_depth++;

                uint32_t next_depth = cur_depth + 1;
                uint32_t u = vertex_ordering_[next_depth];
                uint32_t bn = get_backward_neighbors(next_depth)[0];
                uint32_t bns_depth = get_backward_neighbors_depth(next_depth)[0]; 
                uint32_t key = embedding_[bns_depth];
                local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);

                if(cur_depth == max_depth - 2){
                    uint32_t temp_count = num_local_candidates_[next_depth];
                    uint32_t *temp_buffer = local_candidates_[next_depth];    
                    if(temp_count == 0) {
                        valid_nc = false;
                        break;
                    }
                    prob = prob / (double)temp_count;

                    std::uniform_int_distribution<> d(0,temp_count - 1);
                    int random = d(gen_nc);

                    uint32_t temp_v = temp_buffer[random];
                    if(!visited_[temp_v]){
                        if(catalog_->iso){
                            auto iter = catalog_->res_map_iso.find(embedding_[0]);
                            if(iter == catalog_->res_map_iso.end()) catalog_->res_map_iso[embedding_[0]] = 1/prob;
                            else iter->second += 1/prob;

                            auto iter_succ = catalog_->success_map_iso.find(embedding_[0]);
                            if(iter_succ == catalog_->success_map_iso.end()) catalog_->success_map_iso[embedding_[0]] = 1;
                            else iter_succ->second += 1;
                        }
                        else{
                            Edge ptb = {embedding_[0],embedding_[1]};
                            auto iter = catalog_->res_map.find(ptb);
                            if(iter == catalog_->res_map.end()) catalog_->res_map[ptb] = 1/prob;
                            else iter->second += 1/prob;

                            //#success
                            auto iter_succ = catalog_->success_map.find(ptb);
                            if(iter_succ == catalog_->success_map.end()) catalog_->success_map[ptb] = 1;
                            else iter_succ->second += 1;
                        }
                    }
                    
                    valid_nc = false;
                    break;
                }
                else{
                    cur_depth += 1;
                }
           }
        }   
    }
    
    //reset the "visited"
    for (uint32_t j = start_depth; j < ext_depth; ++j) {
        visited_[embedding_[j]] = false;
    }

    //reset the "embedding"
    for(uint32_t j = start_depth; j < ext_depth; ++j) {
        embedding_[j] = 0;
    }   

    return false;
}

uint64_t leapfrogtriejoin::execute_rwdw_final(){
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = catalog_->combine_core;
    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];

    uint32_t start_depth;
    if(catalog_->iso) start_depth = 1;
    else start_depth = 2;

    std::random_device rd;
    std:mt19937 gen(rd());
    std::discrete_distribution<> d_big_ori(catalog_->weight_vec_big.begin(), catalog_->weight_vec_big.end());

    while (true){
        if(g_exit)
            return count_;

        bool valid = true;
        prob = 1;
        if(catalog_->iso){
            uint32_t random_ptb = d_big_ori(gen);
            embedding_[0] = catalog_->iso_vec[random_ptb];
            auto iter = catalog_->freq_map_iso.find(embedding_[0]);
            if(iter == catalog_->freq_map_iso.end()) catalog_->freq_map_iso[embedding_[0]] = 1;
            else iter->second++; 
        }
        else{
            uint32_t random_ptb = d_big_ori(gen);
            embedding_[0] = catalog_->node_small[random_ptb];
            embedding_[1] = catalog_->node_big[random_ptb];
            pEdge ptb = {embedding_[0],embedding_[1]};
            auto iter = catalog_->freq_map.find(ptb);
            if(iter == catalog_->freq_map.end()) catalog_->freq_map[ptb] = 1;
            else iter->second ++;
        }
        for (uint32_t j = 0; j < start_depth; ++j) {
            uint32_t u = vertex_ordering_[j];
            visited_[embedding_[j]] = true;
            idx_embedding_[j] = catalog_->get_candidate_index(u, embedding_[j]);
        }
        uint32_t cur_depth = start_depth;
        uint32_t ext_depth = start_depth;

        if(start_depth >= num_core_vertex_){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }

            if(num_local_candidates_[start_depth] > 0) {
                bool flag = rwcw_non_core_final(start_depth,max_depth);
            }
        }
        else{
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);
            while(valid){
                while (cur_depth <= max_depth - 2){
                    if(num_local_candidates_[cur_depth] == 0) {
                        valid = false;
                        break;
                    }
                    prob = prob / (double)num_local_candidates_[cur_depth];

                    uint32_t u = vertex_ordering_[cur_depth];
                    std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                    int random = d(gen);

                    uint32_t encoded_id = local_candidates_[cur_depth][random];
                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        valid = false;
                        break;
                    }
                    visited_[v] = true;
                    ext_depth++;
                    embedding_[cur_depth] = v;

                    num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);
                    if(cur_depth == max_depth -2){
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        if(loop_count == 0) {
                            valid = false;
                            break;
                        }
                        prob = prob / (double)num_local_candidates_[cur_depth + 1];

                        std::uniform_int_distribution<> d(0,loop_count - 1);
                        int random = d(gen);

                        uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[random]];
                        if(visited_[temp_v]){
                            valid = false;
                            break;
                        }
                        embedding_[cur_depth + 1] = temp_v;
                        ext_depth++;
                        if(catalog_->iso){
                            auto iter = catalog_->res_map_iso.find(embedding_[0]);
                            if(iter == catalog_->res_map_iso.end()) catalog_->res_map_iso[embedding_[0]] = 1/prob;
                            else iter->second += 1/prob;

                            auto iter_succ = catalog_->success_map_iso.find(embedding_[0]);
                            if(iter_succ == catalog_->success_map_iso.end()) catalog_->success_map_iso[embedding_[0]] = 1;
                            else iter_succ->second++;
                        }
                        else{
                            Edge ptb = {embedding_[0],embedding_[1]};
                            auto iter = catalog_->res_map.find(ptb);
                            if(iter == catalog_->res_map.end()) catalog_->res_map[ptb] = 1/prob;
                            else iter->second += 1/prob;

                            auto iter_succ = catalog_->success_map.find(ptb);
                            if(iter_succ == catalog_->success_map.end()) catalog_->success_map[ptb] = 1;
                            else iter_succ->second++;
                        }
                        valid = false;
                        break;
                    }
                    else if (cur_depth == core_vertex_depth - 2){
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];
                        if(num_local_candidates_[next_depth] == 0) {
                            valid = false;
                            break;
                        }

                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }

                        bool exit = rwcw_non_core_final(next_depth, max_depth);
                        if(exit || !exit) {
                            valid = false;
                            break;
                        }
                    }
                    else{
                        cur_depth += 1;
                        if(cur_depth > max_depth - 2) valid = false;
                    }
                }
            }
        }

        for(int k = 0; k < ext_depth;++k){
            visited_[embedding_[k]] = false;
        }

        for(int k = 0; k < ext_depth;++k){
            embedding_[k] = 0;
        }

        catalog_->total_sp_num++;

        if (g_exit)
            return count_;
    }

    return count_;
}

bool leapfrogtriejoin::rwcw_non_core_epl(uint32_t start_depth, uint32_t max_depth){
    uint32_t ext_depth = start_depth;

    std::random_device rd;
    std:mt19937 gen_nc(rd());

    std::uniform_int_distribution<> d(0,num_local_candidates_[start_depth] - 1);
    int random = d(gen_nc);

    uint32_t v = si_buffer_[random];
    if(visited_[v]){
        return false;
    }
    embedding_[start_depth] = v;
    ext_depth++;

    prob = prob / (double)num_local_candidates_[start_depth];

    uint32_t cur_depth = start_depth + 1;
    uint32_t u = vertex_ordering_[cur_depth];
    uint32_t bn = get_backward_neighbors(cur_depth)[0];
    uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
    uint32_t key = embedding_[bns_depth];
    local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);

    bool valid_nc = true;

    if (cur_depth == max_depth - 1) {
        uint32_t temp_count = num_local_candidates_[cur_depth];
        uint32_t* temp_buffer = local_candidates_[cur_depth];
        if(temp_count > 0) {

            prob = prob / (double)temp_count;
            
            std::uniform_int_distribution<> d(0,temp_count - 1);
            int random = d(gen_nc);

            uint32_t temp_v = temp_buffer[random];
            if(!visited_[temp_v]){
                embedding_[cur_depth] = temp_v;

                uint32_t start_node = embedding_[0];

                auto iter_res = catalog_->res_node.find(start_node);
                if(iter_res == catalog_->res_node.end()) catalog_->res_node[start_node] = 1/prob;
                else iter_res->second += 1/prob;

                auto iter_succ = catalog_->res_success.find(start_node);
                if(iter_succ == catalog_->res_success.end()) catalog_->res_success[start_node] = 1;
                else iter_succ->second++;
            }
        }
    }
    else{
        while (valid_nc)
        {
           while (cur_depth <= max_depth - 2){
                if(num_local_candidates_[cur_depth] == 0) {
                    valid_nc = false;
                    break;
                }
                prob = prob / (double)num_local_candidates_[cur_depth];

                std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                int random = d(gen_nc);

                uint32_t v = local_candidates_[cur_depth][random];
                
                if (visited_[v]) {
                    valid_nc = false;
                    break;
                }

                visited_[v] = true;
                embedding_[cur_depth] = v;
                ext_depth++;

                uint32_t next_depth = cur_depth + 1;
                uint32_t u = vertex_ordering_[next_depth];
                uint32_t bn = get_backward_neighbors(next_depth)[0];
                uint32_t bns_depth = get_backward_neighbors_depth(next_depth)[0]; 
                uint32_t key = embedding_[bns_depth];
                local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);

                if(cur_depth == max_depth - 2){
                    uint32_t temp_count = num_local_candidates_[next_depth];
                    uint32_t *temp_buffer = local_candidates_[next_depth];    
                    if(temp_count == 0) {
                        valid_nc = false;
                        break;
                    }
                    prob = prob / (double)temp_count;

                    std::uniform_int_distribution<> d(0,temp_count - 1);
                    int random = d(gen_nc);

                    uint32_t temp_v = temp_buffer[random];
                    if(!visited_[temp_v]){
                        embedding_[next_depth] = temp_v;
                        uint32_t start_node = embedding_[0];
                        auto iter_res = catalog_->res_node.find(start_node);
                        if(iter_res == catalog_->res_node.end()) catalog_->res_node[start_node] = 1/prob;
                        else iter_res->second += 1/prob;

                        auto iter_succ = catalog_->res_success.find(start_node);
                        if(iter_succ == catalog_->res_success.end()) catalog_->res_success[start_node] = 1;
                        else iter_succ->second++;
                    }
                    
                    valid_nc = false;
                    break;
                }
                else{
                    cur_depth += 1;
                }
           }
        }   
    }
    
    //reset the "visited"
    for (uint32_t j = start_depth; j < ext_depth; ++j) {
        visited_[embedding_[j]] = false;
    }

    //reset the "embedding"
    for(uint32_t j = start_depth; j < ext_depth; ++j) {
        embedding_[j] = 0;
    }   

    return false;    
}

uint64_t leapfrogtriejoin::execute_rwdw_epl(){
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = catalog_->combine_core;

    intermediate_result_count_[0] = input_->get_size();

    uint32_t init_size = input_->get_size();

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    
    std::random_device rd;
    std:mt19937 gen(rd());
    std::uniform_int_distribution<> d(0,init_size - 1);

    while(true){
        bool valid = true;
        int random = d(gen);
        memcpy(embedding_, input_->get_tuple(random), sizeof(uint32_t) * start_depth);

        auto iter_freq = catalog_->freq_sp_exp.find(embedding_[0]);
        if(iter_freq == catalog_->freq_sp_exp.end()) catalog_->freq_sp_exp[embedding_[0]] = 1;
        else iter_freq->second++;

        set_idx(start_depth);

        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }

        prob = 1 / (double)init_size;

        uint32_t cur_depth = start_depth;
        uint32_t ext_depth = start_depth;

        if(start_depth == max_depth - 1){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
            int random = d(gen);

            prob /= (double)num_local_candidates_[cur_depth];
            embedding_[start_depth] = local_candidates_[start_depth][random];
        
            uint32_t start_node = embedding_[0];

            auto iter_res = catalog_->res_node.find(start_node);
            if(iter_res == catalog_->res_node.end()) catalog_->res_node[start_node] = 1/prob;
            else iter_res->second += 1/prob;

            auto iter_succ = catalog_->res_success.find(start_node);
            if(iter_succ == catalog_->res_success.end()) catalog_->res_success[start_node] = 1;
            else iter_succ->second++;
        }
        else if(start_depth >= num_core_vertex_){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }
            if(num_local_candidates_[start_depth] > 0) bool exit = rwcw_non_core_epl(start_depth, max_depth);
        }
        else{
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);
            while(valid){
                while (cur_depth <= max_depth - 2){
                    if(num_local_candidates_[cur_depth] == 0) {
                        valid = false;
                        break;
                    }
                    prob = prob / (double)num_local_candidates_[cur_depth];

                    uint32_t u = vertex_ordering_[cur_depth];
                    std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                    int random = d(gen);

                    uint32_t encoded_id = local_candidates_[cur_depth][random];
                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        valid = false;
                        break;
                    }
                    visited_[v] = true;
                    ext_depth++;
                    embedding_[cur_depth] = v;

                    num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                    if(cur_depth == max_depth - 2){
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        if(loop_count == 0) {
                            valid = false;
                            break;
                        }
                        prob = prob / (double)num_local_candidates_[cur_depth + 1];

                        std::uniform_int_distribution<> d(0,loop_count - 1);
                        int random = d(gen);

                        uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[random]];
                        if(!visited_[temp_v]){
                            embedding_[cur_depth + 1] = temp_v;
                            uint32_t start_node = embedding_[0];
                            auto iter_res = catalog_->res_node.find(start_node);
                            if(iter_res == catalog_->res_node.end()) catalog_->res_node[start_node] = 1/prob;
                            else iter_res->second += 1/prob;

                            auto iter_succ = catalog_->res_success.find(start_node);
                            if(iter_succ == catalog_->res_success.end()) catalog_->res_success[start_node] = 1;
                            else iter_succ->second++;

                            valid = false;
                            break;
                        }

                        valid = false;
                        break;
                    }
                    else if(cur_depth == core_vertex_depth - 2){
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];
                        if(num_local_candidates_[next_depth] == 0) {
                            valid = false;
                            break;
                        }

                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }

                        bool exit = rwcw_non_core_epl(next_depth, max_depth);
                        
                        valid = false;
                        break;
                    }
                    else{
                        cur_depth += 1;
                        if(cur_depth > max_depth - 2) valid = false;
                    }
                }
            }
        }
        for(int k = 0; k < ext_depth;++k){
            visited_[embedding_[k]] = false;
        }

        for(int k = 0; k < ext_depth;++k){
            embedding_[k] = 0;
        }

        catalog_->total_sp_num++;
        
        if (g_exit)
            return count_;
    }
    return count_;
}

uint32_t leapfrogtriejoin::compute_local_candidates_noncore_ptb_nc(uint32_t depth){

    uint32_t exist_count = 0,lc_count = 0;

    uint32_t disconnect_node = vertex_ordering_[catalog_->end_depth];
    uint32_t bn = vertex_ordering_[catalog_->begin_depth];
    uint32_t u = vertex_ordering_[catalog_->end_depth];
    uint32_t key = idx_embedding_[catalog_->begin_depth];
    uint32_t *lc_exist = catalog_->get_non_core_relation_children(bn, u, key, exist_count);

    if(catalog_->query_graph_->getVertexDegree(disconnect_node) > 0){

        uint32_t* candidates = catalog_->get_candidates(disconnect_node);
        uint32_t num_candidates = catalog_->get_num_candidates(disconnect_node);

        for(uint32_t i = 0; i < num_candidates; ++i){
            if(!exist_check(candidates[i],lc_exist,exist_count)){
                local_candidates_[depth][lc_count] = candidates[i];
                lc_count++;
            }
        }

    }
    else{
        uint32_t label = catalog_->query_graph_->getVertexLabel(disconnect_node);
        uint32_t cnt;
        const ui* last_v = catalog_->data_graph_->getVerticesByLabel(label,cnt);

        for(uint32_t i = 0; i < cnt; ++i){
            if(!exist_check(last_v[i],lc_exist,exist_count)){
                local_candidates_[depth][lc_count] = last_v[i];
                lc_count++;
            }
        }
    }

    return lc_count;
}

bool leapfrogtriejoin::rwd_non_core(uint32_t start_depth, uint32_t max_depth){
    uint32_t ext_depth = start_depth;

    std::random_device rd;
    std:mt19937 gen_nc(rd());

    std::uniform_int_distribution<> d(0,num_local_candidates_[start_depth] - 1);
    int random = d(gen_nc);

    uint32_t v = si_buffer_[random];
    if(visited_[v]){
        return false;
    }
    embedding_[start_depth] = v;
    ext_depth++;

    if(start_depth != catalog_->end_depth) prob = prob / (double)num_local_candidates_[start_depth];

    uint32_t cur_depth = start_depth + 1;

    if(catalog_->end_depth == cur_depth){
        num_local_candidates_[cur_depth] = compute_local_candidates_noncore_ptb_nc(catalog_->end_depth);
    }
    else{
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];
        local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);
    }

    bool valid_nc = true;
    
    if (cur_depth == max_depth - 1) {
        uint32_t temp_count = num_local_candidates_[cur_depth];
        uint32_t* temp_buffer = local_candidates_[cur_depth];
        if(temp_count == 0) {
            valid_nc = false;
        }
        else{
            if(catalog_->end_depth != cur_depth) prob = prob / (double)temp_count;

            std::uniform_int_distribution<> d(0,temp_count - 1);
            int random = d(gen_nc);

            uint32_t temp_v = temp_buffer[random];
            if(!visited_[temp_v]){
                embedding_[cur_depth] = temp_v;
                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                auto iter = this->catalog_->res_map.find(ptb);
                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = (uint64_t)(1 / prob);
                else iter->second += (uint64_t)(1 / prob);

                auto iter_success = this->catalog_->success_map.find(ptb);
                if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                else iter_success->second += 1;
            }
        }
    }
    else{
        while (valid_nc)
        {
            while (cur_depth <= max_depth - 2){
                if(num_local_candidates_[cur_depth] == 0) {
                    valid_nc = false;
                    break;
                }
                if(cur_depth != catalog_->end_depth) prob = prob / (double)num_local_candidates_[cur_depth];
                
                std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                int random = d(gen_nc);

                uint32_t v = local_candidates_[cur_depth][random];
                
                if (visited_[v]) {
                    valid_nc = false;
                    break;
                }

                visited_[v] = true;
                embedding_[cur_depth] = v;
                ext_depth++;

                uint32_t next_depth = cur_depth + 1;

                if(this->catalog_->end_depth == next_depth){
                    num_local_candidates_[cur_depth + 1] = compute_local_candidates_noncore_ptb_nc(catalog_->end_depth);
                } 
                else{
                    uint32_t u = vertex_ordering_[next_depth];
                    uint32_t bn = get_backward_neighbors(next_depth)[0];
                    uint32_t bns_depth = get_backward_neighbors_depth(next_depth)[0]; 
                    uint32_t key = embedding_[bns_depth];
                    local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);
                }

                if(cur_depth == max_depth - 2){
                    uint32_t temp_count = num_local_candidates_[next_depth];
                    uint32_t *temp_buffer = local_candidates_[next_depth];    
                    if(temp_count == 0) {
                        valid_nc = false;
                        break;
                    }
                    else{
                        if(cur_depth + 1 != catalog_->end_depth) prob = prob / (double)temp_count;

                        std::uniform_int_distribution<> d(0,temp_count - 1);
                        int random = d(gen_nc);

                        uint32_t temp_v = temp_buffer[random];
                        if(!visited_[temp_v]){
                            embedding_[next_depth] = temp_v;
                            pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                            auto iter = this->catalog_->res_map.find(ptb);
                            if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = (uint64_t)(1 / prob);
                            else iter->second += (uint64_t)(1 / prob);

                            auto iter_success = this->catalog_->success_map.find(ptb);
                            if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                            else iter_success->second += 1;
                        }
                        
                        valid_nc = false;
                        break;
                    }
                }
                else{
                    cur_depth += 1;
                } 
            }
        }
    }

    //reset the "visited"
    for (uint32_t j = start_depth; j < ext_depth; ++j) {
        visited_[embedding_[j]] = false;
    }

    //reset the "embedding"
    for(uint32_t j = start_depth; j < ext_depth; ++j) {
        embedding_[j] = 0;
    }    

    return false;
}

uint32_t leapfrogtriejoin::compute_local_candidates_core_ptb_nc(uint32_t depth) {

    uint32_t disconnect_node = vertex_ordering_[catalog_->end_depth];
    uint32_t lc_count = 0;
    uint32_t bn =  vertex_ordering_[catalog_->begin_depth];
    uint32_t u = vertex_ordering_[catalog_->end_depth];
    uint32_t key = idx_embedding_[catalog_->begin_depth];
    uint32_t exist_count = 0;
    uint32_t *lc_exist = catalog_->get_core_relation_children(bn, u, key, exist_count);

    if(catalog_->query_graph_->getVertexDegree(disconnect_node) > 0){
        uint32_t num_candidates = catalog_->get_num_candidates(disconnect_node);
        for(uint32_t i = 0; i < num_candidates; ++i){
            if(!exist_check(i,lc_exist,exist_count)){
                local_candidates_[depth][lc_count] = i;
                lc_count++;
            }
        }
    }

    if (lc_count == 0){
        EXIT_EMPTY_SET:
#ifdef COLLECT_FAIL_SI_STATISTICS
        fail_si_count_ += 1;
#endif

        return 0;
    }

    return lc_count;
}

uint64_t leapfrogtriejoin::execute_rwd(){
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = catalog_->combine_core;

    intermediate_result_count_[0] = input_->get_size();

    uint32_t init_size = input_->get_size();

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    
    std::random_device rd;
    std:mt19937 gen(rd());
    std::uniform_int_distribution<> d(0,init_size - 1);


    while(true){

        bool valid = true;
        int random = d(gen);
        memcpy(embedding_, input_->get_tuple(random), sizeof(uint32_t) * start_depth);

        set_idx(start_depth);

        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }

        prob = 1 / (double)init_size;
        uint32_t cur_depth = start_depth;
        uint32_t ext_depth = start_depth;

        if(start_depth >= num_core_vertex_){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            if(start_depth == catalog_->end_depth) num_local_candidates_[cur_depth] = compute_local_candidates_noncore_ptb_nc(catalog_->end_depth);
            else local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }

            if(num_local_candidates_[start_depth] > 0) bool flag = rwd_non_core(start_depth,max_depth);
        
        }
        else{
            if(cur_depth == catalog_->end_depth) num_local_candidates_[cur_depth] = compute_local_candidates_core_ptb_nc(cur_depth);
            else  num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);
            while(valid){
                while (cur_depth <= max_depth - 2){
                    if(num_local_candidates_[cur_depth] == 0) {
                        valid = false;
                        break;
                    }
                    if(cur_depth != catalog_->end_depth) prob = prob / (double)num_local_candidates_[cur_depth];

                    uint32_t u = vertex_ordering_[cur_depth];
                    std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                    int random = d(gen);

                    uint32_t encoded_id = local_candidates_[cur_depth][random];
                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        valid = false;
                        break;
                    }
                    
                    visited_[v] = true;
                    ext_depth++;
                    embedding_[cur_depth] = v;

                    if(cur_depth + 1 == this->catalog_->end_depth) num_local_candidates_[cur_depth + 1] = compute_local_candidates_core_ptb_nc(cur_depth + 1);
                    else num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                    if(cur_depth == max_depth -2){
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        if(loop_count == 0) {
                            valid = false;
                            break;
                        }
                        if(cur_depth + 1 != catalog_->end_depth) prob = prob / (double)num_local_candidates_[cur_depth + 1];

                        std::uniform_int_distribution<> d(0,loop_count - 1);
                        int random = d(gen);

                        uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[random]];
                        if(visited_[temp_v]){
                            valid = false;
                            break;
                        }

                        embedding_[cur_depth + 1] = temp_v;
                        ext_depth++;       

                        pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                        auto iter = this->catalog_->res_map.find(ptb);
                        if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = (uint64_t)(1 / prob);
                        else iter->second += (uint64_t)(1 / prob);

                        auto iter_success = this->catalog_->success_map.find(ptb);
                        if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                        else iter_success->second += 1;

                        valid = false;
                        break;
                    }
                    else if(cur_depth == core_vertex_depth - 2){
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];
                        if(num_local_candidates_[next_depth] == 0) {
                            valid = false;
                            break;
                        }
                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }
                        bool flag = rwd_non_core(next_depth,max_depth);

                        valid = false;
                        break;
                    }
                    else{
                        cur_depth += 1;
                        if(cur_depth > max_depth - 2) valid = false;
                    }
                }
            }
        }

        for(int k = 0; k < ext_depth;++k){
            visited_[embedding_[k]] = false;
        }

        for(int k = 0; k < ext_depth;++k){
            embedding_[k] = 0;
        }

        catalog_->total_sp_num++;
        
        if (g_exit)
            return count_;
    }


    return count_;
}

void leapfrogtriejoin::rwch_non_core_sampling(uint32_t start_depth, uint32_t max_depth){
    uint32_t init_size = num_local_candidates_[start_depth];
    std::random_device rd;
    std:mt19937 gen_nc(rd());
    std::uniform_int_distribution<> d(0,init_size - 1);
    
    for(int i = 0; i < catalog_->sampling_num; ++i){
        this->prob = 1.0;
        bool valid_nc = true;
        int random = d(gen_nc);
        uint32_t v = local_candidates_[start_depth][random];
        
        uint32_t exp_depth = start_depth;

        if (visited_[v]) {
            valid_nc = false;
            continue;
        }
        visited_[v] = true;
        embedding_[start_depth] = v;
        exp_depth++;
        prob = prob / (double)init_size;

        //compute ptb candidates
        if(end_node_bn.test(start_depth)) num_ptb_candidates_[start_depth] = pre_compute_ptb_candidates_nc(start_depth);
        else num_ptb_candidates_[start_depth] = 0;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        //early termination
        if(end_node_bn.test(start_depth) && num_ptb_candidates_[start_depth] == 0){
            if(!this->begin_node_first || cur_depth != catalog_->begin_depth){
                visited_[embedding_[start_depth]] = false;
                continue;
            }
        }

        if(cur_depth == catalog_->end_depth) {
            uint32_t last_cnt = 0;
            uint32_t last_depth = 0;
            for(int k = cur_depth - 1; k >= 0; k--){
                if(end_node_bn.test(k)){
                    last_depth = k;
                    last_cnt = num_ptb_candidates_[k];
                    break;
                }
            }
            for(int k = 0; k < last_cnt;++k){
                local_candidates_[cur_depth][k] = ptb_candidates_[last_depth][k];
            }
            num_local_candidates_[cur_depth] = last_cnt;
        }
        else local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);

        if(cur_depth == max_depth - 1){
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            if(temp_count == 0) {
                valid_nc = false;
            }
            else{
                prob = prob / temp_count;

                std::uniform_int_distribution<> d(0,temp_count - 1);
                int random = d(gen_nc);

                uint32_t temp_v = temp_buffer[random];

                if(visited_[temp_v]) {
                    valid_nc = false;
                }
                else{
                    embedding_[cur_depth] = temp_v;
                    pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};

                    auto iter = this->res_tmp.find(ptb);
                    if(iter == this->res_tmp.end()) this->res_tmp[ptb] = (1 / prob);
                    else iter->second += (1 / prob);
                    
                    auto iter_success = this->catalog_->success_map.find(ptb);
                    if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                    else iter_success->second++;
                }
            }
        }
        else{
            while(valid_nc){
                while (cur_depth <= max_depth - 2){
                    if(num_local_candidates_[cur_depth] == 0) {
                        //store failed result
                        valid_nc = false;
                        break;
                    }
                    
                    prob = prob / num_local_candidates_[cur_depth];
                    
                    std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                    int random = d(gen_nc);

                    uint32_t v = local_candidates_[cur_depth][random];

                    if (visited_[v]) {
                        valid_nc = false;
                        break;
                    }
                    embedding_[cur_depth] = v;
                    visited_[v] = true;
                    exp_depth++;

                    if(end_node_bn.test(cur_depth)) num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates_nc(cur_depth);
                    else num_ptb_candidates_[cur_depth] = 0;
                    
                    uint32_t next_depth = cur_depth + 1;
                    u = vertex_ordering_[next_depth];
                    bn = get_backward_neighbors(next_depth)[0];
                    bns_depth = get_backward_neighbors_depth(next_depth)[0]; 
                    key = embedding_[bns_depth];

                    if(next_depth == catalog_->end_depth){
                        uint32_t last_cnt = 0;
                        uint32_t last_depth = 0;
                        for(int k = cur_depth; k >= 0; k--){
                            if(end_node_bn.test(k)){
                                last_depth = k;
                                last_cnt = num_ptb_candidates_[k];
                                break;
                            }
                        }
                        for(int k = 0; k < last_cnt;++k){
                            local_candidates_[next_depth][k] = ptb_candidates_[last_depth][k];
                        }
                        num_local_candidates_[next_depth] = last_cnt;
                    } 
                    else local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);


                    //early termination
                    if(end_node_bn.test(cur_depth) && num_ptb_candidates_[cur_depth] == 0){
                        if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                           num_local_candidates_[next_depth] = 0;
                    }

                    if(cur_depth == max_depth - 2){
                        uint32_t temp_count = num_local_candidates_[next_depth];
                        uint32_t *temp_buffer = local_candidates_[next_depth];    
                        if(temp_count == 0) {
                            valid_nc = false;
                            break;
                        }else{
                            prob = prob / (double)temp_count;

                            std::uniform_int_distribution<> d(0,temp_count - 1);
                            int random = d(gen_nc);

                            uint32_t temp_v = temp_buffer[random];

                            if(visited_[temp_v]){
                                valid_nc = false;
                                break;
                            }
                            else{
                                embedding_[next_depth] = temp_v;
                                //store success sample
                                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};

                                auto iter = this->res_tmp.find(ptb);
                                if(iter == this->res_tmp.end()) this->res_tmp[ptb] = (1 / prob);
                                else iter->second += (1 / prob);

                                auto iter_success = this->catalog_->success_map.find(ptb);
                                if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                                else iter_success->second++;

                                valid_nc = false;
                                break;
                            }
                        }
                    }
                    else{
                        cur_depth += 1;
                    }
                }
            }
        }
    
        //reset the "visited"
        for (uint32_t j = start_depth; j < exp_depth; ++j) {
            visited_[embedding_[j]] = false;
        }

        //reset the "embedding"
        for(uint32_t j = start_depth; j < exp_depth; ++j) {
            embedding_[j] = 0;
        }
    }

}

bool leapfrogtriejoin::rwch_non_core(uint32_t start_depth, uint32_t max_depth){
    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {
        uint32_t v = si_buffer_[idx_[start_depth]];
        if (visited_[v]) {
            continue;
        }
        visited_[v] = true;
        embedding_[start_depth] = v;

        //compute ptb candidates
        if(end_node_bn.test(start_depth)) num_ptb_candidates_[start_depth] = pre_compute_ptb_candidates_nc(start_depth);
        else num_ptb_candidates_[start_depth] = 0;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        if(cur_depth == catalog_->end_depth) {
            uint32_t last_cnt = 0;
            uint32_t last_depth = 0;
            for(int k = cur_depth - 1; k >= 0; k--){
                if(end_node_bn.test(k)){
                    last_depth = k;
                    last_cnt = num_ptb_candidates_[k];
                    break;
                }
            }
            for(int k = 0; k < last_cnt;++k){
                local_candidates_[cur_depth][k] = ptb_candidates_[last_depth][k];
            }
            num_local_candidates_[cur_depth] = last_cnt;
        }
        else local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);

        //early termination
        if(end_node_bn.test(start_depth) && num_ptb_candidates_[start_depth] == 0){
            if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                num_local_candidates_[cur_depth] = 0;
        }

        if (cur_depth == max_depth - 1) {
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                uint32_t temp_v = temp_buffer[j];
                if (visited_[temp_buffer[j]]) {
                    continue;
                }
                embedding_[cur_depth] = temp_v;
                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                auto iter = this->catalog_->res_map.find(ptb);
                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                else iter->second++;
            }
            visited_[v] = false;
        }
        else if(num_local_candidates_[cur_depth] > 0){
            rwch_non_core_sampling(cur_depth, max_depth);
            for(auto iter = res_tmp.begin(); iter != res_tmp.end(); ++iter){
                auto iter_res = this->catalog_->res_map.find(iter->first);
                if(iter_res == this->catalog_->res_map.end()) this->catalog_->res_map[iter->first] = iter->second / catalog_->sampling_num;
                else iter_res->second += iter->second / catalog_->sampling_num;
            }
            res_tmp.clear();
        }
        visited_[embedding_[start_depth]] = false;
    }
    return false;
}

bool leapfrogtriejoin::rwch_non_core_single(uint32_t start_depth, uint32_t max_depth){
    std::random_device rd;
    std:mt19937 gen_nc(rd());

    std::uniform_int_distribution<> d(0,num_local_candidates_[start_depth] - 1);
    int random = d(gen_nc);
    uint32_t exp_depth = start_depth;

    uint32_t v = si_buffer_[random];

    if(visited_[v]){
        return false;
    }

    visited_[v] = true;
    embedding_[start_depth] = v;
    exp_depth++;

    prob = prob / num_local_candidates_[start_depth];

    //compute ptb candidates
    if(end_node_bn.test(start_depth)) num_ptb_candidates_[start_depth] = pre_compute_ptb_candidates_nc(start_depth);
    else num_ptb_candidates_[start_depth] = 0;

    uint32_t cur_depth = start_depth + 1;
    uint32_t u = vertex_ordering_[cur_depth];
    uint32_t bn = get_backward_neighbors(cur_depth)[0];
    uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
    uint32_t key = embedding_[bns_depth];
    
    if(cur_depth == catalog_->end_depth) {
        uint32_t last_cnt = 0;
        uint32_t last_depth = 0;
        for(int k = cur_depth - 1; k >= 0; k--){
            if(end_node_bn.test(k)){
                last_depth = k;
                last_cnt = num_ptb_candidates_[k];
                break;
            }
        }
        for(int k = 0; k < last_cnt;++k){
            local_candidates_[cur_depth][k] = ptb_candidates_[last_depth][k];
        }
        num_local_candidates_[cur_depth] = last_cnt;
    }
    else local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);
        
    //early termination
    if(end_node_bn.test(start_depth) && num_ptb_candidates_[start_depth] == 0){
        if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
            num_local_candidates_[cur_depth] = 0;
    }

    bool valid_nc = true;
    if (cur_depth == max_depth - 1) {
        uint32_t temp_count = num_local_candidates_[cur_depth];
        uint32_t* temp_buffer = local_candidates_[cur_depth];
        if(temp_count == 0) {
            valid_nc = false;
        }
        else{
            std::uniform_int_distribution<> d(0,temp_count - 1);
            int random = d(gen_nc);

            uint32_t temp_v = temp_buffer[random];
            prob = prob / temp_count;

            if(!visited_[temp_v]){
                embedding_[cur_depth] = temp_v;
                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};

                auto iter = this->res_tmp.find(ptb);
                if(iter == this->res_tmp.end()) this->res_tmp[ptb] = (uint64_t)(1 / prob);
                else iter->second += (uint64_t)(1 / prob);
                
                auto iter_success = this->catalog_->success_map.find(ptb);
                if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                else iter_success->second++;
            }
        }
    }
    else{
        while (valid_nc){
            while (cur_depth <= max_depth - 2){
                if(num_local_candidates_[cur_depth] == 0) {
                    valid_nc = false;
                    break;
                }
                prob = prob / num_local_candidates_[cur_depth];

                std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                int random = d(gen_nc);

                uint32_t v = local_candidates_[cur_depth][random];
                
                if (visited_[v]) {
                    valid_nc = false;
                    break;
                }

                visited_[v] = true;
                embedding_[cur_depth] = v;
                exp_depth++;

                if(end_node_bn.test(cur_depth)) num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates_nc(cur_depth);
                else num_ptb_candidates_[cur_depth] = 0;

                uint32_t next_depth = cur_depth + 1;
                u = vertex_ordering_[next_depth];
                bn = get_backward_neighbors(next_depth)[0];
                bns_depth = get_backward_neighbors_depth(next_depth)[0]; 
                key = embedding_[bns_depth];

                if(next_depth == catalog_->end_depth){
                    uint32_t last_cnt = 0;
                    uint32_t last_depth = 0;
                    for(int k = cur_depth; k >= 0; k--){
                        if(end_node_bn.test(k)){
                            last_depth = k;
                            last_cnt = num_ptb_candidates_[k];
                            break;
                        }
                    }
                    for(int k = 0; k < last_cnt;++k){
                        local_candidates_[next_depth][k] = ptb_candidates_[last_depth][k];
                    }
                    num_local_candidates_[next_depth] = last_cnt;
                } 
                else local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);
                
                //early termination
                if(end_node_bn.test(cur_depth) && num_ptb_candidates_[cur_depth] == 0){
                    if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                        num_local_candidates_[next_depth] = 0;
                }

                if(cur_depth == max_depth - 2){
                    uint32_t temp_count = num_local_candidates_[next_depth];
                    uint32_t *temp_buffer = local_candidates_[next_depth];    
                    if(temp_count == 0) {
                        valid_nc = false;
                        break;
                    }
                    else{
                        prob = prob / temp_count;

                        std::uniform_int_distribution<> d(0,temp_count - 1);
                        int random = d(gen_nc);

                        uint32_t temp_v = temp_buffer[random];

                        if(!visited_[temp_v]){
                            embedding_[next_depth] = temp_v;
                            pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};

                            auto iter = this->res_tmp.find(ptb);
                            if(iter == this->res_tmp.end()) this->res_tmp[ptb] = (uint64_t)(1 / prob);
                            else iter->second += (uint64_t)(1 / prob);
                        
                            auto iter_success = this->catalog_->success_map.find(ptb);
                            if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                            else iter_success->second++;
                        }
                        
                        valid_nc = false;
                        break;
                    }
                }
                else{
                    cur_depth += 1;
                }

            }
        }
    }

    //reset the "visited"
    for (uint32_t j = start_depth; j < exp_depth; ++j) {
        visited_[embedding_[j]] = false;
    }

    //reset the "embedding"
    for(uint32_t j = start_depth; j < exp_depth; ++j) {
        embedding_[j] = 0;
    }    

    return false;
}

bool leapfrogtriejoin::rwch_core(uint32_t start_depth, uint32_t max_depth){
    std::random_device rd;
    std:mt19937 gen(rd());

    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {
        uint32_t u = vertex_ordering_[start_depth];
        uint32_t encoded_id = local_candidates_[start_depth][idx_[start_depth]];
        uint32_t v = catalog_->candidate_sets_[u][encoded_id];

        if (visited_[v]) {
            continue;
        }        

        visited_[v] = true;
        idx_embedding_[start_depth] = encoded_id;
        embedding_[start_depth] = v;

        //compute ptb candidates
        if(end_node_bn.test(start_depth)) {
            if(catalog_->end_depth < num_core_vertex_) num_ptb_candidates_[start_depth] = pre_compute_ptb_candidates(start_depth);
            else num_ptb_candidates_[start_depth] = pre_compute_ptb_candidates_nc(start_depth);
        }
        else num_ptb_candidates_[start_depth] = 0;

        uint32_t cur_depth = start_depth + 1;

        if(cur_depth == this->catalog_->end_depth) {
            uint32_t last_cnt = 0;
            uint32_t last_depth = 0;
            for(int k = cur_depth; k >= 0; k--){
                if(end_node_bn.test(k)){
                    last_depth = k;
                    last_cnt = num_ptb_candidates_[k];
                    break;
                }
            }
            for(int k = 0; k < last_cnt;++k){
                local_candidates_[cur_depth][k] = ptb_candidates_[last_depth][k];
            }
            num_local_candidates_[cur_depth] = last_cnt;
        }
        else num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

        //early termination
        if(end_node_bn.test(start_depth) && num_ptb_candidates_[start_depth] == 0){
            if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                num_local_candidates_[cur_depth] = 0;
        }

        if(start_depth == max_depth - 2 && num_local_candidates_[cur_depth] > 0){
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];

            for(uint32_t j = 0; j < temp_count; ++j){
                uint32_t temp_v = temp_buffer[j];
                if(!visited_[temp_v]){
                    embedding_[cur_depth] = temp_v;
                    pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                    auto iter = this->catalog_->res_map.find(ptb);
                    if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                    else iter->second++;
                }
            }
        }
        else if(start_depth == num_core_vertex_ - 2 && num_local_candidates_[cur_depth] > 0){
            uint32_t next_depth = start_depth + 1;
            uint32_t next_vertex = vertex_ordering_[next_depth];

            for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                uint32_t temp_idx = local_candidates_[next_depth][j];
                si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
            }
            for(ui i = 0; i < catalog_->sampling_num; ++i){
                this -> prob = 1;
                bool flag = rwch_non_core_single(next_depth,max_depth); //only one-time sampling
            }
            for(auto iter = res_tmp.begin(); iter != res_tmp.end(); ++iter){
                auto iter_res = this->catalog_->res_map.find(iter->first);
                if(iter_res == this->catalog_->res_map.end()) this->catalog_->res_map[iter->first] = iter->second / catalog_ -> sampling_num;
                else iter_res->second += iter->second / catalog_ -> sampling_num;
            }
            res_tmp.clear();
        }
        else if(num_local_candidates_[cur_depth] > 0 && start_depth < num_core_vertex_ - 2){
            for(ui i = 0; i < catalog_->sampling_num; ++i){
                bool valid = true;
                this -> prob = 1;
                cur_depth = start_depth + 1;
                uint32_t ext_core_depth = cur_depth;

                while(valid){
                    while (cur_depth <= max_depth - 2){
                        if(num_local_candidates_[cur_depth] == 0) {
                            valid = false;
                            break;
                        }
                        uint32_t u = vertex_ordering_[cur_depth];
                        std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                        int random = d(gen);

                        uint32_t encoded_id = local_candidates_[cur_depth][random];
                        uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                        if(visited_[v]){
                            valid = false;
                            break;
                        }
                        visited_[v] = true;
                        prob = prob / num_local_candidates_[cur_depth];
                        
                        idx_embedding_[cur_depth] = encoded_id;
                        embedding_[cur_depth] = v;
                        ext_core_depth++;

                        //compute ptb candidates
                        if(end_node_bn.test(cur_depth)){
                            if(catalog_->end_depth < num_core_vertex_) num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates(cur_depth);
                            else num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates_nc(cur_depth);
                        }
                        else num_ptb_candidates_[cur_depth] = 0;

                        if(cur_depth + 1 == this->catalog_->end_depth) {
                            uint32_t last_cnt = 0;
                            uint32_t last_depth = 0;
                            for(int k = cur_depth; k >= 0; k--){
                                if(end_node_bn.test(k)){
                                    last_depth = k;
                                    last_cnt = num_ptb_candidates_[k];
                                    break;
                                }
                            }
                            for(int k = 0; k < last_cnt;++k){
                                local_candidates_[cur_depth + 1][k] = ptb_candidates_[last_depth][k];
                            }
                            num_local_candidates_[cur_depth + 1] = last_cnt;
                        }
                        else num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                        //early termination
                        if(end_node_bn.test(cur_depth) && num_ptb_candidates_[cur_depth] == 0){
                            if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                            num_local_candidates_[cur_depth + 1] = 0;
                        }

                        if(cur_depth == max_depth - 2){
                            uint32_t temp_count = num_local_candidates_[cur_depth + 1];
                            uint32_t* temp_buffer = local_candidates_[cur_depth + 1];
                            uint32_t uu = vertex_ordering_[cur_depth + 1];
                            if(temp_count > 0) {
                                std::uniform_int_distribution<> d(0,temp_count - 1);
                                int random = d(gen);
                                uint32_t temp_v = catalog_->candidate_sets_[uu][temp_buffer[random]];

                                if(!visited_[temp_v]){

                                    prob = prob / temp_count;

                                    embedding_[cur_depth + 1] = temp_v;
                                    pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};

                                    auto iter = this->res_tmp.find(ptb);
                                    if(iter == this->res_tmp.end()) this->res_tmp[ptb] = (1 / prob);
                                    else iter->second += (1 / prob);

                                    auto iter_success = this->catalog_->success_map.find(ptb);
                                    if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                                    else iter_success->second++;
                                }
                            }
                            valid = false;
                            break;
                        }
                        else if(cur_depth == num_core_vertex_ - 2){
                            uint32_t next_depth = cur_depth + 1;
                            uint32_t next_vertex = vertex_ordering_[next_depth];
                            if(num_local_candidates_[next_depth] == 0) {
                                valid = false;
                                break;
                            }

                            for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                                uint32_t temp_idx = local_candidates_[next_depth][j];
                                si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                            }

                            bool flag = rwch_non_core_single(next_depth,max_depth); //only one-time sampling

                            valid = false;
                            break;

                        }
                        else{
                            cur_depth += 1;
                            if(cur_depth > max_depth - 2) valid = false;
                        }
                    
                    }
                }
                for(int k = start_depth + 1; k < ext_core_depth;++k){
                    visited_[embedding_[k]] = 0;
                }

                for(int k = start_depth + 1; k < ext_core_depth;++k){
                    embedding_[k] = 0;
                }
            }
            for(auto iter = res_tmp.begin(); iter != res_tmp.end(); ++iter){
                auto iter_res = this->catalog_->res_map.find(iter->first);
                if(iter_res == this->catalog_->res_map.end()) this->catalog_->res_map[iter->first] = iter->second / catalog_->sampling_num;
                else iter_res->second += iter->second / catalog_->sampling_num;
            }
            res_tmp.clear();
        }

        visited_[embedding_[start_depth]] = false;
        embedding_[start_depth] = 0;
    }
    return false;
}

uint64_t leapfrogtriejoin::execute_rwch(){

    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = num_core_vertex_;
    uint32_t enum_depth = catalog_-> max_depth;
    bool core_sp = (catalog_->max_depth < core_vertex_depth);

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    intermediate_result_count_[0] = input_->get_size();

    for(int i = 0; i < catalog_->end_depth;++i){
        if(end_node_bn.test(i)){
            if(i == catalog_->begin_depth) begin_node_first = true;
            else begin_node_first = false;
            break;
        }
    }

    for (uint64_t i = 0; i < input_->get_size(); ++i) {

        memcpy(embedding_, input_->get_tuple(i), sizeof(uint32_t) * start_depth);

        set_idx(start_depth);

        uint32_t check_v = embedding_[0];

        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }

        //init look-ahead candidates
        if(end_node_bn.test(start_depth - 1)){
            if(catalog_->end_depth < core_vertex_depth) num_ptb_candidates_[start_depth - 1] = pre_compute_ptb_candidates(start_depth - 1);
            else num_ptb_candidates_[start_depth - 1] = pre_compute_ptb_candidates_nc(start_depth - 1);
        }
        else{
            num_ptb_candidates_[start_depth - 1] = 0;
        }
    
        if(end_node_bn.test(start_depth - 1) && num_ptb_candidates_[start_depth - 1] == 0 && (start_depth - 1 != catalog_->begin_depth)){
            visited_[embedding_[0]] = false;
            continue;
        }

        if (start_depth >= num_core_vertex_) {
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }

            bool exit = rwch_non_core(start_depth, max_depth);

            if (g_exit)
                return count_;
        }
        else{
            uint32_t cur_depth = start_depth;
            idx_[cur_depth] = 0;
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    uint32_t u = vertex_ordering_[cur_depth];

                    uint32_t encoded_id = local_candidates_[cur_depth][idx_[cur_depth]++];

                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        continue;
                    }

                    visited_[v] = true;
                    embedding_[cur_depth] = v;

                    //compute ptb candidates
                    if(end_node_bn.test(cur_depth)){
                        if(catalog_->end_depth < core_vertex_depth) num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates(cur_depth);
                        else num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates_nc(cur_depth);
                    }
                    else num_ptb_candidates_[cur_depth] = 0;

                    if(cur_depth + 1 == this->catalog_->end_depth) {
                        uint32_t last_cnt = 0;
                        uint32_t last_depth = 0;
                        for(int k = cur_depth; k >= 0; k--){
                            if(end_node_bn.test(k)){
                                last_depth = k;
                                last_cnt = num_ptb_candidates_[k];
                                break;
                            }
                        }
                        for(int k = 0; k < last_cnt;++k){
                            local_candidates_[cur_depth + 1][k] = ptb_candidates_[last_depth][k];
                        }
                        num_local_candidates_[cur_depth + 1] = last_cnt;
                    }
                    else num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                    //early termination
                    if(end_node_bn.test(cur_depth) && num_ptb_candidates_[cur_depth] == 0){
                        if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                           num_local_candidates_[cur_depth + 1] = 0;
                    }

                    if (cur_depth == max_depth - 2) {
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        uint32_t uu = vertex_ordering_[cur_depth + 1];
                        for (uint32_t j = 0; j < loop_count; ++j) {
                            uint32_t temp_v = catalog_->candidate_sets_[uu][temp_buffer[j]];
                            if (visited_[temp_v]) {
                                continue;
                            }
                            if (count_ < std::numeric_limits<uint64_t>::max()) {
                                embedding_[cur_depth + 1] = temp_v;

                                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                                auto iter = this->catalog_->res_map.find(ptb);
                                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                                else iter->second++;

                                auto iter_success = this->catalog_->success_map.find(ptb);
                                if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                                else iter_success->second++;
                            }
                        }
                        visited_[v] = false;
                    }
                    else if(cur_depth == enum_depth - 2 && core_sp){
                        uint32_t next_depth = cur_depth + 1;
                        //core-sampling
                        if(num_local_candidates_[next_depth] > 0) {
                            bool flag = rwch_core(next_depth, max_depth);
                        }
                        visited_[v] = false;
                    }
                    else if(cur_depth == core_vertex_depth - 2 && !core_sp){
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];

                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }
                        bool exit = rwch_non_core(next_depth, max_depth);
                        visited_[v] = false;
                    }
                    else{
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
EXIT_TO_MAIN_PROCESS:
                    if (g_exit)
                        return count_;   
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else{
                    visited_[embedding_[cur_depth]] = false;
                }
            }
        }
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = false;
        }
    }

    return count_;
}

bool leapfrogtriejoin::enumerate_non_core_results_rwc(uint32_t start_depth, uint32_t max_depth){
    std::random_device rd;
    std:mt19937 gen_nc(rd());

    std::uniform_int_distribution<> d(0,num_local_candidates_[start_depth] - 1);
    int random = d(gen_nc);
    uint32_t exp_depth = start_depth;

    uint32_t v = si_buffer_[random];
    if(visited_[v]){
        return false;
    }

    visited_[v] = true;
    embedding_[start_depth] = v;
    exp_depth++;

    prob = prob / (double)num_local_candidates_[start_depth];

    uint32_t cur_depth = start_depth + 1;
    uint32_t u = vertex_ordering_[cur_depth];
    uint32_t bn = get_backward_neighbors(cur_depth)[0];
    uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
    uint32_t key = embedding_[bns_depth];
    if(catalog_->end_depth == cur_depth) num_local_candidates_[cur_depth] =  compute_local_candidates_noncore_ptb(bn,u,key,catalog_->end_depth);
    else local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);

    intermediate_result_count_[cur_depth] += num_local_candidates_[cur_depth];
    bool valid_nc = true;

    if (cur_depth == max_depth - 1) {
        uint32_t temp_count = num_local_candidates_[cur_depth];
        uint32_t* temp_buffer = local_candidates_[cur_depth];
        if(temp_count == 0) {
            valid_nc = false;
        }
        else{
            prob = prob / (double)temp_count;

            std::uniform_int_distribution<> d(0,temp_count - 1);
            int random = d(gen_nc);

            uint32_t temp_v = temp_buffer[random];
            if(!visited_[temp_v]){
                embedding_[cur_depth] = temp_v;

                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                auto iter = this->catalog_->res_map.find(ptb);
                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = (uint64_t)(1 / prob);
                else iter->second += (uint64_t)(1 / prob);

                auto iter_success = this->catalog_->success_map.find(ptb);
                if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                else iter_success->second += 1;
               
            }
        }
    }
    else{
        while (valid_nc)
        {
            while (cur_depth <= max_depth - 2){
                if(num_local_candidates_[cur_depth] == 0) {
                    valid_nc = false;
                    break;
                }
                prob = prob / (double)num_local_candidates_[cur_depth];

                std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                int random = d(gen_nc);

                uint32_t v = local_candidates_[cur_depth][random];
                
                if (visited_[v]) {
                    valid_nc = false;
                    break;
                }

                visited_[v] = true;
                embedding_[cur_depth] = v;
                exp_depth++;

                uint32_t next_depth = cur_depth + 1;
                u = vertex_ordering_[next_depth];
                bn = get_backward_neighbors(next_depth)[0];
                bns_depth = get_backward_neighbors_depth(next_depth)[0]; 
                key = embedding_[bns_depth];
                if(this->catalog_->end_depth == next_depth)num_local_candidates_[cur_depth + 1] = compute_local_candidates_noncore_ptb(bn,u,key,catalog_->end_depth);
                else local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);

                if(cur_depth == max_depth - 2){
                    uint32_t temp_count = num_local_candidates_[next_depth];
                    uint32_t *temp_buffer = local_candidates_[next_depth];    
                    if(temp_count == 0) {
                        valid_nc = false;
                        break;
                    }
                    else{
                        prob = prob / (double)temp_count;

                        std::uniform_int_distribution<> d(0,temp_count - 1);
                        int random = d(gen_nc);

                        uint32_t temp_v = temp_buffer[random];
                        if(!visited_[temp_v]){
                            embedding_[next_depth] = temp_v;
                            pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                            auto iter = this->catalog_->res_map.find(ptb);
                            if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = (uint64_t)(1 / prob);
                            else iter->second += (uint64_t)(1 / prob);

                            auto iter_success = this->catalog_->success_map.find(ptb);
                            if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                            else iter_success->second += 1;
                        }
                        
                        valid_nc = false;
                        break;
                    }
                }
                else{
                    cur_depth += 1;
                }
            }
        }   
    }
    
    //reset the "visited"
    for (uint32_t j = start_depth; j < exp_depth; ++j) {
        visited_[embedding_[j]] = false;
    }

    //reset the "embedding"
    for(uint32_t j = start_depth; j < exp_depth; ++j) {
        embedding_[j] = 0;
    }    

    return false;    
}

uint64_t leapfrogtriejoin::execute_rwc(){
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = num_core_vertex_;

    intermediate_result_count_[0] = input_->get_size();

    uint32_t init_size = input_->get_size();

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    
    std::random_device rd;
    std:mt19937 gen(rd());
    std::uniform_int_distribution<> d(0,init_size - 1);

    while(true){
        bool valid = true;
        int random = d(gen);
        memcpy(embedding_, input_->get_tuple(random), sizeof(uint32_t) * start_depth);

        set_idx(start_depth);

        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }

        prob = 1 / (double)init_size;
        uint32_t cur_depth = start_depth;
        uint32_t ext_depth = start_depth;

        if(start_depth >= num_core_vertex_){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }

            if(num_local_candidates_[start_depth] > 0) {
                bool flag = enumerate_non_core_results_rwc(start_depth,max_depth);
            }
        }
        else{
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

            while(valid){
                while (cur_depth <= max_depth - 2){
                    if(num_local_candidates_[cur_depth] == 0) {
                        valid = false;
                        break;
                    }
                    prob = prob / (double)num_local_candidates_[cur_depth];
                    

                    uint32_t u = vertex_ordering_[cur_depth];
                    std::uniform_int_distribution<> d(0,num_local_candidates_[cur_depth] - 1);
                    int random = d(gen);

                    uint32_t encoded_id = local_candidates_[cur_depth][random];
                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        valid = false;
                        break;
                    }
                    visited_[v] = true;
                    ext_depth++;
                    embedding_[cur_depth] = v;

                    //end_depth
                    if(cur_depth + 1 == this->catalog_->end_depth)  num_local_candidates_[cur_depth + 1] = compute_local_candidates_core_ptb(cur_depth + 1);
                    else num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                    if(cur_depth == max_depth - 2){
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        if(loop_count == 0) {
                            valid = false;
                            break;
                        }
                        prob = prob / (double)num_local_candidates_[cur_depth + 1];

                        std::uniform_int_distribution<> d(0,loop_count - 1);
                        int random = d(gen);

                        uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[random]];
                        if(visited_[temp_v]){
                            valid = false;
                            break;
                        }
                        embedding_[cur_depth + 1] = temp_v;
                        ext_depth++;

                        pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                        auto iter = this->catalog_->res_map.find(ptb);
                        if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = (uint64_t)(1 / prob);
                        else iter->second += (uint64_t)(1 / prob);

                        auto iter_success = this->catalog_->success_map.find(ptb);
                        if(iter_success == this->catalog_->success_map.end()) this->catalog_->success_map[ptb] = 1;
                        else iter_success->second += 1;

                        valid = false;
                        break;
                    }
                    else if(cur_depth == core_vertex_depth - 2){
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];
                        if(num_local_candidates_[next_depth] == 0) {
                            valid = false;
                            break;
                        }

                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }

                        bool exit = enumerate_non_core_results_rwc(next_depth, max_depth);
                        if(exit || !exit) {
                            valid = false;
                            break;
                        }
                    }
                    else{
                        cur_depth += 1;
                        if(cur_depth > max_depth - 2) valid = false;
                    }
                }
            }
        }

        for(int k = 0; k < ext_depth;++k){
            visited_[embedding_[k]] = false;
        }

        for(int k = 0; k < ext_depth;++k){
            embedding_[k] = 0;
        }

        catalog_->total_sp_num++;
        
        if (g_exit)
            return count_;

    }

    return count_;
}

bool leapfrogtriejoin::enumerate_non_core_results_pe(uint32_t start_depth, uint32_t max_depth){
    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {

        uint32_t v = si_buffer_[idx_[start_depth]];
        if (visited_[v]) {
            continue;
        }
        visited_[v] = true;

        embedding_[start_depth] = v;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);

        if (cur_depth == max_depth - 1) {
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                if (visited_[temp_buffer[j]]) {
                    continue;
                }

                embedding_[cur_depth] = temp_buffer[j];

                probe();
            }
            visited_[v] = false;
        }
        else{
            idx_[cur_depth] = 0;
            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    v = local_candidates_[cur_depth][idx_[cur_depth]++];

                    if (visited_[v]) {
                        continue;
                    }

                    visited_[v] = true;

                    embedding_[cur_depth] = v;
                    uint32_t next_depth = cur_depth + 1;
                    u = vertex_ordering_[next_depth];

                    bn = get_backward_neighbors(next_depth)[0];
                    bns_depth = get_backward_neighbors_depth(next_depth)[0];

                    key = embedding_[bns_depth];
                    local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);

                    if (cur_depth == max_depth - 2) {
                        uint32_t temp_count = num_local_candidates_[next_depth];
                        uint32_t *temp_buffer = local_candidates_[next_depth];
                        for (uint32_t j = 0; j < temp_count; ++j) {

                            if (visited_[temp_buffer[j]]) {
                                continue;
                            }
                            embedding_[next_depth] = temp_buffer[j];

                            probe();
                        }
                        visited_[v] = false;
                    }
                    else{
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
                    if (g_exit)
                        return true;
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else{
                    visited_[embedding_[cur_depth]] = false;
                }
            }
            visited_[embedding_[start_depth]] = false;
        }
    }

    return false;
}

bool leapfrogtriejoin::enumerate_non_core_results_probe_mat(uint32_t start_depth, uint32_t max_depth){
    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {

        uint32_t v = si_buffer_[idx_[start_depth]];
        if (visited_[v]) {
            continue;
        }
        visited_[v] = true;

        embedding_[start_depth] = v;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);

        if (cur_depth == max_depth - 1) {
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                if (visited_[temp_buffer[j]]) {
                    continue;
                }

                embedding_[cur_depth] = temp_buffer[j];

                for(auto depth:catalog_->cld_pe){
                    catalog_->probe_instances.emplace_back(embedding_[depth]);
                }
                catalog_->probe_cnt++;
            }
            visited_[v] = false;
        }
        else{
            idx_[cur_depth] = 0;
            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    v = local_candidates_[cur_depth][idx_[cur_depth]++];

                    if (visited_[v]) {
                        continue;
                    }

                    visited_[v] = true;

                    embedding_[cur_depth] = v;
                    uint32_t next_depth = cur_depth + 1;
                    u = vertex_ordering_[next_depth];

                    bn = get_backward_neighbors(next_depth)[0];
                    bns_depth = get_backward_neighbors_depth(next_depth)[0];

                    key = embedding_[bns_depth];
                    local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);

                    if (cur_depth == max_depth - 2) {
                        uint32_t temp_count = num_local_candidates_[next_depth];
                        uint32_t *temp_buffer = local_candidates_[next_depth];
                        for (uint32_t j = 0; j < temp_count; ++j) {

                            if (visited_[temp_buffer[j]]) {
                                continue;
                            }
                            embedding_[next_depth] = temp_buffer[j];

                            for(auto depth:catalog_->cld_pe){
                                catalog_->probe_instances.emplace_back(embedding_[depth]);
                            }
                            catalog_->probe_cnt++;

                        }
                        visited_[v] = false;
                    }
                    else{
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
                    if (g_exit)
                        return true;
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else{
                    visited_[embedding_[cur_depth]] = false;
                }
            }
            visited_[embedding_[start_depth]] = false;
        }
    }

    return false;
}

uint64_t leapfrogtriejoin::execute_mat_probe(){
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = num_core_vertex_;

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    intermediate_result_count_[0] = input_->get_size();

    for (uint64_t i = 0; i < input_->get_size(); ++i) {
        memcpy(embedding_, input_->get_tuple(i), sizeof(uint32_t) * start_depth);
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }
        if(start_depth == max_depth - 1){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            uint32_t *temp_buffer = local_candidates_[start_depth];
            uint32_t loop_count = num_local_candidates_[start_depth];

            for (uint32_t j = 0; j < loop_count; ++j) {
                uint32_t temp_v = local_candidates_[start_depth][j];
                embedding_[start_depth] = temp_v;

                for(auto depth:catalog_->cld_pe){
                    catalog_->probe_instances.emplace_back(embedding_[depth]);
                }
                catalog_->probe_cnt++;
            }
        }
        else if(start_depth >= num_core_vertex_){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }

            bool exit = enumerate_non_core_results_probe_mat(start_depth, max_depth);

            if (exit) g_exit = true;

            if (g_exit)
                return count_;
        }
        else{
            uint32_t cur_depth = start_depth;
            idx_[cur_depth] = 0;
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    uint32_t u = vertex_ordering_[cur_depth];
                    uint32_t encoded_id = local_candidates_[cur_depth][idx_[cur_depth]++];

                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];
                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        continue;
                    }
                    visited_[v] = true;

                    embedding_[cur_depth] = v;

                    num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                    if (cur_depth == max_depth - 2) {
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        for (uint32_t j = 0; j < loop_count; ++j) {
                            uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[j]];
                            if (visited_[temp_v]) {
                                continue;
                            }
                            embedding_[cur_depth + 1] = temp_v;

                            for(auto depth:catalog_->cld_pe){
                                catalog_->probe_instances.emplace_back(embedding_[depth]);
                            }
                            catalog_->probe_cnt++;
                        }
                        visited_[v] = false;
                    }
                    else if(cur_depth == core_vertex_depth - 2){
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];
                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }
                        bool exit = enumerate_non_core_results_probe_mat(cur_depth, max_depth);
                        visited_[v] = false;

                        if (exit) g_exit = true;

                        if (g_exit)
                            return count_;
                    }
                    else{
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
EXIT_TO_MAIN_PROCESS:
                    if (g_exit)
                        return count_;
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else {
                    visited_[embedding_[cur_depth]] = false;
                }
            }
        }
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = false;
        }
    }
    return count_;
}

void leapfrogtriejoin::probe(){
    //probe
    uint32_t head = embedding_[catalog_->end_depth];

    for (uint64_t k = 0; k < catalog_->probe_cnt; ++k) {
        uint64_t start_current = k*catalog_->probe_depth;
        uint64_t end_current = (k+1)*catalog_->probe_depth;

        uint32_t tail = catalog_->probe_instances[end_current - 1];
        if(!catalog_->data_graph_->checkEdgeExistence(head,tail)){
            bool valid = true;
            for(ui m = 0; m < num_vertex_;++m){
                auto iter = std::find(catalog_->probe_instances.begin() + start_current,catalog_->probe_instances.begin() + end_current,embedding_[m]);
                if(iter != catalog_->probe_instances.begin() + end_current) {
                    valid = false;
                    break;
                }
            }
            if(valid) {
                pEdge ptb = {head,tail};
                auto iter = this->catalog_->res_map.find(ptb);
                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                else iter->second++;                            
            }
        }
    }
}

bool leapfrogtriejoin::enumerate_non_core_results_iee(uint32_t start_depth, uint32_t max_depth){
    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {

        uint32_t v = si_buffer_[idx_[start_depth]];

        if (visited_[v]) {
            continue;
        }

        visited_[v] = true;

        embedding_[start_depth] = v;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);
        
        if (cur_depth == max_depth - 1) {
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                uint32_t temp_v = temp_buffer[j];
                if (visited_[temp_buffer[j]]) {
                    continue;
                }
                embedding_[cur_depth] = temp_v;
                pEdge ptb = {embedding_[0],embedding_[1]};
                auto iter = this->catalog_->res_map.find(ptb);
                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                else iter->second++;
            }
            visited_[v] = false;
        }
        else{
            idx_[cur_depth] = 0;
            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    v = local_candidates_[cur_depth][idx_[cur_depth]++];
                    if (visited_[v]) {
                        continue;
                    }
                    visited_[v] = true;
                    embedding_[cur_depth] = v;
                    uint32_t next_depth = cur_depth + 1;
                    u = vertex_ordering_[next_depth];

                    bn = get_backward_neighbors(next_depth)[0];
                    bns_depth = get_backward_neighbors_depth(next_depth)[0];

                    key = embedding_[bns_depth];
                    local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);
                    
                    if (cur_depth == max_depth - 2) {
                         uint32_t temp_count = num_local_candidates_[next_depth];
                        uint32_t *temp_buffer = local_candidates_[next_depth];
                        for (uint32_t j = 0; j < temp_count; ++j) {
                            uint32_t temp_v = temp_buffer[j];
                            if (visited_[temp_v]) {
                                continue;
                            }
                            pEdge ptb = {embedding_[0],embedding_[1]};
                            auto iter = this->catalog_->res_map.find(ptb);
                            if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                            else iter->second++;  
                        }
                        visited_[v] = false;
                    }
                    else{
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
                    if (g_exit)
                        return true;
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else{
                    visited_[embedding_[cur_depth]] = false;
                }
            }
            visited_[embedding_[start_depth]] = false;
        }
    }
    return false;
}

uint64_t leapfrogtriejoin::execute_iee(){
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = catalog_->combine_core;

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    intermediate_result_count_[0] = input_->get_size();
    uint32_t start_depth = 2;

    for (auto &pair:catalog_->candidate_pairs) {
        embedding_[0] = pair.first;
        embedding_[1] = pair.second;

        for (uint32_t j = 0; j < start_depth; ++j) {
            uint32_t u = vertex_ordering_[j];
            visited_[embedding_[j]] = true;
            idx_embedding_[j] = catalog_->get_candidate_index(u, embedding_[j]);
        }

        if (start_depth >= num_core_vertex_) {
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);
            bool exit = enumerate_non_core_results_iee(start_depth, max_depth);
            if (exit) g_exit = true;

            if (g_exit)
                return count_;
        }
        else{
            uint32_t cur_depth = start_depth;
            idx_[cur_depth] = 0;
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    uint32_t u = vertex_ordering_[cur_depth];

                    uint32_t encoded_id = local_candidates_[cur_depth][idx_[cur_depth]++];

                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        continue;
                    }
                    visited_[v] = true;

                    embedding_[cur_depth] = v;

                    num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                    if (cur_depth == max_depth - 2) {
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        for (uint32_t j = 0; j < loop_count; ++j) {
                            uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[j]];
                            if (visited_[temp_v]) {
                                continue;
                            }
                            pEdge ptb = {embedding_[0],embedding_[1]};
                            auto iter = this->catalog_->res_map.find(ptb);
                            if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                            else iter->second++;
                        }
                        visited_[v] = false;
                    }
                    else if (cur_depth == core_vertex_depth - 2) {
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];
                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }

                        bool exit = enumerate_non_core_results_iee(next_depth, max_depth);
                        if (exit) g_exit = true;

                        visited_[v] = false;
                    }
                    else {
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
EXIT_TO_MAIN_PROCESS:
                    if (g_exit)
                        return count_;
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else{
                    visited_[embedding_[cur_depth]] = false;
                }
            }
        }
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = false;
        }
    }

    return count_;
}

uint64_t leapfrogtriejoin::execute_pe(){
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = num_core_vertex_;

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    
    for (uint64_t i = 0; i < input_->get_size(); ++i) {
        memcpy(embedding_, input_->get_tuple(i), sizeof(uint32_t) * start_depth);

        set_idx(start_depth);

        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }

        if(start_depth == max_depth - 1){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            uint32_t *temp_buffer = local_candidates_[start_depth];
            uint32_t loop_count = num_local_candidates_[start_depth];

            for (uint32_t j = 0; j < loop_count; ++j) {
                uint32_t temp_v = local_candidates_[start_depth][j];
                embedding_[start_depth] = temp_v;
                
                probe();
            }
        }
        else if(start_depth >= num_core_vertex_){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }

            bool exit = enumerate_non_core_results_pe(start_depth, max_depth);
            if (exit) g_exit = true;

            if (g_exit)
                return count_;
        }
        else{
            uint32_t cur_depth = start_depth;
            idx_[cur_depth] = 0;
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    uint32_t u = vertex_ordering_[cur_depth];
                    uint32_t encoded_id = local_candidates_[cur_depth][idx_[cur_depth]++];

                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];
                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        continue;
                    }
                    visited_[v] = true;

                    embedding_[cur_depth] = v;

                    num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);
                    if (cur_depth == max_depth - 2) {
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        for (uint32_t j = 0; j < loop_count; ++j) {
                            uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[j]];
                            if (visited_[temp_v]) {
                                continue;
                            }
                            embedding_[cur_depth + 1] = temp_v;

                            probe();
                        }
                        visited_[v] = false;
                    }
                    else if(cur_depth == core_vertex_depth - 2){
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];
                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }
                        bool exit = enumerate_non_core_results_pe(cur_depth, max_depth);
                        visited_[v] = false;

                        if (exit) g_exit = true;

                        if (g_exit)
                            return count_;
                    }
                    else{
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
EXIT_TO_MAIN_PROCESS:
                    if (g_exit)
                        return count_;
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else {
                    visited_[embedding_[cur_depth]] = false;
                }
            }
        }
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = false;
        }
    }
    return count_;
}

bool leapfrogtriejoin::enumerate_non_core_results_drq_single(uint32_t start_depth, uint32_t max_depth){
    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {

        uint32_t v = si_buffer_[idx_[start_depth]];
#ifndef HOMOMORPHISM
        if (visited_[v]) {
            continue;
        }
        visited_[v] = true;
#endif
        embedding_[start_depth] = v;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);

        if (cur_depth == max_depth - 1) {
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {

                if (visited_[temp_buffer[j]]) {
                    continue;
                }

                embedding_[cur_depth] = temp_buffer[j];

                uint32_t start_node = embedding_[catalog_->end_depth];


                auto iter_res = catalog_->res_map_single.find(start_node);
                if(iter_res == catalog_->res_map_single.end()) catalog_->res_map_single[start_node] = 1;
                else iter_res->second++;

                count_ += 1;

                if(!catalog_->common_label_vec.empty()){
                    auto iter = catalog_->common_label_map.find(start_node);
                    if(iter == catalog_->common_label_map.end()) {
                        catalog_->common_label_map[start_node].resize(catalog_->common_label_vec.size());
                        for(int k = 0; k < catalog_->common_label_vec.size();++k){
                            for(auto &depth: catalog_->cl_depth_vec[k]){
                                uint32_t cl_node = embedding_[depth];
                                auto iter_cl = catalog_->common_label_map[start_node][k].find(cl_node);
                                if(iter_cl == catalog_->common_label_map[start_node][k].end()) catalog_->common_label_map[start_node][k][cl_node] = 1;
                                else iter_cl->second++;
                            }
                        }
                    }
                    else{
                        for(int k = 0; k < catalog_->common_label_vec.size();++k){
                            for(auto &depth: catalog_->cl_depth_vec[k]){
                                uint32_t cl_node = embedding_[depth];
                                auto iter_cl = iter->second[k].find(cl_node);
                                if(iter_cl ==  iter->second[k].end())iter->second[k][cl_node] = 1;
                                else iter_cl->second++;
                            }
                        }
                    }
                }
            }
            visited_[v] = false;
        }
        else{
            idx_[cur_depth] = 0;
            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    v = local_candidates_[cur_depth][idx_[cur_depth]++];

                    if (visited_[v]) {
                        continue;
                    }

                    visited_[v] = true;

                    embedding_[cur_depth] = v;
                    uint32_t next_depth = cur_depth + 1;
                    u = vertex_ordering_[next_depth];

                    bn = get_backward_neighbors(next_depth)[0];
                    bns_depth = get_backward_neighbors_depth(next_depth)[0];

                    key = embedding_[bns_depth];
                    local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);

                    if (cur_depth == max_depth - 2) {
                        uint32_t temp_count = num_local_candidates_[next_depth];
                        uint32_t *temp_buffer = local_candidates_[next_depth];
                        for (uint32_t j = 0; j < temp_count; ++j) {

                            if (visited_[temp_buffer[j]]) {
                                continue;
                            }

                            embedding_[next_depth] = temp_buffer[j];

                            uint32_t start_node = embedding_[catalog_->end_depth];

                            auto iter_res = catalog_->res_map_single.find(start_node);
                            if(iter_res == catalog_->res_map_single.end()) catalog_->res_map_single[start_node] = 1;
                            else iter_res->second++;

                            count_ += 1;

                            if(!catalog_->common_label_vec.empty()){
                                auto iter = catalog_->common_label_map.find(start_node);
                                if(iter == catalog_->common_label_map.end()) {
                                    catalog_->common_label_map[start_node].resize(catalog_->common_label_vec.size());
                                    for(int k = 0; k < catalog_->common_label_vec.size();++k){
                                        for(auto &depth: catalog_->cl_depth_vec[k]){
                                            uint32_t cl_node = embedding_[depth];
                                            auto iter_cl = catalog_->common_label_map[start_node][k].find(cl_node);
                                            if(iter_cl == catalog_->common_label_map[start_node][k].end()) catalog_->common_label_map[start_node][k][cl_node] = 1;
                                            else iter_cl->second++;
                                        }
                                    }
                                }
                                else{
                                    for(int k = 0; k < catalog_->common_label_vec.size();++k){
                                        for(auto &depth: catalog_->cl_depth_vec[k]){
                                            uint32_t cl_node = embedding_[depth];
                                            auto iter_cl = iter->second[k].find(cl_node);
                                            if(iter_cl ==  iter->second[k].end())iter->second[k][cl_node] = 1;
                                            else iter_cl->second++;
                                        }
                                    }
                                }
                            }
                        }

                        visited_[v] = false;

                    }else {
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }

                    if (g_exit)
                        return true;
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else{
                    visited_[embedding_[cur_depth]] = false;
                }
            }
            visited_[embedding_[start_depth]] = false;
        }
    }

    return false;
}

uint64_t leapfrogtriejoin::execute(){
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = num_core_vertex_;

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    intermediate_result_count_[0] = input_->get_size();


    for (uint64_t i = 0; i < input_->get_size(); ++i) {

        memcpy(embedding_, input_->get_tuple(i), sizeof(uint32_t) * start_depth);

#ifdef COLLECT_INVALID_PR_STATISTICS
        enter_result_count_[start_depth - 1] = count_;
#endif

#if RELATION_STRUCTURE == 0
        set_idx(start_depth);
#endif
        uint32_t check_v = embedding_[0];
#ifndef HOMOMORPHISM
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }
#endif

        if(start_depth == max_depth - 1){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            uint32_t *temp_buffer = local_candidates_[start_depth];
            uint32_t loop_count = num_local_candidates_[start_depth];

            for (uint32_t j = 0; j < loop_count; ++j) {
                uint32_t temp_v = local_candidates_[start_depth][j];
                embedding_[start_depth] = temp_v;
                
                uint32_t start_node = embedding_[catalog_->end_depth];

                auto iter_res = catalog_->res_map_single.find(start_node);
                if(iter_res == catalog_->res_map_single.end()) catalog_->res_map_single[start_node] = 1;
                else iter_res->second++;

                count_ += 1;

                if(!catalog_->common_label_vec.empty()){
                    auto iter = catalog_->common_label_map.find(start_node);
                    if(iter == catalog_->common_label_map.end()) {
                        catalog_->common_label_map[start_node].resize(catalog_->common_label_vec.size());
                        for(int k = 0; k < catalog_->common_label_vec.size();++k){
                            for(auto &depth: catalog_->cl_depth_vec[k]){
                                uint32_t cl_node = embedding_[depth];
                                auto iter_cl = catalog_->common_label_map[start_node][k].find(cl_node);
                                if(iter_cl == catalog_->common_label_map[start_node][k].end()) catalog_->common_label_map[start_node][k][cl_node] = 1;
                                else iter_cl->second++;
                            }
                        }
                    }
                    else{
                        for(int k = 0; k < catalog_->common_label_vec.size();++k){
                            for(auto &depth: catalog_->cl_depth_vec[k]){
                                uint32_t cl_node = embedding_[depth];
                                auto iter_cl = iter->second[k].find(cl_node);
                                if(iter_cl ==  iter->second[k].end())iter->second[k][cl_node] = 1;
                                else iter_cl->second++;
                            }
                        }
                    }
                }
            }
        }
        else if(start_depth >= num_core_vertex_){
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }

            bool exit = enumerate_non_core_results_drq_single(start_depth, max_depth);

            if (exit) g_exit = true;

            if (g_exit)
                return count_;
        }
        else{
            uint32_t cur_depth = start_depth;
            idx_[cur_depth] = 0;
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    uint32_t u = vertex_ordering_[cur_depth];
                    uint32_t encoded_id = local_candidates_[cur_depth][idx_[cur_depth]++];

                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];
                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
                        continue;
                    }
                    visited_[v] = true;

                    embedding_[cur_depth] = v;

                    num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                    if (cur_depth == max_depth - 2) {
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        for (uint32_t j = 0; j < loop_count; ++j) {
                            uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[j]];
                            if (visited_[temp_v]) {
                                continue;
                            }
                            embedding_[cur_depth + 1] = temp_v;

                            uint32_t start_node = embedding_[catalog_->end_depth];

                            auto iter_res = catalog_->res_map_single.find(start_node);
                            if(iter_res == catalog_->res_map_single.end()) catalog_->res_map_single[start_node] = 1;
                            else iter_res->second++;

                            count_ += 1;

                            if(!catalog_->common_label_vec.empty()){
                                auto iter = catalog_->common_label_map.find(start_node);
                                if(iter == catalog_->common_label_map.end()) {
                                    catalog_->common_label_map[start_node].resize(catalog_->common_label_vec.size());
                                    for(int k = 0; k < catalog_->common_label_vec.size();++k){
                                        for(auto &depth: catalog_->cl_depth_vec[k]){
                                            uint32_t cl_node = embedding_[depth];
                                            auto iter_cl = catalog_->common_label_map[start_node][k].find(cl_node);
                                            if(iter_cl == catalog_->common_label_map[start_node][k].end()) catalog_->common_label_map[start_node][k][cl_node] = 1;
                                            else iter_cl->second++;
                                        }
                                    }
                                }
                                else{
                                    for(int k = 0; k < catalog_->common_label_vec.size();++k){
                                        for(auto &depth: catalog_->cl_depth_vec[k]){
                                            uint32_t cl_node = embedding_[depth];
                                            auto iter_cl = iter->second[k].find(cl_node);
                                            if(iter_cl ==  iter->second[k].end())iter->second[k][cl_node] = 1;
                                            else iter_cl->second++;
                                        }
                                    }
                                }
                            }
                        
                        }
                        visited_[v] = false;
                    }
                    else if(cur_depth == core_vertex_depth - 2){
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];
                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }
                        bool exit = enumerate_non_core_results_drq_single(cur_depth, max_depth);
                        visited_[v] = false;

                        if (exit) g_exit = true;

                        if (g_exit)
                            return count_;
                    }
                    else{
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
EXIT_TO_MAIN_PROCESS:
                    if (g_exit)
                        return count_;
                }
                cur_depth -= 1;
                if (cur_depth < start_depth) {
                    break;
                }else {
                    visited_[embedding_[cur_depth]] = false;
                }
            }
        }

#ifndef HOMOMORPHISM
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = false;
        }
#endif

    }
    return count_;
}

uint64_t leapfrogtriejoin::execute_je(){
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = num_core_vertex_;

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    intermediate_result_count_[0] = input_->get_size();

    for (uint64_t i = 0; i < input_->get_size(); ++i) {
        memcpy(embedding_, input_->get_tuple(i), sizeof(uint32_t) * start_depth);

#ifdef COLLECT_INVALID_PR_STATISTICS
        enter_result_count_[start_depth - 1] = count_;
#endif

#if RELATION_STRUCTURE == 0
        set_idx(start_depth);
#endif
        uint32_t check_v = embedding_[0];
#ifndef HOMOMORPHISM
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }
#endif

#ifdef FAILING_SET_PRUNING
        for (uint32_t j = 0; j < start_depth; ++j) {
            reverse_embedding_[embedding_[j]] = vertex_ordering_[j];
        }
#endif

        if (start_depth >= num_core_vertex_) {
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

#if RELATION_STRUCTURE == 0
            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }
#endif

#ifdef FAILING_SET_PRUNING
            if (num_local_candidates_[start_depth] == 0) {
                vec_failing_set_[start_depth - 1] = ancestors_[vertex_ordering_[start_depth]];
            } else {
                vec_failing_set_[start_depth - 1].reset();
            }
#endif

            bool exit = enumerate_non_core_results_je(start_depth, max_depth);
            if (exit) g_exit = true;
#ifdef COLLECT_INVALID_PR_STATISTICS
            if(enter_result_count_[start_depth - 1] == count_)
                invalid_leaf_pr_count_ += 1;
#endif

            if (g_exit)
                return count_;
        }
        else {
            uint32_t cur_depth = start_depth;
            idx_[cur_depth] = 0;
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

#ifdef FAILING_SET_PRUNING
            if (num_local_candidates_[cur_depth] == 0) {
                vec_failing_set_[cur_depth - 1] = ancestors_[vertex_ordering_[cur_depth]];
            } else {
                vec_failing_set_[cur_depth - 1].reset();
            }
#endif

            intermediate_result_count_[cur_depth] += num_local_candidates_[cur_depth];

            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    
                    //uint32_t check_num = num_local_candidates_[cur_depth];
                    uint32_t prev_v = embedding_[cur_depth - 1];
                    uint32_t num_can_curr = num_local_candidates_[cur_depth];
                    uint32_t u = vertex_ordering_[cur_depth];
#if RELATION_STRUCTURE == 0
                    uint32_t encoded_id = local_candidates_[cur_depth][idx_[cur_depth]++];

                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                    idx_embedding_[cur_depth] = encoded_id;
                    
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                    uint32_t v = local_candidates_[cur_depth][idx_[cur_depth]++];
#endif

#ifndef HOMOMORPHISM
                    if (visited_[v]) {
#ifdef COLLECT_FAIL_SI_STATISTICS
                        iso_conflict_count_ += 1;
#endif
#ifdef FAILING_SET_PRUNING
                        vec_failing_set_[cur_depth] = ancestors_[u];
                        vec_failing_set_[cur_depth] |= ancestors_[reverse_embedding_[v]];
                        vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
#endif
                        continue;
                    }
                       
                    visited_[v] = true;

#ifdef FAILING_SET_PRUNING
                    reverse_embedding_[v] = u;
#endif

#endif
                    embedding_[cur_depth] = v;

#ifdef COLLECT_INVALID_PR_STATISTICS
                    enter_result_count_[cur_depth] = count_;
#endif
                    if(cur_depth + 1 == this->catalog_->end_depth)  num_local_candidates_[cur_depth + 1] = compute_local_candidates_core_ptb(cur_depth + 1);
                    else num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

#ifdef FAILING_SET_PRUNING
                    if (num_local_candidates_[cur_depth + 1] == 0) {
                        vec_failing_set_[cur_depth] = ancestors_[vertex_ordering_[cur_depth + 1]];
                    } else {
                        vec_failing_set_[cur_depth].reset();
                    }
#endif

                    intermediate_result_count_[cur_depth + 1] += num_local_candidates_[cur_depth + 1];

                    if (cur_depth == max_depth - 2) {

#if defined(HOMOMORPHISM) && defined(OUTPUT_OPTIMIZATION) && RELATION_STRUCTURE == 0
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];

                        for (uint32_t x = 0; x < loop_count; ++x) {
                            last_level_candidates_[x] = temp_buffer[x];

#ifdef ENABLE_OUTPUT
                            if (count_ < OUTPUT_RESULT_NUM_LIMIT) {
                                embedding_[cur_depth + 1] = temp_buffer[x];
                                memcpy(current_position_, embedding_, num_vertex_ * sizeof(uint32_t));
                                current_position_ += num_vertex_;
                            }
                            else {
                                return true;
                            }
#endif

                            count_ += 1;
                        }
#else
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        uint32_t uu = vertex_ordering_[cur_depth + 1];
                        for (uint32_t j = 0; j < loop_count; ++j) {
#if RELATION_STRUCTURE == 0
                            //uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[j]];
                            uint32_t temp_v = catalog_->candidate_sets_[uu][temp_buffer[j]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                            uint32_t temp_v = temp_buffer[j];
#endif
                            
#ifndef HOMOMORPHISM
                            if (visited_[temp_v]) {
#ifdef COLLECT_FAIL_SI_STATISTICS
                                iso_conflict_count_ += 1;
#endif
#ifdef FAILING_SET_PRUNING
                                //vec_failing_set_[cur_depth + 1] = ancestors_[last_vertex];
                                vec_failing_set_[cur_depth + 1] = ancestors_[uu];
                                vec_failing_set_[cur_depth + 1] |= ancestors_[reverse_embedding_[temp_v]];
                                vec_failing_set_[cur_depth] |= vec_failing_set_[cur_depth + 1];
#endif
                                continue;
                            }
#endif

                            
#ifdef FAILING_SET_PRUNING
                            vec_failing_set_[cur_depth + 1].set();
                            vec_failing_set_[cur_depth] |= vec_failing_set_[cur_depth + 1];
#endif

#ifndef COUNT_RESULTS
                            last_level_candidates_[j] = temp_v;
#endif

#ifdef ENABLE_OUTPUT
                            if (count_ < std::numeric_limits<uint64_t>::max()) {
                                embedding_[cur_depth + 1] = temp_v;

                                //map stats instead of materialze
                                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                                auto iter = this->catalog_->res_map.find(ptb);
                                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                                else iter->second++;

                            }
                            else {
                                
                                return true;
                            }
#endif

                            count_ += 1;

                            
                        }
#endif
#ifdef FAILING_SET_PRUNING
                        if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                            vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                            idx_[cur_depth] = num_local_candidates_[cur_depth];
                        } else {
                            vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                        }

                        reverse_embedding_.erase(v);
#endif

#ifndef HOMOMORPHISM
                        visited_[v] = false;
#endif

#ifdef COLLECT_INVALID_PR_STATISTICS
                        if (count_ == enter_result_count_[cur_depth]) {
                            invalid_core_pr_count_ += 1;
                        }
#endif

                        if (count_ >= output_count_limit_) {
                            return count_;
                        }
                    }
                    else if (cur_depth == core_vertex_depth - 2) {
                        uint32_t next_depth = cur_depth + 1;
#if RELATION_STRUCTURE == 0
                        uint32_t next_vertex = vertex_ordering_[next_depth];

                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            // local_candidates_[next_depth][j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }
#endif

                        bool exit = enumerate_non_core_results_je(next_depth, max_depth);
                        if (exit) g_exit = true;
                        
#ifdef FAILING_SET_PRUNING
                        if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                            vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                            idx_[cur_depth] = num_local_candidates_[cur_depth];
                        } else {
                            vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                        }

                        reverse_embedding_.erase(v);
#endif

#ifdef COLLECT_INVALID_PR_STATISTICS
                        if (count_ == enter_result_count_[cur_depth]) {
                            invalid_core_pr_count_ += 1;
                        }
#endif

#ifndef HOMOMORPHISM
                        visited_[v] = false;
#endif
                    }
                    else {
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
EXIT_TO_MAIN_PROCESS:
                    if (g_exit)
                        return count_;
                }

                cur_depth -= 1;
#ifdef COLLECT_INVALID_PR_STATISTICS
                if (count_ == enter_result_count_[cur_depth]) {
                    invalid_core_pr_count_ += 1;
                }
#endif

                if (cur_depth < start_depth) {
                    break;
                } else {

#ifdef FAILING_SET_PRUNING
                    if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                        vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                        idx_[cur_depth] = num_local_candidates_[cur_depth];
                    } else {
                        vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                    }

                    reverse_embedding_.erase(embedding_[cur_depth]);
#endif

#ifndef HOMOMORPHISM
                    visited_[embedding_[cur_depth]] = false;
#endif
                }
            }
        }


#ifndef HOMOMORPHISM
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = false;
        }
#endif

#ifdef FAILING_SET_PRUNING
        reverse_embedding_.clear();
#endif

    }

    return count_;
}

uint32_t leapfrogtriejoin::pre_compute_ptb_candidates(uint32_t depth){
    uint32_t last_cnt = 0, final_cnt = 0;
    uint32_t last_depth = 0;
    for(int i = depth - 1; i >= 0; i--){
        if(end_node_bn.test(i)){
            last_depth = i;
            last_cnt = num_ptb_candidates_[i];
            break;
        }
    }

    if(last_cnt == 0) {
        uint32_t u = vertex_ordering_[catalog_->end_depth];
        
        uint32_t* bns = get_backward_neighbors(depth);
        uint32_t* bns_depth = get_backward_neighbors_depth(depth);

        uint32_t lc_count = 0;
        uint32_t bn = vertex_ordering_[depth];
        uint32_t key = idx_embedding_[depth];
        uint32_t *lc = catalog_->get_core_relation_children(bn, u, key, lc_count);
        final_cnt = lc_count;

        for(int k = 0; k < final_cnt;++k){
            ptb_candidates_[depth][k] = lc[k];
        }
    }
    else{
        if(last_depth == catalog_->begin_depth && begin_node_first){
            //set difference
            uint32_t u = vertex_ordering_[catalog_->end_depth];
            uint32_t bn = vertex_ordering_[depth];
            uint32_t key = idx_embedding_[depth];
            uint32_t lc_count = 0, cur_cnt = 0;
            uint32_t *lc = catalog_->get_core_relation_children(bn, u, key, lc_count);

            size_t i = 0, j = 0;
            while (i < lc_count && j < last_cnt) {
                if (lc[i] < ptb_candidates_[last_depth][j]) {
                    si_buffer_precomputed[cur_cnt++] = lc[i];
                    ++i;
                } else if (lc[i] > ptb_candidates_[last_depth][j]) {
                    ++j;
                } else {
                    ++i;
                    ++j;
                }
            }

            // Add remaining elements of lc
            while (i < lc_count) {
                si_buffer_precomputed[cur_cnt++] = lc[i];
                ++i;
            }

            //remove visited
            for(int k = 0; k < cur_cnt;++k){
                uint32_t encoded_id = si_buffer_precomputed[k];
                uint32_t v = catalog_->candidate_sets_[u][encoded_id];
                if(!visited_[v]){
                    ptb_candidates_[depth][final_cnt++] = si_buffer_precomputed[k];
                }
            }        
        }
        else if(depth != catalog_->begin_depth){
            //set intersection
            uint32_t u = vertex_ordering_[catalog_->end_depth];
            uint32_t* ptb_prev = ptb_candidates_[last_depth];
            uint32_t ptb_prev_count = num_ptb_candidates_[last_depth];

            uint32_t lc_count = 0, cur_cnt = 0;
            uint32_t bn = vertex_ordering_[depth];
            uint32_t key = idx_embedding_[depth];
            uint32_t *lc = catalog_->get_core_relation_children(bn, u, key, lc_count);

            ComputeSetIntersection::ComputeCandidates(ptb_prev, ptb_prev_count, lc, lc_count, si_buffer_precomputed, cur_cnt);

            for(int k = 0; k < cur_cnt;++k){
                uint32_t encoded_id = si_buffer_precomputed[k];
                uint32_t v = catalog_->candidate_sets_[u][encoded_id];
                if(!visited_[v]){
                    ptb_candidates_[depth][final_cnt++] = si_buffer_precomputed[k];
                }
            }
        }
        else if(depth == catalog_->begin_depth){
            //set difference
            uint32_t u = vertex_ordering_[catalog_->end_depth];
            uint32_t bn = vertex_ordering_[depth];
            uint32_t key = idx_embedding_[depth];
            uint32_t lc_count = 0, cur_cnt = 0;
            uint32_t *lc = catalog_->get_core_relation_children(bn, u, key, lc_count);

            size_t i = 0, j = 0;
            while (i < last_cnt && j < lc_count) {
                if (ptb_candidates_[last_depth][i] < lc[j]) {
                    si_buffer_precomputed[cur_cnt++] = ptb_candidates_[last_depth][i];
                    ++i;
                } else if (ptb_candidates_[last_depth][i] > lc[j]) {
                    ++j;
                } else {
                    ++i;
                    ++j;
                }
            }

            // Add remaining elements of lc1
            while (i < last_cnt) {
                si_buffer_precomputed[cur_cnt++] = ptb_candidates_[last_depth][i];
                ++i;
            }

            //remove visited
            for(int k = 0; k < cur_cnt;++k){
                uint32_t encoded_id = si_buffer_precomputed[k];
                uint32_t v = catalog_->candidate_sets_[u][encoded_id];
                if(!visited_[v]){
                    ptb_candidates_[depth][final_cnt++] = si_buffer_precomputed[k];
                }
            }
        }
    }
    
    return final_cnt;
}

uint32_t leapfrogtriejoin::pre_compute_ptb_candidates_nc(uint32_t depth){
    uint32_t last_cnt = 0, final_cnt = 0;
    uint32_t last_depth = 0;
    for(int i = depth - 1; i >= 0; i--){
        if(end_node_bn.test(i)){
            last_depth = i;
            last_cnt = num_ptb_candidates_[i];
            break;
        }
    }

    if(last_cnt == 0) {
        uint32_t u = vertex_ordering_[catalog_->end_depth];
        
        uint32_t* bns = get_backward_neighbors(depth);
        uint32_t* bns_depth = get_backward_neighbors_depth(depth);

        uint32_t lc_count = 0;
        uint32_t bn = vertex_ordering_[depth];
        uint32_t key = embedding_[depth];
        uint32_t *lc = catalog_->get_non_core_relation_children(bn, u, key, lc_count);
        final_cnt = lc_count;

        for(int k = 0; k < final_cnt;++k){
            ptb_candidates_[depth][k] = lc[k];
        }
    }
    else if(last_cnt > 0 && depth == catalog_->begin_depth){
        //set difference
        uint32_t u = vertex_ordering_[catalog_->end_depth];
        uint32_t bn = vertex_ordering_[depth];
        uint32_t key = embedding_[depth];
        uint32_t lc_count = 0,cur_cnt = 0;
        uint32_t *lc = catalog_->get_non_core_relation_children(bn, u, key, lc_count);

        size_t i = 0, j = 0;
        while (i < last_cnt && j < lc_count) {
            if (ptb_candidates_[last_depth][i] < lc[j]) {
                si_buffer_precomputed[cur_cnt++] = ptb_candidates_[last_depth][i];
                ++i;
            } else if (ptb_candidates_[last_depth][i] > lc[j]) {
                ++j;
            } else {
                ++i;
                ++j;
            }
        }

        // Add remaining elements of lc1
        while (i < last_cnt) {
            si_buffer_precomputed[cur_cnt++] = ptb_candidates_[last_depth][i];
            ++i;
        }

        //remove visited
        for(int k = 0; k < cur_cnt;++k){
            uint32_t v = si_buffer_precomputed[k];
            if(!visited_[v]){
                ptb_candidates_[depth][final_cnt++] = si_buffer_precomputed[k];
            }
        }
    }
    else if(last_cnt > 0 && last_depth == catalog_->begin_depth){

        //set difference
        uint32_t u = vertex_ordering_[catalog_->end_depth];
        uint32_t bn = vertex_ordering_[depth];
        uint32_t key = embedding_[depth];
        uint32_t lc_count = 0,cur_cnt = 0;
        uint32_t *lc = catalog_->get_non_core_relation_children(bn, u, key, lc_count);

        size_t i = 0, j = 0;
        while (i < lc_count && j < last_cnt) {
            if (lc[i] < ptb_candidates_[last_depth][j]) {
                si_buffer_precomputed[cur_cnt++] = lc[i];
                ++i;
            } else if (lc[i] > ptb_candidates_[last_depth][j]) {
                ++j;
            } else {
                ++i;
                ++j;
            }
        }

        // Add remaining elements of lc1
        while (i < lc_count) {
            si_buffer_precomputed[cur_cnt++] = lc[i];
            ++i;
        }

        //remove visited
        for(int k = 0; k < cur_cnt;++k){
            uint32_t v = si_buffer_precomputed[k];
            if(!visited_[v]){
                ptb_candidates_[depth][final_cnt++] = si_buffer_precomputed[k];
            }
        }
    }
    return final_cnt;
}

uint64_t leapfrogtriejoin::execute_lae() {
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = num_core_vertex_;

    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    intermediate_result_count_[0] = input_->get_size();

    for(int i = 0; i < catalog_->end_depth;++i){
        if(end_node_bn.test(i)){
            if(i == catalog_->begin_depth) begin_node_first = true;
            else begin_node_first = false;
            break;
        }
    }

    for (uint64_t i = 0; i < input_->get_size(); ++i) {
        memcpy(embedding_, input_->get_tuple(i), sizeof(uint32_t) * start_depth);

#if RELATION_STRUCTURE == 0
        set_idx(start_depth);
#endif
        uint32_t check_v = embedding_[0];
#ifndef HOMOMORPHISM
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }
#endif

        //init_ptb_candidates
        if(end_node_bn.test(start_depth - 1)){
            if(catalog_->end_depth < core_vertex_depth) num_ptb_candidates_[start_depth - 1] = pre_compute_ptb_candidates(start_depth - 1);
            else num_ptb_candidates_[start_depth - 1] = pre_compute_ptb_candidates_nc(start_depth - 1);
        }
        else{
            num_ptb_candidates_[start_depth - 1] = 0;
        }

        if(end_node_bn.test(start_depth - 1) && num_ptb_candidates_[start_depth - 1] == 0 && (start_depth - 1 != catalog_->begin_depth)){
            visited_[embedding_[0]] = false;
            continue;
        }

        if (start_depth >= num_core_vertex_) {
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);

#if RELATION_STRUCTURE == 0
            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }
#endif

            bool exit = enumerate_non_core_results_lae(start_depth, max_depth);
            if (exit) g_exit = true;

            if (g_exit)
                return count_;
        }
        else {
            uint32_t cur_depth = start_depth;
            idx_[cur_depth] = 0;
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);

            intermediate_result_count_[cur_depth] += num_local_candidates_[cur_depth];

            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    
                    uint32_t u = vertex_ordering_[cur_depth];
#if RELATION_STRUCTURE == 0
                    uint32_t encoded_id = local_candidates_[cur_depth][idx_[cur_depth]++];

                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];

                    idx_embedding_[cur_depth] = encoded_id;
#endif

#ifndef HOMOMORPHISM
                    if (visited_[v]) {
                        continue;
                    }                   
                    visited_[v] = true;
#endif
                    embedding_[cur_depth] = v;

                    //compute ptb candidates
                    if(end_node_bn.test(cur_depth)){
                        if(catalog_->end_depth < core_vertex_depth) num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates(cur_depth);
                        else num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates_nc(cur_depth);
                    }
                    else num_ptb_candidates_[cur_depth] = 0;
                    
                    if(cur_depth + 1 == this->catalog_->end_depth) {
                        uint32_t last_cnt = 0;
                        uint32_t last_depth = 0;
                        for(int k = cur_depth; k >= 0; k--){
                            if(end_node_bn.test(k)){
                                last_depth = k;
                                last_cnt = num_ptb_candidates_[k];
                                break;
                            }
                        }
                        for(int k = 0; k < last_cnt;++k){
                            local_candidates_[cur_depth + 1][k] = ptb_candidates_[last_depth][k];
                        }
                        num_local_candidates_[cur_depth + 1] = last_cnt;
                    }
                    else num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);

                    //early termination
                    if(end_node_bn.test(cur_depth) && num_ptb_candidates_[cur_depth] == 0){
                        if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                           num_local_candidates_[cur_depth + 1] = 0;
                    } 
                    intermediate_result_count_[cur_depth + 1] += num_local_candidates_[cur_depth + 1];

                    if (cur_depth == max_depth - 2) {

#if defined(HOMOMORPHISM) && defined(OUTPUT_OPTIMIZATION) && RELATION_STRUCTURE == 0
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];

                        for (uint32_t x = 0; x < loop_count; ++x) {
                            last_level_candidates_[x] = temp_buffer[x];

#ifdef ENABLE_OUTPUT
                            if (count_ < OUTPUT_RESULT_NUM_LIMIT) {
                                embedding_[cur_depth + 1] = temp_buffer[x];
                                memcpy(current_position_, embedding_, num_vertex_ * sizeof(uint32_t));
                                current_position_ += num_vertex_;
                            }
                            else {
                                return true;
                            }
#endif

                            count_ += 1;
                        }
#else
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        uint32_t uu = vertex_ordering_[cur_depth + 1];
                        for (uint32_t j = 0; j < loop_count; ++j) {

                            uint32_t temp_v = catalog_->candidate_sets_[uu][temp_buffer[j]];
                            
#ifndef HOMOMORPHISM
                            if (visited_[temp_v]) {
#ifdef COLLECT_FAIL_SI_STATISTICS
                                iso_conflict_count_ += 1;
#endif
                                continue;
                            }
#endif

#ifdef ENABLE_OUTPUT
                            if (count_ < std::numeric_limits<uint64_t>::max()) {
                                embedding_[cur_depth + 1] = temp_v;

                                //map stats instead of materialze
                                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                                auto iter = this->catalog_->res_map.find(ptb);
                                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                                else iter->second++;

                            }
                            else {
                                catalog_->exceed_limit = true;
                                return true;
                            }
#endif

                            count_ += 1;

                        }
#endif

#ifndef HOMOMORPHISM
                        visited_[v] = false;
#endif

                        if (count_ >= output_count_limit_) {
                            return count_;
                        }
                    }
                    else if (cur_depth == core_vertex_depth - 2) {
                        uint32_t next_depth = cur_depth + 1;
#if RELATION_STRUCTURE == 0
                        uint32_t next_vertex = vertex_ordering_[next_depth];

                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            // local_candidates_[next_depth][j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }
#endif

                        bool exit = enumerate_non_core_results_lae(next_depth, max_depth);
                        if (exit) g_exit = true;

#ifndef HOMOMORPHISM
                        visited_[v] = false;
#endif
                    }
                    else {
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
EXIT_TO_MAIN_PROCESS:
                    if (g_exit)
                        return count_;
                }

                cur_depth -= 1;

                if (cur_depth < start_depth) {
                    break;
                } else {

#ifndef HOMOMORPHISM
                    visited_[embedding_[cur_depth]] = false;
#endif
                }
            }
        }
#ifndef HOMOMORPHISM
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = false;
        }
#endif
    }

    return count_;
}

uint32_t leapfrogtriejoin::compute_local_candidates_noncore_ptb(uint32_t bn, uint32_t u, uint32_t key,uint32_t depth){

    uint32_t exist_count = 0,ref_count = 0,lc_count = 0;

    uint32_t* candidate  = catalog_->get_non_core_relation_children(bn, u, key, exist_count);

    this->lc_exist_count = exist_count;

    uint32_t ui = vertex_ordering_[catalog_->begin_depth];
    uint32_t uj = vertex_ordering_[catalog_->end_depth];

    uint32_t key_vi = embedding_[catalog_->begin_depth];
    uint32_t* ref = catalog_->get_non_core_relation_children(ui, uj, key_vi, ref_count);

    for(int i = 0; i < exist_count; ++i){
        if(!exist_check(candidate[i],ref,ref_count)){            
            local_candidates_[depth][lc_count] = candidate[i];
            lc_count++;
        }
    }
    return lc_count;
}

bool leapfrogtriejoin::enumerate_non_core_results_lae(uint32_t start_depth, uint32_t max_depth) {

    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {
#if RELATION_STRUCTURE == 0
        uint32_t v = si_buffer_[idx_[start_depth]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
        uint32_t v = local_candidates_[start_depth][idx_[start_depth]];
#endif

#ifndef HOMOMORPHISM
        if (visited_[v]) {
            continue;
        }
        visited_[v] = true;
#endif

        embedding_[start_depth] = v;

        //compute ptb candidates
        if(end_node_bn.test(start_depth)) num_ptb_candidates_[start_depth] = pre_compute_ptb_candidates_nc(start_depth);
        else num_ptb_candidates_[start_depth] = 0;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        if(cur_depth == catalog_->end_depth) {
            uint32_t last_cnt = 0;
            uint32_t last_depth = 0;
            for(int k = cur_depth - 1; k >= 0; k--){
                if(end_node_bn.test(k)){
                    last_depth = k;
                    last_cnt = num_ptb_candidates_[k];
                    break;
                }
            }
            for(int k = 0; k < last_cnt;++k){
                local_candidates_[cur_depth][k] = ptb_candidates_[last_depth][k];
            }
            num_local_candidates_[cur_depth] = last_cnt;
        }
        else local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);
        
        intermediate_result_count_[cur_depth] += num_local_candidates_[cur_depth];

        //early termination
        if(end_node_bn.test(start_depth) && num_ptb_candidates_[start_depth] == 0){
            if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                num_local_candidates_[cur_depth] = 0;
        }

        if (cur_depth == max_depth - 1) {
#if defined(HOMOMORPHISM) && defined(OUTPUT_OPTIMIZATION)

#ifdef ENABLE_OUTPUT
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                if (count_ < OUTPUT_RESULT_NUM_LIMIT) {
                    embedding_[cur_depth] = temp_buffer[j];
                    memcpy(current_position_, embedding_, num_vertex_ * sizeof(uint32_t));
                    current_position_ += num_vertex_;
                }
                else {
                    return true;
                }
                count_ += 1;
            }

#else
            count_ += num_local_candidates_[cur_depth];
#endif

#else
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                uint32_t temp_v = temp_buffer[j];
               
#ifndef HOMOMORPHISM
                if (visited_[temp_buffer[j]]) {
#ifdef COLLECT_FAIL_SI_STATISTICS
                    iso_conflict_count_ += 1;
#endif
                    continue;
                }
#endif

#ifdef ENABLE_OUTPUT
                if (count_ < std::numeric_limits<uint64_t>::max()) {

                    embedding_[cur_depth] = temp_v;
                    //map stats instead of materialze
                    pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                    auto iter = this->catalog_->res_map.find(ptb);
                    if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                    else iter->second++;
                }
                else {
                    catalog_->exceed_limit = true;
                    return true;
                }
#endif
                count_ += 1;
            }
#endif

#ifndef HOMOMORPHISM
            visited_[v] = false;
#endif
            // if (count_ >= output_count_limit_) {
            //     return true;
            // }
        }
        else {
            idx_[cur_depth] = 0;
            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    v = local_candidates_[cur_depth][idx_[cur_depth]++];
                    
#ifndef HOMOMORPHISM
                    if (visited_[v]) {
                        continue;
                    }
                    visited_[v] = true;
#endif
                    embedding_[cur_depth] = v;

                    if(end_node_bn.test(cur_depth)) num_ptb_candidates_[cur_depth] = pre_compute_ptb_candidates_nc(cur_depth);
                    else num_ptb_candidates_[cur_depth] = 0;

                    uint32_t next_depth = cur_depth + 1;
                    u = vertex_ordering_[next_depth];

                    bn = get_backward_neighbors(next_depth)[0];
                    bns_depth = get_backward_neighbors_depth(next_depth)[0];

                    key = embedding_[bns_depth];

                    if(next_depth == catalog_->end_depth){
                        uint32_t last_cnt = 0;
                        uint32_t last_depth = 0;
                        for(int k = cur_depth; k >= 0; k--){
                            if(end_node_bn.test(k)){
                                last_depth = k;
                                last_cnt = num_ptb_candidates_[k];
                                break;
                            }
                        }
                        for(int k = 0; k < last_cnt;++k){
                            local_candidates_[next_depth][k] = ptb_candidates_[last_depth][k];
                        }
                        num_local_candidates_[next_depth] = last_cnt;
                    } 
                    else local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);
                    
                    intermediate_result_count_[next_depth] += num_local_candidates_[next_depth];

                    //early termination
                    if(end_node_bn.test(cur_depth) && num_ptb_candidates_[cur_depth] == 0){
                        if(!this->begin_node_first || cur_depth != catalog_->begin_depth)
                           num_local_candidates_[next_depth] = 0;
                    }

                    if (cur_depth == max_depth - 2) {
                        uint32_t temp_count = num_local_candidates_[next_depth];
                        uint32_t *temp_buffer = local_candidates_[next_depth];
                        for (uint32_t j = 0; j < temp_count; ++j) {
                            uint32_t temp_v = temp_buffer[j];
                            
#ifndef HOMOMORPHISM
                            if (visited_[temp_buffer[j]]) {
                                continue;
                            }
#endif

#ifdef ENABLE_OUTPUT
                            if (count_ < std::numeric_limits<uint64_t>::max()) {
                                embedding_[cur_depth + 1] = temp_buffer[j];
                                //map stats instead of materialze
                                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                                auto iter = this->catalog_->res_map.find(ptb);
                                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                                else iter->second++;
                            }
                            else {
                                catalog_->exceed_limit = true;
                                return true;
                            }
#endif

                            count_ += 1;
                        }

#ifndef HOMOMORPHISM
                        visited_[v] = false;
#endif
                        if (count_ >= output_count_limit_) {
                            return true;
                        }
                    } else {
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
                    
                    if (g_exit)
                        return true;
                }

                cur_depth -= 1;

                if (cur_depth <= start_depth) {
                    break;
                } else {

#ifndef HOMOMORPHISM
                    visited_[embedding_[cur_depth]] = false;
#endif
                }

            }

#ifndef HOMOMORPHISM
            visited_[embedding_[start_depth]] = false;
#endif
        }
    }

    return false;
}

bool leapfrogtriejoin::enumerate_non_core_results_je(uint32_t start_depth, uint32_t max_depth) {

    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {
#if RELATION_STRUCTURE == 0
        uint32_t v = si_buffer_[idx_[start_depth]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
        uint32_t v = local_candidates_[start_depth][idx_[start_depth]];
#endif

#ifndef HOMOMORPHISM
        if (visited_[v]) {
#ifdef COLLECT_FAIL_SI_STATISTICS
            iso_conflict_count_ += 1;
#endif
#ifdef FAILING_SET_PRUNING
            vec_failing_set_[start_depth] = ancestors_[vertex_ordering_[start_depth]];
            vec_failing_set_[start_depth] |= ancestors_[reverse_embedding_[v]];
            vec_failing_set_[start_depth - 1] |= vec_failing_set_[start_depth];
#endif

            continue;
        }

#ifdef FAILING_SET_PRUNING
        reverse_embedding_[v] = vertex_ordering_[start_depth];
#endif
        visited_[v] = true;
#endif

#ifdef COLLECT_INVALID_PR_STATISTICS
        enter_result_count_[start_depth] = count_;
#endif

        embedding_[start_depth] = v;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        if(cur_depth == catalog_->end_depth) num_local_candidates_[cur_depth] = compute_local_candidates_noncore_ptb(bn,u,key,catalog_->end_depth);
        else local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);
        
        intermediate_result_count_[cur_depth] += num_local_candidates_[cur_depth];

#ifdef FAILING_SET_PRUNING
        if (num_local_candidates_[cur_depth] == 0) {
            vec_failing_set_[cur_depth - 1] = ancestors_[vertex_ordering_[cur_depth]];
        } else {
            vec_failing_set_[cur_depth - 1].reset();
        }
#endif


        if (cur_depth == max_depth - 1) {
#if defined(HOMOMORPHISM) && defined(OUTPUT_OPTIMIZATION)

#ifdef ENABLE_OUTPUT
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                if (count_ < OUTPUT_RESULT_NUM_LIMIT) {
                    embedding_[cur_depth] = temp_buffer[j];
                    memcpy(current_position_, embedding_, num_vertex_ * sizeof(uint32_t));
                    current_position_ += num_vertex_;
                }
                else {
                    return true;
                }
                count_ += 1;
            }

#else
            count_ += num_local_candidates_[cur_depth];
#endif

#else
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                uint32_t temp_v = temp_buffer[j];
               
#ifndef HOMOMORPHISM
                if (visited_[temp_buffer[j]]) {
#ifdef COLLECT_FAIL_SI_STATISTICS
                    iso_conflict_count_ += 1;
#endif
#ifdef FAILING_SET_PRUNING
                    vec_failing_set_[cur_depth] = ancestors_[vertex_ordering_[cur_depth]];
                    vec_failing_set_[cur_depth] |= ancestors_[reverse_embedding_[temp_buffer[j]]];
                    vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
#endif
                    continue;
                }
#endif

#ifdef FAILING_SET_PRUNING
                vec_failing_set_[cur_depth].set();
                vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
#endif

#ifndef COUNT_RESULTS
                last_level_candidates_[j] = temp_buffer[j];
#endif

#ifdef ENABLE_OUTPUT
                if (count_ < std::numeric_limits<uint64_t>::max()) {

                    embedding_[cur_depth] = temp_v;
                    //map stats instead of materialze
                    pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                    auto iter = this->catalog_->res_map.find(ptb);
                    if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                    else iter->second++;
                }
                else {
                    catalog_->exceed_limit = true;
                    return true;
                }
#endif
                count_ += 1;
            }
#endif

#ifdef FAILING_SET_PRUNING
            if (!vec_failing_set_[cur_depth - 1].test(vertex_ordering_[cur_depth - 1])) {
                vec_failing_set_[cur_depth - 2] = vec_failing_set_[cur_depth - 1];
                idx_[cur_depth - 1] = num_local_candidates_[cur_depth - 1];
            } else {
                vec_failing_set_[cur_depth - 2] |= vec_failing_set_[cur_depth - 1];
            }

            reverse_embedding_.erase(v);
#endif


#ifndef HOMOMORPHISM
            visited_[v] = false;
#endif

        }
        else {
            idx_[cur_depth] = 0;
            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    v = local_candidates_[cur_depth][idx_[cur_depth]++];
                    
#ifndef HOMOMORPHISM
                    if (visited_[v]) {
#ifdef COLLECT_FAIL_SI_STATISTICS
                        iso_conflict_count_ += 1;
#endif
#ifdef FAILING_SET_PRUNING
                        vec_failing_set_[cur_depth] = ancestors_[vertex_ordering_[cur_depth]];
                        vec_failing_set_[cur_depth] |= ancestors_[reverse_embedding_[v]];
                        vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
#endif
                        continue;
                    }

#ifdef FAILING_SET_PRUNING
                    reverse_embedding_[v] = vertex_ordering_[cur_depth];
#endif
                    visited_[v] = true;
#endif

#ifdef COLLECT_INVALID_PR_STATISTICS
                    enter_result_count_[cur_depth] = count_;
#endif

                    embedding_[cur_depth] = v;
                    uint32_t next_depth = cur_depth + 1;
                    u = vertex_ordering_[next_depth];

                    bn = get_backward_neighbors(next_depth)[0];
                    bns_depth = get_backward_neighbors_depth(next_depth)[0];

                    key = embedding_[bns_depth];

                    if(next_depth == catalog_->end_depth) num_local_candidates_[next_depth] = compute_local_candidates_noncore_ptb(bn,u,key,catalog_->end_depth);
                    else local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,num_local_candidates_[next_depth]);
                    
                    intermediate_result_count_[next_depth] += num_local_candidates_[next_depth];

#ifdef FAILING_SET_PRUNING
                    if (num_local_candidates_[next_depth] == 0) {
                        vec_failing_set_[next_depth - 1] = ancestors_[vertex_ordering_[next_depth]];
                    } else {
                        vec_failing_set_[next_depth - 1].reset();
                    }
#endif
                    if (cur_depth == max_depth - 2) {
                        uint32_t temp_count = num_local_candidates_[next_depth];
                        uint32_t *temp_buffer = local_candidates_[next_depth];
                        for (uint32_t j = 0; j < temp_count; ++j) {
                            uint32_t temp_v = temp_buffer[j];
                            
#ifndef HOMOMORPHISM
                            if (visited_[temp_buffer[j]]) {
#ifdef COLLECT_FAIL_SI_STATISTICS
                                iso_conflict_count_ += 1;
#endif
#ifdef FAILING_SET_PRUNING
                                vec_failing_set_[next_depth] = ancestors_[vertex_ordering_[next_depth]];
                                vec_failing_set_[next_depth] |= ancestors_[reverse_embedding_[temp_buffer[j]]];
                                vec_failing_set_[next_depth - 1] |= vec_failing_set_[next_depth];
#endif
                                continue;
                            }
#endif

#ifdef FAILING_SET_PRUNING
                            vec_failing_set_[next_depth].set();
                            vec_failing_set_[next_depth - 1] |= vec_failing_set_[next_depth];
#endif

#ifndef COUNT_RESULTS
                            last_level_candidates_[j] = temp_buffer[j];
#endif

#ifdef ENABLE_OUTPUT
                            if (count_ < std::numeric_limits<uint64_t>::max()) {
                                embedding_[cur_depth + 1] = temp_buffer[j];

                                //map stats instead of materialze
                                pEdge ptb = {embedding_[this->catalog_->begin_depth],embedding_[this->catalog_->end_depth]};
                                auto iter = this->catalog_->res_map.find(ptb);
                                if(iter == this->catalog_->res_map.end()) this->catalog_->res_map[ptb] = 1;
                                else iter->second++;
                            }
                            else {
                                catalog_->exceed_limit = true;
                                return true;
                            }
#endif

                            count_ += 1;
                        }

#ifdef FAILING_SET_PRUNING
                        if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                            vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                            idx_[cur_depth] = num_local_candidates_[cur_depth];
                        } else {
                            vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                        }

                        reverse_embedding_.erase(v);
#endif

#ifdef COLLECT_INVALID_PR_STATISTICS
                        if (enter_result_count_[cur_depth] == count_) {
                            invalid_leaf_pr_count_ += 1;
                        }
#endif

#ifndef HOMOMORPHISM
                        visited_[v] = false;
#endif
                        if (count_ >= output_count_limit_) {
                            return true;
                        }
                    } else {
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
                    
                    if (g_exit)
                        return true;
                }

                cur_depth -= 1;

                if (cur_depth <= start_depth) {
                    break;
                } else {
#ifdef COLLECT_INVALID_PR_STATISTICS
                    if (enter_result_count_[cur_depth] == count_) {
                        invalid_leaf_pr_count_ += 1;
                    }
#endif

#ifndef HOMOMORPHISM
                    visited_[embedding_[cur_depth]] = false;
#endif

#ifdef FAILING_SET_PRUNING
                    reverse_embedding_.erase(embedding_[cur_depth]);
                    if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                        vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                        idx_[cur_depth] = num_local_candidates_[cur_depth];
                    } else {
                        vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                    }
#endif
                }

            }
#ifdef COLLECT_INVALID_PR_STATISTICS
            if (count_ == enter_result_count_[start_depth]) {
                if (start_depth == num_core_vertex_ - 1) {
                    invalid_core_pr_count_ += 1;
                }
                else {
                    invalid_leaf_pr_count_ += 1;
                }
            }
#endif

#ifndef HOMOMORPHISM
            visited_[embedding_[start_depth]] = false;
#endif

#ifdef FAILING_SET_PRUNING
            reverse_embedding_.erase(embedding_[start_depth]);
            if (!vec_failing_set_[start_depth].test(vertex_ordering_[start_depth])) {
                vec_failing_set_[start_depth - 1] = vec_failing_set_[start_depth];
                idx_[start_depth] = num_local_candidates_[start_depth];
            } else {
                vec_failing_set_[start_depth - 1] |= vec_failing_set_[start_depth];
            }
#endif
        }
    }

    return false;
}

uint32_t leapfrogtriejoin::compute_local_candidates_core_ptb(uint32_t depth) {
    uint32_t u = vertex_ordering_[depth];
    uint32_t bn_count = get_backward_neighbors_count(depth);
    uint32_t* bns = get_backward_neighbors(depth);
    uint32_t* bns_depth = get_backward_neighbors_depth(depth);
    uint32_t lc_count = 0;
#ifdef INTERSECTION_CACHE
    bool recompute = false;
#endif
    if (bn_count == 1) {
        uint32_t bn = bns[0];
        uint32_t exist_count = 0,ref_count = 0;
#if RELATION_STRUCTURE == 0
        uint32_t key = idx_embedding_[bns_depth[0]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
        uint32_t key = embedding_[bns_depth[0]];
#endif
        uint32_t* candidate = catalog_->get_core_relation_children(bn, u, key, exist_count);

        this->lc_exist_count = exist_count;

        uint32_t ui = vertex_ordering_[catalog_->begin_depth];
        uint32_t uj = vertex_ordering_[catalog_->end_depth];

        uint32_t key_vi = idx_embedding_[catalog_->begin_depth];
        uint32_t* ref = catalog_->get_core_relation_children(ui, uj, key_vi, ref_count);

        for(int i = 0; i < exist_count; ++i){
            if(!exist_check(candidate[i],ref,ref_count)){
                local_candidates_[depth][lc_count] = candidate[i];
                lc_count++;
            }
        }
    }
    else {
#ifndef INTERSECTION_CACHE
        uint32_t bn1 = bns[0];
        uint32_t bn2 = bns[1];

#if RELATION_STRUCTURE == 0
        uint32_t key1 = idx_embedding_[bns_depth[0]];
        uint32_t key2 = idx_embedding_[bns_depth[1]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
        uint32_t key1 = embedding_[bns_depth[0]];
        uint32_t key2 = embedding_[bns_depth[1]];
#endif

#ifndef SPARSE_BITMAP
        uint32_t lc_count1 = 0;
        uint32_t lc_count2 = 0;

        uint32_t *lc1 = catalog_->get_core_relation_children(bn1, u, key1, lc_count1);
        uint32_t *lc2 = catalog_->get_core_relation_children(bn2, u, key2, lc_count2);
#ifdef COLLECT_SI_STATISTICS
        set_intersection_cost_[depth] += std::max(lc_count1, lc_count2);
        set_intersection_count_[depth] += 1;
#endif
        ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, local_candidates_[depth],
                                                  lc_count);

#else
#if RELATION_STRUCTURE == 0
        BSRSet& bsr_set1 = catalog_->bsr_relations_[bn1][u].bsrs[key1];
        BSRSet& bsr_set2 = catalog_->bsr_relations_[bn2][u].bsrs[key2];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
        BSRSet& bsr_set1 = catalog_->bsr_relations_[bn1][u].hash_bsrs_[key1];
        BSRSet& bsr_set2 = catalog_->bsr_relations_[bn2][u].hash_bsrs_[key2];
#endif

#ifdef COLLECT_SI_STATISTICS
        set_intersection_cost_[depth] += std::max(bsr_set1.size_, bsr_set2.size_);
        set_intersection_count_[depth] += 1;
#endif

        cached_bsr_sets_[depth].size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                     bsr_set2.base_, bsr_set2.states_, bsr_set2.size_,
                                                                     cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_);
#endif

        for (uint32_t i = 2; i < bn_count; ++i) {
            uint32_t bn = bns[i];

#if RELATION_STRUCTURE == 0
            uint32_t key = idx_embedding_[bns_depth[i]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            uint32_t key = embedding_[bns_depth[i]];
#endif

#ifndef SPARSE_BITMAP
            lc_count1 = 0;
            lc1 = catalog_->get_core_relation_children(bn, u, key, lc_count1);

            uint32_t temp_lc_count = 0;

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(lc_count1, lc_count);
            set_intersection_count_[depth] += 1;
#endif

            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, local_candidates_[depth], lc_count,
                                                      si_buffer_, temp_lc_count);
            std::swap(local_candidates_[depth], si_buffer_);
            lc_count = temp_lc_count;

            if (lc_count == 0) {
                goto EXIT_EMPTY_SET;
            }

#else
#if RELATION_STRUCTURE == 0
            BSRSet& bsr_set = catalog_->bsr_relations_[bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            BSRSet& bsr_set = catalog_->bsr_relations_[bn][u].hash_bsrs_[key];
#endif

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(bsr_set.size_, cached_bsr_sets_[depth].size_);
            set_intersection_count_[depth] += 1;
#endif

            buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                  cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                  cached_bsr_sets_[depth].size_, buffer_bsr_set1_.base_,
                                                                  buffer_bsr_set1_.states_);

            std::swap(cached_bsr_sets_[depth], buffer_bsr_set1_);
#endif
        }
        this->lc_exist_count = lc_count;
        
        uint32_t bn =  vertex_ordering_[catalog_->begin_depth];
        uint32_t u = vertex_ordering_[catalog_->end_depth];
        uint32_t lc_exist_count = 0;
        uint32_t key = idx_embedding_[catalog_->begin_depth];
        uint32_t *lc_exist = catalog_->get_core_relation_children(bn, u, key, lc_exist_count);

        uint32_t temp_lc_count = 0;
        
        for(int i = 0; i < lc_count; ++i){
            if(!exist_check(local_candidates_[depth][i],lc_exist,lc_exist_count)){
                si_buffer_[temp_lc_count] = local_candidates_[depth][i];
                temp_lc_count++;
            }
        }

        std::swap(local_candidates_[depth], si_buffer_);
        lc_count = temp_lc_count;

#ifdef SPARSE_BITMAP
        lc_count = (uint32_t)offline_bsr_trans_uint(cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                    cached_bsr_sets_[depth].size_, (int*) local_candidates_[depth]);
#endif

#else

        if (intersection_cache_type_[depth] == INTERSECTION_CACHE_TYPE::none || !is_intersection_cache_valid(depth)) {
            uint32_t cached_bn_count = cached_bns_offset_[depth + 1] - cached_bns_offset_[depth];
            uint32_t *cached_bns = cached_bns_ + cached_bns_offset_[depth];
            uint32_t *cached_bns_depth = cached_bns_depth_ + cached_bns_offset_[depth];
            uint32_t cached_bn1 = cached_bns[0];
            uint32_t cached_bn2 = cached_bns[1];

#if RELATION_STRUCTURE == 0
            uint32_t key1 = idx_embedding_[cached_bns_depth[0]];
            uint32_t key2 = idx_embedding_[cached_bns_depth[1]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            uint32_t key1 = embedding_[cached_bns_depth[0]];
            uint32_t key2 = embedding_[cached_bns_depth[1]];
#endif

#ifndef SPARSE_BITMAP
            uint32_t cached_candidates_count = 0;
            uint32_t lc_count1 = 0;

            uint32_t lc_count2 = 0;

            uint32_t *lc1 = catalog_->get_core_relation_children(cached_bn1, u, key1, lc_count1);
            uint32_t *lc2 = catalog_->get_core_relation_children(cached_bn2, u, key2, lc_count2);

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(lc_count1, lc_count2);
            set_intersection_count_[depth] += 1;
#endif
            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, cached_candidates_[depth],
                                                      cached_candidates_count);

#else
#if RELATION_STRUCTURE == 0
            BSRSet& bsr_set1 = catalog_->bsr_relations_[cached_bn1][u].bsrs[key1];
            BSRSet& bsr_set2 = catalog_->bsr_relations_[cached_bn2][u].bsrs[key2];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            BSRSet& bsr_set1 = catalog_->bsr_relations_[cached_bn1][u].hash_bsrs_[key1];
            BSRSet& bsr_set2 = catalog_->bsr_relations_[cached_bn2][u].hash_bsrs_[key2];
#endif

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(bsr_set1.size_, bsr_set2.size_);
            set_intersection_count_[depth] += 1;
#endif

            cached_bsr_sets_[depth].size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                        bsr_set2.base_, bsr_set2.states_, bsr_set2.size_,
                                                                        cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_);
#endif
            for (uint32_t i = 2; i < cached_bn_count; ++i) {
                uint32_t cached_bn = cached_bns[i];

#if RELATION_STRUCTURE == 0
                uint32_t key = idx_embedding_[cached_bns_depth[i]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                uint32_t key = embedding_[cached_bns_depth[i]];
#endif

#ifndef SPARSE_BITMAP
                lc_count1 = 0;
                lc1 = catalog_->get_core_relation_children(cached_bn, u, key, lc_count1);
#ifdef COLLECT_SI_STATISTICS
                set_intersection_cost_[depth] += std::max(lc_count1, cached_candidates_count);
                set_intersection_count_[depth] += 1;
#endif
                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, cached_candidates_[depth], cached_candidates_count,
                                                          si_buffer_, cached_candidates_count);
                std::swap(cached_candidates_[depth], si_buffer_);
#else
#if RELATION_STRUCTURE == 0
                BSRSet& bsr_set = catalog_->bsr_relations_[cached_bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                BSRSet& bsr_set = catalog_->bsr_relations_[cached_bn][u].hash_bsrs_[key];
#endif

#ifdef COLLECT_SI_STATISTICS
                set_intersection_cost_[depth] += std::max(bsr_set.size_, cached_bsr_sets_[depth].size_);
                set_intersection_count_[depth] += 1;
#endif

                buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                      cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                      cached_bsr_sets_[depth].size_, buffer_bsr_set1_.base_,
                                                                      buffer_bsr_set1_.states_);

                std::swap(cached_bsr_sets_[depth], buffer_bsr_set1_);
#endif
            }

#ifndef SPARSE_BITMAP
            cached_candidates_count_[depth] = cached_candidates_count;
#endif
            recompute = true;
        }
        if (intersection_cache_type_[depth] != INTERSECTION_CACHE_TYPE::partial) {
#ifndef SPARSE_BITMAP
            lc_count = cached_candidates_count_[depth];
#else
            lc_count = num_local_candidates_[depth];
#endif
            if (recompute) {
#ifndef SPARSE_BITMAP
                std::swap(cached_candidates_[depth], local_candidates_[depth]);
#else
                lc_count = (uint32_t)offline_bsr_trans_uint(cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                cached_bsr_sets_[depth].size_, (int*) local_candidates_[depth]);
#endif
            }
        }
        else {
            uint32_t no_cached_bn_count = no_cached_bns_offset_[depth + 1] - no_cached_bns_offset_[depth];
            uint32_t* no_cached_bns = no_cached_bns_ + no_cached_bns_offset_[depth];
            uint32_t* no_cached_bns_depth = no_cached_bns_depth_ + no_cached_bns_offset_[depth];

            uint32_t no_cached_bn = no_cached_bns[0];

#if RELATION_STRUCTURE == 0
            uint32_t key = idx_embedding_[no_cached_bns_depth[0]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            uint32_t key = embedding_[no_cached_bns_depth[0]];
#endif

#ifndef SPARSE_BITMAP
            uint32_t lc_count1 = 0;
            uint32_t* lc1 = catalog_->get_core_relation_children(no_cached_bn, u, key, lc_count1);

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(lc_count1, cached_candidates_count_[depth]);
            set_intersection_count_[depth] += 1;
#endif
            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, cached_candidates_[depth], cached_candidates_count_[depth],
                                                      local_candidates_[depth], lc_count);

#else
#if RELATION_STRUCTURE == 0
            BSRSet& bsr_set1 = catalog_->bsr_relations_[no_cached_bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            BSRSet& bsr_set1 = catalog_->bsr_relations_[no_cached_bn][u].hash_bsrs_[key];
#endif

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(bsr_set1.size_, cached_bsr_sets_[depth].size_);
            set_intersection_count_[depth] += 1;
#endif
            buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                 cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_, cached_bsr_sets_[depth].size_,
                                                                 buffer_bsr_set1_.base_, buffer_bsr_set1_.states_);
#endif
            for (uint32_t i = 1; i < no_cached_bn_count; ++i) {
                no_cached_bn = no_cached_bns[i];

#if RELATION_STRUCTURE == 0
                key = idx_embedding_[no_cached_bns_depth[i]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                key = embedding_[no_cached_bns_depth[i]];
#endif

#ifndef SPARSE_BITMAP
                lc_count1 = 0;
                lc1 = catalog_->get_core_relation_children(no_cached_bn, u, key, lc_count1);

#ifdef COLLECT_SI_STATISTICS
                set_intersection_cost_[depth] += std::max(lc_count1, lc_count);
                set_intersection_count_[depth] += 1;
#endif
                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, local_candidates_[depth], lc_count,
                                                          si_buffer_, lc_count);
                std::swap(local_candidates_[depth], si_buffer_);
#else
#if RELATION_STRUCTURE == 0
                BSRSet& bsr_set = catalog_->bsr_relations_[no_cached_bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                BSRSet& bsr_set = catalog_->bsr_relations_[no_cached_bn][u].hash_bsrs_[key];
#endif

#ifdef COLLECT_SI_STATISTICS
                set_intersection_cost_[depth] += std::max(bsr_set.size_, buffer_bsr_set1_.size_);
                set_intersection_count_[depth] += 1;
#endif
                buffer_bsr_set2_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                     buffer_bsr_set1_.base_, buffer_bsr_set1_.states_, buffer_bsr_set1_.size_,
                                                                     buffer_bsr_set2_.base_, buffer_bsr_set2_.states_);

                std::swap(buffer_bsr_set1_, buffer_bsr_set2_);
#endif
            }

#ifdef SPARSE_BITMAP
            lc_count = (uint32_t)offline_bsr_trans_uint(buffer_bsr_set1_.base_, buffer_bsr_set1_.states_,
                                                       buffer_bsr_set1_.size_, (int*) local_candidates_[depth]);
#endif
        }
#endif
    }

    if (lc_count == 0){
        EXIT_EMPTY_SET:
#ifdef COLLECT_FAIL_SI_STATISTICS
        fail_si_count_ += 1;
#endif

        return 0;
    }

    return lc_count;
}

uint32_t leapfrogtriejoin::compute_local_candidates(uint32_t depth) {
    uint32_t u = vertex_ordering_[depth];
    uint32_t bn_count = get_backward_neighbors_count(depth);
    uint32_t* bns = get_backward_neighbors(depth);
    uint32_t* bns_depth = get_backward_neighbors_depth(depth);
    uint32_t lc_count = 0;
#ifdef INTERSECTION_CACHE
    bool recompute = false;
#endif
    if (bn_count == 1) {
        uint32_t bn = bns[0];

#if RELATION_STRUCTURE == 0
        uint32_t key = idx_embedding_[bns_depth[0]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
        uint32_t key = embedding_[bns_depth[0]];
#endif
        local_candidates_[depth] = catalog_->get_core_relation_children(bn, u, key, lc_count);
    }
    else {
#ifndef INTERSECTION_CACHE
        uint32_t bn1 = bns[0];
        uint32_t bn2 = bns[1];

#if RELATION_STRUCTURE == 0
        uint32_t key1 = idx_embedding_[bns_depth[0]];
        uint32_t key2 = idx_embedding_[bns_depth[1]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
        uint32_t key1 = embedding_[bns_depth[0]];
        uint32_t key2 = embedding_[bns_depth[1]];
#endif

#ifndef SPARSE_BITMAP
        uint32_t lc_count1 = 0;
        uint32_t lc_count2 = 0;

        uint32_t *lc1 = catalog_->get_core_relation_children(bn1, u, key1, lc_count1);
        uint32_t *lc2 = catalog_->get_core_relation_children(bn2, u, key2, lc_count2);
#ifdef COLLECT_SI_STATISTICS
        set_intersection_cost_[depth] += std::max(lc_count1, lc_count2);
        set_intersection_count_[depth] += 1;
#endif
        ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, local_candidates_[depth],
                                                  lc_count);

#else
#if RELATION_STRUCTURE == 0
        BSRSet& bsr_set1 = catalog_->bsr_relations_[bn1][u].bsrs[key1];
        BSRSet& bsr_set2 = catalog_->bsr_relations_[bn2][u].bsrs[key2];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
        BSRSet& bsr_set1 = catalog_->bsr_relations_[bn1][u].hash_bsrs_[key1];
        BSRSet& bsr_set2 = catalog_->bsr_relations_[bn2][u].hash_bsrs_[key2];
#endif

#ifdef COLLECT_SI_STATISTICS
        set_intersection_cost_[depth] += std::max(bsr_set1.size_, bsr_set2.size_);
        set_intersection_count_[depth] += 1;
#endif

        cached_bsr_sets_[depth].size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                     bsr_set2.base_, bsr_set2.states_, bsr_set2.size_,
                                                                     cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_);
#endif

        for (uint32_t i = 2; i < bn_count; ++i) {
            uint32_t bn = bns[i];

#if RELATION_STRUCTURE == 0
            uint32_t key = idx_embedding_[bns_depth[i]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            uint32_t key = embedding_[bns_depth[i]];
#endif

#ifndef SPARSE_BITMAP
            lc_count1 = 0;
            lc1 = catalog_->get_core_relation_children(bn, u, key, lc_count1);

            uint32_t temp_lc_count = 0;

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(lc_count1, lc_count);
            set_intersection_count_[depth] += 1;
#endif

            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, local_candidates_[depth], lc_count,
                                                      si_buffer_, temp_lc_count);
            std::swap(local_candidates_[depth], si_buffer_);
            lc_count = temp_lc_count;

            if (lc_count == 0) {
                goto EXIT_EMPTY_SET;
            }
#else
#if RELATION_STRUCTURE == 0
            BSRSet& bsr_set = catalog_->bsr_relations_[bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            BSRSet& bsr_set = catalog_->bsr_relations_[bn][u].hash_bsrs_[key];
#endif

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(bsr_set.size_, cached_bsr_sets_[depth].size_);
            set_intersection_count_[depth] += 1;
#endif

            buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                  cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                  cached_bsr_sets_[depth].size_, buffer_bsr_set1_.base_,
                                                                  buffer_bsr_set1_.states_);

            std::swap(cached_bsr_sets_[depth], buffer_bsr_set1_);
#endif
        }

#ifdef SPARSE_BITMAP
        lc_count = (uint32_t)offline_bsr_trans_uint(cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                    cached_bsr_sets_[depth].size_, (int*) local_candidates_[depth]);
#endif

#else

        if (intersection_cache_type_[depth] == INTERSECTION_CACHE_TYPE::none || !is_intersection_cache_valid(depth)) {
            uint32_t cached_bn_count = cached_bns_offset_[depth + 1] - cached_bns_offset_[depth];
            uint32_t *cached_bns = cached_bns_ + cached_bns_offset_[depth];
            uint32_t *cached_bns_depth = cached_bns_depth_ + cached_bns_offset_[depth];
            uint32_t cached_bn1 = cached_bns[0];
            uint32_t cached_bn2 = cached_bns[1];

#if RELATION_STRUCTURE == 0
            uint32_t key1 = idx_embedding_[cached_bns_depth[0]];
            uint32_t key2 = idx_embedding_[cached_bns_depth[1]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            uint32_t key1 = embedding_[cached_bns_depth[0]];
            uint32_t key2 = embedding_[cached_bns_depth[1]];
#endif

#ifndef SPARSE_BITMAP
            uint32_t cached_candidates_count = 0;
            uint32_t lc_count1 = 0;

            uint32_t lc_count2 = 0;

            uint32_t *lc1 = catalog_->get_core_relation_children(cached_bn1, u, key1, lc_count1);
            uint32_t *lc2 = catalog_->get_core_relation_children(cached_bn2, u, key2, lc_count2);

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(lc_count1, lc_count2);
            set_intersection_count_[depth] += 1;
#endif
            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, cached_candidates_[depth],
                                                      cached_candidates_count);

#else
#if RELATION_STRUCTURE == 0
            BSRSet& bsr_set1 = catalog_->bsr_relations_[cached_bn1][u].bsrs[key1];
            BSRSet& bsr_set2 = catalog_->bsr_relations_[cached_bn2][u].bsrs[key2];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            BSRSet& bsr_set1 = catalog_->bsr_relations_[cached_bn1][u].hash_bsrs_[key1];
            BSRSet& bsr_set2 = catalog_->bsr_relations_[cached_bn2][u].hash_bsrs_[key2];
#endif

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(bsr_set1.size_, bsr_set2.size_);
            set_intersection_count_[depth] += 1;
#endif

            cached_bsr_sets_[depth].size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                        bsr_set2.base_, bsr_set2.states_, bsr_set2.size_,
                                                                        cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_);
#endif
            for (uint32_t i = 2; i < cached_bn_count; ++i) {
                uint32_t cached_bn = cached_bns[i];

#if RELATION_STRUCTURE == 0
                uint32_t key = idx_embedding_[cached_bns_depth[i]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                uint32_t key = embedding_[cached_bns_depth[i]];
#endif

#ifndef SPARSE_BITMAP
                lc_count1 = 0;
                lc1 = catalog_->get_core_relation_children(cached_bn, u, key, lc_count1);
#ifdef COLLECT_SI_STATISTICS
                set_intersection_cost_[depth] += std::max(lc_count1, cached_candidates_count);
                set_intersection_count_[depth] += 1;
#endif
                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, cached_candidates_[depth], cached_candidates_count,
                                                          si_buffer_, cached_candidates_count);
                std::swap(cached_candidates_[depth], si_buffer_);
#else
#if RELATION_STRUCTURE == 0
                BSRSet& bsr_set = catalog_->bsr_relations_[cached_bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                BSRSet& bsr_set = catalog_->bsr_relations_[cached_bn][u].hash_bsrs_[key];
#endif

#ifdef COLLECT_SI_STATISTICS
                set_intersection_cost_[depth] += std::max(bsr_set.size_, cached_bsr_sets_[depth].size_);
                set_intersection_count_[depth] += 1;
#endif

                buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                      cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                      cached_bsr_sets_[depth].size_, buffer_bsr_set1_.base_,
                                                                      buffer_bsr_set1_.states_);

                std::swap(cached_bsr_sets_[depth], buffer_bsr_set1_);
#endif
            }

#ifndef SPARSE_BITMAP
            cached_candidates_count_[depth] = cached_candidates_count;
#endif
            recompute = true;
        }
        if (intersection_cache_type_[depth] != INTERSECTION_CACHE_TYPE::partial) {
#ifndef SPARSE_BITMAP
            lc_count = cached_candidates_count_[depth];
#else
            lc_count = num_local_candidates_[depth];
#endif
            if (recompute) {
#ifndef SPARSE_BITMAP
                std::swap(cached_candidates_[depth], local_candidates_[depth]);
#else
                lc_count = (uint32_t)offline_bsr_trans_uint(cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                cached_bsr_sets_[depth].size_, (int*) local_candidates_[depth]);
#endif
            }
        }
        else {
            uint32_t no_cached_bn_count = no_cached_bns_offset_[depth + 1] - no_cached_bns_offset_[depth];
            uint32_t* no_cached_bns = no_cached_bns_ + no_cached_bns_offset_[depth];
            uint32_t* no_cached_bns_depth = no_cached_bns_depth_ + no_cached_bns_offset_[depth];

            uint32_t no_cached_bn = no_cached_bns[0];

#if RELATION_STRUCTURE == 0
            uint32_t key = idx_embedding_[no_cached_bns_depth[0]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            uint32_t key = embedding_[no_cached_bns_depth[0]];
#endif

#ifndef SPARSE_BITMAP
            uint32_t lc_count1 = 0;
            uint32_t* lc1 = catalog_->get_core_relation_children(no_cached_bn, u, key, lc_count1);

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(lc_count1, cached_candidates_count_[depth]);
            set_intersection_count_[depth] += 1;
#endif
            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, cached_candidates_[depth], cached_candidates_count_[depth],
                                                      local_candidates_[depth], lc_count);

#else
#if RELATION_STRUCTURE == 0
            BSRSet& bsr_set1 = catalog_->bsr_relations_[no_cached_bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
            BSRSet& bsr_set1 = catalog_->bsr_relations_[no_cached_bn][u].hash_bsrs_[key];
#endif

#ifdef COLLECT_SI_STATISTICS
            set_intersection_cost_[depth] += std::max(bsr_set1.size_, cached_bsr_sets_[depth].size_);
            set_intersection_count_[depth] += 1;
#endif
            buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                 cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_, cached_bsr_sets_[depth].size_,
                                                                 buffer_bsr_set1_.base_, buffer_bsr_set1_.states_);
#endif
            for (uint32_t i = 1; i < no_cached_bn_count; ++i) {
                no_cached_bn = no_cached_bns[i];

#if RELATION_STRUCTURE == 0
                key = idx_embedding_[no_cached_bns_depth[i]];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                key = embedding_[no_cached_bns_depth[i]];
#endif

#ifndef SPARSE_BITMAP
                lc_count1 = 0;
                lc1 = catalog_->get_core_relation_children(no_cached_bn, u, key, lc_count1);

#ifdef COLLECT_SI_STATISTICS
                set_intersection_cost_[depth] += std::max(lc_count1, lc_count);
                set_intersection_count_[depth] += 1;
#endif
                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, local_candidates_[depth], lc_count,
                                                          si_buffer_, lc_count);
                std::swap(local_candidates_[depth], si_buffer_);
#else
#if RELATION_STRUCTURE == 0
                BSRSet& bsr_set = catalog_->bsr_relations_[no_cached_bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1 || RELATION_STRUCTURE == 2
                BSRSet& bsr_set = catalog_->bsr_relations_[no_cached_bn][u].hash_bsrs_[key];
#endif

#ifdef COLLECT_SI_STATISTICS
                set_intersection_cost_[depth] += std::max(bsr_set.size_, buffer_bsr_set1_.size_);
                set_intersection_count_[depth] += 1;
#endif
                buffer_bsr_set2_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                     buffer_bsr_set1_.base_, buffer_bsr_set1_.states_, buffer_bsr_set1_.size_,
                                                                     buffer_bsr_set2_.base_, buffer_bsr_set2_.states_);

                std::swap(buffer_bsr_set1_, buffer_bsr_set2_);
#endif
            }

#ifdef SPARSE_BITMAP
            lc_count = (uint32_t)offline_bsr_trans_uint(buffer_bsr_set1_.base_, buffer_bsr_set1_.states_,
                                                       buffer_bsr_set1_.size_, (int*) local_candidates_[depth]);
#endif
        }
#endif
    }

    if (lc_count == 0){
        EXIT_EMPTY_SET:
#ifdef COLLECT_FAIL_SI_STATISTICS
        fail_si_count_ += 1;
#endif

        return 0;
    }

    return lc_count;
}

bool leapfrogtriejoin::exist_check(uint32_t u, uint32_t *candidates, uint32_t count){
    int begin = 0;
    int end = count - 1;
    while (begin <= end) {
        int mid = begin + ((end - begin) >> 1);
        if (candidates[mid] == u) {
            return true;
        }
        else if (candidates[mid] > u)
            end = mid - 1;
        else
            begin = mid + 1;
    }

    return false;
}

void leapfrogtriejoin::initialize() {
    fail_si_count_ = 0;
    iso_conflict_count_ = 0;
    invalid_core_pr_count_ = 0;
    invalid_leaf_pr_count_ = 0;
    enter_result_count_ = new uint64_t[num_vertex_];
    enter_ptb_count_ = new uint64_t[num_vertex_];
    std::fill(enter_result_count_, enter_result_count_ + num_vertex_, 0);
    std::fill(enter_ptb_count_, enter_ptb_count_ + num_vertex_, 0);

    idx_ = new uint32_t[num_vertex_];
    num_local_candidates_ = new uint32_t[num_vertex_];
    num_ptb_candidates_ = new uint32_t[num_vertex_];
    embedding_ = new uint32_t[num_vertex_];
    si_buffer_ = new uint32_t[catalog_->max_num_candidates_per_vertex_];
    si_buffer_precomputed = new uint32_t[catalog_->max_num_candidates_per_vertex_];
    last_level_candidates_ = new uint32_t[catalog_->max_num_candidates_per_vertex_];

    bns_offset_ = new uint32_t[num_vertex_ + 1];
    fns_offset_ = new uint32_t[num_vertex_ + 1];
    bns_ = new uint32_t[num_vertex_ * num_vertex_];
    bns_depth_ = new uint32_t[num_vertex_ * num_vertex_];
    fns_ = new uint32_t[num_vertex_ * num_vertex_];
    fns_depth_ = new uint32_t[num_vertex_ * num_vertex_];

    set_intersection_cost_ = new uint64_t[num_vertex_];
    set_intersection_count_ = new uint64_t[num_vertex_];
    intermediate_result_count_ = new uint64_t[num_vertex_];
    std::fill(set_intersection_cost_, set_intersection_cost_ + num_vertex_, 0);
    std::fill(set_intersection_count_, set_intersection_count_ + num_vertex_, 0);
    std::fill(intermediate_result_count_, intermediate_result_count_ + num_vertex_, 0);

#ifdef INTERSECTION_CACHE
    cached_bns_offset_ = new uint32_t[num_vertex_ + 1];
    cached_bns_ = new uint32_t[num_vertex_ * num_vertex_];
    cached_bns_depth_ = new uint32_t[num_vertex_ * num_vertex_];
    no_cached_bns_offset_ = new uint32_t[num_vertex_ + 1];
    no_cached_bns_ = new uint32_t[num_vertex_ * num_vertex_];
    no_cached_bns_depth_ = new uint32_t[num_vertex_ * num_vertex_];
    intersection_cache_type_ = new INTERSECTION_CACHE_TYPE[num_vertex_];
    mapping_to_cached_bns_ = new uint32_t[num_vertex_ * num_vertex_];
    std::fill(mapping_to_cached_bns_, mapping_to_cached_bns_ + num_vertex_ * num_vertex_, catalog_->max_data_vertex_id_ + 1);
#endif

#ifndef HOMOMORPHISM
    visited_ = new bool[catalog_->max_data_vertex_id_ + 1];
    memset(visited_, 0, sizeof(bool) * (catalog_->max_data_vertex_id_ + 1));
#endif

#if RELATION_STRUCTURE == 0
    idx_embedding_ = new uint32_t[num_vertex_];
#endif

#ifdef SPARSE_BITMAP
    cached_bsr_sets_ = new BSRSet[num_vertex_];
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        cached_bsr_sets_[i].states_ = new int[catalog_->max_num_candidates_per_vertex_];
        cached_bsr_sets_[i].base_ = new int[catalog_->max_num_candidates_per_vertex_];
    }
    buffer_bsr_set1_.states_ = new int[catalog_->max_num_candidates_per_vertex_];
    buffer_bsr_set1_.base_ = new int[catalog_->max_num_candidates_per_vertex_];
    buffer_bsr_set2_.states_ = new int[catalog_->max_num_candidates_per_vertex_];
    buffer_bsr_set2_.base_ = new int[catalog_->max_num_candidates_per_vertex_];
#endif

    if (output_ != nullptr) {
        materialize_ = true;
    }

#ifdef ENABLE_OUTPUT
    current_position_ = output_buffer;
#endif
}

void leapfrogtriejoin::clear() {
    delete[] idx_;
    delete[] num_local_candidates_;
    delete[] num_ptb_candidates_;
    delete[] embedding_;
    delete[] si_buffer_;
    delete[] si_buffer_precomputed;
    delete[] last_level_candidates_;
    delete[] enter_result_count_;
    delete[] enter_ptb_count_;

    if(catalog_->alg == "JointEnum"|| catalog_->alg == "RWC" ){
        for (uint32_t i = 0; i < num_vertex_; ++i) {
            uint32_t bn_count = bns_offset_[i + 1] - bns_offset_[i];
            if (bn_count > 1 || i == catalog_->end_depth) {
                delete[] local_candidates_[i];
            }
        }
    }
    else if(catalog_->alg == "LookAheadEnum"|| catalog_->alg == "RWCH"){
        for (uint32_t i = 0; i < num_vertex_; ++i) {
            uint32_t bn_count = bns_offset_[i + 1] - bns_offset_[i];
            if (bn_count > 1 || i == catalog_->end_depth) {
                delete[] local_candidates_[i];
            }
            if (i < catalog_->end_depth && end_node_bn.test(i)) {
                delete[] ptb_candidates_[i];
            }
        }
    }
    else{
        for (uint32_t i = 0; i < num_vertex_; ++i) {
            uint32_t bn_count = bns_offset_[i + 1] - bns_offset_[i];
            if (bn_count > 1) {
                delete[] local_candidates_[i];
            }
        }
    }

    delete[] local_candidates_;
    delete[] ptb_candidates_;

    delete[] bns_offset_;
    delete[] fns_offset_;
    delete[] bns_;
    delete[] fns_;
    delete[] bns_depth_;
    delete[] fns_depth_;

    delete[] set_intersection_cost_;
    delete[] intermediate_result_count_;
    delete[] set_intersection_count_;

#ifdef INTERSECTION_CACHE
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        delete[] cached_candidates_[i];
    }

    delete[] cached_candidates_;
    delete[] cached_candidates_count_;
    delete[] cached_bns_offset_;
    delete[] cached_bns_;
    delete[] cached_bns_depth_;
    delete[] no_cached_bns_offset_;
    delete[] no_cached_bns_;
    delete[] no_cached_bns_depth_;
    delete[] intersection_cache_type_;
    delete[] mapping_to_cached_bns_;
#endif

#ifndef HOMOMORPHISM
    delete[] visited_;
#endif

#ifdef SPARSE_BITMAP
    delete[] cached_bsr_sets_;
#endif

#if RELATION_STRUCTURE == 0
    delete[] idx_embedding_;
#endif
}

#if RELATION_STRUCTURE == 0
void leapfrogtriejoin::set_idx(uint32_t depth) {
    for (uint32_t i = 0; i < depth; ++i) {
        uint32_t u = vertex_ordering_[i];
        if (get_forward_neighbors_count(i) != 0) {
            idx_embedding_[i] = catalog_->get_candidate_index(u, embedding_[i]);
        }
    }
}
#endif

void leapfrogtriejoin::initialize_bn_fn() {
    uint32_t bn_offset = 0;
    uint32_t fn_offset = 0;
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        uint32_t u = vertex_ordering_[i];
        bns_offset_[i] = bn_offset;
        fns_offset_[i] = fn_offset;
        for (uint32_t j = 0; j < num_vertex_; ++j) {
            uint32_t v = vertex_ordering_[j];
            if (catalog_->is_adjacent(u, v)) {
                if (j < i) {
                    bns_[bn_offset] = v;
                    bns_depth_[bn_offset++] = j;
                    if(i == catalog_->end_depth) end_node_bn.set(j);
                }
                else if (j > i) {
                    fns_[fn_offset] = v;
                    fns_depth_[fn_offset++] = j;
                }
            }
        }
    }
    end_node_bn.set(catalog_->begin_depth);
    bns_offset_[num_vertex_] = bn_offset;
    fns_offset_[num_vertex_] = fn_offset;

    local_candidates_ = new uint32_t*[num_vertex_];
    ptb_candidates_ = new uint32_t*[num_vertex_];

    if(catalog_->alg == "JointEnum" || catalog_->alg == "RWC" ){
        for (uint32_t i = 0; i < num_vertex_; ++i) {
            uint32_t bn_count = bns_offset_[i + 1] - bns_offset_[i];

            if ((i >= input_->get_tuple_length() && bn_count > 1) || (i == catalog_-> end_depth && catalog_-> end_depth > 0)) {
                local_candidates_[i] = new uint32_t[catalog_->max_num_candidates_per_vertex_];
            }
            else {
                local_candidates_[i] = nullptr;
            }
        }
    }
    else if(catalog_->alg == "LookAheadEnum" || catalog_->alg == "RWCH"){
        for (uint32_t i = 0; i < num_vertex_; ++i) {
            uint32_t bn_count = bns_offset_[i + 1] - bns_offset_[i];

            if ((i >= input_->get_tuple_length() && bn_count > 1) || (i == catalog_-> end_depth && catalog_-> end_depth > 0)) {
                local_candidates_[i] = new uint32_t[catalog_->max_num_candidates_per_vertex_];
            }
            else {
                local_candidates_[i] = nullptr;
            }

            if (i < catalog_->end_depth && end_node_bn.test(i)) {
                ptb_candidates_[i] = new uint32_t[catalog_->max_num_candidates_per_vertex_];
            }
            else {
                ptb_candidates_[i] = nullptr;
            }
        }
    }
    else{
        for (uint32_t i = 0; i < num_vertex_; ++i) {
            uint32_t bn_count = bns_offset_[i + 1] - bns_offset_[i];

            if (i >= input_->get_tuple_length() && bn_count > 1) {
                local_candidates_[i] = new uint32_t[catalog_->max_num_candidates_per_vertex_];
            }
            else {
                local_candidates_[i] = nullptr;
            }
        }
    }

}

#ifdef INTERSECTION_CACHE

void leapfrogtriejoin::initialize_intersection_cache() {
    uint32_t cached_bn_offset = 0;
    uint32_t no_cached_bn_offset = 0;

    for (uint32_t i = 0; i < num_vertex_; ++i) {
        cached_bns_offset_[i] = cached_bn_offset;
        no_cached_bns_offset_[i] = no_cached_bn_offset;
        intersection_cache_type_[i] = INTERSECTION_CACHE_TYPE::none;

        // Have at least two backward neighbors.
        if (bns_offset_[i + 1] - bns_offset_[i] >= 2) {
            for (uint32_t j = bns_offset_[i]; j < bns_offset_[i + 1]; ++j) {
                if (bns_depth_[j] < i - 1) {
                    cached_bns_depth_[cached_bn_offset] = bns_depth_[j];
                    cached_bns_[cached_bn_offset++] = bns_[j];
                } else {
                    no_cached_bns_depth_[no_cached_bn_offset] = bns_depth_[j];
                    no_cached_bns_[no_cached_bn_offset++] = bns_[j];
                }
            }

            // If only the depth of the first backward neighbor is lower than i - 1, then
            // intersection cache is none. Otherwise, it is partial or full.
            if (cached_bn_offset - cached_bns_offset_[i] == 1) {
                no_cached_bn_offset = no_cached_bns_offset_[i];
                cached_bns_depth_[cached_bn_offset] = no_cached_bns_depth_[no_cached_bn_offset];
                cached_bns_[cached_bn_offset++] = no_cached_bns_[no_cached_bn_offset];
            } else {
                if (no_cached_bn_offset - no_cached_bns_offset_[i] > 0) {
                    intersection_cache_type_[i] = INTERSECTION_CACHE_TYPE::partial;
                } else {
                    intersection_cache_type_[i] = INTERSECTION_CACHE_TYPE::full;
                }
            }
        }
    }
    cached_bns_offset_[num_vertex_] = cached_bn_offset;
    no_cached_bns_offset_[num_vertex_] = no_cached_bn_offset;

    cached_candidates_count_ = new uint32_t[num_vertex_];
    cached_candidates_ = new uint32_t*[num_vertex_];
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        if (i >= input_->get_tuple_length() && bns_offset_[i + 1] - bns_offset_[i] >= 2) {
            cached_candidates_[i] = new uint32_t[catalog_->max_num_candidates_per_vertex_];
        }
        else {
            cached_candidates_[i] = nullptr;
        }
    }
}

bool leapfrogtriejoin::is_intersection_cache_valid(uint32_t depth) {
    bool valid = true;
    for (uint32_t i = cached_bns_offset_[depth]; i < cached_bns_offset_[depth + 1]; ++i) {
        if (mapping_to_cached_bns_[i] != embedding_[cached_bns_depth_[i]]) {
            valid = false;
            mapping_to_cached_bns_[i] = embedding_[cached_bns_depth_[i]];
        }
    }
    return valid;
}

#endif

#ifdef FAILING_SET_PRUNING

void leapfrogtriejoin::initialize_failing_set_pruning() {
    ancestors_.resize(num_vertex_);
    vec_failing_set_.resize(num_vertex_);
    reverse_embedding_.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    compute_ancestors();
}

void leapfrogtriejoin::compute_ancestors() {
    // Compute the ancestor in the top-down order.
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        uint32_t u = vertex_ordering_[i];
        ancestors_[u].set(u);

        for (uint32_t j = bns_offset_[i]; j < bns_offset_[i + 1]; ++j) {
            uint32_t u_bn = bns_[j];
            ancestors_[u] |= ancestors_[u_bn];
        }

        if(i == catalog_->end_depth && catalog_->end_depth > 0) ancestors_[u] |= ancestors_[vertex_ordering_[catalog_->begin_depth]];
    }
}

#endif

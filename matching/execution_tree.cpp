#include "execution_tree.h"

void execution_tree::import_plan(const std::string plan) {
    std::stringstream ss(plan);
    cereal::JSONInputArchive archive(ss);
    // Get operator num.
    uint32_t num_operators;
    archive(num_operators);
    assert(num_operators > 0);

    // Get operators.
    for (uint32_t i = 0; i < num_operators; ++i) {
        join_operator_type operator_type;
        archive(operator_type);
        if (operator_type == join_operator_type::no_partition_hash_join) {
            no_partition_hash_join_node* node = new no_partition_hash_join_node();
            archive(*node);
            tree_nodes_.push_back(node);
        }
        else if (operator_type == join_operator_type::leap_frog_trie_join) {
            leap_frog_trie_join_node* node = new leap_frog_trie_join_node();
            archive(*node);
            tree_nodes_.push_back(node);
        }
    }
}

std::string execution_tree::export_plan() {
    std::stringstream ss;
    {
        cereal::JSONOutputArchive archive(ss);
        auto num_operators = static_cast<uint32_t>(tree_nodes_.size());
        archive(cereal::make_nvp("Num of Operators", num_operators));
        for (auto node : tree_nodes_) {
            join_operator_type operator_type = node->get_operator_type();
            archive(cereal::make_nvp("Operator type", operator_type));
            if (operator_type == join_operator_type::no_partition_hash_join) {
                archive(cereal::make_nvp("No Partition Hash Join", *(no_partition_hash_join_node *) node));
            } else if (operator_type == join_operator_type::leap_frog_trie_join) {
                archive(cereal::make_nvp("Leap Frog Trie Join", *(leap_frog_trie_join_node *) node));
            }
        }
    }
    return ss.str();
}

uint64_t execution_tree::execute(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_nc(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_nc(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_je(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_je(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_lae(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_lae(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_probe_mat(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_probe_mat(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_pe(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_pe(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_iee(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_iee(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_rwc(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_rwc(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_rwch(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_rwch(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_rwd(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_rwd(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_rwdw_epl(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_rwdw_epl(catalog, relations_, output_count_limit);
    }
    return count;
}

uint64_t execution_tree::execute_rwdw_final(catalog &catalog, uint64_t output_count_limit) {
    uint64_t count = 0;
    for (auto node : tree_nodes_) {
        count = node->execute_rwdw_final(catalog, relations_, output_count_limit);
    }
    return count;
}
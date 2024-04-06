#pragma once
#include <stdexcept>

template <typename key_type, typename value_type, typename weight_type>
typename graph::Graph<key_type, value_type, weight_type>& graph::Graph<key_type, value_type, weight_type>::operator=(const Graph& other) {
    if (this != &other) {
        m_graph_map = other.m_graph_map;
    }
    return *this;
}

template <typename key_type, typename value_type, typename weight_type>
typename graph::Graph<key_type, value_type, weight_type>& graph::Graph<key_type, value_type, weight_type>::operator=(Graph&& other) noexcept {
    if (this != &other) {
        m_graph_map = std::move(other.m_graph_map);
    }
    return *this;
}

template <typename key_type, typename value_type, typename weight_type>
typename graph::Graph<key_type, value_type, weight_type>::Node& graph::Graph<key_type, value_type, weight_type>::Node::operator=(const Node& other) {
    if (this != &other) {
        m_value = other.m_value;
        m_edge = other.m_edge;
    }
    return *this;
}

template <typename key_type, typename value_type, typename weight_type>
typename graph::Graph<key_type, value_type, weight_type>::Node& graph::Graph<key_type, value_type, weight_type>::Node::operator=(Node&& other) noexcept {
    if (this != &other) {
        m_value = std::move(other.m_value);
        m_edge = std::move(other.m_edge);
    }
    return *this;
}

template <typename key_type, typename value_type, typename weight_type>
void graph::Graph<key_type, value_type, weight_type>::Node::swap(Node& obj) noexcept {
    std::swap(m_edge, obj.edge());
    std::swap(m_value, obj.value());
}

template <typename key_type, typename value_type, typename weight_type>//
size_t graph::Graph<key_type, value_type, weight_type>::degree_in(const key_type& key) const {
    size_t count = 0;
    if (m_graph_map.count(key) == 0)
        throw std::runtime_error("There is no node with this key");

    for (const auto& node_pair : m_graph_map) {
        count += node_pair.second.m_edge.count(key);
    }
    return count;
}

template <typename key_type, typename value_type, typename weight_type>//
size_t graph::Graph<key_type, value_type, weight_type>::degree_out(const key_type& key) const {
    if (m_graph_map.count(key) == 0)
        throw std::runtime_error("There is no node with this key");
    return m_graph_map.at(key).m_edge.size();
}

template <typename key_type, typename value_type, typename weight_type>//
bool graph::Graph<key_type, value_type, weight_type>::loop(const key_type& key) const {
    if (m_graph_map.count(key) == 0)
        throw std::runtime_error("There is no node with this key");
    return m_graph_map.at(key).m_edge.count(key) > 0;
}
template <typename key_type, typename value_type, typename weight_type>//
std::pair<typename graph::Graph<key_type, value_type, weight_type>::Node::iterator, bool> graph::Graph<key_type, value_type, weight_type>::insert_edge(std::pair<key_type, key_type> keys, const weight_type& weight) {
    if (m_graph_map.count(keys.first) == 0 || m_graph_map.count(keys.second) == 0)
        throw std::runtime_error("One of the keys is missing");
    return m_graph_map.at(keys.first).m_edge.insert({keys.second, weight});
}

template <typename key_type, typename value_type, typename weight_type>//
std::pair<typename graph::Graph<key_type, value_type, weight_type>::Node::iterator, bool> graph::Graph<key_type, value_type, weight_type>::insert_or_assign_edge(std::pair<key_type, key_type> keys, const weight_type& weight) {
    if (m_graph_map.count(keys.first) == 0 || m_graph_map.count(keys.second) == 0)
        throw std::runtime_error("One of the keys is missing");

    return m_graph_map.at(keys.first).m_edge.insert_or_assign(keys.second, weight);
}

template <typename key_type, typename value_type, typename weight_type>
value_type& graph::Graph<key_type, value_type, weight_type>::at(const key_type& key) {
    iterator it = m_nodes.find(key);
    if (it == m_nodes.end()) {
        throw std::runtime_error("Cannot find key");
    }
    return (*it).second.value();
}


template <typename key_type, typename value_type, typename weight_type>
bool graph::Graph<key_type, value_type, weight_type>::remove_edge(const key_type& key1, const key_type& key2) {
    auto it = m_graph_map.find(key1);
    if (it != m_graph_map.end()) {
        // Узел, из которого удаляем ребро, найден
        auto& edges = it->second.m_edge;
        auto edge_it = edges.find(key2);
        if (edge_it != edges.end()) {
            edges.erase(edge_it);
            return true; // Успешно удалено
        }
    }
    return false; // Исходное ребро не найдено
}

template <typename key_type, typename value_type, typename weight_type>
bool graph::Graph<key_type, value_type, weight_type>::remove_node(const key_type& key) {
    auto it = m_graph_map.find(key);

    if (it != m_graph_map.end()) {
        for (auto& node : m_graph_map) {
            node.second.m_edge.erase(key);
        }

        m_graph_map.erase(key);

        return true;  // Нода успешно удалена
    }

    return false;  // Нода не найдена
}

template <typename key_type, typename value_type, typename weight_type>
void graph::Graph<key_type, value_type, weight_type>::insert_node(const std::initializer_list<std::pair<key_type, value_type>>& list) {
    for (const auto& element : list) {
        m_graph_map.insert(element);
    }
}

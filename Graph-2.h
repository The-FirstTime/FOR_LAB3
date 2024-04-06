#pragma once
#include <unordered_map>

namespace graph {
    template <typename key_type, typename value_type, typename weight_type>
    class Graph {
    public:
        class Node;

        Graph() = default;
        Graph(const Graph& other) : m_graph_map(other.m_graph_map) {}
        Graph(Graph&& other) noexcept : m_graph_map(std::move(other.m_graph_map)) {}

        Graph& operator=(const Graph& other);
        Graph& operator=(Graph&& other) noexcept;

        using iterator = typename std::unordered_map<key_type, Node>::iterator;
        using const_iterator = typename std::unordered_map<key_type, Node>::const_iterator;

        iterator begin() noexcept { return m_graph_map.begin(); }
        iterator end() noexcept { return m_graph_map.end(); }
        iterator find(const key_type& key) { return m_graph_map.find(key); }
        const_iterator begin() const noexcept { return m_graph_map.begin(); }
        const_iterator end() const noexcept { return m_graph_map.end(); }
        const_iterator cbegin() const noexcept { return m_graph_map.cbegin(); }
        const_iterator cend() const noexcept { return m_graph_map.cend(); }

        bool empty() const noexcept { return m_graph_map.empty(); }
        size_t size() const noexcept { return m_graph_map.size(); }
        void clear() noexcept { m_graph_map.clear(); }
        void swap(Graph& obj) noexcept { std::swap(m_graph_map, obj.m_graph_map); }

        size_t degree_in(const key_type& key) const;
        size_t degree_out(const key_type& key) const;
        bool loop(const key_type& key) const;

        std::pair<iterator, bool> insert_node(const key_type& key, const value_type& val) { return m_graph_map.insert({ key, Node(val) }); }
        std::pair<iterator, bool> insert_or_assign_node(const key_type& key, const value_type& val) { return m_graph_map.insert_or_assign(key, Node(val)); }
        std::pair<typename Node::iterator, bool> insert_edge(std::pair<key_type, key_type> keys, const weight_type& weight);//да
        std::pair<typename Node::iterator, bool> insert_or_assign_edge(std::pair<key_type, key_type> keys, const weight_type& weight);//да

        value_type& operator[](const key_type& key) { return m_graph_map[key].value(); }
        value_type& at(const key_type& key);

        bool remove_edge(const key_type& key1, const key_type& key2);
        bool remove_node(const key_type& key);

        void insert_node(const std::initializer_list<std::pair<key_type, value_type>>& list);


    private:
        std::unordered_map<key_type, Node> m_graph_map;
    };
    template <typename key_type, typename value_type, typename weight_type>
    class Graph<key_type, value_type, weight_type>::Node {
    public:
        friend class Graph;

        using iterator = typename std::unordered_map<key_type, weight_type>::iterator;
        using const_iterator = typename std::unordered_map<key_type, weight_type>::const_iterator;

        Node() = default;
        Node(value_type value, const std::unordered_map<key_type, weight_type>& edge = {}) :m_value(value), m_edge(edge) {}
        Node(const Node& other) : m_value(other.m_value), m_edge(other.m_edge) {}
        Node(Node&& other) noexcept : m_value(std::move(other.m_value)), m_edge(std::move(other.m_edge)) {}

        Node& operator=(const Node& other);
        Node& operator= (Node&& other) noexcept;

        bool empty() const noexcept { return m_edge.empty(); }
        size_t size() const noexcept { return m_edge.size(); }
        value_type& value() noexcept { return m_value; }
        const value_type& value() const noexcept { return m_value; }
        void clear() noexcept { m_edge.clear(); }
        void swap(Node& obj) noexcept;

        iterator begin() noexcept { return m_edge.begin(); }//коммит
        iterator end() noexcept { return m_edge.end(); }
        const_iterator begin() const noexcept { return m_edge.begin(); }
        const_iterator end() const noexcept { return m_edge.end(); }
        const_iterator cbegin() const noexcept { return m_edge.cbegin(); }
        const_iterator cend() const noexcept { return m_edge.cend(); }

    private:
        value_type m_value;
        std::unordered_map<key_type, weight_type> m_edge;
    };
}
#include "Graph.hpp"
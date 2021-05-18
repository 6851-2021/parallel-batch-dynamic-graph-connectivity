#pragma once

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "external/parlay/sequence.h"

#include "graph.h"
#include "parallel_euler_tour_tree/include/euler_tour_tree.hpp"
#include "utilities/include/hash.hpp"
// #include "utilities/include/utils.h"

// RESOLVED: I need to get cilk+ with gcc.
// RESOLVED: Move utilities into a proper include path
// RESOLVED: fix the include structure such that eulertree is in include path.
// RESOLVED: Write the make file for this to compiler.
// POSTPONED: Add the pbbs sequence library and related things here./
// TODO: A runner main file to do some benchmark things
// TODO: Get the above to compiler properly. 
// TODO: Update the documentation on all the of functions in headers


namespace detail {

    typedef int8_t Level;

    enum class EdgeType {
        // Edge is in the spanning forest of the graph.
        kNonTree,
        // Edge is not in the spanning forest of the graph.
        kTree,
    };

    // A struct that contains information about a particular edge.
    struct EdgeInfo {
        Level level;
        EdgeType type;
    };

}  // namespace detail

/** This class represents an undirected graph that can undergo efficient edge
 *  insertions, edge deletions, and connectivity queries. Multiple edges between
 *  a pair of vertices are supported.
 */

// TODO: this hackily allows me to use the sequence. Fix the namespacing later.

    
template<typename T>
using parlaysequence = parlay::sequence<T>;

namespace batchDynamicConnectivity {
    using UndirectedEdge = dynamicGraph::UndirectedEdge;
    using UndirectedEdgeHash = dynamicGraph::UndirectedEdgeHash;
    using BatchDynamicET = parallel_euler_tour_tree::EulerTourTree;

    class BatchDynamicConnectivity {
    public:
        /** Initializes an empty graph with a fixed number of vertices.
         *
         *  Efficiency: \f$ O(n \log n ) \f$ where \f$ n \f$ is the number of vertices
         *  in the graph.
         *
         *  @param[in] num_vertices Number of vertices in the graph.
         */
        explicit BatchDynamicConnectivity(int num_vertices);

        explicit BatchDynamicConnectivity(int num_vertices, const parlaysequence<UndirectedEdge> &se);

        /** Deallocates the data structure. */
        ~BatchDynamicConnectivity();

        /** The default constructor is invalid because the number of vertices in the
         *  graph must be known. */
        BatchDynamicConnectivity() = delete;

        /** Copy constructor not implemented. */
        BatchDynamicConnectivity(const BatchDynamicConnectivity &other) = delete;

        /** Copy assignment not implemented. */
        BatchDynamicConnectivity &operator=(const BatchDynamicConnectivity &other) = delete;

        /** Move constructor. */
        BatchDynamicConnectivity(BatchDynamicConnectivity &&other)

        noexcept;

        /** Move assignment not implemented. */
        BatchDynamicConnectivity &operator=(BatchDynamicConnectivity &&other)

        noexcept;

        // TODO:make the API use sequences for everything.
        /** Returns true if vertices \p u and \p v are connected in the graph.
         *
         *  Efficiency: logarithmic in the size of the graph.
         *
         *  @param[in] u Vertex.
         *  @param[in] v Vertex.
         *  @returns True if \p u and \p v are connected, false if they are not.
         */
        parlaysequence<char> BatchConnected(parlaysequence <std::pair<Vertex, Vertex>> suv) const;

        /** Returns true if edge \p edge is in the graph.
         *
         *  Efficiency: constant on average.
         *
         *  @param[in] edge Edge.
         *  @returns True if \p edge is in the graph, false if it is not.
         */
        bool HasEdge(const UndirectedEdge &edge) const;

        /** Returns the number of vertices in `v`'s connected component.
         *
         * Efficiency: logarithmic in the size of the graph.
         *
         * @param[in] v Vertex.
         * @returns The number of vertices in \p v's connected component.
         */
        int64_t GetSizeOfConnectedComponent(Vertex v) const;

        /** Adds an edge to the graph.
         *
         *  The edge must not already be in the graph and must not be a self-loop edge.
         *
         *  Efficiency: \f$ O\left( \log^2 n \right) \f$ amortized where \f$ n \f$ is
         *  the number of vertices in the graph.
         *
         *  @param[in] edge Edge to be added.
         */
        void BatchAddEdges(const parlaysequence <UndirectedEdge> &se);

        /** Deletes an edge from the graph.
         *
         *  An exception will be thrown if the edge is not in the graph.
         *
         *  Efficiency: \f$ O\left( \log^2 n \right) \f$ amortized where \f$ n \f$ is
         *  the number of vertices in the graph.
         *
         *  @param[in] edge Edge to be deleted.
         */
        void BatchDeleteEdges(const parlaysequence <UndirectedEdge> &se);


        parlaysequence <Vertex> BatchFindRepr(const parlaysequence <Vertex> &sv);

    private:

        const int64_t num_vertices_;
        const int64_t max_level_;

        // `spanning_forests_[i]` stores F_i, the spanning forest for the i-th
        // subgraph. In particular, `spanning_forests[0]` is a spanning forest for the
        // whole graph.

        // TODO: Turn this into a sequence
        // TODO: Turn dynamic forest to use parallel Euler tour trees. Convert to ParallelDynamicForest
        parlaysequence <BatchDynamicET*> parallel_spanning_forests_;

        // TODO: fix this so that the non_tree_adjacency_lists_ is now proper.

        // `adjacency_lists_by_level_[i][v]` contains the vertices connected to vertex
        // v by level-i non-tree edges.
        // TODO: make this concurrent map
        parlaysequence <parlaysequence <std::unordered_set < Vertex>>>
        non_tree_adjacency_lists_;

        // TODO: use a concurrent map here.
        // All edges in the graph.
        std::unordered_map <UndirectedEdge, detail::EdgeInfo, UndirectedEdgeHash>
                edges_;

        void AddNonTreeEdge(const UndirectedEdge &edge);

        void BatchAddNonTreeEdge(const parlaysequence <UndirectedEdge> &se);

        void AddTreeEdge(const UndirectedEdge &edge);

        void BatchAddTreeEdge(const parlaysequence <UndirectedEdge> &se);

        void AddEdgeToAdjacencyList(const UndirectedEdge &edge, detail::Level level);

        void BatchUpdateAdjacencyList(const parlaysequence <std::pair<UndirectedEdge, detail::Level>> &sel);

        void DeleteEdgeFromAdjacencyList(
                const UndirectedEdge &edge, detail::Level level);

        void BatchDeleteEdgesInAdjacencyList(const parlaysequence <std::pair<UndirectedEdge, detail::Level>> &sel);

        void ReplaceTreeEdge(const UndirectedEdge &edge, detail::Level level);
    };

    parlaysequence<std::unordered_set < Vertex>> generateVertexLayer(int numVertices, int max_level_){
        auto vtxLayer = parlaysequence<std::unordered_set < Vertex>>(numVertices);

        parallel_for(int i = 0; i < max_level_; i++){
            auto vtxset = std::unordered_set<Vertex> ();
            vtxLayer[i] = vtxset;
        }

        return vtxLayer;
    }

    BatchDynamicConnectivity::BatchDynamicConnectivity(int numVertices)
            : num_vertices_(numVertices), max_level_(log2(numVertices)) {

        parallel_spanning_forests_ = parlaysequence<BatchDynamicET*> (max_level_);
        
        parallel_for(int i = 0; i < max_level_; ++i){
            BatchDynamicET* ET = new BatchDynamicET{numVertices};
            parallel_spanning_forests_[i] = ET;
        }



        non_tree_adjacency_lists_ = parlaysequence<parlaysequence<std::unordered_set < Vertex>>>(max_level_);

        parallel_for(int i = 0; i < max_level_; ++i){
            auto vtxLayer = generateVertexLayer(numVertices, max_level_);
            non_tree_adjacency_lists_[i] = vtxLayer;
        } 

        edges_ = std::unordered_map <UndirectedEdge, detail::EdgeInfo, UndirectedEdgeHash>();
    }


}




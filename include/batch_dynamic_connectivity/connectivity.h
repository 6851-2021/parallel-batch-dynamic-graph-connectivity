#pragma once

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "external/parlay/sequence.h"
#include "external/parlay/primitives.h"
#include <map>
#include "graph.h"
#include "parallel_euler_tour_tree/include/euler_tour_tree.hpp"
#include "utilities/include/hash.hpp"

#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>

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

typedef int64_t Vertex;


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

    using treeSet = std::unordered_set<UndirectedEdge, UndirectedEdgeHash>;

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
        ~BatchDynamicConnectivity() = default;

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

        void printStructure();
        void printLevel(int8_t level);
        void printNonTreeEdges();
        void printNonTreeEdgesForLevel();

    private:

        const int64_t num_vertices_;


        // TODO: convert this to int8_t
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


        treeSet getSpanningTree(const parlaysequence <UndirectedEdge> &se);

        std::pair<int, int>* edgeBatchToPairArray(parlaysequence <UndirectedEdge> &se);

        template<typename Rank, typename Parent>
        treeSet constructTree(Rank& r, Parent& p, const parlaysequence<UndirectedEdge>& se);
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


    BatchDynamicConnectivity::BatchDynamicConnectivity(int numVertices, const parlaysequence <UndirectedEdge> &se)
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

        BatchAddEdges(se);
    }

    parlaysequence<char> BatchDynamicConnectivity::BatchConnected(parlaysequence <std::pair<Vertex, Vertex>> suv) const {
        parlaysequence<char> s(suv.size(), 0);

        BatchDynamicET* pMaxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];
        
        parallel_for(int i=0; i < suv.size(); i++){
            s[i] = pMaxLevelEulerTree->IsConnected(suv[i].first, suv[i].second);
        }
        return s;
    }

    parlaysequence <Vertex> BatchDynamicConnectivity::BatchFindRepr(const parlaysequence <Vertex> &sv) {
        auto pMaxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];

        return parlay::map(sv, [&](Vertex v){return (int64_t) pMaxLevelEulerTree->getRepresentative(v);});
    }

    template<typename Rank, typename Parent>
    treeSet BatchDynamicConnectivity::constructTree(Rank& r, Parent& p, const parlaysequence<UndirectedEdge>& se){
        boost::disjoint_sets<Rank, Parent> dsu(r, p);
        treeSet tree;
        
        // BUG POSSIBLE: check whether this causes a problem.
        for (auto v: se) {dsu.make_set(v.first); dsu.make_set(v.second);}

        for(auto v: se){
            Vertex first = v.first;
            Vertex second = v.second;
            //TODO is there a race condition here if we paralize this? How can we resolve that
            
            if(dsu.find_set(first) != dsu.find_set(second)){
                tree.insert(v);
                dsu.link(first, second);
            }
        }
        return tree;

    } 

    // TODO: add parallel DSU structure to implement this
    treeSet BatchDynamicConnectivity::getSpanningTree(const parlaysequence <UndirectedEdge> &se){
        //I am assuming the interface in
        //https://github.com/ParAlg/gbbs/blob/master/gbbs/union_find.h?fbclid=IwAR0U_Nbe1SpQF7mbmN0CEGLyF-5v362oy1q-9eQLvjQz916jhfTH69bMx9s
        // could be worth paralelizing this

        typedef std::map<Vertex, size_t> rank_t;
        typedef std::map<Vertex, Vertex> parent_t;

        rank_t rank_map;
        parent_t parent_map;

        boost::associative_property_map<rank_t>   rank_pmap(rank_map);
        boost::associative_property_map<parent_t> parent_pmap(parent_map);

        return constructTree(rank_pmap, parent_pmap, se);
    }
        
    std::pair<int, int>* BatchDynamicConnectivity::edgeBatchToPairArray(parlaysequence <UndirectedEdge> &se){
        std::pair<int, int>* array = new std::pair<int, int>[se.size()];
        parallel_for(int i=0; i<se.size(); i++){
            array[i].first =   se[i].first;
            array[i].second =   se[i].second;
        }
        return array;
    }


    void BatchDynamicConnectivity::BatchAddEdges(const parlaysequence <UndirectedEdge> &se) {
        auto maxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];

        parlaysequence <UndirectedEdge> auxiliaryEdges = parlay::map(se, [&](UndirectedEdge e) {
             return UndirectedEdge((Vertex) maxLevelEulerTree->getRepresentative(e.first),
                                   (Vertex) maxLevelEulerTree->getRepresentative(e.second));
        });


        auto tree = getSpanningTree(se);

        parlaysequence<UndirectedEdge> treeEdges;
        parlaysequence<UndirectedEdge> nonTreeEdges;

        // BUG POSSIBLE: we do a bad cast and make sure that max_level_ never overflows.
        
        parallel_for(int i=0; i < se.size(); i++){
            if(tree.count(se[i])){
                treeEdges.push_back(se[i]);
                detail::EdgeInfo ei = {(detail::Level) (max_level_ - 1), detail::EdgeType::kTree};
                edges_[se[i]] = ei;
            } else {
                nonTreeEdges.push_back(se[i]);
                detail::EdgeInfo ei = {(detail::Level) (max_level_ - 1), detail::EdgeType::kNonTree};
                edges_[se[i]] = ei;
            }
        }
        
        // add tree edges
        maxLevelEulerTree->BatchLink(edgeBatchToPairArray(treeEdges), treeEdges.size());

        // add to adjacancy list
        parallel_for(int i = 0; i < nonTreeEdges.size(); i++){
            non_tree_adjacency_lists_[max_level_ - 1][nonTreeEdges[i].first].insert(nonTreeEdges[i].second);
            non_tree_adjacency_lists_[max_level_ - 1][nonTreeEdges[i].second].insert(nonTreeEdges[i].first);
        }
    }


    // void BatchDynamicConnectivity::printLevel(int8_t level){
    //     std::cout << "### Printing level: " << level << "in the graph" << std::endl;
    //     for (auto e: edges_) {
    //         if (e.second.level == level && e.second.type == detail::EdgeType::kTree){
    //             std::cout << "<" << e.first.first << ", " << e.first.second << ">" << " : tree" << std::endl;
    //         }
    //         if (e.second.level == level && e.second.type == detail::EdgeType::kNonTree){
    //             std::cout << "<" << e.first << ", " << e.second << ">" << " : not tree \n" << std::endl;
    //         }
    //     }
    //     std::cout << "### end of level" << level << "\n";
    // }
    // void BatchDynamicConnectivity::printStructure(){
    //     for (int i = 0; i < max_level_; i++){
    //         printLevel(i);
    //     }
    // }
    // void BatchDynamicConnectivity::printNonTreeEdges(){
    //     std::cout << "#####:  printing non tree edges\n";
    //     for (auto e: edges_) {
    //         if (e.second.type == detail::EdgeType::kTree){
    //             std::cout << "<" << e.first.first << ", " << e.first.second << ">" << std::endl;
    //         }
    //     }
    //     std::cout << "### end of tree edges" << "\n";
    // }
    // void BatchDynamicConnectivity::printNonTreeEdgesForLevel(){
    //     std::cout << "#####:  printing non tree edges\n";
    //     for (auto e: edges_) {
    //         if (e.second.type == detail::EdgeType::kNonTree){
    //             std::cout << "<" << e.first.first << ", " << e.first.second << ">" << std::endl;
    //         }
    //     }
    //     std::cout << "### end of non tree edges" << "\n";
    // }
}




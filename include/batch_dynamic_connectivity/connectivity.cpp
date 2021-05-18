//
// Created by sualehasif on 5/18/21.
//

#include "connectivity.h"

// RESOLVED: add the parallel constructor for edges
// TODO: figure out a way to do BatchQuery for representatives.

namespace batchDynamicConnectivity {

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

    // sequence <Vertex> BatchDynamicConnectivity::BatchFindRepr(const sequence <Vertex> &sv) {
    //     auto pMaxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];

    //     return sv.map([&maxLevelEulerTree](vertex v){return maxLevelEulerTree.getRepresentative(v);});

    //     // return trees.at(maxLevel).BatchFindRepr(suv);
    // }

    // // TODO: add parallel DSU structure to implement this
    // std::unordered_set<UndirectedEdge> getSpanningTreeIndices(sequence <UndirectedEdge> &se){
    // //returns an unordered set of the indices in the sequence of the
    // //edges forming a spanning tree
    //     throw 10; // implement me
    // }
    // std::pair<int, int>* edgeBatchToPairArray(sequence <UndirectedEdge> &se){
    //     std::pair<int, int>* array = new std::pair<int, int>[se.size()];
    //     parallel_for(int i=0; i<se.size(); i++){
    //         array[i].first =   se[i].first;
    //         array[i].second =   se[i].second;
    //     }
    //     return array;
    // }

    // void BatchDynamicConnectivity::BatchAddEdges(const sequence <UndirectedEdge> &se) {
    //     auto maxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];
    //     sequence <UndirectedEdge> auxiliaryEdges = se.map([](UndirectedEdge e) {
    //          return UndirectedEdge(maxLevelEulerTree.getRepresentative(e.first),
    //                                maxLevelEulerTree.getRepresentative(e.second))
    //     });
    //     auto treeIndices = getSpanningTreeIndices(sequence <UndirectedEdge> &se);
    //     sequence<UndirectedEdge> treeEdges;
    //     sequence<UndirectedEdge> nonTreeEdges;
    //     parallel_for(int i=0; i < se.size(); i++){
    //         if(treeIndices.contains(i))
    //             treeEdges.push_back(se[i]);
    //         else
    //             nonTreeEdges.push_back(se[i]);
    //     }

    //     //set level
    //     parallel_for(int i=0; i < se.size(); i++){
    //         edges_[se[i]] = max_level_ - 1;
    //     }
        
    //     // add tree edges
    //     maxLevelEulerTree.BatchLink(edgeBatchToPairArray(treeEdges));

    //     // add to adjacancy list
    //     parallel_for(int i = 0; i < nonTreeEdges.size(); i++){
    //         non_tree_adjacency_lists_[max_level_ - 1]
    //     }
        


    //     // trees.at(maxLevel).BatchInsert(treeEdges);
    //     // graph.insert(maxLevel, nonTreeEdges)
    // }

    // void BatchDynamicConnectivity::BatchDeleteEdges(const sequence <UndirectedEdge> &se) {
    //     //TODO: split se into tree and non tree edges
    //     // graph.BatchDelete(nonTreeEdges);
    //     // int64_t minLevel = treeEdges.map([](UndirectedEdge e) { return e.level }).min()
    //     // for (int i = minLevel; i <= maxLevel; i++) {
    //     //     auto tree = trees.at(i);
    //     //     // If tree does not support deleting non existent edges do the following first:
    //     //     //auto toDelete = treeEdges.filter([int i](UndirectedEdge e){return e.level <= i});
    //     //     tree.BatchDeleteEdges(treeEdges);
    //     // }
    //     // sequence <vertex> components = treeEdges.map([](UndirectedEdge e) { return e.first }); //C in paper
    //     // components.merge(treeEdges.map([](UndirectedEdge e) { return e.second }));
    //     // components = sort(components); // I think I saw radix sort somewhere in seq
    //     // components = removeDuplicatesFromSortedSeq(components);
    //     // sequence <UndirectedEdge> promotedEdges; // empty sequence, S in paper
    //     // for (int i = minLevel; i <= maxLevel; i++) {
    //     //     compononents = parallelLevelSearch(components, promotedEdges, i);
    //     // }

    //     return;
    // }

    // sequence <vertex>
    // parallelLevelSearch(sequence <vertex> components, sequence <UndirectedEdge> promotedEdges, int i) {
    // tree = trees.at(i);
    // tree.BatchAddEdges(promotedEdges);
    // sequence < vertex > componentsToConsider = components.filter(
    // [int i, EulerTourTree
    // tree](auto
    // v){ return tree.subtreeSize(v) } <= 1 << (i - 1));
    // sequence < vertex > largeComponents = components.filter([int i, EulerTourTree
    // tree](auto
    // v){ return tree.subtreeSize(v) } > 1 << (i - 1));
    // sequence <UndirectedEdge> R;
    // while (!componentsToConsider.empty()) {
    //     for (vertex v: componentsToConsider) {
    //         //Move all the edges of small components down a level
    //         for (UndirectedEdge e : tree.edges(v)) {  // edges loops over the edges of the tree
    //             if e.level == i{
    //                         e.level = i - 1; // this is probably not nescissary as trees do not internally use edges
    //                         trees.at(i-1).AddEdge(e);
    //                 }
    //         }
    //         R.merge(componentSearch(i, v));
    //     }
    //     //code below could be refactored into a helper with insert
    //     sequence <UndirectedEdge> auxiliaryEdges = R.map([](UndirectedEdge e) {
    //         return UndirectedEdge(trees.at(maxLevel).FindRepr(e.first),
    //                               trees.at(maxLevel).FindRepr(e.second))
    //     });
    //     auto treeEdgeIndices = FindSpanningTree(auxiliaryEdges)
    //     auto treeEdges = getElements(treeEdgeIndices, R) // gets the indexed
    //     auto nonTreeEdges = getOtherElements(treeEdgeIndices, R) //gets non indexed elements
    //     tree.BatchInsert(treeEdges);
    //     for (auto e : treeEdges) {
    //         G.remove(e); // should store index of edge in adjacency in edge object for efficiency
    //     }
    //     promotedEdges.merge(treeEdges);
    //     componentsToConsider = componentsToConsider.map([EulerTourTree
    //     tree](Vertex
    //     v){
    //         return tree.Rep(v)
    //     });
    //     componentsToConsider = sort(componentsToConsider);
    //     componentsToConsider = removeDuplicatesFromSortedSeq(componentsToConsider);
    //     sequence < vertex > largeComponents.merge(componentsToConsider.filter([int i, EulerTourTree
    //     tree](auto
    //     v){ return tree.subtreeSize(v) } > 1 << (i - 1));)
    //     sequence < vertex > componentsToConsider = componentsToConsider.filter(
    //     [int i, EulerTourTree
    //     tree](auto
    //     v){ return tree.subtreeSize(v) } <= 1 << (i - 1));
    // }
    // return largeComponents;

    //     sequence<vertex> v;
    //     return v;
    // }

    // sequence <UndirectedEdge> componentSearch(int i, Vertex v) {
    //     sequence <UndirectedEdge> R;
    //     for (auto u : v.subtree()) {
    //         for (auto e : graph.getEdges(u, i)) { // gets non tree edge incident to u at level i
    //             if (tree.at(i).Repr(e.second) != v) {
    //                 R.pushBack(v)
    //                 return R
    //             }
    //         }
    //     }
    // }

    int64_t log2(int64_t n) {
        int64_t lg = 0;
        int64_t power = 1;
        while (power < n) {
            power = power << 1;
            lg += 1;
        }
        return lg;
    }
}


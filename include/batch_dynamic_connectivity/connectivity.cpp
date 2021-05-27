//
// Created by sualehasif on 5/18/21.
//

#include "connectivity.h"

// RESOLVED: add the parallel constructor for edges
// TODO: figure out a way to do BatchQuery for representatives.
// TODO: something about MST or DSU
// TODO: figure our a way to loop over subtree (preferably edges)
// TODO: get a remove duplicates function from parlay. 
// TODO: 

namespace batchDynamicConnectivity {


    void BatchDynamicConnectivity::BatchDeleteEdges(const sequence <UndirectedEdge> &se) {
        //deletes a badge of edges, must all have been inserted before


        //TODO: split se into tree and non tree edges
        // delete edges from adjacency list
        sequence<UndirectedEdge> treeEdges; //a list of the tree edges in the batch
        auto min_tree_edge_level = max_level_;

        //the loop deletes the edges from the adjacency list, and computes the above
        parallel_for(int i=0; i< se.size(); i++){
            auto level = edges_[se[i]];
            auto u = se[i].first;
            auto v = se[i].second;
            if (non_tree_adjacency_lists_[level][u].contains(v)){
                non_tree_adjacency_lists_[level][u].remove(v);
                non_tree_adjacency_lists_[level][v].remove(u);
            } else {
                treeEdges.push_back(se[i]);
                if (level < min_tree_edge_level ){
                    min_tree_edge_level = level;
                }
            }
        }
        //deletes the edges from each forest level
        for(int level = min_tree_edge_level; level < max_level_; i++){
            auto levelEulerTree = parallel_spanning_forests_[level];
            toDelete = treeEdges.filter([level, &edges_](UndirectedEdge e){return edges_[e] <= level});
            levelEulerTree.BatchCut(edgeBatchToPairArray(toDelete), toDelete.size());    
        }
        //components should end up as a list of represenytatives of the endpoints of eges with no duplicates
        //It starts as just the endpoints, then when trying to find replacement in each level
        //we collapse into representatives in that level

        auto components = treeEdges.map([](UndirectedEdge e) { return e.first });
        //TODO: merge properly
        components.merge(treeEdges.map([](UndirectedEdge e) { return e.second }));
        components = removeDuplicates(components);
        

        sequence <UndirectedEdge> promotedEdges; //edges promoted so far.
        for(int i = minLevel; i <= maxLevel; i++){
            compononents = parallelLevelSearch(components, promotedEdges, i);
        }
        // empty sequence, S in paper
        // int64_t minLevel = treeEdges.map([](UndirectedEdge e) { return e.level }).min()
        // for (int i = minLevel; i <= maxLevel; i++) {
        //     auto tree = trees.at(i);
        //     // If tree does not support deleting non existent edges do the following first:
        //     //auto toDelete = treeEdges.filter([int i](UndirectedEdge e){return e.level <= i});
        //     tree.BatchDeleteEdges(treeEdges);
        // }
        // sequence <vertex> components = treeEdges.map([](UndirectedEdge e) { return e.first }); //C in paper
        // components.merge(treeEdges.map([](UndirectedEdge e) { return e.second }));
        // components = sort(components); // I think I saw radix sort somewhere in seq
        // components = removeDuplicatesFromSortedSeq(components);
        // sequence <UndirectedEdge> promotedEdges; // empty sequence, S in paper
        // for (int i = minLevel; i <= maxLevel; i++) {
        //     compononents = parallelLevelSearch(components, promotedEdges, i);
        // }
        return;
    }

sequence<Vertex> BatchDynamicConnectivity::BatchFindRepr(const sequence<Vertex>& sv){
    return trees.at(maxLevel).BatchFindRepr(suv);
}


void BatchDynamicConnectivity::BatchAddEdges(const sequence<UndirectedEdge>& se){
    //The code for finding spanning tree can be a helper for this and del
    //set the levels of edges to be maxLevel
    sequence<UndirectedEdge> auxiliaryEdges = se.map([](UndirectedEdge e){
        return UndirectedEdge(trees.at(maxLevel).FindRepr(e.first),
         trees.at(maxLevel).FindRepr(e.second))});
          //represent the components connected by the given edges
    auto treeEdgeIndices = FindSpanningTree(auxiliaryEdges)
    auto treeEdges = getElements(treeEdgeIndices, se) // gets the indexed
    auto nonTreeEdges = getOtherElements(treeEdgeIndices, se) //gets non indexed elements
    trees.at(maxLevel).BatchInsert(treeEdges);
    graph.insert(maxLevel, nonTreeEdges)
}

        sequence <UndirectedEdge> R;

        while (componentsToConsider.size() != 0) {
            sequence <vertex> edgesToDropLevel;
            for (vertex v: componentsToConsider) {
                //Move all the edges of small components down a level
                for (UndirectedEdge e : tree.edges(v)) {  // edges loops over the edges of the tree
                    if edges_[e] == level{
                            edges_[e] = level - 1;
                            edgesToDropLevel.push_back(e);
                    }
                }
                R.push_back(componentSearch(i, v));
            }
            parallel_spanning_forests_[level - 1].BatchLink(edgeBatchToPairArray(edgesToDropLevel), edgesToDropLevel.size());


            sequence <UndirectedEdge> auxiliaryEdges = R.map([](UndirectedEdge e) {
                return UndirectedEdge(maxLevelEulerTree.getRepresentative(e.first),
                                       maxLevelEulerTree.getRepresentative(e.second))
            });
            auto promIndices = getSpanningTreeIndices(sequence <UndirectedEdge> &se);
            sequence<UndirectedEdge> newPromotedEdges;
            sequence<UndirectedEdge> notPromotedEdges;
            parallel_for(int i=0; i < se.size(); i++){
                if(promIndices.contains(i))
                    newPromotedEdges.push_back(se[i]);
                else
                    notPromotedEdges.push_back(se[i]);
            }
        
            levelEulerTree.BatchLink(edgeBatchToPairArray(newPromotedEdges), newPromotedEdges.size());
            parallel_for (int i = 0; i < newPromotedEdges.size(); i++) {
                level = edges_[newPromotedEdges[i]];
                auto u = newPromotedEdges[i].first;
                auto v = newPromotedEdges[i].second;
                non_tree_adjacency_lists_[level][u].remove(v);
                non_tree_adjacency_lists_[level][v].remove(u);
                promotedEdges.push_back(newPromotedEdges[i]);
            }


    UndirectedEdge componentSearch(int level, Vertex v) {
        //looks for an out edge in a component that connests it with a different component
        auto levelEulerTree = parallel_spanning_forests_[level];        
        //TODO
        for(auto u: levelEulerTree.subtree(v)){
            for(auto w : non_tree_adjacency_lists_[level][u]){
                if(levelEulerTree.getRepresentative(u) != levelEulerTree.getRepresentative(v)){
                    return UndirectedEdge(u, w)
                }
            }
            R.merge(componentSearch(i, v));
        }
    //code below could be refactored into a helper with insert
    sequence<UndirectedEdge> auxiliaryEdges = R.map([](UndirectedEdge e){
        return UndirectedEdge(trees.at(maxLevel).FindRepr(e.first),
         trees.at(maxLevel).FindRepr(e.second))});
    auto treeEdgeIndices = FindSpanningTree(auxiliaryEdges)
    auto treeEdges = getElements(treeEdgeIndices, R) // gets the indexed
    auto nonTreeEdges = getOtherElements(treeEdgeIndices, R) //gets non indexed elements
    tree.BatchInsert(treeEdges);
    for (auto e : treeEdges){
        G.remove(e); // should store index of edge in adjacency in edge object for efficiency
    }
    promotedEdges.merge(treeEdges);
    componentsToConsider = componentsToConsider.map([EulerTourTree tree](Vertex v){
        return tree.Rep(v)});
    componentsToConsider = sort(componentsToConsider);
    componentsToConsider = removeDuplicatesFromSortedSeq(componentsToConsider);
    sequence<vertex> largeComponents.merge(componentsToConsider.filter([int i, EulerTourTree tree](auto v){return tree.subtreeSize(v)} > 1 << (i - 1));)
    sequence<vertex> componentsToConsider = componentsToConsider.filter(
        [int i, EulerTourTree tree](auto v){return tree.subtreeSize(v)} <= 1 << (i - 1));
    }
}

int64_t log2(int64_t n){
    int64_t lg = 0;
    int64_t power = 1;
    while (power < n){
        power = power << 1;
        lg += 1;
    }
    return lg;
}
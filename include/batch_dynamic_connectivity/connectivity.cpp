//
// Created by sualehasif on 5/18/21.
//

#include "connectivity.h"

// TODO: add the parallel constructor for edges
// TODO:

namespace batchDynamicConnectivity {
    BatchDynamicConnectivity::BatchDynamicConnectivity(int64_t numVertices)
            : num_vertices_(numVertices), max_level_(log2(numVertices)) {
        
        // parallel initialization for the spanning forests using parallel euler tour trees.
        parallel_spanning_forests_ = sequence<BatchDynamicET*>(max_level_, 
            [&numVertices](){
                auto ET = new BatchDynamicET{numVertices};
                return &ET;
            }
        );

        non_tree_adjacency_lists_ = sequence<sequence<std::unordered_set < Vertex>>> (max_level_, 
            [&numVertices](){ 
                auto vertex_sequence = new sequence<std::unordered_set < Vertex>> (numVertices, 
                    [](){
                        return std::unordered_set< Vertex>();
                    }
                );
                return vertex_sequence;
            }
        );

        edges_ = std::unordered_map <UndirectedEdge, detail::EdgeInfo, UndirectedEdgeHash>();
    }

    BatchDynamicConnectivity::BatchDynamicConnectivity(int64_t numVertices, const sequence <UndirectedEdge> &se)
            : num_vertices_(numVertices), max_level_(log2(numVertices)) {
        
        // parallel initialization for the spanning forests using parallel euler tour trees.
        parallel_spanning_forests_ = sequence<BatchDynamicET*>(max_level_, 
            [&numVertices](){
                auto ET = new BatchDynamicET{numVertices};
                return &ET;
            }
        );


        non_tree_adjacency_lists_ = sequence<sequence<std::unordered_set < Vertex>>> (max_level_, 
            [&numVertices](){ 
                auto vertex_sequence = new sequence<std::unordered_set < Vertex>> (numVertices, 
                    [](){
                        return std::unordered_set< Vertex>();
                    }
                );
                return vertex_sequence;
            }
        );

        edges_ = std::unordered_map <UndirectedEdge, detail::EdgeInfo, UndirectedEdgeHash>();

        BatchAddEdges(se);
    }

    sequence<char> BatchDynamicConnectivity::BatchConnected(sequence <std::pair<Vertex, Vertex>> suv) const {
        sequence<char> s(suv.size(), 0);

        BatchDynamicET* pMaxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];
        
        parallel_for(int i=0; i < suv.size(); i++){
            s[i] = pMaxLevelEulerTree->IsConnected(std::get<0>(suv[i]), std::get<1>(suv[i]));
        }
        return s;
    }

    sequence <Vertex> BatchDynamicConnectivity::BatchFindRepr(const sequence <Vertex> &sv) {
        auto pMaxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];

        return sv.map([&maxLevelEulerTree](vertex v){return maxLevelEulerTree.getRepresentative(v);});

        // return trees.at(maxLevel).BatchFindRepr(suv);
    }

    // TODO: add parallel DSU structure to implement this
    std::unordered_set<UndirectedEdge> getSpanningTreeIndices(sequence <UndirectedEdge> &se){
        //I am assuming the interface in
        //https://github.com/ParAlg/gbbs/blob/master/gbbs/union_find.h?fbclid=IwAR0U_Nbe1SpQF7mbmN0CEGLyF-5v362oy1q-9eQLvjQz916jhfTH69bMx9s
        // could be worth paralelizing this
        UnionFind unionFind (num_vertices_);
        std::unordered_set<UndirectedEdge> tree;
        for(int i=0; i<se.size(); i++){
            vertex first = se[i].first;
            vertex second = se[i].second;
            //TODO is there a race condition here if we paralize this? How can we resolve that
            if(unionFind.find(first) != unionFind.find(first)){
                tree.insert(se[i]);
                unionFind.link(first, second);
            }
        }
        return tree;
    }
    std::pair<int, int>* edgeBatchToPairArray(sequence <UndirectedEdge> &se){
        std::pair<int, int>* array = new std::pair<int, int>[se.size()];
        parallel_for(int i=0; i<se.size(); i++){
            array[i].first =   se[i].first;
            array[i].second =   se[i].second;
        }
        return array;
    }

    void BatchDynamicConnectivity::BatchAddEdges(const sequence <UndirectedEdge> &se) {
        auto maxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];
        sequence <UndirectedEdge> auxiliaryEdges = se.map([](UndirectedEdge e) {
             return UndirectedEdge(maxLevelEulerTree.getRepresentative(e.first),
                                   maxLevelEulerTree.getRepresentative(e.second))
        });
        auto treeIndices = getSpanningTreeIndices(sequence <UndirectedEdge> &se);
        sequence<UndirectedEdge> treeEdges;
        sequence<UndirectedEdge> nonTreeEdges;
        parallel_for(int i=0; i < se.size(); i++){
            if(treeIndices.contains(i))
                treeEdges.push_back(se[i]);
            else
                nonTreeEdges.push_back(se[i]);
        }

        //set level
        parallel_for(int i=0; i < se.size(); i++){
            edges_[se[i]] = max_level_ - 1;
        }
        
        // add tree edges
        maxLevelEulerTree.BatchLink(edgeBatchToPairArray(treeEdges), treeEdges.size());

        // add to adjacancy list
        parallel_for(int i = 0; i < nonTreeEdges.size(); i++){
            non_tree_adjacency_lists_[max_level_ - 1][nonTreeEdges[i].first].add(nonTreeEdges[i].second);
            non_tree_adjacency_lists_[max_level_ - 1][nonTreeEdges[i].second].add(nonTreeEdges[i].first);
        }
    }

//TODO implement semisort or any sort
void removeDuplicates(sequence<Vertex>& seq){
    semisort(sequence<Vertex> seq);
    sequence<Vertex> newSeq;
    newSeq.push_back(seq[0]);
    parallel_for(int i =1; i < seq.size(); i++){
        if(seq[i] != seq[i-1]){
            newSeq.push_back(seq[i]);
        }
    }
    return newSeq;
}

    void BatchDynamicConnectivity::BatchDeleteEdges(const sequence <UndirectedEdge> &se) {
        //TODO: split se into tree and non tree edges
        // delete edges from adjacency list
        sequence<UndirectedEdge> treeEdges;
        auto min_tree_edge_level = max_level_;
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
        for(int level = min_tree_edge_level; level < max_level_; i++){
            auto levelEulerTree = parallel_spanning_forests_[level];
            toDelete = treeEdges.filter([level, &edges_](UndirectedEdge e){return edges_[e] <= level});
            levelEulerTree.BatchCut(edgeBatchToPairArray(toDelete), toDelete.size());    
        }
        auto components = treeEdges.map([](UndirectedEdge e) { return e.first });
        //TODO: merge properly
        components.merge(treeEdges.map([](UndirectedEdge e) { return e.second }));
        components = removeDuplicates(components);
        sequence <UndirectedEdge> promotedEdges;
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

    sequence <vertex>
    parallelLevelSearch(sequence <vertex>& components, sequence <UndirectedEdge>& promotedEdges, int level) {
        auto levelEulerTree = parallel_spanning_forests_[level];
        levelEulerTree.BatchLink(promotedEdges);
        levelEulerTree.BatchAddEdges(edgeBatchToPairArray(promotedEdges), promotedEdges.size());
        components = components.map([&levelEulerTree](auto v){levelEulerTree.getRepresentative(v)});
        components = removeDuplicates(components);

        sequence < vertex > componentsToConsider; 
        sequence < vertex > largeComponents;
        parallel_for(int i = 0; i < components.size(); i++){
            if(levelEulerTree.componentSize(components[i]) <= 1 << (level - 1) ){
                componentsToConsider.push_back(components[i]);
            } else{
                largeComponents.push_back(components[i]);
            }
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

            componentsToConsider = 
            componentsToConsider.map(
                [levelEulerTree](Vertex v){return levelEulerTree.getRepresentative(v)}
                );
        
            componentsToConsider = removeDuplicates(componentsToConsider);

            sequence <vertex> newComponentsToConsider;
            parallel_for (int i = 0; i < componentsToConsider.size(); i++){
                if( levelEulerTree.componentSize(componentsToConsider[i])  <= 1 << (level - 1)){
                    newComponentsToConsider.push_back(componentsToConsider[i]);
                } else {
                    largeComponents.push_back(componentsToConsider[i]);
                }
            }
            componentsToConsider = newComponentsToConsider;
        }
        return largeComponents;
    }

    UndirectedEdge componentSearch(int level, Vertex v) {

        auto levelEulerTree = parallel_spanning_forests_[level];        
        //TODO
        for(auto u: levelEulerTree.subtree(v)){
            for(auto w : non_tree_adjacency_lists_[level][u]){
                if(levelEulerTree.getRepresentative(u) != levelEulerTree.getRepresentative(v)){
                    return UndirectedEdge(u, w)
                }
            }
        }
    }

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


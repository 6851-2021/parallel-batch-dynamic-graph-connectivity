#include "connectivity.hpp"

// Below is C-like pseudocode to be converted into C++

BatchDynamicConnectivity::BatchDynamicConnectivity(int64_t numVertices){
    int64_t maxLevel = log2(numVertices);
    Graph graph =  Graph();

    //populate the sequence below with max_height trees with num vertices nodes
    sequence<EulerTourTree::EulerTourTree> trees;
    //methods needed for Tree:
    // -iterate over edges in a subtree  \\ .edges() \\ 
    // -iterate over nodes in a subtree \\ subtree()  \\

    //methods needed for Graph():
    // - iterate over edges incident tp u at level i \\ getEdges(u, i) \\
    // - delete a batch edges   \\BatchDelete(nonTreeEdges) \\
    // - insert a btach of edges at a level \\insert(level, nonTreeEdges)\\ 
}

sequence<char> BatchDynamicConnectivity::BatchConnected(sequence<std::pair<Vertex, Vertex>> suv){
    return trees.at(maxLevel).BatchConnected(suv);
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




void BatchDynamicConnectivity::BatchDeleteEdges(const sequence<UndirectedEdge>& se){
    //TODO: split se into tree and non tree edges
    graph.BatchDelete(nonTreeEdges);
    int64_t minLevel = treeEdges.map([](UndirectedEdge e){return e.level}).min()
    for(int i = minLevel; i <=  maxLevel; i++ ){
        auto tree = trees.at(i);
        // If tree does not support deleting non existent edges do the following first:
        //auto toDelete = treeEdges.filter([int i](UndirectedEdge e){return e.level <= i});
        tree.BatchDeleteEdges(treeEdges); 
    }
    sequence<vertex> components = treeEdges.map([](UndirectedEdge e){return e.first}); //C in paper
    components.merge(treeEdges.map([](UndirectedEdge e){return e.second}));
    components = sort(components); // I think I saw radix sort somewhere in seq
    components = removeDuplicatesFromSortedSeq(components);
    sequence<UndirectedEdge> promotedEdges; // empty sequence, S in paper
    for(int i = minLevel; i <=  maxLevel; i++ ){
        compononents = parallelLevelSearch(components, promotedEdges, i);
    }
}

sequence<vertex> parallelLevelSearch(sequence<vertex> components, sequence<UndirectedEdge> promotedEdges, int i){
    tree = trees.at(i);
    tree.BatchAddEdges(promotedEdges);
    sequence<vertex> componentsToConsider = components.filter(
        [int i, EulerTourTree tree](auto v){return tree.subtreeSize(v)} <= 1 << (i - 1));
    sequence<vertex> largeComponents = components.filter([int i, EulerTourTree tree](auto v){return tree.subtreeSize(v)} > 1 << (i - 1));
    sequence<UndirectedEdge> R;
    while (! componentsToConsider.empty()){
        for (vertex v: componentsToConsider){
            //Move all the edges of small components down a level
            for (UndirectedEdge e : tree.edges(v)){  // edges loops over the edges of the tree
                if e.level == i{
                    e.level = i - 1; // this is probably not nescissary as trees do not internally use edges
                    trees.at(i-1).AddEdge(e);
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
    return largeComponents;
}


sequence<UndirectedEdge> componentSearch(int i, Vertex v){
    sequence<UndirectedEdge> R;
    for (auto u : v.subtree()){
        for(auto e : graph.getEdges(u, i)){ // gets non tree edge incident to u at level i
            if(tree.at(i).Repr(e.second) != v){
                R.pushBack(v)
                return R
            }
        }
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
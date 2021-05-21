//
// Created by sualehasif on 5/18/21.
//

#include "batch_dynamic_connectivity/connectivity.h"

#include <cassert>
#include <random>
#include <utility>
#include <iostream>
// #include <gtest/gtest.h>
// #include "catch.hpp"

using BatchDynamicConnectivity = batchDynamicConnectivity::BatchDynamicConnectivity;

using UndirectedEdge = dynamicGraph::UndirectedEdge;
using UndirectedEdgeHash = dynamicGraph::UndirectedEdgeHash;
using BatchDynamicET = parallel_euler_tour_tree::EulerTourTree;

using treeSet = std::unordered_set<UndirectedEdge, UndirectedEdgeHash>;



bool Test1(){
    parlaysequence<UndirectedEdge> edges;
    edges.push_back(UndirectedEdge(0, 1));
    edges.push_back(UndirectedEdge(1, 2));
    edges.push_back(UndirectedEdge(3, 4));
    
    BatchDynamicConnectivity x (5, edges);
    
    parlaysequence<std::pair<Vertex, Vertex>> queries;
    parlaysequence<char> expectedOut;
    
    queries.push_back(std::make_pair(0,1));

    expectedOut.push_back(true);

    auto result = x.BatchConnected(queries);
    for(int i=0; i < queries.size(); i++){
        if (result[i] != expectedOut[i])
            return false;
    }
    return true;
    
}


bool Test2(){
    parlaysequence<UndirectedEdge> edges;
    edges.push_back(UndirectedEdge(0, 1));
    edges.push_back(UndirectedEdge(1, 2));
    edges.push_back(UndirectedEdge(3, 4));
    
    BatchDynamicConnectivity x (5, edges);
    
    parlaysequence<std::pair<Vertex, Vertex>> queries;
    parlaysequence<char> expectedOut;
    
    queries.push_back(std::make_pair(0,1));
    queries.push_back(std::make_pair(1,2));
    queries.push_back(std::make_pair(0,2));
    queries.push_back(std::make_pair(3,4));
    queries.push_back(std::make_pair(0,3));
    queries.push_back(std::make_pair(0,4));
    queries.push_back(std::make_pair(1,3));
    queries.push_back(std::make_pair(1,4));

    expectedOut.push_back(true);
    expectedOut.push_back(true);
    expectedOut.push_back(true);
    expectedOut.push_back(true);
    expectedOut.push_back(false);
    expectedOut.push_back(false);
    expectedOut.push_back(false);
    expectedOut.push_back(false);

    auto result = x.BatchConnected(queries);
    for(int i=0; i < queries.size(); i++){
        if (result[i] != expectedOut[i])
            return false;
    }
    return true;
    
}

bool Test3(){
    parlaysequence<UndirectedEdge> edges;
     //connected component with first four edges
    for(int i = 0; i < 5; i++){
        for(int j = 0; j < i; j++){
            edges.push_back(UndirectedEdge(i, j));
        }
    }
    edges.push_back(UndirectedEdge(5, 6));
    edges.push_back(UndirectedEdge(6, 7));
    
    BatchDynamicConnectivity x (10, edges);
    
    parlaysequence<std::pair<Vertex, Vertex>> queries;
    parlaysequence<char> expectedOut;
    
    for(long i = 0; i < 10; i++){
        for(long j = 0; j < i; j++){
            queries.push_back(std::make_pair(i, j));
            expectedOut.push_back((i == j) || ((i < 5) && (j < 5)) || ((5 <= i) && (5 <= j) && (i < 8) && (j < 8)));//|| ((i < 10) && (j < 10))
        }
    }
    
    auto result = x.BatchConnected(queries);
    for(int i=0; i < queries.size(); i++){
        if (result[i] != expectedOut[i])
            return false;
    }
    return true;
}


bool Test6(){
    parlaysequence<UndirectedEdge> edges;
     //connected component with first four edges
    for(int i = 0; i < 5; i++){
        for(int j = 0; j < i; j++){
            edges.push_back(UndirectedEdge(i, j));
        }
    }
    edges.push_back(UndirectedEdge(5, 6));
    edges.push_back(UndirectedEdge(6, 7));
    edges.push_back(UndirectedEdge(7, 8));
    edges.push_back(UndirectedEdge(8, 9));
    edges.push_back(UndirectedEdge(9, 10));
    
    BatchDynamicConnectivity x (11, edges);
    
    parlaysequence<std::pair<Vertex, Vertex>> deletes;
    deletes.push_back(0, 1);
    deletes.push_back(1, 2);
    deletes.push_back(3, 4);
    deletes.push_back(4, 1);
    deletes.push_back(7, 8);
    
    x.BatchDeleteEdges(deletes);

    parlaysequence<std::pair<Vertex, Vertex>> queries;
    parlaysequence<char> expectedOut;
    
    for(long i = 0; i < 10; i++){
        for(long j = 0; j < i; j++){
            queries.push_back(std::make_pair(i, j));
            expectedOut.push_back((i == j) || ((i < 5) && (j < 5)) || ((5 <= i) && (5 <= j) && (i < 8) && (j < 8)) || ((8 <= i) && (8 <= j) && (i < 11) && (j < 11)) );
        }
    }
    
    auto result = x.BatchConnected(queries);
    for(int i=0; i < queries.size(); i++){
        if (result[i] != expectedOut[i])
            return false;
    }
    return true;
}

int main(int argc, char **argv) {
    // ::testing::InitGoogleTest(&argc, argv);
    // return RUN_ALL_TESTS();
    if (Test1()){
        std::cout <<"Test 1: passed\n";
    } else{
        std::cout << "Test 1: failed\n";
    }

    if (Test2()){
        std::cout <<"Test 2: passed\n";
    } else{
        std::cout << "Test 2: failed\n";
    }

    if (Test3()){
        std::cout <<"Test 3: passed\n";
    } else{
        std::cout << "Test 3: failed\n";
    }


}
        
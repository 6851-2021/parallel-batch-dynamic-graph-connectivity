//
// Created by sualehasif on 5/18/21.
//

#include "batch_dynamic_connectivity/connectivity.h"

#include <cassert>
#include <random>
#include <utility>
#include <iostream>
#include <gtest/gtest.h>

using BatchDynamicConnectivity = batchDynamicConnectivity::BatchDynamicConnectivity;


TEST(BatchDynamicConnectivity, BatchConnected) {
    sequence<UndirectedEdge> edges;
    //connected component with first four edges
    for(int i = 0; i < 5; i++){
        for(int j = 0; j < i; j++){
            edges.push_back(UndirectedEdge(i, j));
        }
    }
    edges.push_back(UndirectedEdge(5, 6));
    edges.push_back(UndirectedEdge(6, 7));
    BatchDynamicConnectivity x (10, edges);
    sequence<std::pair<Vertex, Vertex>>> querries;
    sequence<char> expectedOut;
    for(int i =0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            querries.push_back(std::pair<Vertex, Vertex>>(i, j));
            expectedOut.push_back((i == j) || (i < 5) && (j < 5) || (i < 8) && (j < 8) || (i < 10) && (j < 10));
        }
    }
    EXPECT_EQ(BatchConnected(querries), expectedOut );    
}


int main(){
    std::cout << "this is a start" << std::endl;
}
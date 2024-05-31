//
// Created by dstein on 01.09.23.
//
#include <stdexcept>
#include <iostream>
#include "andres/graph/complete-digraph.hxx"


inline void test(const bool& pred) {
    if(!pred) throw std::runtime_error("Test failed.");
}

void testDigraph(size_t numberOfElements = 50){
    typedef andres::graph::CompleteDigraph CDG;

    CDG graph(numberOfElements);

    test(graph.numberOfEdges() == numberOfElements*(numberOfElements - 1));
    test(graph.numberOfVertices() == numberOfElements);

    size_t currentEdge = 0;
    for (size_t i = 0; i < graph.numberOfVertices(); ++i){
        for (size_t j = 0; j < graph.numberOfVertices(); ++j){
            if (i == j) continue;

            std::pair<bool, size_t > foundEdge = graph.findEdge(i, j);

            test(foundEdge.first);

            test(foundEdge.second == currentEdge);

            ++currentEdge;
        }
    }

    size_t i = 0;
    size_t j = 1;
    for (size_t edgeIndex = 0; edgeIndex < graph.numberOfEdges(); ++edgeIndex){

        test(graph.vertexOfEdge(edgeIndex, 0) == i);
        test(graph.vertexOfEdge(edgeIndex, 1) == j);

        if (j == numberOfElements - 1){
            j = 0;
            ++i;
        } else if ((j == i-1) && (i != numberOfElements - 1)){
            j = i + 1;
        } else {
            ++j;
        }
    }

}


int main() {
    testDigraph(5);
    return 0;
}
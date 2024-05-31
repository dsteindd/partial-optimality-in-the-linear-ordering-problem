//
// Created by dstein on 01.09.23.
//
#include <stdexcept>
#include <iostream>
#include "dscc/ground-truth/ground-truth-order-problem.hxx"
#include "andres/ilp/linear-ordering-ilp.hxx"
#include "andres/graph/complete-digraph.hxx"
#include "andres/ilp/gurobi.hxx"


inline void test(const bool &pred) {
    if (!pred) throw std::runtime_error("Test failed.");
}

std::vector<size_t>
getAscendingOrder(size_t numberOfElements) {
    std::vector<size_t> order(numberOfElements);
    for (size_t i = 0; i < numberOfElements; ++i) {
        order[i] = i;
    }
    return order;
}

std::vector<size_t>
getDescendingOrder(size_t numberOfElements) {
    std::vector<size_t> order(numberOfElements);
    for (size_t i = 0; i < numberOfElements; ++i) {
        order[i] = numberOfElements - i - 1;
    }
    return order;
}

std::vector<size_t>
orderToEdgeLabels(
        andres::graph::CompleteDigraph const &digraph,
        std::vector<size_t> const &order
) {
    std::vector<size_t> edgeLabels(order.size() * (order.size() - 1));

    for (size_t i = 0; i < digraph.numberOfVertices(); ++i) {
        for (size_t j = 0; j < digraph.numberOfVertices(); ++j) {
            if (i == j)
                continue;

            size_t edge = digraph.findEdge(i, j).second;

            if (order[i] < order[j]) {
                edgeLabels[edge] = 1;
            } else {
                edgeLabels[edge] = 0;
            }
        }
    }

    return edgeLabels;
}


template<class T>
void testSolver(size_t numberOfElements = 30, double sigma0 = 0.0, double sigma1 = 0.4) {
    typedef GroundTruthOrderProblem<T, std::normal_distribution<T>> GTP;
    typedef andres::graph::CompleteDigraph CompleteDigraph;

    std::vector<size_t> order = getAscendingOrder(numberOfElements);

    // setup
    GTP problem = gaussianGroundTruthProblem<T>(
            0.0,
            sigma0,
            sigma1,
            order
    );

    CompleteDigraph digraph(order.size());

    std::vector<double> edgeCosts(digraph.numberOfEdges());
    std::vector<size_t> edgeLabels(digraph.numberOfEdges(), 0);

    for (size_t i = 0; i < order.size(); ++i) {
        for (size_t j = 0; j < order.size(); ++j) {
            if (i == j) continue;

            size_t edge = digraph.findEdge(i, j).second;
            edgeCosts[edge] = problem.cost(i, j);
        }
    }

    andres::graph::ordering::ilp<andres::ilp::Gurobi>(digraph, edgeCosts, edgeLabels, edgeLabels);

    // test totality constraints on found solution
    for (size_t i = 0; i < digraph.numberOfVertices(); ++i) {
        for (size_t j = 0; j < digraph.numberOfVertices(); ++j) {
            if (i == j)
                continue;

            size_t edge = digraph.findEdge(i, j).second;
            size_t inverseEdge = digraph.findEdge(j, i).second;

            test(edgeLabels[edge] + edgeLabels[inverseEdge] == 1);
        }
    }

    // test transitivity constraints on found solution
    for (size_t edge = 0; edge < digraph.numberOfEdges(); ++edge) {
        size_t v0 = digraph.vertexOfEdge(edge, 0);
        size_t v1 = digraph.vertexOfEdge(edge, 1);

        for (size_t i = 0; i < digraph.numberOfVertices(); ++i) {
            if (i == v0 || i == v1)
                continue;

            size_t edge1 = digraph.findEdge(v0, i).second;
            size_t edge2 = digraph.findEdge(i, v1).second;

            test(edgeLabels[edge1] + edgeLabels[edge2] - edgeLabels[edge] <= 1);
        }
    }

    // test that solution matches ground truth solution
    for (size_t edge = 0; edge < digraph.numberOfEdges(); ++edge){
        size_t v0 = digraph.vertexOfEdge(edge, 0);
        size_t v1 = digraph.vertexOfEdge(edge, 1);

        if (order[v0] < order[v1]){
            test(edgeLabels[edge] == 1);
        } else {
            test(edgeLabels[edge] == 0);
        }
    }
}

int main() {
    testSolver<double>(30, 0.0, 0.0);
    return 0;
}
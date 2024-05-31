#include "andres/graph/digraph.hxx"

#ifndef PARTIAL_OPTIMALITY_LINEAR_ORDERING_RESIDUAL_GRAPH_HXX
#define PARTIAL_OPTIMALITY_LINEAR_ORDERING_RESIDUAL_GRAPH_HXX

typedef andres::graph::Digraph<> Digraph;

template<class R>
inline
Digraph
buildResidualGraph(
        Digraph const &flowGraph,
        std::vector<R> const &capacities,
        std::vector<R> const &flows
) {
    Digraph residualGraph(flowGraph.numberOfVertices());

    for (size_t edgeIndex = 0; edgeIndex < flowGraph.numberOfEdges(); edgeIndex++) {
        size_t u = flowGraph.vertexOfEdge(edgeIndex, 0);
        size_t v = flowGraph.vertexOfEdge(edgeIndex, 1);

        if (capacities[edgeIndex] - flows[edgeIndex] > 0) {
            residualGraph.insertEdge(u, v);
        }

        if (flows[edgeIndex] > 0) {
            residualGraph.insertEdge(v, u);
        }
    }

    return residualGraph;
}


#endif //PARTIAL_OPTIMALITY_LINEAR_ORDERING_RESIDUAL_GRAPH_HXX

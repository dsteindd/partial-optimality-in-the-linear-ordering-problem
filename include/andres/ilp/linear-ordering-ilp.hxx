
#pragma once
#ifndef PARTIAL_OPTIMALITY_LINEAR_ORDERING_LINEAR_ORDERING_ILP_HXX
#define PARTIAL_OPTIMALITY_LINEAR_ORDERING_LINEAR_ORDERING_ILP_HXX

#include <array>
#include <andres/graph/linear-ordering/linear-ordering-problem.hxx>
#include <andres/graph/complete-digraph.hxx>
#include "set"

namespace andres::graph::ordering {

    template<typename ILP, typename ECA, typename ELA>
    inline
    void
    ilp(
            CompleteDigraph const &graph,
            ECA const &edgeCosts,
            ELA const &inputLabels,
            ELA &outputLabels,
            size_t numberOfIterations = std::numeric_limits<size_t>::max()
    ) {
        ILP ilp;
        std::array<size_t , 3> transitivityConstraintVariables{};
        std::array<double, 3> transitivityConstraintsCoefficients{};

        std::array<size_t , 2> totalityConstraintVariables{};
        std::array<double, 2> totalityConstraintCoefficients{};

        auto addTransitivityInequalities = [&]() {
            size_t addedConstraints = 0;

            // add transitivity here
            for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge) {
                if (ilp.label(edge) < .5) {
                    transitivityConstraintVariables[2] = edge;

                    auto const v0 = graph.vertexOfEdge(edge, 0);
                    auto const v1 = graph.vertexOfEdge(edge, 1);

                    for (size_t i = 0; i < graph.numberOfVertices(); ++i) {
                        if (i == v0 || i == v1) {
                            continue;
                        }

                        transitivityConstraintVariables[0] = graph.findEdge(v0, i).second;
                        transitivityConstraintVariables[1] = graph.findEdge(i, v1).second;

                        if (ilp.label(transitivityConstraintVariables[0]) > .5 && ilp.label(transitivityConstraintVariables[1]) > .5) {
                            transitivityConstraintsCoefficients[0] = 1.0;
                            transitivityConstraintsCoefficients[1] = 1.0;
                            transitivityConstraintsCoefficients[2] = -1.0;

                            ilp.addConstraint(transitivityConstraintVariables.begin(), transitivityConstraintVariables.end(), transitivityConstraintsCoefficients.begin(), -1, 1);

                            ++addedConstraints;
                        }
                    }
                }
            }

            return addedConstraints;
        };

        auto addTotalityInequalities = [&]() {
            size_t addedConstraints = 0;

            for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge){
                auto const v0 = graph.vertexOfEdge(edge, 0);
                auto const v1 = graph.vertexOfEdge(edge, 1);

                auto const inverseEdge = graph.findEdge(v1, v0).second;

                if ((ilp.label(edge) < .5) == (ilp.label(inverseEdge) < .5)){
                    totalityConstraintVariables[0] = edge;
                    totalityConstraintVariables[1] = inverseEdge;
                    totalityConstraintCoefficients[0] = 1;
                    totalityConstraintCoefficients[1] = 1;

                    ilp.addConstraint(totalityConstraintVariables.begin(), totalityConstraintVariables.end(), totalityConstraintCoefficients.begin(), 1, 1);

                    ++addedConstraints;
                }
            }

            return addedConstraints;
        };

        ilp.initModel(graph.numberOfEdges(), edgeCosts.data());
        ilp.setStart(inputLabels.begin());

        for (size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i) {
            ilp.optimize();

            if ((addTransitivityInequalities() == 0) && (addTotalityInequalities() == 0)) {
                break;
            }
        }

        for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge) {
            outputLabels[edge] = (ilp.label(edge) < .5) ? 0 : 1;
        }
    }

    template<typename ILP, typename ECA, typename ELA>
    inline
    void
    ilp(
            CompleteDigraph const &graph,
            ECA const &edgeCosts,
            ELA const &inputLabels,
            ELA &outputLabels,
            std::set<std::pair<size_t , size_t >> const & oneEdges,
            size_t numberOfIterations = std::numeric_limits<size_t>::max()
    ) {
        ILP ilp;
//        ilp.setVerbosity(true);
        std::array<size_t , 3> transitivityConstraintVariables{};
        std::array<double, 3> transitivityConstraintsCoefficients{};

        std::array<size_t , 2> totalityConstraintVariables{};
        std::array<double, 2> totalityConstraintCoefficients{};

        std::array<size_t , 1> variableConstraintVariables{};
        std::array<double , 1> variableConstraintCoefficients{};

        auto addTransitivityInequalities = [&]() {
            size_t addedConstraints = 0;

            // add transitivity here
            for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge) {
                if (ilp.label(edge) < .5) {
                    transitivityConstraintVariables[2] = edge;

                    auto const v0 = graph.vertexOfEdge(edge, 0);
                    auto const v1 = graph.vertexOfEdge(edge, 1);

                    for (size_t i = 0; i < graph.numberOfVertices(); ++i) {
                        if (i == v0 || i == v1) {
                            continue;
                        }

                        transitivityConstraintVariables[0] = graph.findEdge(v0, i).second;
                        transitivityConstraintVariables[1] = graph.findEdge(i, v1).second;

                        if (ilp.label(transitivityConstraintVariables[0]) > .5 && ilp.label(transitivityConstraintVariables[1]) > .5) {
                            transitivityConstraintsCoefficients[0] = 1.0;
                            transitivityConstraintsCoefficients[1] = 1.0;
                            transitivityConstraintsCoefficients[2] = -1.0;

                            ilp.addConstraint(
                                    transitivityConstraintVariables.begin(),
                                    transitivityConstraintVariables.end(),
                                    transitivityConstraintsCoefficients.begin(),
                                    -std::numeric_limits<double>::infinity(),
                                    1
                                    );

                            ++addedConstraints;
                        }
                    }
                }
            }

            return addedConstraints;
        };

        auto addTotalityInequalities = [&]() {
            size_t addedConstraints = 0;

            for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge){
                auto const v0 = graph.vertexOfEdge(edge, 0);
                auto const v1 = graph.vertexOfEdge(edge, 1);

                auto const inverseEdge = graph.findEdge(v1, v0).second;

                if ((ilp.label(edge) < .5) == (ilp.label(inverseEdge) < .5)){
                    totalityConstraintVariables[0] = edge;
                    totalityConstraintVariables[1] = inverseEdge;
                    totalityConstraintCoefficients[0] = 1;
                    totalityConstraintCoefficients[1] = 1;

                    ilp.addConstraint(totalityConstraintVariables.begin(), totalityConstraintVariables.end(), totalityConstraintCoefficients.begin(), 1, 1);

                    ++addedConstraints;
                }
            }

            return addedConstraints;
        };

        ilp.initModel(graph.numberOfEdges(), edgeCosts.data());

        for (auto oneEdge : oneEdges){
            size_t v0 = oneEdge.first;
            size_t v1 = oneEdge.second;

            auto edgeFound = graph.findEdge(v0, v1);

            if (!edgeFound.first){
                throw std::runtime_error("One Edge supplied was not part of the graph.");
            }

            size_t edge = edgeFound.second;

            variableConstraintVariables[0] = edge;
            variableConstraintCoefficients[0] = 1.0;

            ilp.addConstraint(
                    variableConstraintVariables.begin(),
                    variableConstraintVariables.end(),
                    variableConstraintCoefficients.begin(),
                    1,
                    1
            );

            auto invertedEdgeFound = graph.findEdge(v1, v0);
            if (!invertedEdgeFound.first){
                throw std::runtime_error("Inverted One Edge supplied was not part of the graph.");
            }

            size_t invertedEdge = invertedEdgeFound.second;

            variableConstraintVariables[0] = invertedEdge;
            variableConstraintCoefficients[0] = 1.0;

            ilp.addConstraint(
                    variableConstraintVariables.begin(),
                    variableConstraintVariables.end(),
                    variableConstraintCoefficients.begin(),
                    0,
                    0
            );

        }



        ilp.setStart(inputLabels.begin());

        for (size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i) {
            ilp.optimize();

            if ((addTransitivityInequalities() == 0) && (addTotalityInequalities() == 0)) {
                break;
            }
        }

        for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge) {
            outputLabels[edge] = (ilp.label(edge) < .5) ? 0 : 1;
        }
    }

}


#endif //PARTIAL_OPTIMALITY_LINEAR_ORDERING_LINEAR_ORDERING_ILP_HXX

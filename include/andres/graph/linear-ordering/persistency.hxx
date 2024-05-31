#ifndef PARTIAL_OPTIMALITY_LINEAR_ORDERING_PERSISTENCY_HXX
#define PARTIAL_OPTIMALITY_LINEAR_ORDERING_PERSISTENCY_HXX

#include "linear-ordering-problem.hxx"
#include "andres/graph/subgraph.hxx"
#include "andres/graph/components.hxx"
#include "andres/graph/graph.hxx"
#include "andres/graph/adjacency.hxx"
#include "andres/graph/max-flow.hxx"
#include "andres/graph/residual-graph.hxx"
#include "andres/graph/components.hxx"
#include <numeric>


template<class R>
class Persistency {
    typedef R Rational;
    typedef LinearOrderingProblem<Rational> ProblemView;

public:
    template<class PROBLEM_VIEW>
    explicit Persistency(PROBLEM_VIEW &problemView);

    void checkTranspositionProposition();

    std::pair<bool, std::vector<size_t >> findHeadingSubProblem();

    std::pair<bool, size_t> findHeadingElement();

    std::pair<bool, size_t> findTrailingElement();

    std::pair<bool, std::pair<size_t, size_t> > findHeadingPair();

    std::pair<bool, std::pair<size_t, size_t> > findTrailingPair();

    void checkGroupingPropositionSinglePair();

    void checkShiftingProposition();

    // todo: Better name!
    void checkABIndependentSubsets();

    void checkElementBetweenCondition();

    // getter
    std::set<std::pair<size_t, size_t>> oneEdges() const;


private:
    ProblemView problemView_;

    std::set<std::pair<size_t, size_t> > oneEdges_;

    Rational positivePart(Rational const &value) const;

    Rational negativePart(Rational const &value) const;

    void transitiveClosureOneEdges(size_t const, size_t const);

    std::pair<bool, std::set<std::pair<size_t, size_t >>> findGroupingPairs(
            const andres::graph::Digraph<> &,
            const std::vector<Rational> &,
            std::set<std::pair<size_t, size_t>> &
    ) const;

    void updateGroupingSinglePairGraph(Digraph &, std::vector<Rational> &) const;

    std::vector<char> computeCutSet(Digraph &, size_t) const;

    Rational largeConstant_;
};

template<class R>
inline
std::vector<char>
Persistency<R>::computeCutSet(
        Digraph &residualGraph,
        size_t source
) const {
    std::vector<char> visited(residualGraph.numberOfVertices());
    std::queue<std::size_t> queue;
    queue.push(source);
    visited[source] = 1;
    while (!queue.empty()) {
        std::size_t w = queue.front();
        queue.pop();
        for (auto it = residualGraph.adjacenciesFromVertexBegin(w);
             it != residualGraph.adjacenciesFromVertexEnd(w); ++it) {
            if (!visited[it->vertex()]) {
                queue.push(it->vertex());
                visited[it->vertex()] = 1;
            }
        }
    }
    return visited;
}


template<class R>
void
Persistency<R>::checkElementBetweenCondition() {


    size_t tripleCount = 0;

    for (size_t a = 0; a < problemView_.numberOfElements(); ++a){
        for (size_t b = 0; b < problemView_.numberOfElements(); ++b){
            for (size_t c = 0; c < problemView_.numberOfElements(); ++c){
                if (a == b || b == c || a == c)
                    continue;

                // check whether b in between a and c

                Rational bacBound = 0;
                Rational acbBound = 0;
                Rational bcaBound = 0;
                Rational cabBound = 0;

                // compute bacBound
                bacBound+= problemView_.deltaCost(a, b);

                for (size_t d = 0; d < problemView_.numberOfElements(); ++d){
                    if (d == a || d == b || d == c)
                        continue;

                    bacBound+= positivePart(problemView_.deltaCost(d, b) + problemView_.deltaCost(a, d));
                }

                // compute acb
                acbBound+= problemView_.deltaCost(b, c);

                for (size_t d = 0; d < problemView_.numberOfElements(); ++d){
                    if (d == a || d == b || d == c)
                        continue;

                    acbBound+= positivePart(problemView_.deltaCost(d, c) + problemView_.deltaCost(b, d));
                }

                // compute bca
                bcaBound+= problemView_.deltaCost(c, b);

                for (size_t d = 0; d < problemView_.numberOfElements(); ++d){
                    if (d == a || d == b || d == c)
                        continue;

                    bcaBound+= positivePart(problemView_.deltaCost(d, b) + problemView_.deltaCost(c, d));
                }


                // compute cab
                cabBound+= problemView_.deltaCost(b, a);

                for (size_t d = 0; d < problemView_.numberOfElements(); ++d){
                    if (d == a || d == b || d == c)
                        continue;

                    cabBound+= positivePart(problemView_.deltaCost(b, d) + problemView_.deltaCost(d, a));
                }

                if ((bacBound <= 0) &&
                        (acbBound <= 0) &&
                        (bcaBound <= 0) &&
                        (cabBound <= 0)){
                    tripleCount += 1;
                }

            }
        }
    }
}


template<class R>
std::pair<bool, std::set<std::pair<size_t, size_t >>>
Persistency<R>::findGroupingPairs(
        const andres::graph::Digraph<> &graph,
        const std::vector<Rational> &edgeWeights,
        std::set<std::pair<size_t, size_t>> &alreadyComputedMaxFlow
) const {
    for (size_t i = 0; i < problemView_.numberOfElements(); ++i) {
        for (size_t j = 0; j < problemView_.numberOfElements(); ++j) {
            if ((i == j) || (alreadyComputedMaxFlow.find(std::make_pair(i, j)) != alreadyComputedMaxFlow.end()))
                continue;

            if ((std::find(oneEdges_.begin(), oneEdges_.end(), std::make_pair(i, j)) != oneEdges_.end())
                || (std::find(oneEdges_.begin(), oneEdges_.end(), std::make_pair(j, i)) != oneEdges_.end())) {
                continue;
            }

            Rational lhs = negativePart(problemView_.deltaCost(i, j));

            // rhs requires calculation of a min-st-cut problem for each (i, j) -> this is necessary since the algorithm needs to restart as soon as a one-edge is found
            andres::graph::MaxFlowPushRelabel<andres::graph::Digraph<>, Rational> maxFlowPushRelabel(
                    graph,
                    edgeWeights.begin(),
                    i, j
            );

            Rational rhs = maxFlowPushRelabel.maxFlow();

            std::vector<Rational> flows(graph.numberOfEdges());
            for (size_t edgeIndex = 0; edgeIndex < graph.numberOfEdges(); ++edgeIndex) {
                flows[edgeIndex] = maxFlowPushRelabel.flow(edgeIndex);
            }

            Digraph residualGraph = buildResidualGraph(graph, edgeWeights, flows);
            std::vector<char> visited = computeCutSet(residualGraph, i);

            Rational rhsFromComponents = 0;

            std::set<std::pair<size_t, size_t>> additionalOneEdges;

            bool isValid = true;

            // test that the constraint for R is fulfilled, otherwise continue with next iteration
            for (auto oneEdge: oneEdges_) {
                size_t a = oneEdge.first;
                size_t b = oneEdge.second;

                if ((visited[b] == 1) && (visited[a] == 0)) {
                    alreadyComputedMaxFlow.insert(std::make_pair(i, j));
                    isValid = false;
                }
            }

            if (!isValid) {
                continue;
            }

            if (visited[j] == 1) {
                throw std::runtime_error("Sink in cut-set.");
            }

            for (size_t from = 0; from < residualGraph.numberOfVertices(); ++from) {
                for (size_t to = 0; to < residualGraph.numberOfVertices(); ++to) {
                    if ((visited[from] == 0) || (visited[to] == 1))
                        continue;

                    auto edgeFound = graph.findEdge(from, to);

                    if (edgeFound.first) {
                        rhsFromComponents += edgeWeights[edgeFound.second];
                    }

                    Rational leftHandSide = negativePart(problemView_.deltaCost(from, to));
                    if ((leftHandSide >= rhs)
                        && (oneEdges_.find({to, from}) == oneEdges_.end())
                        && (oneEdges_.find(std::make_pair(from, to)) == oneEdges_.end())) {
                        additionalOneEdges.insert({from, to});
                    }
                }
            }

            if (rhs != rhsFromComponents) {
                throw std::runtime_error("MaxFlow Value and computed cut-set value do not match.");
            }

            alreadyComputedMaxFlow.insert(std::make_pair(i, j));

            if (!additionalOneEdges.empty()) {
                return std::make_pair(true, additionalOneEdges);
            }
        }
    }
    return std::make_pair(false, std::set<std::pair<size_t, size_t>>());
}

template<class R>
std::set<std::pair<size_t, size_t >>
Persistency<R>::oneEdges() const {
    return oneEdges_;
}

template<class R>
inline
void
Persistency<R>::updateGroupingSinglePairGraph(Digraph &graph, std::vector<Rational> &edgeWeights) const {
    for (auto oneEdge: oneEdges_) {
        size_t index1 = oneEdge.first;
        size_t index2 = oneEdge.second;

        auto foundEdge = graph.findEdge(index2, index1);
        if (!foundEdge.first) {
            size_t edgeIndex = graph.insertEdge(index2, index1);
            assert(edgeIndex == edgeWeights.size());
            edgeWeights.push_back(largeConstant_);
        } else {
            size_t edgeIndex = foundEdge.second;
            edgeWeights[edgeIndex] = largeConstant_;
        }
    }
}

template<class R>
void Persistency<R>::checkABIndependentSubsets() {
    size_t fulfilled = 0;
    size_t const total = problemView_.numberOfElements() * (problemView_.numberOfElements() - 1) / 2;

    for (size_t i = 0; i < problemView_.numberOfElements(); i++) {
        for (size_t j = 0; j < problemView_.numberOfElements(); ++j) {

            bool isValid = true;

            for (size_t k = 0; k < problemView_.numberOfElements(); ++k) {
                if (i == k || j == k) continue;

                if (problemView_.cost(i, k) + problemView_.cost(k, j) > 0 ||
                    problemView_.cost(j, k) + problemView_.cost(k, i) > 0) {
                    isValid = false;
                    break;
                }
            }

            if (isValid) {
                fulfilled++;
            }
        }
    }
}

template<class R>
void Persistency<R>::checkShiftingProposition() {
    size_t fulfilledPairs = 0;
    size_t total = problemView_.numberOfElements() * (problemView_.numberOfElements() - 1);

    for (size_t i = 0; i < problemView_.numberOfElements(); ++i) {
        for (size_t j = 0; j < problemView_.numberOfElements(); ++j) {
            if (i == j) continue;

            Rational lhs = problemView_.deltaCost(i, j);

            Rational rhs = 0;
            bool valid = true;
            for (size_t k = 0; k < problemView_.numberOfElements(); ++k) {
                if (k == i || k == j) continue;

                Rational jkCost = problemView_.deltaCost(j, k);
                if (jkCost > 0) {
                    valid = false;
                    break;
                }

                rhs -= negativePart(jkCost);
            }

            if (valid && (lhs <= rhs)) {
                fulfilledPairs += 1;
            }
        }
    }
}

template<class R>
inline
void Persistency<R>::checkGroupingPropositionSinglePair() {
    andres::graph::Digraph<> graph(problemView_.numberOfElements());
    for (size_t i = 0; i < problemView_.numberOfElements(); ++i) {
        for (size_t j = 0; j < problemView_.numberOfElements(); ++j) {
            if (i == j) continue;

            if ((problemView_.deltaCost(i, j) > 0)
                || (std::find(oneEdges_.begin(), oneEdges_.end(), std::make_pair(j, i)) != oneEdges_.end())) {
                graph.insertEdge(i, j);
            }
        }
    }

    std::vector<Rational> edgeWeights(graph.numberOfEdges());
    for (size_t i = 0; i < problemView_.numberOfElements(); ++i) {
        for (size_t j = 0; j < problemView_.numberOfElements(); ++j) {
            if (i == j) continue;

            auto foundEdge = graph.findEdge(i, j);

            if (!foundEdge.first)
                continue;

            size_t edgeIndex = foundEdge.second;

            edgeWeights[edgeIndex] = positivePart(problemView_.deltaCost(i, j));

            if (std::find(oneEdges_.begin(), oneEdges_.end(), std::make_pair(j, i)) != oneEdges_.end()) {
                edgeWeights[edgeIndex] = largeConstant_;
            }
        }
    }

    bool foundPairs = true;

    std::set<std::pair<size_t, size_t>> alreadyComputedMaxFlows;

    while (foundPairs) {

        auto pairs = findGroupingPairs(graph, edgeWeights, alreadyComputedMaxFlows);

        foundPairs = pairs.first;

        if (foundPairs) {

            auto p = pairs.second;

            for (auto pair: p) {
                size_t i = pair.first;
                size_t j = pair.second;

                if (oneEdges_.find({i, j}) == oneEdges_.end()) {
                    oneEdges_.insert({i, j});
                    transitiveClosureOneEdges(i, j);
                }

                updateGroupingSinglePairGraph(graph, edgeWeights);
            }
        }
    }
}


struct DegreeComparator {

    explicit DegreeComparator(andres::graph::Digraph<> const &digraph) :
            digraph_(digraph) {}

    inline bool operator()(const size_t &v1, const size_t &v2) {
        return digraph_.numberOfEdgesFromVertex(v1) < digraph_.numberOfEdgesFromVertex(v2);
    }

private:
    andres::graph::Digraph<> const &digraph_;
};

template<class R>
std::pair<bool, std::vector<size_t >>
Persistency<R>::findHeadingSubProblem() {
    andres::graph::Digraph<> graph(problemView_.numberOfElements());
    for (size_t i = 0; i < problemView_.numberOfElements(); ++i) {
        for (size_t j = 0; j < problemView_.numberOfElements(); ++j) {
            if (i == j)
                continue;

            if (problemView_.deltaCost(i, j) > 0) {
                graph.insertEdge(i, j);
            }
        }
    }

    std::vector<size_t> vertices(graph.numberOfVertices());
    std::iota(vertices.begin(), vertices.end(), 0);
    std::sort(vertices.begin(), vertices.end(), DegreeComparator(graph));

    for (auto v: vertices) {
        std::vector<size_t> component;
        std::vector<char> visited(graph.numberOfVertices());
        std::queue<std::size_t> queue;
        queue.push(v);
        component.push_back(v);
        visited[v] = 1;
        while (!queue.empty()) {
            std::size_t w = queue.front();
            queue.pop();
            for (auto it = graph.adjacenciesFromVertexBegin(w);
                 it != graph.adjacenciesFromVertexEnd(w); ++it) {
                if (!visited[it->vertex()]) {
                    queue.push(it->vertex());
                    component.push_back(it->vertex());
                    visited[it->vertex()] = 1;
                }
            }
        }
        if (component.size() != graph.numberOfVertices()) {
            return std::make_pair(true, component);
        }
    }

    return std::make_pair(false, std::vector<size_t>());
}


template<class R>
std::pair<bool, size_t>
Persistency<R>::findHeadingElement() {

    if (problemView_.numberOfElements() == 1) {
        return std::make_pair(false, -1);
    }

    for (size_t a = 0; a < problemView_.numberOfElements(); ++a) {
        bool fulfilled = true;

        for (size_t d = 0; d < problemView_.numberOfElements(); ++d) {
            if (a == d)
                continue;

            Rational sum = 0;

            for (size_t dp = 0; dp < problemView_.numberOfElements(); ++dp) {
                if (dp == d)
                    continue;

                sum += problemView_.deltaCost(dp, d);
            }

            if (sum > 0) {
                fulfilled = false;
                break;
            }
        }
    }

    return std::make_pair(false, -1);
}

template<class R>
inline
std::pair<bool, std::pair<size_t, size_t> >
Persistency<R>::findHeadingPair() {
    if (problemView_.numberOfElements() == 2) {
        return std::make_pair(false, std::make_pair(-1, -1));
    }

    for (size_t a = 0; a < problemView_.numberOfElements(); ++a) {
        for (size_t b = 0; b < problemView_.numberOfElements(); ++b) {
            if (a == b)
                continue;

            Rational sum1 = 0;
            Rational sum2 = 0;

            for (size_t d = 0; d < problemView_.numberOfElements(); ++d) {
                if (d == a || d == b)
                    continue;

                Rational innerSum = 0;
                for (size_t dp = 0; dp < problemView_.numberOfElements(); ++dp) {
                    if (dp == a || dp == b || dp == d)
                        continue;

                    innerSum += problemView_.deltaCost(d, dp);
                }

                sum1 += std::max({
                                        problemView_.deltaCost(a, d),
                                        problemView_.deltaCost(a, d) + problemView_.deltaCost(b, d),
                                        innerSum
                                });

                sum2 += std::max({
                                         problemView_.deltaCost(b, d),
                                         problemView_.deltaCost(a, d) + problemView_.deltaCost(b, d),
                                         innerSum
                });
            }

            if ((sum1 <= 0) && (sum2 <= 0)) {
                return std::make_pair(true, std::make_pair(a, b));
            }
        }
    }

    return std::make_pair(false, std::make_pair(-1, -1));
}


template<class R>
inline
std::pair<bool, std::pair<size_t, size_t> >
Persistency<R>::findTrailingPair() {
    if (problemView_.numberOfElements() == 2) {
        return std::make_pair(false, std::make_pair(-1, -1));
    }

    for (size_t a = 0; a < problemView_.numberOfElements(); ++a) {
        for (size_t b = 0; b < problemView_.numberOfElements(); ++b) {
            if (a == b)
                continue;

            Rational sum1 = 0;
            Rational sum2 = 0;

            for (size_t d = 0; d < problemView_.numberOfElements(); ++d) {
                if (d == a || d == b)
                    continue;

                Rational innerSum = 0;
                for (size_t dp = 0; dp < problemView_.numberOfElements(); ++dp) {
                    if (dp == a || dp == b || dp == d)
                        continue;

                    innerSum += problemView_.deltaCost(dp, d);
                }

                sum1 += std::max({
                                         problemView_.deltaCost(d, a),
                                         problemView_.deltaCost(d, a) + problemView_.deltaCost(d, b),
                                         innerSum
                                 });

                sum2 += std::max({
                                         problemView_.deltaCost(d, b),
                                         problemView_.deltaCost(d, b) + problemView_.deltaCost(d,b),
                                         innerSum
                                 });
            }

            if ((sum1 <= 0) && (sum2 <= 0)) {
                return std::make_pair(true, std::make_pair(a, b));
            }
        }
    }

    return std::make_pair(false, std::make_pair(-1, -1));
}

template<class R>
std::pair<bool, size_t>
Persistency<R>::findTrailingElement() {
    if (problemView_.numberOfElements() == 1) {
        return std::make_pair(false, -1);
    }


    for (size_t a = 0; a < problemView_.numberOfElements(); ++a) {
        bool isFulfilled = true;

        for (size_t d = 0; d < problemView_.numberOfElements(); ++d) {
            if (d == a)
                continue;

            Rational sum = 0;
            for (size_t dp = 0; dp < problemView_.numberOfElements(); ++dp) {
                if (dp == d)
                    continue;

                sum += problemView_.deltaCost(d, dp);
            }

            if (sum > 0) {
                isFulfilled = false;
                break;
            }

        }

        if (isFulfilled) {
            return std::make_pair(true, a);
        }
    }

    return std::make_pair(false, -1);


}

template<class R>
template<class PROBLEM_VIEW>
Persistency<R>::Persistency(PROBLEM_VIEW &problemView)
        :
        problemView_(problemView.numberOfElements()),
        oneEdges_() {
    for (size_t i = 0; i < problemView_.numberOfElements(); ++i) {
        for (size_t j = 0; j < problemView_.numberOfElements(); ++j) {
            if (i == j) continue;

            problemView_.cost(i, j) = problemView.cost(i, j);
            largeConstant_ += std::abs(problemView_.cost(i, j)) * 100;
        }
    }


}


template<class R>
void Persistency<R>::checkTranspositionProposition() {

    for (size_t i = 0; i < problemView_.numberOfElements(); ++i) {
        for (size_t j = 0; j < problemView_.numberOfElements(); ++j) {
            if (std::find(oneEdges_.begin(), oneEdges_.end(), std::make_pair(i, j)) != oneEdges_.end()) {
                continue;
            }

            if (i == j) continue;

            Rational lhs = problemView_.deltaCost(i, j);

            Rational rhs = 0;
            for (size_t k = 0; k < problemView_.numberOfElements(); ++k) {
                if (k == i || k == j) continue;

                rhs -= positivePart(problemView_.deltaCost(i, k) + problemView_.deltaCost(k, j));
            }

            if (lhs <= rhs) {
                oneEdges_.insert({i, j});
                transitiveClosureOneEdges(i, j);
            }
        }
    }
}

template<class R>
inline
typename Persistency<R>::Rational
Persistency<R>::positivePart(const Persistency::Rational &value) const {
    if (value >= 0) return value;

    return 0;
}

template<class R>
inline
typename Persistency<R>::Rational
Persistency<R>::negativePart(const Persistency::Rational &value) const {
    if (value <= 0) return -value;

    return 0;
}

template<class R>
inline
void Persistency<R>::transitiveClosureOneEdges(size_t const i, size_t const j) {
    // we assume that oneEdges is already transitively closed

    for (auto oneEdge: oneEdges_) {
        size_t index1 = oneEdge.first;
        size_t index2 = oneEdge.second;

        // (i, j) (index1, index2) => (i, index2)
        // (index1, index2), (i, j) => (index1, j)
        if (j == index1) {
            oneEdges_.insert({i, index2});
        }
        if (i == index2) {
            oneEdges_.insert({index1, j});
        }
    }
}

#endif //PARTIAL_OPTIMALITY_LINEAR_ORDERING_PERSISTENCY_HXX

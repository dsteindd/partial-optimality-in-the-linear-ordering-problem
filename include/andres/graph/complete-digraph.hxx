#ifndef PARTIAL_OPTIMALITY_LINEAR_ORDERING_COMPLETE_DIGRAPH_HXX
#define PARTIAL_OPTIMALITY_LINEAR_ORDERING_COMPLETE_DIGRAPH_HXX

#include "cassert"

namespace andres {
    namespace graph {
        class CompleteDigraph {
        public:
            explicit CompleteDigraph(size_t);

            [[nodiscard]] size_t numberOfVertices() const;

            [[nodiscard]] size_t numberOfEdges() const;

            std::pair<bool, size_t> findEdge(size_t const, size_t const) const;
            size_t vertexOfEdge(size_t, size_t const) const;


        private:
            size_t numberOfElements_;
        };

        inline
        CompleteDigraph::CompleteDigraph(const size_t numberOfElements) :
                numberOfElements_(numberOfElements) {};

        inline
        size_t
        CompleteDigraph::numberOfVertices() const {
            return numberOfElements_;
        }

        inline
        size_t
        CompleteDigraph::numberOfEdges() const {
            return numberOfElements_ * (numberOfElements_ - 1);
        }

        inline
        std::pair<bool, size_t > CompleteDigraph::findEdge(size_t const i, size_t const j) const {
            if (i >= numberOfElements_ || j >= numberOfElements_) return std::make_pair(false, 0);

            size_t edgeIndex = i*(numberOfElements_ - 1) + j - (i < j);

            return std::make_pair(true, edgeIndex);
        }

        inline
        size_t CompleteDigraph::vertexOfEdge(size_t edge, size_t const j) const {
            assert(edge < numberOfEdges());
            assert(j < 2);

            size_t vertex0 = edge / (numberOfElements_ - 1);
            if (j == 0) return vertex0;

            edge = edge - vertex0 * (numberOfElements_ - 1);

            // vertex1 - (vertex0 < vertex1) = edge - vertex
            return edge + (edge >= vertex0);
        }
    }
}


#endif //PARTIAL_OPTIMALITY_LINEAR_ORDERING_COMPLETE_DIGRAPH_HXX

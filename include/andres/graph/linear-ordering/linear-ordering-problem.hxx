#pragma once
#ifndef PARTIAL_OPTIMALITY_LINEAR_ORDERING_LINEAR_ORDERING_PROBLEM_HXX
#define PARTIAL_OPTIMALITY_LINEAR_ORDERING_LINEAR_ORDERING_PROBLEM_HXX

#include <cstdio>
#include <vector>
#include <cassert>
#include <limits>
#include <numeric>
#include <algorithm>

template<class T>
class LinearOrderingProblem {

    typedef T Rational;

public:


    explicit LinearOrderingProblem(const size_t &numberOfElements);

    template<typename PROBLEM_VIEW,
            typename =std::enable_if_t<!std::is_integral_v<PROBLEM_VIEW>>>
    explicit LinearOrderingProblem(PROBLEM_VIEW const &);


    Rational &cost(const size_t, const size_t);
    Rational cost(const size_t, const size_t) const;
    size_t numberOfElements() const;
    Rational deltaCost(const size_t, const size_t) const;
    void merge(size_t, size_t);

private:
    size_t numberOfElements_;
    std::vector<Rational> costs_;

    size_t indexOfPair(const size_t, const size_t) const;
};

template<class T>
template<typename PROBLEM_VIEW, typename>
inline
LinearOrderingProblem<T>::LinearOrderingProblem(const PROBLEM_VIEW & problemView):
        LinearOrderingProblem(problemView.numberOfElements())
{
    for (size_t i = 0; i < problemView.numberOfElements(); ++i){
        for (size_t j = 0; j < problemView.numberOfElements(); ++j){
            if (i == j)
                continue;

            costs_[indexOfPair(i, j)] = problemView.cost(i, j);
        }
    }
}

template<class T>
typename LinearOrderingProblem<T>::Rational LinearOrderingProblem<T>::deltaCost(const size_t i, const size_t j) const {
    return cost(i, j) - cost(j, i);
}

template<class T>
inline
LinearOrderingProblem<T>::LinearOrderingProblem(size_t const & numberOfElements):
        numberOfElements_(numberOfElements),
        costs_(numberOfElements_ * (numberOfElements_ - 1)) {}

template<class T>
inline
typename LinearOrderingProblem<T>::Rational &
LinearOrderingProblem<T>::cost(const size_t i, const size_t j) {
    return costs_[indexOfPair(i, j)];
}

template<class T>
typename LinearOrderingProblem<T>::Rational
LinearOrderingProblem<T>::cost(const size_t i, const size_t j) const {
    return costs_[indexOfPair(i, j)];
}

template<class T>
inline
size_t LinearOrderingProblem<T>::indexOfPair(const size_t i, const size_t j) const {
    assert(i != j);

    if (i < j){
        return j*(j-1) + i;
    } else {
        return i*i + j;
    }
}

template<class T>
inline
size_t LinearOrderingProblem<T>::numberOfElements() const {
    return numberOfElements_;
}

template<class T>
inline
void
LinearOrderingProblem<T>::merge(size_t i, size_t j) {
    if (i > j){
        std::swap(i, j);
    }

    size_t lastElement = numberOfElements_ - 1;

    if (j != lastElement){
        // swap j with last element
        for (size_t k = 0; k < numberOfElements_; ++k){
            if (k == j || k == lastElement)
                continue;
            std::swap(costs_[indexOfPair(k, j)], costs_[indexOfPair(k, lastElement)]);
            std::swap(costs_[indexOfPair(j, k)], costs_[indexOfPair(lastElement, k)]);
        }
    }

    // regard i as the merge element
    // cost[i, k] = cost[j, k] + cost[i, k] = cost[lastElement, k] + cost[i, k]
    for (size_t k = 0; k < numberOfElements_; ++k){
        if (k == i || k == lastElement || k == j)
            continue;
        costs_[indexOfPair(i, k)] = costs_[indexOfPair(lastElement, k)] + costs_[indexOfPair(i, k)];
        costs_[indexOfPair(k, i)] = costs_[indexOfPair(k, lastElement)] + costs_[indexOfPair(k, i)];

    }

    numberOfElements_--;
    costs_.resize(numberOfElements_*(numberOfElements_ - 1));

}

template<class T>
inline
LinearOrderingProblem<T>
project(LinearOrderingProblem<T> const & problem, std::vector<size_t> const & indices){
    LinearOrderingProblem<T> projectedProblem(indices.size());

    for (size_t i = 0; i < indices.size(); ++i){
        for (size_t j = 0; j < indices.size(); ++j){
            if (i == j)
                continue;
            projectedProblem.cost(i, j) = problem.cost(indices[i], indices[j]);
        }
    }

    return projectedProblem;
}


template<class T>
inline
std::tuple<LinearOrderingProblem<T>, size_t, std::vector<size_t >>
merge(LinearOrderingProblem<T> const & problem, std::vector<size_t> const & indices) {
    LinearOrderingProblem<T> mergedProblem(problem);

    // maps from mergedProblem label to former label
    std::vector<size_t > labels(mergedProblem.numberOfElements());
    std::iota(labels.begin(), labels.end(), 0);

    size_t n1 = indices[0];

    for (size_t i = 1; i < indices.size(); ++i) {
        // in merged problem
        n1 = std::find(labels.begin(), labels.end(), n1) - labels.begin();
        size_t n2 = std::find(labels.begin(), labels.end(), indices[i]) - labels.begin();
        if (n2 < n1){
            std::swap(n2, n1);
        }
        labels[n2] = mergedProblem.numberOfElements() - 1;

        mergedProblem.merge(n1, n2);
    }

    labels.resize(mergedProblem.numberOfElements());

    return std::make_tuple(mergedProblem, n1, labels);
}



#endif //PARTIAL_OPTIMALITY_LINEAR_ORDERING_LINEAR_ORDERING_PROBLEM_HXX

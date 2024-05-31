//
// Created by dstein on 01.09.23.
//
#include <stdexcept>
#include <numeric>
#include <bits/stdc++.h>
#include "andres/graph/linear-ordering/linear-ordering-problem.hxx"

inline void test(const bool &pred) {
    if (!pred) throw std::runtime_error("Test failed.");
}

template<class T>
void testCostSetting(size_t numberOfElements = 4) {
    typedef LinearOrderingProblem<T> LOP;

    // setup
    LOP problem(numberOfElements);

    for (size_t i = 0; i < numberOfElements; ++i) {
        for (size_t j = 0; j < i; ++j) {
            problem.cost(i, j) = static_cast<T>(i * i + j * j * j);
        }

        for (size_t j = i + 1; j < numberOfElements; ++j) {
            problem.cost(i, j) = static_cast<T>(i * i * i + j * j);
        }
    }

    // assert
    for (size_t i = 0; i < numberOfElements; ++i) {
        for (size_t j = 0; j < i; ++j) {
            test(problem.cost(i, j) == static_cast<T>(i * i + j * j * j));
        }

        for (size_t j = i + 1; j < numberOfElements; ++j) {
            test(problem.cost(i, j) == static_cast<T>(i * i * i + j * j));
        }
    }
}

template<class T>
void testProjection(size_t numberOfElements = 5) {
    LinearOrderingProblem<T> problem(numberOfElements);

    for (size_t i = 0; i < problem.numberOfElements(); ++i) {
        for (size_t j = 0; j < problem.numberOfElements(); ++j) {
            if (i == j)
                continue;
            problem.cost(i, j) = i + j * j;
        }
    }

    {
        std::vector<size_t> indices(numberOfElements / 2, 0);
        std::iota(indices.begin(), indices.end(), 0);


        LinearOrderingProblem<T> projectedProblem = project(problem, indices);

        for (size_t i = 0; i < projectedProblem.numberOfElements(); ++i) {
            for (size_t j = 0; j < projectedProblem.numberOfElements(); ++j) {
                if (i == j)
                    continue;

                test(projectedProblem.cost(i, j) == problem.cost(indices[i], indices[j]));
            }
        }
    }

    {
        std::vector<size_t> indices(numberOfElements / 2);
        std::iota(indices.begin(), indices.end(), numberOfElements / 2);

        LinearOrderingProblem<T> projectedProblem = project(problem, indices);

        for (size_t i = 0; i < projectedProblem.numberOfElements(); ++i) {
            for (size_t j = 0; j < projectedProblem.numberOfElements(); ++j) {
                if (i == j)
                    continue;

                test(projectedProblem.cost(i, j) == problem.cost(indices[i], indices[j]));
            }
        }
    }
}

template<class T>
void testMerging(size_t numberOfElements = 5) {
    LinearOrderingProblem<T> problem(numberOfElements);

    for (size_t i = 0; i < problem.numberOfElements(); ++i) {
        for (size_t j = 0; j < problem.numberOfElements(); ++j) {
            if (i == j)
                continue;
            problem.cost(i, j) = i + j * j;
        }
    }

    {
        std::vector<size_t> indices(numberOfElements / 2, 0);
        std::iota(indices.begin(), indices.end(), 0);

        auto mergeResult = merge(problem, indices);

        LinearOrderingProblem<T> mergedProblem = std::get<0>(mergeResult);
        size_t mergedElement = std::get<1>(mergeResult);
        std::vector<size_t> labels = std::get<2>(mergeResult);

        // test that the labels indeed do not contain merged indices anymore
        for (size_t index: indices) {
            if (index == labels[mergedElement])
                continue;
            test(std::find(labels.begin(), labels.end(), index) == labels.end());
        }

        test(mergedProblem.numberOfElements() == numberOfElements - indices.size() + 1);

        for (size_t i = 0; i < mergedProblem.numberOfElements(); ++i) {
            for (size_t j = 0; j < mergedProblem.numberOfElements(); ++j) {
                if (i == j)
                    continue;

                if (i != mergedElement && j != mergedElement) {
                    test(problem.cost(labels[i], labels[j]) = mergedProblem.cost(i, j));
                } else if (i == mergedElement) {
                    // mergeCost[i, j] = sum_indices problem.cost[labels[k], labels[j]]
                    T cost = 0;
                    for (auto k: indices) {
                        if (k == j)
                            continue;
                        cost += problem.cost(k, labels[j]);
                    }
                    test(cost == mergedProblem.cost(i, j));

                } else if (j == mergedElement) {
                    // cost[labels[i], j] = sum indices cost[labels[i], j]
                    T cost = 0;
                    for (auto k: indices) {
                        if (k == i)
                            continue;
                        cost += problem.cost(labels[i], k);
                    }
                    test(cost == mergedProblem.cost(i, j));
                }

            }
        }
    }

    {
        std::vector<size_t> indices(numberOfElements / 3, 0);
        std::iota(indices.begin(), indices.end(), 2);

        auto mergeResult = merge(problem, indices);

        LinearOrderingProblem<T> mergedProblem = std::get<0>(mergeResult);
        size_t mergedElement = std::get<1>(mergeResult);
        std::vector<size_t> labels = std::get<2>(mergeResult);

        // test that the labels indeed do not contain merged indices anymore
        for (size_t index: indices) {
            if (index == labels[mergedElement])
                continue;
            test(std::find(labels.begin(), labels.end(), index) == labels.end());
        }

        test(mergedProblem.numberOfElements() == numberOfElements - indices.size() + 1);

        for (size_t i = 0; i < mergedProblem.numberOfElements(); ++i) {
            for (size_t j = 0; j < mergedProblem.numberOfElements(); ++j) {
                if (i == j)
                    continue;

                if (i != mergedElement && j != mergedElement) {
                    test(problem.cost(labels[i], labels[j]) = mergedProblem.cost(i, j));
                } else if (i == mergedElement) {
                    // mergeCost[i, j] = sum_indices problem.cost[labels[k], labels[j]]
                    T cost = 0;
                    for (auto k: indices) {
                        if (k == j)
                            continue;
                        cost += problem.cost(k, labels[j]);
                    }
                    test(cost == mergedProblem.cost(i, j));
                } else if (j == mergedElement) {
                    // cost[labels[i], j] = sum indices cost[labels[i], j]
                    T cost = 0;
                    for (auto k: indices) {
                        if (k == i)
                            continue;
                        cost += problem.cost(labels[i], k);
                    }
                    test(cost == mergedProblem.cost(i, j));
                }

            }
        }
    }
}

template<class T>
inline
void
testMergingThreeElements(){


    {
        LinearOrderingProblem<T> problem(3);

        problem.cost(0, 1) = 1;
        problem.cost(0, 2) = 2;
        problem.cost(1, 2) = 3;
        problem.cost(1, 0) = 4;
        problem.cost(2, 0) = 5;
        problem.cost(2, 1) = 6;

        problem.merge(1, 2);
        test(problem.numberOfElements() == 2);
        test(problem.cost(0, 1) == 1 + 2);
        test(problem.cost(1, 0) == 4 + 5);
    }
    {
        LinearOrderingProblem<T> problem(3);

        problem.cost(0, 1) = 1;
        problem.cost(0, 2) = 2;
        problem.cost(1, 2) = 3;
        problem.cost(1, 0) = 4;
        problem.cost(2, 0) = 5;
        problem.cost(2, 1) = 6;

        problem.merge(0, 2);
        test(problem.numberOfElements() == 2);
        test(problem.cost(0, 1) == 1 + 6);
        test(problem.cost(1, 0) == 4 + 3);
    }
    {
        LinearOrderingProblem<T> problem(3);

        problem.cost(0, 1) = 1;
        problem.cost(0, 2) = 2;
        problem.cost(1, 2) = 3;
        problem.cost(1, 0) = 4;
        problem.cost(2, 0) = 5;
        problem.cost(2, 1) = 6;

        problem.merge(1, 2);
        test(problem.numberOfElements() == 2);
        test(problem.cost(0, 1) == 1 + 2);
        test(problem.cost(1, 0) == 4 + 5);
    }
}

int main() {
    testCostSetting<double>(100);
    testCostSetting<long>(100);
    testProjection<double>(100);
    testProjection<long>(100);
    testMerging<double>(100);
    testMerging<long>(100);
    testMergingThreeElements<double>();
    testMergingThreeElements<long>();
    return 0;
}
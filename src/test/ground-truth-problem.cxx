//
// Created by dstein on 01.09.23.
//
#include <stdexcept>
#include <iostream>
#include "dscc/ground-truth/ground-truth-order-problem.hxx"

inline void test(const bool& pred) {
    if(!pred) throw std::runtime_error("Test failed.");
}

std::vector<size_t >
getAscendingOrder(size_t  numberOfElements){
    std::vector<size_t > order(numberOfElements);
    for (size_t i = 0; i < numberOfElements; ++i){
        order[i] = i;
    }
    return order;
}

std::vector<size_t >
getDescendingOrder(size_t  numberOfElements){
    std::vector<size_t > order(numberOfElements);
    for (size_t i = 0; i < numberOfElements; ++i){
        order[i] = numberOfElements - i - 1;
    }
    return order;
}




template<class T>
void testCostSetting(size_t numberOfElements = 50){
    typedef GroundTruthOrderProblem<T, std::normal_distribution<T>> GTP;

    {
        std::vector<size_t > order = getAscendingOrder(numberOfElements);

        // setup
        GTP problem = gaussianGroundTruthProblem<T>(
                0.0,
                0.0,
                0.0,
                order
        );

        for (size_t i = 0; i < problem.numberOfElements(); ++i) {
            for (size_t j = 0; j < problem.numberOfElements(); ++j){
                if (i == j)
                    continue;
                if (order[i] < order[j])
                    test(problem.cost(i, j) == -1);
                else{
                    test(problem.cost(i, j) == 1);
                }
            }
        }
    }

    {
        std::vector<size_t > order = getDescendingOrder(numberOfElements);

        // setup
        GTP problem = gaussianGroundTruthProblem<T>(
                0.0,
                0.0,
                0.0,
                order
        );

        for (size_t i = 0; i < problem.numberOfElements(); ++i) {
            for (size_t j = 0; j < problem.numberOfElements(); ++j){
                if (i == j) continue;
                if (order[i] < order[j])
                    test(problem.cost(i, j) == -1);
                else
                    test(problem.cost(i, j) == 1);
            }
        }
    }

    {
        std::vector<size_t > order = getAscendingOrder(numberOfElements);

        // setup
        GTP problem = gaussianGroundTruthProblem<T>(
                0.0,
                0.0,
                0.0,
                order
        );

        for (size_t i = 0; i < problem.numberOfElements(); ++i) {
            for (size_t j = 0; j < problem.numberOfElements(); ++j){
                if (i == j) continue;
                if (order[i] < order[j])
                    test(problem.cost(i, j) == -1);
                else
                    test(problem.cost(i, j) == 1);
            }
        }
    }

    {
        std::vector<size_t > order = getAscendingOrder(numberOfElements);

        // setup
        GTP problem = gaussianGroundTruthProblem<T>(
                1.0,
                0.0,
                0.0,
                order
        );

        for (size_t i = 0; i < problem.numberOfElements(); ++i) {
            for (size_t j = 0; j < problem.numberOfElements(); ++j){
                if (i == j) continue;
                test(problem.cost(i, j) == 0);
            }
        }
    }

    {
        std::vector<size_t > order = getDescendingOrder(numberOfElements);

        // setup
        GTP problem = gaussianGroundTruthProblem<T>(
                1.0,
                0.0,
                0.0,
                order
        );

        for (size_t i = 0; i < problem.numberOfElements(); ++i) {
            for (size_t j = 0; j < problem.numberOfElements(); ++j){
                if (i == j) continue;
                test(problem.cost(i, j) == 0);
            }
        }
    }


    {
        std::vector<size_t > order = getAscendingOrder(numberOfElements);

        // setup
        GTP problem = gaussianGroundTruthProblem<T>(
                0.0,
                0.0001,
                0.3,
                order
        );

        for (size_t i = 0; i < problem.numberOfElements(); ++i) {
            for (size_t j = 0; j < problem.numberOfElements(); ++j){
                if (i == j) continue;

                if (order[i] < order[j])
                    test(problem.cost(i, j) < 0);
                else
                    test(problem.cost(i, j) > 0);

            }
        }
    }
}

bool isClose(double x, double y, double tol = 1e-6){
    return std::abs(x - y) < tol;
}

void
testCostScaling(size_t numberOfElements, double alpha = 0.1, double sigma = 0.2, double factor = 5){
    GroundTruthOrderProblem<double, std::normal_distribution<double>> problem1 = gaussianGroundTruthProblem<double>(numberOfElements, alpha, sigma, sigma);

    GroundTruthOrderProblem<double, std::normal_distribution<double>> problem2 = gaussianGroundTruthProblem<double>(numberOfElements, alpha, factor*sigma, factor*sigma);

    for (size_t i = 0; i < numberOfElements; ++i) {
        for (size_t j = 0; j < numberOfElements; ++j){
            if (i == j)
                continue;

            if (i < j){
                double mu = -1 + alpha;
                test(isClose(factor*(problem1.cost(i, j) - mu), problem2.cost(i, j) - mu));
            }
            else {
                double mu = 1 - alpha;
                test(isClose(factor*(problem1.cost(i, j) - mu), problem2.cost(i, j) - mu));
            }
        }
    }
}


int main() {
    testCostSetting<double>();
    testCostScaling(100, 0, 0.2, 5);
    testCostScaling(100, 0.5, 1.0, 5);
    return 0;
}
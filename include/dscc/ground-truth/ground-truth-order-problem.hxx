//
// Created by dstein on 08.09.23.
//

#ifndef PARTIAL_OPTIMALITY_LINEAR_ORDERING_GROUND_TRUTH_ORDER_PROBLEM_HXX
#define PARTIAL_OPTIMALITY_LINEAR_ORDERING_GROUND_TRUTH_ORDER_PROBLEM_HXX

#include "random"
#include "cassert"

template<class R, class D = std::normal_distribution<R>>
class GroundTruthOrderProblem {
    typedef R Rational;
    typedef D Distribution;

public:
    GroundTruthOrderProblem(Distribution &, Distribution &, std::vector<size_t> const &, size_t = 42);

    size_t numberOfElements() const;

    Rational &cost(size_t, size_t);

    Rational cost(size_t, size_t) const;

private:
    size_t numberOfElements_;
    std::default_random_engine engine_;
    Distribution &intraDistribution_;
    Distribution &interDistribution_;
    std::vector<Rational> costs_;

    size_t indexOfPair(size_t, size_t) const;
};


template<class R, class D>
inline
GroundTruthOrderProblem<R, D>::GroundTruthOrderProblem(
        Distribution &interDistribution,
        Distribution &intraDistribution,
        const std::vector<size_t> &order,
        size_t const seed
):
        numberOfElements_(order.size()),
        interDistribution_(interDistribution),
        intraDistribution_(intraDistribution),
        costs_(order.size() * (order.size() - 1)),
        engine_(seed) {

    for (size_t i = 0; i < numberOfElements_; ++i) {
        for (size_t j = 0; j < numberOfElements_; ++j) {
            if (i == j)
                continue;

            double cost;
            if (order[i] < order[j]) {
                cost = intraDistribution_(engine_);
            } else {
                cost = interDistribution_(engine_);
            }

            if (std::is_integral<R>()) {
                costs_[indexOfPair(i, j)] = std::round(cost);
            } else {
                costs_[indexOfPair(i, j)] = cost;
            }
        }
    }

}

template<class R, class D>
inline
typename GroundTruthOrderProblem<R, D>::Rational &
GroundTruthOrderProblem<R, D>::cost(const size_t i, const size_t j) {
    return costs_[indexOfPair(i, j)];
}

template<class R, class D>
inline
typename GroundTruthOrderProblem<R, D>::Rational
GroundTruthOrderProblem<R, D>::cost(const size_t i, const size_t j) const {
    return costs_[indexOfPair(i, j)];
}


template<class R, class D>
size_t GroundTruthOrderProblem<R, D>::indexOfPair(const size_t i, const size_t j) const {
    assert(i != j);

    if (i < j) {
        return j * (j - 1) + i;
    } else {
        return i * i + j;
    }
}

template<class R, class D>
inline
size_t GroundTruthOrderProblem<R, D>::numberOfElements() const {
    return numberOfElements_;
}


template<class R>
inline
GroundTruthOrderProblem<R, std::normal_distribution<R>> gaussianGroundTruthProblem(
        double const &alpha,
        double const &sigma0,
        double const &sigma1,
        const std::vector<size_t> &order,
        size_t seed = 42
) {
    typedef GroundTruthOrderProblem<R, std::normal_distribution<double>> GTP;

    assert(alpha >= 0);
    assert(alpha <= 1);
    assert(sigma0 >= 0);
    assert(sigma1 >= 0);

    R sigma = sigma0 + alpha * (sigma1 - sigma0);
    R muInter = 1 - alpha;
    R muIntra = -1 + alpha;

    std::normal_distribution<R> interDistribution(muInter, sigma);
    std::normal_distribution<R> intraDistribution(muIntra, sigma);

    GTP problem = GTP(interDistribution, intraDistribution, order, seed);

    return problem;
}


template<class R>
inline
GroundTruthOrderProblem<R, std::normal_distribution<double>> gaussianGroundTruthProblem(
        size_t const numberOfElements,
        double const &alpha,
        double const &sigma0,
        double const &sigma1,
        double const &multiplier = 1,
        size_t seed = 42
) {
    typedef GroundTruthOrderProblem<R, std::normal_distribution<double>> GTP;

    assert(alpha >= 0);
    assert(alpha <= 1);
    assert(sigma0 >= 0);
    assert(sigma1 >= 0);

    double sigma = sigma0 + alpha * (sigma1 - sigma0);
    double muInter = 1 - alpha;
    double muIntra = -1 + alpha;

    sigma *= multiplier;
    muInter *= multiplier;
    muIntra *= multiplier;

    std::normal_distribution<double> interDistribution(muInter, sigma);
    std::normal_distribution<double> intraDistribution(muIntra, sigma);

    // standard order
    std::vector<size_t> order(numberOfElements);
    std::iota(order.begin(), order.end(), 0);

    GTP problem = GTP(interDistribution, intraDistribution, order, seed);

    return problem;
}


#endif //PARTIAL_OPTIMALITY_LINEAR_ORDERING_GROUND_TRUTH_ORDER_PROBLEM_HXX

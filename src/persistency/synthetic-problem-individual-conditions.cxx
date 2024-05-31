#include "dscc/ground-truth/ground-truth-order-problem.hxx"
#include "andres/graph/linear-ordering/persistency.hxx"
#include "iostream"
#include "fstream"
#include "sstream"
#include "chrono"
#include "filesystem"
#include "argparse.hxx"

struct PersistencyPerformance {
public:
    PersistencyPerformance() = default;

    size_t initialVars;
    size_t remainingVars;
    double eliminatedVarsFraction;
};

typedef GroundTruthOrderProblem<long, std::normal_distribution<double>> GTP;

template<class R>
PersistencyPerformance testHeadingSubproblemCondition(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<R> persistency(currentProblem);

        auto foundHeadingSubProblem = persistency.findHeadingSubProblem();
        if (foundHeadingSubProblem.first) {
            // find inverse set of variables
            auto indices = foundHeadingSubProblem.second;
            std::vector<size_t> indicesBack;
            for (size_t i = 0; i < currentProblem.numberOfElements(); ++i) {
                if (std::find(indices.begin(), indices.end(), i) == indices.end()) {
                    indicesBack.push_back(i);
                }
            }

            problemQueue.push(project(currentProblem, indices));
            problemQueue.push(project(currentProblem, indicesBack));
        } else {
            irreducibleProblems.push_back(currentProblem);
        }
    }

    size_t remainingVariables = 0;
    for (auto &irreducibleProblem: irreducibleProblems) {
        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1);
    }

    double eliminatedVariableFraction =
            static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}


template<class R>
PersistencyPerformance testHeadingElementCondition(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<R> persistency(currentProblem);

        auto foundHeadingSubElement = persistency.findHeadingElement();
        if (foundHeadingSubElement.first) {
            auto index = foundHeadingSubElement.second;
            std::vector<size_t> otherIndices;
            for (size_t i = 0; i < currentProblem.numberOfElements(); ++i) {
                if (i != index) {
                    otherIndices.push_back(i);
                }
            }
            problemQueue.push(project(currentProblem, {index}));
            problemQueue.push(project(currentProblem, otherIndices));
        } else {

            irreducibleProblems.push_back(currentProblem);
        }
    }


    size_t remainingVariables = 0;
    for (auto &irreducibleProblem: irreducibleProblems) {
        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1);
    }

    double eliminatedVariableFraction =
            static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}

template<class R>
PersistencyPerformance testTrailingElementCondition(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<R> persistency(currentProblem);

        auto foundHeadingSubElement = persistency.findTrailingElement();
        if (foundHeadingSubElement.first) {
            auto index = foundHeadingSubElement.second;
            std::vector<size_t> otherIndices;
            for (size_t i = 0; i < currentProblem.numberOfElements(); ++i) {
                if (i != index) {
                    otherIndices.push_back(i);
                }
            }
            problemQueue.push(project(currentProblem, {index}));
            problemQueue.push(project(currentProblem, otherIndices));
        } else {

            irreducibleProblems.push_back(currentProblem);
        }
    }

    size_t remainingVariables = 0;
    for (auto &irreducibleProblem: irreducibleProblems) {
        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1);
    }

    double eliminatedVariableFraction =
            static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}

template<class R>
PersistencyPerformance testHeadingPairCondition(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<R> persistency(currentProblem);

        auto foundHeadingPair = persistency.findHeadingPair();
        if (foundHeadingPair.first) {
            size_t a = foundHeadingPair.second.first;
            size_t b = foundHeadingPair.second.second;
            std::vector<size_t> otherIndices;
            for (size_t i = 0; i < currentProblem.numberOfElements(); ++i) {
                if (i != a && i != b) {
                    otherIndices.push_back(i);
                }
            }
            problemQueue.push(project(currentProblem, {a, b}));
            problemQueue.push(project(currentProblem, otherIndices));
        } else {
            irreducibleProblems.push_back(currentProblem);
        }
    }

    size_t remainingVariables = 0;
    for (auto &irreducibleProblem: irreducibleProblems) {
        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1);
    }

    double eliminatedVariableFraction =
            static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}


template<class R>
PersistencyPerformance testTrailingPairCondition(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<R> persistency(currentProblem);

        auto foundHeadingPair = persistency.findTrailingPair();
        if (foundHeadingPair.first) {
            size_t a = foundHeadingPair.second.first;
            size_t b = foundHeadingPair.second.second;
            std::vector<size_t> otherIndices;
            for (size_t i = 0; i < currentProblem.numberOfElements(); ++i) {
                if (i != a && i != b) {
                    otherIndices.push_back(i);
                }
            }
            problemQueue.push(project(currentProblem, {a, b}));
            problemQueue.push(project(currentProblem, otherIndices));
        } else {
            irreducibleProblems.push_back(currentProblem);
        }
    }

    size_t remainingVariables = 0;
    for (auto &irreducibleProblem: irreducibleProblems) {
        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1);
    }

    double eliminatedVariableFraction =
            static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}

template<class R>
PersistencyPerformance testTranspositionCondition(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    Persistency<R> persistency(problem);
    persistency.checkTranspositionProposition();

    size_t eliminatedVars = 2 * persistency.oneEdges().size();

    size_t remainingVariables = problem.numberOfElements() * (problem.numberOfElements() - 1) - eliminatedVars;

    double eliminatedVariableFraction =
            static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}

template<class R>
PersistencyPerformance testGroupingSinglePairCondition(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    Persistency<R> persistency(problem);
    persistency.checkGroupingPropositionSinglePair();

    size_t eliminatedVars = 2 * persistency.oneEdges().size();

    size_t remainingVariables = problem.numberOfElements() * (problem.numberOfElements() - 1) - eliminatedVars;

    double eliminatedVariableFraction =
            static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}


void testHeadingSubproblemCondition(
        std::ofstream &out,
        size_t numberOfElements,
        double alpha,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
) {
    std::vector<size_t> seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);

    std::cout << "Collecting statistics for numberOfElements=" << numberOfElements << " and " << " alpha=" << alpha
              << std::endl;

    for (auto seed: seeds) {

        GroundTruthOrderProblem<long, std::normal_distribution<double>> groundTruthProblem = gaussianGroundTruthProblem<long>(
                numberOfElements,
                alpha,
                sigma0,
                sigma1,
                multiplier,
                seed);

        LinearOrderingProblem<double> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto performance = testHeadingSubproblemCondition(problem);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds >(end - start).count();

        out << seed << ","
            << numberOfElements << ","
            << alpha << ","
            << sigma0 << ","
            << sigma1 << ","
            << performance.initialVars << ","
            << performance.remainingVars << ","
            << performance.eliminatedVarsFraction << ","
            << duration
            << std::endl;
    }
}

void testHeadingElementCondition(
        std::ofstream &out,
        size_t numberOfElements,
        double alpha,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
) {
    std::vector<size_t> seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);

    std::cout << "Collecting statistics for numberOfElements=" << numberOfElements << " and " << " alpha=" << alpha
              << std::endl;

    for (auto seed: seeds) {

        GTP groundTruthProblem = gaussianGroundTruthProblem<long>(numberOfElements, alpha, sigma0, sigma1,
                                                                  multiplier, seed);

        LinearOrderingProblem<long> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto performance = testHeadingElementCondition(problem);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        out << seed << ","
            << numberOfElements << ","
            << alpha << ","
            << sigma0 << ","
            << sigma1 << ","
            << performance.initialVars << ","
            << performance.remainingVars << ","
            << performance.eliminatedVarsFraction << ","
            << duration
            << std::endl;
    }
}

void testTrailingElementCondition(
        std::ofstream &out,
        size_t numberOfElements,
        double alpha,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
) {
    std::vector<size_t> seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);

    std::cout << "Collecting statistics for numberOfElements=" << numberOfElements << " and " << " alpha=" << alpha
              << std::endl;

    for (auto seed: seeds) {

        GTP groundTruthProblem = gaussianGroundTruthProblem<long>(numberOfElements, alpha, sigma0, sigma1,
                                                                  multiplier, seed);

        LinearOrderingProblem<long> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto performance = testTrailingElementCondition(problem);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        out << seed << ","
            << numberOfElements << ","
            << alpha << ","
            << sigma0 << ","
            << sigma1 << ","
            << performance.initialVars << ","
            << performance.remainingVars << ","
            << performance.eliminatedVarsFraction << ","
            << duration
            << std::endl;
    }
}

void testHeadingPairCondition(
        std::ofstream &out,
        size_t numberOfElements,
        double alpha,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
) {
    std::vector<size_t> seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);

    std::cout << "Collecting statistics for numberOfElements=" << numberOfElements << " and " << " alpha=" << alpha
              << std::endl;

    for (auto seed: seeds) {

        GTP groundTruthProblem = gaussianGroundTruthProblem<long>(numberOfElements, alpha, sigma0, sigma1,
                                                                  multiplier, seed);

        LinearOrderingProblem<long> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto performance = testHeadingPairCondition(problem);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        out << seed << ","
            << numberOfElements << ","
            << alpha << ","
            << sigma0 << ","
            << sigma1 << ","
            << performance.initialVars << ","
            << performance.remainingVars << ","
            << performance.eliminatedVarsFraction << ","
            << duration
            << std::endl;
    }
}


void testTrailingPairCondition(
        std::ofstream &out,
        size_t numberOfElements,
        double alpha,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
) {
    std::vector<size_t> seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);

    std::cout << "Collecting statistics for numberOfElements=" << numberOfElements << " and " << " alpha=" << alpha
              << std::endl;

    for (auto seed: seeds) {

        GTP groundTruthProblem = gaussianGroundTruthProblem<long>(numberOfElements, alpha, sigma0, sigma1,
                                                                  multiplier, seed);

        LinearOrderingProblem<long> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto performance = testTrailingPairCondition(problem);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        out << seed << ","
            << numberOfElements << ","
            << alpha << ","
            << sigma0 << ","
            << sigma1 << ","
            << performance.initialVars << ","
            << performance.remainingVars << ","
            << performance.eliminatedVarsFraction << ","
            << duration
            << std::endl;
    }
}

void testTranspositionCondition(
        std::ofstream &out,
        size_t numberOfElements,
        double alpha,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
) {
    std::vector<size_t> seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);

    std::cout << "Collecting statistics for numberOfElements=" << numberOfElements << " and " << " alpha=" << alpha
              << std::endl;

    for (auto seed: seeds) {

        GTP groundTruthProblem = gaussianGroundTruthProblem<long>(numberOfElements, alpha, sigma0, sigma1,
                                                                  multiplier, seed);

        LinearOrderingProblem<long> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto performance = testTranspositionCondition(problem);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        out << seed << ","
            << numberOfElements << ","
            << alpha << ","
            << sigma0 << ","
            << sigma1 << ","
            << performance.initialVars << ","
            << performance.remainingVars << ","
            << performance.eliminatedVarsFraction << ","
            << duration
            << std::endl;
    }
}

void testGroupingSinglePairCondition(
        std::ofstream &out,
        size_t numberOfElements,
        double alpha,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
) {
    std::vector<size_t> seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);

    std::cout << "Collecting statistics for numberOfElements=" << numberOfElements << " and " << " alpha=" << alpha
              << std::endl;

    for (auto seed: seeds) {

        GTP groundTruthProblem = gaussianGroundTruthProblem<long>(
                numberOfElements,
                alpha,
                sigma0,
                sigma1,
                multiplier,
                seed);

        LinearOrderingProblem<long> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto performance = testGroupingSinglePairCondition(problem);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        out << seed << ","
            << numberOfElements << ","
            << alpha << ","
            << sigma0 << ","
            << sigma1 << ","
            << performance.initialVars << ","
            << performance.remainingVars << ","
            << performance.eliminatedVarsFraction << ","
            << duration
            << std::endl;
    }
}

void runDifferentAlphas(
        std::filesystem::path const & outDirectory,
        std::vector<size_t> const & numbersOfElements = {20, 40, 60, 80, 100, 140, 200}
        ) {
    if (!std::filesystem::is_directory(outDirectory)) {
        std::filesystem::create_directories(outDirectory);
    }
    size_t seeds = 20;
    double sigma0 = 0.2;
    double sigma1 = 0.2;
    double multiplier = 1000;

    std::vector<double> alphas;
    size_t steps = 101;
    for (size_t step = 0; step < steps; step++) {
        alphas.push_back(static_cast<double>(step) / (static_cast<double>(steps) - 1));
    }


    for (auto numberOfElements: numbersOfElements) {

        std::filesystem::path directory = outDirectory / "criterion-heading-subproblem";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "numElements=" << numberOfElements << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);


        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto alpha: alphas) {
            testHeadingSubproblemCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }

        file.close();
    }

    for (auto numberOfElements: numbersOfElements) {

        std::filesystem::path directory = outDirectory / "criterion-heading-element";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "numElements=" << numberOfElements << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto alpha: alphas) {
            testHeadingElementCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
    }

    for (auto numberOfElements: numbersOfElements) {
        std::filesystem::path directory = outDirectory / "criterion-trailing-element";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "numElements=" << numberOfElements << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto alpha: alphas){
            testTrailingElementCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
    }

    for (auto numberOfElements: numbersOfElements) {
        std::filesystem::path directory = outDirectory / "criterion-heading-pair";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "numElements=" << numberOfElements << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto alpha: alphas){
            testHeadingPairCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
    }

    for (auto numberOfElements: numbersOfElements) {
        std::filesystem::path directory = outDirectory / "criterion-transposition";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "numElements=" << numberOfElements << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto alpha: alphas){
            testTranspositionCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
    }

    for (auto numberOfElements: numbersOfElements) {
        std::filesystem::path directory = outDirectory / "criterion-grouping-pair";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "numElements=" << numberOfElements << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto alpha: alphas){
            testGroupingSinglePairCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
    }

}


void runDifferentNumberOfElements(
        std::filesystem::path const & outDirectory,
        size_t max = 100,
        std::vector<double> const & alphas = {0.4, 0.65, 0.7, 1.0}
        ) {
    if (!std::filesystem::is_directory(outDirectory)) {
        std::filesystem::create_directories(outDirectory);
    }

    size_t seeds = 20;
    double sigma0 = 0.2;
    double sigma1 = 0.2;
    double multiplier = 1000;

    size_t increment = 5;
    size_t steps = max / increment;

    std::vector<size_t > numbersOfElements;
    for (size_t step = 1; step < steps + 1; step++){
        numbersOfElements.push_back(step*increment);
    }


    for (auto alpha: alphas) {

        std::filesystem::path directory = outDirectory / "criterion-heading-subproblem";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);


        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto numberOfElements: numbersOfElements) {
            testHeadingSubproblemCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }

        file.close();
    }

    for (auto alpha: alphas) {

        std::filesystem::path directory = outDirectory / "criterion-heading-element";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto numberOfElements: numbersOfElements) {
            testHeadingElementCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
        file.close();
    }

    for (auto alpha: alphas) {
        std::filesystem::path directory = outDirectory / "criterion-trailing-element";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto numberOfElements: numbersOfElements){
            testTrailingElementCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
        file.close();
    }

    for (auto alpha: alphas) {
        std::filesystem::path directory = outDirectory / "criterion-heading-pair";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto numberOfElements: numbersOfElements){
            testHeadingPairCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
        file.close();
    }

    for (auto alpha: alphas) {
        std::filesystem::path directory = outDirectory / "criterion-trailing-pair";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto numberOfElements: numbersOfElements){
            testTrailingPairCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
        file.close();
    }

    for (auto alpha: alphas) {
        std::filesystem::path directory = outDirectory / "criterion-transposition";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto numberOfElements: numbersOfElements){
            testTranspositionCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
        file.close();
    }

    for (auto alpha: alphas) {
        std::filesystem::path directory = outDirectory / "criterion-grouping-pair";
        if (!std::filesystem::is_directory(directory)) {
            std::filesystem::create_directories(directory);
        }

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(directory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds"
             << std::endl;

        for (auto numberOfElements: numbersOfElements){
            testGroupingSinglePairCondition(
                    file,
                    numberOfElements,
                    alpha,
                    sigma0,
                    sigma1,
                    multiplier,
                    seeds
            );
        }
        file.close();
    }
}

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("synthetic-problem-individual-conditions");

    int max = 0;

    // subparser individual
    argparse::ArgumentParser numElements("n");
    numElements.add_description("Run different numberOfElements.");
    numElements.add_argument("--maxElements")
            .default_value(200)
            .store_into(max);
    numElements.add_argument("--alphas")
            .default_value(std::vector<double>{0.65, 0.7, 1.0})
            .nargs(0, 4)
            .scan<'g', double>();

    argparse::ArgumentParser alphas("a");
    alphas.add_description("Run different alphas.");
    alphas.add_argument("--ns")
            .default_value(std::vector<size_t>{20, 40, 60, 80, 100})
            .nargs(0, 10)
            .scan<'i', size_t>();

    program.add_subparser(numElements);
    program.add_subparser(alphas);

    try {
        program.parse_args(argc, argv);

        if (program.is_subcommand_used("n")){
            std::cout << "Running with multiple alphas" << std::endl;
            runDifferentNumberOfElements(
                    "../results/synthetic-individual/alphas",
                    static_cast<size_t>(max),
                    numElements.get<std::vector<double>>("--alphas")
            );
        } else if (program.is_subcommand_used("a")){
            std::cout << "Running with multiple numbers of elements." << std::endl;
            runDifferentAlphas(
                    "../results/synthetic-individual/ns",
                    alphas.get<std::vector<size_t>>("--ns")
            );
        } else {
            std::cout << program.help().str() << std::endl;
            return 0;
        }
    } catch (const std::exception & err){
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    return 0;
}
#include "dscc/lolip2010/lolip-reader.hxx"
#include "andres/graph/linear-ordering/persistency.hxx"
#include <filesystem>
#include <chrono>
#include <andres/ilp/linear-ordering-ilp.hxx>
#include <andres/graph/complete-digraph.hxx>
#include "andres/ilp/gurobi.hxx"
#include "argparse.hxx"

template<class T>
T base_name(T const &path, T const &delims = "/\\") {
    return path.substr(path.find_last_of(delims) + 1);
}

template<class T>
T remove_extension(T const &filename) {
    typename T::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != T::npos ? filename.substr(0, p) : filename;
}

template<class T>
std::string format(T const &number, std::streamsize precision) {
    std::ostringstream ss;
    ss << std::fixed;
    ss.precision(precision);
    ss << number;
    return ss.str();
}

struct PersistencyPerformance {
public:
    PersistencyPerformance() = default;

    size_t initialVars;
    size_t remainingVars;
    double eliminatedVarsFraction;
};

template<class R>
PersistencyPerformance testIndependentSubProblems(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<long> persistency(currentProblem);

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
PersistencyPerformance testHeadingElement(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<long> persistency(currentProblem);

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
PersistencyPerformance testTrailingElement(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<long> persistency(currentProblem);

        auto foundTrailingSubElement = persistency.findTrailingElement();
        if (foundTrailingSubElement.first) {
            auto index = foundTrailingSubElement.second;
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
PersistencyPerformance testHeadingPair(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<long> persistency(currentProblem);

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
PersistencyPerformance testTrailingPair(LinearOrderingProblem<R> const &problem) {


    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

    std::queue<LinearOrderingProblem<R>> problemQueue;
    problemQueue.push(problem);

    std::vector<LinearOrderingProblem<R>> irreducibleProblems;

    while (!problemQueue.empty()) {
        auto currentProblem = problemQueue.front();
        problemQueue.pop();

        Persistency<long> persistency(currentProblem);

        auto foundTrailingPair = persistency.findTrailingPair();
        if (foundTrailingPair.first) {
            size_t a = foundTrailingPair.second.first;
            size_t b = foundTrailingPair.second.second;
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
PersistencyPerformance testTransposition(LinearOrderingProblem<R> const &problem) {

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements * (numberOfElements - 1);

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

    Persistency<R> persistency(problem);
    persistency.checkGroupingPropositionSinglePair();

    size_t eliminatedVars = 2 * persistency.oneEdges().size();

    size_t remainingVariables = problem.numberOfElements() * (problem.numberOfElements() - 1) - eliminatedVars;

    double eliminatedVariableFraction =
            static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}


void testLolipProblemIndependentSubProblems(
        std::filesystem::path const &pathToFile,
        std::ofstream &out
) {
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    auto start = std::chrono::high_resolution_clock::now();

    auto performance = testIndependentSubProblems(problem);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Independent Subproblems: " << performance.eliminatedVarsFraction << std::endl;

    out << problemName << ","
        << problem.numberOfElements() << ","
        << performance.remainingVars << ","
        << performance.initialVars << ","
        << format(performance.eliminatedVarsFraction, 4) << ","
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << std::endl;
}

void testLolipProblemHeadingElement(
        std::filesystem::path const &pathToFile,
        std::ofstream &out
) {
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    auto start = std::chrono::high_resolution_clock::now();

    auto performance = testHeadingElement(problem);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Heading Element: " << performance.eliminatedVarsFraction << std::endl;

    out << problemName << ","
        << problem.numberOfElements() << ","
        << performance.remainingVars << ","
        << performance.initialVars << ","
        << format(performance.eliminatedVarsFraction, 4) << ","
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << std::endl;
}

void testLolipProblemTrailingElement(
        std::filesystem::path const &pathToFile,
        std::ofstream &out
) {
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    auto start = std::chrono::high_resolution_clock::now();

    auto performance = testTrailingElement(problem);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Trailing Element: " << performance.eliminatedVarsFraction << std::endl;

    out << problemName << ","
        << problem.numberOfElements() << ","
        << performance.remainingVars << ","
        << performance.initialVars << ","
        << format(performance.eliminatedVarsFraction, 4) << ","
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << std::endl;
}

void testLolipProblemHeadingPairElement(
        std::filesystem::path const &pathToFile,
        std::ofstream &out
) {
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    auto start = std::chrono::high_resolution_clock::now();

    auto performance = testHeadingPair(problem);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Heading Pair: " << performance.eliminatedVarsFraction << std::endl;

    out << problemName << ","
        << problem.numberOfElements() << ","
        << performance.remainingVars << ","
        << performance.initialVars << ","
        << format(performance.eliminatedVarsFraction, 4) << ","
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << std::endl;
}

void testLolipProblemTrailingPairElement(
        std::filesystem::path const &pathToFile,
        std::ofstream &out
) {
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    auto start = std::chrono::high_resolution_clock::now();

    auto performance = testTrailingPair(problem);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Trailing Pair: " << performance.eliminatedVarsFraction << std::endl;

    out << problemName << ","
        << problem.numberOfElements() << ","
        << performance.remainingVars << ","
        << performance.initialVars << ","
        << format(performance.eliminatedVarsFraction, 4) << ","
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << std::endl;
}

void testLolipProblemTransposition(
        std::filesystem::path const &pathToFile,
        std::ofstream &out
) {
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    auto start = std::chrono::high_resolution_clock::now();

    auto performance = testTransposition(problem);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Transposition: " << performance.eliminatedVarsFraction << std::endl;

    out << problemName << ","
        << problem.numberOfElements() << ","
        << performance.remainingVars << ","
        << performance.initialVars << ","
        << format(performance.eliminatedVarsFraction, 4) << ","
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << std::endl;
}

void testLolipProblemGroupingPairs(
        std::filesystem::path const &pathToFile,
        std::ofstream &out
) {
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    auto start = std::chrono::high_resolution_clock::now();

    auto performance = testGroupingSinglePairCondition(problem);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Grouping-Pair: " << performance.eliminatedVarsFraction << std::endl;

    out << problemName << ","
        << problem.numberOfElements() << ","
        << performance.remainingVars << ","
        << performance.initialVars << ","
        << format(performance.eliminatedVarsFraction, 4) << ","
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << std::endl;
}


std::ofstream getFile(
        std::filesystem::path const & directory,
        std::filesystem::path const & outFileName
        ){
    if (!std::filesystem::is_directory(directory)) {
        std::filesystem::create_directories(directory);
    }

    std::filesystem::path outFileNameIndependentSubProblems = directory / outFileName;
    std::ofstream file(outFileNameIndependentSubProblems);
    file << "problem,n,remainingVars,initialVars,fractionEliminated,durationNanoseconds" << std::endl;

    return file;
}


void lolibProblemDirectory(
        std::filesystem::path const &path,
        std::filesystem::path const &outDirectory
) {

    if (!std::filesystem::is_directory(outDirectory)) {
        std::filesystem::create_directories(outDirectory);
    }

    auto problemClass = base_name<std::string>(path);

    auto outFileName = problemClass + ".csv";

    auto ofIndependenSubProblems = getFile(outDirectory / "criterion-heading-subproblem", outFileName);
    auto ofHeadingElement = getFile(outDirectory / "criterion-heading-element", outFileName);
    auto ofTrailingElement = getFile(outDirectory / "criterion-trailing-element", outFileName);
    auto ofHeadingPair = getFile(outDirectory / "criterion-heading-pair", outFileName);
    auto ofTransposition = getFile(outDirectory / "criterion-transposition", outFileName);
    auto ofGroupingSinglePair = getFile(outDirectory / "criterion-grouping-pair", outFileName);
    auto ofTrailingPair = getFile(outDirectory / "criterion-trailing-pair", outFileName);

    for (const auto &entry: std::filesystem::directory_iterator(path)) {
        std::cout << base_name(entry.path().string()) << std::endl << std::string (100, '=') << std::endl;
        testLolipProblemIndependentSubProblems(entry.path(), ofIndependenSubProblems);
        testLolipProblemHeadingElement(entry.path(), ofHeadingElement);
        testLolipProblemTrailingElement(entry.path(), ofTrailingElement);
        testLolipProblemHeadingPairElement(entry.path(), ofHeadingPair);
        testLolipProblemTrailingPairElement(entry.path(), ofTrailingPair);
        testLolipProblemTransposition(entry.path(), ofTransposition);
        testLolipProblemGroupingPairs(entry.path(), ofGroupingSinglePair);

        std::cout << std::endl;
    }
}

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("lolib-individual");

    program.add_description("Run individual persistency conditions.");
    program.add_argument("--problem")
            .default_value(std::string{"IO"})
            .choices("IO", "SGB", "Spec");


    try {
        program.parse_args(argc, argv);

        std::cout << "Running individual persistency conditions." << std::endl;
        std::cout << "Problem Class: " << program.get<std::string>("--problem") << std::endl;
        lolibProblemDirectory(
                "../data/lolib_2010/" + program.get<std::string>("--problem"),
                "../results/lolib/individual/" + program.get<std::string>("--problem")
        );
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }


    return 0;
}
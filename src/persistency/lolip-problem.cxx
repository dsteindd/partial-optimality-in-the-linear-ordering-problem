#include "dscc/lolip2010/lolip-reader.hxx"
#include "andres/graph/linear-ordering/persistency.hxx"
#include <filesystem>
#include <chrono>
#include <andres/ilp/linear-ordering-ilp.hxx>
#include <andres/graph/complete-digraph.hxx>
#include "andres/ilp/gurobi.hxx"
#include "argparse.hxx"

template<class T>
T base_name(T const & path, T const & delims = "/\\")
{
    return path.substr(path.find_last_of(delims) + 1);
}
template<class T>
T remove_extension(T const & filename)
{
    typename T::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != T::npos ? filename.substr(0, p) : filename;
}

template<class T>
std::string format(T const & number, std::streamsize precision){
    std::ostringstream ss;
    ss << std::fixed;
    ss.precision(precision);
    ss << number;
    return ss.str();
}

struct PersistencyPerformance{
public:
    PersistencyPerformance()= default;

    size_t initialVars;
    size_t remainingVars;
    double eliminatedVarsFraction;
};

template<class T>
struct PersistencySolution{
public:
    typedef std::vector<LinearOrderingProblem<T>> SolutionContainer;
    typedef std::pair<size_t, size_t > Edge;
    typedef std::set<Edge > EdgeSet;
    typedef std::vector<EdgeSet> EdgeSetContainer;

    PersistencySolution() = default;

    SolutionContainer subproblems;
    EdgeSetContainer oneEdgeSets;
};

template<class T>
struct SolverPerformance{
public:
    SolverPerformance() = default;

    double solutionValue;
    long solveTime;
};

template<class R>
PersistencyPerformance testPersistencyForProblem(LinearOrderingProblem<R> const & problem){

    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements*(numberOfElements - 1);

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
            auto foundHeadingSubElement = persistency.findHeadingElement();
            if (foundHeadingSubElement.first){
                auto index = foundHeadingSubElement.second;
                std::vector<size_t > otherIndices;
                for (size_t i = 0; i < currentProblem.numberOfElements(); ++i){
                    if (i != index){
                        otherIndices.push_back(i);
                    }
                }
                problemQueue.push(project(currentProblem, {index}));
                problemQueue.push(project(currentProblem, otherIndices));
            } else {
                auto foundTrailingSubElement = persistency.findTrailingElement();
                if (foundTrailingSubElement.first){
                    auto index = foundTrailingSubElement.second;
                    std::vector<size_t > otherIndices;
                    for (size_t i = 0; i < currentProblem.numberOfElements(); ++i){
                        if (i != index){
                            otherIndices.push_back(i);
                        }
                    }
                    problemQueue.push(project(currentProblem, {index}));
                    problemQueue.push(project(currentProblem, otherIndices));
                } else {
                    auto foundHeadingPair = persistency.findHeadingPair();
                    if (foundHeadingPair.first){
                        size_t a = foundHeadingPair.second.first;
                        size_t b = foundHeadingPair.second.second;
                        std::vector<size_t > otherIndices;
                        for (size_t i = 0; i < currentProblem.numberOfElements(); ++i){
                            if (i != a && i != b){
                                otherIndices.push_back(i);
                            }
                        }
                        problemQueue.push(project(currentProblem, {a, b}));
                        problemQueue.push(project(currentProblem, otherIndices));
                    } else {
                        auto foundTrailingPair = persistency.findTrailingPair();
                        if (foundTrailingPair.first){
                            size_t a = foundTrailingPair.second.first;
                            size_t b = foundTrailingPair.second.second;
                            std::vector<size_t > otherIndices;
                            for (size_t i = 0; i < currentProblem.numberOfElements(); ++i){
                                if (i != a && i != b){
                                    otherIndices.push_back(i);
                                }
                            }
                            problemQueue.push(project(currentProblem, {a, b}));
                            problemQueue.push(project(currentProblem, otherIndices));
                        } else {
                            irreducibleProblems.push_back(currentProblem);
                        }
                    }
                }
            }
        }
    }

    size_t remainingVariables = 0;
    for (auto &irreducibleProblem: irreducibleProblems) {
        Persistency<long> persistency(irreducibleProblem);

        persistency.checkTranspositionProposition();

        persistency.checkGroupingPropositionSinglePair();

        size_t additionalEliminatedVariables = persistency.oneEdges().size()*2;

        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1) - additionalEliminatedVariables;
    }

    double eliminatedVariableFraction = static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}

template<class R>
std::pair<PersistencyPerformance, PersistencySolution<R>> applyPersistencyForPreprocessing(LinearOrderingProblem<R> const & problem){
    typedef std::pair<size_t, size_t > Edge;
    typedef std::set<Edge > EdgeSet;
    typedef std::vector<EdgeSet> EdgeSetContainer;


    size_t numberOfElements = problem.numberOfElements();
    size_t initialVars = numberOfElements*(numberOfElements - 1);

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
            auto foundHeadingSubElement = persistency.findHeadingElement();
            if (foundHeadingSubElement.first){
                auto index = foundHeadingSubElement.second;
                std::vector<size_t > otherIndices;
                for (size_t i = 0; i < currentProblem.numberOfElements(); ++i){
                    if (i != index){
                        otherIndices.push_back(i);
                    }
                }
                problemQueue.push(project(currentProblem, {index}));
                problemQueue.push(project(currentProblem, otherIndices));
            } else {
                auto foundTrailingSubElement = persistency.findTrailingElement();
                if (foundTrailingSubElement.first){
                    auto index = foundTrailingSubElement.second;
                    std::vector<size_t > otherIndices;
                    for (size_t i = 0; i < currentProblem.numberOfElements(); ++i){
                        if (i != index){
                            otherIndices.push_back(i);
                        }
                    }
                    problemQueue.push(project(currentProblem, {index}));
                    problemQueue.push(project(currentProblem, otherIndices));
                } else {
                    auto foundHeadingPair = persistency.findHeadingPair();
                    if (foundHeadingPair.first){
                        size_t a = foundHeadingPair.second.first;
                        size_t b = foundHeadingPair.second.second;
                        std::vector<size_t > otherIndices;
                        for (size_t i = 0; i < currentProblem.numberOfElements(); ++i){
                            if (i != a && i != b){
                                otherIndices.push_back(i);
                            }
                        }
                        problemQueue.push(project(currentProblem, {a, b}));
                        problemQueue.push(project(currentProblem, otherIndices));
                    } else {
                        auto foundTrailingPair = persistency.findTrailingPair();
                        if (foundTrailingPair.first){
                            size_t a = foundTrailingPair.second.first;
                            size_t b = foundTrailingPair.second.second;
                            std::vector<size_t > otherIndices;
                            for (size_t i = 0; i < currentProblem.numberOfElements(); ++i){
                                if (i != a && i != b){
                                    otherIndices.push_back(i);
                                }
                            }
                            problemQueue.push(project(currentProblem, {a, b}));
                            problemQueue.push(project(currentProblem, otherIndices));
                        } else {
                            irreducibleProblems.push_back(currentProblem);
                        }
                    }
                }
            }
        }
    }

    EdgeSetContainer oneEdgeSets;

    size_t remainingVariables = 0;
    for (auto &irreducibleProblem: irreducibleProblems) {
        Persistency<long> persistency(irreducibleProblem);

        persistency.checkTranspositionProposition();

        persistency.checkGroupingPropositionSinglePair();

        oneEdgeSets.push_back(persistency.oneEdges());

        size_t additionalEliminatedVariables = persistency.oneEdges().size()*2;

        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1) - additionalEliminatedVariables;
    }

    double eliminatedVariableFraction = static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {
        {initialVars, remainingVariables, eliminatedVariableFraction},
        {irreducibleProblems, oneEdgeSets}
        };
}

template<class R>
SolverPerformance<R> solveOptimalWithConstraints(LinearOrderingProblem<R> const & problem, std::set<std::pair<size_t , size_t>> const & oneEdgeSet){
    typedef andres::graph::CompleteDigraph Graph;

    Graph graph(problem.numberOfElements());
    std::vector<double> edgeWeights(graph.numberOfEdges());
    std::vector<size_t > labels(graph.numberOfEdges(), 0);

    for (size_t i = 0; i < problem.numberOfElements(); ++i){
        for (size_t j = 0; j < problem.numberOfElements(); ++j){
            if (i == j)
                continue;

            auto edge = graph.findEdge(i, j);

            edgeWeights[edge.second] = problem.cost(i, j);
        }
    }

    auto start = std::chrono::high_resolution_clock ::now();

    if (problem.numberOfElements()*(problem.numberOfElements() - 1) != 2*oneEdgeSet.size()){
        andres::graph::ordering::ilp<andres::ilp::Gurobi>(graph, edgeWeights, labels, labels, oneEdgeSet);
    }

    auto end = std::chrono::high_resolution_clock ::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    double solverValue = 0;
    for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge){
        solverValue += static_cast<double>(labels[edge])*edgeWeights[edge];
    }
    return {solverValue, duration};
}


void testPreprocessingWithOptimalSolver(
        std::filesystem::path const & pathToFile,
        std::ofstream & outPreprocessing
        ) {
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    std::cout << "Applying Persistency for problem " << problemName << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto pair = applyPersistencyForPreprocessing(problem);

    auto end = std::chrono::high_resolution_clock::now();

    auto performance = pair.first;
    auto subInstancesWithConstraints = pair.second;
    auto subproblems = subInstancesWithConstraints.subproblems;
    auto oneEdgeSets = subInstancesWithConstraints.oneEdgeSets;


    std::cout << problemName << ": " << performance.eliminatedVarsFraction << std::endl;

    outPreprocessing << problemName << ","
                     << problem.numberOfElements() << ","
                     << performance.remainingVars << ","
                     << performance.initialVars << ","
                     << format(performance.eliminatedVarsFraction, 4) << ","
                     << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << ",";

    std::cout << "Solving subproblems for problem " << problemName << std::endl;

    long duration = 0;

    for (size_t numberOfProblem = 0; numberOfProblem < subproblems.size(); ++numberOfProblem) {
        auto performance = solveOptimalWithConstraints(subproblems[numberOfProblem], oneEdgeSets[numberOfProblem]);
        duration += performance.solveTime;
    }

    outPreprocessing << duration << ",";

    // solve initial problem optimaly

    auto performanceOptimalInitial = solveOptimalWithConstraints(problem, {});

    outPreprocessing << performanceOptimalInitial.solveTime << std::endl;
}


void testLolipProblem(
        std::filesystem::path const & pathToFile,
        std::ofstream & out
){
    auto problemName = base_name<std::string>(pathToFile);

    LinearOrderingProblem<long> problem = readDataset(pathToFile);

    auto start = std::chrono::high_resolution_clock::now();

    auto performance = testPersistencyForProblem(problem);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << problemName << ": " << performance.eliminatedVarsFraction << std::endl;

    out << problemName << ","
        << problem.numberOfElements() << ","
        << performance.remainingVars << ","
        << performance.initialVars << ","
        << format(performance.eliminatedVarsFraction, 4) << ","
        << std::chrono::duration_cast<std::chrono::nanoseconds >(end - start).count()
        << std::endl;
}

void lolibProblemDirectory(
        std::filesystem::path const & directory,
        std::filesystem::path const & outDirectory,
        std::string const & outFile
        ){
    if (!std::filesystem::is_directory(outDirectory)){
        std::filesystem::create_directories(outDirectory);
    }

    std::filesystem::path outFileName = outDirectory / outFile;
    std::ofstream file(outFileName);

    file << "problem,n,remainingVars,initialVars,fractionEliminated,durationNanoseconds" << std::endl;

    for (const auto & entry : std::filesystem::directory_iterator(directory)){
//        std::cout << entry.path() << std::endl;
        testLolipProblem(entry.path(), file);
    }
}

void lolibProblemDirectoryWithPreprocessing(
        std::filesystem::path const & directory,
        std::filesystem::path const & outDirectory,
        std::string const & outFile
){
    if (!std::filesystem::is_directory(outDirectory)){
        std::filesystem::create_directories(outDirectory);
    }

    std::filesystem::path outFileName = outDirectory / outFile;
    std::ofstream file(outFileName);

    file << "problem,n,remainingVars,initialVars,fractionEliminated,durationPreprocessingNanoseconds,optimalSolverRuntimeNanoseconds,optimalSolutionOnlyDuration" << std::endl;

    for (const auto & entry : std::filesystem::directory_iterator(directory)){
//        std::cout << entry.path() << std::endl;
        testPreprocessingWithOptimalSolver(entry.path(), file);
    }
}

void lolibDirectories(
        std::filesystem::path const & path,
        std::filesystem::path const & outDirectory
        ){

    for (const auto & entry : std::filesystem::directory_iterator(path)){
        if (entry.is_directory()){
            lolibProblemDirectory(
                    entry.path(),
                    outDirectory,
                    base_name<std::string>(entry.path())
                    );
        }

    }
}

int main(int argc, char* argv[]){

    argparse::ArgumentParser program("lolib");

    int max = 0;

    // subparser individual
    argparse::ArgumentParser persistencyOnly("p");
    persistencyOnly.add_description("Run only persistency");
    persistencyOnly.add_argument("--problem")
            .default_value(std::string{"IO"})
            .choices("IO", "SGB", "Spec");

    argparse::ArgumentParser persistencyAndIlp("pi");
    persistencyAndIlp.add_description("Run persistency and ILP Solver.");
    persistencyAndIlp.add_argument("--problem")
            .default_value(std::string{"IO"})
            .choices("IO", "Spec");




    program.add_subparser(persistencyOnly);
    program.add_subparser(persistencyAndIlp);

    try {
        program.parse_args(argc, argv);

        if (program.is_subcommand_used("p")){
            std::cout << "Running persistency only." << std::endl;
            std::cout << "Problem Class: " << persistencyOnly.get<std::string>("--problem") << std::endl;
            lolibProblemDirectory(
                    "../data/lolib_2010/" + persistencyOnly.get<std::string>("--problem"),
                    "../results/lolib/persistency",
                    persistencyOnly.get<std::string>("--problem") + ".csv"
            );
        } else if (program.is_subcommand_used("pi")){
            std::cout << "Running persistency as presolve and then ILP solver." << std::endl;
            std::cout << "Problem Class: " << persistencyAndIlp.get<std::string>("--problem") << std::endl;
            lolibProblemDirectoryWithPreprocessing(
                    "../data/lolib_2010/" + persistencyAndIlp.get<std::string>("--problem"),
                    "../results/lolib/persistency-and-ilp",
                    persistencyAndIlp.get<std::string>("--problem") + ".csv"
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
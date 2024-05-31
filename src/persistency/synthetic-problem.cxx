#include "dscc/ground-truth/ground-truth-order-problem.hxx"
#include "andres/graph/linear-ordering/persistency.hxx"
#include "iostream"
#include "fstream"
#include "sstream"
#include "chrono"
#include "filesystem"
#include <andres/ilp/linear-ordering-ilp.hxx>
#include <andres/graph/complete-digraph.hxx>
#include "andres/ilp/gurobi.hxx"
#include "argparse.hxx"

struct PersistencyPerformance{
public:
    PersistencyPerformance()= default;

    size_t initialVars;
    size_t remainingVars;
    double eliminatedVarsFraction;
};

typedef GroundTruthOrderProblem<long, std::normal_distribution<double>> GTP;

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
        Persistency<R> persistency(irreducibleProblem);

        persistency.checkTranspositionProposition();

//        // todo: seems to be too slow
        persistency.checkGroupingPropositionSinglePair();

        // actually double the oneEdges are eliminated because if xab = 1 then xba = 0
        size_t additionalEliminatedVariables = persistency.oneEdges().size()*2;

        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1) - additionalEliminatedVariables;
    }

    double eliminatedVariableFraction = static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {initialVars, remainingVariables, eliminatedVariableFraction};
}

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

template<class R>
std::pair<PersistencyPerformance, PersistencySolution<R>> applyPersistencyForProblem(LinearOrderingProblem<R> const & problem){
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
        Persistency<R> persistency(irreducibleProblem);

        persistency.checkTranspositionProposition();

//        // todo: seems to be too slow
        persistency.checkGroupingPropositionSinglePair();

        // actually double the oneEdges are eliminated because if xab = 1 then xba = 0
        size_t additionalEliminatedVariables = persistency.oneEdges().size()*2;
        oneEdgeSets.push_back(persistency.oneEdges());

        remainingVariables += irreducibleProblem.numberOfElements() * (irreducibleProblem.numberOfElements() - 1) - additionalEliminatedVariables;
    }

    double eliminatedVariableFraction = static_cast<double>(initialVars - remainingVariables) / static_cast<double>(initialVars);

    return {
        {initialVars, remainingVariables, eliminatedVariableFraction},
        {irreducibleProblems, oneEdgeSets}
    };
}

void testGroundTruthProblem(
        std::ofstream & out,
        size_t numberOfElements,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double alpha = 0.2,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
){
    std::vector<size_t > seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);


    std::cout << "Collecting statistics for numberOfElements " << numberOfElements << " and alpha=" << alpha << std::endl;

    for (auto seed : seeds){

        GTP groundTruthProblem = gaussianGroundTruthProblem<long>(numberOfElements, alpha, sigma0, sigma1, multiplier, seed);

        LinearOrderingProblem<long> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto performance = testPersistencyForProblem(problem);
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


template<class T>
struct SolverPerformance{
public:
    SolverPerformance() = default;

    double solutionValue;
    long solveTime;
};

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

void testGroundTruthProblemWithPreprocessing(
        std::ofstream & out,
        size_t numberOfElements,
        double sigma0 = 0.2,
        double sigma1 = 0.4,
        double alpha = 0.2,
        double multiplier = 1000,
        size_t numberOfSeeds = 50
){
    std::vector<size_t > seeds(numberOfSeeds);
    std::iota(seeds.begin(), seeds.end(), 42);


    std::cout << "Collecting statistics for numberOfElements " << numberOfElements << " and alpha=" << alpha << std::endl;

    for (auto seed : seeds){

        GTP groundTruthProblem = gaussianGroundTruthProblem<long>(numberOfElements, alpha, sigma0, sigma1, multiplier, seed);

        LinearOrderingProblem<long> problem(groundTruthProblem);

        auto start = std::chrono::high_resolution_clock::now();
        auto pair = applyPersistencyForProblem(problem);
        auto end = std::chrono::high_resolution_clock::now();

        auto performance = pair.first;
        auto subproblemsWithConstraints = pair.second;
        auto subproblems = subproblemsWithConstraints.subproblems;
        auto oneEdgeSets = subproblemsWithConstraints.oneEdgeSets;

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        long durationOptimalAfterPreprocessing = 0;
        for (size_t numProblem = 0; numProblem < subproblems.size(); ++numProblem){
            auto subPerformance = solveOptimalWithConstraints(subproblems[numProblem], oneEdgeSets[numProblem]);
            durationOptimalAfterPreprocessing += subPerformance.solveTime;
        }

        auto initialPerformance = solveOptimalWithConstraints(problem, {});

        out << seed << ","
            << numberOfElements << ","
            << alpha << ","
            << sigma0 << ","
            << sigma1 << ","
            << performance.initialVars << ","
            << performance.remainingVars << ","
            << performance.eliminatedVarsFraction << ","
            << duration << ","
            << durationOptimalAfterPreprocessing << ","
            << initialPerformance.solveTime
            << std::endl;
    }
}

void runDifferentAlphas(
        std::filesystem::path const & outDirectory,
        std::vector<size_t> const & numbersOfElements = {20, 40, 60, 80, 100, 150, 200}
        ){
    if (!std::filesystem::is_directory(outDirectory)){
        std::filesystem::create_directories(outDirectory);
    }

    size_t seeds = 20;
    double sigma0 = 0.2;
    double sigma1 = 0.2;
    double multiplier = 1000;

    std::vector<double> alphas;
    size_t steps = 101;
    for (size_t step = 0; step < steps; step++){
        alphas.push_back(static_cast<double>(step)/static_cast<double>(steps - 1));
    }

    for (auto numberOfElements: numbersOfElements){
        std::stringstream ss;
        ss << "numElements=" << numberOfElements << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(outDirectory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds" << std::endl;


        for (auto alpha: alphas){
            testGroundTruthProblem(
                    file,
                    numberOfElements,
                    sigma0,
                    sigma1,
                    alpha,
                    multiplier,
                    seeds
            );
        }
    }
}


void runDifferentAlphasWithPreprocessingAndOptimalSolver(
        std::filesystem::path const & outDirectory,
        std::vector<size_t> const & numbersOfElements = {20, 40}
        ){
    if (!std::filesystem::is_directory(outDirectory)){
        std::filesystem::create_directories(outDirectory);
    }

    size_t seeds = 20;
    double sigma0 = 0.2;
    double sigma1 = 0.2;
    double multiplier = 1000;

    std::vector<double> alphas;
    size_t steps = 51;
    for (size_t step = 0; step < steps; step++){
        alphas.push_back(static_cast<double>(step)/static_cast<double>(steps - 1));
    }

    for (auto numberOfElements: numbersOfElements){
        std::stringstream ss;
        ss << "numElements=" << numberOfElements << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(outDirectory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated"
             << ",durationPreprocessingNanoseconds,optimalSolverRuntimeNanoseconds,optimalSolutionOnlyDuration" << std::endl;


        for (auto alpha: alphas){
            testGroundTruthProblemWithPreprocessing(
                    file,
                    numberOfElements,
                    sigma0,
                    sigma1,
                    alpha,
                    multiplier,
                    seeds
            );
        }
    }
}

void runDifferentNumberOfElementsWithPreprocessingAndOptimalSolver(
        std::filesystem::path const & outDirectory,
        size_t const & maxElements = 100,
        std::vector<double> const & alphas = {1.0}
        ){
    if (!std::filesystem::is_directory(outDirectory)){
        std::filesystem::create_directories(outDirectory);
    }

    size_t seeds = 20;
    double sigma0 = 0.2;
    double sigma1 = 0.2;
    double multiplier = 1000;

    size_t increment = 5;
    size_t steps = maxElements / increment;

    std::vector<size_t > numbersOfElements;
    for (size_t step = 1; step < steps + 1; step++){
        numbersOfElements.push_back(step*increment);
    }

    for (auto alpha: alphas){

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(outDirectory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated"
            << ",durationPreprocessingNanoseconds,optimalSolverRuntimeNanoseconds,optimalSolutionOnlyDuration" << std::endl;

        for (auto numberOfElements: numbersOfElements)
            testGroundTruthProblemWithPreprocessing(
                    file,
                    numberOfElements,
                    sigma0,
                    sigma1,
                    alpha,
                    multiplier,
                    seeds
            );
    }
}

void runDifferentNumberOfElements(
        std::filesystem::path const & outDirectory,
        size_t const & maxElements = 200,
        std::vector<double> const & alphas = {0.65, 0.7, 1.0}
){

    if (!std::filesystem::is_directory(outDirectory)){
        std::filesystem::create_directories(outDirectory);
    }

    size_t seeds = 20;
    double sigma0 = 0.2;
    double sigma1 = 0.2;
    double multiplier = 1000;

    size_t increment = 5;
    size_t steps = maxElements / increment;

    std::vector<size_t > numbersOfElements;
    for (size_t step = 1; step < steps + 1; step++){
        numbersOfElements.push_back(step*increment);
    }

    for (auto alpha: alphas){

        std::stringstream ss;
        ss << "alpha=" << alpha << ".csv";
        std::string fileName = ss.str();
        std::ofstream file(outDirectory / fileName);

        file << "seed,n,alpha,sigma0,sigma1,initialVariables,remainingVars,fractionEliminated,durationNanoseconds" << std::endl;

        for (auto numberOfElements: numbersOfElements)
            testGroundTruthProblem(
                    file,
                    numberOfElements,
                    sigma0,
                    sigma1,
                    alpha,
                    multiplier,
                    seeds
            );
    }
}

int main(int argc, char *argv[]){
    argparse::ArgumentParser program("synthetic-problem");

    int max = 0;

    // subparser individual
    argparse::ArgumentParser numElements("n");
    numElements.add_description("Run different numberOfElements. Optionally set --maxElements");
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

    argparse::ArgumentParser ilp("presolve-and-ilp-n");
    ilp.add_description("Run with Presolve and ILP solver at the end.");
    ilp.add_argument("--maxElements")
        .default_value(100)
        .store_into(max);
    ilp.add_argument("--alphas")
            .default_value(std::vector<double>{0.65, 0.7, 1.0})
            .nargs(0, 4)
            .scan<'g', double>();

    argparse::ArgumentParser ilpAlphas("presolve-and-ilp-a");
    ilpAlphas.add_description("Run with Presolve and ILP solver at the end.");
    ilpAlphas.add_argument("--ns")
            .default_value(std::vector<size_t>{20, 40})
            .nargs(0, 10)
            .scan<'i', size_t>();

    program.add_subparser(numElements);
    program.add_subparser(alphas);
    program.add_subparser(ilp);
    program.add_subparser(ilpAlphas);

    try {
        program.parse_args(argc, argv);

        if (program.is_subcommand_used("n")){
            std::cout << "Running with multiple alphas" << std::endl;
            runDifferentNumberOfElements(
                    "../results/synthetic/alphas",
                    static_cast<size_t>(max),
                    numElements.get<std::vector<double>>("--alphas")
            );
        } else if (program.is_subcommand_used("a")){
            std::cout << "Running with multiple numbers of elements." << std::endl;
            runDifferentAlphas(
                    "../results/synthetic/ns",
                    alphas.get<std::vector<size_t>>("--ns")
            );
        } else if (program.is_subcommand_used("presolve-and-ilp-n")){
            std::cout << "Running with different number of elements and presolve plus ilp" << std::endl;
            runDifferentNumberOfElementsWithPreprocessingAndOptimalSolver(
                    "../results/synthetic/ilp-solver-alphas",
                    static_cast<size_t>(max),
                    ilp.get<std::vector<double>>("--alphas")
            );
        } else if (program.is_subcommand_used("presolve-and-ilp-a")){
            std::cout << "Running with different number of elements and presolve plus ilp" << std::endl;

            runDifferentAlphasWithPreprocessingAndOptimalSolver(
                    "../results/synthetic/ilp-solver-ns",
                    ilpAlphas.get<std::vector<size_t>>("--ns")
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
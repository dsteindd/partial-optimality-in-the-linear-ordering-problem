cmake_minimum_required(VERSION 3.26)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(partial-optimality-linear-ordering)
enable_testing()
include_directories(include)

find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

add_executable(lolip-problem src/persistency/lolip-problem.cxx)
target_include_directories(lolip-problem PUBLIC include)
target_link_libraries(lolip-problem ${GUROBI_LIBRARIES})

add_executable(lolip-problem-individual src/persistency/lolib-problem-individual-conditions.cxx)
target_include_directories(lolip-problem-individual PUBLIC include)
target_link_libraries(lolip-problem-individual ${GUROBI_LIBRARIES})

add_executable(synthetic-problem src/persistency/synthetic-problem.cxx)
target_include_directories(synthetic-problem PUBLIC include)
target_link_libraries(synthetic-problem ${GUROBI_LIBRARIES})

add_executable(synthetic-problem-individual-conditions src/persistency/synthetic-problem-individual-conditions.cxx)
target_include_directories(synthetic-problem-individual-conditions PUBLIC include)
target_link_libraries(synthetic-problem-individual-conditions ${GUROBI_LIBRARIES})

add_executable(lop-problem-tests src/test/lop-problem.cxx)
target_include_directories(lop-problem-tests PUBLIC include)
add_test(lop-problem-tests lop-problem-tests)

add_executable(lolip-reading-tests src/test/lolip-problem-reading.cxx)
target_include_directories(lolip-reading-tests PUBLIC include)
add_test(lolip-reading-tests lolip-reading-test)

add_executable(ground-truth-problem src/test/ground-truth-problem.cxx)
target_include_directories(ground-truth-problem PUBLIC include)
add_test(ground-truth-problem ground-truth-problem)

add_executable(ilp-solver src/test/ilp-solver.cxx)
target_include_directories(ilp-solver PUBLIC include)
target_link_libraries(ilp-solver ${GUROBI_LIBRARIES})
add_test(ilp-solver ilp-solver)

add_executable(digraph-tests src/test/digraph.cxx
        include/argparse.hxx)
target_include_directories(digraph-tests PUBLIC include)
add_test(digraph-tests digraph-tests)


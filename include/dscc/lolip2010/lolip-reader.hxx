//
// Created by dstein on 01.09.23.
//
#pragma once

#ifndef PARTIAL_OPTIMALITY_LINEAR_ORDERING_LOLIP_READER_HXX
#define PARTIAL_OPTIMALITY_LINEAR_ORDERING_LOLIP_READER_HXX

#include <sstream>
#include "andres/graph/linear-ordering/linear-ordering-problem.hxx"
#include "string"
#include "iostream"
#include "fstream"

template<class T>
inline T castString(std::string const & value);

template<> inline float castString<float>(std::string const & value){
    return std::stof(value);
}

template<> inline double castString<double>(std::string const & value){
    return std::stod(value);
}

template<> inline int castString<int>(std::string const & value){
    return std::stoi(value);
}

template<> inline long castString<long>(std::string const & value){
    return std::stol(value);
}


template<class T = long>
inline
LinearOrderingProblem<T> readDataset(std::ifstream & in){
    std::string line;

    // read first line
    getline(in, line);

    size_t numberOfElements = std::stoi(line);
    LinearOrderingProblem<T> problem(numberOfElements);

    size_t i = 0;

    while (getline(in, line)){
        size_t j = 0;
        std::string value;
        std::stringstream lineStream(line);

        while (lineStream >> value){
            if (i != j) {
                problem.cost(i, j) = castString<T>(value);
            }

            ++j;
        }
        ++i;
    }

    return problem;
}

template<class T = long>
inline
LinearOrderingProblem<T> readDataset(std::string const & fileName){
    std::ifstream fileStream(fileName);
    if (fileStream.is_open()){
        return readDataset<T>(fileStream);
    }
    fileStream.close();
}


#endif //PARTIAL_OPTIMALITY_LINEAR_ORDERING_LOLIP_READER_HXX

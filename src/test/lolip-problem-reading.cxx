//
// Created by dstein on 01.09.23.
//
//
// Created by dstein on 01.09.23.
//
#include <stdexcept>
#include "dscc/lolip2010/lolip-reader.hxx"

inline void test(const bool& pred) {
    if(!pred) throw std::runtime_error("Test failed.");
}

template<class T = long>
void testLolipProblemReading(){
    std::string fileName = "/home/dstein/Downloads/lolib_2010/IO/N-be75eec";

    readDataset<long>(fileName);
}


int main() {
    testLolipProblemReading();
    return 0;
}
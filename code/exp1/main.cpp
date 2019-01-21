#include <iostream>
#include "OpenMPExamples.h"
#include <omp.h>
#include<stdio.h>

int main() {
    omp_set_num_threads(4);
    std::cout << "Hello, World!" << std::endl;
    OpenMPExamples openMPExamples;
    openMPExamples.sample5();
    openMPExamples.sample6();
    return 0;
}
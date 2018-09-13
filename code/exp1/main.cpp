#include <iostream>
#include "OpenMPExamples.h"
#include <omp.h>

int main() {
    omp_set_num_threads(4);
    std::cout << "Hello, World!" << std::endl;
    OpenMPExamples openMPExamples;
    //openMPExamples.sample1();
    openMPExamples.sample4();
    return 0;
}
#include <iostream>
#include "OpenMPExamples.h"
#include <omp.h>
#include<stdio.h>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char** argv) {
    int num_threads;
    if(argc>1) {
        num_threads = stoi(argv[1]);
    }

    omp_set_num_threads(num_threads);
    cout << "Number of threads : " << num_threads << endl;
    OpenMPExamples openMPExamples;
    //openMPExamples.sample7();
    //openMPExamples.sample8();
    openMPExamples.sample11(num_threads);
    return 0;
}
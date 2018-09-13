//
// Created by vibhatha on 9/12/18.
//
#include "OpenMPExamples.h"
#include<iostream>
#include<stdio.h>
#include<omp.h>

using namespace std;

int calculate(int n, int i){
    int x = 0;
    x = n*i;
    for (int j = 0; j < x; ++j) {
        for (int k = 0; k < 10; ++k) {
            x = k *j + (k-j)*n + i + (i-n*j)*100;
        }
    }

    return x;
}

OpenMPExamples::OpenMPExamples() {

}

void OpenMPExamples::sample1() {
    int i, j, k = 0;
#pragma omp parallel for
    for (int i = 0; i < 5; ++i) {
        i++;
        std::cout << "Threads " << omp_get_num_threads() << "\n";
    }

}

void OpenMPExamples::sample2() {
    omp_set_num_threads(16); // OPTIONAL — Can also use
    // OMP_NUM_THREADS environment variable
    #pragma omp parallel
    {
        printf("hello, world!\n"); // Execute in parallel
    } // Implicit barrier/join

}

void OpenMPExamples::sample3() {
    clock_t start_time = clock();
    int i=0;
    int n=800000000;
    int *a = new int[n];
    #pragma omp parallel for
    for(i=0; i <n; i++){
        a[i] = i + n*3;
        //printf("%d ",a[i]);
    }
    clock_t end_time = clock();
    cout << "\nTime Taken for Parallel Loops : " << (end_time - start_time)/CLOCKS_PER_SEC << " s" << endl;

    clock_t start_time2 = clock();

    for(i=0; i <n; i++){
        a[i] = i + n*3;
        //printf("%d ",a[i]);
    }
    clock_t end_time2 = clock();
    cout << "\nTime Taken for Sequential Loops : " << (end_time2 - start_time2)/CLOCKS_PER_SEC << " s" << endl;
}

void OpenMPExamples::sample4() {
    clock_t start_time = clock();
    int i=0;
    int j =0;
    int n=50000;
    int *a = new int[n];
#pragma omp parallel for private(i,j) shared(a)
    for(i=0; i <n; i++){
        for (int j = 0; j < n; ++j) {
            a[i] = calculate(n,i);
        }

        //printf("%d ",a[i]);
    }
    clock_t end_time = clock();
    cout << "\nTime Taken for Parallel Loops : " << (end_time - start_time)/CLOCKS_PER_SEC << " s" << endl;

    clock_t start_time2 = clock();

    for(i=0; i <n; i++){
        for (int j = 0; j < n; ++j) {
            a[i] = calculate(n,i);
        }
        //printf("%d ",a[i]);
    }
    clock_t end_time2 = clock();
    cout << "\nTime Taken for Sequential Loops : " << (end_time2 - start_time2)/CLOCKS_PER_SEC << " s" << endl;
}

void OpenMPExamples::sample5() {

}

void OpenMPExamples::sample6() {

}

void OpenMPExamples::sample7() {

}

void OpenMPExamples::sample8() {

}

void OpenMPExamples::sample9() {

}

void OpenMPExamples::sample10() {

}



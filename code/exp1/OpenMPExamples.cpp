//
// Created by vibhatha on 9/12/18.
//
#include "OpenMPExamples.h"
#include<iostream>
#include<stdio.h>
#include<omp.h>
#include <math.h>       /* sin */

#define PI 3.14159265

using namespace std;

int calculate(int n, int i) {
    int x = 0;
    x = n * i;
    for (int j = 0; j < x; ++j) {
        for (int k = 0; k < 10; ++k) {
            x = k * j + (k - j) * n + i + (i - n * j) * 100;
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
    omp_set_num_threads(32); // OPTIONAL — Can also use
    // OMP_NUM_THREADS environment variable
#pragma omp parallel
    {
        printf("hello, world!\n"); // Execute in parallel
    } // Implicit barrier/join

}

void OpenMPExamples::sample3() {
    clock_t start_time = clock();
    double start_time_w = omp_get_wtime();
    int i = 0;
    int n = 800000;
    int *a = new int[n];
    int *b = new int[n];
    omp_set_num_threads(1);
    int num_threads = omp_get_num_threads();
    int sum = 0;
    #pragma omp parallel for
    for (i = 0; i < n; i++) {
        a[i] = i + n * 3;
        sum+= a[i];
        //printf("%d ",a[i]);
    }
    double end_time_w = omp_get_wtime();
    clock_t end_time = clock();
    cout << "Answer : " << sum << endl;
    cout << "\nTime Taken for Parallel Loops : " << (end_time - start_time) / CLOCKS_PER_SEC << " s, threads = " << num_threads << ", Ans = " << endl;
    cout << "\nWall Time Taken for Parallel Loops : " << (end_time_w - start_time_w) / CLOCKS_PER_SEC << " s, threads = " << num_threads << endl;

    clock_t start_time2 = clock();
    double start_time_w_2 = omp_get_wtime();
    int sum1 = 0;
    int j=0;
    for (j = 0; j < n; j++) {
        b[j] = j + n * 3;
        //a[i] = sin (a[i]*PI/180);
        sum1+= b[i];
        //printf("%d ",a[i]);
    }
    double end_time_w_2 = omp_get_wtime();
    clock_t end_time2 = clock();
    cout << "Answer : " << sum1 << endl;
    cout << "\nTime Taken for Sequential Loops : " << (end_time2 - start_time2) / CLOCKS_PER_SEC << " s" << endl;
    cout << "\nWall Time Taken for Sequential Loops : " << (end_time_w_2 - start_time_w_2) / CLOCKS_PER_SEC << " s" << endl;
    delete a;
    delete b;
}

void OpenMPExamples::sample4() {
    omp_set_dynamic(0);
    omp_set_num_threads(8);
    #pragma omp parallel
    printf("Num of Threads Inside Parallel Region : %d \n", omp_get_num_threads());


    clock_t start_time = clock();
    int i = 0;
    int j = 0;
    int n = 100000;
    cout << "Num of Threads Outside Parallel Region : " << omp_get_num_threads() << endl;
    int *a = new int[n];


#pragma omp parallel for private(i, j) shared(a)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; ++j) {
            a[i] = calculate(n, i);
        }

        //printf("%d ",a[i]);
    }
    clock_t end_time = clock();
    cout << "\nTime Taken for Parallel Loops : " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;

    clock_t start_time2 = clock();

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; ++j) {
            a[i] = calculate(n, i);
        }
        //printf("%d ",a[i]);
    }
    clock_t end_time2 = clock();
    cout << "\nTime Taken for Sequential Loops : " << (end_time2 - start_time2) / CLOCKS_PER_SEC << " s" << endl;
}
// TODO : incomplete examples
void OpenMPExamples::sample5() {
    omp_set_dynamic(0);
    omp_set_num_threads(8);
    int N = 100000;
    int *a = new int[N];
    int *b = new int[N];
    int *c = new int[N];
    int *d = new int[N];
    int i = 0;
    int sum = 0;
    for (int j = 0; j < N; ++j) {
        a[j] = j * 10 + 1;
        b[j] = j * 20 + 2;
        c[j] = j * 30 + 3;
        d[j] = j * 40 + 4;
    }
#pragma omp parallel default(none) shared(N, a, b, c, d, sum) private(i)
    {
#pragma  omp fr nowait
        for (i = 0; i < N; ++i) {
            a[i] += b[i];
        }
#pragma omp for nowait
        for (i = 0; i < N; i++) {
            c[i] += d[i];
        }
#pragma omp barrier

#pragma omp for nowait reduction(+:sum)
        for (i = 0; i < N; i++) {
            sum += a[i] + c[i];
        }
    }
    printf("Parallel Sum : %d \n", sum);

}

//TODO : incomplete examples 
void OpenMPExamples::sample6() {
    int N = 100000;
    int *a = new int[N];
    int *b = new int[N];
    int *c = new int[N];
    int *d = new int[N];
    int i = 0;
    int sum = 0;
    for (int j = 0; j < N; ++j) {
        a[j] = j * 10 + 1;
        b[j] = j * 20 + 2;
        c[j] = j * 30 + 3;
        d[j] = j * 40 + 4;
    }


    for (i = 0; i < N; ++i) {
        a[i] += b[i];
    }

    for (i = 0; i < N; i++) {
        c[i] += d[i];
    }


    for (i = 0; i < N; i++) {
        sum += a[i] + c[i];
    }

    printf("Sequential Sum : %d \n", sum);

}

void OpenMPExamples::sample7() {
    static long num_steps = 500000;
    double step;
    int i; double x,pi,sum = 0;
    step = 1.0 / (double) num_steps;
    double t1 = omp_get_wtime();
    for (i = 0; i < num_steps ; i++) {
        x = (i+0.5) * step;
        sum = sum + 4.0/ (1.0 + x*x);
    }
    pi = step * sum;
    double t2 = omp_get_wtime();
    cout << "Sequential : PI Est : " << pi << ", Time Taken : " << (t2-t1) << endl;
}

void OpenMPExamples::sample8() {
    omp_set_num_threads(16);
    static long num_steps = 500000;
    double step;
    int i; double x,pi,sum = 0;
    step = 1.0 / (double) num_steps;
    double t1 = omp_get_wtime();
    #pragma omp parallel for shared(step, num_steps) private(i,x) reduction ( + : sum )
    for (i = 0; i < num_steps ; i++) {
        x = (i+0.5) * step;
        sum = sum + 4.0/ (1.0 + x*x);
    }
    pi = step * sum;
    double t2 = omp_get_wtime();
    cout << "Parallel : PI Est : " << pi << ", Time Taken : " << (t2-t1) << endl;
}

void OpenMPExamples::sample9() {
    int itr = 100;
    int n = 100000;
    const int d = 3;
    double** x;
    double alpha = 0.001;
    double* w = new double [d];
    double* y = new double[n];

    //initialize
    x = new double *[n];
    for (int i = 0; i < n; ++i) {
        x[i] = new double[d];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            x[i][j] = i * j;
        }
    }

    for (int k = 0; k < d; ++k) {
        w[k] = 0;
    }

    for (int l = 0; l < n; ++l) {
        if(l%2==0) {
            y[l] = 1;
        } else {
            y[l] = -1;
        }
    }

    double *x1 = new double [d];
    double *x1y1 = new double [d];
    double y1 = 0;
    double start_time = omp_get_wtime();
    for (int m = 0; m < itr; ++m) {
        for (int i = 0; i < n; ++i) {
            x1 = x[i];

            y1 = y[i];

            for (int j = 0; j < d; ++j) {
                x1y1[j] = x1[j] * y1;
            }

            for (int k = 0; k < d; ++k) {
                w[k] = w[k] - alpha * x1y1[k];
            }
           // printVector(w, d);
        }
       // cout << "---------------------------" << endl;
    }
    double end_time = omp_get_wtime();
    cout << " Weight : ";
    for (int i1 = 0; i1 < d; ++i1) {
        cout << w[i1] << " ";
    }
    cout << endl;
    cout << "Sequential : Time Taken : " << (end_time - start_time) << " s" << endl;
}

void OpenMPExamples::sample10() {
    omp_set_num_threads(2);
    int itr = 100;
    int n = 100000;
    const int d = 3;
    double** x;
    double alpha = 0.001;
    double* w = new double [d];
    double* y = new double[n];

    //initialize
    x = new double *[n];
    for (int i = 0; i < n; ++i) {
        x[i] = new double[d];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            x[i][j] = i * j;
        }
    }

    for (int k = 0; k < d; ++k) {
        w[k] = 0;
    }

    for (int l = 0; l < n; ++l) {
        if(l%2==0) {
            y[l] = 1;
        } else {
            y[l] = -1;
        }
    }

    double *x1 = new double [d];
    double *x1y1 = new double [d];
    double y1 = 0;
    double start_time = omp_get_wtime();

    for (int m = 0; m < itr; ++m) {
        for (int i = 0; i < n; ++i) {
            x1 = x[i];

            y1 = y[i];

            dotProduct(x1, y1, x1y1, d);
            /*for (int j = 0; j < d; ++j) {
                x1y1[j] = x1[j] * y1;
            }*/

            for (int k = 0; k < d; ++k) {
                w[k] = w[k] - alpha * x1y1[k];
            }
            // printVector(w, d);
        }
        // cout << "---------------------------" << endl;
    }
    double end_time = omp_get_wtime();
    cout << " Weight : ";
    for (int i1 = 0; i1 < d; ++i1) {
        cout << w[i1] << " ";
    }
    cout << endl;
    cout << "Parallel : Time Taken : " << (end_time - start_time) << " s" << endl;
}

void OpenMPExamples::printVector(double *v, int size) {
    for (int i = 0; i < size; ++i) {
        cout << v[i] << " ";
    }
    cout << endl;
}

void OpenMPExamples::sample11(int num_threads) {
    omp_set_num_threads(num_threads);
    int d = 254;
    int par = 2;
    int samples = 270000;
    int itz = samples/par;
    double *a = new double [d];
    double b = 10;
    double start_time = omp_get_wtime();
    int id = 0;
    int i,j,k;
#pragma omp parallel private(i,j,k) shared(itz)
 {
     #pragma omp barrier
     for (j = 0; j < 1000 ; ++j) {
        #pragma omp barrier
         for (k = 0; k < itz; ++k) {
            #pragma omp for
             for (i = 0; i < d; ++i) {
                 a[i] = i*b;
             }
        #pragma omp barrier
         }
     #pragma omp barrier
     }

 }
    double end_time = omp_get_wtime();

    cout << "Time Taken : " << (end_time - start_time) << " s " << endl;
}


void OpenMPExamples::dotProduct(double *a, double b, double *ans, int size) {
    #pragma omp parallel for private(ans, a, b) shared(size)
    for (int j = 0; j < size; ++j) {
        ans[j] = a[j] * b;
    }
}



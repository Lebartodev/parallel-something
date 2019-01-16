//#include <iostream>
//#include "math.h"
//#include "stdlib.h"
//#include "stdio.h"
//#include "time.h"
//#include <ctime>
//#include "mpi.h"
//
//static double T = 200;
//static int Nx = 256;
//static int Nt = 10000;
//static double hx = 1.0 / Nx;
//static double ht = 1.0 / Nt;
//static double **v;
//static double **vnext;//p+1 слой
//static double a = 1;
//static int procNum;
//static int rank;
//
//double getV(int i, int j);
//
//double **allocate_solution2D() {
//    double **res = (double **) malloc((Nt + 1) * sizeof(double *));
//    for (int i = 0; i < Nt; i++)
//        res[i] = (double *) malloc((Nx + 1) * sizeof(double));
//    return res;
//}
//
//double rfunc(int i, int j) {
//    double x = 1.0 * (i + 1) / Nx;
//    double y = 1.0 * (j + 1) / Nx;
//    return 1.0 * x * y * (1.1 - x) * (1.1 - y);
//}
//
//double calculateVnext(int i, int j) {
//    return v[i][j] + a * ht * ((getV(i - 1, j) - 2 * getV(i, j) + getV(i + 1, j)) / (hx * hx) +
//                               (getV(i, j - 1) - 2 * getV(i, j) + getV(i, j + 1)) / (hx * hx));
//}
//
//double getV(int i, int j) {
//    if (i < 0)
//        i = Nx - 1;
//    if (i > Nx - 1)
//        i = 0;
//
//    if (j < 0)
//        j = Nx - 1;
//    if (j > Nx - 1)
//        i = 0;
//
//    return v[i][j];
//
//}
//
//
//void difference_scheme() {
//    int portion = 200;
//    int downLimit = (Nx / procNum) * rank;
//    int upLimit = (Nx / procNum) * (rank + 1);
//    double *recvLeft = (double *) malloc(Nx * sizeof(double));
//    double *recvRight = (double *) malloc(Nx * sizeof(double));
//    MPI_Status stat;
//    for (int p = 0; p < Nt - 1; p++) {
//
//        double *sendLeft = v[downLimit];
//        double *sendRight = v[upLimit];
//        if (rank == 0) {
//            MPI_Send(sendLeft, Nx, MPI_DOUBLE, procNum - 1, 1, MPI_COMM_WORLD);
//        } else {
//            MPI_Send(sendLeft, Nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
//        }
//
//        if (rank == procNum - 1) {
//            MPI_Send(sendRight, Nx, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
//        } else {
//            MPI_Send(sendRight, Nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
//        }
//
//        if (rank == 0) {
//            MPI_Recv(recvLeft, Nx, MPI_DOUBLE, procNum - 1, 1, MPI_COMM_WORLD, &stat);
//        } else {
//            MPI_Recv(recvLeft, Nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &stat);
//        }
//
//        if (rank == procNum - 1) {
//            MPI_Recv(recvRight, Nx, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
//        } else {
//            MPI_Recv(recvRight, Nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &stat);
//        }
//
//        if (downLimit != 0) {
//            v[downLimit] = recvLeft;
//        } else {
//            v[Nx - 1] = recvLeft;
//        }
//        if (upLimit != Nx) {
//            v[upLimit] = recvRight;
//        } else {
//            v[0] = recvLeft;
//        }
//
//
//        for (int i = downLimit; i < upLimit; i++) {
//
//            for (int j = 0; j < Nx; j++) {
//                vnext[i][j] = calculateVnext(i, j);
//            }
//        }
//
//        for (int i = downLimit; i < upLimit; i++) {
//            for (int j = 0; j < Nx; j++) {
//                v[i][j] = vnext[i][j];
//            }
//        }
//    }
//
//
//}
//
//int main(int argc, char **argv) {
//    v = allocate_solution2D();
//    vnext = allocate_solution2D();
//
//    for (int m = 0; m < Nx; m++) {
//        for (int n = 0; n < Nx; n++) {
//            v[m][n] = rfunc(m, n);
//        }
//    }
//
////    if (ht > (1.0 / 2.0) * hx * hx) {
////        std::cout << "Egor ";
////        std::cout << " (1.0 / 2.0) * hx * hx =  " << (1.0 / 2.0) * hx * hx;
////        std::cout << "ht = " << ht;
////        throw MATH_ERREXCEPT;
////    }
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    double time = MPI_Wtime();
//
//
//    difference_scheme();
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    time -= MPI_Wtime();
//    if (!rank) {
//        printf("%d\t\t%.8f\n", Nx, -time);
//    }
//
//    MPI_Finalize();
//    return 0;
//}
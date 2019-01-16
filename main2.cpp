//#include <iostream>
//#include "math.h"
//#include "stdlib.h"
//#include "stdio.h"
//#include "time.h"
//
//static double T = 200;
//static int Nx = 100;
//static int Nt = 2000;
//
//
//double **allocate_solution2D() {
//    double **res = (double **) malloc((Nt + 1) * sizeof(double *));
//    for (int i = 0; i < Nt; i++)
//        res[i] = (double *) malloc((Nx + 1) * sizeof(double));
//    return res;
//}
//
//double rfunc(int i, int j, int Nx, int Ny) {
//    double x = 1.0 * (i + 1) / Nx;
//    double y = 1.0 * (j + 1) / Ny;
//    return 1.0 * x * y * (1.1 - x) * (1.1 - y);
//}
//
//void difference_scheme() {
//    double hx = 100.0 / Nx;
//    double ht = 200.0 / Nt;
//
//    if (ht > (1.0 / 2.0) * hx * hx) {
//        std::cout << "Pizdec";
//        throw MATH_ERREXCEPT;
//    }
//
//    double **v = allocate_solution2D();
//    double **vnext = allocate_solution2D();//p+1
//    double **w = allocate_solution2D();
//    double **wnext = allocate_solution2D();//p+1
//
//    for (int m = 0; m < Nx; m++) {
//        for (int n = 0; n < Nx; n++) {
//            v[m][n] = rfunc(m, n, Nx, Nx);
//        }
//    }
//    for (int i = 0; i < Nx; i++) {
//        for (int j = 1; j < Nx; j++) {
//            std::cout << v[i][j] << " ";
//
//        }
//        std::cout << std::endl;
//
//    }
//    std::cout << std::endl;
//    std::cout << std::endl;
//
//
//    for (int p = 0; p < Nt - 1; p++) {
//        for (int n = 0; n < Nx; n++) {
//            vnext[0][n] = v[0][n] + (v[1][n] - 2 * v[0][n] + v[Nx - 1][n]) * (ht / hx / hx);
//            vnext[Nx - 1][n] = v[Nx - 1][n] + (v[0][n] - 2 * v[Nx - 1][n] + v[Nx - 2][n]) * (ht / hx / hx);
//        }
//        for (int i = 1; i < Nx - 1; i++) {
//            for (int j = 0; j < Nx; j++) {
//                vnext[i][j] = v[i][j] + (v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) * (ht / hx / hx);
//            }
//        }
//        for (int i = 0; i < Nx; i++) {
//            for (int j = 0; j < Nx; j++) {
//                w[i][j] = vnext[i][j];
//            }
//        }
//        for (int m = 0; m < Nx; m++) {
//            wnext[m][0] = w[m][0] + (w[m][1] - 2 * w[m][0] + w[m][Nx - 1]) * (ht / hx / hx);
//            wnext[m][Nx - 1] = w[m][Nx - 1] + (w[m][0] - 2 * w[m][Nx - 1] + w[m][Nx - 2]) * (ht / hx / hx);
//        }
//        for (int i = 0; i < Nx; i++) {
//            for (int j = 1; j < Nx - 1; j++) {
//                wnext[i][j] = w[i][j] + (w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) * (ht / hx / hx);
//            }
//        }
//
//
//        for (int i = 0; i < Nx; i++) {
//            for (int j = 0; j < Nx; j++) {
//                v[i][j] = vnext[i][j];
//            }
//        }
//
//        for (int i = 0; i < Nx; i++) {
//            for (int j = 0; j < Nx; j++) {
//                std::cout << wnext[i][j] << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//        std::cout << std::endl;
//        std::cout << std::endl;
//    }
//    for (int i = 0; i < Nx; i++) {
//        for (int j = 0; j < Nx; j++) {
//            std::cout << wnext[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//
//
//}
//
//int main(int argc, char **argv) {
//    difference_scheme();
//
//    return 0;
//}
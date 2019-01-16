#include <iostream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "time.h"
#include <ctime>

static double T = 200;
static int Nx = 50;
static int Nt = 100;
static double hx = 100.0 / Nx;
static double ht = 200.0 / Nt;
static double **v;
static double **vnext;//p+1 слой
static double a = 1;

double getV(int i, int j);

double **allocate_solution2D() {
    double **res = (double **) malloc((Nt + 1) * sizeof(double *));
    for (int i = 0; i < Nt; i++)
        res[i] = (double *) malloc((Nx + 1) * sizeof(double));
    return res;
}

double rfunc(int i, int j) {
    double x = 1.0 * (i + 1) / Nx;
    double y = 1.0 * (j + 1) / Nx;
    return 1.0 * x * y * (1.1 - x) * (1.1 - y);
}

double calculateVnext(int i, int j) {
    return v[i][j] + a * ht * ((getV(i - 1, j) - 2 * getV(i, j) + getV(i + 1, j)) / (hx * hx) +
                               (getV(i, j - 1) - 2 * getV(i, j) + getV(i, j + 1)) / (hx * hx));
}

double getV(int i, int j) {
    if (i < 0)
        i = Nx - 1;
    if (i > Nx - 1)
        i = 0;

    if (j < 0)
        j = Nx - 1;
    if (j > Nx - 1)
        i = 0;

    return v[i][j];

}


void difference_scheme() {

    if (ht > (1.0 / 2.0) * hx * hx) {
        std::cout << "Egor ";
        std::cout << " (1.0 / 2.0) * hx * hx =  "<< (1.0 / 2.0) * hx * hx;
        std::cout << "ht = "<<ht;
        throw MATH_ERREXCEPT;
    }


    for (int m = 0; m < Nx; m++) {
        for (int n = 0; n < Nx; n++) {
            v[m][n] = rfunc(m, n);
        }
    }
    for (int i = 0; i < Nx; i++) {
        for (int j = 1; j < Nx; j++) {
            std::cout << v[i][j] << " ";

        }
        std::cout << std::endl;

    }
    std::cout << std::endl;
    std::cout << std::endl;

    for (int p = 0; p < Nt - 1; p++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Nx; j++) {
                vnext[i][j] = calculateVnext(i, j);
            }
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Nx; j++) {
                v[i][j] = vnext[i][j];
            }
        }
    }

    for (int i = 0; i < Nx; i++) {
        for (int j = 1; j < Nx; j++) {
            std::cout << v[i][j] << " ";

        }
        std::cout << std::endl;

    }
    std::cout << std::endl;
    std::cout << std::endl;


}

int main(int argc, char **argv) {
    unsigned int start_time = clock(); // начальное время


    v = allocate_solution2D();
    vnext = allocate_solution2D();
    difference_scheme();
    unsigned int end_time = clock(); // конечное время
    unsigned int search_time = end_time - start_time; // искомое время
    std::cout << search_time << std::endl;

    return 0;
}
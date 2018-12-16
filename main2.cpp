#include <iostream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "time.h"

static double T = 200;
static int Nx = 20;
static int Nt = 20;


double **allocate_solution2D() {
    double **res = (double **) malloc((Nt + 1) * sizeof(double *));
    for (int i = 0; i < Nt; i++)
        res[i] = (double *) malloc((Nx + 1) * sizeof(double));
    return res;
}

double rfunc(int i, int Nx) {
    double x = 1.0 * i / Nx;
    return -x * (x - 1.0);
}

void difference_scheme() {
    double hx = 100.0 / Nx;
    double ht = 20.0 / Nt;

    if (ht > (1.0 / 2.0) * hx * hx) {
        std::cout << "Pizdec";
        throw MATH_ERREXCEPT;
    }

    double **v = allocate_solution2D();
    double **w = allocate_solution2D();
    for (int i = 0; i < Nx; i++)
        v[0][i] = rfunc(i, Nx);


    for (int p = 0; p < Nt - 1; p++) {
        v[p + 1][0] = v[p][0] + (v[p][1] - 2 * v[p][0] + v[p][Nx - 1]) * (ht / hx / hx);
        v[p + 1][Nx - 1] = v[p][Nx - 1] + (v[p][0] - 2 * v[p][Nx - 1] + v[p][Nx - 2]) * (ht / hx / hx);

        for (int j = 1; j < Nx - 1; j++) {
            v[p + 1][j] = v[p][j] + (v[p][j + 1] - 2 * v[p][j] + v[p][j - 1]) * (ht / hx / hx);
        }


        for (int j = 0; j < Nx; j++) {
            w[p][j] = v[p + 1][j];
        }
        w[p + 1][0] = w[p][0] + (w[p][1] - 2 * w[p][0] + w[p][Nx - 1]) * (ht / hx / hx);
        w[p + 1][Nx - 1] = w[p][Nx - 1] + (w[p][0] - 2 * w[p][Nx - 1] + w[p][Nx - 2]) * (ht / hx / hx);
        for (int j = 1; j < Nx - 1; j++) {
            w[p + 1][j] = w[p][j] + (w[p][j + 1] - 2 * w[p][j] + w[p][j - 1]) * (ht / hx / hx);
        }


        for (int j = 0; j < Nx; j++) {
            std::cout << w[p + 1][j] << " ";
        }
        std::cout << std::endl;


    }

}

int main(int argc, char **argv) {
    difference_scheme();

    return 0;
}
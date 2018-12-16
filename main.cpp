//#include "math.h"
//#include "stdlib.h"
//#include "stdio.h"
//#include "mpi.h"
//#include "time.h"
//
//static double* bufA = 0;
//static double* bufB = 0;
//
//double max(double x, double y){
//    return x>y?x:y;
//}
//
//int max(int x, int y){
//    return x>y?x:y;
//}
//
//void init_buffers(int N){
//    bufA = (double*)malloc(N*sizeof(double));
//    bufB = (double*)malloc(N*sizeof(double));
//}
//
//double** allocate_solution2D(int Nx, int Ny){
//    double** res = (double**)malloc((Nx+1)*sizeof(double*));
//    for(int i=0; i<=Nx; i++)
//        res[i] = (double*)malloc((Ny+1)*sizeof(double));
//    return res;
//}
//
//
//double rfunc(int i, int j, int Nx, int Ny){
//    double x = 1.0*i/Nx;
//    double y = 1.0*j/Ny;
//    return 100*x*y*(x-1.0)*(y-1.0);
//}
//
//void init_condition(double** u, int Nx, int Ny){
//    for(int i=0; i<=Nx; i++)
//        for(int j=0; j<=Ny; j++)
//            u[i][j]= rfunc(i,j,Nx,Ny);
//}
//
//
//
//void difference_scheme (double** u, int Nx, int Ny, int Nt=1){
//    double hx = 1.0/Nx;
//    double hy = 1.0/Ny;
//    double ht = 1.0/Nt;
//    double gx = ht/hx/hx;
//    double gy = ht/hy/hy;
//    double fi;
//    for(int k=0; k<Nt; k++){
////************** x-direct  begin ******************************************************************************
//        for(int j=0; j<=Ny; j++){
//            bufA[0] = 0.0;
//            bufB[0] = u[0][j];
//            for(int i=1; i<Nx; i++)
//            {
//                fi = (1+2*gx)-gx*bufA[i-1];
//                bufA[i] = gx/fi;
//                bufB[i]= (u[i][j]+gx*bufB[i-1])/fi;
//            }
//            // back substitution
//            for(int i=Nx;i>0;i--){
//                u[i-1][j]=bufA[i-1]*u[i][j]+bufB[i-1];
//            }
//        }
////************** x-direct  end ********************************************************************************
////************** y-direct  begin ******************************************************************************
//        for(int i=0; i<=Nx; i++){
//            bufA[0] = 0.0;
//            bufB[0] = u[i][0];
//            for(int j=1; j<Ny; j++)
//            {
//                fi = (1+2*gy)-gy*bufA[j-1];
//                bufA[j] = gy/fi;
//                bufB[j]= (u[i][j]+gy*bufB[j-1])/fi;
//            }
//            for(int j=Ny;j>0;j--){
//                u[i][j-1]=bufA[j-1]*u[i][j]+bufB[j-1];
//            }
//        }
////************** y-direct  end ******************************************************************************
//    }
//}
//
//
//
//int main(int argc, char** argv){
//    //int Nx = (1<<6)-1;
//    //int Ny = (1<<6)-1;
//    //printf("%d\n", (-1)%50);
//    int Nx = 1023;
//    int Ny = Nx;
//    int Nt = 1023;
//    int proc, rank;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &proc);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    if(rank){
//        return 0;
//    }
//    double** u = allocate_solution2D(Nx, Ny);
//    init_condition(u, Nx, Ny);
//    init_buffers(max(Nx, Ny));
//    double time = MPI_Wtime();
//    difference_scheme (u, Nx, Ny, Nt);
//    time -= MPI_Wtime();
//    if(!rank){
//        printf("%d\t\t%.8f\n", Nx, -time);
//    }
//    MPI_Finalize();
//    return 0;
//}
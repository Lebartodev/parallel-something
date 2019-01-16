#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "mpi.h"
#include "time.h"

static double** pbufA = 0;
static double** pbufB = 0;

double max(double x, double y){
    return x>y?x:y;
}


void init_buffers(int N, int rank, int proc){
    pbufA = (double**)malloc((2*(proc-rank)+1)*sizeof(double*));
    pbufB = (double**)malloc((2*(proc-rank)+1)*sizeof(double*));
    for(int i=0; i<2*(proc-rank)+1; i++){
        pbufA[i] = (double*)malloc(N*sizeof(double));
        pbufB[i] = (double*)malloc(N*sizeof(double));
    }
}

double** allocate_solution(int Nx, int Ny, int proc){
    int lx = (Nx+1)/proc;
    double** res = (double**)malloc(lx*sizeof(double*));
    for(int i=0; i<lx; i++)
        res[i] = (double*)malloc((Ny+1)*sizeof(double));
    return res;
}

double rfunc(int i, int j, int Nx, int Ny, int rank, int proc){
    i += (rank-1)*(Nx+1)/proc;
    double x = 1.0*i/Nx;
    double y = 1.0*j/Ny;
    return 100*x*y*(1.0-x)*(1.0-y);
}

void init_condition(double** u, int Nx, int Ny, int rank, int proc){
    int lx = (Nx+1)/proc;
    for(int i=0; i<lx; i++)
        for(int j=0; j<=Ny; j++){
            //u[i][j]= 1.0;
            u[i][j]= rfunc(i,j,Nx, Ny, rank, proc);
        }
}
void difference_scheme (double** u, int Nx, int Ny, int rank, int proc, int Nt=1){
    int lx = (Nx+1)/proc;
    double hx = 1.0/Nx;
    double hy = 1.0/Ny;
    double ht = 1.0/Nt;
    double gx = ht/hx/hx;
    double gy = ht/hy/hy;
    double fi;
    double* bufA;
    double* bufB;
    int j;
    int buf_n = 2*(proc-rank)+1;
    double* buf = (double*)malloc(2*sizeof(double));
    MPI_Status stat;

    for(int k=0; k<Nt; k++){
//************** x-direct  begin ******************************************************************************
        //for(int j=1; j<Ny; j++){
        for(int tact=1-proc; tact<=Ny+proc-1; tact++){
            if((tact+proc-rank>=0)&&(tact+proc-rank<=Ny)){
                j=tact+proc-rank;
                bufA=pbufA[j%buf_n];
                bufB=pbufB[j%buf_n];
                if(rank==1){
                    bufA[0] = 0.0;
                    bufB[0] = u[0][j];
                } else {
                    MPI_Recv(buf, 2, MPI_DOUBLE, rank-2, j, MPI_COMM_WORLD, &stat);
                    fi = (1+2*gx)-gx*buf[0];
                    bufA[0] = gx/fi;
                    bufB[0]= (u[0][j]+gx*buf[1])/fi;
                }
                for(int i=1; i<lx; i++)
                {
                    fi = (1+2*gx)-gx*bufA[i-1];
                    bufA[i] = gx/fi;
                    bufB[i]= (u[i][j]+gx*bufB[i-1])/fi;
                }
                if(rank!=proc){
                    buf[0] = bufA[lx-1];
                    buf[1] = bufB[lx-1];
                    MPI_Send(buf, 2, MPI_DOUBLE, rank, j, MPI_COMM_WORLD);
                }
            }
            if((tact-proc+rank>=0)&&(tact-proc+rank<=Ny)){
                // back substitution
                j=tact-proc+rank;
                bufA=pbufA[j%buf_n];
                bufB=pbufB[j%buf_n];
                if(rank==proc){
                } else {
                    MPI_Recv(buf, 1, MPI_DOUBLE, rank, j, MPI_COMM_WORLD, &stat);
                    u[lx-1][j] = bufA[lx-1]*buf[0]+bufB[lx-1];
                }

                for(int i=lx-1;i>0;i--){
                    u[i-1][j]=bufA[i-1]*u[i][j]+bufB[i-1];
                }
                if(rank!=1){
                    buf[0]=u[0][j];
                    MPI_Send(buf, 1, MPI_DOUBLE, rank-2, j, MPI_COMM_WORLD);
                }
            }
        }
//************** x-direct  end ********************************************************************************
//************** y-direct  begin ******************************************************************************
        bufA = pbufA[0];
        bufB = pbufB[0];
        for(int i=0; i<lx; i++){
            bufA[0] = 0.0;
            bufB[0] = u[i][0];
            for(int j=1; j<Ny; j++)
            {
                fi = (1+2*gy)-gy*bufA[j-1];
                bufA[j] = gy/fi;
                bufB[j]= (u[i][j]+gy*bufB[j-1])/fi;
            }
            // back substitution
            for(int j=Ny;j>0;j--){
                u[i][j-1]=bufA[j-1]*u[i][j]+bufB[j-1];
            }
        }
//************** y-direct  end ******************************************************************************
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(buf);
}

int main(int argc, char** argv){
    //int Nx = (1<<10)-1;
    //int Ny = (1<<10)-1;
    int Nx = 255;
    int Ny = Nx;
    int Nt = 2048;
    int proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if((Nx+1)%proc){
        printf(" Incompatible sizes for %d processors", proc);
        return 0;
    }
    //printf("hello %d \n", rank);
    double** u = allocate_solution(Nx, Ny, proc);
    init_condition(u, Nx, Ny, rank+1, proc);
    init_buffers(max(Nx, Ny)+1, rank+1, proc);
    double time = MPI_Wtime();
    difference_scheme(u, Nx, Ny, rank+1, proc, Nt);
    MPI_Barrier(MPI_COMM_WORLD);
    time -= MPI_Wtime();
    if(!rank){
        printf("%d\t\t%.8f\n", Nx, -time);
    }
    MPI_Finalize();
}

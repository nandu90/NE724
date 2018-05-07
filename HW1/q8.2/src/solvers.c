/***************************************************************************

Author: nsaini
Created: 2018-03-11

***************************************************************************/


#include "solvers.h"
#include "common.h"
#include "memory.h"

void solveSystem(double **vand, double *rhs, double *soln, int size)
{
    int i,j,k;
    //------------------------------------------------------------------------//
    //Set up variables for LAPACK
    double *A;
    allocator1(&A, size*size);

    k=0;
    for(i=0; i<size; i++)
    {
	for(j=0; j<size; j++)
	{
	    A[k] = vand[j][i];  //Be carefule over here. Store Column-wise for correct result
	    k++;
	}
    }

    int M = size; //Rows of matrix
    int N = size;        //Columns of matrix
    int NRHS = 1;
    int LDA = M;
    int *IPIV;
    iallocator1(&IPIV, M);

    double *B;
    allocator1(&B, M);

    for(i=0; i<size; i++)
    {
	B[i] = rhs[i];
    }

    int LDB = M;
    int INFO;

    /*int LWORK = -1;
    double WORK[1];
    char TRANS = 'N';*/
    
    //------------------------------------------------------------------------//
    //Solve the system to get the coefficients
   
    dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
    

    //Assign to the soln array
    for(i=0; i<size; i++)
    {
	soln[i] = B[i];
    }
    

    


    deallocator1(&A, size*size);
    ideallocator1(&IPIV, M);
    deallocator1(&B, M);
}

/***************************************************************************

Author: nsaini
Created: 2018-04-29

***************************************************************************/



#include "common.h"

#include "memory.h"
#include "solvers.h"



char* getexepath()
{
  static char cwd[1024];
  char *err = getcwd(cwd,sizeof(cwd));
  if(err == NULL)
    {
      printf("Error getting the current working directory\n.Exiting...");
      exit(1);
    }
  else
    {
      return cwd;
    }
}

char* concat(char s1[], char s2[])
{
    char* result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the null-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

void applyBC(double *phi, int tcells)
{
    //left
    double phiB = 2.0;
    
    phi[0] = 2.0*phiB - phi[1];

    //right
    phi[tcells-1] = phi[tcells-2];
}


void output(double*x, double *phi, int iter, int tcells)
{
    int i;
    
     //Print the output
    FILE *out1;
    char* path1;
    path1 = concat(getexepath(),"/output/out_00");
    char buffer0[12];
    snprintf(buffer0,12,"%d",iter);
    path1 = concat(path1,buffer0);
    path1 = concat(path1,".txt");
    out1 = fopen(path1,"w");
    free(path1);

    for(i=1; i<tcells-1; i++)
    {
	fprintf(out1,"%.4e %.4e\n",x[i], phi[i]);
    }
    
    fclose(out1);
}


int main(int argc, char **argv)
{
    int nprocs, myrank;
    int master;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    master = 0;

    if(myrank == master)
    {
	printf("----------------------------------------------------------------\n");
	printf("--------------------------And so it begins----------------------\n");
	printf("----------------------------------------------------------------\n");
    }
    
    double time1 = MPI_Wtime();

    //------------------------------------------------------------------------//
    //Create necessary output directories//
    //Let only the master create directory
    if(myrank == master)
    {
	char* path;
	path = concat(getexepath(), "/output");
	mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	free(path);
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    int i,j;
    //------------------------------------------------------------------------//

    int n = 20;

    int cells = n-1;

    double len = 1.0;

    double deltax = len/cells;

    int totaln = n+2;

    int tcells = cells + 2;

    double *phi;
    allocator1(&phi, tcells);
    

    double *phinew;
    allocator1(&phinew, tcells);

    double *phihalf;
    allocator1(&phihalf, tcells);

    for(i=0; i<tcells; i++)
    {
	phi[i] = 1.0;
    }

    double *phipast;
    allocator1(&phipast, tcells);
    for(i=0; i<tcells; i++)
    {
	phipast[i] = 1.0;
    }
    

    double *x;
    allocator1(&x, totaln);

    x[0] = 0.0-deltax/2.0;

    for(i=1; i<totaln; i++)
    {
	x[i] = x[i-1] + deltax;
    }

    double Cv = atof(argv[1]);

   
    
    printf("Selected Courant Number is %.4f\n",Cv);

    double u = 1.0;

    double deltat = Cv*deltax/u;

    printf("Corresponding time step size is %.4f\n",deltat);

    int iter = 0;


    double tol = 1e-8;

    double res, prevres;

    double error;

    double **A;
    allocator2(&A, cells, cells);

    double *vector;
    allocator1(&vector, cells);

    double *B;
    allocator1(&B, cells);

    //------------------------------------------------------------------------//
    //Matrix coefficients are constant for scheme 2. Therefore stored outside
    for(i=0; i<cells; i++)
    {
	for(j=0; j<cells; j++)
	{
	    if(i == j)
	    {
		A[i][j] = (3.0 + 2.0*Cv);
	    }
	    else if(j == i-1)
	    {
		A[i][j] = -2.0*Cv;		
	    }
	}
    }

    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //For printing output
    
    //------------------------------------------------------------------------//


    for(iter=0; iter<1000; iter++)
    {
	
  
	//Construct the RHS Vector
	for(i=0; i<cells; i++)
	{
	    if(i==0)
	    {
		B[i] = 4.0*phi[i+1] - phipast[i+1] + 2.0*Cv*phi[i];
	    }
	    else
	    {
		B[i] = 4.0*phi[i+1] - phipast[i+1];
	    }
	}
	
	solveSystem(A, B, vector, cells);
	
	//Assign solved values back
	for(i=0; i<cells; i++)
	{
	    phinew[i+1] = vector[i];  
	}
	


	applyBC(phinew, tcells);

	error = 0.0;
	    
	for(i=1; i<tcells-1; i++)
	{
	    error += pow(phinew[i] - phi[i],2.0);
	}

	res = sqrt(error);
	
		
	if(iter > 0)
	{
	    printf("Residual is %.4f\n",res);
	    if(fabs(res) < tol)
	    {
		printf("Solution has converged in %d iterations\n", iter+1);
		break;
	    }
	}

	for(i=0; i<tcells; i++)
	{
	    phipast[i] = phi[i];
	    phi[i] = phinew[i];
	}

	prevres = res;

	output(x, phi, iter, tcells);
    }
	

    //------------------------------------------------------------------------//
   
    //------------------------------------------------------------------------//

    
    
    double time2 = MPI_Wtime();
    double secs = time2-time1;
    if(myrank == master)
    {
	printf("Total run time: %.6f secs\n",secs);
    }
    
    if(myrank == master)
    {
      printf("----------------------------------------------------------------\n");
      printf("---------------------That's all folks!--------------------------\n");
      printf("----------------------------------------------------------------\n");
    }
    MPI_Finalize();
}
      

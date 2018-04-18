/***************************************************************************

Author: nsaini
Created: 2018-04-12

***************************************************************************/





#include "common.h"
#include "fileIO.h"
#include "correlations.h"
#include "generalFunc.h"
#include "memory.h"
#include "rhs.h"



int main(int argc, char **argv)
{
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
    //Read control file
    control();
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Read the Node Geometry data file
    struct nodeData *nData = (struct nodeData*)malloc(nodes*sizeof(struct nodeData));
    callocator2(&components, nodes); 
    if(myrank == master)printf("Read the node data file\n\n");
    geomData(nData);
    

    //------------------------------------------------------------------------//
    //Construct the Coefficient matrix for Mass-momentum equation
    
    double a1;
    double b1;
    double a2;
    double b2;
    double ac;
    double bc;

    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Time loop will come here
    double t=0.0;
    double deltat = iniDeltat;
    int iNewton;

    //------------------------------------------------------------------------//
    //Solution variables
    double m1dot = iniM1dot;
    double m2dot = (iniMcdot - iniM1dot)/(nloops - 1);
    double mcdot = iniMcdot;

    double m1old = m1dot;
    double m2old = m2dot;
    double mcold = mcdot;
    double deltaPcore = 0.0;

    double mu = iniMu;
    double rho = iniRho;
    double rhoold = rho;

    double RCP = 100.0;
    double volFlowRate=0.0;
    double volFlowRateOld = 0.0;
    double error;
    int iter = 0;
    //------------------------------------------------------------------------//

    printf("m1dot m2dot mcdot mu rho = %.4e %.4e %.4e %.4e %.4e\n",m1dot, m2dot, mcdot, mu, rho);
    
    while(fabs(volFlowRate - targetMdot) > 1e-5)//(t <= finalTime)//(fabs(mcdot) > 1e-8)
    {
	iter++;
	double mcCheck = mcold;
	double change;
	//Newton's iteration loop
	printf("----------------------------Time Step = %.4e----------------------------\n",t);
	for(iNewton=0; iNewton<maxNewton; iNewton++)
	{
	    loopTerms(&a1, &b1, deltat, nData, m1dot, m1old, rho, mu, rhoold);
	    b1 += RCP;
	    //printf("a1 and b1 outside = %.4f %.4f\n", a1, b1);
	    loopTerms(&a2, &b2, deltat, nData, m2dot, m2old, rho, mu, rhoold);
	    b2 += RCP;
	    //printf("a2 and b2 outside = %.4f %.4f\n", a2, b2);
	    coreTerms(&ac, &bc, deltat, nData, mcdot, mcold, rho, mu, rhoold);
	    //printf("ac and bc outside = %.4f %.4f\n", ac, bc);
	    
	    //Solution from Cramer's loop for next iteration
	    m1dot = (b1*a2 + (nloops-1)*b1*ac - (nloops-1)*ac*b2 + a2*bc)/(a1*a2 + (nloops-1)*a1*ac + a2*ac);
	    deltaPcore = b1 - a1*m1dot;
	    m2dot = (b2 - deltaPcore)/a2;
	    mcdot = m1dot + (nloops-1)*m2dot;

	    //printf("m dot values = %.4f %.4f %.4f %.4f\n",m1dot,m2dot,mcdot,deltaPcore);
	    
	    //exit(1);
	    change = (mcCheck - mcdot)/mcCheck;
	    printf("Newton Iteration = %d\n", iNewton+1);
	    printf("Core Mass Flow rate = %.4e\n", mcdot);
	    printf("Relative change in mass core flow rate = %.4e\n", change);
	    printf("Pressure drop across the core = %.4e\n",deltaPcore);
	    mcCheck = mcdot;

	    if(fabs(change) < 1e-8)
	    {
		printf("Newton's Iteration converged in %d iterations\n",iNewton+1);
		break;
	    }
	    
	}

	volFlowRate = mcdot/nData[0].Ax;
	printf("Volumetric flow Rate in the core = %.2e\n",volFlowRate);
	printf("Difference between target and required value  = %.2e\n",targetMdot-volFlowRate);
	if(fabs(volFlowRate - targetMdot) > 1e-5)
	{
	    /*if(volFlowRate < targetMdot)
	    {
		printf("Ramping Up RCP DeltaP\n");
		RCP += 0.1*fabs(volFlowRate - targetMdot);
	    }
	    else
	    {
		printf("Ramping Down RCP DeltaP\n");
		RCP -= 0.1*fabs(volFlowRate - targetMdot); 
		}*/
	    error = targetMdot - volFlowRate;
	    RCP += 0.1*error + 0.1*(volFlowRate - volFlowRateOld)*deltat;
	    
	}
	else
	{
	    printf("\n\n");
	    printf("Steady State Acheived in %d time steps at time %.4f secs\n",iter,t);
	    printf("Rated DeltaP of the pump = %.4f",RCP);
	}

	printf("\n\n");

	mcold = mcdot;
	m1old = m1dot;
	m2old = m2dot;
	t += deltat;
	volFlowRateOld = volFlowRate;
    }
    
    
    //------------------------------------------------------------------------//


    
    //------------------------------------------------------------------------//
    //Dellocators
    free(nData);
    cdeallocator2(&components,nodes); 
    //------------------------------------------------------------------------//

   
    double time2 = MPI_Wtime();
    double secs = time2-time1;
    if(myrank == master)
    {
	printf("Total run time: %.6f secs\n",secs);
/*fprintf(out,"Total run time: %.6f secs\n",secs);
  fclose(out);*/
    }
    
    if(myrank == master)
    {
      printf("----------------------------------------------------------------\n");
      printf("---------------------That's all folks!--------------------------\n");
      printf("----------------------------------------------------------------\n");
    }
    MPI_Finalize();
}
      

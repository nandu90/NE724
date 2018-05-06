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

    //------------------------------------------------------------------------//
    //Read problem data
    problemData();
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Construct the Coefficient matrix for Mass-momentum equation
    
    double a1;
    double b1;
    double a2;
    double b2;
    double ac;
    double bc;

    
    //------------------------------------------------------------------------//

    plugme = 0;

    //------------------------------------------------------------------------//
    //Time loop will come here
    double t=0.0;
    double deltat = iniDeltat*1.0/3600.0;
    int iNewton;

    //------------------------------------------------------------------------//
    //Loop indexes
    int i;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Solution variables
    double m1dot = iniM1dot;
    double m2dot = (iniMcdot - iniM1dot)/(nloops - 1);
    double mcdot = iniMcdot;

    double m1old = m1dot;
    double m2old = m2dot;
    double mcold = mcdot;
    double deltaPcore = 0.0;

    
    double rhosys = iniRho;

    double RCP = 100.0*144.0;
    double volFlowRate=0.0;
    double volFlowRateOld = 0.0;
    double error = 100.0;
    double RCPtol = 1e-8;
    int iter = 0;

    //Temperature Related
    double *T;
    allocator1(&T, 26);
    for(i=0; i<26; i++)
    {
	T[i] = Tsat;
    }

    double *rho;
    allocator1(&rho, 26);
    for(i=0; i<26; i++)
    {
	rho[i] = rhofromT(Tsat);;
    }

    double *mu;
    allocator1(&mu, 26);
    for(i=0; i<26; i++)
    {
	mu[i] = mufromT(Tsat);;
    }

    double *u;
    allocator1(&u, 26);
    for(i=0; i<26; i++)
    {
	u[i] = ufromT(T[i]);
    }

    double *Tclad;
    allocator1(&Tclad, 4);
   
    double mcCheck;
    double change;
    
    //------------------------------------------------------------------------//

    //printf("m1dot m2dot mcdot mu rho = %.4e %.4e %.4e %.4e %.4e %.4e\n",m1dot, m2dot, mcdot, mu, rhosys, nloops);
    
    while(fabs(error) > RCPtol)//(t <= finalTime)//(fabs(mcdot) > 1e-8)
    {
	iter++;
	mcCheck = mcold;
	//Newton's iteration loop
	printf("----------------------------Time Step = %.4e----------------------------\n",t);
	for(iNewton=0; iNewton<maxNewton; iNewton++)
	{
	    loopTerms(&a1, &b1, deltat, nData, m1dot, m1old, rho, mu, rhosys, 1, 0, 0);
	    b1 += RCP;
	    //printf("a1 and b1 outside = %.4f %.4f\n", a1, b1);
	    loopTerms(&a2, &b2, deltat, nData, m2dot, m2old, rho, mu, rhosys, 2, 0, 0);
	    b2 += RCP;
	    //printf("a2 and b2 outside = %.4f %.4f\n", a2, b2);
	    coreTerms(&ac, &bc, deltat, nData, mcdot, mcold, rho, mu, rhosys);
	    //printf("ac and bc outside = %.4f %.4f\n", ac, bc);
	    
	    //Solution from Cramer's loop for next iteration
	    m1dot = (b1*a2 + (nloops-1.0)*b1*ac - (nloops-1.0)*ac*b2 + a2*bc)/(a1*a2 + (nloops-1.0)*a1*ac + a2*ac);
	    deltaPcore = b1 - a1*m1dot;
	    m2dot = (b2 - deltaPcore)/a2;
	    mcdot = m1dot + (nloops-1.0)*m2dot;

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
	error = targetMdot - volFlowRate;
	printf("Volumetric flow Rate in the core = %.2e\n",volFlowRate);
	printf("Difference between target and required value  = %.2e\n",error);
	if(fabs(error) > 1e-8)
	{
	    RCP += 0.1*error + 0.1*(volFlowRate - volFlowRateOld)*deltat;
	    printf("RCP value is %.4e\n",RCP);
	    
	}
	else
	{
	    printf("\n\n");
	    printf("Steady State Acheived in %d time steps at time %.4f secs\n",iter,t);
	    printf("Mass flow rate of Loop1 %.4e lbm/hr \n",m1dot);
	    printf("Mass flow rate of Loop2 %.4e lbm/hr \n",m2dot);
	    printf("Mass flow rate of Core %.4e lbm/hr \n",mcdot);
	    printf("Volumetric flow rate of Core %.4e lbm/hr-ft^2 \n",mcdot/nData[0].Ax);
	    printf("Rated DeltaP of the pump = %.4f psi",RCP/144.0);
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
    printf("required inlet temp is %.4f\n",inletTemp);
    printf("Entering Stage 2\n");
    double eta = 0.0;
    deltat = 1.0/3600.0;
    double ramp = 0.01*power;
    t = 0.0;
    double Qdot;
    iter = 0;
    double massRateError = 100.0;
    double inTempError = 100.0;
    
    while(fabs(massRateError) > 1e-8)// || fabs(inTempError) > 1e-8)
    {
	iter++;
	Qdot = min(power,(t+deltat) * ramp);
	printf("Power is %.4f MW\n",Qdot/(1.0E6*3.412141633));
	
	//------------------------------------------------------------------------//
	//Solve the mass momentum equations
	mcCheck = mcold;
	for(iNewton=0; iNewton<maxNewton; iNewton++)
	{
	    loopTerms(&a1, &b1, deltat, nData, m1dot, m1old, rho, mu, rhosys, 1, 0, 0);
	    b1 += RCP;
	    loopTerms(&a2, &b2, deltat, nData, m2dot, m2old, rho, mu, rhosys, 2, 0, 0);
	    b2 += RCP;
	    coreTerms(&ac, &bc, deltat, nData, mcdot, mcold, rho, mu, rhosys);
	    
	    //Solution from Cramer's loop for next iteration
	    m1dot = (b1*a2 + (nloops-1.0)*b1*ac - (nloops-1.0)*ac*b2 + a2*bc)/(a1*a2 + (nloops-1.0)*a1*ac + a2*ac);
	    deltaPcore = b1 - a1*m1dot;
	    m2dot = (b2 - deltaPcore)/a2;
	    mcdot = m1dot + (nloops-1.0)*m2dot;

	    
	    //exit(1);
	    change = (mcCheck - mcdot)/mcCheck;
	    /*printf("Newton Iteration = %d\n", iNewton+1);
	    printf("Core Mass Flow rate = %.4e\n", mcdot);
	    printf("Relative change in mass core flow rate = %.4e\n", change);
	    printf("Pressure drop across the core = %.4e\n",deltaPcore);*/
	    mcCheck = mcdot;

	    if(fabs(change) < 1e-8)
	    {
		printf("Newton's Iteration converged in %d iterations\n",iNewton+1);
		break;
	    }
	    
	}
	
	volFlowRate = mcdot/nData[0].Ax;
	printf("Volumetric flow rate of Core %.4e lbm/hr-ft^2 \n",volFlowRate);

	//------------------------------------------------------------------------//
	//Solve the heat equation
	//Adjust the eta value after the core has reached full power
	/*if(fabs(massRateError) < 1e-8)
	{
	    inTempError = inletTemp - T[0];
	    printf("Node 0 temp is %.4f, required is %.4f and error is %.4f and eta value is %.4e\n",T[0], inletTemp, inTempError,eta);
	    if(fabs(inTempError) > 1e-8)
	    {
		eta = 2.32E-04;
		//eta = -0.001 * fabs(inTempError) *  inTempError/fabs(inTempError);
	    }
	    }*/
	solveTemp(&m1dot, &m2dot, &mcdot, nData,  deltat, T, rho, mu,u, eta, Qdot, iter, 0, 0);
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Adjust the RCP value once the core reaches full power
	if(fabs(power - Qdot) <1e-8)// && fabs(massRateError) > 1e-8)
	{
	    volFlowRate = mcdot/nData[0].Ax;
	    massRateError = targetMdot - volFlowRate;

	    if(fabs(massRateError) > 1e-8)
	    {
		RCP += 0.001*massRateError + 0.001*(volFlowRate - volFlowRateOld)*deltat;
		printf("RCP value is %.4e\n",RCP/144.0);
		printf("Mass flow rate error is %.4e and Inlet Temperature Error is %.4e\n",massRateError, inTempError);
	    }
	    else
	    {
		printf("Steady State Acheived in %d time steps at time %.4f secs\n",iter,t);
		printf("Mass flow rate of Loop1 %.4e lbm/hr \n",m1dot);
		printf("Mass flow rate of Loop2 %.4e lbm/hr \n",m2dot);
		printf("Mass flow rate of Core %.4e lbm/hr \n",mcdot);
		printf("Volumetric flow rate of Core %.4e lbm/hr-ft^2 \n",mcdot/nData[0].Ax);
		printf("Rated DeltaP of the pump = %.4f psi\n",RCP/144.0);
		printf("Mass flow rate error is %.4e and Inlet Temperature Error is %.4e\n",massRateError, inTempError);
	    }
	}
	//------------------------------------------------------------------------//

	
	if(fabs(massRateError) < 1e-8)// && fabs(inTempError) < 1e-8)
	{
	    for(i=0 ;i <26; i++)
	    {
		printf("Temperature in node %d is %.4f\n",i,T[i]);
	    }

	    getCladTemp(mcdot, T, Tclad, mu, nData, Qdot);
	    for(i=0 ; i<4; i++)
	    {
		printf("Node %d Clad temperature is %.4f and coolant temperature is %.4f\n",i,Tclad[i], T[i]);
	    }
	}

	
	mcold = mcdot;
	m1old = m1dot;
	m2old = m2dot;
	t += deltat;
	volFlowRateOld = volFlowRate;

	printf("\n\n");
    }
    
    
    
    //------------------------------------------------------------------------//
    plugme = 0;
    //File operations
    FILE *Qout;
    Qout = fopen("data.txt","w");
    if(Qout == NULL)
    {
	printf("Error opening power.txt\n");
	exit(1);
    }
    
    //Time realted variables   
    deltat = 0.01/3600.0;
    double totalTime = 0.2/3600.0;
    t = 0.0;
    double trip = 2.0/3600.0;
    //double toperate = 365.0*24.0*60.0*60.0;   //Operating time in secs
    //double tsincetrip;
    iter = 0;

    //Pump trip model
    //double beta = 18.35;
    //double deltaP0 = RCP;

    double Tsatsys = 652.744;
    //double maxclad;

    int nblock = 500;

    double origm1 = m1dot;
    double origm2 = m2dot;
    double origmc = mcdot;
    
    while(t<totalTime)
    {
	iter++;
	printf("Iteration number %d\n",iter);
	//------------------------------------------------------------------------//
	//Determine power value
	if(t<=trip)
	{
	    Qdot = power;
	}
	else
	{
	    Qdot = power;
	    /*tsincetrip = (t - trip)*3600.0;  //Convert to secs
	    Qdot = pow(tsincetrip + 10.0,-0.2);
	    Qdot += -pow(toperate + tsincetrip + 10.0,-0.2);
	    Qdot += 0.87*pow(toperate + tsincetrip + 2.0E7, -0.2);
	    Qdot += -0.87*pow(tsincetrip + 2.0E7, -0.2);
	    Qdot = Qdot*0.1*power;*/
	}
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Solve the mass momentum equations
	mcCheck = mcold;
	for(iNewton=0; iNewton<maxNewton; iNewton++)
	{
	    loopTerms(&a1, &b1, deltat, nData, m1dot, m1old, rho, mu, rhosys, 1, 1, nblock);
	    b1 += RCP;
	    loopTerms(&a2, &b2, deltat, nData, m2dot, m2old, rho, mu, rhosys, 2, 1, nblock);
	    b2 += RCP;
	    coreTerms(&ac, &bc, deltat, nData, mcdot, mcold, rho, mu, rhosys);
	    
	    //Solution from Cramer's loop for next iteration
	    m1dot = (b1*a2 + (nloops-1.0)*b1*ac - (nloops-1.0)*ac*b2 + a2*bc)/(a1*a2 + (nloops-1.0)*a1*ac + a2*ac);
	    deltaPcore = b1 - a1*m1dot;
	    m2dot = (b2 - deltaPcore)/a2;
	    mcdot = m1dot + (nloops-1.0)*m2dot;

	    change = (mcCheck - mcdot)/mcCheck;
	    mcCheck = mcdot;

	    if(fabs(change) < 1e-8)
	    {
		printf("Newton's Iteration converged in %d iterations\n",iNewton+1);
		break;
	    }
	    
	}
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Solve the temperature equation
	solveTemp(&m1dot, &m2dot, &mcdot, nData,  deltat, T, rho, mu,u, eta, Qdot, iter, 1, nblock);
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Get the corresponding temperature of cladding
	getCladTemp(mcdot, T, Tclad, mu, nData, Qdot);
	printf("Cladding and Coolant Temperatures at time %.4f secs\n",t*3600.0);
	for(i=0; i<4; i++)
	{
	    printf("Node %d CladTemp is %.4f and CoolantTemp is %.4f SatTemp is %.4f\n",i,Tclad[i],T[i],Tsatsys);
	}
	printf("Volumetric flow rate of Core %.4e lbm/hr-ft^2 \n",mcdot/nData[0].Ax);
	printf("Rated DeltaP of the pump = %.4f psi\n",RCP/144.0);

	printf("Loop 1: Orig = %.4e and new = %.4e\n", origm1, m1dot);
	printf("Loop 2: Orig = %.4e and new = %.4e\n", origm2, m2dot);
	printf("Core: Orig = %.4e and new = %.4e\n", origmc, mcdot);
	printf("The flow rate is reduced to %.2e%%\n", (origm1 - m1dot)*100.0/origm1);
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Check if max temp exceeds Tsat
	/*maxclad = max(Tclad[0], max(Tclad[1], max(Tclad[2],Tclad[3])));

	if(maxclad > Tsatsys || T[3] > Tsatsys)
	{
	    printf("Max Clad %.4f or Core exit Temp %.4f exceeds Saturation Temp %.4f\n",maxclad,T[3],Tsatsys);
	    printf("Please increment beta %.4f\n",beta);
	    exit(1);
	    }*/
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Print the output
	fprintf(Qout,"%.4e %.4e ",t*3600.0,Qdot/(1.0E6*3.412141633));
	fprintf(Qout,"%.4e %.4e %.4e %.4e ",Tclad[0], Tclad[1], Tclad[2], Tclad[3]);
	fprintf(Qout,"%.4e %.4e ",T[0], T[3]);
	fprintf(Qout,"%.4e %.4e %.4e",mcdot, m1dot, m2dot);
	fprintf(Qout,"\n");
	//------------------------------------------------------------------------//
	
	mcold = mcdot;
	m1old = m1dot;
	m2old = m2dot;
	t += deltat;
	volFlowRateOld = volFlowRate;

	printf("\n\n");
    }

    fclose(Qout);
    //------------------------------------------------------------------------//

    
    
    //------------------------------------------------------------------------//
    //Dellocators
    free(nData);
    cdeallocator2(&components,nodes);

    deallocator1(&T, 26);
    deallocator1(&rho, 26);
    deallocator1(&mu, 26);
    deallocator1(&u, 26);
    deallocator1(&Tclad, 4);
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
      

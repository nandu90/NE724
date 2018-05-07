/***************************************************************************

Author: nsaini
Created: 2018-05-03

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "rhs.h"
#include "correlations.h"
#include "solvers.h"

void getCladTemp(double mcdot, double *T, double *Tclad, double *mu, struct nodeData *nData, double Qdot)
{
    //------------------------------------------------------------------------//
    //Loop indexes and temporary variables
    int i;
    double ktemp, cptemp;
    double Pr, Re;
    double h;
    double dia, area;
    double qdot;
    double harea;
    //------------------------------------------------------------------------//
    double factor[4] = {0.586, 1.414, 1.414, 0.586};

    //------------------------------------------------------------------------//
    //Loop in the core nodes
    for(i=0; i<4; i++)
    {
	ktemp = kfromT(T[i]);
	cptemp = cpfromT(T[i]);

	Pr = mu[i]*cptemp/ktemp;

	dia = nData[i].De;
	area = nData[i].Ax;
	Re = fabs(dia*mcdot/(area*mu[i]));

	//Dittus-Boelter
	h = 0.023*pow(Re,0.8)*pow(Pr,0.4)*ktemp/dia;

	//Total heat transfer in this node
	qdot = factor[i]*Qdot/4.0;

	//Total heat transfer area of the node
	harea = PI*dRod*nData[i].len * nRodFuel;

	//Clad temperature
	Tclad[i] = qdot/(h*harea) + T[i];
	
    }
}

void solveTemp(double *m1dot, double *m2dot, double *mcdot, struct nodeData *nData, double deltat, double *T, double *rho, double *mu, double *u, double eta, double Qdot, int iter, int plug, int nblock)
{
    double m1 = *m1dot;
    double m2 = *m2dot;
    double mc = *mcdot;

    double m2bar = m2*(nloops - 1.0);

    //------------------------------------------------------------------------//
    //Loop indexes
    int i;
    double a,b,c;
    double vol;
    
    double c1, c2;
    double a1, a2;
    double area;
    //------------------------------------------------------------------------//

    
    
    //------------------------------------------------------------------------//
    //Assemble the matrix
    // 0-3 : Core
    // 4-13: Loop 1
    // 14-23: Loop 2,3,4
    // 24 : Lower Plenum
    // 25 : Upper Plenum

    double **A;
    allocator2(&A, 26, 26);

    //------------------------------------------------------------------------//
    //For the core
    for(i=0; i<4; i++)
    {
	vol = nData[i].len * nData[i].Ax;
	
	a = -(mc/2.0)*(1.0 + fabs(mc)/mc);
	b = vol*rho[i]/deltat + fabs(mc);
	c = (mc/2.0)*(1.0 - fabs(mc)/mc);

	if(i == 0)
	{
	    A[i][i] += b;
	    A[i][i+1] += c;
	    A[i][24] += a;    //Goes to Lower Plenum
	}
	else if(i == 3)
	{
	    A[i][i] += b;
	    A[i][25] += c;
	    A[i][i-1] += a;   //Goes to upper Plenum
 	}
	else
	{
	    A[i][i] += b;
	    A[i][i+1] += c;
	    A[i][i-1] += a;
	}
	//printf("%.2e %.2e %.2e\n",A[i][i-1],A[i][i],A[i][i+1]);
    }

    //For Loop 1: Hot Leg to downcomer
    for(i=4; i<14; i++)
    {
	area = nData[i+1].Ax;

	if(i>=5 && i<=11 && plug == 1)
	{
	    //area = area/nSGTubes;
	    //area = area*((double)(nSGTubes - nblock));
	}
	
	vol = nData[i+1].len * area;

	a = -(m1/2.0)*(1.0 + fabs(m1)/m1);
	b = vol*rho[i]/deltat + fabs(m1);
	c = (m1/2.0)*(1.0 - fabs(m1)/m1);

        if(i == 4)
	{
	    A[i][i] += b;
	    A[i][i+1] += c;
	    A[i][25] += a;    //Goes to Upper Plenum
	}
	else if(i == 13)
	{
	    A[i][i] += b;
	    A[i][24] += c;    //Goes to Lower Plenum
	    A[i][i-1] += a;   
 	}
	else
	{
	    A[i][i] += b;
	    A[i][i+1] += c;
	    A[i][i-1] += a;
	}
    }

    //For Loop 2,3,4: Hot Leg to downcomer
    for(i=14; i<24; i++)                             //DOUBT should i take m2bar or m2dot
    {
	vol = nData[i-9].len * nData[i-9].Ax;

	a = -(m2/2.0)*(1.0 + fabs(m2)/m2);
	b = vol*rho[i]/deltat + fabs(m2);
	c = (m2/2.0)*(1.0 - fabs(m2)/m2);

        if(i == 14)
	{
	    A[i][i] += b;
	    A[i][i+1] += c;
	    A[i][25] += a;    //Goes to Upper Plenum
	}
	else if(i == 23)
	{
	    A[i][i] += b;
	    A[i][24] += c;    //Goes to Lower Plenum
	    A[i][i-1] += a;   
 	}
	else
	{
	    A[i][i] += b;
	    A[i][i+1] += c;
	    A[i][i-1] += a;
	}
    }

    //Upper Plenum
    vol = UPVol;
    
    a = -(mc/2.0)*(1.0 + fabs(mc)/mc);
    b = vol*rho[25]/deltat + 0.5*(fabs(mc) + fabs(m1) + fabs(m2bar));
    c1 = (m1/2.0)*(1.0 - fabs(m1)/m1);
    c2 = (m2bar/2.0)*(1.0 - fabs(m2bar)/m2bar);
    A[25][25] = b;
    A[25][3] = a;
    A[25][4] = c1;
    A[25][14] = c2;

    //Lower Plenum
    vol = LPVol;
    
    a1 = -(m1/2.0)*(1.0 + fabs(m1)/m1);
    a2 = -(m2bar/2.0)*(1.0 + fabs(m2bar)/m2bar);
    b = vol*rho[24]/deltat + 0.5*(fabs(mc) + fabs(m1) + fabs(m2bar));
    c = (mc/2.0)*(1.0 - fabs(mc)/mc);
    A[24][24] = b;
    A[24][0] = c;
    A[24][13] = a1;
    A[24][23] = a2;

    /*int j;
    if(iter == 4 && plug == 1)
    {
	for(i=0; i<26; i++)
	{
	    printf("%d ",i);
	    for(j=0; j<26; j++)
	    {
		if(fabs(A[i][j]) < 1e-10)
		{
		    printf("* ");
		}
		else
		{
		    printf("%.1e ",A[i][j]);
		}
	    }
	    printf("\n");
	}
	}*/
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //RHS for Steam Generators -  Heat transfer term
    double Re, Pr;
    double ktemp, cptemp;
    double dia;
    
    double *B;
    double *hcoeff;
    double *UA;
    double *qdot;
    allocator1(&B, 26);
    allocator1(&hcoeff, 26);
    allocator1(&UA, 26);
    allocator1(&qdot, 26);

    //Loop 1 Steam Generators
    for(i=5; i<12; i++)
    {
	ktemp = kfromT(T[i]);
	cptemp = cpfromT(T[i]);
	Pr = mu[i]*cptemp/ktemp;

	dia = nData[i+1].De;
	area = nData[i+1].Ax;
	if(plug == 1)
	{
	    //area = area/nSGTubes;
	    //area = area*((double)(nSGTubes - nblock));
	    
	}
	Re = fabs(dia*m1/(area*mu[i]));

	//Dittus-Boelter
	hcoeff[i] = 0.023*pow(Re,0.8)*pow(Pr,0.3)*ktemp/dia;
	
	if(plug == 1)
	{
	    UA[i] = (nSGTubes-nblock)*PI*dia*nData[i+1].len/((1.0/hcoeff[i]) + eta);
	}
	else
	{
	    UA[i] = nSGTubes*PI*dia*nData[i+1].len/((1.0/hcoeff[i]) + eta);
	}
	qdot[i] = UA[i]*(Tsat - T[i]);

	/*if(iter == 5 && plug == 1)
	{
	    printf("%.4e %.4e %.4e %.4e\n",Pr,mu[i],cptemp, T[i]);
	    }*/
    }
    //if(iter == 5 && plug == 1)exit(1);

    //Loop 2,3,4 Steam generators
    for(i=15; i<22; i++)
    {
	ktemp = kfromT(T[i]);
	cptemp = cpfromT(T[i]);
	Pr = fabs(mu[i]*cptemp/ktemp);

	dia = nData[i-9].De;
	area = nData[i-9].Ax;
	Re = fabs(dia*m2/(area*mu[i]));

	//Dittus-Boelter
	hcoeff[i] = 0.023*pow(Re,0.8)*pow(Pr,0.3)*ktemp/dia;
	UA[i] = nSGTubes*PI*dia*nData[i-9].len/((1.0/hcoeff[i]) + eta);
	qdot[i] = UA[i]*(Tsat - T[i]);                 //Doubt here
    }
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Call solver for heat conduction in the core here
    
    
    double factor[4] = {0.586, 1.414, 1.414, 0.586};
    for(i=0; i<4; i++)
    {
	qdot[i] = factor[i]*Qdot/4.0;
	//printf("Power for %d is %.4e\n",i,qdot[i]);
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Assemble RHS now
    for(i=0; i<26; i++)
    {
	if(i>=0 && i<=3) //core
	{
	    vol = nData[i].len*nData[i].Ax;
	}
	else if(i>=4 && i<=13)
	{
	    vol = nData[i+1].len*nData[i+1].Ax;
	}
	else if(i>=14 && i<=23)
	{
	    vol = nData[i-9].len*nData[i-9].Ax;
	}
	else if(i==24)
	{
	    vol = LPVol;
	}
	else
	{
	    vol = UPVol;
	}

	B[i] = qdot[i] + vol*rho[i]*u[i]/deltat;

	/*if(iter == 5 && plug == 1)
	{
	    printf("%d RHS is %.4e and qdot is %.4e\n",i, B[i],qdot[i]);
	    }*/
	
    }
    //if(iter == 5 && plug == 1)exit(1);
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Solve for Internal energy value at the next time step
    double *unext;
    allocator1(&unext , 26);
    solveSystem(A, B, unext, 26);
    /*for(i=0; i<26; i++)
    {
	printf("Old: %.4e and New: %.4e\n",u[i], unext[i]);
	}*/

    for(i=0; i<26; i++)
    {
	u[i] = unext[i];
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Compute the temperature from internal energy
    for(i=0; i<26; i++)
    {
	if(plug == 1)printf("Old Temperature was %.2f ",T[i]);
	T[i] = TfromU(u[i]);
	if(plug == 1)printf("and New Temperature is %.2f\n",T[i]);
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Compute the density and viscosity from the new temperature values
    for(i=0; i<26; i++)
    {
	rho[i] = rhofromT(T[i]);
	mu[i] = mufromT(T[i]);
    }
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator2(&A, 26, 26);
    deallocator1(&B, 26);
    deallocator1(&hcoeff, 26);
    deallocator1(&UA, 26);
    deallocator1(&qdot, 26);
    deallocator1(&unext, 26);
    //------------------------------------------------------------------------//

    
}

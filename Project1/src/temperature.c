/***************************************************************************

Author: nsaini
Created: 2018-05-03

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "rhs.h"

void solveTemp(double *m1dot, double *m2dot, double *mcdot, struct nodeData *nData, double deltat)
{
    double m1 = *m1dot;
    double m2 = *m2dot;
    double mc = *mcdot;

    double m2bar = m2*(nloops - 1.0);

    //------------------------------------------------------------------------//
    //Loop indexes
    int i,j;
    double a,b,c;
    double vol;
    double rho =1.0;   //Will have to do something about density
    double c1, c2;
    double a1, a2;
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
	b = vol*rho/deltat + fabs(mc);
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
	vol = nData[i+1].len * nData[i+1].Ax;

	a = -(m1/2.0)*(1.0 + fabs(m1)/m1);
	b = vol*rho/deltat + fabs(m1);
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
    for(i=14; i<24; i++)
    {
	vol = nData[i-9].len * nData[i-9].Ax;

	a = -(m2bar/2.0)*(1.0 + fabs(m2bar)/m2bar);
	b = vol*rho/deltat + fabs(m2bar);
	c = (m2bar/2.0)*(1.0 - fabs(m2bar)/m2bar);

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
    vol = nData[4].len * nData[4].Ax;
    a = -(mc/2.0)*(1.0 + fabs(mc)/mc);
    b = vol*rho/deltat + 0.5*(fabs(mc) + fabs(m1) + fabs(m2bar));
    c1 = (m1/2.0)*(1.0 - fabs(m1)/m1);
    c2 = (m2bar/2.0)*(1.0 - fabs(m2bar)/m2bar);
    A[25][25] = b;
    A[25][3] = a;
    A[25][4] = c1;
    A[25][14] = c2;

    //Lower Plenum
    vol = nData[14].len * nData[14].Ax;
    a1 = -(m1/2.0)*(1.0 - fabs(m1)/m1);
    a2 = -(m2bar/2.0)*(1.0 - fabs(m2bar)/m2bar);
    b = vol*rho/deltat + 0.5*(fabs(mc) + fabs(m1) + fabs(m2bar));
    c = (mc/2.0)*(1.0 + fabs(mc)/mc);
    A[24][24] = b;
    A[24][0] = c;
    A[24][13] = a1;
    A[24][23] = a2;

    for(i=0; i<26; i++)
    {
	for(j=0; j<26; j++)
	{
	    if(fabs(A[i][j]) < 1e-10)
	    {
		printf("- ");
	    }
	    else
	    {
		printf("%.1e ",A[i][j]);
	    }
	}
	printf("\n");
    }
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator2(&A, 26, 26);
    //------------------------------------------------------------------------//

    
}

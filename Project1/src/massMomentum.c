/***************************************************************************

Author: nsaini
Created: 2018-04-15

***************************************************************************/


#include "common.h"
#include "correlations.h"
#include "rhs.h"

void getgeom(double *len, double *area, double *dia, double *deltaH, int index, struct nodeData *nData, double *eqLD, double *n)
{
    (*len) = nData[index].len;
    (*area) = nData[index].Ax;
    (*dia) = nData[index].De;
    (*deltaH) = nData[index].DeltaH;

    if(index == 5)                    //Hot Leg
    {
	(*eqLD) = HLEqLD;
	*n = 1.0;
    }
    else if(index >= 6 && index <=12) //Steam Generators
    {
	*n = 6633.0;
	(*eqLD) = (*len)/(*dia) + BendEqLD;
    }
    else if(index == 13)              //Cold Leg
    {
	(*eqLD) = CLEqLD;
	*n = 1.0;
    }
    else if(index == 14)                             //Downcomer
    {
	(*eqLD) = (*len)/(*dia);
	*n = 1.0;
    }
}

void loopTerms(double *a, double *b, double deltat, struct nodeData *nData, double mdot, double mold, double *rhoarr, double *muarr, double rhosys, int code)
{
    //------------------------------------------------------------------------//
    double len, area, dia, deltaH;
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    double term4 = 0.0;
    
    double F;
    double kin, kout;
    double eqLD;
    double n;

    double rho, mu;
    //Assemble Loop : Hot leg to Downcomer
    //5    : Hot leg
    //6-12 : Steam Generators
    //13   : Cold Leg
    //14   : Downcomer
    int index;
    for(index=5; index < 15; index++)
    {
	if(code == 1)
	{
	    rho = rhoarr[index - 1];
	    mu = muarr[index -1];
	}
	else
	{
	    rho = rhoarr[index + 9];
	    mu = muarr[index + 9];
	}
	getgeom(&len,&area,&dia, &deltaH,  index,nData, &eqLD, &n);
	term1 += len/(area*deltat);
	
	//Get the friction factor
	F = frictionFactor(mdot, mu, dia, area,index);
	term2 += (F *(eqLD))*pow(1.0/area,2.0);
	//term2 += (F *len/dia)*pow(1.0/area,2.0);
	
	if(index == 5)
	{
	    kin = HLInLossCoeff;
	    kout = HLOutLossCoeff;
	}
	else if(index >=6 && index <=12)
	{
	    kin = BundInLossCoeff;
	    kout = BundOutLossCoeff;
	}
	else if(index == 13)
	{
	    kin = CLInLossCoeff;
	    kout = CLOutLossCoeff;
	}
	else
	{
	    kin = 0.0;
	    kout = 0.0;
	}
	term3 += kin*(pow(1.0/area,2.0));
	term3 += kout*(pow(1.0/area,2.0));

	term4 += rho*GRAVITY*deltaH;

	
    }
    term1 = term1/GC;
    term2 = term2/(2.0*rhosys*GC);
    term3 = term3/(2.0*rhosys*GC);
    term4 = term4/GC;

    
    *a = term1 + 2.0*(term2 + term3)*fabs(mdot);
    *b = term1*mold + (term2 + term3)*mdot*fabs(mdot) - term4;

}


void coreTerms(double *a, double *b, double deltat, struct nodeData *nData, double mdot, double mold, double *rhoarr, double *muarr, double rhosys)
{
    double len, area, dia, deltaH;
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    double term4 = 0.0;
    
    double F;
    double kin, kout;

    double eqLD;
    double n;

    int coreIndex;

    int index[6] = {0,1,2,3,4,15};
    int i;

    double rho, mu;

    for(coreIndex =0; coreIndex<4; coreIndex++)
    {
	mu = muarr[coreIndex];
	
	getgeom(&len,&area,&dia, &deltaH,  coreIndex,nData, &eqLD, &n);
	
	
	term1 += len/(area*deltat);
	
	F = frictionFactor(mdot, mu, dia, area,coreIndex);

	term2 += (F *(len/dia))*pow(1.0/area,2.0);
	
		
    }

    
    //Account for inlet and outlet coefficients - Plenums
    kin = cInLossCoeff;
    kout = cOutLossCoeff;
    
    term3 += kin*(pow(1.0/area,2.0));
    term3 += kout*(pow(1.0/area,2.0));
    
    //Add spacer grids
    term3 += ((double)nSpacerGrid)*gridLossCoeff*(pow(1.0/area,2.0));

    //------------------------------------------------------------------------//
    //Get the buoyancy term
    
    
    for(i=0; i<6; i++)
    {
	if(i < 4)
	{
	    rho = rhoarr[i];
	}
	else if(i == 4)
	{
	    rho = rhoarr[25];
	}
	else if(i == 5)
	{
	    rho = rhoarr[24];
	}
	
	getgeom(&len,&area,&dia,&deltaH,index[i],nData, &eqLD, &n);
	term4 += rho*GRAVITY*deltaH;

    }


    term1 = term1/GC;
    term2 = term2/(2.0*rhosys*GC);
    term3 = term3/(2.0*rhosys*GC);
    term4 = term4/GC;
    
    *a = term1 + 2.0*(term2 + term3)*fabs(mdot);
    *b = term1*mold + (term2 + term3)*mdot*fabs(mdot) - term4;

}
    

/***************************************************************************

Author: nsaini
Created: 2018-04-15

***************************************************************************/


#include "common.h"
#include "correlations.h"
#include "rhs.h"

void getgeom(double *len, double *area, double *dia, double *deltaH, int index, struct nodeData *nData)
{
    *len = nData[index].len;
    *area = nData[index].Ax;
    *dia = nData[index].De;
    *deltaH = nData[index].DeltaH;
}

void loopTerms(double *a, double *b, double deltat, struct nodeData *nData, double mdot, double mold, double rho, double mu, double rhoold)
{
    //------------------------------------------------------------------------//
    double len, area, dia, deltaH;
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    double term4 = 0.0;
    
    double F;
    double kin, kout;
    //Assemble Loop : Hot leg to Downcomer
    //5    : Hot leg
    //6-12 : Steam Generators
    //13   : Cold Leg
    //14   : Downcomer
    int index;
    for(index=5; index < 15; index++)
    {
	getgeom(&len,&area,&dia, &deltaH,  index,nData);
	term1 += len/(area*deltat);
	
	//Get the friction factor
	F = frictionFactor(mdot, mu, dia, area,index);
	term2 += (F*len/(dia))*pow(1/area,2.0);
	
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
	term3 += kin*(pow(1/area,2.0));
	term3 += kout*(pow(1/area,2.0));

	term4 += rhoold*GRAVITY*deltaH;
	
    }
    term2 = term2/(2.0*rho);
    term3 = term3/(2.0*rho);

    //printf("values %.4f %.4f %.4f %.4f\n",term1, term2, term3, term4);
    
    
    *a = term1 + 2.0*(term2 + term3)*fabs(mdot);
    *b = term1*mold + (term2 + term3)*mdot*fabs(mdot) - term4;

    //printf(" a and b = %.4f %.4f\n", *a, *b);
    //printf("Buoyancy term in loop = %.4e\n",term4);

    
}


void coreTerms(double *a, double *b, double deltat, struct nodeData *nData, double mdot, double mold, double rho, double mu, double rhoold)
{
    double len, area, dia, deltaH;
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    double term4 = 0.0;
    
    double F;
    double kin, kout;

    getgeom(&len,&area,&dia, &deltaH,  0,nData);

    term1 = len/(area*deltat);

    F = frictionFactor(mdot, mu, dia, area,0);

    term2 = F*len/(dia*area*area*2.0*rho);

    kin = cInLossCoeff;
    kout = cOutLossCoeff;

    term3 = (kin + kout)/(2.0*rho*area*area);

    *a = term1 + 2.0*(term2 + term3)*fabs(mdot);
    

    //------------------------------------------------------------------------//
    //Get the buoyancy term
    int index[6] = {0,1,2,3,4,15};
    int i;
    for(i=0; i<6; i++)
    {
	getgeom(&len,&area,&dia,&deltaH,index[i],nData);
	term4 += rhoold*GRAVITY*deltaH;
    }

    *b = term1*mold + (term2 + term3)*mdot*fabs(mdot) - term4;
    
    //printf("Buoyancy term in core = %.4e\n",term4);
}
    

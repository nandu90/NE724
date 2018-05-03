/***************************************************************************

Author: nsaini
Created: 2018-04-15

***************************************************************************/

#include "common.h"
#include "correlations.h"

double frictionFactor(double mdot, double mu, double D, double area, int index)
{
    double f;
    double Re;
    Re = fabs(D*mdot/(area*mu));

    //printf("Reynolds number in %s is %.2e\n",components[index],Re);

    if(Re <= 500.0)
    {
	f = 0.128;
    }
    else if(Re > 500.0 && Re <= 2300.0)
    {
	f = 64.0/Re;
    }
    else if(Re > 2300.0 && Re <= 4200.0)
    {
	f = (6.0526E-6)*Re + 0.01388;
    }
    else if(Re > 4200.0 && Re <= 30000.0)
    {
	f = 0.3164*pow(Re, -0.25);
    }
    else
    {
	f = 0.184*pow(Re, -0.2);
    }

    return f;
}

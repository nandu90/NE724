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

    printf("Reynolds number in %s is %.2e\n",components[index],Re);

    double f1, f2;
    if(Re >= 640.0 && Re <= 2300.0)
    {
	f = 64.0/Re;
    }
    else if(Re >=4200.0)
    {
	f = 0.316/(pow(Re,0.25)); //Darcy's friction factor: http://www.kolumbus.fi/jukka.kiijarvi/clunowa/fluid_mechanics/pdf_articles/darcy_friction_factor.pdf
    }
    else if(Re > 2300.0 && Re <4200.0)
    {
	//Interpolate for values between 2300 and 4200
	f1 = 64.0/2300.0;
	f2 = 0.316/(pow(4200.0,0.25));

	f = f1 + (Re-2300.0)*(f2-f1)/(4200.0-2300.0);
    }
    else
    {
	f = 0.1;
    }

    return f;
}

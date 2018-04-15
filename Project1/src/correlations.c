/***************************************************************************

Author: nsaini
Created: 2018-04-15

***************************************************************************/

#include "common.h"
#include "correlations.h"

double frictionFactor(double mdot, double mu, double Pw)
{
    double f;
    double Re;
    Re = 4.0*mdot/(Pw*mu);

    double f1, f2;
    if(Re <= 2300.0)
    {
	f = 64.0/Re;
    }
    else if(Re >=4200.0)
    {
	f = 0.316/(pow(Re,0.25)); //Darcy's friction factor: http://www.kolumbus.fi/jukka.kiijarvi/clunowa/fluid_mechanics/pdf_articles/darcy_friction_factor.pdf
    }
    else
    {
	//Interpolate for values between 2300 and 4200
	f1 = 64.0/2300.0;
	f2 = 0.316/(pow(4200.0,0.25));

	f = f1 + (Re-2300.0)*(f2-f1)/(4200.0-2300.0);
    }

    return f;
}

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

double ufromT(double T)
{
    double a,b,c,d;
    double value;
    a = 1.536E-05;
    b = -0.02423;
    c = 13.89;
    d = -2324.0;

    value = a*T*T*T + b*T*T + c*T + d;

    return value;
}

double TfromU(double u)
{
    double a,b,c,d;
    double value;
    a = -2.229E-06;
    b = 0.002936;
    c = -0.3845;
    d = 252.3;

    value = a*u*u*u + b*u*u + c*u + d;

    if(plugme == 1)value = u;
    return value;
}

double rhofromT(double T)
{
    double a,b,c,d;
    double value;
    a = -1.849E-06;
    b = 0.00292;
    c = -1.596;
    d = 349.0;

    value = a*T*T*T + b*T*T + c*T + d;

    return value;
}


double mufromT(double T)
{
    double a,b,c,d;
    double value;
    a = -2.711E-12;
    b = 4.571E-09;
    c = -2.714E-06;
    d = 0.0006234;

    value = a*T*T*T + b*T*T + c*T + d;

    value = value*3600.0;
    
    return value;
}

double kfromT(double T)
{
    double a,b,c,d;
    double value;
    a = -2.883E-10;
    b = -2.248E-06;
    c = 0.00181;
    d = 0.3153;

    value = a*T*T*T + b*T*T + c*T + d;

    value = value*0.5781759824;
    
    return value;
}

double cpfromT(double T)
{
    double a,b,c,d;
    double value;
    a = 2.12E-07;
    b = -0.0003301;
    c = 0.1726;
    d = -29.13;

    value = a*T*T*T + b*T*T + c*T + d;

    value = 1.5;
    return value;
}

double cvfromT(double T)
{
    //double a,b,c,d,e;
    double value;
    /*a = 2.499E-10;
    b = -5.491E-07;
    c = 0.0004529;
    d = -0.1665;
    e = 23.77;

    value = a*T*T*T*T + b*T*T*T + c*T*T + d*T + e;*/

    value = 0.7307638889;

    return value;
}

/***************************************************************************

Author: nsaini
Created: 2018-05-05

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "rhs.h"

void conduction(double *T, double Qdot)
{
    //Volumetric heat generation
    double fuelarea = nRodFuel * PI* (pow(dPellet,2.0)/4.0)*12.0;
    double qvol = Qdot/fuelarea;

    //------------------------------------------------------------------------//
    //Deltax for conduction
    //Have atleast 10 cells in the gap
    double gap = (dRod/2.0) - cladThic - dPellet;
    double deltax = gap/10.0;

    int ncells = dRod/deltax;

    printf("%.4e %.4e %d",gap, deltax, ncells);

    exit(1);
    
    
}

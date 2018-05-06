/***************************************************************************

Author: nsaini
Created: 2018-04-16

***************************************************************************/


#ifndef RHS_H
#define RHS_H

void getgeom(double *, double *, double *, double *, int, struct nodeData *, double *, double *);

void loopTerms(double *, double *, double, struct nodeData *, double, double, double *, double *, double, int);

void coreTerms(double *, double *, double, struct nodeData *, double, double, double *, double *, double);

void solveTemp(double *, double *, double *, struct nodeData *, double, double *, double * , double *, double *, double, double, int);

void getCladTemp(double , double *, double *, double *, struct nodeData *, double );
#endif

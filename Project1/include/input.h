/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#ifndef INPUT_H
#define INPUT_H

//------------------------------------------------------------------------//
//The prblem data is hardcoded for now. No need going through the pain of
//reading everything.

//REACTOR
double power; //in watts
double cHeight;    //in inches (why!!! :( )
int nRodLoc;
int nRodFuel;
double dRod;
double dPellet;
double cladThic;
double gapConduct;
double cladConduct;
double rodPitch;
int nSpacerGrid;
double gridLossCoeff;
double massFlux;
double inletTemp;
double outletTemp;
double pressure;
double cInLossCoeff;
double cOutLossCoeff;

//UPPER PLENUM
double UPLen;
double UPDia;
double UPVol;

//HOT LEGS
int nHL;
double HLLen;
double HLDia;
double HLEqLD;
double HLInLossCoeff;
double HLOutLossCoeff;

//STEAM GENERATORS
int nSG;
double secPressure;
int nSGTubes;
double SGTubeInDia;
double SGTubeLen;
double BendEqLD;
double BundInLossCoeff;
double BundOutLossCoeff;

//COLD LEGS
int nCL;
double CLLen;
double CLDia;
double CLInLossCoeff;
double CLOutLossCoeff;
double CLEqLD;

//DOWNCOMER
double DCInDia;
double DCOutDia;
double DCLen;

//LOWER PLENUM
double LPLen;
double LPDia;
double LPVol;


//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Inputs
double iniP;
double iniMu;
double iniRho;
double iniMcdot;
double iniM1dot;

double finalTime;
double maxNewton;
double iniDeltat;
//------------------------------------------------------------------------//







#endif

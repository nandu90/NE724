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
double power = 3411.0E6; //in watts
double cHeight = 144.0;    //in inches (why!!! :( )
int nRodLoc = 55777;
int nRodFuel = 50952;
double dRod = 0.374;
double dPellet = 0.3225;
double cladThic = 0.0225;
double gapConduct = 1000.0;
double cladConduct = 9.6;
double rodPitch = 0.496;
int nSpacerGrid = 8;
double gridLossCoeff = 0.5;
double massFlux = 2.48E6;
double inletTemp = 552.0;
double outletTemp = 616.0;
double pressure = 2250;
double cInLossCoeff = 4.25;
double cOutLossCoeff = 4.25;

//UPPER PLENUM
double UPLen = 1.5;
double UPDia = 158;
double UPVol = 1373.7;

//HOT LEGS
int nHL = 4;
double HLLen = 20.0;
double HLDia = 29.0;
double HLEqLD = 20.0;
double HLInLossCoeff = 0.5;
double HLOutLossCoeff = 1.0;

//STEAM GENERATORS
int nSG = 4;
double secPressure = 1000.0;
int nSGTubes = 6633;
double SGTubeInDia = 0.6075;
double SGTubeLen = 66.8;
double BendEqLD = 55.0;
double BundInLossCoeff = 0.5;
double BundOutLossCoeff = 1.0;

//COLD LEGS


//All input variables to be read from the control file are specified here







#endif

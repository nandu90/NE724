/***************************************************************************

Author: nsaini
Created: 2018-04-15

***************************************************************************/


#include "common.h"
#include "fileIO.h"

void geomData(struct nodeData* nData)
{
    if(myrank == master)
    {
	FILE *datafile;
	datafile = fopen("geomData.txt","r");
	if(datafile == NULL)
	{
	    printf("Error opening Geom Data file!\n");
	    exit(1);
	}

	char* line = NULL;
	ssize_t size;
	size_t len = 0;
	
	char delim [1] = " ";
	char* word = NULL;

	int index = 0;

	printf("Number of nodes = %d\n\n",nodes);
	while((size = getline(&line, &len, datafile)) != -1)
	{
	    word = strtok(line, delim);
	    while(word != NULL)
	    {		
		word = strtok(NULL,delim);
		
		word = strtok(NULL,delim);
		nData[index].len = atof(word);
		
		word = strtok(NULL,delim);
		nData[index].De = atof(word);
		
		word = strtok(NULL,delim);
		nData[index].Ax = atof(word);
		
		word = strtok(NULL,delim);
		nData[index].Ph = atof(word);
		
		word = strtok(NULL,delim);
		nData[index].DeltaH = atof(word);
		
		index++;
		word = strtok(NULL,delim);
	    }
	}

	fclose(datafile);
	
	free(word);
	free(line);
	
	//Check
	int i;
	for(i=0; i<nodes; i++)
	{
	    printf("%d %.2f %.2f %.2f %.2f %.2f\n\n",i, nData[i].len, nData[i].De, nData[i].Ax, nData[i].Ph, nData[i].DeltaH);
	}
    }
}


void problemData()
{
     power = 3411.0E6; //in watts
     cHeight = 144.0;    //in inches (why!!! :( )
     nRodLoc = 55777;
     nRodFuel = 50952;
     dRod = 0.374;
     dPellet = 0.3225;
     cladThic = 0.0225;
     gapConduct = 1000.0;
     cladConduct = 9.6;
     rodPitch = 0.496;
     nSpacerGrid = 8;
     gridLossCoeff = 0.5;
     massFlux = 2.48E6;
     inletTemp = 552.0;
     outletTemp = 616.0;
     pressure = 2250;
     cInLossCoeff = 4.25;
     cOutLossCoeff = 4.25;
    
//UPPER PLENUM
     UPLen = 1.5;
     UPDia = 158;
     UPVol = 1373.7;
    
//HOT LEGS
     nHL = 4;
     HLLen = 20.0;
     HLDia = 29.0;
     HLEqLD = 20.0;
     HLInLossCoeff = 0.5;
     HLOutLossCoeff = 1.0;
    
//STEAM GENERATORS
     nSG = 4;
     secPressure = 1000.0;
     nSGTubes = 6633;
     SGTubeInDia = 0.6075;
     SGTubeLen = 66.8;
     BendEqLD = 55.0;
     BundInLossCoeff = 0.5;
     BundOutLossCoeff = 1.0;
    
//COLD LEGS
     nCL = 4;
     CLLen = 40.0;
     CLDia = 27.5;
     CLInLossCoeff = 0.5;
     CLOutLossCoeff = 4.6;
     CLEqLD = 18.0;
    
//DOWNCOMER
     DCInDia = 173.0;
     DCOutDia = 158.0;
     DCLen = 18.4;
    
//LOWER PLENUM
     LPLen = 7.2;
     LPDia = 173.0;
     LPVol = 784.4;
}

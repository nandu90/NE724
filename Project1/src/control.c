/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#include "fileIO.h"
#include "common.h"

void control()
{
  if(myrank == master)
    {
	printf("Reading the input file\n\n");
	
	FILE *controlfile;
	controlfile = fopen("control.txt","r");
	if(controlfile == NULL)
	{
	    printf("Error opening control file!\n");
	    exit(1);
	}

	char* line = NULL;
	ssize_t size;
	size_t len = 0;
	
	char delim [1] = " ";
	char* word = NULL;


	while((size = getline(&line, &len, controlfile)) != -1)
	{
	    word = strtok(line, delim);
	    while(word != NULL)
	    {
		if(strcmp(word,"Nodes") == 0)
		{
		    word = strtok(NULL,delim);
		    nodes = atoi(word);
		}
		if(strcmp(word,"Number_of_loops") == 0)
		{
		    word = strtok(NULL,delim);
		    nloops = atof(word);
		}
		if(strcmp(word,"System_pressure") == 0)
		{
		    word = strtok(NULL,delim);
		    iniP = atof(word);
		}
		if(strcmp(word,"Density") == 0)
		{
		    word = strtok(NULL,delim);
		    iniRho = atof(word);
		}
		if(strcmp(word,"Viscosity") == 0)
		{
		    word = strtok(NULL,delim);
		    iniMu = atof(word);
		}
		if(strcmp(word,"Core_Flow_Rate") == 0)
		{
		    word = strtok(NULL,delim);
		    iniMcdot = atof(word);
		}
		if(strcmp(word,"Loop1_Flow_Rate") == 0)
		{
		    word = strtok(NULL,delim);
		    iniM1dot = atof(word);
		}
		if(strcmp(word,"Totaltime") == 0)
		{
		    word = strtok(NULL,delim);
		    finalTime = atof(word);
		}
		if(strcmp(word,"Max_Newton_Iterations") == 0)
		{
		    word = strtok(NULL,delim);
		    maxNewton = atof(word);
		}
		if(strcmp(word,"Time_Step_Size") == 0)
		{
		    word = strtok(NULL,delim);
		    iniDeltat = atof(word);
		}
		if(strcmp(word,"Target_Mass_Flow_Rate") == 0)
		{
		    word = strtok(NULL,delim);
		    targetMdot = atof(word);
		}
		word = strtok(NULL,delim);
	    }
	}

	fclose(controlfile);
	
	free(word);
	free(line);
    }
    
}

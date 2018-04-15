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
		word = strtok(NULL,delim);
	    }
	}

	fclose(controlfile);
	
	free(word);
	free(line);
    }
    
}

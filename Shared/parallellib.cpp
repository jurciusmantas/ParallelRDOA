#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
#include "RDOAlib.h"
#include "populationlib.h"

void initPopulationStructToMPI(MPI_Datatype* population_dt, int numX)
{
	MPI_Type_contiguous(numX + 1, MPI_DOUBLE, &(*population_dt));
  	MPI_Type_commit(&(*population_dt));
}

void calculateDistancesAsync(int numDP, int numProcs, int id, double** distances, double** demandPoints)
{
    std::cout << id << ": calculateDistances START" << std::endl;

	//Allocate memory and fill 0's
	*distances = (double*) calloc(numDP * numDP, sizeof(double));;

	int procChunkSizeForCalc = numDP / numProcs;
	int offsetForCalc = procChunkSizeForCalc * id;

		for (int iIters = offsetForCalc; iIters < offsetForCalc + procChunkSizeForCalc; iIters++)
		{
			for (int jIters = iIters + 1; jIters < numDP; jIters++)
			{
				double distance = HaversineDistance(demandPoints[iIters], demandPoints[jIters]);
				(*distances)[(iIters * numDP) + jIters] = distance;
			}
		}
		
	MPI_Allgather(*(distances)+(offsetForCalc * numDP), procChunkSizeForCalc * numDP, 
		MPI_DOUBLE, *distances, procChunkSizeForCalc * numDP, MPI_DOUBLE, MPI_COMM_WORLD);
}

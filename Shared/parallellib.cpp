#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>
#include "RDOAlib.h"
#include "populationlib.h"

void initPopulationStructToMPI(MPI_Datatype* population_dt, int numX)
{
	MPI_Aint displacements[2]  = { offsetof(populationItem, solution), offsetof(populationItem, locations) };
  	int block_lengths[2]  = { 1, numX };
  	MPI_Datatype types[2] = { MPI_DOUBLE, MPI_INT };

  	int createStructRes = MPI_Type_create_struct(2, block_lengths, displacements, types, &(*population_dt));
	if (createStructRes == MPI_SUCCESS)
		std::cout << "initPopulationStructToMPI - SUCCESS" << std::endl;
	else
	{
		std::cout << "initPopulationStructToMPI - error : " << createStructRes << std::endl;
		abort();
	}
  	MPI_Type_commit(&(*population_dt));
}

void calculateDistancesAsync(int numDP, int numProcs, int id, double** distances, double** demandPoints)
{
    std::cout << id << ": calculateDistances START" << std::endl;

	//Allocate memory and fill 0's
	*distances = (double*) calloc(numDP * numDP, sizeof(double));;

	int procChunkSizeForCalc = numDP / numProcs;
	int offsetForCalc = procChunkSizeForCalc * id;

	#pragma omp parallel
	{
		//Calculation
		#pragma omp for schedule(dynamic)
		for (int iIters = offsetForCalc; iIters < offsetForCalc + procChunkSizeForCalc; iIters++){
			for (int jIters = iIters + 1; jIters < numDP; jIters++){
				double distance = HaversineDistance(demandPoints[iIters], demandPoints[jIters]);
				(*distances)[(iIters * numDP) + jIters] = distance;
			}
		}
	}

	std::cout << id << ": got to calculateDistances barrier" << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
		
	MPI_Allgather(*(distances)+(offsetForCalc * numDP), procChunkSizeForCalc * numDP, 
		MPI_DOUBLE, *distances, procChunkSizeForCalc * numDP, MPI_DOUBLE, MPI_COMM_WORLD);
}
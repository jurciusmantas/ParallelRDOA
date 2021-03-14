#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <sys/time.h>
#include "parallelRDOAlib.h"

double getTime() {
   struct timeval laikas;
   gettimeofday(&laikas, NULL);
   double rez = (double)laikas.tv_sec+(double)laikas.tv_usec/1000000;
   return rez;
}

double HaversineDistance(double* a, double* b) {
   double dlon = fabs(a[0] - b[0]);
   double dlat = fabs(a[1] - b[1]);
   double aa = pow((sin((double)dlon/(double)2*0.01745)),2) + cos(a[0]*0.01745) * cos(b[0]*0.01745) * pow((sin((double)dlat/(double)2*0.01745)),2);
   double c = 2 * atan2(sqrt(aa), sqrt(1-aa));
   double d = 6371 * c; 
   return d;
}

void loadDemandPoints(int numDP, double*** demandPoints) {
    std::cout << "loadDemandPoints START" << std::endl;

	FILE *f;
	f = fopen("../demandPoints.dat", "r");
    *demandPoints = new double*[numDP];
	for (int i = 0; i < numDP; i++) {
		(*demandPoints)[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &(*demandPoints)[i][0], &(*demandPoints)[i][1], &(*demandPoints)[i][2]);
	}
	fclose(f);

    std::cout << "loadDemandPoints END" << std::endl;
}

int rouletteWheel(double* probabilities, int count)
{
    double rndNumber = rand() / (double) RAND_MAX;
    double offset = 0.0;
    int pick = 0;

    for (int i = 0; i < count; i++)
    {
        offset += probabilities[i];
        if (rndNumber < offset)
        {
            pick = i;
            break;
        }
    }

    return pick;
}

int locationAvailable(int location, int numX, int *X, int* bestX)
{
    for (int i = 0; i < numX; i++)
    {
        if (X[i] == location || bestX[i] == location) 
            return 0;
    }

    return 1;
}

#pragma region Generate solution

void generateSolution(int numX, int numCL, int* X, int* bestX, int* ranks, double** distances, int gen_solution_mode)
{
    std::cout << "generateSolution START" << std::endl;

    //New seed on every call
    srand((unsigned)time(0));

    int changed = 0;
    do
    {
        for (int i = 0; i < numX; i ++)
        {
            //Probability for change is 1/s
            //where s is numX
            int probabilityForChange = rand() % numX;
            if (probabilityForChange != 1)
                continue;

            int rankSum = 0, locationProbabilitiesCount = 0, locationChanged = 0;
            for (int j = 0; j < numCL; j++)
            {
                if (locationAvailable(j, numX, X, bestX) == 1)
                    rankSum += ranks[j];
            }

            double* locationProbabilities = new double[numCL];
            for (int j = 0; j < numCL; j++)
            {
                if (locationAvailable(j, numX, X, bestX) == 0)
                    continue;
                
                double locationProbability = 0;
                //RDOA
                if (gen_solution_mode == 0)
                {
                    locationProbabilities[locationProbabilitiesCount] = ranks[j] / (double)rankSum;
                }

                //RDOA-D
                if (gen_solution_mode == 1)
                {
                    double probabilityDenominator = 0;
                    for (int z = 0; z < numCL; z++)
                        probabilityDenominator += ranks[z] / distances[i][j];

                    locationProbabilities[locationProbabilitiesCount] = ranks[j] / (distances[i][j] * probabilityDenominator);
                }

                locationProbabilitiesCount++;
            }

            do
            {
                int pickLocation = rouletteWheel(locationProbabilities, locationProbabilitiesCount);
                if (X[i] != pickLocation)
                {
                    X[i] = pickLocation;
                    locationChanged = 1;
                    changed = 1;
                }
            }
            while (locationChanged == 0);
        }
    }
    while (changed == 0);

    std::cout << "generateSolution END" << std::endl;
}

void randomSolution(int numCL, int numX, int* X) {
	int unique;
	for (int i=0; i<numX; i++) {
		do {
			unique = 1;
			X[i] = (int)((double)rand()/RAND_MAX * numCL);
			for (int j=0; j<i; j++)
				if (X[j] == X[i]) {
					unique = 0;
					break;
				}		
		} while (unique == 0);
	}
}

#pragma endregion

double evaluateSolution(int numX, int numDP, int numPF, int numF, int* X, double** demandPoints, double** distances, int eval_solution_mode) 
{
	double result = 0;
    int bestPF;
    int bestCX;
    double d;

    for (int i = 0; i < numDP; i++) 
    {
        // Nearest from current facility and preexisting
        bestPF = 1e5;
        for (int j = 0; j < numPF; j++) 
        {
            d = distances[i][j];
            if (d < bestPF)
                bestPF = d;
        }

        // Nearest from current facility and solution
        bestCX = 1e5;
        for (int j = 0; j < numX; j++) 
        {
            d = distances[i][X[j]];
            if (d < bestCX) 
                bestCX = d;
        }

        //Binary rule
        if (eval_solution_mode == 0)
        {
            // Attraction is 1/(1 + distance), so smaller distance the better
            if (bestCX < bestPF)
                result += demandPoints[i][2];
            else if (bestCX == bestPF) 
                //result += demandPoints[i][2] * ( 1.0 / (numF + 1)); // Fixed proportion - equal for every firm?
                result += 0.3 * demandPoints[i][2]; // Fixed proportion - 0.3?
        }

        //PartialyBinaryRule
        if (eval_solution_mode == 1)
        {
            double attractionCurrent = 1.0 / (1 + bestCX);
            double attractionPreexisting = 1.0 / (1 + bestPF);

            result += demandPoints[i][2] * (attractionCurrent / (attractionCurrent + (numF * attractionPreexisting)));
        }
    }

	return result;
}

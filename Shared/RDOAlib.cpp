#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <sys/time.h>
#include "RDOAlib.h"

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

void readConfig(int** params, int count)
{
    std::cout << "readConfig START" << std::endl;

    FILE *f;
	f = fopen("../config.dat", "r");

    for (int i = 0; i < count; i++)
        fscanf(f, "%d", &(*params)[i]);

    fclose(f);
    std::cout << "readConfig END" << std::endl;
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

void calculateDistances(int numDP, double*** distances, double** demandPoints)
{
    std::cout << "calculateDistances START" << std::endl;

	//Memory
	*distances = new double*[numDP];
	for (int i = 0; i < numDP; i++) 
		(*distances)[i] = new double[numDP];

	//Calculation
	for (int iIters = 0; iIters < numDP - 1; iIters++){
		(*distances)[iIters][iIters] = 0;

		for (int jIters=iIters + 1; jIters<numDP; jIters++){
			double distance = HaversineDistance(demandPoints[iIters], demandPoints[jIters]);
			(*distances)[iIters][jIters] = distance;
			(*distances)[jIters][iIters] = distance;
		}
	}

	(*distances)[numDP - 1][numDP - 1] = 0;

    std::cout << "calculateDistances END" << std::endl;
}

int rouletteWheel(double* probabilities, int count)
{
    double rndNumber = (double)rand() / RAND_MAX;
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

bool locationAvailable(int location, int numX, int *X, int* bestX)
{
    for (int i = 0; i < numX; i++)
    {
        if (X[i] == location || bestX[i] == location) 
            return false;
    }

    return true;
}

double GetDistance(void* distances, int distancesDim, int i, int j, int numDP)
{
    double result;

    if (distancesDim == 1)
    {
        double* castedDistances = (double*)distances;
        result = castedDistances[i * numDP + j];
        if (result == 0)
            result = castedDistances[j * numDP + i];
    }

    if (distancesDim == 2)
    {
        double** castedDistances = (double**)distances;
        result = castedDistances[i][j];
        if (result == 0)
            result = castedDistances[j][i];
    }

    return result;
}

void updateRanks(int* ranks, int* X, int* bestX, int numCL, int numX, bool success)
{
    bool zeroExists = false;

    for (int i = 0; i < numX; i++)
    {
        if (success)
        {
            ranks[X[i]]++;

            //Search if X contains bestX[i]
            bool contains = false;
            for (int j = 0; j < numX; j++)
            {
                if (X[j] == bestX[i])
                {
                    contains = true;
                    break;
                }
            }

            if (!contains)
            {
                ranks[bestX[i]]--;
                if (ranks[bestX[i]] == 0)
                    zeroExists = true;
            }
        }

        if (!success)
        {
            //Search if bestX contains X[i]
            int contains = 0;
            for (int j = 0; j < numX; j++)
            {
                if (X[i] == bestX[j])
                {
                    contains = true;
                    break;
                }
            }

            if (!contains)
            {
                ranks[X[i]]--;
                if (ranks[X[i]] == 0)
                    zeroExists = true;
            }
        }
    }

    if (zeroExists)
        for(int i = 0; i < numCL; i++)
            ranks[i]++;
}

#pragma region Generate solution

void generateSolutionRDOA(int numX, int numDP, int numCL, int* X, int* bestX, int* ranks, void* distances, int distancesDim)
{
    int changed = 0;
    do
    {
        for (int i = 0; i < numX; i ++)
        {
            //Probability for change is 1/s
            //where s is numX
            double probabilityForChange = rand() / (double) RAND_MAX;
            if (probabilityForChange > (1.0 / numX))
                continue;

            int rankSum = 0;
            for (int j = 0; j < numCL; j++)
            {
                if (locationAvailable(j, numX, X, bestX))
                    rankSum += ranks[j];
            }

            double* locationProbabilities = new double[numCL];
            for (int j = 0; j < numCL; j++)
            {
                if (!locationAvailable(j, numX, X, bestX))
                {
                    locationProbabilities[j] = 0;
                    continue;
                }
                
                locationProbabilities[j] = ranks[j] / (double)rankSum;
            }

            int pickLocation = rouletteWheel(locationProbabilities, numCL);
            X[i] = pickLocation;
            changed = 1;
        }
    }
    while (changed == 0);
}

void generateSolutionRDOAD(int numX, int numDP, int numCL, int* X, int* bestX, int* ranks, void* distances, int distancesDim)
{
    int changed = 0;
    do
    {
        for (int i = 0; i < numX; i ++)
        {
            //Probability for change is 1/s
            //where s is numX
            double probabilityForChange = rand() / (double) RAND_MAX;
            if (probabilityForChange > (1.0 / numX))
                continue;

            double* locationProbabilities = new double[numCL];
            for (int j = 0; j < numCL; j++)
            {
                if (!locationAvailable(j, numX, X, bestX))
                {
                    locationProbabilities[j] = 0;
                    continue;
                }
                
                double probabilityDenominator = 0.0;
                for (int z = 0; z < numCL; z++)
                {
                    if (locationAvailable(z, numX, X, bestX))
                        probabilityDenominator += ranks[z] / GetDistance(distances, distancesDim, z, X[i], numDP);
                }

                locationProbabilities[j] = ranks[j] / (GetDistance(distances, distancesDim, j, X[i], numDP) * probabilityDenominator);
            }

            int pickLocation = rouletteWheel(locationProbabilities, numCL);
            X[i] = pickLocation;
            changed = 1;
        }
    }
    while (changed == 0);
}

void generateSolution(int numX, int numDP, int numCL, int* X, int* bestX, int* ranks, void* distances, int distancesDim, int gen_solution_mode)
{
    if (gen_solution_mode == 0)
        generateSolutionRDOA(numX, numDP, numCL, X, bestX, ranks, distances, distancesDim);
    else if (gen_solution_mode == 1)
        generateSolutionRDOAD(numX, numDP, numCL, X, bestX, ranks, distances, distancesDim);
    else
    {
        std::cout << "generateSolution mode = " << gen_solution_mode << " not supported" << std::endl;
        abort();
    }
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

#pragma region Evaluate solution

double evaluateSolutionBinary(int numX, int numDP, int numPF, int numF, int* X, double** demandPoints, void* distances, int distancesDim)
{
    double result = 0;
    double d;

    for (int i = 0; i < numDP; i++) 
    {
        // Attraction is 1/(1 + distance), so the smaller distance the better

        double bestPF = 1e5;
        // Nearest from current facility from preexisting firms
        for (int j = 0; j < numF * numPF; j++)
        {
            d = GetDistance(distances, distancesDim, i, j, numDP);
            if (d < bestPF)
                bestPF = d;
        }

        double bestCX = 1e5;
        // Nearest from current facility and solution
        for (int j = 0; j < numX; j++) 
        {
            d = GetDistance(distances, distancesDim, i, X[j], numDP);
            if (d < bestCX) 
                bestCX = d;
        }

        if (bestCX < bestPF)
            result += demandPoints[i][2];
        else if (bestCX == bestPF)
            result += 0.3 * demandPoints[i][2]; // Fixed proportion - 0.3?
    }

	return result;
}

double evaluateSolutionPartialyBinary(int numX, int numDP, int numPF, int numF, int* X, double** demandPoints, void* distances, int distancesDim)
{
    double result = 0;
    int* bestPFs;
    double d;

    for (int i = 0; i < numDP; i++) 
    {
        bestPFs = new int[numF];

        // Nearest from current facility and preexisting firms
        for (int preexistingIndex = 0; preexistingIndex < numF; preexistingIndex++)
        {
            double bestPF = 1e5;
            for (int j = 0; j < numPF; j++)
            {
                int secondIndex = preexistingIndex + (j * numF);
                d = GetDistance(distances, distancesDim, i, secondIndex, numDP);
                if (d < bestPF)
                    bestPF = d;
            }
            bestPFs[preexistingIndex] = bestPF;
        }

        // Nearest from current facility and solution
        double bestCX = 1e5;
        for (int j = 0; j < numX; j++) 
        {
            d = GetDistance(distances, distancesDim, i, X[j], numDP);
            if (d < bestCX) 
                bestCX = d;
        }

        /* PartialyBinaryRule */
        double attractionCurrent = 1.0 / (1.0 + bestCX);
        double sumAttractionPreexisting = 0.0;
        for (int j = 0; j < numF; j++)
            sumAttractionPreexisting += 1.0 / (1.0 + bestPFs[j]);

        result += demandPoints[i][2] * (attractionCurrent / (attractionCurrent + sumAttractionPreexisting));

        delete bestPFs;
    }

	return result;
}

double evaluateSolution(int numX, int numDP, int numPF, int numF, int* X, double** demandPoints, void* distances, int distancesDim, int eval_solution_mode) 
{
	if (eval_solution_mode == 0)
        return evaluateSolutionBinary(numX, numDP, numPF, numF, X, demandPoints, distances, distancesDim);
    else if (eval_solution_mode == 1)
        return evaluateSolutionPartialyBinary(numX, numDP, numPF, numF, X, demandPoints, distances, distancesDim);
    else
    {
        std::cout << "evaluateSolution mode = " << eval_solution_mode << " not supported" << std::endl;
        abort();
    }
}

#pragma endregion

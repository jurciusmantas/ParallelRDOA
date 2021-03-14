#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../Shared/parallelRDOAlib.h"

/*  
    GenSolution
    0 - RDOA, 1 - RDOA-D
    EvalSolution
    0 - BinaryRule, 1 - PartialyBinaryRule 
*/
#define GEN_SOLUTION 0
#define EVAL_SOLUTION 0
#define ITERS 10000

using namespace std;

int numDP = 10000;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 5;          // Esanciu objektu skaicius (preexisting facilities)
int numF  = 3;          // Esanciu imoniu skaicius (firms)
int numCL = 25;         // Kandidatu naujiems objektams skaicius (candidate locations)
int numX  = 3;          // Nauju objektu skaicius

double **demandPoints, **distances;
int *X, *bestX, *ranks;

//=============================================================================

double getTime();
void loadDemandPoints();
void calculateDistances();
void generateSolution();
int locationAvailable(int location);
void randomSolution();
double HaversineDistance(double* a, double* b);
double evaluateSolution();
void updateRanks(int success);

//=============================================================================

int main() {
    double ts_start = getTime();
	loadDemandPoints();
	calculateDistances();
	
    X = new int[numX];
	bestX = new int[numX];
	double bestU = -1;
    
    //Init ranks
    ranks = new int[numCL];
    for (int i = 0; i < numCL; i++)
        ranks[i] = i + 1;

    randomSolution();
    double u = evaluateSolution();
    bestU = u;
    for (int i=0; i<numX; i++) 
        bestX[i] = X[i];
	
    for (int iters = 0; iters < ITERS; iters++) {
        printf("iteration - %d \n", iters);
        generateSolution();
        u = evaluateSolution();

        if (u > bestU) 
        {
            updateRanks(1);
            bestU = u;

            /*  Note:
                Shouldn't we generate the solution BEFORE assigning X = X'?
                If we do it affter assigning, the locations that were in X but weren't in X'
                become available to pick.
            */
            for (int i=0; i<numX; i++) 
                bestX[i] = X[i];
        }
        else
            updateRanks(0);
    }

    // Write results
    ofstream resultsFile;
    resultsFile.open("results.txt", ios_base::app);
	for (int i=0; i<numX; i++) 
        resultsFile << bestX[i] << " ";

	resultsFile << "(" << bestU << "), " << getTime() - ts_start << endl;
}

#pragma region Demand points

void loadDemandPoints() {
    printf("loadDemandPoints START\n");

	FILE *f;
	f = fopen("../demandPoints.dat", "r");
	demandPoints = new double*[numDP];
	for (int i=0; i<numDP; i++) {
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);

    printf("loadDemandPoints END\n");
}

void calculateDistances(){
    printf("calculateDistances START\n");

	//Memory
	distances = new double*[numDP];
	for (int i = 0; i < numDP; i++) 
		distances[i] = new double[numDP];

	//Calculation
	for (int iIters = 0; iIters < numDP - 1; iIters++){
		distances[iIters][iIters] = 0;

		for (int jIters=iIters + 1; jIters<numDP; jIters++){
			double distance = HaversineDistance(demandPoints[iIters], demandPoints[jIters]);
			distances[iIters][jIters] = distance;
			distances[jIters][iIters] = distance;
		}
	}

	distances[numDP - 1][numDP - 1] = 0;

    printf("calculateDistances END\n");
}

double HaversineDistance(double* a, double* b) {
   double dlon = fabs(a[0] - b[0]);
   double dlat = fabs(a[1] - b[1]);
   double aa = pow((sin((double)dlon/(double)2*0.01745)),2) + cos(a[0]*0.01745) * cos(b[0]*0.01745) * pow((sin((double)dlat/(double)2*0.01745)),2);
   double c = 2 * atan2(sqrt(aa), sqrt(1-aa));
   double d = 6371 * c; 
   return d;
}

#pragma endregion

#pragma region Generate solution

void randomSolution() {
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

void generateSolution()
{
    printf("generateSolution START\n");

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
                if (locationAvailable(j) == 1)
                    rankSum += ranks[j];
            }

            double* locationProbabilities = new double[numCL];
            for (int j = 0; j < numCL; j++)
            {
                if (locationAvailable(j) == 0)
                    continue;
                
                double locationProbability = 0;
                //RDOA
                if (GEN_SOLUTION == 0)
                {
                    locationProbabilities[locationProbabilitiesCount] = ranks[j] / (double)rankSum;
                }

                //RDOA-D
                if (GEN_SOLUTION == 1)
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

    printf("generateSolution END\n");
}

int locationAvailable(int location)
{
    for (int i = 0; i < numX; i++)
    {
        if (X[i] == location || bestX[i] == location) 
            return 0;
    }

    return 1;
}

#pragma endregion

double evaluateSolution() 
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
        if (EVAL_SOLUTION == 0)
        {
            // Attraction is 1/(1 + distance), so smaller distance the better
            if (bestCX < bestPF)
                result += demandPoints[i][2];
            else if (bestCX == bestPF) 
                //result += demandPoints[i][2] * ( 1.0 / (numF + 1)); // Fixed proportion - equal for every firm?
                result += 0.3 * demandPoints[i][2]; // Fixed proportion - 0.3?
        }

        //PartialyBinaryRule
        if (EVAL_SOLUTION == 1)
        {
            double attractionCurrent = 1.0 / (1 + bestCX);
            double attractionPreexisting = 1.0 / (1 + bestPF);

            result += demandPoints[i][2] * (attractionCurrent / (attractionCurrent + (numF * attractionPreexisting)));
        }
    }

	return result;
}

void updateRanks(int success)
{
    int zeroExists = 0;

    for (int i = 0; i < numX; i++)
    {
        if (success == 1)
        {
            ranks[X[i]]++;

            //Search if X contains bestX[i]
            int contains = 0;
            for (int j = 0; j < numX; j++)
            {
                if (X[j] == bestX[i])
                {
                    contains = 1;
                    break;
                }
            }

            if (contains == 0)
            {
                ranks[bestX[i]]--;
                if (ranks[bestX[i]] == 0)
                    zeroExists = 1;
            }
        }

        if (success == 0)
        {
            //Search if bestX contains X[i]
            int contains = 0;
            for (int j = 0; j < numX; j++)
            {
                if (X[i] == bestX[j])
                {
                    contains = 1;
                    break;
                }
            }

            if (contains == 0)
            {
                ranks[X[i]]--;
                if (ranks[X[i]] == 0)
                    zeroExists = 1;
            }
        }
    }

    if (zeroExists)
        for(int i = 0; i < numCL; i++)
            ranks[i]++;
}

double getTime() {
   struct timeval laikas;
   gettimeofday(&laikas, NULL);
   double rez = (double)laikas.tv_sec+(double)laikas.tv_usec/1000000;
   return rez;
}

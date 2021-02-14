#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

/*  
    GenSolution
    0 - RDOA, 1 - RDOA-D
    EvalSolution
    0 - BinaryRule, 1 - PartialyBinaryRule 
*/
#DEFINE GEN_SOLUTION = 0;
#DEFINE EVAL_SOLUTION = 0;

using namespace std;

int numDP = 10000;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 5;          // Esanciu objektu skaicius (preexisting facilities)
int numF  = 3;          // Esanciu imoniu skaicius (firms)
int numCL = 50;         // Kandidatu naujiems objektams skaicius (candidate locations)
int numX  = 3;          // Nauju objektu skaicius

double **demandPoints, **distances;
int *X, *bestX, *ranks;
long rankSum;

//=============================================================================

double getTime();
void loadDemandPoints();
void calculateDistances();
void generateSolution();
int locationAvailable(int location);
void randomSolution(int *X);
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X);
void updateRanks();

//=============================================================================

int main() {
	loadDemandPoints();
	calculateDistances();
	
    X = new int[numX];
	bestX = new int[numX];
	double bestU = -1;
    
    //Init ranks and rankSum
    rankSum = 0;
    ranks = new int[numDP];
    for (int i = 0; i < numDP; i++)
    {
        ranks[i] = i + 1;
        rankSum += ranks[i];
    }

    randomSolution(X);
    double u = evaluateSolution(X);
    bestU = u;
    bestX = X;
	
    for (int iters = 0; iters < 10000; iters++) {
        generateSolution();
        u = evaluateSolution();

        if (u > bestU) 
        {
            bestU = u;
            updateRanks();

            /*  Note:
                Shouldn't we generate the solution BEFORE assigning X = X'?
                If we do it affter assigning, the locations that were in X but weren't in X'
                become available to pick.
             */
            for (int i=0; i<numX; i++) 
                bestX[i] = X[i];
        }
    }

	
	double tf = getTime();     // Skaiciavimu pabaigos laikas

	cout << "Geriausias sprendinys: ";
	for (int i=0; i<numX; i++) cout << bestX[i] << " ";
	cout << "(" << bestU << ")" << endl << "Skaiciavimo trukme: " << tf-ts_start << endl;
}

#pragma region Demand points

void loadDemandPoints() {
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double*[numDP];
	for (int i=0; i<numDP; i++) {
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

void calculateDistances(){
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

void randomSolution(int *X) {
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
    //New seed on every call (?)
    srand((unsigned)time(0));

    int changed = 0;
    do
    {
        for (int i = 0; i < numX; i ++)
        {
            //Probability for change is 1/s
            //where s is numX
            int probability = rand() % numX;
            if (probability != 1)
                continue;

            /*  Notes to ask:
                1)  Is the probability to choose the location correct? The article states:
                    'a single facility will be changed in average'
                    What if none facilities will be changed? (--> thats why do/while, but is it optimal?)
                
                2)  If codeflow manages to come here (past the if above, do we MUST change this location?
                    There is a possibility that the below cycle will not pick any location.
                    So, do/while is neccassary here too?
            */
            for (int j = 0; j < numDP; j++) 
            {
                if (locationAvailable(j) == 0)
                    continue;

                double locationProbability = 0;
                //RDOA
                if (GEN_SOLUTION == 0)
                    locationProbability = ranks[j] / rankSum;

                //RDOA-D
                if (GEN_SOLUTION == 1)
                {
                    double probabilityDenominator = 0;
                    for (int z = 0; z < numDP; z++)
                        probabilityDenominator += ranks[z] / distances[i][j];

                    locationProbability = ranks[j] / (distances[i][j] * probabilityDenominator);
                }

                if (rand() % 100 < (locationProbability * 100))
                {
                    X[i] = j;
                    changed = 1;
                    break;
                }
            }
        }
    }
    while (changed == 0);
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

    for (int i = 0; i < numDP; i++) 
    {
        int bestPF;
        int bestX;
        double d;
        bestPF = 1e5;

        // Nearest from current facility and preexisting
        for (int j = 0; j < numPF; j++) 
        {
            d = distances[i][j];
            if (d < bestPF)
                bestPF = d;
        }

        // Nearest from current facility and solution
        bestX = 1e5;
        for (int j = 0; j < numX; j++) 
        {
            d = distances[i][X[j]];
            if (d < bestX) 
                bestX = d;
        }

        //Binary rule
        if (EVAL_SOLUTION == 0)
        {
            // Attraction is 1/(1 + distance), so smaller distance the better
            if (bestX < bestPF)
                result += demandPoints[i][2];
            else if (bestX == bestPF) 
                result += demandPoints[i][2] * ( 1 / (numF + 1)); // Fixed proportion - equal for every firm?
        }

        //PartialyBinaryRule
        if (EVAL_SOLUTION == 1)
        {
            double attractionCurrent = 1 / (1 + bestX);
            double attractionPreexisting = 1 / (1 + bestPF);

            result += demandPoints[i][2] * (attractionCurrent / (attractionCurrent + (numF * attractionPreexisting));
        }

    }

	return result;
}

double getTime() {
   struct timeval laikas;
   gettimeofday(&laikas, NULL);
   double rez = (double)laikas.tv_sec+(double)laikas.tv_usec/1000000;
   return rez;
}

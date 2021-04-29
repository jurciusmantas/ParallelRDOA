/* Population struct */
typedef struct pop_item
{
    double solution;
    int* locations;
}
populationItem;

/* Initialize population */
void initPopulation(populationItem** population, int populationSize, int numX);

/* Search for solution in population */
populationItem search(populationItem* population, int populationSize, int* X, int numX);

/* Insert (if neccessary) */
void insert(populationItem* population, int* locations, int numX, double solution, int* itemsInPopulation, int populationSize, int* popRanks = NULL);

double* popItemToArray(populationItem item, int numX);
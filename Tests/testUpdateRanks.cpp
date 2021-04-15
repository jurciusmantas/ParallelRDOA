#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "../Shared/RDOAlib.h"

using namespace std;

bool testBetter();
bool testNotBetter();

int main() {
    cout << "Update ranks tests: " << endl;

    bool testBetterResult = testBetter();
    if (!testBetterResult)
        cout << "testBetter - FALSE" << endl;

    bool testNotBetterResut = testNotBetter();
    if (!testNotBetterResut)
        cout << "testNotBetter - FALSE" << endl;

    system("pause");
}

bool testBetter(){
    int ranks[5] = { 1, 2, 3, 4, 5 };
    int x1[] = { 1, 2 };
    int bestX1[] = { 2, 3 };

    updateRanks(ranks, x1, bestX1, 5, 2, true);
    if (ranks[0] == 1 && //same 
        ranks[1] == 3 && //++ 
        ranks[2] == 4 && //++ 
        ranks[3] == 3 && //--
        ranks[4] == 5)   //same
        cout << "testBetter - 1 - OK" << endl;
    else
    {
        cout << "testBetter - 1 - FALSE" << endl;
        return false;
    }

    for (int i = 0; i < 5; i++)
        ranks[i] = i;

    int x2[] = { 0, 1 };
    int bestX2[] = { 3, 4 };

    updateRanks(ranks, x2, bestX2, 5, 2, true);
    if (ranks[0] == 1 && //++ 
        ranks[1] == 2 && //++ 
        ranks[2] == 2 && //same 
        ranks[3] == 2 && //--
        ranks[4] == 3)   //--
        cout << "testBetter - 2 - OK" << endl;
    else
    {
        cout << "testBetter - 2 - FALSE" << endl;
        return false;
    }

    for (int i = 0; i < 5; i++)
        ranks[i] = i;

    int x3[] = { 0, 4 };
    int bestX3[] = { 1, 3 };

    // Test zeroExists
    updateRanks(ranks, x3, bestX3, 5, 2, true);
    if (ranks[0] == 2 && //++ ++
        ranks[1] == 1 && //-- ++
        ranks[2] == 3 && //   ++ 
        ranks[3] == 3 && //-- ++
        ranks[4] == 6)   //++ ++
        cout << "testBetter - 3 - OK" << endl;
    else
    {
        cout << "testBetter - 3 - FALSE" << endl;
        return false;
    }

    return true;
}

bool testNotBetter(){
    int ranks[5] = { 1, 2, 3, 4, 5 };
    int x1[] = { 2, 3 };
    int bestX1[] = { 1, 2 };

    updateRanks(ranks, x1, bestX1, 5, 2, false);
    if (ranks[0] == 1 && //same 
        ranks[1] == 2 && //same 
        ranks[2] == 3 && //same
        ranks[3] == 3 && //--
        ranks[4] == 5)   //same
        cout << "testNotBetter - 1 - OK" << endl;
    else
    {
        cout << "testNotBetter - 1 - FALSE" << endl;
        return false;
    }

    for (int i = 0; i < 5; i++)
        ranks[i] = i + 1;

    int x2[] = { 3, 4 };
    int bestX2[] = { 0, 1 };

    updateRanks(ranks, x2, bestX2, 5, 2, false);
    if (ranks[0] == 1 && //same 
        ranks[1] == 2 && //same 
        ranks[2] == 3 && //same 
        ranks[3] == 3 && //--
        ranks[4] == 4)   //--
        cout << "testNotBetter - 2 - OK" << endl;
    else
    {
        cout << "testNotBetter - 2 - FALSE" << endl;
        return false;
    }

    for (int i = 0; i < 5; i++)
        ranks[i] = i + 1;

    int x3[] = { 0, 4 };
    int bestX3[] = { 1, 3 };

    // Test zeroExists
    updateRanks(ranks, x3, bestX3, 5, 2, false);
    if (ranks[0] == 1 && //-- ++
        ranks[1] == 3 && //   ++
        ranks[2] == 4 && //   ++ 
        ranks[3] == 5 && //   ++
        ranks[4] == 5)   //-- ++
        cout << "testNotBetter - 3 - OK" << endl;
    else
    {
        cout << "testNotBetter - 3 - FALSE" << endl;
        return false;
    }

    return true;
}
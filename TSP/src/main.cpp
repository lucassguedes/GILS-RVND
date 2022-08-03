#include "readData.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdio>
#include <iomanip>
#include <ctime>
#include <queue>
#include <exception>
#include "utils.hpp"

#define _2_OPT_MIN_SIZE 3
#define REINSERTION_SIZE 1
#define _OR_OPT2_SIZE 2
#define _OR_OPT3_SIZE 3

#define MAX_ALPHA_PERCENT 26

/**
 * @file main.cpp
 *
 * @brief Implementation of GILS-RVND metaheuristic to solve TSP.
 *
 * @author Lucas Guedes
 * Contact: lucassilva@eng.ci.ufpb.br
 */

#define ORIGIN_VERTEX 1
 
double **matrizAdj;
int dimension;

double globabBestCost;

enum Neighborhood
{
    SWAP,
    _2_OPT,
    REINSERTION,
    OR_OPT2,
    OR_OPT3,
}; //Values to be used in RVND

std::vector<Neighborhood> g_NL = {SWAP, _2_OPT, REINSERTION, OR_OPT2, OR_OPT3};

/*Begin: Functions prototypes*/
void doubleBridge(std::vector<int> &s, int i, int j);
void reinsertion(std::vector<int>&s, const int i, const int j, const int subRouteSize);
void findBestNeighbor(const Neighborhood neighborhood, std::vector<int> &s, double &bestCost);
void initializeCandidateList(std::vector<int> &candidateList);
void RVND(std::vector<int> &s, double &currentCost);
std::vector<int> Perturb(std::vector<int> s);
std::vector<int> construction(const double alpha);
std::vector<int> GILS_RVND(const int Imax, const int Iils, std::vector<double> R);
/*End: Functions prototypes*/

int main(int argc, char **argv)
{
    readData(argc, argv, &dimension, &matrizAdj);

    srand(time(NULL));
    double currentCost;


    clock_t beginTime = clock();

    /*Begin: Setting parameters*/
    const int Imax = 50;
    int Iils = (dimension >= 150) ? dimension / 2 : dimension;
    std::vector<double> R(MAX_ALPHA_PERCENT);

    for(int i = 0; i < MAX_ALPHA_PERCENT; i++)
        R[i] = (double)i/100;

    /*End: Setting parameters*/

    std::vector<int> s = GILS_RVND(Imax, Iils, R);


    clock_t endTime = clock();

    std::cout << "TIME: " << (double)(endTime - beginTime) / CLOCKS_PER_SEC << std::endl;
    std::cout << "COST: " << getSolutionCost(matrizAdj, s) << std::endl;

    /*Checking the solution*/
    // if (isValidSolution(s, dimension))
    // {
    //     std::cout << "A solução encontrada é VÁLIDA!\n";
    // }
    // else
    // {
    //     std::cout << "A solução encontrada é INVÁLIDA!\n";
    // }

    // showVector("Final solution: ", s);/*uncomment this to see the final solution*/

    return 0;
}

std::vector<int> GILS_RVND(const int Imax, const int Iils, std::vector<double> R){
    std::vector<int>final_s,s,s1;
    double alpha, cost, cost1;
    int iterILS;
    globabBestCost = INT_MAX;

    srand(time(NULL));

    for(size_t i = 0; i < Imax; i++){
        alpha = R[rand() % R.size()];
        s = construction(alpha);

        s1 = s;

        cost = cost1 = getSolutionCost(matrizAdj, s1);

        iterILS = 0;
        while(iterILS < Iils){
        
            RVND(s, cost);
            if(cost < cost1){
                s1 = s;
                cost1 = cost;
                iterILS = 0;
            }
            s = Perturb(s1);

            cost = getSolutionCost(matrizAdj, s);

            
            iterILS++;

        }   

        if(cost1 < globabBestCost){
            final_s = s1;
            globabBestCost = cost1;
        }

    
    }

    return final_s;

}

std::vector<int> construction(const double alpha)
{
    std::vector<int> s = {ORIGIN_VERTEX};
    std::vector<int> candidateList;
    int sizeOfRCL, random_index;

    initializeCandidateList(candidateList);

    int root = ORIGIN_VERTEX;

    while (!candidateList.empty())
    {
        /*Sorting the candidate list according to their distances to current root.*/
        std::sort(candidateList.begin(), candidateList.end(), [root, s](int a, int b)
                  { return matrizAdj[root][a] < matrizAdj[root][b]; });

        sizeOfRCL = alpha * candidateList.size();

        if (sizeOfRCL == 0)
        {
            random_index = 0;
        }
        else
        {
            random_index = rand() % sizeOfRCL;
        }

        root = candidateList[random_index];
        s.push_back(root);
        candidateList.erase(candidateList.begin() + random_index);
    }
    s.push_back(ORIGIN_VERTEX);

    return s;
}

void RVND(std::vector<int> &s, double &currentCost)
{
    std::vector<int> s1 = s;
    std::vector<Neighborhood> NL = g_NL;
    Neighborhood randomNeigboorhood;
    double rvndCost = currentCost;
    int randomIndex;

    while (!NL.empty())
    {
        if (NL.size() == 1)
            randomIndex = 0;
        else
            randomIndex = rand() % NL.size();

        randomNeigboorhood = NL[randomIndex];

        findBestNeighbor(randomNeigboorhood, s1, rvndCost);

        if (rvndCost < currentCost)
        {
            s = s1;
            currentCost = rvndCost;
            NL = g_NL;
        }
        else
        {
            NL.erase(NL.begin() + randomIndex);
        }
    }
}

std::vector<int> Perturb(std::vector<int> s){
    int i,j;

    i = rand() % (s.size() - 4) + 1; 
    j = rand() % (s.size() - 4) + 1;

    if(i == j)
    {
        if(i + 2 <= s.size() - 1 - 2){
            j += 2;
        }else{  
            j -= 2;
        }
    }


    doubleBridge(s, i,j);

    return s;
}

void doubleBridge(std::vector<int> &s, int i, int j){
  int amount = 2;
  
  std::swap_ranges(s.begin() + i, s.begin() + i + amount, s.begin() + j);
}

void findBestNeighbor(const Neighborhood neighborhood, std::vector<int> &s, double &bestCost)
{
    double delta, deltaAval;
    double reusableCosts[2];
    int choosedIndex1, choosedIndex2;
    size_t subRouteSize;

    delta = 0;

    size_t routeSize = s.size();

    switch (neighborhood)
    {
    case SWAP:

        for (int i = 1; i < routeSize - 2; i++)
        {
            reusableCosts[0] = -matrizAdj[s[i - 1]][s[i]];
            reusableCosts[1] = -matrizAdj[s[i - 1]][s[i]] - matrizAdj[s[i]][s[i + 1]];
            for (int j = i + 1; j < routeSize - 1; j++)
            {
                if (abs(j - i) == 1) //are neighbors
                {
                    deltaAval = matrizAdj[s[i]][s[j + 1]] 
                                + matrizAdj[s[i - 1]][s[j]] 
                                - matrizAdj[s[j]][s[j + 1]] 
                                + reusableCosts[0];
                }
                else
                {
                    deltaAval = matrizAdj[s[i - 1]][s[j]] 
                                + matrizAdj[s[j]][s[i + 1]] 
                                + matrizAdj[s[j - 1]][s[i]] 
                                + matrizAdj[s[i]][s[j + 1]] 
                                - matrizAdj[s[j - 1]][s[j]] 
                                - matrizAdj[s[j]][s[j + 1]] 
                                + reusableCosts[1];
                }

                if (deltaAval < delta)
                {
                    choosedIndex1 = i;
                    choosedIndex2 = j;
                    delta = deltaAval;
                }
            }
        }

        if (delta < 0)
        {
            std::swap(s[choosedIndex1], s[choosedIndex2]);
            bestCost += delta;
        }

        break;
    case _2_OPT:

        for(int i = 1; i < routeSize - _2_OPT_MIN_SIZE - 1; i++){
            reusableCosts[0] = - matrizAdj[s[i]][s[i-1]];
            for(int j = i+_2_OPT_MIN_SIZE; j < routeSize - 1; j++){
                
                deltaAval = matrizAdj[s[j]][s[i-1]] + matrizAdj[s[i]][s[j+1]] - matrizAdj[s[j]][s[j+1]] + reusableCosts[0];

                if(deltaAval < delta){
                    choosedIndex1 = i;
                    choosedIndex2 = j;
                    delta = deltaAval;
                }

            }
        }

        if(delta < 0){
            reverse(s.begin()+choosedIndex1,s.begin()+choosedIndex2+1); 
            bestCost += delta; 
        }

        break;
    case REINSERTION:
        subRouteSize = REINSERTION_SIZE;
        for(int i = 1; i < routeSize - subRouteSize - 1; i++){

            reusableCosts[0] = + matrizAdj[s[i-1]][s[i+1]]
                               - matrizAdj[s[i-1]][s[i]]
                               - matrizAdj[s[i]][s[i+1]];

            for(int j = 1; j < routeSize - 1 ; j++){

                if(j >= i + subRouteSize + 1 || j < i){

                    deltaAval = matrizAdj[s[j-1]][s[i]] 
                                + matrizAdj[s[i]][s[j]] 
                                - matrizAdj[s[j-1]][s[j]]
                                + reusableCosts[0];


                    if(deltaAval < delta){
                        choosedIndex1 = i;
                        choosedIndex2 = j;
                        delta = deltaAval;
                    }

                }
            
                
            }
        }

        if(delta < 0){
            reinsertion(s,choosedIndex1, choosedIndex2, subRouteSize);
            bestCost += delta;
        }
    
        break;

    case OR_OPT2:
        subRouteSize = _OR_OPT2_SIZE;
        for(int i = 1; i < routeSize - subRouteSize - 1; i++){

            reusableCosts[0] = + matrizAdj[s[i-1]][s[i+subRouteSize]]
                               - matrizAdj[s[i-1]][s[i]]
                               - matrizAdj[s[i + subRouteSize - 1]][s[i+subRouteSize]];

            for(int j = 1; j < routeSize - 1 ; j++){

                if(j >= i + subRouteSize + 1 || j < i){

                    deltaAval = matrizAdj[s[j-1]][s[i]] 
                                + matrizAdj[s[i + subRouteSize - 1]][s[j]] 
                                - matrizAdj[s[j-1]][s[j]]
                                + reusableCosts[0];


                    if(deltaAval < delta){
                        choosedIndex1 = i;
                        choosedIndex2 = j;
                        delta = deltaAval;
                    }

                }
            
                
            }
        }

        if(delta < 0){
            reinsertion(s,choosedIndex1, choosedIndex2, subRouteSize);
            bestCost += delta;
        }


        break;

    case OR_OPT3:
        subRouteSize = _OR_OPT3_SIZE;
        for(int i = 1; i < routeSize - subRouteSize - 1; i++){

            reusableCosts[0] = + matrizAdj[s[i-1]][s[i+subRouteSize]]
                               - matrizAdj[s[i-1]][s[i]]
                               - matrizAdj[s[i + subRouteSize - 1]][s[i+subRouteSize]];

            for(int j = 1; j < routeSize - 1 ; j++){

                if(j >= i + subRouteSize + 1 || j < i){

                    deltaAval = matrizAdj[s[j-1]][s[i]] 
                                + matrizAdj[s[i + subRouteSize - 1]][s[j]] 
                                - matrizAdj[s[j-1]][s[j]]
                                + reusableCosts[0];


                    if(deltaAval < delta){
                        choosedIndex1 = i;
                        choosedIndex2 = j;
                        delta = deltaAval;
                    }

                }
            
                
            }
        }

        if(delta < 0){
            reinsertion(s,choosedIndex1, choosedIndex2, subRouteSize);
            bestCost += delta;
        }
        break;
    }
}

void reinsertion(std::vector<int>&s, int i, int j, const int subRouteSize){

    std::vector<int>subRoute(subRouteSize);
    std::copy(s.begin()+i,s.begin()+i+subRouteSize,subRoute.begin());
    s.erase(s.begin()+i,s.begin()+i+subRouteSize);
    if(j > i){
        j -= (subRouteSize);
    }
    s.insert(s.begin()+j, subRoute.begin(), subRoute.end());

}

void initializeCandidateList(std::vector<int> &candidateList)
{
    for (int i = 2; i <= dimension; i++)
    {
        candidateList.push_back(i);
    }
}
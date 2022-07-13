#include<iostream>
#include "../include/readData.h"
#include <vector>
#include <algorithm>

double **matrix;
int dimension;

double global_best_cost;

#define ORIGIN_VERTEX 1

enum Neighborhood
{
    SWAP,
    _2_OPT,
    REINSERTION,
    OR_OPT2,
    OR_OPT3,
};//Values to be used in RVND

std::vector<Neighborhood> g_NL = {SWAP, _2_OPT, REINSERTION, OR_OPT2, OR_OPT3};

void initialize_candidate_list(std::vector<int> &candidate_list);
std::vector<int> construction(const double alpha);

int main(int argc, char **argv){
    std::vector<int> s;

    readData(argc, argv, &dimension, &matrix);

    s = construction(0.20);


    for(auto k : s)
    {
        std::cout << k << " ";
    }
    std::cout << std::endl;

    return 0;
}

std::vector<int> construction(const double alpha)
{
    std::vector<int> s = {ORIGIN_VERTEX};
    std::vector<int> candidate_list;
    int size_of_rcl, random_index;

    initialize_candidate_list(candidate_list);

    int root = ORIGIN_VERTEX;

    while(!candidate_list.empty())
    {
        /*Sorting the candidate list according to their distances to current root*/
        std::sort(candidate_list.begin(), candidate_list.end(), [root, s](int a, int b){return matrix[root][a] < matrix[root][b];});

        size_of_rcl = alpha * candidate_list.size();
        random_index = (size_of_rcl == 0) ? 0 : rand() % size_of_rcl;

        root = candidate_list[random_index];
        s.push_back(root);
        candidate_list.erase(candidate_list.begin() + random_index);
    }
    s.push_back(ORIGIN_VERTEX);

    return s;
}

void initialize_candidate_list(std::vector<int> &candidate_list)
{
    for(size_t i = 2; i <= dimension; i++)
    {
        candidate_list.push_back(i);
    }
}
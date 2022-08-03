#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <iostream>
#include <vector>
#include <algorithm>

#define ORIGIN_VERTEX 1

void showVector(std::string msg, std::vector<int> v)
{
    std::cout << msg << std::endl;
    for (auto k : v)
        std::cout << k << " ";
    std::cout << std::endl;
}

bool isValidSolution(std::vector<int> s, const int dimension)
{
    size_t solutionSize = s.size();

    if (s[0] != ORIGIN_VERTEX || s[solutionSize - 1] != ORIGIN_VERTEX)
        return false;

    for (int i = 1; i < dimension; i++)
    {
        if (std::find(s.begin(), s.end(), i) == s.end())
            return false;
    }

    return true;
}

double getSolutionCost(double **matrix, std::vector<int> s)
{
    size_t solSize = s.size();
    double totalCost = 0;
    for (int i = 0; i < solSize - 1; i++)
    {
        totalCost += matrix[s[i]][s[i + 1]];
    }

    return totalCost;
}

#endif
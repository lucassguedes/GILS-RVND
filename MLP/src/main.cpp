#include<iostream>
#include "../include/readData.h"
#include <vector>
#include <algorithm>

double **matrix;
int dimension;

double global_best_cost;

#define ORIGIN_VERTEX 1

struct SubSeqInfo{
    int begin_vertex;
    int end_vertex;

    int w;
    double t;
    double c;
};

class ReoptData
{
    public: 
    std::vector< std::vector<int> > w;
    std::vector< std::vector<double> > t;
    std::vector< std::vector<double> > c;

    ReoptData(int dimension)
    {
        this->create_reopt_data(dimension);
    }

    ReoptData(std::vector<int>s,int dimension)
    {
        this->create_reopt_data(dimension);
        this->fill_reopt_data(s,dimension);
    }

    void create_reopt_data(int dimension)
    {
        this->w = std::vector<std::vector<int> > (dimension+1,std::vector<int>(dimension+1,0));
        this->t = std::vector<std::vector<double> > (dimension+1,std::vector<double>(dimension+1,0));
        this->c = std::vector<std::vector<double> > (dimension+1,std::vector<double>(dimension+1,0));
    }

    void fill_reopt_data(std::vector<int>s, int dimension)
    {
        for(int i = 0; i <= dimension; i++)
            this->w[i][i] = (i > 0);

        for(int i = 0; i < dimension; i++)
        {
            for(int j = i+1; j <= dimension; j++)
            {
                this->w[i][j] = this->w[i][j-1] + this->w[j][j];
                this->t[i][j] = this->t[i][j-1] + matrix[s[j-1]][s[j]] + this->t[j][j];
                this->c[i][j] = this->c[i][j-1] + this->w[j][j]*(this->t[i][j-1] + matrix[s[j-1]][s[j]]) + this->c[j][j];
            }
        }

        for(int i = dimension; i > 1 ; i--)
        {
            for(int j = i-1; j >= 0 ; j--)
            {
                this->w[i][j] = this->w[i][j+1] + this->w[j][j];
                this->t[i][j] = this->t[i][j+1] + matrix[s[j+1]][s[j]] + this->t[j][j];
                this->c[i][j] = this->c[i][j+1] + this->w[j][j]*(this->t[i][j+1] + matrix[s[j+1]][s[j]]) + this->c[j][j];
            }
        }
    }
};

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

double get_total_cost(std::vector<int>s);
double get_latency(std::vector<int>s);


void test_cost_c(std::vector<int>s, ReoptData reopt)
{
    std::vector<int>s1;
    size_t seqsize = s.size();
    double cost1, cost2;
    for(size_t i = 0; i < seqsize; i++)
    {
        for(size_t j = i+1; j < seqsize; j++)
        {   
            s1.clear();
            s1.resize((j+1) - i);
            std::copy(s.begin()+i, s.begin()+j+1, s1.begin());

            std::reverse(s1.begin(), s1.end());
            cost1 = get_total_cost(s1); 
            // cost1 = get_latency(s1);
            cost2 = reopt.c[j][i];
            // cost2 = reopt.t[i][j];

            if(cost1 != cost2)
                std::cout << "DIFERENTES - i = " << i << ", j = " << j <<std::endl;   
            else 
                std::cout << "IGUAIS - i = " << i << ", j = " << j <<std::endl;
        }
    }
}

int main(int argc, char **argv){
    std::vector<int> s;

    readData(argc, argv, &dimension, &matrix);

    clock_t begin_time = clock();

    s = construction(0.20);

    ReoptData reopt = ReoptData(s,dimension);

    test_cost_c(s,reopt);


    clock_t end_time = clock();

    std::cout << "COST: " << get_total_cost(s) << std::endl;
    std::cout << "TIME: " << (double)(end_time - begin_time) / CLOCKS_PER_SEC << std::endl;

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


double get_total_cost(std::vector<int>s)
{
    double total_cost = 0;
    double latency = 0;
    size_t seqsize = s.size();

    for(size_t i = 0; i < seqsize - 1; i++)
    {
        latency += matrix[s[i]][s[i+1]];
        total_cost += latency;
    }

    return total_cost;
}

double get_latency(std::vector<int>s)
{
    double latency = 0;
    size_t seqsize = s.size();

    for(size_t i = 0; i < seqsize - 1; i++)
    {
        latency += matrix[s[i]][s[i+1]];
    }

    return latency;
}
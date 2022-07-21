#include<iostream>
#include "../include/readData.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <ctime>

double **matrix;
int dimension;

double global_best_cost;

#define ORIGIN_VERTEX 1

class Candidate{
    public:
    int vertex;
    double distance_to_origin;

    Candidate(int vertex, double distance_to_origin){
        this->vertex = vertex;
        this->distance_to_origin = distance_to_origin;
    }
};

class SubSeqInfo{
    public:
    int begin_vertex;
    int end_vertex;

    int w;
    double t;
    double c;

    SubSeqInfo(int begin_vertex, int end_vertex, int w, double t, double c){
        this->begin_vertex = begin_vertex;
        this->end_vertex = end_vertex;

        this->w = w;
        this->t = t;
        this->c = c;
    }

    SubSeqInfo(){}
};

SubSeqInfo merge_subseq(std::vector<int> &s, SubSeqInfo &ss1, SubSeqInfo &ss2)
{
    SubSeqInfo result;

    result.begin_vertex = ss1.begin_vertex;
    result.end_vertex = ss2.end_vertex;

    result.w = ss1.w + ss2.w;
    result.t = ss1.t + matrix[s[ss1.end_vertex]][s[ss2.begin_vertex]] + ss2.t;
    result.c = ss1.c + ss2.w * (ss1.t + matrix[s[ss1.end_vertex]][s[ss2.begin_vertex]]) + ss2.c;

    return result;
}

class ReoptData
{
    public: 
    std::vector< std::vector<int> > w;
    std::vector< std::vector<double> > t;
    std::vector< std::vector<double> > c;

    ReoptData()  
    { 
        this->create_reopt_data();
    }

    ReoptData(std::vector<int>&s)
    {
        this->create_reopt_data();
        this->fill_reopt_data(s);
    }

    void create_reopt_data()
    {
        this->w = std::vector<std::vector<int> > (dimension+1,std::vector<int>(dimension+1,0));
        this->t = std::vector<std::vector<double> > (dimension+1,std::vector<double>(dimension+1,0));
        this->c = std::vector<std::vector<double> > (dimension+1,std::vector<double>(dimension+1,0));
    }

    void fill_reopt_data(std::vector<int>s)
    {
        for(int i = 0; i <= dimension; i++)
            this->w[i][i] = (i > 0);

        for(int i = 0, k = dimension; i < dimension; i++, k--)
        {
            for(int j = i+1, u = k - 1; j <= dimension; j++, u--)
            {
                this->w[i][j] = this->w[i][j-1] + this->w[j][j];
                this->t[i][j] = this->t[i][j-1] + matrix[s[j-1]][s[j]] + this->t[j][j];
                this->c[i][j] = this->c[i][j-1] + this->w[j][j]*(this->t[i][j-1] + matrix[s[j-1]][s[j]]) + this->c[j][j];

                this->w[k][u] = this->w[k][u+1] + this->w[u][u];
                this->t[k][u] = this->t[k][u+1] + matrix[s[u+1]][s[u]] + this->t[u][u];
                this->c[k][u] = this->c[k][u+1] + this->w[u][u]*(this->t[k][u+1] + matrix[s[u+1]][s[u]]) + this->c[u][u];
            }
        }

        // for(int i = dimension; i > 1 ; i--)
        // {
        //     for(int j = i-1; j >= 0 ; j--)
        //     {
        //         this->w[i][j] = this->w[i][j+1] + this->w[j][j];
        //         this->t[i][j] = this->t[i][j+1] + matrix[s[j+1]][s[j]] + this->t[j][j];
        //         this->c[i][j] = this->c[i][j+1] + this->w[j][j]*(this->t[i][j+1] + matrix[s[j+1]][s[j]]) + this->c[j][j];
        //     }
        // }
    }

    void update_reopt_data(std::vector<int> &s, int best_i, int best_j)
    {
        for(int i = best_i; i < dimension; i++)
        {
            for(int j = 0; j < best_j; j++)
            {
                this->w[i][j] = this->w[i][j-1] + this->w[j][j];
                this->t[i][j] = this->t[i][j-1] + matrix[s[j-1]][s[j]] + this->t[j][j];
                this->c[i][j] = this->c[i][j-1] + this->w[j][j]*(this->t[i][j-1] + matrix[s[j-1]][s[j]]) + this->c[j][j];
            }
        }

        for(int i = 0;  i < best_j; i++)
        {
            for(int j = best_i; j <= dimension; j++)
            {
                this->w[i][j] = this->w[i][j-1] + this->w[j][j];
                this->t[i][j] = this->t[i][j-1] + matrix[s[j-1]][s[j]] + this->t[j][j];
                this->c[i][j] = this->c[i][j-1] + this->w[j][j]*(this->t[i][j-1] + matrix[s[j-1]][s[j]]) + this->c[j][j];
            }
        }

    }
};

#define _2_OPT_MIN_SIZE 3
#define REINSERTION_SIZE 1
#define _OR_OPT2_SIZE 2
#define _OR_OPT3_SIZE 3

#define MAX_ALPHA_PERCENT 26

enum Neighborhood
{
    SWAP,
    _2_OPT,
    REINSERTION,
    OR_OPT2,
    OR_OPT3,
};//Values to be used in RVND

// std::vector<Neighborhood> g_NL = {SWAP, _2_OPT, REINSERTION, OR_OPT2, OR_OPT3};

ReoptData reopt;

void initialize_candidate_list(std::vector<Candidate> &candidate_list, int root);
std::vector<int> construction(const double alpha);
void RVND(std::vector<int>&s, double &current_cost);
void find_best_neighbor(const Neighborhood neighborhood, std::vector<int>&s, double &rvnd_cost, int &best_i, int &best_j);
void reinsertion(std::vector<int>&s, int i, int j, const int subroute_size);
std::vector<int> Perturb(std::vector<int> s);
void double_bridge(std::vector<int>&s, int i, int j);
std::vector<int> GILS_RVND(const int Imax, const int Iils, std::vector<double>R);

double get_total_cost(std::vector<int>s);
double get_latency(std::vector<int>s);

int main(int argc, char **argv){
    std::vector<int> s;

    readData(argc, argv, &dimension, &matrix);

    clock_t begin_time = clock();


    const int Imax = 50;
    int Iils = (dimension >= 150) ? dimension / 2 : dimension;
    std::vector<double> R(MAX_ALPHA_PERCENT);

    for(int i = 0; i < MAX_ALPHA_PERCENT; i++)
        R[i] = (double)i/100;


    s = GILS_RVND(Imax, Iils, R);


    // s = construction(0.20);

    // double current_cost = get_total_cost(s);

    // RVND(s, current_cost); 

    // reopt = ReoptData(s, dimension);

    // test_or_opt2(s);

    clock_t end_time = clock();
    for(auto k : s)
        std::cout << k  << " ";
    std::cout << std::endl;
    std::cout << "COST: " << get_total_cost(s) << std::endl;
    std::cout << "TIME: " << (double)(end_time - begin_time) / CLOCKS_PER_SEC << std::endl;

    return 0;
}

std::vector<int> GILS_RVND(const int Imax, const int Iils, std::vector<double>R)
{
    std::vector<int> final_s, s, s1;
    double alpha, cost, cost1;
    int iterILS;
    global_best_cost = std::numeric_limits<double>::max();

    srand(time(NULL));

    for(size_t i = 0; i < Imax; i++)
    {
        alpha = R[rand() % R.size()];
        s = construction(alpha);

        s1 = s;

        cost = cost1 = get_total_cost(s1);


        iterILS = 0;
        while(iterILS < Iils)
        {
            RVND(s, cost);

            if(cost < cost1)
            {
                s1 = s;
                cost1 = cost;
                iterILS = 0;
            }

            s = Perturb(s1);

            cost = get_total_cost(s);

            iterILS++;
        }

        if(cost1 < global_best_cost)
        {
            final_s = s1;
            global_best_cost = cost1;
        }
    }

    return final_s;
}

void RVND(std::vector<int>&s, double &current_cost)
{
    std::vector<int> s1 = s;
    std::vector<Neighborhood> NL = {SWAP, _2_OPT, REINSERTION, OR_OPT2, OR_OPT3};
    Neighborhood random_neigborhood;
    int random_index;
    int best_i, best_j;

    reopt = ReoptData(s1); /*Initialize re-optimization data structures*/
    
    double rvnd_cost = current_cost = reopt.c[0][dimension];

    while(!NL.empty())
    {
        if(NL.size() == 1)
            random_index = 0;
        else
            random_index = rand() % NL.size();

        // random_neigborhood = NL[random_index];

        find_best_neighbor(NL[random_index], s1, rvnd_cost, best_i, best_j);

        if(rvnd_cost < current_cost)
        {
            s = s1;
            current_cost = rvnd_cost;
            NL = {SWAP, _2_OPT, REINSERTION, OR_OPT2, OR_OPT3};
            // reopt.update_reopt_data(s, best_i, best_j);
            reopt = ReoptData(s);
        }else{
            s1 = s;
            rvnd_cost = current_cost;
            NL.erase(NL.begin() + random_index);
        }
    }
}

void find_best_neighbor(const Neighborhood neighborhood, std::vector<int>&s, double &rvnd_cost, int &best_i, int &best_j)
{
    double cost, best_cost_found;
    int choosed_index_1, choosed_index_2;
    size_t subroute_size;

    SubSeqInfo ss1,ss2,ss3,ss4,ss5,ss6,ssr0,ssr1,ssr2,ssr3;

    best_cost_found = rvnd_cost;

    size_t route_size = s.size();

    switch(neighborhood)
    {
        case SWAP:
            for(int i = 1; i < route_size - 2; i++)
            {
                ssr0 = SubSeqInfo(0,i-1,reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]);
                ssr1 = SubSeqInfo(i, i, reopt.w[i][i], reopt.t[i][i], reopt.c[i][i]);
                for(int j = i + 1; j < route_size - 1; j++)
                {
                    if(abs(j - i) == 1) /*are neighbors*/
                    {
                        /*merge here*/
                        ss3 = SubSeqInfo(j,j,reopt.w[j][j], reopt.t[j][j], reopt.c[j][j]);
                        ss4 = SubSeqInfo(j+1, route_size - 1, reopt.w[j+1][route_size-1], reopt.t[j+1][route_size-1], reopt.c[j+1][route_size-1]);
                        ss1 = merge_subseq(s, ssr0, ss3);
                        ss2 = merge_subseq(s, ssr1, ss4);

                        cost = merge_subseq(s, ss1, ss2).c;
                    }else{
                        /*merge here*/
                        ss4 = SubSeqInfo(j,j,reopt.w[j][j], reopt.t[j][j], reopt.c[j][j]);
                        ss5 = SubSeqInfo(i+1, j-1, reopt.w[i+1][j-1], reopt.t[i+1][j-1], reopt.c[i+1][j-1]);
                        ss6 = SubSeqInfo(j+1, route_size-1, reopt.w[j+1][route_size-1], reopt.t[j+1][route_size-1], reopt.c[j+1][route_size-1]);
                        ss1 = merge_subseq(s, ssr0, ss4);
                        ss2 = merge_subseq(s, ss5, ssr1);
                        ss3 = merge_subseq(s, ss1, ss2);

                        cost = merge_subseq(s, ss3, ss6).c;
                    }

                    if(cost < best_cost_found)
                    {
                        choosed_index_1 = i;
                        choosed_index_2 = j;
                        best_cost_found = cost;
                    }

                    // std::vector<int>s1 = s;
                    // std::swap(s1[i], s1[j]);

                    // double f = get_total_cost(s1);

                    // if(f != cost)
                    //     std::cout << "(SWAP) - Diferentes!\n";
                }
            }

            if(best_cost_found < rvnd_cost)
            {
                std::swap(s[choosed_index_1], s[choosed_index_2]);
                rvnd_cost = best_cost_found;
            }
        break;
        case _2_OPT:
            for(int i = 1; i < route_size - _2_OPT_MIN_SIZE - 1; i++)
            {
                ssr0 = SubSeqInfo(0,i-1, reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]);
                for(int j = i+_2_OPT_MIN_SIZE; j < route_size - 1; j++)
                {
                    ss2 = SubSeqInfo(j,i, reopt.w[j][i], reopt.t[j][i], reopt.c[j][i]);
                    ss3 = SubSeqInfo(j+1,route_size-1, reopt.w[j+1][route_size-1], reopt.t[j+1][route_size-1], reopt.c[j+1][route_size-1]);
                    ss1 = merge_subseq(s, ssr0, ss2);
                    cost = merge_subseq(s, ss1, ss3).c;

                    if(cost < best_cost_found)
                    {
                        choosed_index_1 = i;
                        choosed_index_2 = j;
                        best_cost_found = cost;
                    }

                    // std::vector<int>s1 = s;
                    // reverse(s1.begin()+i, s1.begin()+j+1);

                    // double f = get_total_cost(s1);

                    // if(f != cost)
                    //     std::cout << "(2-OPT) - Diferentes!\n";
                }
            }

            if(best_cost_found < rvnd_cost)
            {
                reverse(s.begin()+choosed_index_1, s.begin()+choosed_index_2+1);
                rvnd_cost = best_cost_found;
            }
        break;
        case REINSERTION:
            subroute_size = REINSERTION_SIZE;
            for(int i = 1; i < route_size - subroute_size - 1; i++)
            {
                ssr0 = SubSeqInfo(0, i - 1, reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]);
                ssr1 = SubSeqInfo(i, i+subroute_size-1, reopt.w[i][i+subroute_size-1], reopt.t[i][i+subroute_size-1], reopt.c[i][i+subroute_size-1]);
                ssr2 = SubSeqInfo(i+subroute_size, route_size-1, reopt.w[i+subroute_size][route_size-1], reopt.t[i+subroute_size][route_size-1], reopt.c[i+subroute_size][route_size-1]);
                for(int j = 1; j < route_size - 1; j++)
                {
                    if(j >= i + subroute_size + 1){
                        ss3 = SubSeqInfo(i+subroute_size, j-1, reopt.w[i+subroute_size][j-1], reopt.t[i+subroute_size][j-1], reopt.c[i+subroute_size][j-1]);
                        ss4 = SubSeqInfo(j, route_size-1, reopt.w[j][route_size-1], reopt.t[j][route_size-1], reopt.c[j][route_size-1]);
                        ss1 = merge_subseq(s, ssr0, ss3);
                        ss2 = merge_subseq(s, ssr1, ss4);

                        cost = merge_subseq(s, ss1, ss2).c;
                    }else if(j < i){
                        /*merge here*/
                        ss3 = SubSeqInfo(0, j-1, reopt.w[0][j-1], reopt.t[0][j-1], reopt.c[0][j-1]);
                        ss4 = SubSeqInfo(j, i-1, reopt.w[j][i-1], reopt.t[j][i-1], reopt.c[j][i-1]);
                        ss1 = merge_subseq(s, ss3, ssr1);
                        ss2 = merge_subseq(s, ss4, ssr2);

                        cost = merge_subseq(s, ss1, ss2).c;
                    }else{
                        continue;
                    }

                    if(cost < best_cost_found)
                    {
                        choosed_index_1 = i;
                        choosed_index_2 = j;
                        best_cost_found = cost;
                    }

                    // std::vector<int>s1 = s;
                    // reinsertion(s1, i, j, subroute_size);

                    // double f = get_total_cost(s1);

                    // if(f != cost)
                    //     std::cout << "(REINSERTION) - Diferentes!\n";
                }
            }

            if(best_cost_found < rvnd_cost)
            {
                reinsertion(s, choosed_index_1, choosed_index_2, subroute_size);
                rvnd_cost = best_cost_found;
            }
        break;
        case OR_OPT2:
            subroute_size = _OR_OPT2_SIZE;
            for(int i = 1; i < route_size - subroute_size - 1; i++)
            {
                ssr0 = SubSeqInfo(0, i - 1, reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]);
                ssr1 = SubSeqInfo(i, i+subroute_size-1, reopt.w[i][i+subroute_size-1], reopt.t[i][i+subroute_size-1], reopt.c[i][i+subroute_size-1]);
                ssr2 = SubSeqInfo(i+subroute_size, route_size-1, reopt.w[i+subroute_size][route_size-1], reopt.t[i+subroute_size][route_size-1], reopt.c[i+subroute_size][route_size-1]);
                for(int j = 1; j < route_size - 1; j++)
                {
                    if(j >= i + subroute_size + 1){
                        ss3 = SubSeqInfo(i+subroute_size, j-1, reopt.w[i+subroute_size][j-1], reopt.t[i+subroute_size][j-1], reopt.c[i+subroute_size][j-1]);
                        ss4 = SubSeqInfo(j, route_size-1, reopt.w[j][route_size-1], reopt.t[j][route_size-1], reopt.c[j][route_size-1]);
                        ss1 = merge_subseq(s, ssr0, ss3);
                        ss2 = merge_subseq(s, ssr1, ss4);

                        cost = merge_subseq(s, ss1, ss2).c;
                    }else if(j < i){
                        /*merge here*/
                        ss3 = SubSeqInfo(0, j-1, reopt.w[0][j-1], reopt.t[0][j-1], reopt.c[0][j-1]);
                        ss4 = SubSeqInfo(j, i-1, reopt.w[j][i-1], reopt.t[j][i-1], reopt.c[j][i-1]);
                        ss1 = merge_subseq(s, ss3, ssr1);
                        ss2 = merge_subseq(s, ss4, ssr2);

                        cost = merge_subseq(s, ss1, ss2).c;
                    }else{
                        continue;
                    }

                    if(cost < best_cost_found)
                    {
                        choosed_index_1 = i;
                        choosed_index_2 = j;
                        best_cost_found = cost;
                    }

                    // std::vector<int>s1 = s;
                    // reinsertion(s1, i, j, subroute_size);

                    // double f = get_total_cost(s1);

                    // if(f != cost)
                    //     std::cout << "(REINSERTION) - Diferentes!\n";
                }
            }

            if(best_cost_found < rvnd_cost)
            {
                reinsertion(s, choosed_index_1, choosed_index_2, subroute_size);
                rvnd_cost = best_cost_found;
            }
        break;
        case OR_OPT3:
            subroute_size = _OR_OPT3_SIZE;
            for(int i = 1; i < route_size - subroute_size - 1; i++)
            {
                ssr0 = SubSeqInfo(0, i - 1, reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]);
                ssr1 = SubSeqInfo(i, i+subroute_size-1, reopt.w[i][i+subroute_size-1], reopt.t[i][i+subroute_size-1], reopt.c[i][i+subroute_size-1]);
                ssr2 = SubSeqInfo(i+subroute_size, route_size-1, reopt.w[i+subroute_size][route_size-1], reopt.t[i+subroute_size][route_size-1], reopt.c[i+subroute_size][route_size-1]);
                for(int j = 1; j < route_size - 1; j++)
                {
                    if(j >= i + subroute_size + 1){
                        ss3 = SubSeqInfo(i+subroute_size, j-1, reopt.w[i+subroute_size][j-1], reopt.t[i+subroute_size][j-1], reopt.c[i+subroute_size][j-1]);
                        ss4 = SubSeqInfo(j, route_size-1, reopt.w[j][route_size-1], reopt.t[j][route_size-1], reopt.c[j][route_size-1]);
                        ss1 = merge_subseq(s, ssr0, ss3);
                        ss2 = merge_subseq(s, ssr1, ss4);

                        cost = merge_subseq(s, ss1, ss2).c;
                    }else if(j < i){
                        /*merge here*/
                        ss3 = SubSeqInfo(0, j-1, reopt.w[0][j-1], reopt.t[0][j-1], reopt.c[0][j-1]);
                        ss4 = SubSeqInfo(j, i-1, reopt.w[j][i-1], reopt.t[j][i-1], reopt.c[j][i-1]);
                        ss1 = merge_subseq(s, ss3, ssr1);
                        ss2 = merge_subseq(s, ss4, ssr2);

                        cost = merge_subseq(s, ss1, ss2).c;
                    }else{
                        continue;
                    }

                    if(cost < best_cost_found)
                    {
                        choosed_index_1 = i;
                        choosed_index_2 = j;
                        best_cost_found = cost;
                    }

                    // std::vector<int>s1 = s;
                    // reinsertion(s1, i, j, subroute_size);

                    // double f = get_total_cost(s1);

                    // if(f != cost)
                    //     std::cout << "(REINSERTION) - Diferentes!\n";
                }
            }

            if(best_cost_found < rvnd_cost)
            {
                reinsertion(s, choosed_index_1, choosed_index_2, subroute_size);
                rvnd_cost = best_cost_found;
            }
        break;
    }
    
    best_i = choosed_index_1;
    best_j = choosed_index_2;
}

std::vector<int> Perturb(std::vector<int> s)
{
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

    double_bridge(s, i, j);

    return s;
}

void double_bridge(std::vector<int>&s, int i, int j)
{
    int n = 2;

    std::swap_ranges(s.begin() + i, s.begin() + i + n, s.begin() + j);
}

void reinsertion(std::vector<int>&s, int i, int j, const int subroute_size)
{
    std::vector<int>subroute(subroute_size);
    std::copy(s.begin()+i, s.begin()+i+subroute_size, subroute.begin());
    s.erase(s.begin()+i, s.begin()+i+subroute_size);

    if(j > i)
    {
        j -= (subroute_size);
    }

    s.insert(s.begin()+j, subroute.begin(), subroute.end());
}

bool by_distance_to_root(Candidate a, Candidate b)
{
    return a.distance_to_origin < b.distance_to_origin;
}


std::vector<int> construction(const double alpha)
{
    std::vector<int> s = {ORIGIN_VERTEX};
    std::vector<Candidate> candidate_list;
    int size_of_rcl, random_index;

    int root = ORIGIN_VERTEX;

    initialize_candidate_list(candidate_list, root);

    while(!candidate_list.empty())
    {
        /*Sorting the candidate list according to their distances to current root*/
        std::sort(candidate_list.begin(), candidate_list.end(), by_distance_to_root);

        size_of_rcl = alpha * candidate_list.size();
        random_index = (size_of_rcl == 0) ? 0 : rand() % size_of_rcl;

        root = candidate_list[random_index].vertex;
        s.push_back(root);
        candidate_list.erase(candidate_list.begin() + random_index);

        for(size_t i = 0; i < candidate_list.size(); i++)
        {
            candidate_list[i].distance_to_origin = matrix[candidate_list[i].vertex][root];
        }
    }
    s.push_back(ORIGIN_VERTEX);

    return s;
}

void initialize_candidate_list(std::vector<Candidate> &candidate_list, int root)
{
    for(size_t i = 2; i <= dimension; i++)
    {
        candidate_list.push_back(Candidate(i, matrix[i][root]));
        
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
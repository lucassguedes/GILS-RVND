#include "../include/utils.h"


// void test_or_opt2(std::vector<int>s)
// {
//     size_t route_size = s.size();
//     size_t subroute_size = _OR_OPT2_SIZE;

//     double cost1, cost2;

//     for(int i = 1; i < route_size - subroute_size - 1; i++)
//     {
//         for(int j = 1; j < route_size - 1; j++)
//         {
//             if(j >= i + subroute_size + 1){
//                 SubSeqInfo ss1 = merge_subseq(s, SubSeqInfo(0, i - 1, reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]), SubSeqInfo(i+subroute_size, j-1, reopt.w[i+subroute_size][j-1], reopt.t[i+subroute_size][j-1], reopt.c[i+subroute_size][j-1]));
//                 SubSeqInfo ss2 = merge_subseq(s, SubSeqInfo(i, i+subroute_size-1, reopt.w[i][i+subroute_size-1], reopt.t[i][i+subroute_size-1], reopt.c[i][i+subroute_size-1]), SubSeqInfo(j, route_size-1, reopt.w[j][route_size-1], reopt.t[j][route_size-1], reopt.c[j][route_size-1]));

//                 cost1 = merge_subseq(s, ss1, ss2).c;
//             }else if(j < i){
//                 /*merge here*/
//                 SubSeqInfo ss1 = merge_subseq(s, SubSeqInfo(0, j-1, reopt.w[0][j-1], reopt.t[0][j-1], reopt.c[0][j-1]), SubSeqInfo(i, i+subroute_size-1, reopt.w[i][i+subroute_size-1], reopt.t[i][i+subroute_size-1], reopt.c[i][i+subroute_size-1]));
//                 SubSeqInfo ss2 = merge_subseq(s, SubSeqInfo(j, i-1, reopt.w[j][i-1], reopt.t[j][i-1], reopt.c[j][i-1]), SubSeqInfo(i+subroute_size, route_size-1, reopt.w[i+subroute_size][route_size-1], reopt.t[i+subroute_size][route_size-1], reopt.c[i+subroute_size][route_size-1]));

//                 cost1 = merge_subseq(s, ss1, ss2).c;
//             }else{
//                 continue;
//             }

//             std::vector<int>s1 = s;

//             reinsertion(s1,i, j, subroute_size);

//             cost2 = get_total_cost(s1);

//             if(cost1 != cost2)
//             {
//                 std::cout << "[OR-OPT2] DIFERENTES - i = " << i << ", j = " << j <<std::endl;   
//             }
//             else 
//             {
//                 std::cout << "[OR-OPT2] IGUAIS - i = " << i << ", j = " << j <<std::endl;
//             }
//         }
//     }
// }

// void test_reinsertion(std::vector<int>s)
// {
//     size_t route_size = s.size();
//     size_t subroute_size = REINSERTION_SIZE;

//     double cost1, cost2;

//     for(int i = 1; i < route_size - subroute_size - 1; i++)
//     {
//         for(int j = 1; j < route_size - 1; j++)
//         {
//             if(j >= i + subroute_size + 1){
//                 SubSeqInfo ss1 = merge_subseq(s, SubSeqInfo(0, i - 1, reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]), SubSeqInfo(i+subroute_size, j-1, reopt.w[i+subroute_size][j-1], reopt.t[i+subroute_size][j-1], reopt.c[i+subroute_size][j-1]));
//                 SubSeqInfo ss2 = merge_subseq(s, SubSeqInfo(i, i+subroute_size-1, reopt.w[i][i+subroute_size-1], reopt.t[i][i+subroute_size-1], reopt.c[i][i+subroute_size-1]), SubSeqInfo(j, route_size-1, reopt.w[j][route_size-1], reopt.t[j][route_size-1], reopt.c[j][route_size-1]));

//                 cost1 = merge_subseq(s, ss1, ss2).c;
//             }else if(j < i){
//                 /*merge here*/
//                 SubSeqInfo ss1 = merge_subseq(s, SubSeqInfo(0, j-1, reopt.w[0][j-1], reopt.t[0][j-1], reopt.c[0][j-1]), SubSeqInfo(i, i+subroute_size-1, reopt.w[i][i+subroute_size-1], reopt.t[i][i+subroute_size-1], reopt.c[i][i+subroute_size-1]));
//                 SubSeqInfo ss2 = merge_subseq(s, SubSeqInfo(j, i-1, reopt.w[j][i-1], reopt.t[j][i-1], reopt.c[j][i-1]), SubSeqInfo(i+subroute_size, route_size-1, reopt.w[i+subroute_size][route_size-1], reopt.t[i+subroute_size][route_size-1], reopt.c[i+subroute_size][route_size-1]));

//                 cost1 = merge_subseq(s, ss1, ss2).c;
//             }else{
//                 continue;
//             }

//             std::vector<int>s1 = s;

//             reinsertion(s1,i, j, subroute_size);

//             cost2 = get_total_cost(s1);

//             if(cost1 != cost2)
//             {
//                 std::cout << "[Reinsertion] DIFERENTES - i = " << i << ", j = " << j <<std::endl;   
//             }
//             else 
//             {
//                 std::cout << "[Reinsertion] IGUAIS - i = " << i << ", j = " << j <<std::endl;
//             }
//         }
//     }
// }


// void test_2opt(std::vector<int>s)
// {
//     double cost1, cost2;
//     size_t route_size = s.size();
//     for(int i = 1; i < route_size - _2_OPT_MIN_SIZE - 1; i++)
//     {
//         for(int j = i+_2_OPT_MIN_SIZE; j < route_size - 1; j++)
//         {
//             SubSeqInfo ss1 = merge_subseq(s, SubSeqInfo(0,i-1, reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]), SubSeqInfo(j,i, reopt.w[j][i], reopt.t[j][i], reopt.c[j][i]));
//             cost1 = merge_subseq(s, ss1, SubSeqInfo(j+1,route_size-1, reopt.w[j+1][route_size-1], reopt.t[j+1][route_size-1], reopt.c[j+1][route_size-1])).c;

//             std::vector<int>s1 = s;

//             reverse(s1.begin()+i,s1.begin()+j+1); 

//             cost2 = get_total_cost(s1);

//             if(cost1 != cost2)
//             {
//                 std::cout << "[2-OPT] DIFERENTES - i = " << i << ", j = " << j <<std::endl;   
//             }
//             else 
//             {
//                 std::cout << "[2-OPT] IGUAIS - i = " << i << ", j = " << j <<std::endl;
//             }

//         }
//     }
// }


// void test_swap(std::vector<int>s)
// {
//     int route_size = s.size();
//     double cost1, cost2;
//     for(int i = 1; i < route_size - 2; i++)
//     {
//         for(int j = i + 1; j < route_size - 1; j++)
//         {
//             if(abs(j - i) == 1) /*are neighbors*/
//             {
//                 /*merge here*/
//                 SubSeqInfo ss1 = merge_subseq(s, SubSeqInfo(0,i-1,reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]), SubSeqInfo(j,j,reopt.w[j][j], reopt.t[j][j], reopt.c[j][j]));
//                 SubSeqInfo ss2 = merge_subseq(s, SubSeqInfo(i, i, reopt.w[i][i], reopt.t[i][i], reopt.c[i][i]), SubSeqInfo(j+1, route_size - 1, reopt.w[j+1][route_size-1], reopt.t[j+1][route_size-1], reopt.c[j+1][route_size-1]));

//                 cost1 = merge_subseq(s, ss1, ss2).c;
//             }else{
//                 /*merge here*/
//                 SubSeqInfo ss1 = merge_subseq(s, SubSeqInfo(0, i -1, reopt.w[0][i-1], reopt.t[0][i-1], reopt.c[0][i-1]), SubSeqInfo(j,j,reopt.w[j][j], reopt.t[j][j], reopt.c[j][j]));
//                 SubSeqInfo ss2 = merge_subseq(s, SubSeqInfo(i+1, j-1, reopt.w[i+1][j-1], reopt.t[i+1][j-1], reopt.c[i+1][j-1]), SubSeqInfo(i, i, reopt.w[i][i], reopt.t[i][i], reopt.c[i][i]));
//                 SubSeqInfo ss3 = merge_subseq(s, ss1, ss2);

//                 cost1 = merge_subseq(s, ss3, SubSeqInfo(j+1, route_size-1, reopt.w[j+1][route_size-1], reopt.t[j+1][route_size-1], reopt.c[j+1][route_size-1])).c;
//             }

//             std::vector<int> s1 = s;

//             std::swap(s1[i], s1[j]);

//             cost2 = get_total_cost(s1);

//             if(cost1 != cost2)
//             {
//                 std::cout << "DIFERENTES - i = " << i << ", j = " << j <<std::endl;   
//             }
//             else 
//             {
//                 std::cout << "IGUAIS - i = " << i << ", j = " << j <<std::endl;
//             }
//         }
//     }

// }



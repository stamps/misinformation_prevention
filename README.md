# Reverse prevention sampling for misinformation mitigation in social networks
Code for the paper "Reverse prevention sampling for misinformation mitigation in social networks" in ICDT 2020.

Information:
--------------------------------------------------------
Version 1.0: Implementation of RPS Algorithm for Misinformation Prevention. For more details about RPS, please read our paper: "Michael Simpson, Venkatesh Srinivasan, Alex Thomo, Reverse prevention sampling for misinformation mitigation in social networks, ICDT 2020"

Contact Authors: Michael Simpson (michaelesimp@gmail.com)


Requirements:
--------------------------------------------------------
In order to compile all the tools, it requires GCC 4.7.2 and later (which also support OpenMP 3.1).


Compile:
--------------------------------------------------------
Use `make' command to compile everything

How to use:
--------------------------------------------------------
 g++ -std=c++11 LimitMisinformation_tang_parallel.cpp -fopenmp -O3 sfmt/SFMT.c -o LimitMisinformation_tang_parallel

./LimitMisinformation_tang_parallel wordassociation-2011 k topk weighted

# limit_misinformation
code for misinformation prevention

 g++ -std=c++11 LimitMisinformation_tang_parallel.cpp -fopenmp -O3 sfmt/SFMT.c -o LimitMisinformation_tang_parallel

./LimitMisinformation_tang_parallel wordassociation-2011 k topk weighted

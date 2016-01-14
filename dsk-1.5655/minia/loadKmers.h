
#ifndef LOADKMERS_H
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <algorithm> //for sort
#include <math.h>
#include <unordered_map>
#define LOADKMERS_H

using namespace std;

extern int *Kmerlist; //list of all K values
extern int smallestKmer; // smallest value of kmer asked for computation. largest value is equal to sizeKmer in the functions
extern int totalKmers; //number of k values for computation
extern unordered_map<int,int> kmerlength_map; // to take care of begining of k-mer lengths
int* loadKmers(char *kmerfname);
int reestimate_partitions(int size,uint64_t partition_volume, double * lkmer_counts, long * hash_vals, int * partition_files,int fact);
typedef pair<int,double> clmer;
bool comparator(const clmer& l,const clmer& r);
#endif

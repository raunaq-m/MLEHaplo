#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <inttypes.h>
#include <sys/stat.h> // for S_IRWXU etc
#include <vector>
#include <sys/statvfs.h> // to determine available disk space
#include <dirent.h> // to clear the temp directory
#include <cmath>
#include <unordered_map> // to map small kmer counts which decides the partitions into which k-mers are separated

#ifndef SORTINGCOUNT_H
#define SORTINGCOUNT_H

#include "Bank.h"
#include "Kmer.h"
#include "Utils.h"
#include "OAHash.h"
#include "loadKmers.h"

using namespace std;

typedef uint32_t uint_abundance_t;
void sorting_count(Bank *Sequences, char *prefix, int max_memory, int max_disk_space, bool write_count, int verbose, bool reverse);

#endif

//
//  Bank.cpp
//
//  Created by Guillaume Rizk on 28/11/11.
//

//TEST

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#define __STDC_FORMAT_MACROS    // Added by Raunaq

#include <algorithm>
#include <iostream>
#include <sys/stat.h>
#include <inttypes.h>
       #include <stdio.h>
       #include <stdlib.h>
       #include <string.h>
#include <cmath> // for log2f

#include "Bank.h"
#include "Kmer.h" // Bank (almost) doesn't need Kmer.h, but KmersBuffer certainly does
#include "lut.h"
#include <errno.h>
using namespace std;

off_t fsize(const char *filename) {
    struct stat st; 
    
    if (stat(filename, &st) == 0)
        return st.st_size;
    
    return -1; 
}

// just a macro to open file indexed by i
void Bank::open_stream(int i)
{
    buffered_file[i]->stream = gzopen(buffered_file[i]->fname,"r");
    if (buffered_file[i]->stream == NULL)
    {
        printf("error opening file: %s\n",buffered_file[i]->fname);
        exit(1);
    }
}
// and close it
void Bank::close_stream(int i)
{
    gzclose(buffered_file[i]->stream);
    buffered_file[i]->stream = NULL;
}

// the following functions are adapted from kseq.h by Heng Li (https://github.com/attractivechaos/klib)
inline bool rebuffer(buffered_file_t *bf)
{
    if (bf->eof)
        return false;
    bf->buffer_start = 0;
    bf->buffer_end = gzread(bf->stream, bf->buffer, BUFFER_SIZE);
    if (bf->buffer_end < BUFFER_SIZE)
        bf->eof = 1;
    if (bf->buffer_end == 0) 
        return false;
    return true;
}

inline signed char buffered_getc(buffered_file_t *bf)
{
    if (bf->buffer_start >= bf->buffer_end)
        if (! rebuffer(bf))
            return -1;
    return (signed char) ( bf->buffer[bf->buffer_start++] );
}

#define nearest_power_of_2(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

inline signed int Bank::buffered_gets(buffered_file_t *bf, variable_string_t *s, char *dret, bool append, bool allow_spaces)
{
    if (dret) *dret = 0;
    if (!append)
        s->length = 0;
    if (bf->buffer_start >= bf->buffer_end && bf->eof)
        return -1;
    while (1)
    {
        int i;
        if (bf->buffer_start >= bf->buffer_end)
            if (! rebuffer(bf))
                break;
        if (allow_spaces)
        {
            for (i = bf->buffer_start; i < bf->buffer_end ; i++)
                if (bf->buffer[i] == '\n') 
                    break;
        }
        else
        {
            for (i = bf->buffer_start; i < bf->buffer_end ; i++)
                // isspace() answers yes for ' ', \t, \n, \v, \f, \r
                if (isspace(bf->buffer[i]))
                    break;
        }
        if (s->max - s->length < (i - bf->buffer_start + 1))
        {
            s->max = s->length + (i - bf->buffer_start + 1);
            nearest_power_of_2(s->max);
            s->string = (char*)realloc(s->string,s->max);
        } 
        memcpy(s->string + s->length, bf->buffer + bf->buffer_start, i - bf->buffer_start);
        s->length += i - bf->buffer_start;
        bf->buffer_start = i + 1;
        if (i < bf->buffer_end)
        {
            if (dret)
                *dret = bf->buffer[i];
            break;
        }
    }
    if (s->string == NULL)
    {
        s->max = 256;
        s->string = (char*)calloc(256,1);
    }
    else if ( allow_spaces && s->length > 1 && s->string[s->length-1] == '\r')
        s->length--;
    s->string[s->length]= '\0';
    return s->length;
}

void Bank::rewind_all()
{
    for (int i=0; i<nb_files; i++)
    {
        if (buffered_file[i]->stream != NULL)
        {
            gzclose(buffered_file[i]->stream);
            buffered_file[i]->stream = NULL;
        }
        buffered_file[i]->last_char = buffered_file[i]->eof = buffered_file[i]->buffer_start = buffered_file[i]->buffer_end = 0;
    }
    index_file = 0;
    open_stream(index_file);
}

// THIS READS FASTQ or FASTA, compressed with gzip or not
// no limit on read length, allows multi-line reads
// returns true if a read was successfuly read
//         false if end of file
// adapted from kseq.h by Heng Li (https://github.com/attractivechaos/klib)
bool  Bank::get_next_seq_from_file(char **nseq, char **cheader, int *len, int *hlen, int file_id)
{
    signed char c;
    buffered_file_t *bf = buffered_file[file_id];
    if (bf->last_char == 0)
    {
        while ( (c = buffered_getc(bf)) != -1 && c != '>' && c != '@'); // go to next header
        if (c == -1)
            return false; // eof
        bf->last_char = c;
    }
    read->length = dummy->length = 0;

    if (buffered_gets(bf, header, (char *)&c, false, false) < 0) //ici
        return false; // eof
    if (c != '\n')
        buffered_gets(bf, dummy, NULL, true, true); // read header //dummy instead of header to stop before first space
    
    if (read->string == NULL)
    {
        read->max = 256;
        read->string = (char*) malloc(read->max);
    }
    while ( (c = buffered_getc(bf)) != -1 && c != '>' && c != '+' && c != '@')
    {
        if (c == '\n')
            continue; // empty line
        read->string[read->length++] = c;
        buffered_gets(bf, read, NULL, true, true);
    }
    if (c == '>' || c == '@')
        bf->last_char = c;
    if (read->length + 1 >= read->max)
    {
        read->max = read->length + 2;
        nearest_power_of_2(read->max);
        read->string = (char*) realloc(read->string, read->max);
    }
    read->string[read->length] = '\0';
    if (c == '+') // fastq
    {
        if (dummy->max < read->max) // resize quality to match read length
        {
            dummy->max = read->max;
            dummy->string = (char*)realloc(dummy->string, dummy->max);
        }
        while ( (c = buffered_getc(bf)) != -1 && c != '\n'); // read rest of quality comment
        while (buffered_gets(bf, dummy, NULL, true, true) >= 0 && dummy->length < read->length); // read rest of quality
        bf->last_char = 0;
    }
    *len = read->length;
    *nseq = read->string;
    if (cheader && hlen)
    {
        *cheader = header->string;
        *hlen = header->length;
    }

    return true;
}

// wrapper
bool  Bank::get_next_seq_from_file(char **nseq, int *len, int file_id)
{
    return get_next_seq_from_file(nseq,NULL,len,NULL,file_id);
}

// wrapper
bool Bank::get_next_seq(char **nseq, char **cheader, int *len, int *hlen, int * id_file)
{
    * id_file = index_file;
    bool success = get_next_seq_from_file(nseq,cheader,len,hlen,index_file);
    if (success)
        return true;
    
    // cycle to next file if possible
    if ( index_file < nb_files-1 )
    {
        close_stream(index_file);
	index_file++;
        open_stream(index_file);
	return get_next_seq(nseq,cheader, len,hlen, id_file);
    }
    return false;
}

// wrapper
bool Bank::get_next_seq(char **nseq, char **cheader, int *len, int *hlen)
{
    bool success = get_next_seq_from_file(nseq,cheader,len,hlen,index_file);
    if (success)
        return true;
    
    // cycle to next file if possible
    if ( index_file < nb_files-1 )
    {
        close_stream(index_file);
        index_file++;
        open_stream(index_file);
        return get_next_seq(nseq,cheader, len,hlen);
    }
    return false;
}

//wrapper
bool Bank::get_next_seq(char **nseq, int *len)
{
  return get_next_seq(nseq,NULL,len,NULL);
}
//wrapper
bool Bank::get_next_seq(char **nseq, int *len, int * id_file)
{
    return get_next_seq(nseq,NULL,len,NULL,id_file);
}
//function for subset of get_next_seq
bool Bank::get_next_seq_subset( char **nseq, char **cheader, int *len, int *hlen)
{
	// Only read a subset of the nb_files 5 % of the files to estimate the k-mer counts
	bool success = get_next_seq_from_file(nseq,cheader,len,hlen,index_file);
	
	if(success)
		return true;
		
	int limit = nb_files/20;
	//printf("Limit of files to process %d \n",limit);
	if(index_file < limit -1)
	{
		close_stream(index_file);
		index_file++;
		open_stream(index_file);
		return get_next_seq_subset(nseq,cheader,len,hlen);
	}
	return false;
}

// had to move the Bank(x,x) constructor to an init() to avoid calling a constructor inside the Bank(x) constructor
void Bank::init(char **fname, int nb_files_)
{
    int64_t i;
    nb_files = nb_files_;
    filesizes = 0;

    // open the reads file, don't know if it is a fasta/q file or a list of file names yet
    gzFile tempfile = gzopen(fname[0],"r");
    if (tempfile == NULL)
    {
        char *buffer = (char*)malloc(BUFSIZ);
        char * errorMessage = (char * ) strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s   %s \n",errorMessage,fname[0]);
        free(buffer);
        exit(1);
    }
    char deb=(char)gzgetc(tempfile);

    char **nfname;// [MAX_NB_FILES][TAILLE_NOM];
    nfname = (char**) malloc(sizeof(char*)*MAX_NB_FILES);
    for(int jj=0; jj<MAX_NB_FILES; jj++ )
	nfname [jj] =  (char*) malloc(sizeof(char)*TAILLE_NOM);
      
    if(deb=='>' || deb=='@' || deb==EOF)
    { // file is a fasta/q file
        gzclose(tempfile);
    }
    else // file contains a list of file names
    {
        char* ret;
        gzungetc(deb,tempfile);
        printf("File %s starts with character \"%c\", hence is interpreted as a list of file names\n",fname[0],deb );
        int ii;
        // get the filenames
        for (ii=0; ii<MAX_NB_FILES ; ii++)
        {
            ret = gzgets(tempfile, nfname[ii], BUFFER_SIZE);
            if (ret != NULL) {
                // remove \r \n chars
                char *endline = strchr(nfname[ii], '\n');
                if (endline)
                    *endline='\0';
                endline = strchr(nfname[ii], '\r');
                if (endline)
                    *endline='\0';
            }
            else // no more filenames
                break;
        }
        printf("Reading %i read files\n",ii);
        if(ii==MAX_NB_FILES)
            printf("Warning! using max number of read files (%i)\n",ii);

        nb_files = ii;
        fname = (char **) nfname;
        gzclose(tempfile);

    }
    //initialize the buffers
    buffered_file = (buffered_file_t**) malloc(sizeof(buffered_file_t *)*nb_files);
    for (i=0; i<nb_files; i++)
    {
        buffered_file[i] = (buffered_file_t *)calloc(1, sizeof(buffered_file_t));
        buffered_file[i]->buffer = (unsigned char*) malloc(BUFFER_SIZE); 
        buffered_file[i]->fname = strdup(fname[i]);
    }
    // estimate total size of files
    for (i=0; i<nb_files; i++)
    {
        bool compressed = false;
        uint64_t estimated_filesize;

        if (strstr(fname[i],"gz") == (fname[i]+strlen(fname[i])-2) ) compressed=true;
        if (compressed)
            // crude hack, based on Quip paper reporting compression ratio (~0.3). 
            // gzseek(SEEK_END) isn't supported. need to read whole file otherwise :/
            estimated_filesize = fsize(fname[i]) * 4; 
        else
            estimated_filesize = fsize(fname[i]);

        buffered_file[i]->estimated_filesize = estimated_filesize;
        filesizes += estimated_filesize;
    }
    //    printf("Filesize estimated from reads is %llu\n",filesizes );

    // open each file for reading
    rewind_all(); // initialize the get_next_seq iterator to the first file

    // init read and dummy (for readname and quality)
    read = (variable_string_t*) calloc(1,sizeof(variable_string_t));
    dummy = (variable_string_t*) calloc(1,sizeof(variable_string_t));
    header = (variable_string_t*) calloc(1,sizeof(variable_string_t));

    
    for(int jj=0; jj<MAX_NB_FILES; jj++ )
      free	(nfname [jj]); 
    free(nfname);
}

Bank::Bank(char *fname0)
{
    char *fname[1] = { fname0 };
    init(fname, 1);
}


Bank::Bank(char **fname, int nb_files_)
{
    init(fname,nb_files_);
}

Bank::~Bank(){
    variable_string_t * to_free[3] = {read, dummy, header};
    for (int i = 0; i < 3; i++)
    {
        if (to_free[i])

        {
            if (to_free[i]->string)
                free(to_free[i]->string);
            free(to_free[i]);
        }
    }
    for (int i=0; i<nb_files; i++)
    {
        free(buffered_file[i]->buffer);
        free(buffered_file[i]);
    }
}

void Bank::close()
{
    for (int i=0; i<nb_files; i++)
        gzclose(buffered_file[i]->stream);
}
void Bank::count_kmers_for_small_value(int l, double *lmer_counts)
{
	printf("Reestimating partitions sizes and number of passes based on %d-mers \n",l);
	char *rseq;
	int readlen;
	long total_bins = pow(4,l);
	int kmer_nbits = sizeof(kmer_type)*8;
	rewind_all();
	//long lmer_counts[total_bins];
	char * pt_begin; //used from conversion of reads to binary format
	int idx=0;
	long lmer_mask = (((long)1)<<(l*2)) -1;
	//initialize all lmer_counts to zero
	uint64_t total_counts=0;
	for(int i=0;i<total_bins;i++)
		lmer_counts[i]=0;
		uint64_t multiplication_factor ;
	if(filesizes > 21474836480 || nb_files>20)
	{	
		//printf("Limit of files to read %d \n",nb_files/20);
		while(get_next_seq_subset(&rseq,NULL,&readlen,NULL))
		{
			//count l-mers from the read and store them in lmer_counts
			pt_begin = rseq;	
			while(pt_begin < (rseq + readlen))
			{
				idx=0;
				//skips NN
				while(*pt_begin =='N' && pt_begin < (rseq+readlen))
					pt_begin++;
				// goes to next N or end of seq
				while( (pt_begin[idx] !='N') && ((pt_begin+idx) < (rseq + readlen)) )
				{	
					idx++;
				}
				// we have  a seq begining at pt_begin of size idx, if idx > l, count l-mers and store them in the array
				if(idx>=l)
				{
					//update kmer counts
					int c_pt = 0; int lmer_to_read =idx-l+1;
					long seed_graine = seed_lmer(pt_begin,l); //temporary
					lmer_counts[seed_graine]++;
					total_counts++;
					c_pt++; pt_begin+=l;
					while(c_pt<lmer_to_read)
					{
						seed_graine = (seed_graine*4 + NT2int(pt_begin[0])) & lmer_mask;
						lmer_counts[seed_graine]++;
						total_counts++;
						c_pt++; pt_begin++;
					}
				}else{
					pt_begin +=idx;
				}
			}
		}
		multiplication_factor = 20;
	}else
	{
		while(get_next_seq(&rseq,NULL,&readlen,NULL))
		{
			//count l-mers from the read and store them in lmer_counts
			pt_begin = rseq;	
			while(pt_begin < (rseq + readlen))
			{
				idx=0;
				//skips NN
				while(*pt_begin =='N' && pt_begin < (rseq+readlen))
					pt_begin++;
				// goes to next N or end of seq
				while( (pt_begin[idx] !='N') && ((pt_begin+idx) < (rseq + readlen)) )
				{	
					idx++;
				}
				// we have  a seq begining at pt_begin of size idx, if idx > l, count l-mers and store them in the array
				if(idx>=l)
				{
					//update kmer counts
					int c_pt = 0; int lmer_to_read =idx-l+1;
					long seed_graine = seed_lmer(pt_begin,l); //temporary
					lmer_counts[seed_graine]++;
					total_counts++;
					c_pt++; pt_begin+=l;
					while(c_pt<lmer_to_read)
					{
						seed_graine = (seed_graine*4 + NT2int(pt_begin[0])) & lmer_mask;
						lmer_counts[seed_graine]++;
						total_counts++;
						c_pt++; pt_begin++;
					}
				}else{
					pt_begin +=idx;
				}
			}
		}
		multiplication_factor = 1;
	}
	// reestimate the number of partitions file required in total based on volume per partition limitation 
	uint64_t total_partitions=0;
	//int volume_size = estimate_kmers_volume(sizeKmer);
	//printf("Total counts is %llu and size is %d \n",total_counts,volume_size);
	//printf("Factor for multiplication %llu\n",multiplication_factor);
	if( multiplication_factor ==0 )
		multiplication_factor = 1;
	//printf("Factor for multiplication %llu\n",multiplication_factor);
	for(int i=0;i<total_bins;i++)
	{
	//	printf("%d %lu\n",i,lmer_counts[i]);
		//ensure that each partition file doesn't goes over 0.75 of its volume_per_partitions, due to skew in homopolymer chains
		lmer_counts[i] = lmer_counts[i]*kmer_nbits*multiplication_factor/(1.0*1024*1024*8); // converting counts to MB
//		lmer_counts[i]= ceil(lmer_counts[i]*kmer_nbits/(0.75*volume_per_partition)/1024/1024); 
//		lmer_counts[i] = ceil(log(lmer_counts[i])/log(4));
//		lmer_counts[i]= pow(4,lmer_counts[i]);
//		total_partitions += (uint64_t) lmer_counts[i];
	}
	//return total_partitions;
}


uint64_t Bank::reestimate_partitions(int l)
{
//Depricated May 2014
//reestimate the partitions size based on low kmer counts, as homo polymers will be more frequent in the database
	printf("Re-estimating partitions sizes and number of passes based on %d-mers \n",l);
	char * rseq;
	int readlen;
//	long total_bins = pow(4,l); // total number of bins in which counts of l-mers will be stored
	int kmer_nbits = sizeof(kmer_type)*8; // because a k
	rewind_all();
//	long lmer_counts[total_bins];
	long lmer_counts[4]; // most common lmer counts will be all A's, or C with all A's .... so just compute them for reestimation
	long lmers[4];	
	lmers[0]=0; lmers[1]= (((long)1)<<(2*l-1)); lmers[2] =  (((long)1)<<(2*l)); lmers[3]= lmers[2]|lmers[1];
	char * pt_begin;
	int idx=0;
	long lmer_mask = (((long)1)<<(l*2)) - 1;
	//initialize counts to zero
	for (int i=0;i<4;i++)
		lmer_counts[i]=0;
	while(get_next_seq(&rseq,&readlen))
	{
		//count l-mers from the read and store them in lmer_counts
		pt_begin = rseq;	
		while(pt_begin < (rseq + readlen))
		{
			idx=0;
			//skips NN
			while(*pt_begin =='N' && pt_begin < (rseq+readlen))
				pt_begin++;
			// goes to next N or end of seq
			while( (pt_begin[idx] !='N') && ((pt_begin+idx) < (rseq + readlen)) )
			{	
				idx++;
			}
			// we have  a seq begining at pt_begin of size idx, if idx > l, count l-mers and store them in the array
			if(idx>=l)
			{
				//update kmer counts
				int c_pt = 0; int lmer_to_read =idx-l+1;
				long seed_graine = seed_lmer(pt_begin,l); //temporary
				//lmer_counts[seed_graine]++;
				for(int j=0;j<4;j++)
				{
					if(seed_graine==lmers[j])
						lmer_counts[j]++;
				}
				c_pt++; pt_begin+=l;
				while(c_pt<lmer_to_read)
				{
					seed_graine = (seed_graine*4 + NT2int(pt_begin[0])) & lmer_mask;
					//lmer_counts[seed_graine]++;
					for(int j=0;j<4;j++)
					{
						if(seed_graine==lmers[j])
							lmer_counts[j]++;
					}
					c_pt++; pt_begin++;
				}
			}else{
				pt_begin +=idx;
			}
		}
	}
	// find max value in the list of lmer_counts
	long max=-10, max_pos = -10;
	for(int i=0;i<4;i++)
	{
	//	printf(" %d %ld \n", i, lmer_counts[i]);
		if(max<lmer_counts[i]) 
		{
			max = lmer_counts[i];
			max_pos = i;
		}
	}
	printf("Max count of %d mer is %ld kmers %ld \n", l,max,max_pos);
	max = max*kmer_nbits/1024/1024/8;
	printf("Max count of %d mer is %ld MB %ld \n", l,max,max_pos);
	return max;
}
// estimate the volume of all redundant kmers in the reads, if they were to be stored in 2bits
uint64_t Bank::estimate_kmers_volume(int k)
{
	return estimate_kmers_volume(k,sizeof(kmer_type)*8);
}
uint64_t Bank::estimate_kmers_volume(int k, int kmer_nbits)
{

    char * rseq;
    int readlen;
    int NbRead = 0;
    //int kmer_nbits = std::max(64,(int)pow(2,ceilf(log2f(2*k)))); // Bank assumes that a kmer is stored in the smallest integer type (e.g. uint64_t or uint128_t) // not accurate anymore with _ttmath/_largeint
    //int kmer_nbits = sizeof(kmer_type)*8;
    rewind_all();
    uint64_t total_volume = 0;

    while ( index_file < nb_files ) 
    {
	open_stream(index_file);
	int NbRead = 0;
	uint64_t volume_for_file = 0;
	while (get_next_seq_from_file(&rseq,NULL,&readlen,NULL,index_file))
    	{
        	if (readlen >= k)
      	     		volume_for_file += (readlen-k+1) * (uint64_t) kmer_nbits;
        	if (NbRead++ == 100000)
            		break;
    	}
  	if ( gztell(buffered_file[index_file]->stream) != 0) // empty file
    	{    //return 1;

            volume_for_file = (uint64_t) ( ( (float) volume_for_file ) * ( ( (float)(buffered_file[index_file]->estimated_filesize)) / ((float) gztell(buffered_file[index_file]->stream)) ) );
            total_volume += volume_for_file;

    	}
	close_stream(index_file);
	index_file++;
    }
    total_volume = total_volume / 1024 /1024 /8; // put it in MB
    
    if (total_volume == 0)  // tiny files fix
        total_volume = 1;

    rewind_all();
    return total_volume;
}

// estimate the number of reads
uint64_t Bank::estimate_nb_reads()
{
    char * rseq;
    int readlen;
    int NbRead = 0;
    rewind_all();
    
    uint64_t volume = 0;
    while (get_next_seq(&rseq,&readlen))
    {
        volume += 1;
        if (NbRead++ == 1000)
            break;
    }

    if ( gztell(buffered_file[index_file]->stream) == 0) // empty file
        return 1;

    volume = (volume * filesizes) / gztell(buffered_file[index_file]->stream); // linear extrapolation from the first 1k reads 

    rewind_all();
    return volume;
}

// estimate maximum read length
// from the first 10000 reads of each file
int Bank::estimate_max_readlen()
{
    char * rseq;
    int readlen;
    rewind_all();
    int max_readlen = 0;
    uint64_t volume = 0;

    index_file = 0;
    
    while ( index_file < nb_files )
    {
        int NbRead = 0;
        while (get_next_seq_from_file(&rseq,NULL,&readlen,NULL,index_file))
        {
            max_readlen = max(readlen, max_readlen);
            if (NbRead++ == 10000)
                break;
        }
        index_file++;
    } 

    index_file = 0;
    return max_readlen;
}

void Bank::save_position()
{
    restore_index_file = index_file;
    restore_pos = gztell(buffered_file[index_file]->stream) - buffered_file[index_file]->buffer_end + buffered_file[index_file]->buffer_start;
}

void Bank::load_position()
{
    close_stream(index_file);
    index_file = restore_index_file;
    open_stream(index_file);
    gzseek(buffered_file[index_file]->stream, restore_pos, SEEK_SET);
    buffered_file[index_file]->eof = false;
    rebuffer(buffered_file[index_file]);
}

// BinaryBank: a binary file containing kmers

BinaryBank::BinaryBank(char *given_filename, int given_sizeElement, bool write) : sizeElement(given_sizeElement)
{
    strcpy(filename,given_filename);
    open(write);
    buffer_size_nelem= (WRITE_BUFFER/given_sizeElement);
    buffer = (void *) malloc(given_sizeElement * buffer_size_nelem);
    cpt_buffer=0;
}


BinaryBankConcurrent::BinaryBankConcurrent(char *given_filename, int given_sizeElement, bool write, int given_nthreads) : BinaryBank(given_filename,given_sizeElement,write) 
{
    nthreads = given_nthreads;
    
    //free(buffer); buffer =NULL; //cannot do that
    bufferT = (void **) malloc(sizeof(void*) * nthreads);

    for (int i= 0; i< nthreads; i++)
    {
         ((void ** )bufferT)[i]= (void *) malloc( WRITE_BUFFER);
      //  ((void ** )bufferT)[i]= (void *) malloc(sizeElement* WRITE_BUFFER);

    }
    cpt_buffer_tid = (int  *)malloc(sizeof(int) * nthreads);
    memset (cpt_buffer_tid,0,sizeof(int) * nthreads);
}


void BinaryBankConcurrent::write_element_buffered( void *element, int tid)
{
    write_buffered(element,sizeElement,tid);
}


void BinaryBankConcurrent::write_buffered( void *element, int size, int tid)
{
    write_buffered( element, size, tid, true);
}

void BinaryBankConcurrent::write_buffered( void *element, int size, int tid, bool can_flush)
{
    if(cpt_buffer_tid[tid]>= WRITE_BUFFER -100 && can_flush)
    {
        flush(tid);
    }
    
    char * buf_pt = ((char**) bufferT)[tid];    
    memcpy(buf_pt + cpt_buffer_tid[tid] , element, size);
    
    cpt_buffer_tid[tid]+=size;
    // cpt_buffer_tid[tid]++;
    

}



void BinaryBankConcurrent::flush(int tid)
{
    flockfile(binary_read_file);
    if (!fwrite( ((void **)bufferT)[tid], 1, cpt_buffer_tid[tid], binary_read_file))            
    {
        printf("error: can't fwrite (disk full?)\n");
        funlockfile(binary_read_file);
        exit(1);
    }
    cpt_buffer_tid[tid]=0;
    funlockfile(binary_read_file);

}

int BinaryBankConcurrent::find_error()
{
	return ferror(binary_read_file);
}
//should be called by only one of the threads
void BinaryBankConcurrent::close()
{
    //flush buffer // if close Bank in read mode with data in the readbuffer, will result in error
    for(int ii=0; ii< nthreads; ii++)
    {
        if(cpt_buffer_tid[ii])
        {
//	printf("In close function close 1 \n");
            if (!fwrite(((void **)bufferT)[ii], 1, cpt_buffer_tid[ii], binary_read_file))
          //      if (!fwrite(((void **)bufferT)[ii], sizeElement, cpt_buffer_tid[ii], binary_read_file))

            {
                printf("error: can't fwrite (disk full?)\n");
                exit(1);
            }
        }
//	printf("In close function close 2 \n");
        cpt_buffer_tid[ii]=0;
    }
    try{ 
    	//if(!feof(binary_read_file))
	 fclose(binary_read_file);
	//else 
	 //printf("End of binary file but closing is  problem \n");
    }
    catch (int e) 
    {
	printf("Error number %d \n",e);
	throw;	
    }
//	printf("In close function close 3\n");
}


void BinaryBank::write_element( void *element)
{
  //  flockfile(binary_read_file);
   // fprintf(stderr,"write elem %lli \n",*(int64_t *)element);
    if (!fwrite(element, sizeElement, 1, binary_read_file))
    {
       // funlockfile(binary_read_file);
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
  //  funlockfile(binary_read_file);
}


void BinaryBank::write_element_buffered( void *element)
{
    
    if(cpt_buffer==buffer_size_nelem)
    {
        if (!fwrite(buffer, sizeElement, buffer_size_nelem, binary_read_file))
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
        cpt_buffer=0;
    }
    
    
    //((kmer_type *)buffer)[cpt_buffer]= *((kmer_type *)element);
    memcpy((unsigned char *)buffer + (cpt_buffer * sizeElement), element, sizeElement);
    cpt_buffer++;
    
}



size_t BinaryBank::read_element( void *element)
{
    return fread(element, sizeElement,1, binary_read_file);
}

size_t BinaryBank::read_element_buffered( void *element)
{
    if(cpt_buffer==0)
    {
        if(feof(binary_read_file)) { //printf( "End of file reached 1\n"); 
		return 0; }
	cpt_buffer=fread(buffer, sizeElement,buffer_size_nelem, binary_read_file);
	//printf("number of elements read is %d %d \n",cpt_buffer,buffer_size_nelem);
	if (cpt_buffer==0) return 0;
        //if(feof(binary_read_file)) { printf( "End of file reached 2 \n"); return 0; }
        cpt_init_buffer = cpt_buffer;
    }
    memcpy(element, (unsigned char *)buffer + (cpt_buffer-1) * sizeElement, sizeElement);
    //printf("number of the cpt_buffer is %d \n",cpt_buffer);
    cpt_buffer --;
    return cpt_buffer+1; // nb remaining before read
}

/*size_t BinaryBank::read_element_buffered( void *element,int size) 
{
	if(cpt_buffer==0)
	{
		cpt_buffer=fread(buffer, sizeElement, buffer_size_nelem, binary_read_file);
		if( cpt_buffer==0)	return 0;
	}

	memcpy(element, (unsigned char*)buffer + (cpt_buffer-1)*size,size);
	cpt_buffer--;
	return cpt_buffer+1;
}*/
// used to read/write raw information to the binary file (e.g. kmer count)

void BinaryBank::write( void *element, int size)
{
    if (!fwrite(element, size, 1, binary_read_file))
    {
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
}

size_t BinaryBank::read( void *element, int size)
{
    return fread(element, size,1, binary_read_file);
}


void BinaryBank::rewind_all()
{
    rewind(binary_read_file);
}

void BinaryBank::close()
{
    //flush buffer // if close Bank in read mode with data in the readbuffer, will result in error
    if(cpt_buffer)
    {
    if (!fwrite(buffer, sizeElement, cpt_buffer, binary_read_file))
    {
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
    }
    cpt_buffer=0;
    
    fclose(binary_read_file);
}

void BinaryBank::open(bool write)
{
    binary_read_file = fopen(filename,write?"wb":"rb");
    if( binary_read_file == NULL )
    {
        char *buffer = (char*)malloc(BUFSIZ);
        char * errorMessage = (char * ) strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s  write %i  %s\n",errorMessage,write,filename);
        free(buffer);
        exit(1);
    }

}

off_t BinaryBank::nb_elements()
{
  return fsize(filename)/sizeElement;
}


BinaryBank::~BinaryBank()
{
    if(buffer!=NULL)
    {
        free (buffer); //buffer =NULL;
    }
}


BinaryBankConcurrent::~BinaryBankConcurrent()
{
    
    for (int i= 0; i< nthreads; i++)
    {
        free(((void ** )bufferT)[i]);
        ((void ** )bufferT)[i]=NULL;
    }
    free(bufferT);
}



/////////////class BinaryReads a file containing reads

BinaryReads::~BinaryReads()
{
    free (buffer); buffer = NULL;
}


BinaryReads::BinaryReads(char *given_filename,  bool write)
{
    read_write_buffer_size = BINREADS_BUFFER;
    strcpy(filename,given_filename);
    open(write);
    buffer = (unsigned char *) malloc(read_write_buffer_size*sizeof(unsigned char));
    cpt_buffer = 0;
}


void BinaryReads::rewind_all()
{
    rewind(binary_read_file);
}

void BinaryReads::close()
{
    unsigned int block_size =0;
    //flush buffer
    if(cpt_buffer)
    {
        //printf("close :write block %i \n",cpt_buffer);
        block_size = cpt_buffer;
        fwrite(&block_size, sizeof(unsigned int), 1, binary_read_file); // block header
        if (!fwrite(buffer, 1, cpt_buffer, binary_read_file))
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
    }
    cpt_buffer=0;
    
    fclose(binary_read_file);
}

void BinaryReads::open(bool write)
{
    binary_read_file = fopen(filename,write?"wb":"rb");
    if( binary_read_file == NULL )
    {
        char *buffer = (char*)malloc(BUFSIZ);
        char * errorMessage = (char * ) strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s  write %i  %s\n",errorMessage,write,filename);
        free(buffer);
        exit(1);
    }
    
}



//format is
// 32 bit integer = readlen,  then seq in binary
// then next read..
//32 bit len is overkill but simpler
//also makes buffer then write block with header : size of block to read, with n reads .... will allow large fread when reading this file ...
void BinaryReads::write_read(char * read, int readlen)
{
    int tai = readlen;
    unsigned char rbin;
    char * pt = read;
    unsigned int block_size = 0;
    
 //   printf("write read %i / %i   readlen %i \n",cpt_buffer,read_write_buffer_size,readlen);
    //todo : also flush to disk  sometimes (ie if very large buffer, to create smaller blocks..)
    if(cpt_buffer >= (read_write_buffer_size-readlen) || cpt_buffer > 10000000 )  ////not enough space to store next read   true space is 4 + readlen/4 + rem
        //flush buffer to disk
    {
        
        block_size = cpt_buffer;
        
        //printf("write block %i\n",block_size);
        if(block_size) fwrite(&block_size, sizeof(unsigned int), 1, binary_read_file); // block header
        if (!fwrite(buffer, 1, cpt_buffer, binary_read_file)) // write a block, it ends at end of a read
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
        cpt_buffer=0;
    }
    
    //check if still not enough space in empty buffer : can happen if large read, then enlarge buffer
    if(read_write_buffer_size < readlen)
    {
        read_write_buffer_size = 2*readlen; // too large but ok
        buffer =  (unsigned char *) realloc(buffer,sizeof(unsigned char) * read_write_buffer_size);
    }
    
    memcpy(buffer+cpt_buffer,&readlen,sizeof(int));
    cpt_buffer+= sizeof(int);
    
    //fwrite( (void *) &readlen, sizeof(int), 1, binary_read_file);

    
    for (tai=readlen; tai>=4  ; tai-=4)
    {
        rbin = code4NT(pt);
      //  fwrite((void *) &rbin, 1,1,binary_read_file );
        buffer[cpt_buffer]=rbin; cpt_buffer++;
        pt +=4;
    }
    
    //then remaining
    if(tai)
    {
        rbin = code_n_NT(pt,tai);
       // fwrite( (void *) &rbin,1,1,binary_read_file);
        buffer[cpt_buffer]=rbin; cpt_buffer++;
    }
}



void  compute_kmer_table_from_one_seq(int readlen, char * seq, kmer_type * kmer_table )  //,char * pkmer_table //pour remplissage table loc
{
    kmer_type graine = codeSeed(seq);
    kmer_type graine_revcomp = revcomp(graine);
    kmer_table[0] = min(graine,graine_revcomp);
    seq++;
    for (int i=1; i<readlen-sizeKmer+1; i++)
    {
        graine =   (graine * 4 + NT2int(seq[sizeKmer-1])) & kmerMask   ;
        graine_revcomp =  ((graine_revcomp >> 2) +  ( ((kmer_type) comp_NT[NT2int(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask ;
        kmer_table[i] = min(graine,graine_revcomp);
        seq++;
    }
}

void compute_kmer_table_from_one_seq(int readlen, char * seq, kmer_type * kmer_table,int * kmer_length_table, int sizesmallest)
{
	int diffK = sizeKmer - sizesmallest;
	int smallK_read = 0;
	kmer_type graine;
	int kmers_read = 0;
	
	while(smallK_read<diffK) 
	{
		graine = codeSeed(seq,sizesmallest + smallK_read,kmerMask);
		kmer_table[kmers_read] = graine;
		kmer_length_table[kmers_read] = sizesmallest + smallK_read;
		kmers_read++;
		smallK_read++;
	}
	graine = codeSeed(seq);
	kmer_table[kmers_read] = graine;
	kmer_length_table[kmers_read] = sizeKmer;
	seq++;
	kmers_read++;
	for (int i=1;i<readlen-sizeKmer+1;i++)
	{
        	graine =   (graine * 4 + NT2int(seq[sizeKmer-1])) & kmerMask   ;
		kmer_table[kmers_read] = graine;
		kmer_length_table[kmers_read] = sizeKmer;
		seq++;kmers_read++;
	}
}



////kmers buffer


KmersBuffer::KmersBuffer(BinaryReads *bfile, int  pbuffer_size, int nseq_task )
{
    
    read_write_buffer_size = BINREADS_BUFFER;
    buffer = ( char *) malloc(read_write_buffer_size*sizeof( char));
    cpt_buffer = 0;
    cpt_binSeq_read =0; binSeq_toread =0; 
    cpt_read_count =0; //Added by raunaq, to make sure that from every read we keep also get smaller kmers from the begining of the read
    smallK = sizeKmer-smallestKmer;
    max_read_length = KMERSBUFFER_MAX_READLEN;
    binfile = bfile;
    buffer_size = pbuffer_size;
    kmers_buffer =  (kmer_type *) malloc(sizeof(kmer_type) * buffer_size);
    //kmer_length = (kmer_type *) malloc(sizeof(kmer_type) * buffer_size); // Added by raunaq; so that we have a length of kmer stored along with the kmer
    kmer_length = (int*) malloc(sizeof(int) * buffer_size); // Added by raunaq; so that we have a length of kmer stored along with the kmer
   // binSeq =  (char *) malloc(sizeof(char) * max_read_length); // no need to alloc ram for binse : will points to buffer
    binSeq_extended =  (char *) malloc(sizeof(char) * max_read_length);
    blocksize_toread =0;


    nseq_step = nseq_task;
    binary_read_file = bfile->binary_read_file;
    
}


void KmersBuffer::reset_max_readlen(int read_length)
{

    max_read_length = read_length;

  //  binSeq =  (char *) realloc(binSeq,sizeof(char) * max_read_length);
    binSeq_extended =  (char *) realloc(binSeq_extended,sizeof(char) * max_read_length);

}

 KmersBuffer::~KmersBuffer()
{
    free (kmers_buffer);
    free(buffer);
    free(kmer_length); // Added by Raunaq
    //free(binSeq);
    free(binSeq_extended);
}


//now returns number of kmers read
int KmersBuffer::readkmers()
{
    
    int llen;
    int * len = & llen ;
    unsigned int block_size =0;
    int temp; //check if the buffer size block is set 

    //////reading new block from disk if needed
    // cpt_buffer == blocksize_toread    tells we finished reading previous buffer
    // (binSeq_toread  <= cpt_binSeq_read) tells  we finished reading the last sequence
    if(cpt_buffer == blocksize_toread  && (binSeq_toread  <= cpt_binSeq_read))
    {
        flockfile(binary_read_file);
	if( ! fread(&block_size,sizeof(unsigned int),1, binary_read_file)) //read block header
        {
            funlockfile(binary_read_file);
            return -1; // no more blocks to read
        }
        
        
        if(block_size >= read_write_buffer_size) // block buffer need to be enlarged
        {
            read_write_buffer_size = 2*block_size;
            buffer =  ( char *) realloc(buffer,sizeof( char) * read_write_buffer_size);
        }
        
        temp =fread(buffer,sizeof( char),block_size, binary_read_file); // read a block of sequences into the buffer
        //printf("Block size to read %d total read from file %d \n",block_size,temp);
	funlockfile(binary_read_file);
        cpt_buffer = 0;
        blocksize_toread = block_size;        
    }
    ///////////////////////
    
    //	fprintf(stderr,"Block size to read %d current buffer position %d\n",blocksize_toread,cpt_buffer);
    
    
    
    //now parse the whole block in ram
    
    int i,j;
    int nchar;
    int diffK = sizeKmer-smallestKmer; //Added by Raunaq  
    unsigned char fournt;
    char code2NT[4] = {'A','C','T','G'};
    char temp_seq[512];

    nkmers = 0;
    int nseq_lues = 0;
    //cpt_buffer : how much we have already read in the buffer
    //blocksize_toread : how much there is to read in the buffer
    while(cpt_buffer < blocksize_toread || ( binSeq_toread  > cpt_binSeq_read)) //while work to do
    {
        
        if( binSeq_toread <= cpt_binSeq_read)// read new sequence if needed  //we  will put one sequence into binSeq_extended
        {
	    memcpy(len,buffer+cpt_buffer,sizeof(int)); // the sequence length
            cpt_buffer += sizeof(int);
            nseq_lues ++;
            
            if( (*len) > max_read_length) reset_max_readlen((int)(1.2*(*len))); // resize memory for sequence if needed
            nchar = ((*len)+3)/4; // number of bytes used to encode the sequence in its binary format (4 nt per byte)

            binSeq = buffer + cpt_buffer; // point binseq to correct place //cpt_buffer ==  where we are now in the buffer
            cpt_buffer += nchar;
            if( (*len) >= smallestKmer) 
	    {		
       		// on disk data was encoded with 4 nucleotides per bytes,
            	// here we expand one sequence  to one nucl per byte into binSeq_extended
            	//nucleotides are still encoded in [0-3]
            	j=0;
         	for(i=0; i<nchar; i++)
            	{
                	fournt = binSeq[i]; 
          	      	binSeq_extended[j+3]=fournt & 3; fournt = fournt >> 2;  temp_seq[j+3] = code2NT[binSeq_extended[j+3]];
                	binSeq_extended[j+2]=fournt & 3; fournt = fournt >> 2;  temp_seq[j+2] = code2NT[binSeq_extended[j+2]];
                	binSeq_extended[j+1]=fournt & 3; fournt = fournt >> 2;  temp_seq[j+1] = code2NT[binSeq_extended[j+1]];
                	binSeq_extended[j+0]=fournt & 3;  temp_seq[j+0] = code2NT[binSeq_extended[j+0]];
                	j+=4;
            	}
      	      //printf("Reading a new read from the file %d cpt_buffer %d len\n");
      	      binSeq_toread = *len-sizeKmer+1; // binSeq_toread tells how many kmers there are in this sequence
	    	if((*len)<sizeKmer) 
			smallK = *len-smallestKmer+1;
		else 
			smallK = diffK;
	      //if(binSeq_toread<0)
		//printf(">%d_%d_%d_%d\n%s\n",nseq_lues,*len,binSeq_toread,smallK,temp_seq);
	      //binSeq_toread = *len-smallestKmer+1; // binSeq_toread tells how many kmers there are in this sequence including the last set of smaller size kmers which help in counting kmers for smaller values of k. 
	      //fprintf(stderr,"%d Smallest kmer size, %d length of reading kmers  %d ",smallestKmer, binSeq_toread,diffK);
	      cpt_binSeq_read = 0;  // tells how many kmers we have currently parsed in this sequence
              cpt_read_count = 0;    
	    }
        }
        
        
        {
            // binSeq_extended = beginning of the sequence,
            //  cpt_binSeq_read = how much we have already read in this sequence (when kmers_buffer is full, we can halt parsing kmers (see below) in the middle of a sequence, so this value is not necessarily 0)
            char *seq = binSeq_extended+cpt_binSeq_read;  
            kmer_type graine;
            kmer_type graine_revcomp;
	    if(cpt_read_count<smallK && cpt_binSeq_read==0 )
	    {
		graine = codeSeed_bin(seq,smallestKmer + cpt_read_count);
		if(nkmers>=buffer_size)
		{
			return nkmers;
		}
		kmers_buffer[nkmers] = graine;
		kmer_length[nkmers] = smallestKmer + cpt_read_count;
		nkmers++; cpt_read_count++;
            }// loop to list smaller kmers from the begining of the read 
	    while(cpt_read_count<smallK && cpt_binSeq_read==0)
	    {
		graine = graine*4 + seq[smallestKmer+cpt_read_count-1]; 
		//graine_revcomp = revcomp(graine,smallestKmer + cpt_read_count);
		if(nkmers>=buffer_size)
		{
			printf("Return from kbuff 1\n");
			return nkmers;
		}
		//kmers_buffer[nkmers] = min(graine,graine_revcomp);
		kmers_buffer[nkmers] = graine;
		kmer_length[nkmers] = smallestKmer + cpt_read_count;
		//printf("%s\n",print_kmer(kmers_buffer[nkmers],smallestKmer+cpt_read_count,kmerMask));
		nkmers++; cpt_read_count ++;
	    }
	    if( binSeq_toread - cpt_binSeq_read > 0 )
            //if( binSeq_toread > cpt_binSeq_read )
                // there are still unread kmers in this sequence, here we read the first one,
                // we put it in graine / graine_revcomp  and store it in the kmers_buffer
            {
		graine = codeSeed_bin(seq);
                //graine_revcomp = revcomp(graine);
                
                if(nkmers>=buffer_size)
                {
                   printf("Return from kbuff 2\n"); 
		   return nkmers;
                }
                //kmers_buffer[nkmers] = min(graine,graine_revcomp); //fprintf(stderr,"%u \n",kmers_buffer[nkmers]);
                kmers_buffer[nkmers] = graine;
		kmer_length[nkmers] = sizeKmer;
		//fprintf(stderr,"%d %" PRIu64 "\n",nkmers,kmers_buffer[nkmers]);
                //printf("%s\n",print_kmer(kmers_buffer[nkmers],sizeKmer,kmerMask));
		nkmers++; cpt_binSeq_read ++;
                seq++;
            } 
	   /* else 
	    {
		//fprintf(stderr,"We are here \n");
		graine = codeSeed_bin( seq, binSeq_toread + smallestKmer - 1 - cpt_binSeq_read );
		graine_revcomp = revcomp( graine, binSeq_toread + smallestKmer - 1 - cpt_binSeq_read );
	   	if(nkmers>=buffer_size) 
		{
			return nkmers;
		} 
		kmers_buffer[nkmers] = min(graine,graine_revcomp);
                //fprintf(stderr,"%d %" PRIu64 "\n",nkmers,kmers_buffer[nkmers]);
		nkmers++; cpt_binSeq_read ++;
		seq++;
	   }
          */
            
            while( binSeq_toread > cpt_binSeq_read ) //while there remains kmers to be read in this sequence
            {                
                graine =  (graine * 4 + (seq[sizeKmer-1])) & kmerMask ; //parse next nucleotide to construc the next kmer
                //graine_revcomp =  ((graine_revcomp >> 2) +  ( ((kmer_type) comp_NT[(int)(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask;
                //kmers_buffer[nkmers] = min(graine,graine_revcomp); 
                kmers_buffer[nkmers]  = graine;
		kmer_length[nkmers] = sizeKmer;
		//fprintf(stderr,"Kmer number %d Buffer value %" PRIu64 " Kmermask %" PRIu64 "\n",nkmers,kmers_buffer[nkmers],kmerMask);
                //printf("%s\n",print_kmer(kmers_buffer[nkmers],sizeKmer,kmerMask));
		nkmers ++; cpt_binSeq_read ++; //we store the kmer in the kmers_buffer
		seq++;
                if(nkmers>=buffer_size) //the kmers_buffer is full, we stop
                {                    
		    printf("Return from kbuff 3\n");
                    return nkmers;
                }
            }
	    /*while( binSeq_toread > cpt_binSeq_read )
	    {
		//Looking at end region of a read. This can be implemented using the kmerMask, but it wasn't too clear to me. 
		graine = codeSeed_bin( seq, binSeq_toread + smallestKmer - 1 - cpt_binSeq_read );
		graine_revcomp = revcomp( graine, binSeq_toread + smallestKmer - 1 - cpt_binSeq_read );
                kmers_buffer[nkmers] = min(graine,graine_revcomp);
		//fprintf(stderr,"End region: Total kmers %d Current position %d Buffer value %" PRIu64"\n",binSeq_toread,cpt_binSeq_read,kmers_buffer[nkmers]);
		nkmers++;
		cpt_binSeq_read++; 
		seq++; 
	   	if(nkmers>=buffer_size)
		{
			return nkmers;
		} 
	   }*/
            
        }
    }
    // we stop when we finished one block, or when kmers_buffer is full,
    // it can happen in the middle of a sequence : the next time we call readkmers we will have to continue
    // from where we stopped  in this sequence (counter cpt_binSeq_read tells us that) 

    //while buffer is non empty, we 'expand' a sequence into  binSeq_extended
    //then we parse binSeq_extended to store kmers in the kmers_buffer
    
    //printf("Outside readkmerloop, read %d reads \n",nseq_lues);
    return nkmers;
    
    
}

void  compute_partition_hashing_kmers(uint64_t no_partition,long * lmer_counts,int * partition_kmer_depth,uint32_t *match_partition)
{
	// store the counts in 
	long iter=0,iter_lmer=0;
	while(iter<no_partition)
	{
		for(int j=0;j<lmer_counts[iter_lmer];j++)
		{
			match_partition[iter] = (((uint32_t)j)<<(partition_kmer_depth[iter]*2)) | iter_lmer; //shift the current j by 2 * passes_hash+partitions_hash number and then or it with the current iter_lmer value
			partition_kmer_depth[iter]+=log(lmer_counts[iter_lmer])/log(4);
			iter++;
		}
		iter_lmer++;
	}
}

int compute_min_kmer_depth(uint64_t no_partition,int *partition_kmer_depth)
{
	int min=1000;
	for(int j=0;j<no_partition;j++)
	{
		if(partition_kmer_depth[j]<min)
			min=partition_kmer_depth[j];
	}
	return min;
}

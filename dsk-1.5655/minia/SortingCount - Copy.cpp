#include "SortingCount.h"
#include "inttypes.h"
#include <sys/resource.h> // for getrlimit()
#if OMP
#include "omp.h"
#endif

#define SINGLE_BAR 1
#define __STDC_FORMAT_MACROS    // Added by Raunaq

bool clear_cache = false; // clear file cache from memory (for timing only)

bool hybrid_mode = false;
bool use_hashing = true; // use hashing instead of sorting (better control of memory)
float load_factor = 0.7;
bool use_compressed_reads = true ; // true; // write compressed read file
extern int *Kmerlist;
extern unordered_map<int,int> kmerlength_map; //Map of largest length to smallest length 
bool output_histo;


// main k-mer counting function, shared between minia and dsk
// verbose == 0 : stderr progress bar
// verbose >= 1 : print basic status
// verbose >= 2 : print extra partition information
// write_count == True: include kmer count in results file, in that form:
//           - save kmer count for each kmer in the resulting binary file
//           - the very first four bytes of the result file are the kmer length
void sorting_count(Bank *Sequences, char *prefix, int max_memory, int max_disk_space, bool write_count, int verbose)
{

    // create a temp dir from the prefix
    char temp_dir[1024];
    sprintf(temp_dir,"%s_temp",prefix);

    // clear the temp folder (needs to be done before estimating disk space)
    DIR*            dp;
    struct dirent*  ep;
    char            p_buf[512] = {0};
    dp = opendir(temp_dir);
    while ( (dp != NULL) && ((ep = readdir(dp)) != NULL)) {
        sprintf(p_buf, "%s/%s", temp_dir, ep->d_name);
        remove(p_buf);
    }
    if(dp != NULL)
        closedir(dp);

    if (max_disk_space == 0)
    {
        // default max disk space
        struct statvfs buffer ;
        char current_path[1000];
        getcwd(current_path,sizeof(current_path));
        // int ret =
        statvfs(current_path, &buffer);
        int available = (int)(((double)buffer.f_bavail * (double)buffer.f_bsize) / 1024 / 1024);
	uint32_t tt_new_temp = (uint32_t) (((double)Sequences->filesizes)/(1024*1024));
        printf("Available disk space in %s: %d  %u %llu MB\n",current_path,available,tt_new_temp,Sequences->filesizes); // not working in osx (is that a TODO then?)
        max_disk_space = min((uint32_t)available/2, tt_new_temp);
    } 
    if (max_disk_space <= 0) // still 0?
        max_disk_space = 10000; // = default for osx

    // estimate number of iterations TODO Check if multiplication with totalKmers is actually required or not. It may be just increasing number of partitions for no reason
    //uint64_t volume = totalKmers*Sequences->estimate_kmers_volume(smallestKmer);  //Since there are totalKmers no of kmers and an upper bound can be estimated by using the smallest size of kmer. Added by Raunaq
    uint64_t volume = Sequences->estimate_kmers_volume(smallestKmer);  //Since there are totalKmers no of kmers and an upper bound can be estimated by using the smallest size of kmer. Added by Raunaq
    uint32_t nb_passes = ( volume / max_disk_space ) + 1;
    int passes_hash ;
    
    int nb_threads=1;
    
#if OMP
    use_compressed_reads =true;
    nb_threads = 8;
    max_memory /= nb_threads;
    max_memory = max (max_memory,1);
#endif
    
    // temp bugfix: don't use compressed reads for long reads
    if (Sequences->estimate_max_readlen() > 1000000)
        use_compressed_reads = false;
    
    
    uint64_t volume_per_pass,volume_per_partition;
    uint32_t nb_partitions;
    int partitions_hash;

    // loop to lower the number of partitions below the maximum number of simulatenously open files
    do
    {
        volume_per_pass = volume / nb_passes;
        nb_partitions = ( volume_per_pass * totalKmers / max_memory ) + 1; 
	//printf("volume per pass and total volume %llu %llu \n",volume_per_pass,(unsigned long long)volume);
        // if partitions are hashed instead of sorted, adjust for load factor
        // (as in the worst case, all kmers in the partition are distinct and partition may be slightly bigger due to hash-repartition)
        if (use_hashing)
        {
            nb_partitions = (uint32_t) ceil((float) nb_partitions / load_factor);
            nb_partitions = ((nb_partitions * OAHash::size_entry() ) + sizeof(key_type)-1) / sizeof(key_type); // also adjust for hash overhead
        }

        struct rlimit lim;
        int max_open_files = 1000;
        int err = getrlimit(RLIMIT_NOFILE, &lim);
        if (err == 0)
            max_open_files = lim.rlim_cur / 2;
        if (nb_partitions >= max_open_files)
            nb_passes++;
        else
            break;
    }
    while (1);
    volume_per_partition= volume_per_pass/nb_partitions;
    passes_hash = ceil(log(nb_passes)/log(4));
    partitions_hash = ceil(log(nb_partitions)/log(4));
    int size_for_reestimation = ceil((passes_hash + partitions_hash)*1.8);
    double * lmer_counts = (double * ) malloc(sizeof(long)*pow(4,size_for_reestimation));
    long * lmers_for_hash = (long * ) malloc(sizeof(long)*pow(4,size_for_reestimation));
    int * partitions_for_lmers =(int * ) malloc(sizeof(int)*pow(4,size_for_reestimation));
    Sequences->count_kmers_for_small_value(size_for_reestimation,lmer_counts);
    int temp_partition=reestimate_partitions(size_for_reestimation,volume_per_partition,lmer_counts,lmers_for_hash,partitions_for_lmers);
    unordered_map<long,int> part_hash;
    int total_lmers=pow(4,size_for_reestimation);
    for(int it=0;it<total_lmers;it++)
    {
	pair<long,int> temp_pair(lmers_for_hash[it],partitions_for_lmers[it]);
        part_hash.insert (temp_pair); // Add element to the hash 
    }
    //uint64_t up_passes_size = volume_per_pass;
      	do
	{
		//recompute the number of partitions based on updated partitions estimate
		nb_partitions = ceil(temp_partition*1.0/nb_passes);
		struct rlimit lim;
	        int max_open_files = 1000;
	        int err = getrlimit(RLIMIT_NOFILE, &lim);
	        if (err == 0)
        	    max_open_files = lim.rlim_cur / 2;
	        if (nb_partitions >= max_open_files)
        	    nb_passes++;
	        else
        	    break;
	}while(1);
    	printf("no of partitions before %lu and after %d passes %lu \n",nb_partitions*nb_passes,temp_partition,nb_passes);
    uint64_t total_IO =   volume * 2LL * 1024LL*1024LL   ;// in bytes  +   nb_passes * ( volume / (sizeof(kmer_type)*4) )    ; // in bytes
    uint64_t temp_IO = 0;
    BinaryBankConcurrent * redundant_partitions_file[nb_partitions]; 
    char redundant_filename[nb_partitions][256];
    kmer_type kmer;
    int max_read_length = KMERSBUFFER_MAX_READLEN;
    kmer_type * kmer_table_seq = (kmer_type * ) malloc(sizeof(kmer_type)*max_read_length); ;
    kmer_type * kmer_length_table_seq = (kmer_type * ) malloc(sizeof(kmer_type)*max_read_length);

    BinaryReads *  binread = NULL;
    if(use_compressed_reads)
        binread = new BinaryReads(return_file_name(binary_read_file),true);

    fprintf(stderr,"Sequentially counting ~%llu MB of kmers with %d partition(s) and %d passes using %d thread(s), ~%d MB of memory and ~%d MB of disk space\n", (unsigned long long)volume, nb_partitions,nb_passes, nb_threads, max_memory * nb_threads, max_disk_space);

    STARTWALL(count);

    mkdir(temp_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    // Open totalKmers files to store counts of totalKmers different k's
    BinaryBankConcurrent * SolidKmers[totalKmers];
    for (int s=0;s<totalKmers;s++) 
    {	
		char temp[1024];
		sprintf(temp,"%s.%d",return_file_name(solid_kmers_file),Kmerlist[s]);
		uint64_t exp = (((uint64_t)1)<<(Kmerlist[s]*2))-1;
		SolidKmers[s] = new BinaryBankConcurrent(temp,sizeof(kmer),true,nb_threads);
		//printf("kmer is %d exp is %llu \n",Kmerlist[s],exp);
		//BinaryBankConcurrent * SolidKmers = new BinaryBankConcurrent(return_file_name(solid_kmers_file),sizeof(kmer),true,nb_threads);

	    if (write_count)
	    {
        	// write k-mer nbits as the first 4 bytes; and actual k-mer size as the next 4 bits
      		  uint32_t kmer_nbits = sizeof(kmer) * 8;
   	     	SolidKmers[s]->write_buffered(&kmer_nbits, 4,0);
        	SolidKmers[s]->write_buffered(&Kmerlist[s], 4,0);
        	SolidKmers[s]->flush(0);
        }
   }

    int64_t estimated_NbReads = Sequences->estimate_nb_reads();
    char * rseq;
    int readlen;
    int64_t NbSolid = 0;
    int64_t * NbSolid_omp = (int64_t  *) calloc(nb_threads,sizeof(int64_t));
    //long total_kmers_per_partition[nb_partitions]; //guillaume probably commented it because updating this variable would require synchronization
    long distinct_kmers_per_partition[nb_partitions];
    uint64_t  * histo_count = (uint64_t  *) calloc(10001,sizeof(uint64_t));


#if OMP
    uint64_t  **  histo_count_omp = (uint64_t  **) calloc(nb_threads,sizeof(uint64_t *));
    for(int ii=0;ii<nb_threads;ii++)
    {
        histo_count_omp[ii]= (uint64_t  *) calloc(10001,sizeof(uint64_t));
    }
#endif
    

    
   
    //start by the conversion of the file to binary format

    if(use_compressed_reads)
    {
        char * pt_begin;
        int idx =0 ;
        int64_t NbRead = 0;
        Progress progress_conversion;
       // progress_conversion.timer_mode=1; // to switch to timer mode (show elapsed and estimated remaining time)
        progress_conversion.init(estimated_NbReads,"First step: Converting input file into Binary format");
        
        Sequences->rewind_all();
        while(1)
        {
            if(! Sequences->get_next_seq(&rseq,&readlen)) break; // read  original fasta file
            if(readlen > max_read_length) // realloc kmer_table_seq if needed
            {
                max_read_length = 2*readlen;
                kmer_table_seq = (kmer_type * ) realloc(kmer_table_seq,sizeof(kmer_type)*max_read_length);
            	kmer_length_table_seq = (kmer_type * ) realloc(kmer_length_table_seq,sizeof(kmer_type)*max_read_length);
	    }
            
            pt_begin = rseq;
            //should be ok
            while (pt_begin < (rseq+ readlen))
            {
                idx=0; // start a new read

                //skips NN
                while (*pt_begin =='N' && pt_begin < (rseq+ readlen))
                {
                    pt_begin ++;
                }
                // goes to next N or end of seq
                while ( (pt_begin[idx] !='N') &&  ((pt_begin +idx) < (rseq+ readlen))  )
                {
                    idx++;
                }
                
                //we have a seq beginning at  pt_begin of size idx  ,without any N, will be treated as a read:
                binread->write_read(pt_begin,idx);
		revcomp_sequence(pt_begin,idx); // reverse complement the string 
		binread->write_read(pt_begin,idx); // write reverse complement string 
		revcomp_sequence(pt_begin,idx); // restore the string 

		pt_begin += idx;
            }
            
            // binread->write_read(rseq,readlen);
            
            
            NbRead++;
            if ((NbRead%10000)==0)
            {
                progress_conversion.inc(10000);
            }
        }
	//printf("Number of reads converted to binary %d \n",NbRead);
        progress_conversion.finish();
        binread->close();

    }
    ///fin conversion
    if (clear_cache)
    {
#ifdef OSX
        system("purge");
#else
        system("echo 3 > /proc/sys/vm/drop_caches");
#endif
    }
    
    
    
#if SINGLE_BAR
    Progress progress;
    char message[1000];
    sprintf(message,"Counting kmers");
    progress.timer_mode=1;
    if (verbose == 0 )
        progress.init(total_IO,message);
#endif
    
    //use_compressed_reads=false; // for testing compute_kmer_from_one_seq 
    // how many times we will traverse the whole reads file (has an influence on temp disk space)

   uint64_t iter_partition=0;
    for (uint32_t current_pass = 0; current_pass < nb_passes; current_pass ++)
    {
	// stop computing if all partitions are done Added by Raunaq
        if (iter_partition==temp_partition)
		break;
	if(use_compressed_reads ) //open binary reads for reading
            binread->open(false);
        
        STARTWALL(debpass);
        STARTWALL(debw);
	int initial_value = current_pass*nb_partitions;
        for (uint32_t p=0;p<nb_partitions;p++)
        {
            sprintf(redundant_filename[p],"%s/partition%d.redundant_kmers",temp_dir,p);
            redundant_partitions_file[p] =  new BinaryBankConcurrent (redundant_filename[p],sizeof(kmer_type),true, nb_threads);
            distinct_kmers_per_partition[p]=0;
       	}
	int final_value = ((current_pass+1)*nb_partitions)-1;
	printf("Storing k-mers in partition files between %d and %d \n",initial_value,final_value);
        Sequences->rewind_all();
#if !SINGLE_BAR
        Progress progress;
        progress.timer_mode=1; // to switch to timer mode (show elapsed and estimated remaining time)
        char message[1000];
        sprintf(message,"Pass %d/%d, Step 1: partitioning",current_pass+1,nb_passes);
        if (verbose == 0 )
            progress.init(estimated_NbReads,message);
#endif
     

        
        //current_pass> 0 &&
#if OMP
#pragma omp parallel if(use_compressed_reads)  num_threads(nb_threads)
#endif
        {
            int64_t  nbkmers_written =0;
            int tid =0;
            int64_t NbRead = 0;
            int64_t nread =0;
            int64_t tempread =0;
	    long it_zero_wrt =0;
#if OMP

            tid = omp_get_thread_num();
#endif
            int nreads_in_buffer= 1000;
            KmersBuffer * kbuff =NULL;
            if(use_compressed_reads)
            {
                kbuff = new KmersBuffer (binread, 1000000,  nreads_in_buffer); //buffer size (in nb of kmers), seq per task // the buffer is per thread
                kbuff->binary_read_file = binread->binary_read_file;
            }

            kmer_type * kmer_table ;
            kmer_type * kmer_length_info ; // Added by Raunaq, to store the length of read into the partitions file
	    while(1)
            {

                //read the fasta file
                if(use_compressed_reads) // && current_pass>0
                {
                    nread = kbuff->readkmers();
                    if( nread < 0) break;
                    NbRead+= nread;
                    tempread+= nread;
                }
                else
                {
                    if(! Sequences->get_next_seq(&rseq,&readlen)) break; // read  original fasta file
                    if(readlen > max_read_length) // realloc kmer_table_seq if needed
                    {
                        max_read_length = 2*readlen;
                        kmer_table_seq = (kmer_type * ) realloc(kmer_table_seq,sizeof(kmer_type)*max_read_length);
            		kmer_length_table_seq = (kmer_type * ) realloc(kmer_length_table_seq,sizeof(kmer_type)*max_read_length);
                    }

                }

//                if(use_compressed_reads ) //write compressed read file at first pass //&& current_pass==0
//                    binread->write_read(rseq,readlen);

                int i;
                int nbkmers =readlen-sizeKmer+1;  

                if( use_compressed_reads) //current_pass >0 &&
                {
                    nbkmers = kbuff->nkmers;
                    kmer_table = kbuff->kmers_buffer;
		    kmer_length_info = kbuff->kmer_length;
                } 
                else //old fashion   
                {
                    compute_kmer_table_from_one_seq(readlen,rseq,kmer_table_seq,kmer_length_table_seq,Kmerlist[totalKmers-1]); // Added by Raunaq for computing kmers for all values of k 
                    nbkmers =readlen-Kmerlist[totalKmers-1]+1;  
                    kmer_table = kmer_table_seq;
		    kmer_length_info = kmer_length_table_seq;
                    NbRead++;
                    //printf("Number of kmers read from seq %d \n",nbkmers);
		}
		
                nbkmers_written= 0;
		char  temp_kmer[256];
		int zero;
                //compute the kmers stored in the buffer kmer_table
                for (i=0; i<nbkmers; i++)
                {
                    kmer_type lkmer;
					kmer_type lkmer_length;
                    // kmer = extractKmerFromRead(rseq,i,&graine,&graine_revcomp);

                    lkmer = kmer_table[i];
					lkmer_length = kmer_length_info[i];
		   // zero = code2seq(lkmer,temp_kmer);
					long pass_lkmer = code2first_n_nucleotide(lkmer,size_for_reestimation);
					unordered_map<long,int>::const_iterator got = part_hash.find(pass_lkmer);
					int p;// compute in which partition this kmer falls into
					if(got==part_hash.end())
						continue;
					else
						p = got->second; 
                    // check if this kmer should be included in the current pass
                    if(!(p >= initial_value && p<= final_value))
						continue;


/*		
#ifdef _ttmath
                    (reduced_kmer % nb_partitions).ToInt(p);
#else
                    p = reduced_kmer % nb_partitions;
#endif
*/
					p = p - current_pass*nb_partitions;  
                    nbkmers_written++;

                    redundant_partitions_file[p]->write_element_buffered(&lkmer,tid); // save this kmer to the right partition file
					redundant_partitions_file[p]->write_buffered(&lkmer_length,sizeof(lkmer_length),tid,false); // save the kmer length next to the kmer in the same partition file
		    // total_kmers_per_partition[p]++; // guillaume probably commented it because updating this variable would require synchronization

                }
                //NbRead++;
#if SINGLE_BAR
                if(verbose==0)
                {
                if (nb_threads == 1)
                    progress.inc(nbkmers_written * sizeof(kmer_type));
                else
                    progress.inc(nbkmers_written * sizeof(kmer_type),tid);
                }
#endif
             //   if ((NbRead%10000)==0)
                if(tempread> 10000)
                {
                    tempread -= 10000;
                    if (verbose)
                        fprintf (stderr,"%cPass %d/%d, loop through reads to separate (redundant) kmers into partitions, processed %lluM reads out of %lluM",13,current_pass+1,nb_passes,(unsigned long long)(NbRead/1000/1000),(unsigned long long)(estimated_NbReads/1000/1000));
#if !SINGLE_BAR
                    else
                        if (nb_threads == 1)
                            progress.set(NbRead);
                        else
                            progress.inc(10000,tid);
#endif
                }
            } //end while
           // printf("Count of zero in write is %lu \n",it_zero_wrt);
            if(use_compressed_reads)
                delete kbuff;
        } // end OMP 


        
#if !SINGLE_BAR
        if (verbose == 0)
        {
            if (nb_threads == 1)
             progress.finish();
            else
              progress.finish_threaded();  // here only one thread
            
            sprintf(message,"Pass %d/%d, Step 2: computing kmer count per partition",current_pass+1,nb_passes);
            progress.init(nb_partitions+1,message);
        }
#endif
        
        if (verbose)fprintf(stderr,"\n");

        if (verbose >= 2)
        {
            STOPWALL(debw,"Writing redundant kmers");
        }
        STARTWALL(debtri);
	


            for (uint32_t p=0;p<nb_partitions;p++)
            {	
				redundant_partitions_file[p]->close();
                redundant_partitions_file[p]->open(false);
            }



        // for better timing: clear the file cache, since the partitions may still be in memory, that's unfair to low mem machines
        if (clear_cache)
        {
#ifdef OSX
            system("purge");
#else
            system("echo 3 > /proc/sys/vm/drop_caches");
#endif
        }

        //quick and dirty parall with omp, testing
        //todo if we want omp and histo : separate histo_count tab per thread that needs to be merged at the end
        // TODO to guillaume: remove that todo above, because it is done, right?
        kmer_type lkmer,lkmer_length,lkmer_temp,exp;
	long it_zero=0;
	OAHash * hash;
	int p,s;
#if OMP 
        //omp_set_numthreads(2);  //num_threads(2) //if(!output_histo) num_threads(nb_threads)
#pragma omp parallel for private (p,s,lkmer,lkmer_length,hash,lkmer_temp,exp)  num_threads(nb_threads)
#endif        
        // load, sort each partition to output solid kmers
        for ( p=0;p<nb_partitions;p++)
        {
			char temp_kmer[256];  // bug check code 
			int zero;
			kmer_type lkmer_revcomp; // to store revcomps
				
           	bool use_hashing_for_this_partition = use_hashing;
			if(hybrid_mode)
			{
				if(   (redundant_partitions_file[p]->nb_elements()*sizeof(kmer_type)) <  (max_memory*1024LL*1024LL) )  // Maintain totalKmers hash for each partition file
				{	
					use_hashing_for_this_partition = false;
				}
				else
				{
					use_hashing_for_this_partition = true;
				}
			}
            int tid =0;
			//int s;
			//Computing if hashing should be used or not for this partition
#if OMP
            tid = omp_get_thread_num();
#endif
            //use_hashing_for_this_partition = false;  //to check the vector part of the code
           	if (use_hashing_for_this_partition)
            {
                // hash partition and save to solid file
				
				hash = new OAHash(max_memory*1024LL*1024LL/2); // One hash to store all types of k-mer lengths

				uint64_t nkmers_read=0;
				redundant_partitions_file[p]->read_element_buffered(&lkmer_length);

				while (redundant_partitions_file[p]->read_element_buffered(&lkmer))
				{
			
					if(lkmer_length == Kmerlist[0])  //only add the largest k-mer 
						hash->increment(lkmer,convert_to_int(lkmer_length));
					else
					{
						unordered_map<int,int>::const_iterator got = kmerlength_map.find(convert_to_int(lkmer_length));
						exp = (((kmer_type)1)<<(got->second*2))-1;
						lkmer_temp = lkmer & exp;
						hash->increment(lkmer_temp,got->second);

					}	
					if(!redundant_partitions_file[p]->read_element_buffered(&lkmer_length)) 
					{
						break;
					}
					nkmers_read++;
#if SINGLE_BAR
					if(verbose==0 && nkmers_read==10000)
					{
						if (nb_threads == 1)
							progress.inc(nkmers_read*sizeof(kmer_type));
						else
							progress.inc(nkmers_read*sizeof(kmer_type),tid);
						nkmers_read=0;
					}
#endif
                }
                
                
				if (verbose >= 2)
					 printf("Pass %d/%d partition %d/%d hash load factor: %0.3f\n",current_pass+1,nb_passes,p+1,nb_partitions,hash->load_factor());
                	for( s=0;s<totalKmers;s++) 
					{
						OAHash * temp_ = new OAHash(max_memory*1024LL*1024LL/2);
						hash->start_iterator();
						while (hash->next_iterator())
                				{
							uint_abundance_t abundance = hash->iterator->value;
        	       		 			uint_abundance_t abund_tid = (current_pass+1)*100+p;
							if(output_histo)
							{
							 uint_abundance_t saturated_abundance;
							 saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
#if OMP
							 histo_count_omp[tid][saturated_abundance]++;
#else
					
							 histo_count[saturated_abundance]++;
#endif
							}
							int length_kmer = hash->iterator->length;
							lkmer = hash->iterator->key;
	                    				if (abundance >= nks && abundance <= max_couv && length_kmer == Kmerlist[s])
							{
								//write if lkmer is the smaller of it and its reverse complement
								lkmer_revcomp = revcomp(lkmer,length_kmer);
								if(lkmer < lkmer_revcomp)
								{
								SolidKmers[s]->write_element_buffered(&(hash->iterator->key),tid);
							
								 NbSolid_omp[tid]++;
								if (write_count)
										SolidKmers[s]->write_buffered(&abundance, sizeof(abundance),tid, false);
								}
							}
		                    			distinct_kmers_per_partition[p]++;
							if(s!=totalKmers-1)
							{
								if(length_kmer == Kmerlist[s])
								{
									exp = (((kmer_type)1)<<(Kmerlist[s+1]*2))-1;
									lkmer_temp = lkmer & exp;
									temp_->increment_by_value(lkmer_temp,abundance,Kmerlist[s+1]);
								}else {
									temp_->increment_by_value(lkmer,abundance,length_kmer);
								}
							}
						}
						hash->~OAHash();
						hash = temp_;
					}
				hash->~OAHash();
			//printf("All hashes closed and destroyed \n");
			}
            
			else
			{
				// This part does it in slower fashion
				// sort partition and save to solid file 
        	    //vector < kmer_type > kmers;
				vector < kmer_type > kmers[totalKmers];
                uint64_t nkmers_read=0;
               	//int s=0; 
				
				redundant_partitions_file[p]->read_element_buffered(&lkmer_length);
				while (redundant_partitions_file[p]->read_element_buffered (&lkmer))
				{
    		        for(s=0;s<totalKmers;s++)
					{
						//kmer_type lkmer_temp;
						//kmer_type exp;
						if(lkmer_length<Kmerlist[s])
							continue;
						if(s==0)
							kmers[s].push_back (lkmer);
						else
						{
							exp = (((kmer_type)1)<<(Kmerlist[s]*2))-1;
							lkmer_temp = lkmer & exp; // Converting the kmer to its smaller equivalent in binary 
							kmers[s].push_back (lkmer_temp);
						}
                    }
					nkmers_read++;
					if(!redundant_partitions_file[p]->read_element_buffered(&lkmer_length)) break;  //Added to get the next length of kmer
#if SINGLE_BAR
					if(verbose==0 && nkmers_read==10000)
					{
						if (nb_threads == 1)
							progress.inc(nkmers_read*sizeof(kmer_type));
						else
							progress.inc(nkmers_read*sizeof(kmer_type),tid);
						nkmers_read=0;
					}
#endif
                }
                
                for(s=0;s<totalKmers;s++)
               	{
					sort (kmers[s].begin (), kmers[s].end ());
                
					kmer_type previous_kmer = *(kmers[s].begin ());
					uint_abundance_t abundance = 0;
					for (vector < kmer_type >::iterator it = kmers[s].begin (); it != kmers[s].end ();it++)
					{
						kmer_type current_kmer = *it;
					
						if (current_kmer == previous_kmer)
							abundance++;
						else
						{
							if(output_histo)
							{
									uint_abundance_t saturated_abundance;
									saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
#if OMP
									histo_count_omp[tid][saturated_abundance]++;
#else
									histo_count[saturated_abundance]++;
#endif
					
							}
							if (abundance >= nks  && abundance <= max_couv)
							{
								 NbSolid_omp[tid]++;
								 SolidKmers[s]->write_element_buffered(&previous_kmer,tid);
						
								 if (write_count)
									SolidKmers[s]->write_buffered(&abundance, sizeof(abundance),tid, false);
							}		
								abundance = 1;
							distinct_kmers_per_partition[p]++;
						}
						previous_kmer = current_kmer;
					}
                
                //last kmer
					distinct_kmers_per_partition[p]++;
					if(output_histo)
					{
							uint_abundance_t saturated_abundance;
							saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
#if OMP
							histo_count_omp[tid][saturated_abundance]++;
#else
							histo_count[saturated_abundance]++;
#endif
				
					}
					if (abundance >= nks && abundance <= max_couv)
					{
							NbSolid_omp[tid]++;
							SolidKmers[s]->write_element_buffered(&previous_kmer,tid);
				
							if (write_count)
							   SolidKmers[s]->write_buffered(&abundance, sizeof(abundance),tid, false);
				
					}	
				}
			}
            
            
		//printf("Done writing kmers for all K \n");
            
            	if (verbose >= 1)
                	fprintf(stderr,"%cPass %d/%d, loaded and sorted partition %d/%d, found %lld solid kmers so far",13,current_pass+1,nb_passes,p+1,nb_partitions,(long long)(NbSolid_omp[tid]));
            
		//printf("Done writing kmers for all K %d check 1 \n",p);
            	if (verbose >= 2)
                	printf("\nPass %d/%d partition %d/%d %ld distinct kmers\n",current_pass+1,nb_passes,p+1,nb_partitions,/*total_kmers_per_partition[p],*/distinct_kmers_per_partition[p]);
            
#if !SINGLE_BAR
            	if (verbose == 0 && nb_threads==1)
                	progress.inc(1);
            	else if (verbose == 0 && nb_threads>1)
                	progress.inc(1,tid);
#endif
            
            	//if(redundant_partitions_file[p]->find_error()) {
		//	printf("Error in the binary file \n");
		//}
            	redundant_partitions_file[p]->close();
		
            	remove(redundant_filename[p]);
 
        } // end for partitions

#if OMP
        //merge histo
        if(output_histo)

        {
            for (int cc=1; cc<10001; cc++) {
                uint64_t sum_omp = 0;
                for(int ii=0;ii<nb_threads;ii++)
                {
                    sum_omp += histo_count_omp[ii][cc];
                }
                histo_count[cc] = sum_omp;
            }
        }
#endif
        
#if !SINGLE_BAR
        if (verbose == 0 && nb_threads == 1)
            progress.finish();
        else if (verbose == 0 && nb_threads > 1 )
            progress.finish_threaded();
#endif

        if (verbose) fprintf(stderr,"\n");

        if (verbose >= 2)
        {
            STOPWALL(debtri,"Reading and sorting partitions");
            STOPWALL(debpass,"Pass total");

        }
       
	//printf("Done writing kmers for all K check 4 \n");
        if(use_compressed_reads)
            binread->close();
        
        //delete
            for (uint32_t p=0;p<nb_partitions;p++)
            {
                delete redundant_partitions_file[p] ;
            }
        
    }
	//printf("Done writing kmers for all K check 5 \n");

    //single bar
#if SINGLE_BAR
    if (verbose == 0 && nb_threads == 1)
        progress.finish();
    else if (verbose == 0 && nb_threads > 1 )
        progress.finish_threaded();
#endif
    
    if(output_histo)
    {
        FILE * histo_file = fopen(return_file_name(histo_file_name),"w");
        for (int cc=1; cc<10001; cc++) {
            fprintf(histo_file,"%i\t%llu\n",cc,(unsigned long long)(histo_count[cc]));
        }
        fclose(histo_file);
    }
    free(histo_count);

    NbSolid = NbSolid_omp[0];
#if OMP
    NbSolid=0;
    for(int ii=0;ii<nb_threads;ii++)
    {
        NbSolid += NbSolid_omp[ii];
    }
#endif
   for ( int s=0;s<totalKmers;s++) 
    	SolidKmers[s]->close();
    printf("\nSaved %lld solid kmers\n",(long long)NbSolid);
    rmdir(temp_dir);

    STOPWALL(count,"Counted kmers");
    fprintf(stderr,"\n------------------ Counted kmers and kept those with abundance >=%i,     \n",nks);
} 




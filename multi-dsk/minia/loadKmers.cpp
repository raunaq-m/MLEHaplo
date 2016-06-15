#include "loadKmers.h"
int totalKmers;
int smallestKmer ;
int *Kmerlist;
unordered_map<int,int> kmerlength_map;
int* loadKmers(char *kmerfname) {

	//Load the file containing list of k-mers (assuming kmers entered in decreasing order) 
	string line;	
	int count=0,counter=0;
	ifstream infile(kmerfname);
	if(!infile.is_open()) {
		fprintf(stderr, "%s file doesn't exist or failed to load. Please check if kmer list file exists \n",kmerfname);
	}
	while(!infile.eof()){
		getline(infile,line);
		count++;
	}
	//fprintf(stderr,"No of kmers asked for computation  %d\n",count-1);
	infile.close();
	totalKmers = count-1;
 	Kmerlist = new int [count]; 
	infile.open(kmerfname);
	while(!infile.eof()){
		getline(infile,line);
		Kmerlist[counter++] =atoi(line.c_str()); 
	}
	int largestKmer = Kmerlist[0];
	smallestKmer = Kmerlist[totalKmers-1];
	//fprintf(stderr,"Largest %d, Smallest %d\n",Kmerlist[0],Kmerlist[totalKmers-1]);
	// Map numbers from Kmerlist[totalKmers-1] to Kmerlist[0] to the 
	int iter = 0;
	for( counter=Kmerlist[0];counter>=Kmerlist[totalKmers-1];counter--)
	{
		if (counter < Kmerlist[iter]) 
			iter++; 
		//printf("Mapped %d to %d int \n",counter, Kmerlist[iter]);
		pair<int,int> temp_pair(counter,Kmerlist[iter]);
		kmerlength_map.insert(temp_pair);
	}
	/*for (counter=0;counter<count-1;counter++) {
		fprintf(stderr,"%d kmer value\n",Kmerlist[counter]);
	}*/
	
	return Kmerlist;
}

int reestimate_partitions(int size_of_lmers,uint64_t partition_volume,double * lkmer_counts, long * hash_vals, int * partition_files,int fact)
{
	//Sort the counts using merge sort, and also maintain the positions of sorted numbers
	vector<clmer> sorted_lmers;
	long total_bins = pow(4,size_of_lmers);
	// make pairs of lmercount and lmer 
	for (int i=0;i<total_bins;i++)
		sorted_lmers.push_back(make_pair(i,lkmer_counts[i]*fact)); //Factor 2 is for reverse complements//REVERSECOMPLEMENT 
	// sort the counts of kmers and maintain their index 
	sort(sorted_lmers.begin(),sorted_lmers.end(),comparator);
	double current_vol=0; int part = 0, iter=0, items_in_cur = 0;
	FILE *cTrack = fopen("write_counts.txt","w");
	double min = 100;
	for(vector<clmer>::iterator it=sorted_lmers.begin(); it!=sorted_lmers.end();++it)
	{
		clmer temp = *it;
		//if(current_vol>=partition_volume/2 || items_in_cur > 5000 )
		// if(current_vol>=partition_volume )
		// {
			// current_vol = 0;
			// items_in_cur = 1;
			// part++;
			// hash_vals[iter] = temp.first;
			// partition_files[iter]=part;
		// }
		// else
		// {
		if(temp.second !=0)
			if(temp.second < min)
				min = temp.second;
		if(temp.second ==0)
			temp.second = min;

		// Make sure that a single file doesn't gets too many small k-mer counts
		current_vol +=temp.second;
		if(current_vol <= partition_volume)
		{	
			items_in_cur++;
			hash_vals[iter] = temp.first;
			partition_files[iter]=part;
		}else
		{
			if(items_in_cur == 0) 
			{
				current_vol = 0;
			}else
			{
				current_vol = temp.second;
				part++;
			}
			items_in_cur = 1;
			hash_vals[iter] = temp.first;
			partition_files[iter] = part;
		}
		//}
		fprintf(cTrack,"iterr %d kmer count %f,current vol %f, kmer %lu, partiti    on file %d\n",iter, temp.second,current_vol,hash_vals[iter],partition_files[iter]);
		//printf("iterr %d kmer count %f,current vol %f, kmer %lu, partition file %d\n",iter, temp.second,current_vol,hash_vals[iter],partition_files[iter]);
		iter++;
	}
	// print out lmers their counts and the partition files 
	iter=0;
/*	sort(sorted_lmers.begin(),sorted_lmers.end(),comparator);
	for(vector<clmer>::iterator it=sorted_lmers.begin(); it!=sorted_lmers.end();++it)
	{
		clmer temp= *it;
		printf("kmer count %lu, kmer %lu, partition file %d\n",temp.second,hash_vals[iter],partition_files[iter]);
		iter++;
	}
*/
	fclose(cTrack);
	return part+1;
}

bool comparator(const clmer& l,const clmer& r)
{
	return l.second>r.second;
}

/*int main(int argc, char *argv[]){
	
	int i=10;
	long a1[10],b[10];
	int c[10];
	a1[0]=10;a1[1]=30;a1[2]=50;a1[3]=40; a1[4]=90;a1[5]=70;a1[6]=20;a1[7]=60;a1[8]=80;a1[9]=100;
	for(int j=0;j<i;j++)
		printf("%lu\n",a1[j]);
	Sort(10,a1,b,c);
	for(int j=0;j<i;j++)
		printf("%lu %lu %d\n",a1[j],b[j],c[j]);
	return 0;
}*/

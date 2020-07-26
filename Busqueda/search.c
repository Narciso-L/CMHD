/*Copyright (C) <2016>  <Narciso López-López>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "search.h"



/*loads vocabulary, louds & dimension_size from all dimensions, 
moreover loads the number of levels, total elements per level and elements per level from array_dim_result (all merged dimensions) */
void load_dims(dim **array_dims, uint numDims,char *path){

	uint i = 0, j = 0; 
	char dest[500];						//path's file 
	strcpy(dest, path);
	strcat(dest, "/dims.dat");			//filename
	FILE *dims_file = fopen(dest,"rb");	//open file for reading

	if(dims_file == NULL){				//fiel not found
		perror(dest);
		exit(-1);
	}

	for(i=0;i<numDims;i++){		
		array_dims[i] = (dim*)malloc(sizeof(dim));									//allocates memory for each dimension
		array_dims[i]->bitmap = (bitRankW32Int*)malloc(sizeof(bitRankW32Int));		//allocates memory for bitmap
		load(array_dims[i]->bitmap,dims_file);										//loads bitmap from file
	}

	for(i=0;i<numDims;i++){		
		fread(&array_dims[i]->dimension_size,sizeof(uint),1,dims_file);				//loads dimension_size from file
	}

	uint lenWord;
	while(j<numDims){
		array_dims[j]->vocabulary = (char**)malloc(sizeof(char*)*array_dims[j]->dimension_size);
		for(i=0;i<array_dims[j]->dimension_size;i++){
			fread(&lenWord,1,sizeof(uint),dims_file);								//loads word's size
			array_dims[j]->vocabulary[i] = (char*)malloc(sizeof(char)*(lenWord));
			fread(array_dims[j]->vocabulary[i],lenWord,sizeof(char),dims_file);		//loads vocabulary for each dimension
			//printf("voc: %s\n",array_dims[j]->vocabulary[i]);
		}
		j++;
	}

	//allocates memory for array_dim_result
	array_dims[numDims] = (dim*)malloc(sizeof(dim));

	fread(&array_dims[numDims]->levels,sizeof(uint),1,dims_file);					//loads number of levels
	array_dims[numDims]->total_elements_per_level = (uint*)malloc(sizeof(uint)*(array_dims[numDims]->levels+1));
	fread(array_dims[numDims]->total_elements_per_level, sizeof(uint), array_dims[numDims]->levels+1, dims_file); //loads array total elements

	array_dims[numDims]->elements_per_level = (uint**)malloc(sizeof(uint*)*(array_dims[numDims]->levels+1));	//memory for matrix elements per level

	uint e;
	uint numPos = 1;
	for(i=0;i<array_dims[numDims]->levels+1;i++){  
		array_dims[numDims]->elements_per_level[i] = (uint*)malloc(sizeof(uint)*(numPos));
		FTRep *arrayDac = loadFTFlex(dims_file);										//loads DAC
		for(e = 0; e<numPos; e++) {
			array_dims[numDims]->elements_per_level[i][e] = accessFT(arrayDac, e + 1);	//loads elements per level
		}
		free(arrayDac);
		numPos = array_dims[numDims]->total_elements_per_level[i]; 
		//fread(array_dims[numDims]->elements_per_level[i],array_dims[numDims]->total_elements_per_level[i],sizeof(uint),dims_file);
	}

	fclose(dims_file);															//closes file

	// for(i=1;i<array_dims[numDims]->levels+1;i++){
	// 	for(j=0;j<array_dims[numDims]->total_elements_per_level[i-1];j++){
	// 		printf("%d ",array_dims[numDims]->elements_per_level[i][j] );		//prints elements per level
	// 	}
	// 	printf("\n\n");
	// }

	// printf("array_dims[numDims]->levels %d\n", array_dims[numDims]->levels);	//prints levels

}


//loads array_sum or array_sum2 and bitrank1 & 2
void load_bitmaps(bitmap *bitmap_dims, char *path){

	char dest[500];			//path's file 
	
	bitmap_dims->bitRank1 = (bitRankW32Int*)malloc(sizeof(bitRankW32Int));	//allocates memory for bitrank1
	bitmap_dims->bitRank2 = (bitRankW32Int*)malloc(sizeof(bitRankW32Int));	//allocates memory for bitrank1

	strcpy(dest, path);
	strcat(dest, "/bitmaps.dat");								//filename
	FILE *bitmaps = fopen(dest, "r");							//open file for reading

	if(bitmaps == NULL){										//file not found
		perror(dest);
		exit(-1);
	}

	load(bitmap_dims->bitRank1,bitmaps);						//loads from file bitrank1
	load(bitmap_dims->bitRank2,bitmaps);						//loads from file bitrank2
	fread(&bitmap_dims->regenerate,sizeof(uint),1,bitmaps );	//loads the flag to indicate if you have to use array_sum or array_sum2

	fclose(bitmaps);											//closes file

	// printf("bitmap: %d bits\n", bitmap_dims->bitRank1->n);	//prints the bits from bitrank1
	// printf("bitmap: %d bits\n", bitmap_dims->bitRank2->n);	//prints the bits from bitrank2

	// printf("REGENERATE %d\n",bitmap_dims->regenerate );		//prints the flag regenerate
}



/*------------------------------------------------------------------
 Initilizes the structures used.
 ---------------------------------------------------------------- */
//initializes hash table
void init_hash(unsigned long sizeVoc){
	
	_memMgr = createMemoryManager();	// Creates a Memory Manager
    initialize_hash(sizeVoc);			// Initilization of data structures used by the hashTable

    positionInTH = (unsigned long*) malloc (sizeVoc * sizeof(unsigned long));
	zeroNode = 0;
}


//loads the vocabulary from all dimensions and assigns the positions of the hash table corresponding to LOUDS
void hash_table(dim **array_dims, uint numDims){

	uint total_size = 0;
	uint i = 0;

	for(i=0;i<numDims;i++){		
		total_size += array_dims[i]->dimension_size;	//total size of the dimensions' vocabulary
	}
//	printf("total_size %d\n",total_size );

	init_hash(total_size); 								//initializes compressor for up to N_value words    
	uint z = 0, pos = 0, size = 0, dim=0;

	for(dim = 0; dim<numDims; dim++){	
		for(pos = 0; pos<array_dims[dim]->dimension_size; pos++){
			size = strlen(array_dims[dim]->vocabulary[pos]);
			//looks for a word in the hash table and returns its position in the vocabulary
			z = search ((unsigned char *)array_dims[dim]->vocabulary[pos], size, &addrInTH );	

		    if (z==zeroNode) {
		    	//inserts a new word in a given position of the hashTable (position previously computed)
				insertElement ((unsigned char *) array_dims[dim]->vocabulary[pos], size, &addrInTH);
			//	printf("voc is %s\n",array_dims[dim]->vocabulary[pos] );
			//	printf("size is %d\n",size );
			//	printf("%lu\n",addrInTH );
				hash[addrInTH].size = 0;
				hash[addrInTH].len = size;										
				hash[addrInTH].posInBitVector = select1(array_dims[dim]->bitmap,pos+2);  //stores the word's position from LOUDS in hash table
				hash[addrInTH].dimension = dim+1;										 //stores the dimension
			//	printf("posInBitVector %lu\n",hash[addrInTH].posInBitVector );
				positionInTH[zeroNode] = addrInTH;
			//  printf("%lu\n", positionInTH[zeroNode]);
			//  printf("NumElem %lu\n",NumElem );
			//  printf("zeronode %lu\n", zeroNode);
				zeroNode++;
			}	
		}	
	}
}


//cleans the character \n of each line
void cleanCharacter(char *string){		//clean the newline character of a string that is passed as parameter
	
	char *character;

	character = strchr (string, '\n');	//the character in the string is searched
	if (character)
		*character = '\0';				//it is replaced by end of string
}


//reads the query input file and stores the tags in datawarehouse->tags
void read_query_input_file(char *argv[], data *datawarehouse, uint numQueries, uint numDims){

	char buffer[BUFFER_SIZE];	
	char *word;
	uint i = 0; 

	datawarehouse->tags = (char**)malloc(sizeof(char*)*(numQueries*(numDims)));  

	FILE *input = fopen(argv[3], "r");					//open for reading  

	if(input == NULL){									//file error
		perror(argv[3]);
		exit(1);
	}

	while(fgets(buffer, BUFFER_SIZE, input)!= NULL){	//read line
		word = strtok(buffer, ",");
		cleanCharacter(word);							//cleans \n from each line

		while(word!=NULL){
			datawarehouse->tags[i] = (char*)malloc(sizeof(char)*strlen(word));	//allocates memory for word's size
			strcpy(datawarehouse->tags[i],word);		//stores the tag
		//	printf("%s\n",word );
			word = strtok(NULL, ",");
			i++;
		}
	}

	// for(i=0;i<numQueries*numDims;i++){
	// 	printf("%s\n",datawarehouse->tags[i] );			//prints the tags
	// }

	fclose(input);	
}


//frees the formula's allocated memory
void free_formula(formel *formula, uint levels){

/*	uint i = 0;
	for(i=0;i<levels;i++){
		free(formula->ones[i]);
	}

	free(formula->ones);

	for(i=0;i<levels;i++){
		free(formula->result[i]);
	} */

	free(formula->result);
	free(formula->rank0s);
	free(formula->select0s);
	free(formula);
}


//frees the offset's allocated memory
void free_offset(data *datawarehouse, uint levels){

/*	uint i = 0;
	for(i=0;i<levels;i++){
		free(datawarehouse->offset[i]);
	}
*/
	free(datawarehouse->offset);
}


//frees the ancestor's allocated memory
void free_ancestor(data *datawarehouse, uint levels){

	uint i=0;
	for(i=0;i<levels;i++){
		free(datawarehouse->ancestor[i]);
	}

	free(datawarehouse->ancestor);
}


//frees the datawarehouse's allocated memory
void free_datawarehouse(data *datawarehouse, uint numQueries, uint numDims){

	uint i = 0;
	for(i=0;i<numQueries*(numDims);i++){	//frees the tags
		free(datawarehouse->tags[i]);
	}

	free(datawarehouse->tags);

	// for(i=0;i<levels;i++){
	// 	free(datawarehouse->ancestor[i]);
	// }

	// free(datawarehouse->ancestor);

	// for(i=0;i<levels;i++){
	// 	free(datawarehouse->offset[i]);
	// }

	// free(datawarehouse->offset);
	free(datawarehouse);
}


//frees array_sum or array_sum2 and bitrank1 & 2
void free_bitmap_dims(bitmap *bitmap_dims,uint regenerate){

	destroyBitRankW32Int(bitmap_dims->bitRank1);		//frees bitrank1
	destroyBitRankW32Int(bitmap_dims->bitRank2);		//frees bitrank2

	if(regenerate == 0 ) free(bitmap_dims->array_sum);	//frees array_sum 1 or 2
	else free(bitmap_dims->array_sum2);

	free(bitmap_dims);
}


//frees the dimensions' allocated memory
void free_dims(dim **array_dims, uint numDims){
	
	uint i = 0, j = 0; 	

	//frees array_dim_result (total elements & elements per level)
	free(array_dims[numDims]->total_elements_per_level);

	for(j=0;j<array_dims[numDims]->levels;j++){
		free(array_dims[numDims]->elements_per_level[j]);		
	}
	free(array_dims[numDims]->elements_per_level);

	free(array_dims[numDims]);

 	//frees dimensions' vocabulary and bitmaps 
	for(i=0;i<numDims;i++){

		for(j=0;j<array_dims[i]->dimension_size;j++){
			free(array_dims[i]->vocabulary[j]);		
		}
		free(array_dims[i]->vocabulary);

		destroyBitRankW32Int(array_dims[i]->bitmap);	//frees the dimension's bitmap 

		free(array_dims[i]);
	}

	free(array_dims);
}



/********************************************** MAIN BODY *****************************************************/

int main(int argc, char *argv[]){

	uint numDims = atoi(argv[1]);		//get number of dimensions from command line
	uint numQueries = atoi(argv[2]);	//get number of queries from command line

	if(argc!= 5){
		printf("Input format: %s (filename) %d (N numDims) %d (Q numQueries) %s (query input file) %s (folder with binaries)  \n", argv[0],numDims,numQueries,argv[3],argv[4] );
		return -1;
	}

	if(numDims<2){
		printf("The number of dimensions has to be greater than 1\n");
		return -1;
	}

	
//	printf("\n\n******************* IMPORT DATA ************************\n");

	
	// LOAD DIMS
	dim **array_dims = (dim**)malloc(sizeof(dim*)*numDims+1);

	/*loads vocabulary, louds & dimension_size from all dimensions, 
	moreover loads the number of levels, total elements per level and elements per level from array_dim_result (all merged dimensions) */
	load_dims(array_dims,numDims,argv[argc - 1]);

	//WAREHOUSE
	data *datawarehouse = (data*)malloc(sizeof(data));

	//LOAD BITRANKS 1 & 2
	bitmap *bitmap_dims = (bitmap*)malloc(sizeof(bitmap));
	load_bitmaps(bitmap_dims,argv[argc - 1]);				//loads array_sum or array_sum2 and bitrank1 & 2

	FTRep *array_sumDac; 
	char dest[500];			//path's file 

	if(bitmap_dims->regenerate == 0){	
		strcpy(dest, argv[argc - 1]);
		strcat(dest, "/array_sum.dat");		//filename
		array_sumDac = loadFT(dest);		// loads ARRAY_SUM from DAC
	}
	else{
		strcpy(dest, argv[argc - 1]);
		strcat(dest, "/array_sum2.dat");	//filename
		array_sumDac = loadFT(dest);		// loads ARRAY_SUM 2 from DAC
	}

	

//	printf("\n\n/************************* HASH TABLE **************************/\n");

	//loads the vocabulary from all dimensions and assigns the positions of the hash table corresponding to LOUDS
	hash_table(array_dims,numDims);		



//	printf("\n\n/**************** READ INPUT FILE WITH TAGS AND VALUES ********************/\n");

	//reads the query input file and stores the tags in datawarehouse->tags
	read_query_input_file(argv, datawarehouse, numQueries, numDims);


//	printf("\n\n/********************** SEARCH *****************************\n");

	clock_t begin, end;
	double time_spent = 0.0;

	uint queries = 0;
	unsigned long total_cube = 0;			//stores the whole sum/final result

	while(queries<(numQueries*(numDims))){
	
		unsigned long lookfor = 0; 
		uint size = 0;
		uint i = 0, j = 0;
	
		uint levels = array_dims[numDims]->levels-1;	//assigns levels

		datawarehouse->ancestor = (uint**)malloc(sizeof(uint*)*(levels+1));
		datawarehouse->ancestor[levels] = (uint*)malloc(sizeof(uint)*(numDims));  
	
		uint *reps = (uint*)malloc(sizeof(uint)*numDims);
		uint times = 0;
		uint queries_times = 0;
	
		for(i=0;i<numDims;i++){ 	
			cleanCharacter(datawarehouse->tags[i+queries]);		//cleans \n
			size = strlen(datawarehouse->tags[i+queries]);		//gets the tag's size
		//  printf("size is %d\n",size );
		// 	printf("datawarehouse->tags[i] %s\n", datawarehouse->tags[i+queries]);
			//looks for a word in the hash table and returns its position in the vocabulary
			lookfor = search ((unsigned char *)datawarehouse->tags[i+queries], size, &addrInTH );
		//	printf("%lu\n",lookfor );
		//	printf("%s\n",hash[addrInTH].word );
		//	printf("posInBitVector %lu\n",hash[addrInTH].posInBitVector );

			if(lookfor==0)
				datawarehouse->ancestor[levels][i] = hash[addrInTH].posInBitVector;  //stores the position from LOUDS (last level)
			queries_times++;
		//	printf("queries_times are %d\n",queries_times );
		}
		queries += queries_times;

		begin = clock();						//init clock

		// printf("\n\n/*********************** LOOK FOR THE ANCESTOR **************************/\n");
		
		uint dim = 0;
		j = 0;

		//looks for the parent of each level in each dimension (from penultimate level)
		while(levels>0){
			datawarehouse->ancestor[levels-1]=(uint*)malloc(sizeof(uint)*(numDims));
			times = 1;
			for(j=0;j<numDims;j++){
				datawarehouse->ancestor[levels-1][j] = parent(array_dims[dim]->bitmap,datawarehouse->ancestor[levels][j]);
				dim++;
				if(dim == numDims) dim = 0;
			}
			levels--;	
		}


		// printf("\n\n/********************** REPETITIONS ************************\n");

		levels = array_dims[numDims]->levels;

		i=0; 
		for(j=0;j<numDims;j++){			//obtain the level for each dimension
			times=0;
			for(i=0;i<levels;i++){
				if(datawarehouse->ancestor[i][j]!=0) times++;
				// printf("datawarehouse->ancestor[%d][%d] %d\n",i,j,datawarehouse->ancestor[i][j] );
				// printf("times are %d\n",times );
			}
			reps[j] = times;	
		}

	    // printf("\n\n/****************** SAME LEVEL OR NOT ******************\n");

		uint same_level = 0;
		i=0;
		while(i<numDims-1){				//know if the query will be of the same level or not 
			if(reps[i]==reps[i+1]){
				same_level = 0;
			}
			else{ //different level
				same_level = 1;
			}
			if(same_level == 1) break;
			i++;
		}	

		// printf("same_level %d\n",same_level );

		// printf("\n\n/********************* LOOK FOR THE LOWEST LEVEL AMONG DIMENSIONS *******************\n");

		uint lowest_level = 0;
		i=0;
		while(i<numDims){
			if(reps[i]>=lowest_level){		//finds the lowest level among dimensions
				lowest_level = reps[i];
			}
			i++;
		}	

		// printf("lowest_level %d\n",lowest_level );


		/************************************************************************************************************/
		/************************************************************************************************************/
		/************************************************************************************************************/


		// printf("\n\n/************************** SAME LEVEL QUERIES *******************************/\n");
		if(same_level==0){
		
			levels = array_dims[numDims]->levels-1;

			//CALCULATES THE NUMBER OF ONES BETWEEN THE PREVIOUS AND NEXT 0 TO IMPLEMENT THE FORMULA

			/*
				rank0(pos) = x    *1*
				select0(x)+1 = y  *2*
				Rx = pos-y        *3*

				//Iterative
				(R1*D2)+R2 = X    *4*
				(X*D3)+R3  = X		where Dx is the number of ones in the next dimension  *5*
			*/

			formel *formula = (formel*)malloc(sizeof(formel));
			formula->rank0s = (uint*)malloc(sizeof(uint)*(numDims));	//number of rank0 before the position in louds
			formula->select0s = (uint*)malloc(sizeof(uint)*(numDims));	//select0 to obtain the position in louds 
			formula->ones = (uint**)malloc(sizeof(uint*)*(levels+1));	//number of ones between zeros		
			formula->result = (uint**)malloc(sizeof(uint*)*(levels+1)); //results to obtain the offsets	

			uint rank0 = 0, pos_prev0 = 0, pos_next0 = 0;
			
			dim = 0;
			uint jj = 0;
		
			while(levels>=array_dims[numDims]->levels-reps[0]){ 
				formula->result[levels]=(uint*)malloc(sizeof(uint)*numDims);
				formula->ones[levels]=(uint*)malloc(sizeof(uint)*numDims);
				for(jj=0; jj<numDims; jj++){
		
					formula->rank0s[jj] = datawarehouse->ancestor[levels][jj] - rank(array_dims[dim]->bitmap,datawarehouse->ancestor[levels][jj]-1); // *1*
					//printf("Rank 0's is %d\n",formula->rank0s[jj] );
					formula->select0s[jj] = select0(array_dims[dim]->bitmap,formula->rank0s[jj])+1;	// *2*
					//printf("Position is %d\n",formula->select0s[jj] );

					//CALCULATE THE NUMBER OF 1's BETWEEN PREVIOUS AND NEXT 0 TO USE THE FORMULA
					rank0 = datawarehouse->ancestor[levels][jj] - rank(array_dims[dim]->bitmap,datawarehouse->ancestor[levels][jj]-1);
					//printf("rank0 is %d\n",rank0);
					pos_prev0 = select0(array_dims[dim]->bitmap,rank0);
					//printf("pos_prev0 is %d\n",pos_prev0);
					pos_next0 = select0(array_dims[dim]->bitmap,rank0+1);
					//printf("pos_next0 is %d\n",pos_next0);
					formula->ones[levels][jj] = (pos_next0 - pos_prev0)-1;	//Dx
					formula->result[levels][jj] = datawarehouse->ancestor[levels][jj] - formula->select0s[jj];	// *3* Rx
					dim++;
					//printf("Result[%d][%d] is %d\n",levels,jj, formula->result[levels][jj] );
					//printf("ones[%d][%d] is %d\n",levels,jj, formula->ones[levels][jj]);
					if(dim == numDims) dim = 0;
				}

				levels--;
				if(levels == -1) break;
			}

			// NOW IS THE MOMENT TO CALCULATE THE OFFSET PER LEVEL IN THE BITMAP 2

			levels = array_dims[numDims]->levels-1;
			uint pos = 0 ; uint dims = 0;
			datawarehouse->offset = (uint**)malloc(sizeof(uint*)*(levels+1)); 

			while(levels>=array_dims[numDims]->levels-reps[0]){ 
	
				datawarehouse->offset[levels]=(uint*)malloc(sizeof(uint)*(numDims));
				for(jj=0; jj<numDims; jj=jj+(numDims)){
					//FOR 2 DIMENSIONS
					datawarehouse->offset[levels][pos] = (formula->result[levels][jj]*formula->ones[levels][jj+1])+formula->result[levels][jj+1];	// *4*
					// printf("Result[%d][%d] is %d\n",levels,jj, formula->result[levels][jj] );
					// printf("ones/desplaz[%d][%d] is %d\n",levels,jj+1, formula->ones[levels][jj+1]);
					// printf("OFFSET[%d][%d] is %d\n",levels,jj,datawarehouse->offset[levels][jj] );
					// printf("\n");
					for(dims = 2; dims < numDims; dims++){	//MORE THAN 2 DIMENSIONS
						// printf("ANCestor[%d][%d] %d\n",levels,jj, datawarehouse->ancestor[levels][jj]);
						// printf("Result[%d][%d] is %d\n",levels,jj+1, formula->result[levels][jj+1] );
						// printf("ones/desplaz[%d][%d] is %d\n",levels,jj+1, formula->ones[levels][jj+1]);
						// printf("OFFSET[%d][%d] is %d\n",levels,jj,datawarehouse->offset[levels][jj] );
						datawarehouse->offset[levels][pos] = (datawarehouse->offset[levels][pos]*formula->ones[levels][jj+dims])+formula->result[levels][jj+dims]; // *4* & *5*
					}
				//	printf("datawarehouse->offset[%d][%d] %d \n",levels,pos,datawarehouse->offset[levels][pos] );
					pos++;
				}

				pos = 0;
				levels--;
				if(levels == -1) break;
			}

		//	 printf("\n\n/************** GET THE SON OF EACH POSITION IN EVERY LEVEL ***************/\n\n");

			uint total = 0, position = 0, son = 0;

			levels = 1; dim = 1; jj=0;
			uint levels_max = array_dims[numDims]->levels;
			levels = array_dims[numDims]->levels-reps[0]; 
		
			if(bitmap_dims->regenerate == 0){			//ARRAY_SUM		
				
				if(reps[0] == 1){
					total_cube += accessFT(array_sumDac, datawarehouse->offset[levels][position]+1);
				 // printf("Query is %d\n",accessFT(array_sumDac, datawarehouse->offset[levels][position]+1)); 
				}
				else{
					while(levels<levels_max){
						total = datawarehouse->offset[levels][position] + son + 1;
						//printf("total %d\n",total );
						son = select0(bitmap_dims->bitRank2,total)+1;
						//printf("son is %d\n",son );

						levels++;
						if(levels == levels_max-1){	
							total = son + datawarehouse->offset[levels][position];
						//	printf("Level %u DAC -> %u (%u %u)\n",levels, total+1, son, datawarehouse->offset[levels][position]);
							total_cube += accessFT(array_sumDac, total+1);
						//	printf("Query is %d\n",accessFT(array_sumDac, total+1)); 
						}
					}
				}
				
			}
			else{		//ARRAY_SUM2
				total = 0, position = 0, son = 0;
				levels = array_dims[numDims]->levels-reps[0];
			
				if(reps[0] == 1){
					if(isBitSet(bitmap_dims->bitRank1,datawarehouse->offset[levels][position]) == 0){
					//	printf("Query is %d\n", 0);
					//	break;
					}
					else{
						total_cube += accessFT(array_sumDac, datawarehouse->offset[levels][position]+1);
					//	printf("Query is %d\n",accessFT(array_sumDac, datawarehouse->datawarehouse->offset[levels][position]+1));
					} 
				}

				else{
				
					while(levels<levels_max){	
					//	printf("datawarehouse->offset[levels][position] %d\n",datawarehouse->offset[levels][position] );
					//	printf("Levels -> %u || %u\n", levels, levels_max-1);
						if(levels == levels_max-1) {
							total = son + datawarehouse->offset[levels][position];
							total_cube += accessFT(array_sumDac, total+1);
						//	printf("Level %u DAC -> %u (%u %u)\n",levels, total+1, son, datawarehouse->offset[levels][position]);
						//	printf("Query is %d\n",accessFT(array_sumDac, total+1)); 
							break;
						}

						total = son + datawarehouse->offset[levels][position];
						//printf("Level %u position -> %u (%u)\n",levels, total, datawarehouse->offset[levels][position]);
						if(isBitSet(bitmap_dims->bitRank1,total) == 0) {
						//	printf("Query is bitset %d\n",isBitSet(bitmap_dims->bitRank1,total));
						//	printf("Query is %d\n",0);  
							break;
						}
						total = rank(bitmap_dims->bitRank1,total);
						son = select0(bitmap_dims->bitRank2,total)+1;
						//	printf("son %u || nze %u\n", son, total);
						levels++;

					}
				}
			}

			end = clock();
			time_spent += (double)(end - begin) / (CLOCKS_PER_SEC);
			
			//FREE FORMULA
			free_formula(formula, array_dims[numDims]->levels);
		
			//FREE OFFSET
			free_offset(datawarehouse,array_dims[numDims]->levels);

			//FREE ANCESTOR
			free_ancestor(datawarehouse, array_dims[numDims]->levels);

			free(reps);
		}

		else{ 
		//  printf("\n\n/************************** DISTINCT LEVEL QUERIES *******************************/\n");

		//	printf("\n\n/*************** LOOK FOR SONS AT EACH DIMENSION **************\n");

			son **sons = (son**)malloc(sizeof(son*)*numDims);

			uint x = 0;
			levels = array_dims[numDims]->levels;

			for(x=0;x<numDims;x++){
				sons[x] = (son*)malloc(sizeof(son));
				sons[x]->positions = (uint**)calloc(levels,sizeof(uint*));
			}

			uint total_sum = 0; 		//stores the sum for different level queries 
			uint dim = 0;
			
			while(dim<numDims){			//determine the children of a parent

				uint i = 0, j = 0, k = 0, m = 0;
			/*	printf("reps[dim] %d\n",reps[dim] );
				printf("lowest %d\n",lowest_level );
				printf("levels %d\n",levels );*/
				if(reps[dim]==lowest_level){
					sons[dim]->num_childs = 1;
					sons[dim]->positions[i] = (uint*)malloc(sizeof(uint)*SONS_PER_DIMENSION); 
					sons[dim]->positions[i][j] = datawarehouse->ancestor[levels-1][dim];   
				//	printf("sons[dim]->positions[i][j] %d\n",sons[dim]->positions[i][j] );
				}
				else{
					
					sons[dim]->num_childs = 1;
					sons[dim]->positions[i] = (uint*)calloc(sons[dim]->num_childs,sizeof(uint));
					sons[dim]->positions[i][j] = datawarehouse->ancestor[levels-1][dim];  
			
					for(i=1;i<lowest_level;i++){
						sons[dim]->positions[i] = (uint*)calloc(SONS_PER_DIMENSION,sizeof(uint));
					}
					
					i=0;
					sons[dim]->num_childs = 1;

					//printf("lowest_level-reps[dim] %d\n",lowest_level-reps[dim]);
					while(i<lowest_level-reps[dim]){
						//printf("sons[dim]->num_childs %d\n",sons[dim]->num_childs );
					
						m=0;
						sons[dim]->total_childs = 0;
						for(k=0;k<sons[dim]->num_childs;k++){
							// printf("i %d\n",i );
							// printf("k %d\n",k );
							sons[dim]->first = first_child(array_dims[dim]->bitmap, sons[dim]->positions[i][k]);
							sons[dim]->last = last_child(array_dims[dim]->bitmap, sons[dim]->positions[i][k]);
							sons[dim]->cant_childs = children(sons[dim]->first,sons[dim]->last);
							sons[dim]->total_childs += children(sons[dim]->first,sons[dim]->last);

							for(j=0;j<sons[dim]->cant_childs;j++){
								sons[dim]->positions[i+1][m] = next_child(sons[dim]->first,j);
							//  printf("sons[dim]->positions[i+1][m] %d\n",sons[dim]->positions[i+1][m] );
								m++;
							}
						}
						
						sons[dim]->num_childs = sons[dim]->total_childs;
						sons[dim]->cant_childs = 0;
						i++;
						//printf("\n");
					}
				}	
				dim++;
			}
			
			/*
			 for(i=0;i<numDims;i++){
				for(j=0;j<20;j++){
					//if(sons[i]->positions[levels-reps[i]][j]!=0)
					printf("%d ",sons[i]->positions[lowest_level-reps[i]][j] );
				}
				printf("\n\n");
			}
			*/

			//GET THE SONS TO COMBINE
			uint **combine = (uint**)malloc(sizeof(uint*)*numDims);
			uint *total_elements_per_dimension = (uint*)malloc(sizeof(uint)*numDims);
			uint cont = 0;

			for(i=0;i<numDims;i++){
				combine[i] = (uint*)malloc(sizeof(uint)*SONS_PER_DIMENSION);
				for(j=0;j<SONS_PER_DIMENSION;j++){
					if(sons[i]->positions[lowest_level-reps[i]][j]!=0)
					combine[i][j] = sons[i]->positions[lowest_level-reps[i]][j] ;
					cont++;
					if(sons[i]->positions[lowest_level-reps[i]][j+1]==0 ) break;
				}
				total_elements_per_dimension[i] = cont;
				cont=0;
			}
			/*
			for(i=0;i<numDims;i++){
				for(j=0;j<SONS_PER_DIMENSION;j++){
					if(combine[i][j]!=0)
						printf("%d ",combine[i][j] );
				}
				printf("\n\n");
			}*/
				
			
			// OBTAIN THE TOTAL ELEMENTS FROM ALL DIMENSIONS TO ALLOCATE MEMORY
			uint total = 1;
			for(i=0;i<numDims;i++){
				total *=total_elements_per_dimension[i];
			}
			//printf("total is %d\n", total);

			//printf("\n\n/**************** COMBINE CELLS FOR EACH DIMENSION ************************\n"); 
			uint **result = (uint**)malloc(sizeof(uint*)*total);  
			uint *temp = (uint*)malloc(sizeof(uint)*numDims); 
			uint *pos_dim = (uint*)malloc(sizeof(uint)*numDims); 

			result[0] = (uint*)malloc(sizeof(uint)*numDims); 
			for(i=0;i<numDims;i++){
				temp[i]= combine[i][0];
				pos_dim[i]= 0;
				//printf("%d\n",temp[i] );
			}
			
			uint d = numDims; 
			uint p =0;
			j = 0, i = 0;	

			for(j=0;j<total;j++){
				result[i] = (uint*)malloc(sizeof(uint)*numDims); 
			    for (p = 0; p < numDims; p++) {
			        result[i][p] = temp[p];
			    //  	printf(" %d ", result[i][p] );
			    }
			    // printf("\n\n");
		   	 	i++;
			    for (d = numDims - 1; d >= 0; d--) {
			    	pos_dim[d]++;
			        if (pos_dim[d] >= total_elements_per_dimension[d] ) {
			            if (d == 0) {
			                break;
			            } else {
			            	pos_dim[d]=0;
			                temp[d] = combine[d][0];
			            }
			        } else {
			        	temp[d] = combine[d][pos_dim[d]];
			            break;
			        }
			    }
			}
		
			free(total_elements_per_dimension);
			free(temp);
			free(pos_dim);
			free(combine);

			//SHOWS COMBINATORIAL
			/* for(i=0;i<total;i++){	
			 	for(j=0;j<numDims;j++){
			 		printf("%d ",result[i][j] );
			 	}
			 	printf("\n");	
			 } */
			 	
			begin = clock();												//init clock
			uint levels = lowest_level-1;
			uint **ancestor = (uint**)malloc(sizeof(uint*)*(levels+1));		//allocates memory for ancestor
		
			ancestor[levels] = (uint*)malloc(sizeof(uint)*(total*numDims)); //allocates memory for ancestor last level
	
			uint k = 0;
			for(i = 0; i<total; i++){		//assigns the result of combining cells to the lowest level
				for(j=0;j<numDims;j++){
					ancestor[levels][k] = result[i][j]; 	
				//	printf("%d ",ancestor[levels][k] );
					k++;
				}
			//	printf("\n");	
			}

			for(i=0;i<total;i++){
				free(result[i]);
			}
			free(result);

			//printf("\n\n/*********************** LOOK FOR THE ANCESTOR **************************/\n");
				
			dim = 0, j = 0;
		
			//looks for the parent of each level in each dimension
			while(levels>0){
				ancestor[levels-1]=(uint*)malloc(sizeof(uint)*(total*numDims));		//allocates memory for ancestor levels
				for(j=0; j<total*numDims; j++){
					ancestor[levels-1][j] = parent(array_dims[dim]->bitmap,ancestor[levels][j]);
				//	printf("%d ",ancestor[levels-1][j]  );
					dim++;
					if(dim == numDims) dim = 0;
				}
				levels--;	
			}

			levels = lowest_level-1;
		
			//CALCULATES THE NUMBER OF ONES BETWEEN THE PREVIOUS AND NEXT 0 TO IMPLEMENT THE FORMULA

			/*
				rank0(pos) = x
				select0(x)+1 = y
				Rx = pos-y 

				//Iterative
				(R1*D2)+R2 = X
				(X*D3)+R3  = X		where Dx is the number of ones in the next dimension
			*/

			uint **offset = (uint**)malloc(sizeof(uint*)*(levels+1));
			formel *formula = (formel*)malloc(sizeof(formel));
			formula->rank0s = (uint*)malloc(sizeof(uint)*(total*numDims));		//number of rank0 before the position in louds
			formula->select0s = (uint*)malloc(sizeof(uint)*(total*numDims));	//select0 to obtain the position in louds 
			formula->ones = (uint**)malloc(sizeof(uint*)*(levels+1));			//number of ones between zeros
			formula->result = (uint**)malloc(sizeof(uint*)*(levels+1)); 		//results to obtain the offsets
	
			uint rank0 = 0, pos_prev0 = 0, pos_next0 = 0, jj = 0;
			dim = 0;

			while(levels>=0){ 
				formula->result[levels]=(uint*)malloc(sizeof(uint)*total*numDims);
				formula->ones[levels]=(uint*)malloc(sizeof(uint)*total*numDims);
				for(jj=0;jj<total*numDims;jj++){	
					formula->rank0s[jj] = ancestor[levels][jj] - rank(array_dims[dim]->bitmap,ancestor[levels][jj]-1);
				//	printf("Rank 0's is %d\n",formula->rank0s[jj] );
					formula->select0s[jj] = select0(array_dims[dim]->bitmap,formula->rank0s[jj])+1;
				//	printf("Position is %d\n",formula->select0s[jj] );
					
					//CALCULATE THE NUMBER OF 1's BETWEEN PREVIOUS AND NEXT 0 TO USE THE FORMULA
					rank0 = ancestor[levels][jj] - rank(array_dims[dim]->bitmap,ancestor[levels][jj]-1);
				//	printf("rank0 is %d\n",rank0);
					pos_prev0 = select0(array_dims[dim]->bitmap,rank0);
				//	printf("pos_prev0 is %d\n",pos_prev0);
					pos_next0 = select0(array_dims[dim]->bitmap,rank0+1);
				//	printf("pos_next0 is %d\n",pos_next0);
					formula->ones[levels][jj] = (pos_next0 - pos_prev0)-1;
					formula->result[levels][jj] = ancestor[levels][jj] - formula->select0s[jj];
					dim++;
				//	printf("Result[%d][%d] is %d\n",levels,jj, formula->result[levels][jj] );
				//	printf("ones[%d][%d] is %d\n",levels,jj, formula->ones[levels][jj]);
					if(dim == numDims) dim = 0;
				}

				levels--;
				if(levels == -1) break;
			}

			levels = lowest_level-1;
			uint pos = 0, dims = 0;
			
			while(levels>=0){ 
				offset[levels]=(uint*)malloc(sizeof(uint)*(total));
				for(jj=0; jj<total*numDims; jj=jj+(numDims)){		//FOR 2 DIMENSIONS
					offset[levels][pos] = (formula->result[levels][jj]*formula->ones[levels][jj+1])+formula->result[levels][jj+1];
					// printf("Result[%d][%d] is %d\n",levels,jj, formula->result[levels][jj] );
					// printf("ones/desplaz[%d][%d] is %d\n",levels,jj+1, formula->ones[levels][jj+1]);
					// printf("OFFSET[%d][%d] is %d\n",levels,pos,offset[levels][pos] );
					// printf("\n");
					for(dims=2; dims<numDims; dims++){		//MORE THAN 2 DIMENSIONS
						// printf("Ancestor[%d][%d] %d\n",levels,jj, ancestor[levels][jj]);
						// printf("Result[%d][%d] is %d\n",levels,jj+1, formula->result[levels][jj+1] );
						// printf("ones/desplaz[%d][%d] is %d\n",levels,jj+1, formula->ones[levels][jj+1]);
						// printf("OFFSET[%d][%d] is %d\n",levels,jj,offset[levels][jj] );
						offset[levels][pos] = (offset[levels][pos]*formula->ones[levels][jj+dims])+formula->result[levels][jj+dims];
					}
				//	printf("pos %d\n",pos );
				//	printf("offset[%d][%d] %d \n",levels,pos,offset[levels][pos] );
					pos++;
				}

				pos = 0;
				levels--;
				if(levels == -1) break;
			}

			free(ancestor);

			// printf("\n\n/************** OBTAINS THE SON FROM EACH POSITION AT EACH LEVEL ***************/\n\n");

			uint total_offset = 0, rr2 = 0, son = 0;
			
			levels = 0; 
		
			while(rr2<total){
			
				levels = 0; 
				son=0;
				if(bitmap_dims->regenerate == 0){	//ARRAY_SUM
					if(lowest_level == 1){
							total_cube += accessFT(array_sumDac, offset[levels][rr2]+1);
						//	printf("DISTINCT LEVEL ARRAYSUM Query is %d\n",bitmap_dims->array_sum[offset[levels][rr2]]);
							break;
					}
					else{
						while(levels<lowest_level){   
							total_offset = offset[levels][rr2] + son + 1;
							son = select0(bitmap_dims->bitRank2,total_offset)+1;
							levels++;
							if(levels == lowest_level-1){
								total_offset = son + offset[levels][rr2];
								total_sum += accessFT(array_sumDac, total_offset+1);
							//	printf("Query is %d\n",bitmap_dims->array_sum[total_offset]); 
							}
						}
					}
				}
				else{		//ARRAY_SUM2

					total_offset = 0,  son = 0; //rr2 = 0,
					levels = 0;
					
					if(lowest_level == 1){
						if(isBitSet(bitmap_dims->bitRank1,offset[levels][rr2]) == 0) {//bitmap_dims->array_sum[
						//	printf("Query is %d\n", bitmap_dims->array_sum2[offset[levels][rr2]] );
						//	printf("1Sum -> %u \n", 0);
							break;
						}
						else{
							total_cube += accessFT(array_sumDac, offset[levels][rr2]+1);
						//	printf("2Sum -> %u \n", accessFT(array_sumDac, offset[levels][rr2]+1));
						//	printf("Query is %d\n", accessFT(array_sumDac, offset[levels][rr2]+1));
						} 
					}
					else{						
						while(levels<lowest_level){	
						
							if(levels == lowest_level-1){
								total_offset = son + offset[levels][rr2];
								total_sum += accessFT(array_sumDac, total_offset+1);
							//	printf("3Sum -> %u \n", accessFT(array_sumDac, total_offset+1));
								break;
							}

							total_offset = son + offset[levels][rr2];
							if(isBitSet(bitmap_dims->bitRank1,total_offset) ==0) {
							//	printf("4Sum -> %u \n", 0);
								break;
							}
							total_offset = rank(bitmap_dims->bitRank1,total_offset);
							son = select0(bitmap_dims->bitRank2,total_offset)+1;	
							levels++;
						}
					}
				}
				rr2++;
			} 

			for(i=0;i<levels;i++){
				free(offset[i]);
			}
			free(offset);

			free(sons);

		//	printf("Query is %d\n",total_sum );
			total_cube +=total_sum;
			
			end = clock();
			time_spent += (double)(end - begin) / (CLOCKS_PER_SEC);
		}	//end else different level queries
	}		//end while all queries

	end = clock();

	printf("time = %.7f, results = %lu, us/query = %.7f, us/res=%.7f\n",time_spent , total_cube, (time_spent*1000000) /numQueries, (time_spent*1000) / total_cube );
//	printf("time = %.7f, results = %lu, us/query = %.7f, us/res=%.7f\n",time_spent / 1000, total_cube, (time_spent*1000) /numQueries, (time_spent*1000) / total_cube );


	// printf("\n\n/*********** FREE MEMORY ******************/\n");
	
	//FREE DATAWAREHOUSE
	free_datawarehouse(datawarehouse, numQueries, numDims);

	//FREE BITMAP_DIMS
	if(bitmap_dims->regenerate == 0)
		free_bitmap_dims(bitmap_dims,0);
	else free_bitmap_dims(bitmap_dims,1);

	//FREE DIMENSIONS
	free_dims(array_dims,numDims);

	//FREE HASH TABLE
	freeHashTable();
	

	return EXIT_SUCCESS;
}

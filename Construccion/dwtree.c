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

#include "dwtree.h"


//calculates the time of construction
double get_process_time() {
    struct rusage usage;
    if( 0 == getrusage(RUSAGE_SELF, &usage) ) {
        return (double)(usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) +
               (double)(usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1.0e6;
    }
    return 0;
}


//cleans the character \n of each line
void cleanCharacter(char *string){		//clean the newline character of a string that is passed as parameter
	
	char *character;

	character = strchr (string, '\n');	//the character in the string is searched
	if (character)
		*character = '\0';				//it is replaced by end of string
}


//reads the vocabulary of each dimension to obtain the vocabulary and the dimension size
void read_vocabulary(char *filein, dim *dimension){

	char buffer[BUFFER_SIZE];	
	char *word;
	uint i = 0; 

	dimension->vocabulary = (char**)malloc(sizeof(char*)*MAX_N);	//allocates memory for the number of elements per dimension (vocabulary)

	FILE *input = fopen(filein, "r");		//open for reading

	if(input == NULL){
		perror(filein);
		exit(1);
	}

	while(fgets(buffer, BUFFER_SIZE, input)!= NULL){	//reads each line
		word = strtok(buffer, "\n");
		while(word!= NULL){
			dimension->vocabulary[i] = (char*)malloc(sizeof(char)*strlen(word));
			strcpy(dimension->vocabulary[i],word);
			word = strtok(NULL, "\n");
			i++;
		}
	}

	dimension->dimension_size = i;			//assigns the dimension size
	
	fclose(input);							//close file input
}


//reads the bitmap of each dimension to get the vars: bitmap, levels, total_elements_per_level and elements_per_level
void read_dimension(char *filein, dim *dimension){

	char buffer[BUFFER_SIZE];	
	char *word;
	char *token;
	uint i = 0, j = 0, k = 0, l = 0, m = 0; 
	uint cont = 0, count_ones = 0;

	uint *sourcebitmap = (uint*)malloc(((MAX_N+31)/32)*sizeof(uint));	//allocates memory for the bitmap

	FILE *input = fopen(filein, "r");		//open for reading

	if(input == NULL){
		perror(filein);
		exit(1);
	}
	
	dimension->levels = atoi(fgets(buffer,BUFFER_SIZE,input));		//gets the number of levels

	word = fgets(buffer,BUFFER_SIZE,input);							//reads bitmap
	cleanCharacter(word);
	token = strtok(word," ;,.-"); 									//there are five delimiters here 
	dimension->elements_per_level = (uint**)malloc(sizeof(uint*)*(dimension->levels+1));
	uint *total = (uint*)malloc(sizeof(uint)*dimension->levels);

    while (token != NULL){
    	dimension->elements_per_level[j] = (uint*)calloc(MAX_DIVS_LEVEL,sizeof(uint));	//allocates memory for elements_per_level
		for(k = 0; k<strlen(token); k++){				
			if(token[k] == '1'){						//if bitmap from file has 1 it puts a 1 in the sourcebitmap, else it puts a 0
				cont++;
				bitset(sourcebitmap,m);
				m++;
				count_ones++;
			}
			if (token[k]=='0'){
				dimension->elements_per_level[j][l] = cont;
				l++;
				cont = 0;
				bitclean(sourcebitmap,m);
				m++;
			}	
		}
		l=0;
		total[i] = count_ones;							//gets the number of 1's to assign them to total_elements_per_level
		count_ones = 0;
		j++;
		token = strtok(NULL," ;,.-");
		i++;
    }

    dimension->bitmap = createBitRankW32Int(sourcebitmap,m,1,20);	//creates the bitmap

	dimension->total_elements_per_level = (uint*)malloc(sizeof(uint)*(dimension->levels));	//allocates memory for total_elements_per_level

	for(j=0;j<=dimension->levels;j++){
		dimension->total_elements_per_level[j] = total[j];			//assigns total_elements_per_level in dimension
	//	printf("%d\n", dimension->total_elements_per_level[j] );
	}

	free(total);

	fclose(input);
}


//prints the vocabulary of one dimension
void print_vocabulary(dim *dimension){
	
	uint i = 0;
	for(i=0;i<dimension->dimension_size;i++){
		printf("%s \n",dimension->vocabulary[i]);
	}
}


//prints the bitmap of one dimension
void print_bitmap(bitRankW32Int *bitmap){
	
	uint k = 0;
	uint nBits;
	
	nBits = bitmap->n;							//get the number of bits
	printf("bitmap-> %d bits\n", nBits);
	for (k=0;k<nBits;k++) {						
		printf("%d",bitget(bitmap->data,k));	//prints the bitmap
	}
}


//prints the elements_per_level of one dimension
void print_elements_per_level(dim *dimension){

	uint i = 0, j = 0;

	for(i=1;i<dimension->levels+1;i++){
		for (j = 0; j<dimension->total_elements_per_level[i-1];j++){ //dimension->total_elements_per_level[dimension->levels-2]
			//if (dimension->elements_per_level[i][j]!=0)
				printf("%u ",dimension->elements_per_level[i][j]);
		}
		printf("\n");
	}
}


//prints the final dimension with all the sizes
void print_final_dimension(dim **array_dims, uint numDims){

/*	uint i, j = 0; 

	for(i=1;i<=array_dims[numDims]->levels;i++){ //realmente hay que bajar tantas como niveles tenga CAMBIAR
		for(j=0;j<array_dims[numDims]->total_elements_per_level[i-1];j++){
			//printf("elems per level %d \n", array_dims[numDims]->total_elements_per_level[i-1]);
			if((array_dims[numDims]->elements_per_level[i][j]!= 0))
				printf("%d ",array_dims[numDims]->elements_per_level[i][j]);
		}
	}
	printf("\n"); */

	uint y = 0, z = 0;
	printf("\n");
	for(y=0;y<array_dims[numDims]->levels+1;y++){
		for(z=0;z<array_dims[numDims]->total_elements_per_level[y];z++){
			if(array_dims[numDims]->elements_per_level[y][z]!=0)
			printf("%d ",array_dims[numDims]->elements_per_level[y][z]);
		}
		printf("\n");
	} 

}


//frees the allocated memory
void free_dims(dim **array_dims, uint numDims){
	
	uint x = 0, y = 0; 	

	//frees array_dim_result (total elements & elements per level)
	free(array_dims[numDims]->total_elements_per_level);

	for(y=0;y<array_dims[x]->levels+1;y++){
		free(array_dims[numDims]->elements_per_level[y]);		
	}
	free(array_dims[numDims]->elements_per_level);

	free(array_dims[numDims]);

 	//frees dimensions' vocabulary, bitmaps, total_elements and elements_per_level 
	for(x=0;x<numDims;x++){

		free(array_dims[x]->total_elements_per_level);

		for(y=0;y<array_dims[x]->dimension_size;y++){
			free(array_dims[x]->vocabulary[y]);		
		}
		free(array_dims[x]->vocabulary);

		for(y=0;y<array_dims[x]->levels+1;y++){
			free(array_dims[x]->elements_per_level[y]);		
		}
		free(array_dims[x]->elements_per_level);

		destroyBitRankW32Int(array_dims[x]->bitmap);

		free(array_dims[x]);
	}

	
	free(array_dims);
}


//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
uint calculate_level(dim *array_dim_result,uint level){

	uint total=0; uint i=0;

		for(i=1; i<=level; i++){
			total +=array_dim_result->total_elements_per_level[i];
		//	printf("total is %d\n",total );
		}
		return total;
}


//multiplies all the levels between two dimensions
uint *multiply_dims(dim *array_dims1,dim *array_dims2, uint levels,uint *total_elements_in_level, uint rest){

	uint pos_D1_level_before = 0;	//position in D1 previous level
	uint pos_D2_level_before = 0;	//position in D2 previous level

	uint acum_D1_level_before = array_dims1->elements_per_level[levels-1][pos_D1_level_before];	//accumulator D1 previous level
	uint acum_D2_level_before = array_dims2->elements_per_level[levels-1][pos_D2_level_before];	//accumulator D2 previous level

	//OBTENER EL NUMERO CONCRETO PARA UTILIZAR EL TOTAL MALLOC EN CADA UNO DE LOS NIVELES, EJ: 1 6 35
	*total_elements_in_level = array_dims1->total_elements_per_level[levels-1]*array_dims2->total_elements_per_level[levels-1];
	//printf("total malloc is %d\n",*total_elements_in_level );

	uint temp = 0;
	uint *result = (uint*)calloc(*total_elements_in_level,sizeof(uint));
	uint pos_D1 = 0, pos_D2 = 0;	//position in the last level of D1 and D2
	uint posicion = 0;

	while(pos_D1<array_dims1->total_elements_per_level[levels-1]){
		while(pos_D2<array_dims2->total_elements_per_level[levels-1]){
			
			temp = array_dims1->elements_per_level[levels][pos_D1]*array_dims2->elements_per_level[levels][pos_D2];
			result[posicion] = temp;
			//printf("%d ",result[posicion] );
			if(pos_D2 == acum_D2_level_before-1){
				pos_D2 = pos_D2 - array_dims2->elements_per_level[levels-1][pos_D2_level_before]+1;
				break;
			}
			else pos_D2++;

			posicion++;
		}

		if(pos_D1+1 == acum_D1_level_before){
			pos_D1 = pos_D1+1 - array_dims1->elements_per_level[levels-1][pos_D1_level_before];
			pos_D2_level_before++;
			pos_D2 = acum_D2_level_before;
			acum_D2_level_before += array_dims2->elements_per_level[levels-1][pos_D2_level_before];
		}
		else pos_D1++;

		posicion++;

		if( (pos_D2 == array_dims2->total_elements_per_level[levels-1]) ){  
			pos_D2 = 0;
			pos_D2_level_before = 0;
			pos_D1_level_before++;
			pos_D1 = acum_D1_level_before;
			acum_D1_level_before += array_dims1->elements_per_level[levels-1][pos_D1_level_before];
			acum_D2_level_before = array_dims2->elements_per_level[levels-1][pos_D2_level_before];
		}
	}

	if(rest == 1)
		free(array_dims1->elements_per_level[levels]);

	return result;
}


//gets the final dimension as a result of multiplying all existing dimensions
void all_dimensions(dim **array_dims, uint numDims){

	/********** Multiplication between two dimensions ***********/

	uint last_level=0;

	array_dims[numDims]->levels = array_dims[0]->levels;
	array_dims[numDims]->elements_per_level = (uint**)malloc(sizeof(uint*)*(array_dims[numDims]->levels+1));
	array_dims[numDims]->elements_per_level[0] = (uint*)malloc(sizeof(uint));
 	array_dims[numDims]->elements_per_level[0][0] = 1;

 	array_dims[numDims]->total_elements_per_level = (uint*)malloc(sizeof(uint)*(array_dims[numDims]->levels+1));  

    uint total_elements_in_level = 0;

	for(last_level= array_dims[numDims]->levels;last_level>0;last_level--){ 
		array_dims[numDims]->elements_per_level[last_level] = multiply_dims(array_dims[0],array_dims[1],last_level, &total_elements_in_level,0 );
		array_dims[numDims]->total_elements_per_level[last_level-1] = total_elements_in_level;
		printf("\n");
	}


	/********** Multiplication among three dimensions ***********/
	
	last_level = 0;

	uint actual_dimension = 2;

	while(actual_dimension<numDims){
		for(last_level=array_dims[numDims]->levels;last_level>0;last_level--){ 
			array_dims[numDims]->elements_per_level[last_level] = multiply_dims(array_dims[numDims],array_dims[actual_dimension],last_level,&total_elements_in_level,1);
			array_dims[numDims]->total_elements_per_level[last_level-1] = total_elements_in_level;
		//	printf("\n");
		}
		actual_dimension++;
	//	printf("\n");
	}

	//CALCULATE LAST LEVEL
	uint i = 0, j = 0, total = 1;
	uint levels = array_dims[numDims]->levels;
	
	for(i=levels-1; i<levels; i++){
		for(j=0;j<numDims;j++){
			total*= array_dims[j]->total_elements_per_level[levels];
		}
	}

	array_dims[numDims]->total_elements_per_level[levels] = total;

}


//performs exponentiation recursive
uint power(int x,int y){
 	return y?(y%2?x:1)*power(x*x,y>>1):1;
}


//allocates memory for total elements and elements per level (array_dim_result)
void alloc_array_dim_result(dim **array_dims, uint numDims){

	array_dims[numDims]->levels = array_dims[0]->levels;		//assigns levels to array_dim_result

	//allocates memory for total elements per level
	array_dims[numDims]->total_elements_per_level = (uint*)calloc(array_dims[numDims]->levels+1,sizeof(uint));  

	//allocates memory for elements_per_level
	array_dims[numDims]->elements_per_level = (uint**)malloc(sizeof(uint*)*(array_dims[numDims]->levels+1));

	array_dims[numDims]->elements_per_level[0] = (uint*)calloc(1,sizeof(uint));	
	array_dims[numDims]->elements_per_level[0][0] = 1;				//assigns one element to the level 0 (root)

	uint i = 0; // j = 1;
	uint exponent = power(2,array_dims[numDims]->levels);			//makes exponentiation for the last level, i.e. 2^5 = 32
	uint total_bytes = power(exponent,numDims);						//makes exponentiation for the number of dimensions, i.e. 32^5
	for(i = array_dims[numDims]->levels ; i>=1; i--){
		//printf("total_bytes es %u\n",total_bytes );
		//allocates memory for elements per level until level 1
		array_dims[numDims]->elements_per_level[i] = (uint*)calloc(total_bytes,sizeof(uint));
	//	total_bytes = (total_bytes/(power(2,array_dims[numDims]->levels-j)) )+ (total_bytes/(power(2,array_dims[numDims]->levels*2)) );   
	//	j++;
		total_bytes = (total_bytes/2) + (total_bytes/8);
	}

}


//reads the bitmaps of all dimensions and stores them into matrices to use them in multiDims
void readBitmaps(char **argv, int numDims,bitRankW32Int ***bitmapss){

	uint numLevels;
    uint d = 0;

    char buffer[BUFFER_SIZE];	
	char *word;
	char *token;
	uint k = 0, l = 0, posBitmap = 0; 
	
	for (d = 0; d < numDims; d++) {							//loop for all the number of dimensions
	
		FILE *input = fopen(argv[3+numDims+d], "r");		//open for reading

		if(input == NULL){
			perror(argv[3+numDims+d]);
			exit(1);
		}
		//else printf("%s\n",argv[3+numDims+d] );

		numLevels = atoi(fgets(buffer,BUFFER_SIZE,input));							//gets the number of levels from input
		bitmapss[d] = (bitRankW32Int**)malloc(sizeof(bitRankW32Int*)*numLevels);	//allocates memory for each bitmap

		word = fgets(buffer,BUFFER_SIZE,input);										//reads bitmap from input
		cleanCharacter(word);
		token = strtok(word," ;,.-"); 												//there are five delimiters here

		l=0;
		
        while(token!= NULL){
    		posBitmap = 0;
   			uint *sourcebitmap = (uint*)calloc(sizeof(uint),((MAX_N+31)/32));		//allocates memory for the bitmap
			for(k = 0; k<strlen(token); k++){	
				if(token[k] == '1'){				//if bitmap from file has 1 it puts a 1 in the sourcebitmap, else it puts a 0
					bitset(sourcebitmap,posBitmap);
					posBitmap++;
				}
				else{
					posBitmap++;
				}	
			}
			
			if(l!=0) {
        		bitmapss[d][l-1] = createBitRankW32Int(sourcebitmap,posBitmap,1,20);	//creates the bitrank for each dimension and level
    		}
    		l++;
    		
    		token = strtok(NULL, " ;,.-");
		}
        	
		fclose(input);
	}

}


//multiplies recursively all dimensions and levels to obtain the final structure
void multiDims(uint numDimensions, uint numLevels, uint *starPosition, uint *endPosition, uint level, bitRankW32Int ***bitmaps, dim *array_dim_result){
    
    uint value = 1;
    uint isFinished = 0;
    uint d = 0;

    uint *tmpPositions = (uint *)calloc(numDimensions,sizeof(uint));
    for (d = 0; d < numDimensions; d++){
        tmpPositions[d] = starPosition[d];
    }
  
    while (!isFinished) {
        if (level < numLevels -1){
            uint *nStartPosition = (uint *) calloc(numDimensions,sizeof(uint));
            uint *nEndPosition = (uint *) calloc(numDimensions,sizeof(uint));
            value = 1;
           	for (d = 0; d < numDimensions; d++) {
                nStartPosition[d] = tmpPositions[d] == 0 ? 0 : (select0(bitmaps[d][level],tmpPositions[d]) - (tmpPositions[d] - 1));
                nEndPosition[d] = select0(bitmaps[d][level],tmpPositions[d] + 1) - (tmpPositions[d] + 1);

            //    printf("Dim %u level %u (%u): from %u to %u get %u \n", d, level, tmpPositions[d],
             //       nStartPosition[d], nEndPosition[d], nEndPosition[d] - nStartPosition[d] + 1);
             //   printf("\n");
             //   printf("nEndPosition[d] %u\n",nEndPosition[d] * (nEndPosition[d] - nStartPosition[d] + 1) );
                value *= (nEndPosition[d] - nStartPosition[d] + 1);
            }
          	
            multiDims(numDimensions, numLevels, nStartPosition, nEndPosition, level + 1, bitmaps,array_dim_result);
            free(nStartPosition);
            free(nEndPosition);
        }
        else{
        	value = 1;
            for (d = 0; d < numDimensions; d++) {
                value *= select0(bitmaps[d][level],tmpPositions[d] + 1 )  -
                         (tmpPositions[d] == 0 ? 0 : select0(bitmaps[d][level],tmpPositions[d]) + 1 );
            }      
        //    printf("%u ", value);
         
        }
		array_dim_result->elements_per_level[level][array_dim_result->total_elements_per_level[level]++] = value;
        
        //Next value
        for (d = numDimensions - 1; d >= 0; d--){
            tmpPositions[d] += 1;
            if (tmpPositions[d] > endPosition[d]){
                if (d == 0){
                    isFinished = 1;
                    break;
                }
                tmpPositions[d] = starPosition[d];
            } 
            else{
                break;
            }
        }
    }
    free(tmpPositions);
}


//fills the last level of total_elements_per_level
void fill_last_level_array_dim_result(dim *array_dim_result){
	uint p = 0;
	for(p = 0; p<array_dim_result->total_elements_per_level[array_dim_result->levels-1]; p++){
		array_dim_result->total_elements_per_level[array_dim_result->levels]+=array_dim_result->elements_per_level[array_dim_result->levels-1][p];
	//	printf("%d ", array_dim_result->elements_per_level[array_dim_result->levels][p]);
		
	}
//	printf("%d ",array_dim_result->total_elements_per_level[array_dim_result->levels] );
}


/*------------------------------------------------------------------
 Initilizes the structures used.
 ---------------------------------------------------------------- */
void initialize2(unsigned long sizeVoc){
	
	_memMgr = createMemoryManager();	// Creates a Memory Manager
    initialize_hash(sizeVoc);			// Initilization of data structures used by the hashTable

    positionInTH = (unsigned long*) malloc (sizeVoc * sizeof(unsigned long));
	zeroNode = 0;
}


//loads the vocabulary from all dimensions and assigns the positions of the hash table corresponding to LOUDS
void hash_table(dim **array_dims, uint numDims){

	uint total_size = 0;
	uint i = 0;

	for(i=0; i<numDims; i++){
		total_size += array_dims[i]->dimension_size;	//total size of the dimensions' vocabulary
	}
	//printf("total_size %d\n",total_size );

	initialize2(total_size); 			//initializes compressor for up to N_value words    
	uint z = 0, pos = 0, size = 0, dim=0;

	for(dim=0;dim<numDims;dim++){	
		for(pos=0;pos<array_dims[dim]->dimension_size;pos++){
			size = strlen(array_dims[dim]->vocabulary[pos]);
			//looks for a word in the hash table and returns its position in the vocabulary
			z = search ((unsigned char *)array_dims[dim]->vocabulary[pos], size, &addrInTH );

		    if (z==zeroNode) {
		    	//inserts a new word in a given position of the hashTable (position previously computed)
				insertElement ((unsigned char *) array_dims[dim]->vocabulary[pos], size, &addrInTH);
				//printf("voc is %s\n",array_dims[dim]->vocabulary[pos] );
				//printf("size is %d\n",size );
			 	//printf("%lu\n",addrInTH );
				hash[addrInTH].size = 0;
				hash[addrInTH].len = size;
				hash[addrInTH].posInBitVector = select1(array_dims[dim]->bitmap,pos+2);  	//stores the word's position from LOUDS in hash table
				hash[addrInTH].dimension = dim+1;											//stores the dimension (future use)
				//printf("posInBitVector %lu\n",hash[addrInTH].posInBitVector );
				positionInTH[zeroNode] = addrInTH;
				//printf("%lu\n", positionInTH[zeroNode]);
				//printf("NumElem %lu\n",NumElem );
				//printf("zeronode %lu\n", zeroNode);
				zeroNode++;
			}	
		}	
	}
}


//reads the input file with tags and values
void read_input_file(char *argv[], data *datawarehouse, uint total, uint numDims){

	char buffer[BUFFER_SIZE];	
	char *word;
	uint i = 0, j = 0, k = 0; 

	datawarehouse->tags = (char**)malloc(sizeof(char*)*(total*(numDims)));  
	datawarehouse->value = (uint*)malloc(sizeof(uint*)*(total));
	
	FILE *input = fopen(argv[2], "r");					//open for reading  

	if(input == NULL){
		perror(argv[2]);
		exit(1);
	}

	while(fgets(buffer, BUFFER_SIZE, input)!= NULL){	//reads each line from input file
		word = strtok(buffer, ",");
		cleanCharacter(word);							//cleans the character \n of each line

		for(i=0;i<numDims;i++){
			datawarehouse->tags[j] = (char*)malloc(sizeof(char)*strlen(word));
			strcpy(datawarehouse->tags[j],word);								//assigns tags from file to datawarehouse->tags
		//	printf("%s\n",datawarehouse->tags[j] );
			word = strtok(NULL, ",");
			j++;
		}
		datawarehouse->value[k] = atoi(word);									//assigns value form file to datawarehouse->value
		k++;
	}

	// for(i=0; i<total; i++){
	// 	printf("%d ",datawarehouse->value[i] );			//prints the values
	// }

	// for(i=0; i<total*numDims; i++){
	// 	printf("%s\n",datawarehouse->tags[i] );			//prints the tags
	// }

	fclose(input);	
}


//creates the bitmaps1 and 2, the bitrank2 and the array_sum
void create_bitmaps(bitmap *bitmap_dims, dim *array_dim_result){

	uint i, j, k;
	uint inc_sourcebitmap = 0;

	//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
	uint total_sum = calculate_level(array_dim_result,array_dim_result->levels);

	//printf("total sum es %d\n",total_sum );
	bitmap_dims->sourcebitmap1 = (uint*)malloc(((total_sum+31)/32)*sizeof(uint) );	//allocates memory for the bitmap
	bitmap_dims->sourcebitmap2 = (uint*)malloc(((total_sum+31)/32)*sizeof(uint) );	//allocates memory for the bitmap2

	//generates bitmaps according to the structure
	for(i=0;i<array_dim_result->levels;i++){
		for(j=0;j<array_dim_result->total_elements_per_level[i];j++){
			//printf("%d \n",array_dim_result->elements_per_level[i][j] );
			for(k=0;k<array_dim_result->elements_per_level[i][j];k++){
				bitset(bitmap_dims->sourcebitmap1,inc_sourcebitmap);		//if segment has 5 positions it puts 5 ones to sourcebitmap1
				if(k == array_dim_result->elements_per_level[i][j]-1){		//if segment has 5 positions it puts 4 ones and 1 zero to sourcebitmap2
					bitclean(bitmap_dims->sourcebitmap2,inc_sourcebitmap);	
				}
				else{
					bitset(bitmap_dims->sourcebitmap2,inc_sourcebitmap);
				}		
				inc_sourcebitmap++;		
			}
		}
	}

	bitmap_dims->bitRank1 = createBitRankW32Int(bitmap_dims->sourcebitmap1,total_sum,1,20);	//creates the bitmap1
	bitmap_dims->bitRank2 = createBitRankW32Int(bitmap_dims->sourcebitmap2,total_sum,1,20);	//creates the bitmap2
	bitmap_dims->array_sum = (uint*)calloc(total_sum,sizeof(uint));							//creates array_sum

	bitmap_dims->remove_leaves = 0;															//initializes remove_leaves

	// uint z = 0;

	/*
	for(z=0; z<total_sum; z++){
	 	printf("%d ",bitmap_dims->array_sum[z] );				//prints array_sum
	}

	// printf("\n\n");

	for(z=0; z<total_sum; z++){
	 	printf("%d ",bitget(bitmap_dims->sourcebitmap1,z) );	//prints sourcebitmap1
	}

	printf("total is %d\n",z );

	for(z=0; z<total_sum; z++){
		printf("%d ",bitget(bitmap_dims->sourcebitmap2,z) );	//prints sourcebitmap2
	}

	printf("total is %d\n",z );
	*/
}


//looks for the ancestor at every level in every dimension
void look_for_ancestor(dim **array_dims, data *datawarehouse, uint total, uint numDims){

	uint levels = array_dims[numDims]->levels-1;

	datawarehouse->ancestor = (uint**)malloc(sizeof(uint*)*(levels+1));
	
	datawarehouse->ancestor[levels] = (uint*)malloc(sizeof(uint)*(total*(numDims)));		
	
	unsigned long lookfor = 0;
	uint size = 0;
	uint i=0;

	for(i=0;i<(total*(numDims));i++){ 
		//	printf("%s\n",datawarehouse->tags[i] );
			size = strlen(datawarehouse->tags[i]);						//gets the tag's size
		//	printf("size is %d\n",size );
		// 	printf("datawarehouse->tags[i] %s\n", datawarehouse->tags[i]);
			//looks for a word in the hash table and returns its position in the vocabulary
			lookfor = search ((unsigned char *)datawarehouse->tags[i], size, &addrInTH );
		//	printf("%lu\n",lookfor );
		//	printf("%s\n",hash[addrInTH].word );
			if(lookfor==0)		
				datawarehouse->ancestor[levels][i] = hash[addrInTH].posInBitVector;  //stores the position from LOUDS (last level)
				
		//	printf("posInBitVector %lu\n",hash[addrInTH].posInBitVector );
		//	printf("dimension %lu\n",hash[addrInTH].dimension );
	}

	//SHOWS THE POSITION IN BIT VECTOR
	//for(i=0;i<(total*(numDims));i++){
	//	printf("%d ",datawarehouse->ancestor[levels][i]);		//prints ancestor's last level 
	//}

	uint dim = 0;
	uint j = 0;

	//looks for the parent of each level in each dimension (from penultimate level)
	while(levels>0){
		datawarehouse->ancestor[levels-1]=(uint*)malloc(sizeof(uint)*(total*(numDims)));
		for(j=0;j<total*(numDims);j++){
			//looks for the parent in louds per dimension
			datawarehouse->ancestor[levels-1][j] = parent(array_dims[dim]->bitmap,datawarehouse->ancestor[levels][j]);	
			//printf("datawarehouse->ancestor[%d][%d] --> %d\n",levels-1,j, datawarehouse->ancestor[levels-1][j] );
			dim++;
			if(dim == numDims) dim = 0;
		}

		levels--;	
	}

}


//calculates the offset of each level using the formula
void calculate_offset(formel *formula, dim **array_dims, data *datawarehouse, uint total, uint numDims){

	//CALCULATE THE NUMBER OF ONES BETWEEN THE PREVIOUS AND NEXT 0 TO IMPLEMENT THE FORMULA

	/*
		rank0(pos) = x 		*1*
		select0(x)+1 = y   	*2*
		Rx = pos-y          *3*

		//Iterative
		(R1*D2)+R2 = X 		*4*
		(X*D3)+R3  = X		where Dx is the number of ones in the next dimension  *5*
	*/

	uint levels = array_dims[0]->levels-1;

	datawarehouse->offset = (uint**)malloc(sizeof(uint*)*(levels+1));	//stores the final offset per level
	formula->rank0s = (uint*)malloc(sizeof(uint)*(total*(numDims)));	//number of rank0 before the position in louds
	formula->select0s = (uint*)malloc(sizeof(uint)*(total*(numDims)));	//select0 to obtain the position in louds 
	formula->ones = (uint**)malloc(sizeof(uint*)*(levels+1));			//number of ones between zeros	
	formula->result = (uint**)malloc(sizeof(uint*)*(levels+1));			//results to obtain the offsets	
			
	uint rank0 = 0, pos_prev0 = 0, pos_next0 = 0;
	
	uint dim = 0;
	uint jj = 0;

	while(levels>=0){
		formula->result[levels]=(uint*)malloc(sizeof(uint)*(total*(numDims)));
		formula->ones[levels]=(uint*)malloc(sizeof(uint)*(total*(numDims)));
		for(jj=0;jj<(total*numDims);jj++){

			formula->rank0s[jj] = datawarehouse->ancestor[levels][jj] - rank(array_dims[dim]->bitmap,datawarehouse->ancestor[levels][jj]-1); // *1*
			//printf("Rank 0's is %d\n",formula->rank0s[jj] );
			formula->select0s[jj] = select0(array_dims[dim]->bitmap,formula->rank0s[jj])+1;		// *2*
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
			//printf("Result[%d][%d] es %d\n",levels,jj, formula->result[levels][jj] );
			//printf("ones[%d][%d] es %d\n",levels,jj, formula->ones[levels][jj]);
			if(dim == numDims) dim = 0;
		}

		levels--;	
		if(levels == -1) break;
	}
	

	// NOW IS THE MOMENT TO CALCULATE THE OFFSET IN THE BITMAP 2

	levels = array_dims[0]->levels-1;
		
	uint pos =0 ; uint dims = 0;
	
	while(levels>=0){

		datawarehouse->offset[levels]=(uint*)malloc(sizeof(uint)*(total));
		for(jj=0;jj<total*(numDims);jj=jj+(numDims)){
			//FOR 2 DIMENSIONS
			datawarehouse->offset[levels][pos] = (formula->result[levels][jj]*formula->ones[levels][jj+1])+formula->result[levels][jj+1]; // *4*
			// printf("\n\n");
			// printf("Result[%d][%d] is %d\n",levels,jj, formula->result[levels][jj] );
			// printf("ones/desplaz[%d][%d] is %d\n",levels,jj+1, formula->ones[levels][jj+1]);
			// printf("OFFSET[%d][%d] is %d\n",levels,jj,datawarehouse->offset[levels][jj] );
			// printf("\n");
			for(dims=2;dims<numDims;dims++){	//MORE THAN 2 DIMENSIONS
				// printf("ANCestor[%d][%d] %d\n",levels,jj, datawarehouse->ancestor[levels][jj]);
				// printf("Result[%d][%d] is %d\n",levels,jj+1, formula->result[levels][jj+1] );
				// printf("ones/desplaz[%d][%d] is %d\n",levels,jj+1, formula->ones[levels][jj+1]);
				// printf("OFFSET[%d][%d] is %d\n",levels,jj,datawarehouse->offset[levels][jj] );
				datawarehouse->offset[levels][pos] = (datawarehouse->offset[levels][pos]*formula->ones[levels][jj+dims])+formula->result[levels][jj+dims]; // *4* && *5*
			}
			
				// printf("offset[%d][%d] %d \n",levels,pos,datawarehouse->offset[levels][pos] );
			pos++;
		}

		pos = 0;
		levels--;
		if(levels == -1) break;
	}
	
}


//fills the leaves of the CMHD with the read input file
void fill_leaves(bitmap *bitmap_dims, data *datawarehouse, dim **array_dims, uint total, uint numDims){

	uint total_offset = 0, rr = 0, son = 0;
	uint levels = 0;  

	while(rr<total){
		levels = 0;
		son = 0;

		while(levels<array_dims[numDims]->levels-1){	
			total_offset = datawarehouse->offset[levels][rr] + son +1;				//offset in level
	//		printf("datawarehouse->offset[%d][%d] %d\n",levels,rr,datawarehouse->offset[levels][rr]  );
	//		printf("total_offset %d\n",total_offset );
			son = select0(bitmap_dims->bitRank2,total_offset)+1;					//selects the i-th zero
	
			levels++;
			if(levels == array_dims[numDims]->levels-1){
				total_offset =  son + datawarehouse->offset[levels][rr];			//offset in level
				bitmap_dims->array_sum[total_offset] = datawarehouse->value[rr]; 	//gets the value to array_sum
	//			printf("Value %d, Position %d\n\n", bitmap_dims->array_sum[total_offset], total_offset);
			}		
		}
		rr++;
	}

	//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
	uint total_sum = calculate_level(array_dims[numDims],array_dims[numDims]->levels);
	printf("total is %d\n",total_sum );
/*	uint z = 0;
	
	for(z = 0; z<total_sum; z++){
 		printf("%d ",bitmap_dims->array_sum[z] );		//prints the leaves from array_sum
	} 
	printf("total is %d\n",z ); 
*/
	
}


//sums cells to obtain the aggregated sum at upper levels
void sum_dw_tree(dim *array_dim_result,bitmap *bitmap_dims){

	uint partial = 0;				//flag to determine if a segment has regions with 0's
	uint dad = 0;					//ancestor at upper level
	uint rank0 = 0;
	bitmap_dims->regenerate = 0;	//flag to regenerate bitmaps and array_sum
	uint i = 0, j = 0, k = 0;
	
	uint levels = array_dim_result->levels;  
	
	while(levels>1){
		//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
		i = calculate_level(array_dim_result,levels-1);	
	//	  printf("I is %d\n",i );
	//	  printf("levels is %d\n",levels );

		for(j=0;j<array_dim_result->total_elements_per_level[levels];j++){
			partial = 0;
			// printf("%d\n", array_dim_result->total_elements_per_level[levels]);
			// printf("%d\n", array_dim_result->elements_per_level[levels][j]);
			for(k=0;k<array_dim_result->elements_per_level[levels-1][j];k++){  	//looking for a segment with all zeros
				// printf("j is %d\n",j );
				// printf("i is: %d\n",i );
				// printf("k is %d\n",k );
				// printf("array_sum: %d\n",bitmap_dims->array_sum[i+k] );
				// printf("i + k are %d\n",i+k );

				partial += bitmap_dims->array_sum[i+k];
			//	if(array_dims[numDims]->elements_per_level[levels][j]== (k-1)) i++;
			//	printf("partial is %d \n",partial );
				if(partial == 0 ) bitmap_dims->regenerate = 1;	//if a segment has 0's, the flag turns into 1
			}
	//		printf("i is: %d\n",i);
			rank0 = i - rank(bitmap_dims->bitRank2,i-1);		//counts the number of zeros until the position
	//		printf("rank0 is %d\n",rank0 );
			dad = select1(bitmap_dims->bitRank1,rank0);			//selects the i-th one (I think that here I can use parent from LOUDS)
	//		printf("dad at position %d\n",dad );
			
	//		printf("partial is out %d\n",partial );
	//		printf("array_sum is before %d\n",bitmap_dims->array_sum[dad] );
			bitmap_dims->array_sum[dad] = partial;
	//		printf("array_sum is after %d\n",bitmap_dims->array_sum[dad] );
			
			i=i+array_dim_result->elements_per_level[levels-1][j];
	//		printf("i is %d\n",i );
			if(i == calculate_level(array_dim_result,levels)) break;
		}
	//	printf("\n\n");
		levels--;	
		if(levels==1) break;
	}

	//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
	int total_sum = calculate_level(array_dim_result,array_dim_result->levels);

	
	uint z = 0;
/*	for(z=0; z<calculate_level(array_dim_result,array_dim_result->levels); z++){
		printf("%d ",bitmap_dims->array_sum[z] );										//prints the array_sum filled
	}
	printf("\n\n");
	printf("%d\n",z );
*/	
	// printf("\n\n");
	
	// printf("total sum %d\n",total_sum );
	for(z=0; z<total_sum; z++){
		if(bitmap_dims->array_sum[z]==0)
			bitclean(bitmap_dims->sourcebitmap1,z);			//modifies sourcebitmap1 according to array_sum
	//	printf("%d ",bitget(bitmap_dims->sourcebitmap1,z) );
	}
	//printf("\n\n");
	//printf("%d\n",z );

}


//removes the zeros of CMHD and generates new bitmaps
void remove_zeros(dim *array_dim_result, bitmap *bitmap_dims,char *argv){

	uint levels = 0;
	uint partial = 0;
	uint partial_sum = 0;
	uint p = 0, i = 0, j = 0, k = 0;
	
	uint *new_total_elements_per_level = (uint*)malloc(sizeof(uint)*array_dim_result->levels+1);	//aux array to store total elements per level
	
	for(p=0;p<array_dim_result->levels+1;p++){
		new_total_elements_per_level[p] = array_dim_result->total_elements_per_level[p];			//fill new total elements per level
		// printf("%d\n",array_dim_result->total_elements_per_level[p] );
	}

	uint **new_elements_per_level = (uint**)malloc(sizeof(uint*)*array_dim_result->levels+1);		//aux matrix to store elements per level
	uint q=0,r=0;

	for(q=0;q<=array_dim_result->levels-1;q++){														//fill new elements per level
		// printf("array_dim_result->levels %d\n",array_dim_result->levels );
		new_elements_per_level[q] = (uint*)calloc(new_total_elements_per_level[q],sizeof(uint));
		
		for(r=0;r<array_dim_result->total_elements_per_level[q];r++){	
			new_elements_per_level[q][r] = array_dim_result->elements_per_level[q][r];
			// printf("new_elements_per_level[q]%d[r]%d  %d\n",q,r,new_elements_per_level[q][r] );
		}
	}

	
	while(levels<=array_dim_result->levels-1){								//cuts the regions with zeros and moves the sums & leaves x positions
		for(i=0;i<array_dim_result->total_elements_per_level[levels];i++){
			//	printf("array_dims[numDims]->total_elements_per_level[%d] %d\n",levels, array_dim_result->total_elements_per_level[levels]);
			//	printf("array_dims[numDims]->elements_per_level[%d][%d] %d\n",levels,i, array_dim_result->elements_per_level[levels][i]);
			partial = 0;
			
			for(j=0;j<array_dim_result->elements_per_level[levels][i];j++){			//looking for a segment with all zeros
				partial += bitmap_dims->array_sum[k+j];
				//	printf("partial is %d\n",partial );
				//	printf("bitmap_dims->array_sum[k+j] is %d\n",bitmap_dims->array_sum[k+j] );
			}

			// printf("\n\n");
			// printf("partial before %d\n",partial );
			if(partial==0){															//segment with zeros
				partial_sum += array_dim_result->elements_per_level[levels][i];
				//	printf("array_dim_result->elements_per_level[levels][i] %d\n",array_dim_result->elements_per_level[levels][i] );
				new_total_elements_per_level[levels-1]= new_total_elements_per_level[levels-1]-1;
				//	printf("new_total_elements_per_level[levels-1] %d\n",new_total_elements_per_level[levels-1] );
				new_elements_per_level[levels][i] = 0;
				if(array_dim_result->levels>=levels){								//calculates the leaves that have been removed
					bitmap_dims->remove_leaves+=partial_sum;
				//	printf("bitmap_dims->remove %d\n",bitmap_dims->remove_leaves );
					partial_sum = 0;
				}
			}
			k +=j;						//increments position
		}
		// printf("\n\n");
		levels++;						//next level
	}

	//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
	bitmap_dims->total_sum = calculate_level(array_dim_result,array_dim_result->levels);
//	printf("total_sum is %d\n",bitmap_dims->total_sum );
	
	bitmap_dims->total_sum -= bitmap_dims->remove_leaves ;		//now the total sum turns into total sum less the leaves removed

	bitmap_dims->array_sum2 = (uint*)malloc(sizeof(uint)*bitmap_dims->total_sum);	//allocates memory for array_sum2
	levels = 0;
	k = 0;
	uint l = 0;
	uint m = 0;

	while(levels<=array_dim_result->levels-1){									//copies the new array_sum to array_sum2 without zero regions
		for(i=0;i<array_dim_result->total_elements_per_level[levels];i++){
			partial = 0;
			for(j=0;j<array_dim_result->elements_per_level[levels][i];j++){
				partial += bitmap_dims->array_sum[k+j];
			}
			if(partial!=0){
				for(l=0;l<array_dim_result->elements_per_level[levels][i];l++){
					bitmap_dims->array_sum2[m]= bitmap_dims->array_sum[k+l];
					// printf("bitmap_dims->array_sum2[m] %d\n",bitmap_dims->array_sum2[m] );
					// printf("bitmap_dims->array_sum[k+l] %d\n",bitmap_dims->array_sum[k+l] );
					m++;
				}
			}
			k +=j;			//increments position
		}
		levels++;			//next level
	}

	//EXPORT ARRAY_SUM2
	
	char dest[500];																	//file path size

	FTRep *arrayDac = createFT(bitmap_dims->array_sum2,bitmap_dims->total_sum);		//creates DAC
	strcpy(dest, argv);																																																																																																																																																												
	strcat(dest, "/array_sum2.dat");
	saveFT(arrayDac, dest);															//saves DAC

/*	q=0,r=0;

	for(q=0;q<=array_dim_result->levels-1;q++){										//prints new elements per level
		for(r=0;r<array_dim_result->total_elements_per_level[q];r++){
			if(new_elements_per_level[q][r]!=0)
			printf("%d ",new_elements_per_level[q][r] );
		}
		printf("\n\n");
	}
*/
	uint z = 0;
/*	for(z=0;z<bitmap_dims->total_sum;z++){											//prints array_sum2
	 	printf("%d ",bitmap_dims->array_sum2[z] );
	}

	printf("z is %d\n", z);
*/
	
	destroyBitRankW32Int(bitmap_dims->bitRank1);									//destroys & frees bitrank1
	destroyBitRankW32Int(bitmap_dims->bitRank2);									//destroys & frees bitrank2


	//printf("\n\n/********************** CREATE NEW SOURCEBITMAPS **********************/\n\n");

	uint inc_sourcebitmapnew2 = 0;

	//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
	//in this particular case, we don't need the leaves in the bitmap1
	uint nodes_remaining = calculate_level(array_dim_result,array_dim_result->levels-1);
	//printf("nodes_remaining %d\n", nodes_remaining);

	bitmap_dims->sourcebitmap1 = (uint*)calloc((((nodes_remaining)+31)/32),sizeof(uint) );		//allocates memory for the bitmap1
	bitmap_dims->sourcebitmap2 = (uint*)malloc(((bitmap_dims->total_sum+31)/32)*sizeof(uint) );	//allocates memory for the bitmap2

	i = 0, j = 0, k =0 , l = 0;

	for(i=0;i<=array_dim_result->levels-1;i++){
		for(j=0;j<array_dim_result->total_elements_per_level[i];j++){
		//	printf("new_elements_per_level[i]%d[j]%d is %d \n",i,j,new_elements_per_level[i][j] );
			if(new_elements_per_level[i][j]!=0){
				for(k=0;k<new_elements_per_level[i][j];k++){
					if(k == new_elements_per_level[i][j]-1){					//if segment has 5 positions it puts 4 ones and 1 zero to sourcebitmap2
						bitclean(bitmap_dims->sourcebitmap2,inc_sourcebitmapnew2);	//puts zero
						inc_sourcebitmapnew2++;
					}
					else{
						bitset(bitmap_dims->sourcebitmap2,inc_sourcebitmapnew2);	//puts one
						inc_sourcebitmapnew2++;
					}		
				}
			}
		}
	}
	
	z=0;

	for(z=0;z<nodes_remaining;z++){
		if(bitmap_dims->array_sum2[z]!=0){
			bitset(bitmap_dims->sourcebitmap1,z);		//puts 1 in sourcebitmap according to array_sum2
		}
	}

/*
	printf("\n\n");
	uint ww = 0, gg = 0;
	printf("BITMAP 2:\n");
	for(ww=0;ww<bitmap_dims->total_sum;ww++){
		gg= bitget(bitmap_dims->sourcebitmap2,ww);		//prints sourcebitmap2
		printf("%d ",gg );
	}
	printf("\n ww is %d\n", ww);

	printf("BITMAP 1:\n");
	printf("total sum is %d\n",bitmap_dims->total_sum );
	for(ww = 0; ww<nodes_remaining; ww++){		
		gg = bitget(bitmap_dims->sourcebitmap1,ww);		//prints sourcebitmap1
		printf("%d ",gg );
	}
	printf("\n nodes_remaining %d\n", ww);
*/

//	printf("/************************ CREATE NEW BITRANKS 1 & 2 ********************************/\n\n");
	
	bitmap_dims->bitRank1 = createBitRankW32Int(bitmap_dims->sourcebitmap1,nodes_remaining,1,20);			//creates the bitmap1 without leaves
	bitmap_dims->bitRank2 = createBitRankW32Int(bitmap_dims->sourcebitmap2,bitmap_dims->total_sum,1,20);	//creates the bitmap2 to navigate

//	printf("bitmap1 has %d\n",bitmap_dims->bitRank1->n );	//shows the number of bits of bitRank1
//	printf("bitmap2 has %d\n",bitmap_dims->bitRank2->n );	//shows the number of bits of bitRank2

	free(new_total_elements_per_level);				//frees new_total_elements_per_level

//	uint z=0;
	for(z=1;z<array_dim_result->levels+1;z++){		//frees new_elements_per_level
		free(new_elements_per_level[z]);
	}
	free(new_elements_per_level);
}



//initializes dimensions (reads vocabulary & dimension), moreover creates the array_dim_result (array_dims[numDims])
void initialize(dim **array_dims,char *argv[],uint numDims){

	uint x = 0;

	for(x=0;x<numDims;x++){
		array_dims[x] = (dim*)malloc(sizeof(dim));
		//reads the vocabulary of each dimension to obtain the vocabulary and the dimension size
		read_vocabulary(argv[x+3],array_dims[x]);			
		//reads the bitmap of each dimension to get the vars: bitmap, levels, total_elements_per_level and elements_per_level
		read_dimension(argv[x+numDims+3],array_dims[x]); 	
	}

	array_dims[numDims] = (dim*)malloc(sizeof(dim));	//array_dim_result
}


//frees the allocated memory for the dw
void free_datawarehouse(data *datawarehouse, uint total, uint numDims){

	free(datawarehouse->value);

	uint i = 0;
	for(i=0;i<total*(numDims);i++){
		free(datawarehouse->tags[i]);		//frees the tags
	}

	free(datawarehouse->tags);

	for(i=0;i<numDims;i++){
		free(datawarehouse->ancestor[i]);	//frees ancestor
	}

	free(datawarehouse->ancestor);

	for(i=0;i<numDims;i++){
		free(datawarehouse->offset[i]);		//frees offset
	}

	free(datawarehouse->offset);
	free(datawarehouse);
}


//frees the allocated memory for the formula
void free_formula(formel *formula, uint numDims){

	uint i = 0;
	for(i=0;i<numDims;i++){
		free(formula->ones[i]);				//frees ones
	}

	free(formula->ones);

	for(i=0;i<numDims;i++){
		free(formula->result[i]);			//frees result
	}

	free(formula->result);					
	free(formula->rank0s);
	free(formula->select0s);
	free(formula);

}


//frees the allocated memory for dimensions
void free_bitmap_dims(bitmap *bitmap_dims){

	destroyBitRankW32Int(bitmap_dims->bitRank1);	//destroys & frees bitRank1
	destroyBitRankW32Int(bitmap_dims->bitRank2);	//destroys & frees bitRank2
	free(bitmap_dims->array_sum);					//frees array_sum
	free(bitmap_dims);
}


//saves array_sum or array_sum2 and bitrank1 & 2
void save_bitmaps(bitmap *bitmap_dims, FILE *f, dim **array_dims,uint numDims){
	uint numOfLeaves = array_dims[numDims]->total_elements_per_level[array_dims[numDims]->levels];
	numOfLeaves -= bitmap_dims->remove_leaves;
//	printf("numOfLeaves %d\n",numOfLeaves );
	uint numOfNodes = bitmap_dims->total_sum - numOfLeaves;
//	printf("numOfNodes %d\n",numOfNodes );

	uint *newBitmap = (uint*)malloc(((numOfNodes+31)/32)*sizeof(uint) );			//allocates memory for the bitmap

	if(bitmap_dims->regenerate == 0){					//ARRAY_SUM
		uint i;
		for (i = 0; i < numOfNodes; i++) {
			if (isBitSet(bitmap_dims->bitRank1, i)) {
				bitset(newBitmap,i);					//puts a 1
			}
			else bitclean(newBitmap,i);					//puts a 0
		}
//		printf("numOfNodes %u of %u (%u) %u\n", numOfNodes, calculate_level(array_dims[numDims],array_dims[numDims]->levels), bitmap_dims->total_sum,bitmap_dims->remove_leaves);

		destroyBitRankW32Int(bitmap_dims->bitRank1);
		bitmap_dims->bitRank1 = createBitRankW32Int(newBitmap,numOfNodes,1,20);		//creates the bitmap1
		free(newBitmap);

//		printf("number of bits %u\n",bitmap_dims->bitRank1->n );					//prints the number of bits from bitrank1
//		printf("number of bits %u\n",bitmap_dims->bitRank2->n );					//prints the number of bits from bitrank2
		save(bitmap_dims->bitRank1, f);
		save(bitmap_dims->bitRank2, f);
	}
	else{												//ARRAY_SUM2
//		printf("number of bits %u\n",bitmap_dims->bitRank1->n );					//prints the number of bits from bitrank1
//		printf("number of bits %u\n",bitmap_dims->bitRank2->n );					//prints the number of bits from bitrank2
		save(bitmap_dims->bitRank1, f);
		save(bitmap_dims->bitRank2, f);
	}

	fwrite(&bitmap_dims->regenerate,sizeof(uint),1,f );	//stores the flag to indicate if you have to use array_sum or array_sum2
}


/*saves vocabulary, louds & dimension_size from all dimensions, 
moreover saves the number of levels, total elements per level and elements per level from array_dim_result (all merged dimensions) */
void save_dims(dim **array_dims, FILE *f, uint numDims, uint total){
	uint i = 0, j = 0;

	//ARRAY_DIMS
	for(i=0;i<numDims;i++){
		save(array_dims[i]->bitmap,f);										//saves bitmap to file
	}

	for(i=0;i<numDims;i++){
		fwrite(&array_dims[i]->dimension_size,sizeof(uint),1,f );			//saves dimension_size to file
	}

	uint lenWord;
	while(j<numDims){
		for(i=0;i<array_dims[j]->dimension_size;i++){
			lenWord = strlen(array_dims[j]->vocabulary[i]);					//gets the word's lenght 
			fwrite(&lenWord,sizeof(uint),1,f);								//saves lenght to file
			fwrite(array_dims[j]->vocabulary[i],sizeof(char),lenWord,f);	//saves vocabulary to file
		}
		j++;
	}

	//ARRAY_DIM_RESULT
	fwrite(&array_dims[numDims]->levels,sizeof(uint),1,f );					//saves levels to file (array_dim_result)
	array_dims[numDims]->total_elements_per_level[array_dims[numDims]->levels] = total;
	fwrite(array_dims[numDims]->total_elements_per_level, sizeof(uint), array_dims[numDims]->levels+1, f);	//saves total_elements_per_level

	FTRep *arrayDac = createFT(array_dims[numDims]->elements_per_level[0],1);	//creates DAC for level 0 //esto sobra
	saveFTFlex(arrayDac,f);
	free(arrayDac);
	
	for(i=0;i<array_dims[numDims]->levels;i++){  
		arrayDac = createFT(array_dims[numDims]->elements_per_level[i],array_dims[numDims]->total_elements_per_level[i]);	//creates DAC for each level
		saveFTFlex(arrayDac,f);
		free(arrayDac);
	}
	
}


//saves ancestor and offset from DW
void save_warehouse(data *datawarehouse, FILE *f, uint total, uint numDims){
	uint k = 0;

	for(k=0;k<numDims;k++){
		fwrite(datawarehouse->ancestor[k],sizeof(uint),total*numDims,f);
	}

	for(k=0;k<numDims;k++){			
		fwrite(datawarehouse->offset[k],sizeof(uint),total*numDims,f);
	}
}

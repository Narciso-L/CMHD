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





//MAIN BODY
int main(int argc, char *argv[]){

	uint numDims = atoi(argv[1]);								//gets the number of dimensions from command line

	if(argc!= (numDims*2)+4){
		printf("Input format: %s (filename) %d (N numDims) tags dim1(vocabulary) ... dimN(vocabulary) dim1bitmap ... dimNbitmap \n", argv[0],numDims );
		return -1;
	}

	if(numDims<2){
		printf("The number of dimensions has to be greater than 1\n");
		return -1;
	}

	
	dim **array_dims = (dim**)malloc(sizeof(dim*)*numDims+1);	//allocates memory for all dimensions +1 (dimension result)

	initialize(array_dims,argv,numDims);						//initializes dimensions (reads vocabulary & dimension), 
																//moreover creates the array_dim_result (array_dims[numDims])

	clock_t begin, end;
	double time_spent = 0.0;

	begin = clock();											//init clock


	printf("\n\n/********************** FEATURES **********************/\n\n");

//	uint dim = 0;
/*
	for(dim=0; dim<numDims; dim++){
		printf("Dimension %d vocabulary size: %d\n",dim+1,array_dims[dim]->dimension_size );	//prints vocabulary size
		printf("Number of levels: %d\n",  array_dims[dim]->levels);								//prints number of levels in the dimension
		printf("\n");
	}
*/

	printf("\n\n/********************** PRINT DIMENSIONS **********************/\n\n");
/*
	for(dim=0; dim<numDims; dim++){
		printf("Dimension %d\n",dim+1 );
		print_elements_per_level(array_dims[dim]);	//prints the elements_per_level of one dimension
		printf("\n");
	}
*/

	printf("\n\n/********************** PRINT BITMAPS **********************/\n\n");
/*
	for(dim=0; dim<numDims; dim++){
		printf("Dimension %d\n",dim+1 );
		print_bitmap(array_dims[dim]->bitmap);	//prints the bitmap 
		printf("\n\n");
	}
*/

	printf("\n\n/********************** ALL BITMAPS **********************/\n\n");

 	alloc_array_dim_result(array_dims,numDims);			//allocates memory for total elements and elements per level (array_dim_result)

	bitRankW32Int ***bitmapss = (bitRankW32Int***)malloc(sizeof(bitRankW32Int**)*numDims);	//allocates memory for dimensions bitmaps

	readBitmaps(argv,numDims,bitmapss);		//reads the bitmaps of all dimensions and stores them into matrices to use them in multiDims

	
    printf("\n\n/********************** FINAL DIMENSION **********************/\n");

    uint d = 0;
	uint *startPositions = (uint*)malloc(sizeof(uint)*numDims);
    uint *endPositions = (uint*)malloc(sizeof(uint)*numDims);
    for (d= 0; d < numDims; d++) {
        startPositions[d] = 0;
        endPositions[d] = 0;
    }

    //multiplies recursively all dimensions and levels to obtain the final structure
    multiDims(numDims, array_dims[0]->levels, startPositions, endPositions, 0, bitmapss, array_dims[numDims]);

    fill_last_level_array_dim_result(array_dims[numDims]);	//fills the last level of total_elements_per_level

	//all_dimensions(array_dims, numDims);
	//print_final_dimension(array_dims,numDims);


	printf("\n\n/********************** HASH TABLE **********************/\n");

	hash_table(array_dims,numDims);		//loads the vocabulary from all dimensions and assigns 
										//the positions of the hash table corresponding to LOUDS


	printf("\n\n/*********** CREATE SOURCEBITMAP 1 & 2 , BITRANK2 and ARRAY_SUM ***********/\n");

	bitmap *bitmap_dims = (bitmap*)malloc(sizeof(bitmap));	//allocates memory for the struct to store bitmaps and arrays of sums
//	bitmap_dims->remove_leaves = 0;
    create_bitmaps(bitmap_dims, array_dims[numDims]);		//creates the bitmaps1 and 2, the bitrank2 and the array_sum


	printf("\n\n/********************** READ INPUT FILE WITH TAGS AND VALUES **********************/\n");

	uint total = array_dims[numDims]->total_elements_per_level[array_dims[numDims]->levels];	//obtains the whole file elements
//	printf("%d\n",total );
	data *datawarehouse = (data*)malloc(sizeof(data));		//allocates memory for the struct to store the tags, parents and offsets per level

	read_input_file(argv, datawarehouse, total, numDims);	//reads the input file with tags and values


	printf("\n\n/********************** LOOK FOR THE ANCESTOR **********************/\n");

	look_for_ancestor(array_dims, datawarehouse, total, numDims);	//looks for the ancestor at every level in every dimension
	
	//FREE HASH TABLE
	//freeHashTable();


	printf("\n\n/********************** CALCULATE OFFSET AT EACH LEVEL **********************/\n");

	formel *formula = (formel*)malloc(sizeof(formel));		//allocates memory for the struct for implementing the formula and obtain offsets

	calculate_offset(formula, array_dims, datawarehouse, total, numDims);	//calculates the offset of each level using the formula

	
	printf("\n\n/********************** FILL THE LEAVES OF THE DW **********************/\n");

	fill_leaves(bitmap_dims, datawarehouse, array_dims, total, numDims);	//fills the leaves of the CMHD with the read input file

	//FREE FORMULA
	//free_formula(formula, numDims);
	//FREE DATAWAREHOUSE
	//free_datawarehouse(datawarehouse, total, numDims);


	printf("\n\n/********************** SUM **********************/\n");

	sum_dw_tree(array_dims[numDims],bitmap_dims);	//sums cells to obtain the aggregated sum at upper levels


	printf("\n\n/********************** REMOVE REGIONS OF ZEROS **********************\n");

	//printf("bitmap_dims->regenerate %d\n",bitmap_dims->regenerate );
	if(bitmap_dims->regenerate == 1){
		remove_zeros(array_dims[numDims],bitmap_dims,argv[argc - 1]);	//removes the zeros of CMHD and generates new bitmaps
		free(bitmap_dims->array_sum2);

		end = clock();
		printf("Elapsed time regenerating: %.7f seconds \n",time_spent = (double)(end - begin) / CLOCKS_PER_SEC);
	}

	end = clock();
	printf("Elapsed time: %.7f seconds \n",time_spent = (double)(end - begin) / CLOCKS_PER_SEC);


	printf("\n\n/********************** EXPORT DATA **********************/\n");

	char dest[500];				//path's file 

	//SAVE & EXPORT ARRAY_SUM

	//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
	bitmap_dims->total_sum = calculate_level(array_dims[numDims],array_dims[numDims]->levels);
	//printf("total sum es %d\n",bitmap_dims->total_sum );

	if(bitmap_dims->regenerate == 0){
		FTRep *arrayDac = createFT(bitmap_dims->array_sum,bitmap_dims->total_sum);	//creates DAC (array_sum)
		strcpy(dest, argv[argc - 1]);																																																																																																																																																												
		strcat(dest, "/array_sum.dat");	
		saveFT(arrayDac, dest);									//saves DAC (array_sum)
	}

	
	//SAVE & EXPORT BITMAPS 1 & 2
	strcpy(dest, argv[argc - 1]);
	strcat(dest, "/bitmaps.dat");
	FILE *bitmaps = fopen(dest, "w+"); 							//opens file for writing
	save_bitmaps(bitmap_dims, bitmaps,array_dims,numDims);		//saves array_sum or array_sum2 and bitrank1 & 2
	fclose(bitmaps);

	
	//SAVE & EXPORT LOUDS DIMS & VOCABULARY & DIMENSION SIZE
	strcpy(dest, argv[argc - 1]);
	strcat(dest, "/dims.dat");
	FILE *dims = fopen(dest,"w+"); 								//opens file for writing
	/*saves vocabulary, louds & dimension_size from all dimensions, moreover saves the number of levels, 
	total elements per level and elements per level from array_dim_result (all merged dimensions) */
	save_dims(array_dims, dims, numDims, total);				
	fclose(dims);


	//EXPORT ANCESTOR & OFFSET
	// FILE *warehouse = fopen("../Busqueda/warehouse.dat","w+"); 
	// save_warehouse(datawarehouse, warehouse, total, numDims);
	// fclose(warehouse);

	
	printf("\n\n/********************** FREE MEMORY **********************/\n");
	
	//FREE DATAWAREHOUSE
	//free_datawarehouse(datawarehouse, total, numDims);

	//FREE FORMULA
	//free_formula(formula, numDims);

	//FREE BITMAP_DIMS
//	free_bitmap_dims(bitmap_dims);

	//FREE DIMENSIONS
//	free_dims(array_dims,numDims);

	//FREE HASH TABLE
//	 freeHashTable();
	

	return EXIT_SUCCESS;

}
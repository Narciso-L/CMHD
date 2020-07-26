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

#include "louds.h"
#include "hash.h"
#include <ctype.h>
#include "dacs.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>


#define BUFFER_SIZE 2048			//line size
#define MAX_DIVS_LEVEL 128			//for elements_per_level per dimension
#define MAX_N 128					//number of elements per dimension (vocabulary, level...)


//struct to store each one of the levels and the size of each matrix
typedef struct dimension {
	char **vocabulary;				//stores the vocabulary of the dimension
	uint dimension_size;			//stores the number of words in the vocabulary
	bitRankW32Int *bitmap;			//stores the bitmap
	uint levels;					//stores the number of levels in the dimension
	uint *total_elements_per_level; //stores an array with all the elements per level, i.e. for dim1 -> 1 2 5		root is not counted
	uint **elements_per_level;		/*stores a matrix that counts the number of ones per level, i.e for dim1 -> 1 	this is root
																												2             L1
																												2 3           L2
																												4 2 2 5 2     L3  */
}dim;


//struct to store bitmaps and arrays of sums
typedef struct bm{
	uint *sourcebitmap1;
	uint *sourcebitmap2;
	bitRankW32Int *bitRank1;	//stores bit array without leaves
	bitRankW32Int *bitRank2;	//stores bit array to navigate
	uint *array_sum;			//stores the leaves and all sums
	uint *array_sum2;			//stores the leaves and all sums when the array has had regions of zeros
	uint regenerate;			//flag to identificate if it has to use array_sum or array_sum2
	uint total_sum;				//total sum of the cells in arrays and bitrank2
	uint remove_leaves;			//leaves that have to be removed

}bitmap;


//struct to store the tags, parents and offsets per level
typedef struct warehouse{
	char **tags;				//stores the tags from input file
	uint *value;				//stores the value of the data cube from input file
	uint **ancestor;			//stores the ancestor per level and dimension
	uint **offset;				//stores the offset per level and dimension
}data;


//struct for implementing the formula and obtain offsets
typedef struct form{
	uint *rank0s;				//number of rank0 before the position in louds
	uint *select0s;				//select0 to obtain the position in louds
	uint **ones;				//number of ones between zeros
	uint **result;				//results to obtain the offsets
}formel;



//calculates the time of construction
double get_process_time();


//cleans the character \n of each line
void cleanCharacter(char *string);


//reads the vocabulary of each dimension to obtain the vocabulary and the dimension size
void read_vocabulary(char *filein, dim *dimension);


//reads the bitmap of each dimension to get the vars: bitmap, levels, total_elements_per_level and elements_per_level
void read_dimension(char *filein, dim *dimension);


//prints the vocabulary of one dimension
void print_vocabulary(dim *dimension);


//prints the bitmap of one dimension
void print_bitmap(bitRankW32Int *bitmap);


//prints the elements_per_level of one dimension
void print_elements_per_level(dim *dimension);


//prints the final dimension with all the sizes
void print_final_dimension(dim **array_dims, uint numDims);


//frees the allocated memory
void free_dims(dim **array_dims, uint numDims);


//calculates the sum from the parameter level, i.e. for level 2, it calculates the total sum of sizes from level 1 and level 2
uint calculate_level(dim *array_dim_result,uint level);


//multiplies all the levels between two dimensions
uint *multiply_dims(dim *array_dims1,dim *array_dims2, uint levels,uint *total_elements_in_level, uint rest);


//gets the final dimension as a result of multiplying all existing dimensions
void all_dimensions(dim **array_dims, uint numDims);


//performs exponentiation recursive
uint power(int x,int y);


//allocates memory for total elements and elements per level (array_dim_result)
void alloc_array_dim_result(dim **array_dims, uint numDims);


//reads the bitmaps of all dimensions and stores them into matrices to use them in multiDims
void readBitmaps(char **argv, int numDimensions,bitRankW32Int ***bitmaps);


//multiplies recursively all dimensions and levels to obtain the final structure
void multiDims(uint numDimensions, uint numLevels, uint *starPosition, uint *endPosition, uint level, bitRankW32Int ***bitmaps, dim *array_dim_result);


//fills the last level of total_elements_per_level
void fill_last_level_array_dim_result(dim *array_dim_result);


//initializes hash table
void initialize2(unsigned long sizeVoc);


//loads the vocabulary from all dimensions and assigns the positions of the hash table corresponding to LOUDS
void hash_table(dim **array_dims, uint numDims);


//reads the input file with tags and values
void read_input_file(char *argv[], data *datawarehouse, uint total, uint numDims);


//creates the bitmaps1 and 2, the bitrank2 and the array_sum
void create_bitmaps(bitmap *bitmap_dims, dim *array_dim_result);


//looks for the ancestor at every level in every dimension
void look_for_ancestor(dim **array_dims, data *datawarehouse, uint total, uint numDims);


//calculates the offset of each level using the formula
void calculate_offset(formel *formula, dim **array_dims, data *datawarehouse, uint total, uint numDims);


//fills the leaves of the CMHD with the read input file
void fill_leaves(bitmap *bitmap_dims, data *datawarehouse, dim **array_dims, uint total, uint numDims);


//sums cells to obtain the aggregated sum at upper levels
void sum_dw_tree(dim *array_dim_result,bitmap *bitmap_dims);


//removes the zeros of CMHD and generates new bitmaps
void remove_zeros(dim *array_dim_result, bitmap *bitmap_dims,char *argv);


//initializes dimensions (reads vocabulary & dimension), moreover creates the array_dim_result (array_dims[numDims])
void initialize(dim **array_dims,char *argv[],uint numDims);


//frees the allocated memory for the dw
void free_datawarehouse(data *datawarehouse, uint total, uint numDims);


//frees the allocated memory for the formula
void free_formula(formel *formula, uint numDims);


//frees the allocated memory for dimensions
void free_bitmap_dims(bitmap *bitmap_dims);


//saves array_sum or array_sum2 and bitrank1 & 2
void save_bitmaps(bitmap *bitmap_dims, FILE *f, dim **array_dims,uint numDims);


/*saves vocabulary, louds & dimension_size from all dimensions, 
moreover saves the number of levels, total elements per level and elements per level from array_dim_result (all merged dimensions) */
void save_dims(dim **array_dims, FILE *f, uint numDims, uint total);


//saves ancestor and offset from dw
void save_warehouse(data *datawarehouse, FILE *f, uint total, uint numDims);


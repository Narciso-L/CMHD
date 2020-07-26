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


//#include "dwtree.h"

#include "louds.h"
#include "hash.h"
#include <ctype.h>
#include "dacs.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>


#define BUFFER_SIZE 2048
#define SONS_PER_DIMENSION 64   	//number of childs per dimension



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
	bitRankW32Int *bitRank1;	//stores bit array without leaves
	bitRankW32Int *bitRank2;	//stores bit array to navigate
	uint *array_sum;			//stores the leaves and all sums
	uint *array_sum2;			//stores the leaves and all sums when the array has had regions of zeros
	uint regenerate;			//flag to identificate if it has to use array_sum or array_sum2
	uint total_sum;				//total sum of the cells in arrays and bitrank2

}bitmap;


//struct to store the tags, parents and offsets per level
typedef struct warehouse{
	char **tags;			//stores the tags from input file
	uint **ancestor;		//stores the ancestor per level and dimension
	uint **offset;			//stores the offset per level and dimension
}data;


//struct for implementing the formula and obtain offsets
typedef struct form{
	uint *rank0s;		//number of rank0 before the position in louds
	uint *select0s;		//select0 to obtain the position in louds 
	uint **ones;		//number of ones between zeros
	uint **result;		//results to obtain the offsets
}formel;


//struct to determine the children of a parent
typedef struct pos{
	uint first;			//first child
	uint last;			//last child
	uint num_childs;	//number of childs all levels
	uint cant_childs;	//number of childs per level
	uint total_childs;	//number of childs in previous levels
	uint **positions;	//child's positions per level & dimension in louds
}son;



/*loads vocabulary, louds & dimension_size from all dimensions, 
moreover loads the number of levels, total elements per level and elements per level from array_dim_result (all merged dimensions) */
void load_dims(dim **array_dims, uint numDims,char *path);


//loads array_sum or array_sum2 and bitrank1 & 2
void load_bitmaps(bitmap *bitmap_dims,char *path);


//initializes hash table
void init_hash(unsigned long sizeVoc);


//loads the vocabulary from all dimensions and assigns the positions of the hash table corresponding to LOUDS
void hash_table(dim **array_dims, uint numDims);


//cleans the character \n of each line
void cleanCharacter(char *string);


//reads the query input file and stores the tags in datawarehouse->tags
void read_query_input_file(char *argv[], data *datawarehouse, uint numQueries, uint numDims);


//frees the formula's allocated memory
void free_formula(formel *formula, uint levels);


//frees the offset's allocated memory
void free_offset(data *datawarehouse, uint levels);


//frees the ancestor's allocated memory
void free_ancestor(data *datawarehouse, uint levels);


//frees the datawarehouse's allocated memory
void free_datawarehouse(data *datawarehouse, uint numQueries, uint numDims);


//frees array_sum or array_sum2 and bitrank1 & 2
void free_bitmap_dims(bitmap *bitmap_dims, uint regenerate);


//frees the dimensions' allocated memory
void free_dims(dim **array_dims, uint numDims);

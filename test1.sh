# Copyright (C) <2016>  <Narciso López-López>

    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.

    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <https://www.gnu.org/licenses/>.



#!/bin/bash

ENCODE_EXE=./Busqueda/search
OUTPUT_FOLDER=./Busqueda/
QUERIES_FOLDER=./Queries/

dataset_dir="$(dirname "$0")"
source $dataset_dir/dataset.sh

#set -e
for i in `seq 0 $(( ndatasets - 1 ))`
do
    dataset=${datasets[$i]}
    numDimensions=${dims[$i]}
    for z in `seq 0 $(( nZeros - 1 ))`
	do
		zero=${zeros[$z]}
	   	echo ""
	    echo "*********************************************************************"
	    echo ""
	    echo "-----------------------$dataset-----Zeros: $zero -------------------"

		################################
		##### BINARY               #####
		################################
		echo "-------------------- Binary------------------------------------------"
		echo "---------- Ultime ----------"
		$ENCODE_EXE $numDimensions $numOfQueries $QUERIES_FOLDER/queries_$dataset"_B"_ultime.txt $OUTPUT_FOLDER/$dataset"_"$zero"_B"  
		#$ENCODE_EXENarciso $numDimensions $numOfQueries $QUERIES_FOLDER/queries_$dataset"_B"_ultime.txt $OUTPUT_FOLDER/Narciso/$dataset"_B_"$zero

		################################
		##### RANDOM                #####
		################################
		echo "--------------------- Random------------"
		echo "---------- Ultime ----------"
		$ENCODE_EXE $numDimensions $numOfQueries $QUERIES_FOLDER/queries_$dataset"_"ultime.txt $OUTPUT_FOLDER/$dataset"_"$zero  
		#$ENCODE_EXENarciso $numDimensions $numOfQueries $QUERIES_FOLDER/queries_$dataset"_"ultime.txt $OUTPUT_FOLDER/Narciso/$dataset"_"$zero
	done
 done

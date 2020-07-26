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

ENCODE_EXE=./Construccion/builddwtree
STRUCTURE_FOLDER=./Structure/ 
DATA_FOLDER=./Data/     
OUTPUT_FOLDER=./Busqueda/

dataset_dir="$(dirname "$0")"
source $dataset_dir/dataset.sh

set -e
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


		########################################
		##### Binary dataset -- CMHD B #####
		########################################
		INPUT_DATA=""
		for d in `seq 1 $numDimensions`
		do
			INPUT_DATA=" $INPUT_DATA $STRUCTURE_FOLDER/$dataset""_B""/dim$d " 
		done

			for d in `seq 1 $numDimensions`
		do
			INPUT_DATA=" $INPUT_DATA $STRUCTURE_FOLDER/$dataset""_B""/bitmap$d " 
		done

		echo ""
		echo ""
		echo "-----------$dataset""_B""------ CMHD ------------"
		mkdir -p $OUTPUT_FOLDER/$dataset"_"$zero"_B"/
		$ENCODE_EXE $numDimensions $DATA_FOLDER/$dataset"_"$zero.txt $INPUT_DATA $OUTPUT_FOLDER/$dataset"_"$zero"_B"

		########################################
		##### RANDOM dataset -- CMHD #####
		########################################
		INPUT_DATA=""
		for d in `seq 1 $numDimensions`
		do
			INPUT_DATA=" $INPUT_DATA $STRUCTURE_FOLDER/$dataset/dim$d " 
		done

			for d in `seq 1 $numDimensions`
		do
			INPUT_DATA=" $INPUT_DATA $STRUCTURE_FOLDER/$dataset/bitmap$d " 
		done


		echo "-----------$dataset------ CMHD ------------"
		mkdir -p $OUTPUT_FOLDER/$dataset"_"$zero/
		$ENCODE_EXE $numDimensions $DATA_FOLDER/$dataset"_"$zero.txt $INPUT_DATA $OUTPUT_FOLDER/$dataset"_"$zero
	done
done

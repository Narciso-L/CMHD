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

OUTPUT_FOLDER=./Busqueda/

dataset_dir="$(dirname "$0")"
source $dataset_dir/dataset.sh

for i in `seq 0 $(( ndatasets - 1 ))`
do
    dataset=${datasets[$i]}

    echo "-----------$dataset""_B""------ FSMD -------------"
	OUTPUT_FOLDER_GRAPH="$OUTPUT_FOLDER/$dataset""_B/"
	du -sb $OUTPUT_FOLDER_GRAPH

	echo "-----------$dataset-------- FSMD ------------------"
	OUTPUT_FOLDER_GRAPH="$OUTPUT_FOLDER/$dataset/"
	du -sb $OUTPUT_FOLDER_GRAPH
done

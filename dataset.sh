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

ndatasets=13
datasets=(3D_8 3D_16 3D_32 4D_8 4D_16 4D_32 5D_8 5D_16 5D_32 6D_8 6D_16 7D_8 8D_8)
dims=(3 3 3 4 4 4 5 5 5 6 6 7 8)

#ndatasets=2
#datasets=(4D_16 4D_32)
#dims=(4 4)


# ZEROS
nZeros=3
zeros=(0 25 50)


# QUERIES
minValue=1
maxValue=10
numOfQueries=10000

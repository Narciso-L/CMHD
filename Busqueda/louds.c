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


//finds position of the first child for node at the position i
uint first_child(bitRankW32Int *bitmap, uint pos){

    return select0(bitmap,rank(bitmap,pos))+1;
}


//finds position of the last child for node at the position i
uint last_child(bitRankW32Int *bitmap, uint pos){

    return select0(bitmap,rank(bitmap,pos)+1)-1;
}


//finds position of the parent for the node at the position i
uint parent(bitRankW32Int *bitmap, uint pos){

	uint rank0 = pos-(rank(bitmap,pos-1));
    return select1(bitmap,rank0);
}


//return number of children for node at the position i
uint children(uint first, uint last){

    return (last - first)+1;
}


//returns position of num-th child for the node at the position i, num >= 0
uint next_child(uint first, uint num){

     return first + num;
}

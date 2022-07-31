# sk_is18
is there a 18 clues in a sudoku solution grid

This file created in April 2021 is a stand alone program preparing a function to be called in C C++ to check whether a given sudoku solution grid can be solved by a 18 clues sudoku.

"blue" an old member of the New Sudoku Players forum supplied long ago such a function for a windows 10 system, but the code is not known.
This code is somehow derived from "blue's" general view on the process to apply. This code is delivered to give a possibility to run a LINUX application doing the job.

The code first version has been tested in April/May 2021.

Benchmarks showed that blue's code was much better. The third version expected to be ready in summer 2022 tries to approach blue's performance.

The Unavoidable sets (UA) generator has been completely rewritten, using a mixed funciton combining a "brute force generator" and the previously seen unavoidable sets. 
The generation is now done in some milliseconds compared to several hundred milliseconds with the previous code.

But the main change is in the expansion of unavoidable sets.

Here, the bands are ordered to have at the top bands with the lowest count of valid one band solutions.

Unavoidable sets used are 
. all bands and stack UAs (internal table)
. all 2 digits UAs
. all 3;4;5 digits UAs for bands 1+2
. all UAs solving 4 boxes in bands 1+2
. all UAs with 2 cells in band 3 with 3 4 5 digits

The program directly expand th UAs within bands 1+2 to get valid solutions for bands 1+2.
If the solution is not valid, a fresh UA is used and stored
If the solution is valid, an attempt to find a 18 adding clues in band 3 is done.

A special process is applied if there is room in bands 1+2 for more clues.

In this third vresion, the task is covered usin six permutations each band/stack being once reordered as band 3.
This is very close to the process applied in fact for the 17 clues  scan looking for all 17.
In this process, most of the redundancy is filtered. The distribution 666 is searched once (out of the 6 permutations) ,666 both in bands and stacks.

Generally speaking, the costly final check that the band 1+2 proosed out of the uas list is postponed to the first proposed 18. 

<<<<<<<<<<<<<<<<<<<<<<<<<<<< Compiling the file
The file names are given in line with g++ 

a .cpp file is separate module for compilation
a .h module can contain some code, but then is attached ot a .cpp module

the main.cpp module is there to handle the command line and to read the entry ile, celling the function

expected .cpp files are

main.cpp 
go-12sol-tables.cpp containing tables specific to the processing of solution grids and ED bands
sk-bitfields.cpp some code for bit field of different sizes (mainly 128 bit fields)
sk_t.cpp code shared in many applications
zhb1b2b.cpp brute force solver limited ot 2 bands
zhn.cpp brute force solver for the entire solution grid

sk-is18.cpp the function main file


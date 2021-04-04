# sk_is18
is there a 18 clues in a sudoku solution grid

This file created in April 2021 is a stand alone program preparing a function to be called in C C++ to check whether a given sudoku solution grid can be solved by a 18 clues sudoku.

"blue" an old member of the New Sudoku Players forum supplied long ago such a function for a windows 10 system, but the code is not known.
This code is somehow derived from "blue's" general view on the process to apply. This code is delivered to give a possibility to run a LINUX application doing the job.

The code is expected to be fullt tested in April/May 2021.

This version has an external code reading a file of the solution grids to test, calling the future function.

As this is a test version, a pre processor option is available to process a file containing solution grids with a known 18. Then, the 18 to search is located in positions 83... of the entry file.

The contains contains many other debugging options, generally driven by pre processor options
The current code does not process solutions grids having one or more bands with a 2 clues valid band (band 29 in the list 1-416). If such a grid is sent, the return code is -1. A specific code will come later for this specific case known to have a very limited number of valid ED solution grids well identified. 

The lay-out of the function (for a given solution grid) is the following :

a) preliminary tasks:
  identify the 3 bands number (0-415) and morph
  expand the 3 bands for all valid bands of the gangster having 3 clues to 7 clues.
    This is done morphing an internal table of all UAs of each band having no subset.
  Sort the band in increasing order of the number of valid 6 clues of the band
  call one of three permutations of the bands where the 18 searched will have the highest number of clues in band 3
    Only the first permutation will search the worst case, the 6;6;6 distribution of clues
    
b) processing one permutation

b.1 preliminary tasks 
    the preliminary task is the collection of unavoidable sets (UAs) for the permutations
    (b1a) . UAs specific to bands 1+2 (excluding uas located in one band)
    (b1b) . UAs having a limited number of clues in bands 3 (basically a socket 2 clues or 3 clues with bands 1+2 or a double socket 2x2)
    
b.2 Selecting valid band1;valid band2 pairs hitting all UAs (b1)
    this is a critical process where the number of pairs can be billions.
    
    The "internal loops" is called a step. In a step 
     all valid band 1 have 2 or more common celss 
     all valid bands 2 have 2 or more common cells
      (in fact, such an index is created at the band expansion in the preliminary tasks)
     An "outer loop" is added usig short uas of (b1a).
     
 b.3 a pair passing this filter is pre checked  against the maximum number of clues expected in band 2  
     the validity of the pairs solution for bands 1+2 is controlled (possible new UA (b1a))
     
 b.4 for each pair passing all controls
     the analysis on UAs still valid in band 3 is done with specific cases 
       depending on the number of clues forced by small residual uas 
    
 The process is halted as soon as possible when the first 18 clues appears
 
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


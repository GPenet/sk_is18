
#define _CRT_SECURE_NO_DEPRECATE

#define TEST_ON
#define HAVEKNOWN 1

//#define COLOIN


#define GUALIMSIZE 18
#define TUASIZE 5000

#define UA32_10 0xffc00000
#define UA64_54 0x3fffffffffffff
//============================================== 
/* all bands seen
maxindex= 983
maxn5= 51516
maxn6= 237762
maxdet5= 261
maxdet6= 2004
*/

#define G17MORESIZE 128

#include <sys/timeb.h>
#include "main.h"  // main and main tables and basic brute force
#include "go_17sol_tables.h"     
#include "Zh1b2b.h"  // brute force 2 bands  
//_________________ brute force handling 1 2 3 bands 
extern SGO sgo;
extern ZHOU    zhou[50],zhou_i;// , zhou_i, zhou_solve;
extern ZH_GLOBAL zh_g;
extern ZH_GLOBAL2 zh_g2;
extern ZHGXN zhgxn;
extern ZHOU2 zhou2[5];
extern ZHOU3 zhou3[10];

extern ZH2B_GLOBAL   zh2b_g;
extern ZH2B5_GLOBAL   zh2b5_g;   // 2_5 digits 2 bands
extern ZH2B_1D_GLOBAL zh2b1d_g;  // one digit 2 bands
extern ZH2GXN zh2gxn;
extern ZH2_3  zh2_3[10];
extern ZH2_4  zh2_4[20];
extern ZH2_5  zh2_5[20];

extern ZH2B zh2b[40], zh2b_i, zh2b_i1;
extern ZH2B5 zh2b5[10]; // solved digit per digit
extern ZH2B_1D zh2b1d[6]; // maxi 6 guesses, one per row
extern ZHONE_GLOBAL   zh1b_g;
extern ZHONE zhone[20];
extern ZHONE zhone_i;
extern ofstream  fout1;
//FINPUT finput;

#include "sk_is18.h" 
GCHK gchk;
TUAS81 tuas81;
TUASB12 tuasb12;
STD_B416 bax[6],baxs[3];// the  bands/stacks at start  
//STD_B1_2 myband1, myband2;
uint64_t  p_cpt2g[70];

uint64_t to_clean[200];

#include "go_17sol_tables.h"
#include "go_18_bands_cpp.h" 
#include "go_chkis18_bs_cpp.h"  
#include "go_chk18_commands_cpp.h"





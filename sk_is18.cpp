
#define _CRT_SECURE_NO_DEPRECATE
/* program organisation
	main is the standard frame including the basic brute force
*/

#define SEARCH17SOL


//#define DEFPHASE -2
#ifdef DEFPHASE
#endif
//#define DEBUGKNOWN
#ifdef DEBUGKNOWN
#endif

//#define DEBUGONE
#ifdef DEBUGONE
#endif

//#define DEBUGINIT
#ifdef DEBUGINIT
#endif

//#define DEBUGEXL
#ifdef DEBUGEXL
#endif

//#define DEBUGL1L2 95
#ifdef DEBUGL1L2
#endif

//#define DEBUGSTEP 154643
#ifdef DEBUGSTEP
#endif

//#define DEBUCLEANB 160656673
#ifdef DEBUCLEANB
#endif
//#define DEBUGPAT4
#ifdef DEBUGPAT4
#endif

//#define NEWPAT4
#ifdef NEWPAT4
#endif
//#define GTEST17_ON 1
/* control of the UAs generation
	UALIMSIZE to test the range 18-20
	TUASIZE table for all uas 
	TUASIZE_ONE max for one call fo the
	  brute force generator
*/
#define XCHUNK64 200
#define YCHUNK64 200
#define XCHUNK128 200
#define YCHUNK128 200
#define XCHUNK256 100
#define YCHUNK256 100

#define UALIMSIZE 22
#define GUALIMSIZE 18
#define TUASIZE 5000
#define TUASIZE_ONE 50

#define UA32_10 0xffc00000
#define UA64_54 0x3fffffffffffff
#define TUA64_12SIZE 2000
#define ADD_N128 4
//============================================== 
/* all bands seen
maxindex= 983
maxn5= 51516
maxn6= 237762
maxdet5= 261
maxdet6= 2004
*/
#define MAXN5 51520
#define MAXN6 237770 
#define MAXN7 820000 
#define MAXSTEP5 5000
#define MAXSTEP6 23000 
#define MAXNIND5 300
#define MAXNIND6 2100
#define MAXNIND7 40000

#define G17MORESIZE 32

#define G17TESTUASGUASLIMITS 1

#include <sys/timeb.h>
#include "main.h"  // main and main tables and basic brute force
#include "go_17sol_tables.h"     
#include "Zh1b2b.h"  // brute force 2 bands  
//_________________ brute force handling 1 2 3 bands 
extern SGO sgo;
extern ZHOU    zhou[50],zhou_i;// , zhou_i, zhou_solve;
extern ZH_GLOBAL zh_g;
extern ZH_GLOBAL2 zh_g2;
extern ZH2B_GLOBAL   zh2b_g;
extern ZH2B5_GLOBAL   zh2b5_g;   // 2_5 digits 2 bands
extern ZH2B_1D_GLOBAL zh2b1d_g;  // one digit 2 bands
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

STD_B416 bax[3];// the 3 bands at start
STD_B1_2 myband1, myband2;
uint64_t  p_cpt2g[70];

// working areas for external loop
//X base band can be band1 or band2   Y is the other band
//   X is split in yes and others
//   X yes is processed with Y
//   then Y is reduced to yes
BI2_32 bi2_b1w[250], bi2_b1yes[250];
BI2_32 bi2_b2w[250], bi2_b2yes[250];
BI2_32 bi2_b1w2[250], bi2_b1yes2[250];
BI2_32 bi2_b2w2[250], bi2_b2yes2[250];
#define MAXEXP7 1200000
VALIDB vab1w[MAXEXP7], vab1yes[MAXEXP7]; 
VALIDB vab2w[MAXEXP7], vab2yes[MAXEXP7]; 
VALIDB vab1w2[MAXEXP7], vab1yes2[MAXEXP7];
VALIDB vab2w2[MAXEXP7], vab2yes2[MAXEXP7];

VALIDB64 vab64b1[MAXEXP7], vab64b2[MAXEXP7];

ZS128 zs128b1[MAXEXP7], zs128b2[MAXEXP7];;

uint64_t to_clean[2000000];

GENUAS_B12 genuasb12;
GEN_BANDES_12 genb12;



#include "go_17sol_tables.h"
#include "go_18_bands_cpp.h" 
#include "go_18_genb12_cpp.h"     
#include "go_chkis18_bs_cpp.h"  
//#include "zh5_cpp.h"
#include "go_chk18_commands_cpp.h"





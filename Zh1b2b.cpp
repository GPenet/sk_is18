/*
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/
#include "main.h"   
#include "Zh1b2b.h"
extern uint64_t  p_cpt2g[70];

// row unsolved 6 bits per digit lookup table
uint64_t zh2b_t_runsolved[9] = { 077 , 077 << 6 , 077 << 12 , 077 << 18 , 077 << 24 ,
(uint64_t)077 << 32 ,(uint64_t)077 << 38 ,(uint64_t)077 << 44 ,(uint64_t)077 << 50 };
uint32_t zh2b_t_runsolvedshift[9] = { 0,6,12,18,24,32,38,44,50 };
//========================= global variable and working areas 
extern ZH_GLOBAL zh_g;
ZH2B zh2b_i, zh2b_i1,zh2b[40] ;
ZH2B_GLOBAL   zh2b_g;   // 2 bands 9 digits
ZH2GXN zh2gxn;
ZH2_3  zh2_3[10];
ZH2_4  zh2_4[20];
ZH2_5  zh2_5[20];

ZH2B5_GLOBAL   zh2b5_g;   // 2_5 digits 2 bands
ZH2B_1D_GLOBAL zh2b1d_g;  // one digit 2 bands

ZH2B_1D zh2b1d[6]; // maxi 6 guesses, one per row
ZH2B5 zh2b5[10]; // solved digit per digit

ZHONE_GLOBAL zh1b_g;
ZHONE zhone[20]; // one band 9 digits 
ZHONE zhone_i; 


//#include "go_17sol_zx_UaCollector_cpp.h"
//GENUAS_1Bx genuas1b;
//============================= ZH_GLOBAL code
ZH2B_GLOBAL::ZH2B_GLOBAL(){
	zsol = 0; // no solution unless required buy the user
	nctlg =  0;

}
void ZH2B_GLOBAL::GetBands(int * g1, int * g2) {
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = puz0[cell];
			fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = puz0[cell+27];
			fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
}

void ZH2B_GLOBAL::InitGangster(int * g0, int * g1) {
	uint64_t col64 = (uint64_t)Zhoucol | ((uint64_t)Zhoucol << 32);
	memcpy(fd_revised, fd_sols[1], sizeof fd_revised);
	for (int i = 0; i < 9; i++, col64 <<= 1) {
		if (g0[i] == g1[i])continue;
		int changes = g0[i] ^ g1[i]; // one added one cleared
		for (int d = 0, bit = 1; d < 9; d++, bit <<= 1) {// check digits
			if (!(changes & bit)) continue;
			if (g0[i] & bit)fd_revised[d] &= ~col64;
			else fd_revised[d] |= col64;
		}		
	}
}


uint64_t ZH2B_GLOBAL::BuildUaret(BF64 * wsol) {
	ua_ret = 0;
	for (int i = 0; i < 9; i++) {
		BF64 w = wsol[i] - fd_sols[0][i];
		ua_ret |= w.bf.u64;
	}
	return ua_ret;
}


//============= zh2B code for uas 2 bands and validity 2 bands puzzle
/*
void ZH2b::Start_nFloors(int floors, BF64 * mypm) {
	cells_unsolved;bf.u64=BIT_SET_2X^ solved_cells;
	for (uint32_t i = 0; i < nd; i++) {
		myfd[i] &= cells_unsolved;
		if (myfd[i].count == 6)return ;
	}
	
}*/
//============================================ ZH2B code

void ZH2B::Init_gang() {//init after zh1b_g InitGangster
	zh2b_g.ndigits = 9;
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_revised, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}

void ZH2B::DebugSol() {//init after zh1b_g getband
	zh2b_g.ndigits = 9;
	memset(this, 0, sizeof zh2b[0]);
	memcpy(FD, zh2b_g.fd_sols[0], sizeof FD);
	ImageCandidats();
}

#define UPDN(I,J)A=FD[I].bf.u32[J];\
Shrink = (TblShrinkMask[A & 0x1FF] | \
TblShrinkMask[ (A>>9) & 0x1FF]<<3 | \
TblShrinkMask[ (A>>18) & 0x1FF]<<6);\
if ((A &=TblComplexMask[Shrink]) ==0)  return 0; \
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF); \
FD[I].bf.u32[1 - J] &= TblMaskSingle[S]; \
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]]; \
CompFD[I].bf.u32[J] = FD[I].bf.u32[J] = A;

#define UPWCL(I,P,Q,R,T,U,V,W,X)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;\
wcl[U]&= cl;wcl[V]&= cl;wcl[W]&= cl;wcl[X]&= cl;

#define UPWCL3(I,P,Q)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;

#define UPWCL4(I,P,Q,R)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;

#define UPWCL5(I,P,Q,R,T)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;


int ZH2B::ApplySingleOrEmptyCells() {
	zh2b_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map = map[2],
		R3 = R2 & Map, R4;// digits 12
	R2 |= R1 & Map; R1 |= Map;

	Map = map[3]; R4 = R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[4];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[5];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[6];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[7];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[8];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		//if (zh2b_g.diag) cout << Char2Xout(R1) << " apply R1 " << endl;
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 9; idig++) {
				if (map[idig] & bit) {// this is the digit
					int cell = From_128_To_81[res];
					Assign(idig, cell, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		R2 &= ~R3;
		R3 &= ~R4;
		uint32_t res, cell;
		if (zh2b_g.diag) 	cout<<Char2Xout(R2) << "R2  count "<< _popcnt64(R2) << endl;
		if (!R2) {
			if (!R3) R3 = R4;
			bitscanforward64(res, R3);
			zh2b_g.guess_xcell = res;
			return 0;
		}
		// try to get 2 cells or more
		while (bitscanforward64(res, R2)) {
			zh2b_g.guess_xcell = res;
			cell = From_128_To_81[res];
			uint64_t mask=R2 & cell_z3x[cell].u64[0];
			if (_popcnt64(mask)) return 0;
			if (_popcnt64(R2) < 3)return 0;
			R2 ^=( uint64_t)1<< res;
		}
		return 0;
	}
}
int ZH2B::FullUpdate() {
	//if (zh2b_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}


//================================ZH2B code
/*
*/


inline int ZH2B::Seta(int digit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell],
		block = TblBoard_Block[cell];
	if (FD[digit].Off(xcell)) return 1; // not valid
	Assign(digit, cell, xcell);
	BF64 *Fd = &FD[digit];
	*Fd &= AssignMask_Digit[cell].u64[0];
	int ddig = 6 * digit;
	if (digit > 4) ddig += 2;// second bloc of 32 bits
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	cells_unsolved.Clear(xcell);
	BF64 * RF = &FD[8];
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}

int ZH2B::Update(){
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink ){
	  Shrink = 0;
	  if (!rows_unsolved.bf.u32[0])goto digit5;

	  {register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
	  if (!(AR & 077))goto digit1;

//=digit 0
	  if (FD[0].bf.u32[0] ==  CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0,0)	if ((AR & 7) != S){
				AR &= 07777777770 | S;	UPWCL(0, 2,4,6,8,10,12,14,16)	}

digit0b:if (FD[0].bf.u32[1] ==CompFD[0].bf.u32[1])goto digit1;
		UPDN(0,1)	if (((AR >> 3) & 7) != S){
				AR &= 07777777707 | (S << 3);	UPWCL(1, 3,5,7,9,11,13,15,17)	}

digit1:	if (!(AR  & 07700))goto digit2;

		if (FD[1].bf.u32[0] ==	CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1,0)	if (((AR >> 6) & 7) != S){
			AR &= 07777777077 | (S << 6);UPWCL(0, 0, 4, 6, 8, 10, 12, 14, 16)	}

digit1b:if (FD[1].bf.u32[1] ==	CompFD[1].bf.u32[1])goto digit2;
		UPDN(1,1)		if (((AR >> 9) & 7) != S){
			AR &= 07777770777 | (S << 9);UPWCL(1, 1, 5, 7, 9, 11, 13, 15, 17)	}

digit2:	if (!(AR  & 0770000))goto digit3;

		if (FD[2].bf.u32[0] ==	CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2,0)	if (((AR >> 12) & 7) != S){
			AR &= 07777707777 | (S << 12);	UPWCL(0, 0, 2, 6, 8, 10, 12, 14, 16)}

digit2b:if (FD[2].bf.u32[1] ==	CompFD[2].bf.u32[1])goto digit3;
		UPDN(2,1)	if (((AR >> 15) & 7) != S){
			AR &= 07777077777 | (S << 15);	UPWCL(1, 1, 3, 7, 9, 11, 13, 15, 17)	}

digit3: if (!(AR & 077000000))goto digit4;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		  UPDN(3,0)	  if (((AR >> 18) & 7) != S){
			  AR &= 07770777777 | (S << 18); UPWCL(0, 0, 2, 4, 8, 10, 12, 14, 16)	 }

digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
		  UPDN(3, 1)if (((AR >> 21) & 7) != S){
			  AR &= 07707777777 | (S << 21); UPWCL(1, 1, 3, 5, 9, 11, 13, 15, 17)	}

digit4:if (!(AR & 07700000000))goto end01234;

		if (FD[4].bf.u32[0] ==	CompFD[4].bf.u32[0])goto digit4b;
		UPDN(4, 0)if (((AR >> 24) & 7) != S){
			AR &= 07077777777 | (S << 24);  UPWCL(0, 0, 2, 4, 6, 10, 12, 14, 16)	}

digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto end01234;
		UPDN(4, 1)if (((AR >> 27) & 7) != S){
			AR &= 0777777777 | (S << 27);  UPWCL(1, 1, 3, 5, 7, 11, 13, 15, 17)	}

end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


digit5:
	  if (!rows_unsolved.bf.u32[1])continue; // second lot  4 digits

	  {register unsigned int  AR = rows_unsolved.bf.u32[1];// valid for digits 5,6,7,8
	  if (!(AR & 077))goto digit6;

	  if (FD[5].bf.u32[0] ==	CompFD[5].bf.u32[0])goto digit5b;
		UPDN(5, 0) if ((AR & 7) != S){
			AR &= 07777777770 | S;	UPWCL(0, 0, 2, 4, 6, 8, 12, 14, 16)		}

digit5b:if (FD[5].bf.u32[1] == CompFD[5].bf.u32[1])goto digit6;
		UPDN(5, 1) 	if (((AR >> 3) & 7) != S){
			AR &= 07777777707 | (S << 3); UPWCL(1, 1, 3, 5, 7, 9, 13, 15, 17)	}

digit6:  if (!(AR & 07700))goto digit7;

		if (FD[6].bf.u32[0] ==  CompFD[6].bf.u32[0])goto digit6b;
		UPDN(6, 0) if (((AR >> 6) & 7) != S){
			AR &= 07777777077 | (S << 6);	UPWCL(0, 0, 2, 4, 6, 8, 10, 14, 16)	  }

digit6b: if (FD[6].bf.u32[1] == CompFD[6].bf.u32[1])goto digit7;
		UPDN(6, 1) if (((AR >> 9) & 7) != S){
			AR &= 07777770777 | (S << 9); UPWCL(1, 1, 3, 5, 7, 9, 11, 15, 17)  }

digit7:   if (!(AR  & 0770000))goto digit8;

		if (FD[7].bf.u32[0] ==			  CompFD[7].bf.u32[0])goto digit7b;
		  UPDN(7, 0)if (((AR >> 12) & 7) != S){
			  AR &= 07777707777 | (S << 12);  UPWCL(0, 0, 2, 4, 6, 8, 10, 12, 16)  }

digit7b:  if (FD[7].bf.u32[1] ==  CompFD[7].bf.u32[1])goto digit8;
		  UPDN(7, 1) if (((AR >> 15) & 7) != S){
			  AR &= 07777077777 | (S << 15);  UPWCL(1, 1, 3, 5, 7, 9, 11, 13, 17)  }

digit8:  if (!(AR  & 077000000))goto end5678;

		  if (FD[8].bf.u32[0] ==  CompFD[8].bf.u32[0])goto digit8b;
		  UPDN(8,0)	  if (((AR >> 18) & 7) != S){
			  AR &= 07770777777 | (S << 18);  UPWCL(0, 0, 2, 4, 6, 8, 10, 12, 14)	  }

digit8b: if (FD[8].bf.u32[1] == CompFD[8].bf.u32[1])goto end5678;
		  UPDN(8,1)	  if (((AR >> 21) & 7) != S){
			  AR &= 07707777777 | (S << 21);  UPWCL(1, 1, 3, 5, 7, 9, 11, 13, 15)  }

end5678:rows_unsolved.bf.u32[1] = AR;
	  }// end of validity for AR
  }// end while
#ifdef DIAG
	cout << "end update cycle" << endl;
	Debug(1);
	ImageCandidats();
#endif
  return 1;
}

char * ZH2B::SetKnown(char * zs) {
	strcpy(zs, &empty_puzzle[27]);
	int tdig[9];// build the table of digits
	register int A = rows_unsolved.bf.u32[0];
	tdig[0] = A & 077; A >>= 6; tdig[1] = A & 077; A >>= 6; tdig[2] = A & 077; A >>= 6;
	tdig[3] = A & 077; A >>= 6; tdig[4] = A & 077;
	A = rows_unsolved.bf.u32[1];
	tdig[5] = A & 077; A >>= 6; tdig[6] = A & 077; A >>= 6;
	tdig[7] = A & 077; A >>= 6; tdig[8] = A & 077;

	for (int idig = 0; idig < 9; idig++) {// one digit
		int arowsj = tdig[idig];// 6 rows
		if (arowsj == 077) continue;
		for (int ib = 0; ib < 2; ib++) {// 3 blocs per digit
			int arows = (arowsj >> (3 * ib)) & 7;
			if (arows == 7) continue; // not assigned
			unsigned int band = FD[idig].bf.u32[ib];
			for (int j = 0; j < 3; j++) if (!(arows & (1 << j))) {
				int row = (band >> TblMult9[j]) & 0x1ff;
				uint32_t  irow;
				bitscanforward(irow, row);
				int	cell = Tblstartblock[ib] + TblMult9[j] + irow;
				zs[cell] = idig + '1';
			}
		}
	}
	return zs;
}
void ZH2B::Init_std_bands() {//init after zh1b_g getband
	zh2b_g.ndigits = 9;
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}
inline void ZH2B::InitTclues(uint32_t * tclues, int n) {
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	for (int icell = 0; icell < n; icell++) {
		int cell = tclues[icell], digit = zh2b_g.puz0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map

		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_init[i];
}






void ZH2B::Debug(int all) {
	cout << "DEBUG  nbsol=" << zh2b_g.nsol << " unsolved=" << Unsolved_Count() << endl;
	//	cout << zh1b_g.out27 << " band1 "<<endl;
	char zi[82];  SetKnown(zi);
	cout << zi << " known rows 1_6 digits " << endl;
	if (!all) return;

	cout << "map per digit bands 2 3" << endl;
	for (int ib = 0; ib < 2; ib++) {
		for (int ir = 0; ir < 3; ir++) {
			for (int idig = 0; idig < 9; idig++) {
				unsigned vf = FD[idig].bf.u32[ib];
				unsigned wr = (vf >> (9 * ir)) & 0x1ff;
				for (int k = 0; k < 9; k++) {
					if (wr & (1 << k))		cout << idig + 1;
					else 		cout << ".";
					if (k == 2 || k == 5) 	cout << " ";
				}
				cout << "  ";
			}
			cout << endl; //end of row
		}
		cout << endl; // end of block
	}
	cout << endl; // end of map per digit

}
int ZH2B::GetAllDigits(int cell) {
	int ir = 0, xcell = C_To128[cell];;
	for (int i = 0; i < 9; i++) if (FD[i].On(xcell))ir |= (1 << i);
	return ir;
}
void ZH2B::ImageCandidats() {
	int dig_cells[81]; for (int i = 0; i < 54; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " << endl << endl;
	for (i = 0; i < 9; i++) {  // attention ici i indice colonne
		lcol[i] = 2;    // 2  mini tous chiffres impos�s
		for (j = 0; j < 6; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 6; i++) { // maintenant indice ligne
		if ((i == 3) || (i == 6)) {
			for (int ix = 0; ix < (tcol + 10); ix++)       cout << (char)'-';
			cout << endl;
		}
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell], 
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < 9; id++)if (digs & (1 << id))
				cout << id + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}

//______________4 box uas collector
void ZH2B::InitBands12(int * g0) {
	zh2b_g.ndigits = 9;
	memcpy(zh2b_g.puz0, g0, sizeof zh2b_g.puz0);
	memset(zh2b_g.fd_sols, 0, sizeof zh2b_g.fd_sols);
	// build sol per digit and pm per digit at start
	for (int i = 0; i < 9; i++) {// i column
		for (int j = 0; j < 3; j++) {//j row band1 or row band 2
			int cell = 9 * j + i, dig = zh2b_g.puz0[cell];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = zh2b_g.puz0[cell + 27];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	//ImageCandidats();

}

uint64_t ZH2B::IsValid(uint32_t * tclues, int n,int onlyone) {
	*this = zh2b[0];
	InitTclues(tclues, n);
	zh2gxn.nua = 0;	zh2gxn.uamin = 100;
	zh2gxn.onlyone = onlyone;
	zh2b_g.go_back = 0;
	ComputeNext();
	return zh2gxn.nua;
}
inline void ZH2B::ComputeNext() {
	if (zh2b_g.go_back) return;
	int ir = FullUpdate();
	//cout << "back full update ir="<<ir	<< endl;		
	//ImageCandidats();
	if (ir == 1)GuessValidB12();
	else if (ir == 2) {// solved 
		zh2gxn.uaw = 0;
		for (int i = 0; i < 9; i++) {
			BF64 w = FD[i] - zh2b_g.fd_sols[0][i];
			zh2gxn.uaw |= w.bf.u64;
		}
		if (zh2gxn.uaw) {
			uint64_t cc = _popcnt64(zh2gxn.uaw);
			if (cc < zh2gxn.uamin)zh2gxn.uamin = cc;
			if (cc > zh2gxn.uamin) return;
			zh2gxn.tua[zh2gxn.nua++] = zh2gxn.uaw;
			//cout << Char2Xout(zh2gxn.uaw) << " uaw "
				//<< cc << " " << zh2gxn.nua << endl;
			if (zh2gxn.onlyone)zh2b_g.go_back = 1;
		}
	}
}
void ZH2B::GuessValidB12() {// 
	if (zh2b_g.go_back) return;
	uint32_t xcell = zh2b_g.guess_xcell, cell, digit;
	cell = From_128_To_81[xcell];
	digit = zh2b_g.puz0[cell];
	uint64_t bit = (uint64_t)1 << xcell;
	// true first if possible
	if (FD[digit].bf.u64 & bit) {
		//cout << digit+1 << cellsFixedData[cell].pt << " try true" << endl;
		//cout << "Guess okr " << digit + 1 << cellsFixedData[cell].pt << endl;
		ZH2B * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNext();
		if (zh2b_g.go_back) return;
	}
	// then false 
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig].bf.u64 & bit) {
			//cout << "Guess nokr " << idig + 1 << cellsFixedData[cell].pt << endl;
			ZH2B * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNext();
			if (zh2b_g.go_back) return;
		}
	}

}


void ZH2B::InitB1245(int * g0) {
	zh2b_g.ndigits = 9;
	memcpy(zh2b_g.puz0, g0, sizeof zh2b_g.puz0);
	memset(zh2b_g.fd_sols, 0, sizeof zh2b_g.fd_sols);
	// build sol per digit and pm per digit at start
	for (int i = 0; i < 9; i++) {// i column
		for (int j = 0; j < 3; j++) {//j row band1 or row band 2
			int cell = 9 * j + i, dig = zh2b_g.puz0[cell];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = zh2b_g.puz0[cell + 27];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// assign box 3 6
	memset(zh2b_g.Digit_cell_Assigned_step, 0, sizeof zh2b_g.Digit_cell_Assigned_step);
	for (int row = 0; row < 6; row++) for (int col = 6; col < 9; col++) {
		int cell = 9 * row + col, xcell = C_To128[cell],
			digit = g0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_step[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_step[i];
	//ImageCandidats();

}
void ZH2B::InitB1346(int * g0) {
	zh2b_g.ndigits = 9;
	memcpy(zh2b_g.puz0, g0, sizeof zh2b_g.puz0);
	memset(zh2b_g.fd_sols, 0, sizeof zh2b_g.fd_sols);
	// build sol per digit and pm per digit at start
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = zh2b_g.puz0[cell];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = zh2b_g.puz0[cell + 27];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// assign box 3 6
	memset(zh2b_g.Digit_cell_Assigned_step, 0, sizeof zh2b_g.Digit_cell_Assigned_step);
	for (int row = 0; row < 6; row++) for (int col = 3; col < 6; col++) {
		int cell = 9 * row + col, xcell = C_To128[cell],
			digit = g0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_step[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_step[i];
	//ImageCandidats();

}
void ZH2B::InitB2356(int * g0) {
	zh2b_g.ndigits = 9;
	memcpy(zh2b_g.puz0, g0, sizeof zh2b_g.puz0);
	memset(zh2b_g.fd_sols, 0, sizeof zh2b_g.fd_sols);
	// build sol per digit and pm per digit at start
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = zh2b_g.puz0[cell];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = zh2b_g.puz0[cell + 27];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// assign box 3 6
	memset(zh2b_g.Digit_cell_Assigned_step, 0, sizeof zh2b_g.Digit_cell_Assigned_step);
	for (int row = 0; row < 6; row++) for (int col = 0; col < 3; col++) {
		int cell = 9 * row + col, xcell = C_To128[cell],
			digit = g0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_step[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_step[i];
	//ImageCandidats();

}
int ZH2B::Do4bGo() {// called in zhou3[1]
	zh2gxn.nua = 0;
	zh2b_g.go_back = 0;
	ComputeNext4box();
	return 0;
}


int ZH2B::FullUpdate4box() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells4box())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}
int ZH2B::ApplySingleOrEmptyCells4box() {
	zh2b_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map = map[2],
		R3 = R2 & Map, R4;// digits 12
	R2 |= R1 & Map; R1 |= Map;

	Map = map[3]; R4 = R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[4];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[5];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[6];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[7];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[8];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		//if (zh2b_g.diag) cout << Char2Xout(R1) << " apply R1 " << endl;
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 9; idig++) {
				if (map[idig] & bit) {// this is the digit
					int cell = From_128_To_81[res];
					Assign(idig, cell, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		R2 &= ~R3;
		R3 &= ~R4;
		R2 &= ~R3;
		if (!R2) R2 = R3;
		zh2gxn.rx = R2;
		int ir = GetNetUaCell4box();
		if (ir < 0)return -1; // dead branch
		zh2b_g.guess_xcell = ir;
		return 0;
	}
}
int ZH2B::GetNetUaCell4box() {
	// setup ok cells
	uint64_t sok = 0; // zh2gxn.unsolved_field;
	for (int i = 0; i < 9; i++) {
		sok |= (FD[i].bf.u64 & zh2gxn.fsol[i]);
	}
	uint64_t sok2 = sok & ~cells_unsolved.bf.u64;
	sok &= cells_unsolved.bf.u64;
	int nu = *zh2gxn.nknownuas, cret;
	//cout << Char2Xout(sok) << " sok next cell" << endl;
	//cout << Char2Xout(sok2) << " sok2 next cell" << endl;
		// first non hit ua to solve
	uint64_t u_or = 0; // zh2gxn.unsolved_field;
	for (int iu = 0; iu < nu; iu++) {
		register uint64_t U = zh2gxn.knownuas[iu];
		if (U &sok2) continue;
		//cout << Char2Xout(U) << "U to fill" << endl;
		U &= sok;
		if (!U)return -1; // dead branch
		u_or |= U;
	}
	//cout << Char2Xout(u_or) << " u_or" << endl;
	if (u_or) {// some ua(s) to fill
		if (u_or & zh2gxn.rx) u_or &= zh2gxn.rx;
		bitscanforward64(cret, u_or);
		return cret;
	}
	// all uas solved, first cell in rx
	bitscanforward64(cret, zh2gxn.rx);
	return cret;
}

inline void ZH2B::ComputeNext4box() {
	if (zh2b_g.go_back) return;
	int ir = FullUpdate4box();
	//cout << "back full ir=" << ir << endl;
	//ImageCandidats();
	if (ir == 1)Guess4box();
	else if (ir == 2) {// solved 
		uint64_t ww = 0;
		for (int i = 0; i < 9; i++) {
			BF64 w = FD[i] - zh2b_g.fd_sols[0][i];
			ww |= w.bf.u64;
		}
		if (ww) {
			uint64_t cc = _popcnt64(ww);
			if (cc > 20) return;
			//if (zh2gxn.nua > 40 && cc > 18)return;
			//if (cc > 22) return;
			if (zh2gxn.nua > 40 && cc > 20)return;
			// check no ua false
			int nu = *zh2gxn.nknownuas;
			//cout << Char2Xout(ww) << " seen nu="<<nu << endl;
			if (nu) {
				uint64_t sok = 0; // zh2gxn.unsolved_field;
				for (int i = 0; i < 9; i++) {
					sok |= (FD[i].bf.u64 & zh2gxn.fsol[i]);
				}
				//cout << Char2Xout(sok) << " sok" << endl;
				for (int iu = 0; iu < nu; iu++) {
					register uint64_t U = zh2gxn.knownuas[iu];
					if (!(U &sok)) return;;
				}
			}

			zh2gxn.tua[zh2gxn.nua++] = ww;
			zh2gxn.knownuas[(*zh2gxn.nknownuas)++]	= ww;
			//cout << Char2Xout(ww) << " added zh2gxn.nua="<< zh2gxn.nua 
			//	<< " cc="<<cc<< endl;
			if (zh2gxn.nua >80)zh2b_g.go_back = 1;
		}
	}
}

void ZH2B::Guess4box() {// 
	if (zh2b_g.go_back) return;
	uint32_t xcell = zh2b_g.guess_xcell, cell, digit;
	cell = From_128_To_81[xcell];
	digit = zh2b_g.puz0[cell];
	uint64_t bit = (uint64_t)1 << xcell;
	// true first if possible
	if (FD[digit].bf.u64 & bit) {
		//cout << "guess ok" <<digit+1<< cellsFixedData[cell].pt << endl;
		ZH2B * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNext4box();
	}
	// then false 
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig].bf.u64 & bit) {
			//cout << "guess nok" << idig + 1 << cellsFixedData[cell].pt << endl;
			ZH2B * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNext4box();
		}
	}

}


//____________ uas collector 

void ZH2GXN::SetupFsol(int * grid0) {
	g0 = grid0;
	memset(fsol, 0, sizeof fsol);
	memset(gangsters, 0, sizeof gangsters);
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		fsol[g0[i]] |= bit;
		gangsters[C_col[i]] |= 1 << g0[i];
	}
}

//________________ 3 digits 
void ZH2_3::GoZ3A(int fl) {
	zh2gxn.nua = 0;
	*zh2gxn.nknownuas = 0;// be sure to start with no ua
	int n = 0;
	uint64_t isfl = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl&bit) {
			isfl |= zh2gxn.fsol[i];
			zh2gxn.maptodigit[n] = i;
			zh2gxn.fsolw[n].bf.u64 = zh2gxn.fsol[i];
			zh2gxn.digit_map[i] = n++;
		}
	}
	cells_unsolved.bf.u64 = isfl;
	zh2gxn.unsolved_field = isfl;
	memset(FD, 0, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}
int ZH2_3::GoZ3(int  fl) {
	if (_popcnt32(fl) != 3) {
		cout << "bug fl not 3 digits" << endl;
		return -1;// not valid fl
	}
	GoZ3A(fl);// start shared with gangsters g2
	// init pm using gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64&bit))continue;
		for (int idig = 0; idig < 3; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 0777777;//3*6 bits
	if(!FullUpdate()) return -1;
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_3::GoZ3G2(int fl, int c1, int d1, int c2, int d2) {
	if (_popcnt32(fl) != 3) {
		cout << "bug fl not 3 digits" << endl;
		return -1;// not valid fl
	}
	GoZ3A(fl);// start shared with gangsters g2

	// init pm using revised gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	int bit12 = (1 << d1)|(1 << d2);
	gx[c1] ^= bit12;// must do remove one add the other
	gx[c2] ^= bit12;

	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64&bit))continue;
		for (int idig = 0; idig < 3; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 0777777;//3*6 bits


	int ir = FullUpdate();
	if (!ir) return -1;
	if (ir == 2) {// immediate ua size 18
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 3; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
		//cout << Char2Xout(w.bf.u64) << " seen " 
		//	<<_popcnt64(w.bf.u64 )<< endl;
		return 1;
	}
	//cout << "after update" << endl;
	//ImageCandidats();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_3::DoZ3Go() {// called in zhou3[1]
	zh2gxn.nua = 0;
	if (*zh2gxn.nknownuas) {
		int ir= GetNetUaCell();
		if(ir<0)return -1; // dead branch
		zh2b_g.guess_xcell =ir; // dead branch
	}
	Guess();
	return 0;
}
inline void ZH2_3::Assign(int rdigit, int cell, int xcell) {
	FD[rdigit] &= AssignMask_Digit[cell].u64[0];
	cells_unsolved.Clear(xcell);
	int ddig = 6 * rdigit;
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
}
int ZH2_3::Seta(int rdigit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell];
	if (FD[rdigit].Off(xcell)) return 1; // not valid
	Assign(rdigit, cell, xcell);
	BF64 *Fd = &FD[rdigit];
	BF64 * RF = &FD[2];
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}
int ZH2_3::FullUpdate() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}
int ZH2_3::ApplySingleOrEmptyCells() {
	zh2b_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map = map[2],
		R3 = R2 & Map;// digits 12
	R2 |= R1 & Map; R1 |= Map;
	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 4; idig++) {
				if (map[idig] & bit) {// this is the digit
					Seta(idig, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		R2 &= ~R3;
		if (!R2) R2 = R3;
		zh2gxn.rx = R2;
		int ir= GetNetUaCell();
		if (ir < 0)return -1; // dead branch
		zh2b_g.guess_xcell=ir;
		return 0;
	}
}
void ZH2_3::Guess() {
	uint32_t xcell = zh2b_g.guess_xcell, cell;
	cell = From_128_To_81[xcell];
	uint64_t bit = (uint64_t)1 << xcell;
	//ImageCandidats();
	int digit = zh2gxn.digit_map[zh2gxn.g0[cell]];
	if (FD[digit].bf.u64 & bit) {
		ZH2_3 * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNext();
	}
	for (int idig = 0; idig < 3; idig++) {
		if (idig == digit) continue;
		if (FD[idig].bf.u64 & bit) {
			ZH2_3 * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNext();
		}
	}
}
int ZH2_3::GetNetUaCell() {
	// setup ok cells
	uint64_t sok = 0; // zh2gxn.unsolved_field;
	for (int i = 0; i < 3; i++) {
		sok |= (FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64);
	}
	uint64_t sok2 = sok & ~cells_unsolved.bf.u64;
	sok &= cells_unsolved.bf.u64;
	int nu = *zh2gxn.nknownuas,cret;
	// first non hit ua to solve
	for (int iu = 0; iu < nu; iu++) {
		register uint64_t U = zh2gxn.knownuas[iu];
		if (U &sok2) continue;
		U &= sok;
		if (!U)return -1; // dead branch
		if (U & zh2gxn.rx) U &= zh2gxn.rx;
		bitscanforward64(cret, U);
		return cret;
	}
	// all uas solved, first cell in rx
	bitscanforward64(cret, zh2gxn.rx);
	return cret;
}
void ZH2_3::ComputeNext() {
	int ir = FullUpdate();
	if (ir == 1)Guess();
	else if (ir == 2) {// solved 
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 3; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		uint64_t cc = _popcnt64(w.bf.u64);
		if (w.bf.u64)	 {
			// check no ua false
			int nu = *zh2gxn.nknownuas;
			if (nu) {
				uint64_t sok = 0; // zh2gxn.unsolved_field;
				for (int i = 0; i < 3; i++) {
					sok |= FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64;
				}
				for (int iu = 0; iu < nu; iu++) {
					register uint64_t U = zh2gxn.knownuas[iu];
					if (!(U &sok)) return;;
				}
			}
			zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
			zh2gxn.knownuas[(*zh2gxn.nknownuas)++]
				= w.bf.u64;
			//cout << Char2Xout(w.bf.u64) << " added" << endl;
		}
	}
}
int ZH2_3::Update() {
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		if (!rows_unsolved.bf.u32[0])break;

		{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
		if (!(AR & 077))goto digit1;

		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL3(0, 2, 46)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL3(1, 3, 5)
		}

	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL3(0, 0, 4)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL3(1, 1, 5)
		}

	digit2:	if (!(AR & 0770000))goto end01234;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL3(0, 0, 2)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto end01234;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL3(1, 1, 3)
		}


	end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


	}// end while

	return 1;
}
void ZH2_3::ImageCandidats() {
	BF64  R2 = FD[0] & FD[1], R1 = FD[0] | FD[1];
	BF64 R3 = R2 & FD[2];
	R2 |= R1 & FD[2];	R1 |= FD[2];

	for (int i = 0; i < 6; i++) { // rows
		if ((i == 3)) {
			for (int ix = 0; ix < 45; ix++)       cout << (char)'-';
			cout << endl;
		}
		for (int j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, xcell = C_To128[cell];
			uint64_t	bit = (uint64_t)1 << xcell;
			if (!(R1.bf.u64&bit)) cout << "-   ";
			else if (R3.bf.u64&bit) cout << "123 ";
			else {
				for (int i = 0; i < 3; i++)
					if (FD[i].bf.u64&bit)cout << i + 1;
				if (R2.bf.u64&bit)cout << "  ";
				else cout << "   ";
			}

		} // end for j
		cout << endl;
	}
	//cout << endl;

}


//________________ 4 digits 
void ZH2_4::GoZ4A(int  fl) {
	zh2gxn.nua = 0;
	int n = 0;
	uint64_t isfl = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl&bit) {
			isfl |= zh2gxn.fsol[i];
			zh2gxn.maptodigit[n] = i;
			zh2gxn.fsolw[n].bf.u64 = zh2gxn.fsol[i];
			zh2gxn.digit_map[i] = n++;
		}
	}
	cells_unsolved.bf.u64 = isfl;
	zh2gxn.unsolved_field = isfl;
	memset(FD, 0, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);

}
int ZH2_4::GoZ4(int  fl) {
	if (_popcnt32(fl) != 4) {
		cout << "bug fl not 4 digits" << endl;
		return -1;// not valid fl
	}
	GoZ4A( fl);
	// init pm using gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64&bit))continue;
		for (int idig = 0; idig < 4; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 077777777;//4*6 bits
	FullUpdate();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;;
}
int ZH2_4::GoZ4G2(int fl, int c1, int d1, int c2, int d2) {
	if (_popcnt32(fl) != 4) {
		cout << "bug fl not 4 digits" << endl;
		return -1;// not valid fl
	}
	GoZ4A(fl);// start shared with gangsters g2

	// init pm using revised gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	int bit12 = (1 << d1) | (1 << d2);
	gx[c1] ^= bit12;// must do remove one add the other
	gx[c2] ^= bit12;

	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64&bit))continue;
		for (int idig = 0; idig < 4; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 077777777;//4*6 bits


	int ir = FullUpdate();
	if (!ir) return -1;
	if (ir == 2) {// immediate ua size 18
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 4; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
		//cout << Char2Xout(w.bf.u64) << " seen " 
		//	<<_popcnt64(w.bf.u64 )<< endl;
		return 1;
	}
	//cout << "after update" << endl;
	//ImageCandidats();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_4::DoZ4Go() {// called in zhou3[1]
	zh2gxn.nua = 0;
	if (*zh2gxn.nknownuas) {
		int ir = GetNetUaCell();
		if (ir < 0)return -1; // dead branch
		zh2b_g.guess_xcell = ir; // dead branch
	}
	Guess();
	return 0;
}
inline void ZH2_4::Assign(int rdigit, int cell, int xcell) {
	FD[rdigit] &= AssignMask_Digit[cell].u64[0];
	cells_unsolved.Clear(xcell);
	int ddig = 6 * rdigit;
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
}
int ZH2_4::Seta(int rdigit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell],
		block = TblBoard_Block[cell];
	if (FD[rdigit].Off(xcell)) return 1; // not valid
	Assign(rdigit, cell, xcell);
	BF64 *Fd = &FD[rdigit];
	BF64 * RF = &FD[3];
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}
int ZH2_4::FullUpdate() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}
int ZH2_4::ApplySingleOrEmptyCells() {
	zh2b_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map = map[2],
		R3 = R2 & Map;// digits 12
	R2 |= R1 & Map; R1 |= Map;
	Map = map[3]; 	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 4; idig++) {
				if (map[idig] & bit) {// this is the digit
					Seta(idig, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		R2 &= ~R3;
		if (!R2) R2 = R3;
		zh2gxn.rx = R2;
		int ir = GetNetUaCell();
		if (ir < 0)return -1; // dead branch
		zh2b_g.guess_xcell = ir;
		return 0;
	}
}
void ZH2_4::Guess() {
	uint32_t xcell = zh2b_g.guess_xcell, cell;
	cell = From_128_To_81[xcell];
	uint64_t bit = (uint64_t)1 << xcell;
	//cout << "guess" << cellsFixedData[cell].pt << endl;
	int digit = zh2gxn.digit_map[zh2gxn.g0[cell]];
	if (FD[digit].bf.u64 & bit) {
		//cout << "guess ok" <<digit+1<< cellsFixedData[cell].pt << endl;
		ZH2_4 * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNext();
	}
	for (int idig = 0; idig < 4; idig++) {
		if (FD[idig].bf.u64 & bit) {
			if (idig == digit) continue;
			//cout << "guess nok" << idig + 1 << cellsFixedData[cell].pt << endl;
			ZH2_4 * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNext();
		}
	}
}
int ZH2_4::GetNetUaCell() {
	// setup ok cells
	uint64_t sok = 0; // zh2gxn.unsolved_field;
	for (int i = 0; i < 4; i++) {
		sok |= (FD[i] & zh2gxn.fsolw[i]).bf.u64;
	}
	uint64_t sok2= sok & ~cells_unsolved.bf.u64;
	sok&= cells_unsolved.bf.u64;
	//cout << Char2Xout(sok) << " sok next cell" << endl;
	//cout << Char2Xout(sok2) << " sok2 next cell" << endl;
	int nu = *zh2gxn.nknownuas, cret;
	// first non hit ua to solve
	for (int iu = 0; iu < nu; iu++) {
		register uint64_t U = zh2gxn.knownuas[iu];
		if (U &sok2) continue;
		//cout << Char2Xout(U) << "U to fill" << endl;
		U &= sok;
		if (!U)return -1; // dead branch

		// must be on possible true if not dead
		//this ua would be a subset
		if (U & zh2gxn.rx) U &= zh2gxn.rx;
		bitscanforward64(cret, U);
		return cret;
	}
	// all uas solved, first cell in rx
	bitscanforward64(cret, zh2gxn.rx);
	return cret;
}
void ZH2_4::ComputeNext() {
	int ir = FullUpdate();
	//cout << "compnext ir=" << ir << endl;
	//ImageCandidats();
	if (ir == 1)Guess();
	else if (ir == 2) {// solved 
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 4; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		uint64_t cc = _popcnt64(w.bf.u64);
		if (w.bf.u64)	 {
			// check no ua false
			//cout << Char2Xout(w.bf.u64) << " seen" << endl;
			int nu = *zh2gxn.nknownuas;
			if (nu) {
				uint64_t sok = 0; // zh2gxn.unsolved_field;
				for (int i = 0; i < 4; i++) {
					sok |= FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64;
				}
				for (int iu = 0; iu < nu; iu++) {
					register uint64_t U = zh2gxn.knownuas[iu];
					if (!(U &sok)) return;;
				}
			}
			zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
			zh2gxn.knownuas[(*zh2gxn.nknownuas)++]
				= w.bf.u64;
			//cout << Char2Xout(w.bf.u64) << " added zh2gxn.nua="<< zh2gxn.nua << endl;
		}
	}
}
int ZH2_4::Update() {
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		if (!rows_unsolved.bf.u32[0])break;

		{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
		if (!(AR & 077))goto digit1;

		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL4(0, 2, 4, 6)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL4(1, 3, 5, 7)
		}

	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL4(0, 0, 4, 6)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL4(1, 1, 5, 7)
		}

	digit2:	if (!(AR & 0770000))goto digit3;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL4(0, 0, 2, 6)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto digit3;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL4(1, 1, 3, 7)
		}

	digit3: if (!(AR & 077000000))goto end01234;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		UPDN(3, 0)	  if (((AR >> 18) & 7) != S) {
			AR &= 07770777777 | (S << 18); UPWCL4(0, 0, 2, 4)
		}

	digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto end01234;
		UPDN(3, 1)if (((AR >> 21) & 7) != S) {
			AR &= 07707777777 | (S << 21); UPWCL4(1, 1, 3, 5)
		}


	end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


	}// end while

	return 1;
}
void ZH2_4::ImageCandidats() {
	BF64  R2 = FD[0] & FD[1], R1 = FD[0] | FD[1];
	BF64 R3 = R2 & FD[2];
	R2 |= R1 & FD[2];	R1 |= FD[2];
	BF64 R4 = R3 & FD[3];
	R3 |= R2 & FD[3]; R2 |= R1 & FD[3];	R1 |= FD[3];

	for (int i = 0; i < 6; i++) { // rows
		if ((i == 3)) {
			for (int ix = 0; ix < 45; ix++)       cout << (char)'-';
			cout << endl;
		}
		for (int j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, xcell = C_To128[cell];
			uint64_t	bit = (uint64_t)1 << xcell;
			if (!(R1.bf.u64&bit)) cout << "-    ";
			else if (R4.bf.u64&bit) cout << "1234 ";
			else {
				for (int i = 0; i < 4; i++)
					if (FD[i].bf.u64&bit)cout << i + 1;
				if (R3.bf.u64&bit)cout << "  ";
				else if (R2.bf.u64&bit)cout << "   ";
				else cout << "    ";
			}

		} // end for j
		cout << endl;
	}
	cout << endl;

}

//_______________________________  5 digits
void ZH2_5::GoZ5A(int  fl) {
	zh2gxn.nua = 0;
	int n = 0;
	uint64_t isfl = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl&bit) {
			isfl |= zh2gxn.fsol[i];
			zh2gxn.maptodigit[n] = i;
			zh2gxn.fsolw[n].bf.u64 = zh2gxn.fsol[i];
			zh2gxn.digit_map[i] = n++;
		}
	}
	cells_unsolved.bf.u64 = isfl;
	zh2gxn.unsolved_field = isfl;
	memset(FD, 0, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}
	int ZH2_5::GoZ5(int  fl) {
		if (_popcnt32(fl) != 5) {
			cout << "bug fl not 5 digits" << endl;
			return -1;// not valid fl
		}
		GoZ5A(fl);
		uint32_t gx[9];
		for (int i = 0; i < 9; i++) {
			gx[i] = zh2gxn.gangsters[i] & fl;
		}
		for (int i = 0; i < 54; i++) {
			int xi = C_To128[i];
			uint64_t bit = (uint64_t)1 << xi;
			if (!(cells_unsolved.bf.u64&bit))continue;
			for (int idig = 0; idig < 5; idig++) {
				uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
				if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
			}
		}
		rows_unsolved.bf.u64 = 07777777777;//5*6 bits
		//FullUpdate();
		uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
		return (int)cc;;
	}
	int ZH2_5::GoZ5G2(int fl, int c1, int d1, int c2, int d2) {
		if (_popcnt32(fl) != 5) {
			cout << "bug fl not 5 digits" << endl;
			return -1;// not valid fl
		}
		GoZ5A(fl);
		uint32_t gx[9];
		for (int i = 0; i < 9; i++) {
			gx[i] = zh2gxn.gangsters[i] & fl;
		}
		int bit12 = (1 << d1) | (1 << d2);
		gx[c1] ^= bit12;// must do remove one add the other
		gx[c2] ^= bit12;
		for (int i = 0; i < 54; i++) {
			int xi = C_To128[i];
			uint64_t bit = (uint64_t)1 << xi;
			if (!(cells_unsolved.bf.u64&bit))continue;
			for (int idig = 0; idig < 5; idig++) {
				uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
				if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
			}
		}
		rows_unsolved.bf.u64 = 07777777777;//5*6 bits


		int ir = FullUpdate();
		if (!ir) return -1;
		if (ir == 2) {// immediate ua  
			BF64 w; w.bf.u64 = 0;
			for (int i = 0; i < 5; i++) {
				w |= FD[i] - zh2gxn.fsolw[i];
			}
			zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
			//cout << Char2Xout(w.bf.u64) << " seen " 
			//	<<_popcnt64(w.bf.u64 )<< endl;
			return 1;
		}
		//cout << "after update" << endl;
		//ImageCandidats();
		uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
		return (int)cc;;
	}



	int ZH2_5::DoZ5Go() {// called in zhou3[1]
		zh2gxn.nua = 0;
		if (*zh2gxn.nknownuas) {
			int ir = GetNetUaCell();
			if (ir < 0)return -1; // dead branch
			zh2b_g.guess_xcell = ir; // dead branch
		}
		Guess();
		return 0;
	}
	inline void ZH2_5::Assign(int rdigit, int cell, int xcell) {
		FD[rdigit] &= AssignMask_Digit[cell].u64[0];
		cells_unsolved.Clear(xcell);
		int ddig = 6 * rdigit;
		rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	}
	int ZH2_5::Seta(int rdigit, int xcell) { // single in cell
		int cell = From_128_To_81[xcell],
			block = TblBoard_Block[cell];
		if (FD[rdigit].Off(xcell)) return 1; // not valid
		Assign(rdigit, cell, xcell);
		BF64 *Fd = &FD[rdigit];
		BF64 * RF = &FD[4];
		for (; RF >= FD; RF--)RF->Clear(xcell);
		Fd->Set(xcell); // restore bit for digit assigned
		return 0;
	}
	int ZH2_5::FullUpdate() {
		while (1) {
			if (!Update()) return 0; // game locked in update
			if (!Unsolved_Count()) return 2;
			if (ApplySingleOrEmptyCells())	return 0; // locked 
			if (zh2b_g.single_applied) 			continue;
			break;
		}
		return 1;
	}
	int ZH2_5::ApplySingleOrEmptyCells() {
		zh2b_g.single_applied = 0;
		uint64_t * map = &FD[0].bf.u64;
		uint64_t unsolved = cells_unsolved.bf.u64;
		register uint64_t R2 = map[0] & map[1],
			R1 = (map[0] | map[1]), Map = map[2],
			R3 = R2 & Map, R4;// digits 12
		R2 |= R1 & Map; R1 |= Map;

		Map = map[3]; R4 = R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
		Map = map[4];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

		if (unsolved & (~R1)) return 1; // locked
		R1 &= ~R2;
		R1 &= unsolved; // these are new singles	
		if (R1) {
			zh2b_g.single_applied = 1;
			while (R1) {// usually a very small number of cells to assign
				uint32_t res;
				if (!bitscanforward64(res, R1)) break;
				uint64_t bit = (uint64_t)1 << res; // switch to the bit value
				R1 &= ~bit;  // clear the bit
				for (int idig = 0; idig < 9; idig++) {
					if (map[idig] & bit) {// this is the digit
						int cell = From_128_To_81[res];
						Assign(idig, cell, res);
						goto nextr1;// finished for that cell
					}
				}
				return 1; //conflict with a previous cell assugn
			nextr1: {}
			}
			return 0;
		}
		else {
			R2 &= ~R3;
			R3 &= ~R4;
			if (!R2) R2 = R3;
			if (!R2) R2 = R4;
			zh2gxn.rx = R2;
			int ir = GetNetUaCell();
			if (ir < 0)return -1; // dead branch
			zh2b_g.guess_xcell = ir;
			return 0;
		}
	}

	int ZH2_5::GetNetUaCell() {
		// setup ok cells
		uint64_t sok = 0; // zh2gxn.unsolved_field;
		for (int i = 0; i < 5; i++) {
			sok |= (FD[i] & zh2gxn.fsolw[i]).bf.u64;
		}
		uint64_t sok2 = sok & ~cells_unsolved.bf.u64;
		sok &= cells_unsolved.bf.u64;
		int nu = *zh2gxn.nknownuas, cret;
		//cout << Char2Xout(sok) << " sok next cell" << endl;
		//cout << Char2Xout(sok2) << " sok2 next cell" << endl;
			// first non hit ua to solve
		for (int iu = 0; iu < nu; iu++) {
			register uint64_t U = zh2gxn.knownuas[iu];
			if (U &sok2) continue;
			//cout << Char2Xout(U) << "U to fill" << endl;
			U &= sok;
			if (!U)return -1; // dead branch
			if (U & zh2gxn.rx) U &= zh2gxn.rx;
			bitscanforward64(cret, U);
			return cret;
		}
		// all uas solved, first cell in rx
		bitscanforward64(cret, zh2gxn.rx);
		return cret;
	}
	void ZH2_5::Guess() {// 
		uint32_t xcell = zh2b_g.guess_xcell, cell;
		cell = From_128_To_81[xcell];
		uint64_t bit = (uint64_t)1 << xcell;
		//cout << "guess" << cellsFixedData[cell].pt << endl;
		int digit = zh2gxn.digit_map[zh2gxn.g0[cell]];
		if (FD[digit].bf.u64 & bit) {
			//cout << "guess ok" <<digit+1<< cellsFixedData[cell].pt << endl;
			ZH2_5 * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(digit, xcell);
			mynext->ComputeNext();
		}
		for (int idig = 0; idig < 5; idig++) {
			if (FD[idig].bf.u64 & bit) {
				if (idig == digit) continue;
				//cout << "guess nok" << idig + 1 << cellsFixedData[cell].pt << endl;
				ZH2_5 * mynext = this + 1; // start next guess
				*mynext = *this;
				mynext->Seta(idig, xcell);
				mynext->ComputeNext();
			}
		}
	}
	void ZH2_5::ComputeNext() {
		//ImageCandidats();
		int ir = FullUpdate();
		//cout << "compnext ir=" << ir << endl;
		//ImageCandidats();
		if (ir == 1)Guess();
		else if (ir == 2) {// solved 
			BF64 w; w.bf.u64 = 0;
			for (int i = 0; i < 5; i++) {
				w |= FD[i] - zh2gxn.fsolw[i];
			}
			uint64_t cc = _popcnt64(w.bf.u64);
			if (w.bf.u64) {
				// check no ua false
				//cout << Char2Xout(w.bf.u64) << " seen" << endl;
				int nu = *zh2gxn.nknownuas;
				if (nu) {
					uint64_t sok = 0; // zh2gxn.unsolved_field;
					for (int i = 0; i < 5; i++) {
						sok |= FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64;
					}
					for (int iu = 0; iu < nu; iu++) {
						register uint64_t U = zh2gxn.knownuas[iu];
						if (!(U &sok)) return;;
					}
				}
				zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
				zh2gxn.knownuas[(*zh2gxn.nknownuas)++]
					= w.bf.u64;
				//cout << Char2Xout(w.bf.u64) << " added zh2gxn.nua="<< zh2gxn.nua
				//	<< " " << cc<< endl;
			}
		}
	}

	int ZH2_5::Update() {
		int Shrink = 1;
		register int S, A;
		register unsigned int cl, *wcl = FD[0].bf.u32;
		while (Shrink) {
			Shrink = 0;
			if (!rows_unsolved.bf.u32[0])break;

			{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
			if (!(AR & 077))goto digit1;

			//=digit 0
			if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
			UPDN(0, 0)	if ((AR & 7) != S) {
				AR &= 07777777770 | S;	UPWCL5(0, 2, 4, 6, 8)
			}

		digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
			UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
				AR &= 07777777707 | (S << 3);	UPWCL5(1, 3, 5, 7, 9)
			}

		digit1:	if (!(AR & 07700))goto digit2;

			if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
			UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
				AR &= 07777777077 | (S << 6); UPWCL5(0, 0, 4, 6, 8)
			}

		digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
			UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
				AR &= 07777770777 | (S << 9); UPWCL5(1, 1, 5, 7, 9)
			}

		digit2:	if (!(AR & 0770000))goto digit3;

			if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
			UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
				AR &= 07777707777 | (S << 12);	UPWCL5(0, 0, 2, 6, 8)
			}

		digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto digit3;
			UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
				AR &= 07777077777 | (S << 15);	UPWCL5(1, 1, 3, 7, 9)
			}

		digit3: if (!(AR & 077000000))goto digit4;

			if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
			UPDN(3, 0)	  if (((AR >> 18) & 7) != S) {
				AR &= 07770777777 | (S << 18); UPWCL5(0, 0, 2, 4, 8)
			}

		digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
			UPDN(3, 1)if (((AR >> 21) & 7) != S) {
				AR &= 07707777777 | (S << 21); UPWCL5(1, 1, 3, 5, 9)
			}

		digit4:if (!(AR & 07700000000))goto end01234;

			if (FD[4].bf.u32[0] == CompFD[4].bf.u32[0])goto digit4b;
			UPDN(4, 0)if (((AR >> 24) & 7) != S) {
				AR &= 07077777777 | (S << 24);  UPWCL5(0, 0, 2, 4, 6)
			}

		digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto end01234;
			UPDN(4, 1)if (((AR >> 27) & 7) != S) {
				AR &= 0777777777 | (S << 27);  UPWCL5(1, 1, 3, 5, 7)
			}

		end01234: rows_unsolved.bf.u32[0] = AR;
			}// end of validity for AR 01234

		}// end while

		return 1;
	}

	void ZH2_5::ImageCandidats() {
		BF64  R2 = FD[0] & FD[1], R1 = FD[0] | FD[1];
		BF64 R3 = R2 & FD[2];
		R2 |= R1 & FD[2];	R1 |= FD[2];
		BF64 R4 = R3 & FD[3];
		R3 |= R2 & FD[3]; R2 |= R1 & FD[3];	R1 |= FD[3];
		BF64 R5 = R4 & FD[4];
		R4 |= R3 & FD[4]; R3 |= R2 & FD[4];
		R2 |= R1 & FD[4];	R1 |= FD[4];

		for (int i = 0; i < 6; i++) { // rows
			if ((i == 3)) {
				for (int ix = 0; ix < 54; ix++)       cout << (char)'-';
				cout << endl;
			}
			for (int j = 0; j < 9; j++) {
				if ((j == 3) || (j == 6))cout << "|";
				int cell = 9 * i + j, xcell = C_To128[cell];
				uint64_t	bit = (uint64_t)1 << xcell;
				if (!(R1.bf.u64&bit)) cout <<   "-     ";
				else if (R5.bf.u64&bit) cout << "12345 ";
				else {
					for (int i = 0; i < 5; i++)
						if (FD[i].bf.u64&bit)cout << i + 1;
					if (R4.bf.u64&bit)cout << "  ";
					else if (R3.bf.u64&bit)cout << "   ";
					else if (R2.bf.u64&bit)cout << "    ";
					else                   cout << "     ";
				}

			} // end for j
			cout << endl;
		}
		cout << endl;

	}


//============================= ZH2B5 uas guas creation

/*
process using as source
	fd_sols[1] is pm if source=0  fd_revised if  source=1
the solution is known and is fd_sols[0]
all digits not in bitfield fl are known (usualy fl 2-4 digits)
internal myfd as initial value  cleaned in zb2b_1d step by step

*/
uint64_t  ZH2B5_GLOBAL::FindUAsInit(int fl, int source) {
	BF64 * mypm = zh2b_g.fd_sols[1];
	if (source) mypm = zh2b_g.fd_revised;
	uint64_t solved_cells = 0;
	uint32_t nd = 0;
	for (int idig = 0, bit = 1; idig < 9; idig++, bit <<= 1) {
		if (fl & bit) {
			fdsw[0][nd] = zh2b_g.fd_sols[0][idig];
			fdsw[2][nd].bf.u64 = (~fdsw[0][nd].bf.u64) & BIT_SET_2X;
			myfd[nd++] = mypm[idig];
		}
		else solved_cells |= zh2b_g.fd_sols[0][idig].bf.u64;
	}
	ndigits = nd;
	cells_unsolved.bf.u64 = BIT_SET_2X ^ solved_cells;
	for (uint32_t i = 0; i < nd; i++) {
		myfd[i] &= cells_unsolved;
		if (myfd[i].Count() == 6)return 0;
	}
	return solved_cells;
	// first cleaning return unsolved cells in bits 	
}
void ZH2B5_GLOBAL::CollectUas5() {
	if (diag)cout << "entry CollectUas5()" << endl;
	nuaf5 = 0;
	memset(&zh2b5[0], 0, sizeof zh2b5[0]);
	memcpy(zh2b5[0].FD, myfd, sizeof myfd);
	memset(zh2b5[0].CompFD, 0, sizeof myfd);
	zh2b5[0].cells_unsolved = cells_unsolved;
	uint32_t t[5] = { 077,07777,0777777,077777777,07777777777 };
	zh2b5[0].rows_unsolved.bf.u32[0] = t[ndigits-1];
	if(diag)zh2b5[0].ImageCandidats();
	//if (1) return;
	zh2b5[0].ComputeNext5();
}
void ZH2B5_GLOBAL::ValidSol5(uint64_t * sol) {//Ua to catch
	if (zh2b5_g.diag) {
		cout << "valid sol nuaf5="<<nuaf5<< endl;
	}
	uint64_t ua = 0, *s0 = &fdsw[0][0].bf.u64;
	for (uint32_t i = 0; i < ndigits; i++) {
		uint64_t w = sol[i] & ~s0[i];// digit ua
		if (!w) {
			if (zh2b5_g.diag) cout << "subset assumed rdigit" << i << endl;
			return; //a subset exists
		}
		ua |= w;
	}
	/* not true can be small gua not a band ua
	if (!modevalid) {// not gua mode must be 1n both bands
		register uint64_t R = BIT_SET_27;
		if (!(ua&R))return;; // ua in band2
		R <<= 32;
		if (!(ua&R))return; // ua in band1
	}
	*/
	uint64_t  cc = (int)_popcnt64(ua);
	if (zh2b5_g.diag)cout << "cc=" << cc << endl;
	if (cc > sizef5) return;
	if (nuaf5 >= 40 && modevalid) {// reasonnable limit for a given step
		int limit = tuaf5[39].bf.u64 >> 59;
		if (cc >= limit)return; // too many uas found here, skip it
		nuaf5 = 39;
	}
	ua |= cc << 59;
	//int ir = genuasb12.AddUA64(&tuaf5[0].bf.u64, nuaf5);
	//uint64_t cc = _popcnt64(genuasb12.ua);
	//genuasb12.ua.bf.u64  |= cc << 59;
	AddUA64(&tuaf5[0].bf.u64, nuaf5,ua);
	if (zh2b5_g.diag)cout << Char2Xout(ua) << " ua found end nuaf5=" << nuaf5 << "cc=" << cc << endl;
}

//======================== ZH2B5  2-5 digits
int ZH2B5::Update5() {
	if(zh2b5_g.ndigits>5)return 0; // force false if bad use
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		if (!rows_unsolved.bf.u32[0])return 1;// solved

		{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
		if (!(AR & 077))goto digit1;

		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL5(0, 2, 4, 6, 8)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL5(1, 3, 5, 7, 9)
		}

	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL5(0, 0, 4, 6, 8)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL5(1, 1, 5, 7, 9)
		}

	digit2:	
		if (zh2b5_g.ndigits < 3) goto end01234;
		if (!(AR & 0770000))goto digit3;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL5(0, 0, 2, 6, 8)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto digit3;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL5(1, 1, 3, 7, 9)
		}

	digit3: 
		if (zh2b5_g.ndigits < 4) goto end01234;
		if (!(AR & 077000000))goto digit4;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		UPDN(3, 0)	  if (((AR >> 18) & 7) != S) {
			AR &= 07770777777 | (S << 18); UPWCL5(0, 0, 2, 4, 8)
		}

	digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
		UPDN(3, 1)if (((AR >> 21) & 7) != S) {
			AR &= 07707777777 | (S << 21); UPWCL5(1, 1, 3, 5, 9)
		}

	digit4:if (zh2b5_g.ndigits < 5) goto end01234;
		if (!(AR & 07700000000))goto end01234;

		if (FD[4].bf.u32[0] == CompFD[4].bf.u32[0])goto digit4b;
		UPDN(4, 0)if (((AR >> 24) & 7) != S) {
			AR &= 07077777777 | (S << 24);  UPWCL5(0, 0, 2, 4, 6)
		}

	digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto end01234;
		UPDN(4, 1)if (((AR >> 27) & 7) != S) {
			AR &= 0777777777 | (S << 27);  UPWCL5(1, 1, 3, 5, 7)
		}

	end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234

	}// end while
	return 1;
}
int ZH2B5::FullUpdate5() {
	while (1) {
		if (!Update5()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (zh2b5_g.ndigits > 3) {
			if (ApplySingleOrEmptyCells5())	return 0; // locked 
			if (zh2b5_g.single_applied) continue;
		}
		break;
	}
	return 1;
}
int ZH2B5::ApplySingleOrEmptyCells5() {
#define NAKED5(X) Map=map[X];R2|=R1&Map;R1|=Map 
	zh2b5_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map;// digits 12
	if (zh2b5_g.ndigits > 2)NAKED5(2);
	if (zh2b5_g.ndigits > 3)NAKED5(3);
	if (zh2b5_g.ndigits > 4)NAKED5(4);
	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // forget solved seen as singles
	if (!R1) return 0;
	zh2b5_g.single_applied = 1;

	while (R1) {// usually a very small number of cells to assign
		uint32_t  res;
		if (!bitscanforward64(res, R1)) break;
		uint64_t bit = (uint64_t)1 << res; // switch to the bit value
		R1 &= ~bit;  // clear the bit
		// call Seta(int digit, int xcell) so find digit
		for (uint32_t idig = 0; idig < zh2b5_g.ndigits; idig++) {
			if (map[idig] & bit) {// this is the digit
				if (FD[idig].Off(res))  return 1; // invalid, gane locked
				Seta5(idig, res);
				goto nextr1;// finished for that cell
			}
		}
		return 1; //conflict with a previous cell assugn
	nextr1: {}
	}
	return 0;// not locked
}
int ZH2B5::Seta5(int digit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell],
		block = TblBoard_Block[cell];
	if (FD[digit].Off(xcell)) return 1; // not valid
	//Assign(digit, cell, xcell);
	BF64 *Fd = &FD[digit];
	*Fd &= AssignMask_Digit[cell].u64[0];
	int ddig = 6 * digit;
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	cells_unsolved.Clear(xcell);
	BF64 * RF = &FD[zh2b5_g.ndigits - 1];//last used digit
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}
void  ZH2B5::Guess5() {// solve in table digit with lowest count
	if (zh2b5_g.diag) {
		cout << oct << rows_unsolved.bf.u64 << dec << " unsolved guess5" << endl;
		ImageCandidats();
	}
	if (!rows_unsolved.bf.u64) {// this is a solution
		//cout << "valid sol" << endl;
		zh2b5_g.ValidSol5(&FD[0].bf.u64);
		return;
	}
	int mincount = 100, digmin,ncount=0;
	for (uint32_t idig = 0; idig < zh2b5_g.ndigits; idig++) {
		int cc = FD[idig].Count();
		if (cc < 7) continue;
		ncount++;
		if (cc < mincount) {
			mincount = cc;
			digmin = idig;
		}
	}
	if (zh2b5_g.diag) {
		cout << oct << rows_unsolved.bf.u64 << dec << " unsolved guess5 for dig"
			<<digmin+1<< endl;
		//ImageCandidats();
	}
	// put in table all digmin valid solutions 
	BF64 tuaw[100], tsolw[100];
	//cout << "call solve 1 digit for digitw=" << digmin+1 << endl;
	int nuaw=zh2b1d_g.Go(zh2b5_g.fdsw[0][digmin], FD[digmin], tsolw, tuaw,
		(rows_unsolved.bf.u32[0]) >> (6 * digmin) & 077);
	rows_unsolved.bf.u32[0] &= ~(077 << (6 * digmin));// clear unsolved
	if (zh2b5_g.diag)cout << "return nuaw=" << nuaw << " digmin="<<digmin+1<< endl;
	//try each digit full solution if not fully true (subset would exist)
	for (int i = 0; i < nuaw; i++) {
		ZH2B5 * nextz = this + 1;// open a new zh2b
		*nextz = *this;
		{// clear the solution for other digits
			register uint64_t R = tsolw[i].bf.u64, nR = ~R;
			if (zh2b5_g.diag)cout<<Char2Xout(R) << "apply for i=" << i << " digmin=" << digmin+1 << endl;
			uint64_t * RF = &nextz->FD[zh2b5_g.ndigits - 1].bf.u64;
			for (; RF >= &nextz->FD->bf.u64; RF--)*RF &= nR;
			nextz->FD[digmin].bf.u64 = R;// restore the solution for digit
		}	
		// if ncount=2 it is a solution
		if (ncount == 2) {
			if (zh2b5_g.diag) nextz->ImageCandidats();
			zh2b5_g.ValidSol5(&nextz->FD[0].bf.u64);
		}
		else nextz->ComputeNext5();// continue solving
	}
}
int ZH2B5::GetAllDigits(int cell) {
	int ir = 0, xcell = C_To128[cell];;
	for (uint32_t i = 0; i < zh2b5_g.ndigits; i++) if (FD[i].On(xcell))ir |= (1 << i);
	return ir;
}
void ZH2B5::ImageCandidats() {
	int dig_cells[81]; for (int i = 0; i < 54; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " << endl << endl;
	for (i = 0; i < 9; i++) {  // attention ici i indice colonne
		lcol[i] = 2;    // 2  mini tous chiffres impos�s
		for (j = 0; j < 6; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 6; i++) { // maintenant indice ligne
		if ((i == 3) || (i == 6)) {
			for (int ix = 0; ix < (tcol + 10); ix++)       cout << (char)'-';
			cout << endl;
		}
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell],
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < 9; id++)if (digs & (1 << id))
				cout << id + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}


//======================== ZH2B1B one digit 2 bands table of solutions

int ZH2B_1D_GLOBAL::Go(BF64 & sol, BF64 & fde, BF64 *tsol, BF64 *tua,int ru) {
	tsolw = tsol;
	tuaw = tua;
	mysol = sol;
	nsolw = 0;
	myandsol.bf.u64 = BIT_SET_2X;
	zh2b1d[0].FD = fde;
	zh2b1d[0].CompFD.bf.u64 = 0;
	zh2b1d[0].ComputeNext(ru);
	return nsolw;
}
int ZH2B_1D::GetSols( int ru) {
	CompFD.bf.u64 = 0;
	FD = zh2b_g.mystart &zh2b_g.myandsol;
	zh2b_g.nsolw = 0;
	ComputeNext(ru);
	return zh2b_g.nsolw;
}
int ZH2B_1D::GetAllSols(BF64 & fde, int ru, BF64 & fdsol) {
	zh2b_g.mystart = FD = fde;
	zh2b_g.mysol = fdsol;
	zh2b_g.myandsol.bf.u64 = BIT_SET_2X;
	zh2b_g.nsolw = 0;
	ComputeNext(ru);
	return zh2b_g.nsolw;
}
void ZH2B_1D::ComputeNext(int ru) {
	if (Update(ru)) {
		if (!ru) {// new sol put it in table
			if (FD != zh2b1d_g.mysol) {
				zh2b1d_g.tuaw[zh2b1d_g.nsolw] = FD- zh2b1d_g.mysol;// partial ua
				zh2b1d_g.tsolw[zh2b1d_g.nsolw++] = FD;// partial solution
			}
		}
		else Guess(ru);
	}
}
int ZH2B_1D::Update(int &ru) {// solve as much as possible 
	int Shrink = 1;
	register int S, A;
	while (Shrink) {
		Shrink = 0;
		if (!(ru & 7)) goto band2;
		A = FD.bf.u32[0];
		if (A == CompFD.bf.u32[0]) goto band2;
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[1] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		CompFD.bf.u32[0] = 
			FD.bf.u32[0] = A;
		ru = (ru & 070) | S;
	band2:;
		if (!(ru & 070)) goto exit;
		A = FD.bf.u32[1];
		if (A == CompFD.bf.u32[1]) goto exit;
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[0] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		CompFD.bf.u32[1] = 
			FD.bf.u32[1] = A;
		ru = (ru & 7) | (S<<3);
	exit:;
	}
	//cout << Char2Xout(FD.bf.u64) << " sortie update ru=" 
	//	<< oct << ru << dec << endl;

	return 1;
}
void ZH2B_1D::Guess(int ru) {//
	int dcell = 0,v,ruw=ru;
	{
		register uint32_t a = FD.bf.u32[0], b = FD.bf.u32[1];
		if (!(ru & 7)) goto guess_b2;
		if ((ru & 070) && (_popcnt32(b) < _popcnt32(a))) goto guess_b2;
		// guess in band1
		ruw = (-ru) &ru;
		v = FD.bf.u32[0] & TblRowUnsolved[ruw];// unknown in last unknown row
		goto gogo;
	guess_b2:	// guess in band2
		dcell = 27;
		ruw=ru >> 3;
		ruw = (-ruw) &ruw;
		v = FD.bf.u32[1] & TblRowUnsolved[ruw];// unknown in last unknown row
		ruw <<= 3; //relocate bit to the right place
	}
gogo:;
	uint32_t cell;
	ru ^= ruw;// clear the row for next steps
	while (bitscanforward(cell, v)) {
		v ^= 1 << cell;
		cell += dcell;
		ZH2B_1D * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(cell);
		mynext->ComputeNext(ru);
	}
}
int ZH2B_1D::IsValid(uint64_t v) {
	CompFD.bf.u64 = 0;
	FD.bf.u64 = v;
	int ru = 077;
	return Update(ru);
}

//=================================== ZHONE


#define UPDN1(I,P,Q,R,T,U,V,W,X)Shrink = (TblShrinkMask[A & 0x1FF] | \
TblShrinkMask[ (A>>9) & 0x1FF]<<3 | \
TblShrinkMask[ (A>>18) & 0x1FF]<<6);\
if ((A &=TblComplexMask[Shrink]) ==0)  return 0; \
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF); \
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]]; \
if ((A>>27) != S){\
cl = ~(A & TblRowMask[S]); \
A=(A&BIT_SET_27) | (S<<27);\
cells_unsolved &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;\
wcl[U]&= cl;wcl[V]&= cl;wcl[W]&= cl;wcl[X]&= cl;}\
CompFD[I] = FD[I] = A
//===================== now ZHone code 

ZHONE_GLOBAL::ZHONE_GLOBAL() {
	zsol  = 0; // no solution unless required buy the user
}
void ZHONE_GLOBAL::SetUp(int * b, uint32_t * t, uint32_t nt) {// source
	band0 = b; tua = t;  nua = nt;
	// build sol per digit and pm per digit at start
	memset(fds, 0, sizeof fds);
	memset(pms, 0, sizeof pms);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = b[cell];
			pms[dig] |= Zhoucol << i; // add candidates in the column
			fds[dig] |= 1 << cell;
		}
	}

}

//============================================ ZHONE
void ZHONE::Init() {//init after zh1b_g getband
	cells_unsolved = BIT_SET_27;
	memcpy(FD, zh1b_g.pms, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	for (int i = 0; i < 9; i++)	FD[i] |= 7 << 27;//set unknown rows
	zh1b_g.ndigits = 9;
}



#define NAKED(X) 	Map=map[X];R3|=R2&Map;R2|=R1&Map;R1|=Map;
int ZHONE::ApplySingleOrEmptyCells() {
	zh1b_g.single_applied = 0;
	uint32_t * map = FD, unsolved = cells_unsolved;
	register int R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map, R3;// digits 12
	Map = map[2]; R3 = R2 & Map; R2 |= R1 & Map; R1 |= Map;
	NAKED(3) NAKED(4) NAKED(5) NAKED(6)	NAKED(7) NAKED(8) // digits 3-9
		if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // forget solved seen as singles
	if (R1) zh1b_g.single_applied = 1;
	else {
		zh1b_g.pairs = R2 & (~R3);
		zh1b_g.triplets = R3;
		return 0;
	}
	while (R1) {// usually a very small number of cells to assign
		uint32_t res;
		if (!bitscanforward(res, R1)) break;
		int cell = res, bit = 1 << cell; // switch to the bit value
		//		cout << "naked cell" << From_128_To_81[res]<< endl;
		R1 &= ~bit;  // clear the bit
		// call Seta(int digit, int xcell) so find digit 
		for (int idig = 0; idig < 9; idig++) {
			if (map[idig] & bit) {// this is the digit
				//				if (FD[idig].Off(res))  return 1; // invalid, gane locked
				//				Seta(idig, res);
				int cell = res;
				Assign(idig, cell);
				goto nextr1;// finished for that cell
			}
		}
		return 1; //conflict with a previous cell assign
	nextr1: {}
	}
	return 0;// not locked 
}
void ZHONE::Seta(int digit, int cell) { // single in cell
	Assign(digit, cell);
	register  int  bit= 1 << cell,cl = ~bit;
	FD[0]&=cl; FD[1] &= cl; FD[2] &= cl; FD[3] &= cl; FD[4] &= cl;
	FD[5] &= cl; FD[6] &= cl; FD[7] &= cl; FD[8] &= cl;
	FD[digit] |= bit;// restore digit
}
int ZHONE::Update() {
	int Shrink = 1;
	register int S, A, cl;
	register uint32_t *wcl = FD;
	while (Shrink) {
		Shrink = 0;
		if ((A = FD[0]) - CompFD[0]) { UPDN1(0, 1, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[1]) - CompFD[1]) { UPDN1(1, 0, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[2]) - CompFD[2]) { UPDN1(2, 0, 1, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[3]) - CompFD[3]) { UPDN1(3, 0, 1, 2, 4, 5, 6, 7, 8); }
		if ((A = FD[4]) - CompFD[4]) { UPDN1(4, 0, 1, 2, 3, 5, 6, 7, 8); }
		if ((A = FD[5]) - CompFD[5]) { UPDN1(5, 0, 1, 2, 3, 4, 6, 7, 8); }
		if ((A = FD[6]) - CompFD[6]) { UPDN1(6, 0, 1, 2, 3, 4, 5, 7, 8); }
		if ((A = FD[7]) - CompFD[7]) { UPDN1(7, 0, 1, 2, 3, 4, 5, 6, 8); }
		if ((A = FD[8]) - CompFD[8]) { UPDN1(8, 0, 1, 2, 3, 4, 5, 6, 7); }
		//	  Debug(1);
	}// end while

	return 1;
}
void ZHONE::Guess() {
	if (!cells_unsolved) {
		if (zh1b_g.zsol && (!zh1b_g.nsol)) SetKnown(zh1b_g.zsol);// store the first solution
		zh1b_g.nsol++;
		if (zh1b_g.nsol > zh1b_g.lim) zh1b_g.go_back = 1;
		return;
	}
	uint32_t res;
	int ndig = 2;
	register int R3 = zh1b_g.pairs;
	if (bitscanforward(res, R3))goto gogo; // first pair is ok to go
	ndig = 3;
	R3 = zh1b_g.triplets;
	bitscanforward(res, R3); // first triplet is ok to go
gogo:
	int cell = res, bit = 1 << cell;
	//cout << "brute force guess  cell " << cellsFixedData[cell].pt << endl;
	for (int idig = 0; idig < 9; idig++) {
		if (FD[idig] & bit) {// one valid digit		
			if (--ndig) {
				ZHONE * mynext = this + 1; // start next guess
				*mynext=*this;
				mynext->Seta(idig, res);
				mynext->ComputeNext();
				if (zh1b_g.go_back) return;
			}
			else {// this is the last do it in the same place
				Seta(idig, res);
				ComputeNext();
				return;
			}
		}
	}
}
int ZHONE::FullUpdate() {
	if (zh1b_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())			return 0; // locked empty cell or conflict singles in cells
		if (zh1b_g.single_applied) {
			//			cout << "after singles" << endl; Debug();
			continue;
		}
		break;
	}
	return 1;
}

//___________________________________________________
char * ZHONE::SetKnown(char * zs) {
	strcpy(zs, &empty_puzzle[27]);
	zs[27] = 0;
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++) {// one digit
		int arows = FD[idig] >> 27;// 3 rows
		if (arows == 7) continue;
		int band = FD[idig];
		for (int j = 0; j < 3; j++) if (!(arows & (1 << j))) {
			int row = (band >> TblMult9[j]) & 0x1ff;
			uint32_t  irow;
			bitscanforward(irow, row);
			int	cell = TblMult9[j] + irow;
			zs[cell] = idig + '1';
		}
	}
	return zs;
}
void ZHONE::Debug(int all) {
	cout << "DEBUG  nbsol=" << zh1b_g.nsol << " unsolved=" << Unsolved_Count() 
		<< " index to zhone="<<this-zhone<< endl;
	//cout << zh1b_g.out27 << " band1 " << endl;
	char zi[82];  SetKnown(zi);
	cout << zi << " known rows 1_3 digits " << endl;
	if (!all) return;

	cout << "map per digit one band" << endl;
	for (int ir = 0; ir < 3; ir++) {
		for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++) {
			int vf = FD[idig], wr = (vf >> (9 * ir)) & 0x1ff;
			for (int k = 0; k < 9; k++) {
				if (wr & (1 << k))		cout << idig + 1;
				else 		cout << ".";
				if (k == 2 || k == 5) 	cout << " ";
			}
			cout << "  ";
		}
		cout << endl; //end of row
	}
	cout << endl; // end of map per digit

}
void ZHONE::DebugDigit(int digit) {
	cout << "DEBUG  digit=" << digit + 1 << endl;
	for (int ir = 0; ir < 3; ir++) {
		int vf = FD[digit], wr = (vf >> (9 * ir)) & 0x1ff;
		for (int k = 0; k < 9; k++) {
			if (wr & (1 << k))		cout << digit + 1;
			else 		cout << ".";
			if (k == 2 || k == 5) 	cout << " ";
		}
		cout << endl; //end of row
	}
	cout << endl; // end of block

}
int ZHONE::GetAllDigits(int cell) {
	int ir = 0, bit = 1 << cell;
	for (int i = 0; i < (int)zh1b_g.ndigits; i++) if (FD[i] & bit)ir |= (1 << i);
	return ir;
}
void ZHONE::ImageCandidats() {
	int dig_cells[27]; for (int i = 0; i < 27; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " << endl << endl;
	for (i = 0; i < 9; i++) {  // column
		lcol[i] = 2;    // 2  mini 
		for (j = 0; j < 3; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 3; i++) { // now row 
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell], 
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < (int)zh1b_g.ndigits; id++)
				if (digs & (1 << id))
					if(zh1b_g.ndigits==9)	cout << id + 1;
					else cout << zh1b_g.digmap [id] + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}

//============== collector mode
void ZHONE::Checkstart() {
	zh1b_g.ndigits = 9;
	memcpy(FD, zh1b_g.fd_sols[1], sizeof FD);
	cells_unsolved = BIT_SET_27 ;
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++)
		FD[idig] |= 7 << 27;
	ImageCandidats();
}
void ZHONE::InitGuess() {// after start floors before guess
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++)
		FD[idig] |= 7 << 27;
	memset(zh1b_g.previous_ua_status, 0,
		sizeof zh1b_g.previous_ua_status); //no previous status
	for (int i = 0; i < 6; i++)
		zh1b_g.upstream_unsolved_cells[i] = BIT_SET_27;
	if (zh1b_g.diag) {
		cout << "end Init guess" << endl;
		ImageCandidats();
		cout << "end Init guess bis" << endl;
		ImageCandidats();
	}
}
int ZHONE::UpdateDigit(int digit) {
	register int S, A = FD[digit], cl;
	if (A  != CompFD[digit]) {
		int Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		if ((A >> 27) != S) {
			cl = ~(A & TblRowMask[S]);
			A = (A&BIT_SET_27) | (S << 27);
			cells_unsolved &= cl;
		}
		FD[digit] = CompFD[digit] = A;
	}// end if
	return 1;
}
void ZHONE::Set2(int cell) { // single in cell
	Assign(0, cell);
	register  int  cl = ~(1 << cell);
	FD[1] &= cl;
}

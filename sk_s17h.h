//#define DODEBUG
//#define DO2X2
#define TESTCLEAN

struct SPOT_E64 {// spots to find band12 valid solutions n clues
	SPOT_E64 * sp;
	uint64_t  all_previous_cells, active_cells;
	uint32_t * start_possibles, n_possibles, ipos, ispot;
	uint64_t * tua;
	uint32_t stack[3], bands[2], missing_clues, nua;
	inline void Copy(SPOT_E64 * old) {
		*this = *old;
		start_possibles += n_possibles;
		ispot++;
		missing_clues--;
		ipos = 0;
		tua += nua;
	}
	inline void AddCellBandStack(int cell, uint32_t * ncb) {
		// if the stack is limit update sn active
		int st = C_stack[cell];
		stack[st]++;
		if (stack[st] > 5) {
			//cout << "stack pleine" << st << endl;
			active_cells &= ~band3xBM[st + 3].u64[0];

		}
		// if the band is limit update sn active
		int b = C_div27[cell];
		bands[b]++;
		if (bands[b] > 5) {// more in mode b 656 only
			//cout << "bande pleine" << b<< endl;
			active_cells &= ~band3xBM[b].u64[0];
		}
	}
	inline void GetOne(uint64_t v) {
		n_possibles = 1;
		bitscanforward64(start_possibles[0], v);
		start_possibles[0] = From_128_To_81[start_possibles[0]];
	}
	inline void AddPossibles(uint64_t v) {
		uint32_t cc;
		while (bitscanforward64(cc, v)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			v ^= (uint64_t)1 << cc;// clear bit
			start_possibles[n_possibles++] = From_128_To_81[cc];
		}
	}
	inline int GetUa(uint64_t v) {
		n_possibles = 0;
		AddPossibles(v);
		return n_possibles;
	}
	uint32_t GetPossibles() {
		//p_cpt2g[22] ++;
		//cout << p_cpt2g[22] << " possibles missing  " << missing_clues << endl;
		// UAs are limited to active cells no empty or single ua
		if (missing_clues < 2) return 0; // minimum to call this process
		uint32_t cells_count[64],
			min_elims = (nua + missing_clues - 1) / missing_clues;
		memset(cells_count, 0, sizeof cells_count);
		uint32_t cc;
		for (uint32_t iua = 0; iua < nua; iua++) {
			register uint64_t Rw = tua[iua] & BIT_SET_2X;
			while (bitscanforward64(cc, Rw)) {// look for  possible cells
				register uint64_t bit2 = (uint64_t)1 << cc;
				Rw ^= (uint64_t)1 << cc;// clear bit
				cells_count[cc]++;
				//if (!cc) cout << "cell 0 pour i=" << iua << " cpt="<< cells_count[cc] << endl;
			}
		}
		//cout << "compte brut critical "<< min_elims << endl;
		//for (int i = 0; i < 64; i++)if (cells_count[i])
		//	cout << From_128_To_81[i] << "\t" << cells_count[i] << endl;

		// collect cells over critical count
		GINT64 tcells[64], temp;
		uint32_t ntcells = 0;
		for (int i = 0; i < 54; i++) {
			register uint32_t my_cell_count = cells_count[C_To128[i]];
			if (my_cell_count >= min_elims) {
				GINT64 & w = tcells[ntcells++];
				w.u32[0] = i;
				w.u32[1] = my_cell_count;
			}
		}
		if (!ntcells) return 0;
		if (ntcells > 1) {// sort in decreasing order
			for (uint32_t i = 0; i < ntcells - 1; i++) {
				for (uint32_t j = i + 1; j < ntcells; j++) {
					if (tcells[i].u64 < tcells[j].u64) {
						temp.u64 = tcells[i].u64;
						tcells[i].u64 = tcells[j].u64;
						tcells[j].u64 = temp.u64;
					}
				}
			}
			//if (ntcells > 64 - missing_clues) ntcells = 64 - missing_clues;
		}
		// load the final table of cells to consider
		//cout << "final count" << endl;
		for (uint32_t i = 0; i < ntcells; i++) {
			start_possibles[i] = tcells[i].u32[0];
			//cout <<i<<"\t"<< tcells[i].u32[0]<<"\t"<< tcells[i].u32[1] <<endl;
		}
		n_possibles = ntcells;
		return ntcells;
	}
	inline int GetLast() {
		n_possibles = 0;
		uint64_t andx =(uint64_t) BIT_SET_2X;
		for (uint32_t iua = 0; iua < nua; iua++) {
			andx &= tua[iua];
			if (!andx)return 0;
		}
		uint32_t cc;
		while (bitscanforward64(cc, andx)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			andx ^= (uint64_t)1 << cc;// clear bit
			start_possibles[n_possibles++] = From_128_To_81[cc];
		}
		return n_possibles;
	}
	inline uint32_t GetPossiblesP() {// using parallel
		//p_cpt2g[24] ++;
		uint32_t lim = (nua + missing_clues - 1) / missing_clues,
			lim2 = lim + 5;
		//uint64_t tval[17];
		uint64_t tval[14];
		memset(tval, 0, sizeof tval);
		tval[1] = tua[0] & BIT_SET_2X;
		uint32_t tend = 1; // last used so far
		for (uint32_t iua = 1; iua < nua; iua++) {
			register uint64_t Rw = tua[iua] & BIT_SET_2X,
				Rmore = tval[tend] & Rw;
			if (tend < lim2 && Rmore) tend++;
			switch (tend) {
				//case 15:tval[15]|=tval[14] & Rw;
				//case 14:tval[14]|=tval[13] & Rw;
			case 13:tval[13] |= tval[12] & Rw;
			case 12:tval[12] |= tval[11] & Rw;
			case 11:tval[11] |= tval[10] & Rw;
			case 10:tval[10] |= tval[9] & Rw;
			case 9:tval[9] |= tval[8] & Rw;
			case 8:tval[8] |= tval[7] & Rw;
			case 7:tval[7] |= tval[6] & Rw;
			case 6:tval[6] |= tval[5] & Rw;
			case 5:tval[5] |= tval[4] & Rw;
			case 4:tval[4] |= tval[3] & Rw;
			case 3:tval[3] |= tval[2] & Rw;
			case 2:tval[2] |= tval[1] & Rw;
			}
			tval[1] |= Rw;
		}
		if (tend < lim)return 0;// nothing to do
		n_possibles = 0;
		for (uint32_t i = lim; i < tend; i++) tval[i] &= ~tval[i + 1];
		for (uint32_t i = tend; i >= lim; i--)AddPossibles(tval[i]);
		//if (p_cpt2g[23] < 10)
		//	cout << p_cpt2g[23] << " possibles 2 missing " << missing_clues << " n="<<n<< endl;
		return n_possibles;
	}
	inline uint32_t GetPossibles2() {// using parallel
		//p_cpt2g[23] ++;
		register uint64_t R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R2) return 0;
		R2 &= ~R3; R3 &= ~R4; R4 &= ~R5;
		n_possibles = 0;
		if (R5)					AddPossibles(R5);
		if (R4)					AddPossibles(R4);
		if (R3)					AddPossibles(R3);
		if (R2)					AddPossibles(R2);
		return n_possibles;
	}
	inline uint32_t GetPossibles3() {// using parallel
		//p_cpt2g[23] ++;
		register uint64_t R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R3) return 0;
		R3 &= ~R4; R4 &= ~R5; R5 &= ~R6;
		n_possibles = 0;
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		if (R4)	AddPossibles(R4);
		if (R3)	AddPossibles(R3);
		return n_possibles;
	}
	inline uint32_t GetPossibles4() {// using parallel
		//p_cpt2g[23] ++;
		register uint64_t R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R4) return 0;
		R4 &= ~R5; R5 &= ~R6;
		n_possibles = 0;
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		if (R4)	AddPossibles(R4);
		return n_possibles;
	}
	inline uint32_t GetPossibles5() {// using parallel
		//p_cpt2g[23] ++;
		register uint64_t R7 = 0, R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R7 |= R6 & R; R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R5) return 0;
		R4 &= ~R5; R5 &= ~R6; R6 &= ~R7;
		n_possibles = 0;
		if (R7)	AddPossibles(R7);
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		return n_possibles;
	}
	void D1() {
		cout << "D1\t" << ispot << "\t" << ipos << endl;
		//cout << Char2Xout(active_cells) << " spot=" << ispot << " pos=" << ipos
		//	<< " npos=" << n_possibles << endl;
	}
	void D2(int all=0) {
		//cout << Char2Xout(active_cells) << " Shrink active"  << endl;
		cout << Char2Xout(all_previous_cells) << " known nuas=" << nua << endl;
		if (!all) return;
		for (uint32_t iua = 0; iua < nua; iua++)
			cout << Char2Xout(tua[iua]) << endl;
	}

	void D3() {
		cout << Char2Xout(all_previous_cells) << " known" << endl;
		cout << "get possibles  nua=" << nua << "  missing_clues " << missing_clues
			<< "n_possibles" << n_possibles << endl;
	}
	inline void Ddead(uint32_t iua) {
		//			cout <<"\t\t"<< ispot <<" "<<ipos<<" dead branch iua=" << iua << endl;
	}
	inline void Dass(uint32_t iua, uint64_t Ru) {
		//			cout<<"\t\t" << ispot << " " << ipos << " assign iua=" << iua
		//				<< " cell=" << start_possibles[0] << endl;
		//			cout << ispot << " " << ipos << " assign iua=" << iua <<endl
		//				<< Char2Xout(Ru) <<" cell="<< start_possibles[0] << endl;
	}
	inline void DNoMore() {
		cout << Char2Xout(all_previous_cells) << "no more uas "
			<< bands[0] << bands[1] << " "
			<< stack[0] << stack[1] << stack[2] << endl;

	}
};
struct MINCOUNT {
	uint32_t mini_bf1, mini_bf2, mini_bf3, mini_triplet,
		critbf, pairs27;
	uint32_t all_used_minis, mincount,minplus;
	inline void SetMincount() {// after direct setting minis
		all_used_minis = mini_bf1 | mini_triplet;
		mini_triplet &= ~mini_bf1;// count only triplets with no pair
		mincount = _popcnt32(all_used_minis) + _popcnt32(mini_bf3);
		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs

		// set up pair + triplet bitfield
		if (mini_triplet) {// must add triplet minirow
			for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
				if (mini_triplet&bit)
					critbf |= field;
		}
		minplus = mincount;
	}
	GINT64_t  Count_per_stack() {
		GINT64 cc; cc.u64 = 0;
		for (int i = 0, st = 0111; i < 3; i++, st <<= 1) {
			cc.u16[i]= _popcnt32(all_used_minis&st) + 
				_popcnt32(mini_bf3&st);
		}
		return cc;
	}
	void Status(const char * lib) {
		cout<<lib << "critical Status mincount ="<<mincount<< " minplus=" <<minplus << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char9out(mini_bf1) << "     minis bf1" << endl;
		cout << Char9out(mini_bf2) << "     minis bf2" << endl;
		cout << Char9out(mini_bf3) << "     minis bf3" << endl;
		cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;

	}
};
struct G17TMORE {// FIFO table of more for bands 1+2
	uint64_t  t[G17MORESIZE];
	int nt, maxt, curt;
	inline void Init() { maxt = G17MORESIZE; nt = 0; }
	inline void Add(uint64_t v) {//add a new more in FIFO 
		if (nt < maxt) {// use next location
			curt = nt;
			t[nt++] = v;
		}
		else {// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			t[curt] = v;
		}

	}
	inline void Add_If_New(uint64_t v) {// check oldest first
		register uint64_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// from curt to 0
		if ((*Rt) == V)return;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) goto add;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if ((*Rt) == V)return;
	add:
		Add(v);
	}
	inline int Check(uint64_t v) {// check oldest first
		if (!nt) return 0;
		register uint64_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if (!((*Rt) & V))return 1;
		return 0;
	}
	void Print(int modegua) {
		register uint64_t *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		{
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return;

		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl) {
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}

	}


	void PrintUasDirect() {
		for (int i = 0; i < nt; i++) {
			register uint64_t w = t[i];
			cout << Char2Xout(w) << endl;
		}
	}

};
struct MORE32 {// FIFO table of more for band b
	uint32_t  t[32];
	int nt, maxt, curt;
	inline void Init() { maxt = 32; nt = 0; }
	inline void Add(uint32_t v) {//add a new more in FIFO 
		if (nt < maxt) {// use next location
			curt = nt;
			t[nt++] = v;
		}
		else {// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			t[curt] = v;
		}

	}

	inline int Check(uint32_t v) {// check oldest first
		if (!nt) return 0;
		register uint32_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if (!((*Rt) & V))return 1;
		return 0;
	}

};
  
// standard first band (or unique band)
struct STD_B416 {
	char band[28];
	int band0[27], i416, gangster[9], map[27], dband;
	uint32_t tua[100], nua;//   maximum 81  
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void Initstd();
	void GetBandTable(int i);
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas()
		;
	void InitC10(int i);
	void InitG12(int i);
	void InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
		, int iband = 1);

	
	void PrintStatus();
};
struct STD_B1_2 :STD_B416 {
	// row solution pattern in digit
	int mini_digs[9], mini_pairs[27],
		revised_g[9];// after false forced in 1/2 minirows
	int  tv_pairs[27], nvpairs; //9 or 27 bits 
	void FillMiniDigsMiniPairs(STD_B1_2 & bb);
	inline void InitRevisedg() {
		memcpy(revised_g, gangster, sizeof gangster);
	}
	int ReviseG_triplet(int imini, int ip, STD_B1_2 * bb);
	uint32_t GetMiniData(int index, uint32_t & bcells, STD_B1_2 *bb);
	void PrintShortStatus();
};
struct STD_B3 :STD_B416 {// data specific to bands 3
	struct GUAs {
		BF128 isguasocket2, isguasocket3, isguasocket2_46;// active i81
		BF128 isguasocketc2, isguasocketc3, isguasocketc2_46;// active i81
		int triplet[9];//same gua3s
		int triplet_imini[81];
		int ua_pair[81], ua_triplet[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81],
			ua2pair27[81],ua2bit[81],ua3bit[81];
	}guas;
	BF128 isguasocketc246;//all linked to a socket 2
	BF128 issocket2x2;
	uint32_t pat2x2[81];
	uint32_t nvalid6, nvb128,idiag,idiag127;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	//_______________ handling mincount
	uint32_t t2[81], nt2 , t3[9], nt3 ,	t46[81], nt46,t246[81],nt246,
		andmiss1,noutmiss1,wactive0;

	//____  handling socket vectors (link i81 to V128b3s[256][1858][27 + 9]    
	uint32_t tsock[MAXSOCKB3], ntsock;// tsock is the pattern to use
	int index2[81], index3[81], index4_6[81], index246[81]; // index to tsock

	struct NOT_EMPTY {// blocs of 128 passing filters
		BF128 v;
		uint64_t index, filler;
	}bnot_empty[1858]; 
	uint32_t  nbnot_empty;// blocs 128 valids passing filters


	MINCOUNT smin;
	GINT64  stack_count;// after BuildPossiblesB3()
	//_______________________


	void InitBand3(int i16, char * ze, BANDMINLEX::PERM & p);
	int IsGua(int i81);
	int IsGua3(int i81);
	int Is2Rows(int i81_1, int i81_2);
	int GetI81_2(int bf) {
		for (int i = 0; i < 81; i++) if (guas.ua_pair[i] == bf) return i;
		return -1;
	}
	int GetI81_3(int bf) {
		for (int i = 0; i < 81; i++) if (guas.ua_triplet[i] == bf) return i;
		return -1;
	}
	void Build_tsock();

	int Clean_valid_bands3A();
	void Clean_valid_bands3B(uint32_t ib3);
	void Clean_valid_b128(uint32_t ib3);
	void ExpandB3();
	void SetUpMincountxy(BF128 & validsockets);
	void SetUpMincountb1b2();
	uint32_t QuickCheckMore();
	void PrintB3Status();
	void DiagExpand(int ib3) ;
};

//================== UA collector 2 bands 

struct GENUAS_B12 {// for uas collection in bands 1+2 using brute force 
	int dig_cells[9][9],
		gangbf[9],// columns for band 3 in bit field
		revised_gangbf[9],// same revised UA2s UA3s ***
		mini_digs[9], mini_pairs[27], // UA2s UA3  ***
		//valid_pairs, //  27 bits valid sockets UA2s ***
		nfloors, limstep, map[9], cptdebug, modemore;
	BF128 valid_sockets;

	//=============== uas collector 
	int limsize, floors;
	uint64_t  tuaold[1000],// previous non hit uas infinal table of uas for bands 1+2
		tua[TUA64_12SIZE]// 
		,tuab1b2[200]	;// collecting bands uas in 2x mode
	uint32_t nuaold, nua, nuab1b2,
		tuamore[500];
	//_______________uamore control
	STD_B1_2 *ba, *bb;
	uint32_t patb, ib, digp,colb, cola;
	uint64_t w0, ua,p12;
	//_____________________ functions collect UAs bands 1+2
	int Initgen();
	void BuildFloorsAndCollectOlds(int fl);
	//int AddUA64(uint64_t * t, uint32_t & nt);
	inline void AddUA(uint64_t v) {
		ua = v; AddUA64(tua, nua, ua);
	}
	inline void AddUACheck(uint64_t v) {
		if (nua >= TUA64_12SIZE) nua = TUA64_12SIZE - 1;
		ua = v; AddUA64(tua, nua, ua);
	}
	int BuilOldUAs(uint32_t r0);
	int CheckOld();
	int CheckMain(uint64_t wua);
	void CollectMore2digits();
	void Collect2digits2_4_cols();
	void CollectMoreTry6_7();
	void EndCollectMoreStep();
	void CollectTriplets();
	void CollectMore2minirows();
	//_____________________ functions collect UA2s UAs3 socket 

	void ProcessSocket2(int i81);
	int DebugUas();
};

#define SIZETGUA 150
#define GUAREDSIZE 100
struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[512];
	int modeb12, go_back, diagmore,diagbug,ip20,
		it16, it16_2, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16, maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[82], tc[6], ntc;
	int gcheck[82], ib2check, ib3check;
	int skip, last;// restart point; last entry in the batch
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2, pband3, pcheck2, pcheck3;
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;
	int cold[9],coldf[9], rowd[6], boxd[6], rowdb3[3], boxdb3[3]; //free digits 
	//_________________ gangster 
	int gangcols[9];// 9 bits fields band3 for each colum (see valid band2)
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int gangb12[9]; // digit bf bands 12 per column
	int   *gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// active cols for a given digit
	//____________structs hosting the 81 GUA entries
	struct SGUA2 {// 81 possible UA2 sockets
		// permanent data
		uint64_t * tua;
		int col1, col2;// columns of the socket
		int i_81,iguan; // index 0_80 for this 
		int i9;// complementary column in minirow
		int id1, id2; // index of digits in gang 27 
		// Current band1+2 data
		int digs, dig1, dig2;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		int gangcols[9];// revised gangster
		uint32_t nua;// nua_start, nua_end;
		void Debug(const char * lib);

	}tsgua2[81];
	struct SGUA3 {// 81 possible UA3 sockets
		// permanent data
		uint64_t * tua,killer;
		int col1;// first columns 0-9 
		int i_81, imini,iguan; // index 0_80 for this 
		int id1, id2, id3; // index of digits in gang 27 
		// Current band1+2 data
		int  dig1, dig2, dig3;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		uint32_t nua;// nua_start, nua_end, nua;
		void Debug(const char * lib);
	}tsgua3[81];
	// __________________________  primary UAs tables and creation of such tables
	uint64_t  // tua3x[3000],// dynamic sub tables
		*ptua2;// pointer to current table cycle search 2/3
	uint32_t  ntua2, ntua3, nua2; // nua2 for cycle search 2/3  
	//================== bands 3 and gangster band 3 analysis
	int nband3;
	int tactive2[81], nactive2, tactive3[81], nactive3;
	int   tcolok[2], ncolok;

	int ngua6_7, c1, c2, band, floors, digp, i81;
	uint64_t wua0, ua;// partial gua ua to check
	uint64_t tuacheck[100], tua_for_check[500];
	uint32_t uadigs_for_check[500], nua_for_check, nua_check;
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues

	GEN_BANDES_12() {
		gang27 = gang[0];
		InitialSockets2Setup();
		InitialSockets3Setup();
	}
	void InitialSockets2Setup();// batch level
	void InitialSockets3Setup();// batch level
	void BuildGang9x3();
	void Build_CheckUAs_Subsets_Table();
	void Build_tuacheck(int fl);
	int Have_tuacheck_subset();
	void SecondSockets2Setup();// band1+2 level
	void SecondSockets2MoreUAs();// band1+2 level
	void GuaCollectMore();
	void SecondSockets3Setup();// band1+2 level
	void GuaCollect(int fl, int diag = 0);
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode = 0);
	void NewBand1(int iw);
	int Band2Check();
	int Band3Check();
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10 = 1);
	int DebugFreshUA(uint64_t ua);
	int Debug17(SGUA2 & w);
	//int FindBand3Unique();//test or  debugging code see the corresponding file
	//================ B creating a catalogue for the 17 search 
	//same as A exchanging bands 2/3


	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
	// debugging code special call
	//int Afterb1b2(int option = 0);
};

struct G17B3HANDLER {
	MINCOUNT smin;
	int known_b3, rknown_b3, active_b3,   ib3, nb3,
		active_sub, ndead, wactive0, nmiss, //ncritical,
		irloop, wua, stack;
	uint32_t *uasb3if, nuasb3if, *uasb3of, nuasb3of , andoutf;
	GINT64 stack_count;
	int diagh;
	// ================== entry in the proces
	void Init(STD_B3 & b3);
	void GoMiss0(STD_B3 & b3 );
	void GoMiss1(STD_B3 & b3);
	void Do_miss1();
	void GoMiss2Init(STD_B3 & b3);
	void GoMiss2(STD_B3 & b3, uint32_t uamin);
	void Init( );
	inline void AddCellMiss1(uint32_t cell, int bit) {
		stack_count.u16[C_stack[cell]]++;
		known_b3 |= bit;
	}
	void AddCell_Miss2(uint32_t * t);
	inline int AddCell_Of(uint32_t cell, int bit) {
		register int s= C_stack[cell];
		if (stack_count.u16[s] > 5) return 0;
		stack_count.u16[s]++;
		if (stack_count.u16[s] > 5) {
			s =~( 07007007 << (3 * s));// mask
			wua &= s;
			wactive0 &= s;
		}
		nmiss--;
		known_b3 |= bit;
		return 1;
	}
	uint32_t IsMultiple(int bf);
	int ShrinkUas1();
	//=============== process critical
	void CriticalAssignCell(int Ru);
	void Critical2pairs();
	void Go_Critical();
	void CriticalLoop();
	//==================== process subcritical no cell added outside the GUAs field
	void SubMini(int M, int mask);
	void Go_Subcritical();
	void Go_SubcriticalMiniRow();
	void Go_SubcriticalMiniRow_End(int stack);
	//===================== process not critical
	void ShrinkUasOfB3();
	void Go_miss1_b3();
	void Go_miss2_b3();
	void Go_miss3_b3();

};

/*
entry 92
maxindex= 983
maxn5= 51516
maxn6= 237762
maxdet5= 261
maxdet6= 2004
*/
struct GUAN {// one occurrence of the active guan 
	uint64_t *tua, killer;
	uint32_t colbf, digsbf, ncol, // 2-4 cols
		nua, nfree;
	int i81;// i81 or pattern (if 4 columns)
	void Enter(uint64_t *t, uint32_t n, uint32_t cbf,
		int32_t dbf, uint32_t ind) {
		i81 = ind;// index 0-80 if gua2 gua3
		tua = t; nua = n; colbf = cbf; digsbf = dbf;
		killer = BIT_SET_2X;
		for (uint32_t i = 0; i < nua; i++)killer &= tua[i];
		ncol = _popcnt32(colbf);
	}
	void Debug1Guan(int i, int all = 0) {
		cout << Char2Xout(killer);
		cout << "kill  i=" << i << "\tcols" << Char9out(colbf);
		cout << "\tdigs" << Char9out(digsbf) << " ncol=" << ncol
			<< " nua=" << nua << " i81=" << i81
			<< endl;
		if (all) {
			for (uint32_t i = 0; i < nua; i++)
				cout << Char2Xout(tua[i]) << endl;
		}
	}

};


struct TU_GUAN {// GUAN process (used GUA all kinds)
	GUAN tguan[256], guanw, tguanb1[160],g2[160];
	uint64_t  guabuf[15000], *pguabuf;
	uint64_t  guabufb1[15000], *pguabufb1;
	uint64_t  guabufr[10000], *pguabufr;
	uint32_t //i81socket2[80], i81socket3[80],
		tsockets2[80], tsockets3[80],// tsockets4[80],
		ntsockets2, ntsockets3, //ntsockets4,
		ntsockets2_2, ntsockets3_2;// , ntsockets4_2;
	//______________________ valid handler
	//uint64_t vv2, vv3, vv2_2, vv3_2;
	uint32_t nguan, nguanb1, nguastepb2, ng2 ;
	void AddGuan(uint64_t *t, uint32_t n, uint32_t cbf,
		int32_t dbf, uint32_t ind) {
		if (nguan < 256) {
			pguabuf = &pguabuf[n];// lock space
			tguan[nguan++].Enter(t, n, cbf, dbf, ind);
		}
		//tguan[nguan - 1].Debug1Guan(nguan - 1);
	}

	void Init(){
		pguabuf = guabuf;// reinit gua buffer use
		nguan = 0;//and guan table
	}
	void DoStepb1();// reduce tables 
	void DoStepb2();// reduce tables again
	//void SetValidVector();
	void SetValidV23();




	void Debug1() {
		for (uint32_t i = 0; i < nguan; i++)
			tguan[i].Debug1Guan(i,1);
	}
	void Debug1_for_Bf12(uint64_t bft);
	void Debug2() {
		for (uint32_t i = 0; i < nguanb1; i++)
			tguanb1[i].Debug1Guan(i, 1);
	}
	void Debug3() {
		cout << "debug tuguan stepb2  total guan" <<	ng2   << endl;
		cout << ntsockets2_2 << " " << ntsockets3_2 << " known sockets 2 3" << endl;
		for (uint32_t i = 0; i < ng2; i++)
			g2[i].Debug1Guan(i, 1);
		cout << endl;
	}
	void Debugsetv23() {

		cout << ntsockets2_2 << " " << ntsockets3_2
			<< " sockets début" << endl;
		cout << ntsockets2 << " " << ntsockets3
			<< " sockets fin" << endl;
		cout << "sockets2 ";
		for (uint32_t i = 0; i < ntsockets2; i++)
			cout << " " << tsockets2[i];
		cout <<endl<< "sockets3 ";
		for (uint32_t i = 0; i < ntsockets3; i++)
			cout << " " << tsockets3[i];
		cout << endl;
	}
};
struct BI2 {//Index 2 valid in band
	uint64_t bf, // 2 cells in bit fiekd
		active; // remaining possible cells
	uint32_t tval[2], // 2 cells in int mode
		istart,iend;// in the VALIB table

 };
struct VALIDB {//  valid band 1+2
	uint64_t bf; // 2 cells in bit fiekd
		//active; // remaining possible cells
	uint32_t tval[4]; // 0-4 cells in int mode
	uint64_t nval;// n cells over the 2
	inline void Enter(uint64_t ebf, uint32_t * tc2) {
		bf = ebf ;
		nval= _popcnt64(bf)-2;
		memcpy(tval, tc2, sizeof tval); 
		//Print();
	}
	int DebugFindKnown17();
	void Print() {// debugging
		cout << Char2Xout(bf) << " valid n="<<nval << endl;
	}
};
struct VALIDB1 {// minimal valid band no index
	uint64_t bf; // 2 cells in bit fiekd
	uint32_t tval[6]; // 0-6 cells in int mode
	uint64_t nval;// n cells over the 2
	void Print() { // debugging
		cout << Char2Xout(bf) << " valid n=" << nval << endl;
	}
};

#define ZST6 800
#define ZST5 500
#define ZST4 300
#define ZST3 50
#define ZST2 5
struct SPLITSTEPB1 {
	VALIDB1 * va0, *vaa;
	uint64_t bf,a_orall,b_orall;
	uint32_t n0, na;
	SPLITSTEPB1();
	void BuildStep( uint32_t  n0e);
	int GetA(uint64_t bfe) {// debugging known 17
		for (uint32_t i = 0; i < na; i++)
			if (bfe == vaa[i].bf)return i;
		return -1;
	}
	int GetB(uint64_t bfe) {// debugging known 17
		for (uint32_t i = 0; i < n0; i++)
			if (bfe == va0[i].bf)return i;
		return -1;
	}
};
struct ZSTEP {
	VALIDB t6[ZST6], t5[ZST5], t4[ZST4], t3[ZST3], t2[ZST2];
	uint64_t nt6, nt5, nt4, nt3,nt2;
	void CopyIf(ZSTEP & o,uint64_t * tua, uint32_t nua) {// must hit all tus
		nt6= nt5= nt4= nt3=nt2=0;// source must hit all uas to be valid
		for (uint64_t i = 0; i < o.nt6; i++) {
			register uint64_t Rua = o.t6[i].bf;
			for (uint32_t j = 0; j < nua; j++) {
				if(!(Rua&tua[j]))goto nextt6;
			}
			t6[nt6++] = o.t6[i];
			nextt6:;
		}
		for (uint64_t i = 0; i < o.nt5; i++) {
			register uint64_t Rua = o.t5[i].bf;
			for (uint32_t j = 0; j < nua; j++) {
				if (!(Rua&tua[j]))goto nextt5;
			}
			t5[nt5++] = o.t5[i];
		nextt5:;
		}

	}
	int DebugFindKnown17();

	void Print() {// debugging
		cout << "print t.." << nt6 << " " << nt5  << endl;
		for (uint64_t i = 0; i < nt6; i++)
			cout << Char2Xout(t6[i].bf) << endl;
		for (uint64_t i = 0; i < nt5; i++)
			cout << Char2Xout(t5[i].bf) << endl;
	}
};
struct ZS64 {// final loop 64 bits 
	uint64_t bf,  v;
};
struct ZS128 {// final loop 64 bits 
	BF128 v;
	uint64_t bf, active;
};
struct ZS256 {// final loop 64 bits 
	BF128 v,v2;
	uint64_t bf, bfm;// bfm not in the step for >256
};
struct ZS384 {// final loop 64 bits 
	BF128 v, v2,v3;
	uint64_t bf, bfm;// bfm not in the step for >256
};
struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	BF128 p17diag;// known 17 pattern for tests
	int loopb1,b3lim, debug17,debug17_check,
		diag, diagbug,debugb3, aigstop, 
		iretb1,doloopnotok,
		npuz, a_17_found_here;
	uint32_t	iband1,iband2, step1count;
	uint64_t breakcount,b2count, totb2;
	G17B3HANDLER g17hh0;
	//______sockets common to  all bands 3  
	BF128 isguasocket2all, isguasocket3all;
#ifdef DO2X2
	BF128 socks2x2;// sockets 2x2  4 cells  2 digits in band3
	int ts2x2[81], nts2x2; // storing active sockets 2x2
	uint64_t tuas2x2[81]; // first ua (usually one) for an active socket 2x2
	int ts2x2_clean[81], nts2x2_clean; // same status for a given XY 
#endif

	//====== data for band expansion
	uint32_t nexp, bnua;
	uint64_t btua[300],start_active, b1cpt[8], b2cpt[8],b1cptdiag;
	BI2 * mybi2t,wbi_1,wbi_2;
	VALIDB * myokt, validb_known,validb_known_b1, validb_known_b2;
	uint32_t nmybi2t, nmyokt,nbi2_1,nbi2_2,
		nvb1,nvb1steps[10],
		nzs1_6,nzs2_6, nzs1_5, nzs2_5;
	ZSTEP zstep_all_b1[250],stepb2,stepb1,zfb1,zfb2;
	//======= status after step 2 in band 2 then 2 uas in band 1
	SPLITSTEPB1 splitstepb1;
	uint64_t tusb1[2000], tusb2_12[1000], tusb2[1000], temptyb1[500];
	uint32_t ntusb1, ntusb2, ntusb2_12, nemptyb1;
	VALIDB1 * tvb1go;// part A or part B
	uint32_t ntvb1go;
	uint64_t  temptyb2[500];
	uint32_t  ntua_128, ntua_256, nemptyb2;
	uint64_t fb1,acb1, fb2, fb12, acb2a,  acb2, acb12;

	//============ vectors 64 128 bits 
	BF128 v128uas,  vc128[54];
	BF128 v256uas, vc256[54];// 128 to 256 uas
	uint64_t v64uas, vc64[54],nmainloop11;
	ZS64 zs64_1_6[MAXSTEP6], zs64_1_5[MAXSTEP5],
		zs64_2_6[MAXSTEP6], zs64_2_5[MAXSTEP5];
	ZS128 zs128_1_6[MAXSTEP6], zs128_1_5[MAXSTEP5],
		zs128_2_6[MAXSTEP6], zs128_2_5[MAXSTEP5];
	ZS256 zs256_1_6[MAXSTEP6], zs256_1_5[MAXSTEP5],
		zs256_2_6[MAXSTEP6], zs256_2_5[MAXSTEP5];
	ZS384 zs384_1_6[MAXSTEP6], zs384_1_5[MAXSTEP5],
		zs384_2_6[MAXSTEP6], zs384_2_5[MAXSTEP5];

	//============================ b12 no more uas to test
	G17TMORE moreuas_AB, moreuas_AB_small, moreuas_AB_big;
	uint64_t wb12bf, wb12active,myua;
	GINT64  stack_count_step, stack_count, stack_countf;
	uint32_t tclues[40],tb3[256],*tcluesxy; 
	int nclues_step, nclues,ntb3;
	BF128 bands_active_pairs, bands_active_triplets,
		valid_vect;
	//============  go band3
	uint32_t free1, free2, free3;

	//int cur_ib;
	uint32_t tcluesb12[20], ncluesb3x;
	uint32_t   nmiss;
	uint32_t uasb3_1[2000], uasb3_2[2000], uas_in[2000],
		nuasb3_1, nuasb3_2, nuas_in, b3_andout;
	MINCOUNT smin;
	MORE32 moreuas_b3, moreuas_b3_small;

	
	
	
	
	
	
	//=====================process for a new band 2 / set of bands 3
	void GoM10();// standard entry
	void GoM10Known();// entry for a known 17
	void GoM10Uas();// collect uas guas
	void GoM10GUas2x2();// collect guas 2 soclets 2 digits
	void GoM10B3s();// expand B3 and vectors B3
	void BuildSocketsB1(uint64_t active);
	void ExpandB1();
	void BuildNoIndexB1();
	void ExpandB2();
	void ExpandB3(uint32_t ib3);
	void ExpandOneBand(int ib);// find 2-6 valid bands no redundant clue
	void DoLoopB1(VALIDB1 * t1, uint32_t nt1);
	void DoStepB1_From_Do_Loop(VALIDB1 * t1, uint32_t nt1);
	void Apply_Band1_Loop();// extract as killed in loop
	void Apply_Band1_Step();// shrink tables
	int Apply_Band2_Step();// shrink tables again
			//________ end of step band 1
	void ShrinkBand2();

					  
					  //___________ outer loop on steps band2 band1

		//___________ extract potential valid bands 1+2 (no more uas)
	void Do64uas();
	void DoChunk64(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb);
	void Do64uas_11(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb);

		//________ same 65 to 128 uas
	void Do128uas();
	void DoChunk128(ZS128 * a, ZS128 * b, uint64_t na, uint64_t nb);
	void Do128uas_11(ZS128 * a, ZS128 * b, uint64_t na, uint64_t nb);
	
	//________ same >128 uas
	void Do256uas();
	void DoChunk256(ZS256 * a, ZS256 * b, uint64_t na, uint64_t nb);
	void Do256uas_11(ZS256 * a, ZS256 * b, uint64_t na, uint64_t nb);

	//_______ processing potential valid bands 1+2
	void CleanAll();
	void B12InTclues();
	inline void AddXClue(uint32_t * t, int & n, uint32_t xcell) {
		uint32_t cell = From_128_To_81[xcell];
		stack_count.u16[C_stack[cell]]++;
		t[n++] = cell;
	}
	int Is_B12_Not_Unique();
	int BuildPossiblesB3();
	void GoPossiblesB3();
	void GoB3(STD_B3 & b);
	void DebugAdd12();

	void Clean_11(uint64_t bf);
	void Clean_11_n(uint64_t bf);
	void Clean_guas_first();
	void Clean_guas_first2();
	void Clean_valid_first();
	void Clean_valid_first2();

	void FinalCheckB3(uint32_t bfb3);
	void Out17(uint32_t bfb3);
	void MergeUasExpand(STD_B3 & b);
	void ExpandBand3(uint32_t *tua, uint32_t nua);
	void Debug_If_Of_b3();




	//================ debugging code
	inline uint32_t k17x(int ix) {	return p17diag.bf.u32[ix];	}
	uint64_t DebugFindKnown17_valids_2bands();
	void Debug_b1b2cpt();
	void DebugGetPuz(const char * p) {
		p17diag.SetAll_0();
		for (int i = 0; i < 81; i++)
			if (p[i] != '.')p17diag.Set_c(i);

		cout <<"this is a "
			<<_popcnt32(p17diag.bf.u32[0])
			<< _popcnt32(p17diag.bf.u32[1])
			<< _popcnt32(p17diag.bf.u32[2]) 
			<<" pattern for the expected puzzle"<< endl;
	}
	//int DebugK17M10();
	void GodebugInit(int mode);
	//int GodebugFindKnown17();
	int GodebugCheckUas(const char * lib);
};
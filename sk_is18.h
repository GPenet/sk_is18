#define ZST6 800
#define ZST5 500
#define ZST4 300
#define ZST3 50
#define ZST2 5
struct MINCOUNT {
	uint32_t mini_bf1, mini_bf2, mini_bf3, mini_triplet,
		critbf, pairs27;
	uint32_t all_used_minis, mincount, minplus;
	void SetMincount() {// after direct setting minis
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
	inline void SetMinplus(int nof, int andw) {
		minplus = mincount;
		if (nof) {
			minplus++;
			if (!andw) minplus++;
		}

	}
	void Status(const char * lib) {
		cout << lib << "critical Status mincount =" << mincount << " minplus=" << minplus << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char9out(mini_bf1) << "     minis bf1" << endl;
		cout << Char9out(mini_bf2) << "     minis bf2" << endl;
		cout << Char9out(mini_bf3) << "     minis bf3" << endl;
		cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;

	}
}smin, sminb2;
//================ Bands 
struct BI2 {//Index 2 valid in band
	uint32_t bf, // 2 cells in bit fiekd
		active; // remaining possible cells
	uint32_t tval[2], // 2 cells in int mode
		istart, iend;// in the VALIB table

};
struct VALIDB {//  valid band 1+2
	uint32_t bf; // 2 cells in bit fiekd
		//active; // remaining possible cells
	uint32_t tval[5]; // 0-5 cells in int mode
	uint32_t nval;// n cells over the 2
	inline void Enter(uint32_t ebf, uint32_t * tc2) {
		bf = ebf;
		nval = _popcnt32(bf) - 2;
		memcpy(tval, tc2, sizeof tval);
	}
};
struct VALIDB64 {// valid band  mode 64 bits
	uint64_t bf; // 2 cells in bit fiekd
	uint32_t tval[7]; // 0-7 cells in int mode
	uint32_t nval;// n cells 
};


struct STD_B416 {
	BI2 * my_bi2;
	VALIDB * my_validb;
	//VALIDB1 * my_validb1;
	uint32_t nbi2, nvalidb;//, nvalidb1;
	char band[28];
	int band0[27], i416, gangster[9], map[27], dband;
	uint32_t tua[100], nua;//   maximum 81  
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas()	;
	void InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
		, int iband = 1);
	void InitExpand(BI2 * bi2, VALIDB * validb) { my_bi2 = bi2; my_validb = validb; }
	void ExpandOneBand();

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

#define MAXSOCKB3 100

struct STD_B3 :STD_B416 {// data specific to bands 3
	struct GUAs {
		BF128 isguasocket2, isguasocket3, isguasocket2_46;// active i81
		int triplet[9];//same gua3s
		int triplet_imini[81];
		int ua_pair[81], ua_triplet[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81],
			ua2pair27[81], ua2bit[81], ua3bit[81];
	}guas;
	BF128 isguasocketc246;//all linked to a socket 2
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	//_______________________

	void EndInitBand3();
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
}myband3;

struct BINDEXN {
	BI2 * t2;
	VALIDB * tvb;
	uint32_t nt2, ntvb;
	inline void Attach(BI2 * t2e, VALIDB * tvbe) { t2 = t2e; tvb = tvbe; }
	void Copy(STD_B1_2 & b);
	void Copy(BINDEXN & b);
	void Copy_no7clues(STD_B1_2 & b);

}bin_b1,bin_b2,bin_b1yes, bin_b2yes,
bin2_b1, bin2_b2, bin2_b1yes, bin2_b2yes;
struct INDEX_XY {
	uint64_t bf,and_g,or_g;// 2 cells common to the lot
	uint64_t ntotvb,ncluesmin;
	struct ITEM {// one of the tables per lot 
		VALIDB64 * tvb;
		uint32_t ntvb,sum_vb;
	}titem[5];
	void Debug() {
		cout << "\tnand=" << _popcnt64(and_g);
		cout << "\tnor=" << _popcnt64(or_g);
		for(int i=0;i<5;i++)
			cout << "\t" << titem[i].ntvb <<";" << titem[i].sum_vb;
	}

}index_xy_b1,index_xy_b2;
BF128 uas_buffer[20000],uasnn_buffer[30000];


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
struct MOREALL {// FIFO table of more for bands 1+2
	BF128 ta[G17MORESIZE];
	uint64_t  t[G17MORESIZE];
	int nt, maxt, curt;
	inline void Init() { maxt = G17MORESIZE; nt = 0; }
	inline void Add(BF128 v) {//add a new more in FIFO 
		if (nt < maxt) {// use next location
			curt = nt;
			ta[nt++] = v;
		}
		else {// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			ta[curt] = v;
		}

	}
	inline int Check(uint64_t v) {// check oldest first
		if (!nt) return 0;
		register uint64_t V = v;
		register BF128 *Rt = &ta[curt], *Rtl;
	loop1:// form curt to 0
		if (!(Rt->bf.u64[0] & V))return 1;
		if (--Rt >= ta)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &ta[curt];
		Rt = &ta[maxt];
		while (--Rt > Rtl)if (!(Rt->bf.u64[0] & V))return 1;
		return 0;
	}
	void PrintUasDirect() {
		for (int i = 0; i < nt; i++) {
			BF128 wa = ta[i];
			cout << Char27out(wa.bf.u32[2])<<"\t";
			cout << Char2Xout(wa.bf.u64[0]) << endl;
		}
	}
};

struct GVCELLS {// vectors cells for guas
	BF128 v21[27][54]; // first 128 uas GUA2
	BF128 v22[27][54]; // next 128 uas GUA2
	BF128 v31[9][54]; // first 128 uas GUA3
	BF128 v32[9][54]; // next 128 uas GUA3
	BF128 vpx[256][54]; // room for 256 atterns 4/6 in band3
}gvcells;
struct GVSTEPS {// base vectors start,loopb1,loopb2
	BF128 v21[27];
	BF128 v22[27];
	BF128 v31[9];
	BF128 v32[9];
	BF128 vpx[256];
}gvs_start,gvs_b1,gvs_b2;

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
	uint64_t  tuaold[1000],// previous non hit uas in final table of uas for bands 1+2
		tua[TUA64_12SIZE]// 
		, tuab1b2[200];// collecting bands uas in 2x mode
	uint32_t nuaold, nua, nuab1b2,
		tuamore[500];
	//_______________uamore control
	STD_B1_2 *ba, *bb;
	uint32_t patb, ib, digp, colb, cola;
	uint64_t w0, ua, p12;
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
struct GEN_BANDES_12 {// encapsulating global data 
	int modeb12, go_back, diagmore, diagbug, ip20,
		it16, it16_2, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[82], tc[6], ntc;
	int gcheck[82], ib2check, ib3check;
	int skip, last;// restart point; last entry in the batch
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2, pband3, pcheck2, pcheck3;
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;
	int cold[9], coldf[9], rowd[6], boxd[6], rowdb3[3], boxdb3[3]; //free digits 
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
		int i_81, iguan; // index 0_80 for this 
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
		uint64_t * tua, killer;
		int col1;// first columns 0-9 
		int i_81, imini, iguan; // index 0_80 for this 
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
	void SetUpGuas2();
	void SetUpGuas3();
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

};

struct TU_GUAN {// GUAN process (used GUA all kinds)
	GUAN tguan[256], guanw ;
	uint64_t  guabuf[30000], *pguabuf;
	uint32_t nguan, nguanb1, nguastepb2, ng2;
	void AddGuan(uint64_t *t, uint32_t n, uint32_t cbf,
		int32_t dbf, uint32_t ind) {
		if (nguan < 256) {
			pguabuf = &pguabuf[n];// lock space
			tguan[nguan++].Enter(t, n, cbf, dbf, ind);
		}
	}

	void Init() {
		pguabuf = guabuf;// reinit gua buffer use
		nguan = 0;//and guan table
	}

}tuguan;

struct G17B3HANDLER {
	MINCOUNT smin;
	int known_b3, rknown_b3, active_b3, ib3, nb3,
		active_sub, ndead, wactive0, nmiss, //ncritical,
		irloop, wua, stack;
	uint32_t *uasb3if, nuasb3if, *uasb3of, nuasb3of, andoutf;
	int diagh;
	// ================== entry in the proces
	void Init();
	void GoMiss0();
	void GoMiss1();
	void Do_miss1();
	void GoMiss2Init();
	void GoMiss2( uint32_t uamin);
	//void Init();
	void AddCell_Miss2(uint32_t * t);
	inline int AddCell_Of(uint32_t cell, int bit) {
		nmiss--;
		known_b3 |= bit;
		return 1;
	}
	uint32_t IsMultiple(int bf);
	//int ShrinkUas1();
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
struct ZS64 {// final loop 64 bits 
	uint64_t bf, v;
};
struct ZS128 {// final loop 64 bits 
	BF128 v;
	uint64_t bf, active;
	void dump() {
		cout << Char2Xout(bf) << " ";
		cout << Char64out(v.bf.u64[0]) ;
		cout << Char64out(v.bf.u64[1])  << endl;
	}
};
struct ZS256 {// final loop 64 bits 
	BF128 v, v2;
	uint64_t bf, bfm;// bfm not in the step for >256
};
struct ZS384 {// final loop 64 bits 
	BF128 v, v2, v3;
	uint64_t bf, bfm;// bfm not in the step for >256
};

struct GCHK {
	//=========== kwon puzzle filters
	BF128 puzknown,puzknown_perm;
	uint32_t kpfilt[4]; //initial statu 3 
	int kn_ir1, kn_ir2;


	int aigstop, start_perm, *tpw, *tsortw,
		a_18_seen;
	uint32_t band_order[3];

	//___________________ studied solution 
	char * ze;// given solution grid
	char * zp;// first 18 if any
	char zsol[82];// solution grid morphed
	int grid0[81];// same 0 based
	STD_B416 * bands_abc[3],bA,bB;

	//_______________ GUAs 
	MINCOUNT smin;
	uint32_t guapats[256], nguapats,// active not gua2s gua3s
		guapats2[27], guapats3[9],
		nclues_b3,wactive0,
		tpatx[128],npatx,
		andmiss1, noutmiss1;
	void GuapatsInit() {
		for (int i = 0,mask=7,ij=0; i < 9; i++,mask<<=3) 
			guapats3[i ] = mask;
		for (int i = 0; i < 27; i++) {
			int j = i / 3, mask = guapats3[j], bit = 1 << i;
			guapats2[i] = mask ^ bit;
			
		}
		/*
		for (int i = 0; i < 27; i++) {
			cout << Char27out(guapats2[i]) << " ";
			cout << Char27out(gua2bit[i]) << " i=" << endl;

		}
		for (int i = 0; i < 9; i++)
			cout << Char27out(guapats3[i]) << endl;

		*/
	}


	/*
	//_____ sockets 2x2
	BF128 socks2x2;// sockets 2x2  4 cells  2 digits in band3
	int ts2x2[81], nts2x2; // storing active sockets 2x2
	uint64_t tuas2x2[81]; // first ua (usually one) for an active socket 2x2
	int ts2x2_clean[81], nts2x2_clean; // same status for a given XY
	*/
	//___________________ external loop 
	int loopb1;
	uint32_t  b2count, breakcount;
	struct EXTL {
		uint64_t noxyes;
		uint32_t bfx, tbfy[20], ntbfy,mode,ratio;
		void Init(uint32_t bf, int mod) { bfx = bf; ntbfy = 0; mode = mod; }
		void Debug() {
			if (mode == 1) {
				cout << Char27out(bfx) << " bfx b1 ntbfy= " << ntbfy << endl;
				for (uint32_t i1 = 0; i1 < ntbfy; i1++) {
					cout << "\t\t" << Char27out(tbfy[i1]) << " bfy" << endl;
				}
			}
			else {
				cout <<"\t\t"<< Char27out(bfx) << " bfx b2 ntbfy= " << ntbfy << endl;
				for (uint32_t i1 = 0; i1 < ntbfy; i1++) {
					cout << Char27out(tbfy[i1]) << " bfy" << endl;
				}

			}

		}
	}extl1[10],extl2[10],extlw,extlr;
	uint32_t nextl1, nextl2;
	uint64_t n_nob1 = 0, ntotb1 = 0, n_nob2 = 0, ntotb2 = 0, minratio;
	uint64_t FindSockets(uint64_t active, uint64_t lim);
	void ExtractMin(uint64_t active, BINDEXN & bin1, BINDEXN & bin2);
	void ExtSplitY(BINDEXN & binw, uint32_t *tbf, uint32_t ntbf, uint32_t & activer);
	void ExtSplitX(BINDEXN & bin1no, BINDEXN & bin1yes,
		uint32_t bf, uint32_t & activer);
	//_______________  loops XY 
	G17TMORE moreuas_12_13, moreuas_14, moreuas_15, moreuas_AB, moreuas_AB_small, moreuas_AB_big;
	MOREALL mall9, mall10, mall11, mall12, mallxx;// uas seen in 18 check per size
	uint64_t tusb1[2000], tusb2_12[1000], tusb2[1000]; 
	uint32_t ntusb1, ntusb2, ntusb2_12; 
	uint32_t ntvb1go;
	//uint32_t   ntua_128, ntua_256;
	uint64_t fb1, acb1, fb2, fb12,  acb2, acb12;
	uint32_t tclues[40], *tcluesxy;// mini 25+band a
	int nclues_step, nclues;
	uint64_t n_to_clean;

	int G3_SplitBi2( int mode,int kill7,
		BINDEXN & binw ,uint32_t ibi2,
		INDEX_XY & indxyw,VALIDB64 * pvb );

	void G3_SplitAll(int mode, BINDEXN & binw,  
		INDEX_XY & indxyw, VALIDB64 * pvb);
	
	//============ main loop process
	BF128 v128uas, vc128[54];
	BF128 v256uas, vc256[54];// 128 to 256 uas
	uint64_t v64uas, vc64[54];
	BF128 v64_192uas, vc64_192[54],bv192;
	BF128 v192_320uas, vc192_320[54],bv320;
	BF128 v320_448uas, vc320_448[54],bv448;
	BF128 v448_576uas, vc448_576[54],bv576;
	BF128 v576_704uas, vc576_704[54],bv704;
	BF128 v704_832uas, vc704_832[54],bv832;
	BF128 v832_960uas, vc832_960[54],bv960;



	uint32_t nua_add, nua_3x3y, nxy_filt1,nvalid,n128add;

	//================== clean process
	uint64_t wb12bf, wb12active, myua;



	//________ clean and valid
	uint32_t uasb3_1[2000], uasb3_2[2000], uas_in[2000],
		nuasb3_1, nuasb3_2, nuas_in, b3_andout;

	int Clean_valid_bands3A();
	void Clean_valid_bands3B();

	
	//======================= band b  index 3 reduction
	//uint32_t  indf;
	

	//=============== band B when band A is locked 
	//BF128  final81_2, final81_3;
	//BANDB sbb;
	//_____  initial infield outfield and more  
	uint32_t btuaif[256], btuaof[3000], tuaif[3000],
		nbif, nbof;
	uint32_t more_of[128], nmoreof, more_if[128], nmoreif;
	uint32_t  mode_ab, myuab, filt32, 
		mincluesb,maxcluesb,	mincluesc, maxcluesc ;
	MORE32 moreuas_b3;
	//==================== current band 3 to process
	uint32_t  ncluesb3;
	int   nmiss;

	//____ start the process for a triplet banda band b bandc 
	void Start(STD_B416 * bandx, int * tsort,int ip);
	void Copy_Check_7clues_needed();
	void Go1_Collect_Uas();
	void Go1GUas2x2();
	void Go2_Ext_Loop();
	void Go2b_Ext_Loop(uint64_t activeloop,uint32_t mode2);


	//___ process sub lots 
	void Go3(BINDEXN & bin1, BINDEXN & bin2);
	void Apply_Band1_Step();
	void Apply_Band1_Step_GUAs();
	int Apply_Band2_Step();
	void Apply_Band2_Step_GUAs();
	int Apply_Band1_2_Step();
	//___________ extract potential valid bands 1+2 (no more uas)
	void Do64uas();
	void DoChunk64(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb);
	void Do64uas_11(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb);
	

	//_______ processing potential valid bands 1+2
	void CleanAll();
	void ExpandB3();
	int Is_B12_Not_Unique();
	void FinalCheckB3(uint32_t bfb3);
	void Out17(uint32_t bfb3);
	void NewUaB3();

};


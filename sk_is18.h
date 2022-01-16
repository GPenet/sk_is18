#define ZST6 800
#define ZST5 500
#define ZST4 300
#define ZST3 50
#define ZST2 5

struct TUAS81 {
	BF128 tall[5000],
		told[3000], // later  has band 3
		t2_2[1000];//band 3 2/3 after 5 clues 
	uint32_t ntall, ntold, nt2_2;
	int Add(BF128 w, uint32_t floor) {
		uint32_t cc = w.Count96(), nfloor = ~floor;
		for (uint32_t iua = 0; iua < ntall; iua++) {
			BF128 wt = tall[iua];
			uint32_t cct = wt.bf.u32[3] >> 16;
			uint32_t floort = wt.bf.u32[3] & 0777;
			wt.bf.u32[3] = 0;
			if (cct < cc) {// look for subset
				if (floort&nfloor) continue;//can not be subset
				if ((wt - w).isEmpty()) return 0; // subset
				continue;
			}
			if (wt == w) return 0; // redundancy
			if (cct == cc)continue; // not sorted here
			// insert here and check super set
			for (uint32_t jua = ntall; jua > iua; jua--)
				tall[jua] = tall[jua - 1];
			tall[iua] = w;// new inserted
			tall[iua].bf.u32[3] = floor | (cc << 16);
			ntall++;
			// is it a subset of a previous entry
			for (iua++; iua < ntall; iua++) {
				if ((w - tall[iua]).isEmpty()) {// we have a subset
					for (uint32_t k = iua + 1; k < ntall; k++)
						tall[k - 1] = tall[k];
					ntall--;
					iua--; //continue same position
				}
			}
			return 2;
		}
		w.bf.u32[3] = floor | (cc << 16);
		tall[ntall++] = w;// new added
		return 1;
	}
	int Add2(BF128 w, uint32_t cc) {
		for (uint32_t iua = 0; iua < ntall; iua++) {
			BF128 wt = tall[iua];
			uint32_t cct = wt.bf.u32[3];
			wt.bf.u32[3] = 0;
			if (cct < cc) {// look for subset
				if ((wt - w).isEmpty()) return 0; // subset
				continue;
			}
			if (wt == w) return 0; // redundancy
			if (cct == cc)continue; // not sorted here
			// insert here and check super set
			for (uint32_t jua = ntall; jua > iua; jua--)
				tall[jua] = tall[jua - 1];
			tall[iua] = w;// new inserted
			tall[iua].bf.u32[3] = cc;
			ntall++;
			// is it a subset of a previous entry
			for (iua++; iua < ntall; iua++) {
				if ((w - tall[iua]).isEmpty()) {// we have a subset
					for (uint32_t k = iua + 1; k < ntall; k++)
						tall[k - 1] = tall[k];
					ntall--;
					iua--; //continue same position
				}
			}
			return 2;
		}
		w.bf.u32[3] = cc;
		tall[ntall++] = w;// new added
		return 1;
	}
	int New(BF128 w, uint32_t floor) {
		told[ntold++] = w;
		return Add(w, floor);
	}
	void SetupOld(uint32_t floors) {
		ntold = 0;
		register uint32_t  F = ~floors;
		for (uint32_t i = 0; i < ntall; i++) {
			BF128 w = tall[i];
			uint32_t fl = w.bf.u16[6];
			if (!(F&fl))told[ntold++] = w;
		}
	}
	void Debug(BF128 *t, uint32_t nt, char * lib) {
		cout << lib << " " << nt << endl;
		for (uint32_t iua = 0; iua < nt; iua++) {
			BF128 wt = t[iua];
			cout << Char27out(wt.bf.u32[0]) << "|";
			cout << Char27out(wt.bf.u32[1]) << "|";
			cout << Char27out(wt.bf.u32[2]) << " digs=";
			cout << Char9out(wt.bf.u16[6])
				<< " i=" << iua
				<< " " << _popcnt64(wt.bf.u64[0])
				<< " " << _popcnt32(wt.bf.u32[2]) << endl;
		}
	}
	void DebugAll() {
		Debug(tall, ntall, "all uas status nuas=");
	}
	void DebugOld() {
		Debug(told, ntold, "old uas status nuas=");
	}

};
struct TUASB12 {
	uint64_t tua[TUA64_12SIZE], t2[3000];
	uint32_t nua, nt2;
	void SortBySize() {
		uint64_t t[50][300], nt[50];
		memset(nt, 0, sizeof nt);
		for (uint32_t iua = 0; iua < nua; iua++) {
			register uint64_t wu = tua[iua] & BIT_SET_2X,
				cc = _popcnt64(wu);
			t[cc][nt[cc]++] = wu;
		}
		nua = 0;
		for (uint32_t i = 0; i < 50; i++)if (nt[i]) {
			uint64_t *to = t[i];
			for (uint32_t j = 0; j < nt[i]; j++)
				tua[nua++] = to[j];
		}
	}
	uint32_t CountT2(uint64_t filter) {
		nt2 = 0;
		register uint64_t F = filter;
		for (uint32_t iua = 0; iua < nua; iua++)
			if (!(tua[iua] & F))nt2++;
		return nt2;
	}
	uint32_t SetupT2(uint64_t filter, uint64_t active) {
		nt2 = 0;
		register uint64_t F = filter, w;
		for (uint32_t iua = 0; iua < nua; iua++) {
			w = tua[iua];
			if (!(w&F))t2[nt2++] = w & active;
		}
		return nt2;
	}
	void GetT2(uint64_t filter) {
		nt2 = 0;
		register uint64_t Fn = ~filter, w;
		for (uint32_t iua = 0; iua < nua; iua++) {
			w = tua[iua] & BIT_SET_2X;
			if (!(w&Fn))t2[nt2++] = w;
		}
	}
	void GetT2(uint64_t f,uint64_t *tui27, uint64_t nui27 ) {
		nt2 = 0;
		register uint64_t Fn = ~f, w;
		for (uint32_t iua = 0; iua < nua; iua++) {
			w = tua[iua] & BIT_SET_2X;
			if (!(w&Fn))t2[nt2++] = w;
		}
		for (uint64_t iua = 0; iua < nui27; iua++) {
			w = tui27[iua] & BIT_SET_2X;
			if (!(w&Fn))t2[nt2++] = w;
		}
	}
	void Dump() {
		cout << "dump tua nua=" << nua << endl;
		for (uint32_t iua = 0; iua < nua; iua++)
			cout << Char2Xout(tua[iua]) << " i=" << iua << " " << (tua[iua] >> 59) << endl;;
	}
	void DumpT2() {
		cout << "dump t2 nua=" << nt2 << endl;
		for (uint32_t iua = 0; iua < nt2; iua++)
			cout << Char2Xout(t2[iua]) << " i=" << iua << endl;;
	}
};

struct CHUNK3B {// storing 64 uas and vectors 3 bands
	uint64_t tu54[64], v0, vc[54], nt;
	uint32_t tu27[64];// index 0-26 or pattern
	inline void Init() {
		nt=0,		v0 = 0;
		memset(vc, 255, sizeof vc);
	}
	inline void Add(uint64_t u54, uint32_t u27) {//add a new ua
		if (nt >= 64) return;// safety should never be
		uint64_t bit = (uint64_t)1 << nt;
		tu27[nt] = u27;
		tu54[nt++] = u54;
		v0 |= bit;
		uint32_t cc54;// build cells vectors
		register  uint64_t Rw = u54;
		while (bitscanforward64(cc54, Rw)) {
			Rw ^= (uint64_t)1 << cc54;// clear bit
			vc[cc54]^=bit;
		}
	}
	inline uint64_t ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		if (!nt) return 0;
		uint64_t w = v0;
		for (uint32_t i = 0; i < ntcells; i++)
			w &= vc[tcells[i]];
		return w;
	}
	void Get1(BF128 * td, uint32_t&  ntd,
		uint64_t bf54, uint64_t ac54) {
		for (uint32_t i = 0; i < nt; i++) {
			register uint64_t U = tu54[i],F= bf54;
			if (!(F&U)) {
				BF128 & w = td[ntd++];
				w.bf.u64[0] = U & ac54;
				w.bf.u32[2] = tu27[i];
			}
		}
	}
	void Get2(BF128 * td, uint32_t&  ntd,
		uint32_t * tc, uint32_t ntc,
		uint64_t ac54) {
		register uint64_t V = v0;
		uint32_t ir;
		for (uint32_t i = 0; i < ntc; i++)
			V &= vc[tc[i]];
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			BF128 & w = td[ntd++];
			w.bf.u64[0] = tu54[ir] & ac54;
			w.bf.u32[2] = tu27[ir];
		}
	}
	void DebugC2(int all=1) {
		for (uint32_t i = 0; i < nt; i++) {
			cout << Char54out(tu54[i]) << " i27=" << tu27[i] << " i=" << i << endl;
		}
		if (!all) return;
		cout << Char64out(v0) << " v0" << endl;
		for(int i=0;i<54;i++) 
			if ((vc[i] & v0) != v0)
				cout << Char64out(vc[i] & v0) << " cell=" <<i<< endl;
	}
	void DebugMore(int all = 1) {
		for (uint32_t i = 0; i < nt; i++) {
			cout << Char54out(tu54[i]) << " ";
			cout << Char27out(tu27[i]) << " i=" << i << endl;
		}
		if (!all) return;
		cout << Char64out(v0) << " v0" << endl;
		for (int i = 0; i < 54; i++)
			if ((vc[i] & v0) != v0)
				cout << Char64out(vc[i] & v0) << " cell=" << i << endl;
	}
};
struct CHUNK1B {//storing 64 uas and vectors band 3s
	uint64_t v0, vc[27],nt;
	uint32_t tua[64];
	inline void Init() {
		nt = 0, v0 = 0;
		memset(vc, 255, sizeof vc);
	}
	inline void Add(uint32_t v) {//add a new ua
		if (nt >= 64) return;// safety should never be
		uint64_t bit = (uint64_t)1 << nt;
		v &= BIT_SET_27;//no extra bit
		tua[nt++] = v;
		v0 |= bit;
		uint32_t cc27;// build cells vectors
		register  uint32_t Rw = v;
		while (bitscanforward(cc27, Rw)) {
			Rw ^= 1 << cc27;// clear bit
			vc[cc27] ^= bit;
		}
	}
	inline uint64_t ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		if (!nt) return 0;
		uint64_t w = v0;
		for (uint32_t i = 0; i < ntcells; i++)
			w &= vc[tcells[i]];
		return w;
	}
	void Debug(int all=1) {
		for (uint32_t i = 0; i < nt; i++) {
			cout << Char27out(tua[i]) << " i=" << i << endl;
		}
		if (!all) return;
		cout << Char64out(v0) << " v0" << endl;
		for (int i = 0; i < 27; i++) 
			if( (vc[i] & v0)!=v0)
			cout << Char64out(vc[i] & v0) << " cell=" << i << endl;
	}
}b3direct;
struct CHUNKS_HANDLER{
	CHUNK3B c2[40],  cmore[20];
	CHUNK1B band3[2], misc[2];// up to 81 band3 up to 128 others
	uint32_t ic2,  icmore, iband3, imisc;
	void Init() {
		ic2=  icmore= iband3= imisc=0;
		c2[0].Init();  cmore[0].Init();
		band3[0].Init(); band3[1].Init();
		misc[0].Init();
	}
	inline int GetC2Count() { return 64 * ic2 +(int) c2[ic2].nt; }
	void Add128(BF128 w);// switch to index 0-26
	void Addc2(uint64_t ua12, uint32_t ua) {
		uint64_t ua54 = (ua12 & BIT_SET_27) |
			((ua12 & BIT_SET_B2) >> 5);
		if (ic2 == 39 && c2[39].nt > 63) return;
		if (c2[ic2].nt > 63)c2[++ic2].Init();
		c2[ic2].Add(ua54, ua);
	}
	void Addband3(uint32_t w) {
		if (iband3 == 1 && band3[1].nt > 63) return;
		if (band3[iband3].nt > 63)band3[++iband3].Init();
		band3[iband3].Add(w);
	}
	void Addmisc(uint32_t w) {
		if (imisc == 1 && misc[1].nt > 63) return;
		if (misc[imisc].nt > 63)misc[++imisc].Init();
		misc[imisc].Add(w);
	}

	void GetMore(BF128 * td, uint32_t&  ntd,
		uint64_t bf54, uint64_t ac54,
		uint32_t * tc, uint32_t ntc ) {
		ntd = 0;
		for (uint32_t i = 0; i <= icmore; i++)
			if (cmore[i].nt < 15)
				cmore[i].Get1(td, ntd, bf54, ac54);
			else cmore[i].Get2(td, ntd, tc, ntc, ac54);

	}
	void Get2(BF128 * td, uint32_t&  ntd,
		uint64_t bf54, uint64_t ac54,
		uint32_t * tc, uint32_t ntc) {
		ntd = 0;
		for (uint32_t i = 0; i <= ic2; i++)
			if (c2[i].nt < 15)
				c2[i].Get1(td, ntd, bf54, ac54);
			else c2[i].Get2(td, ntd, tc, ntc, ac54);

	}
	int C2Count() { return(int) (64 * ic2 + c2[ic2].nt); }
	void DebugAll() {
		cout <<"guas 2 cells)"<<endl;
		for (uint32_t i = 0; i <= ic2; i++)
			c2[i].DebugC2(0);
		cout << "guas more cells)" << endl;
		for (uint32_t i = 0; i <= icmore; i++)
			cmore[i].DebugMore(0);
		cout << "band 3)" << endl;
		for (uint32_t i = 0; i <= iband3; i++)
			band3[i].Debug(0);

	}
	void C2Status(){
		for (uint32_t i = 0; i <= ic2; i++)
			c2[i].DebugC2(0);
	}

}chunkh;
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
		cout << lib << "critical Status mincount =" << mincount << endl;
			//<< " minplus=" << minplus << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		if (mini_bf1)cout << Char9out(mini_bf1) << "     minis bf1" << endl;
		if (mini_bf2)cout << Char9out(mini_bf2) << "     minis bf2" << endl;
		if (mini_bf3)cout << Char9out(mini_bf3) << "     minis bf3" << endl;
		//cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;

	}
};
//================ Bands 

#define myband3 bax[2]

struct STD_B416 {
	char band[28];
	int band0[27], i416, gangster[9], map[27], dband;
	uint32_t tua[82], nua;//   maximum 81  
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas()	;
	void InitBand2_3( char * ze, BANDMINLEX::PERM & p
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


struct MORE64VECT {// FIFO table of more for bands 1+2
	BF128 vect, vc[54];
	uint64_t  t[128];
	uint32_t nt;
	inline void Init() {
		nt = 0;
		memset(&vect, 0, sizeof vect);
		memset(vc, 255, sizeof vc);

	}

	inline void Add(uint64_t v) {//add a new more if <128 
		if (nt < 128) {
			uint32_t cc64;// build cells vectors 
			vect.Set(nt);
			register uint64_t Rw = v;
			while (bitscanforward64(cc64, Rw)) {
				Rw ^= (uint64_t)1 << cc64;// clear bit
				if (cc64 > 26)cc64 -= 5;
				vc[cc64].clearBit(nt);
			}
			t[nt++] = v;
		}
	}
	inline void Add2(uint64_t v, uint64_t ac) {
		uint32_t cc64;// build cells vectors 
		vect.Set(nt);
		register uint64_t Rw = v & ac;
		while (bitscanforward64(cc64, Rw)) {
			Rw ^= (uint64_t)1 << cc64;// clear bit
			if (cc64 > 26)cc64 -= 5;
			vc[cc64].clearBit(nt);
		}
		t[nt++] = v;
	}
	inline int ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		if (!nt) return 0;
		BF128 w = vect;
		for (uint32_t i = 0; i < ntcells; i++)
			w &= vc[tcells[i]];
		return (w.isNotEmpty());
	}
	void Add_to(uint64_t * td, uint32_t & ntd) {
		for (uint32_t i = 0; i < nt; i++)
			td[ntd++] = t[i];
	}


};
struct MOREV2 {// 2 more64vect paired
	MORE64VECT mv1, mv2;
	uint32_t sw12;
	inline void Init() { mv1.Init(); mv2.Init(); sw12 = 0; }
	inline void Add(uint64_t v) {//add a new more if <128 
		if (sw12 && (mv2.nt == 128)) {
			mv1.Init();
			sw12 = 0;// mv1 active
		}
		if ((!sw12) && (mv1.nt == 128)) {
			mv2.Init();
			sw12 = 1;// mv2 active
		}
		if (sw12)mv2.Add(v);
		else mv1.Add(v);
	}
	inline int ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		if (mv1.ApplyXY(tcells, ntcells)) return 1;
		return mv2.ApplyXY(tcells, ntcells);
	}

	void Status(const char * lib) {
		cout << "status for more uas " << lib << endl;
		cout << mv1.nt << " " << mv2.nt << " " << mv1.nt + mv2.nt << endl;
	}
	void Add_to(uint64_t * td, uint32_t & ntd) {
		mv1.Add_to(td, ntd);
		mv2.Add_to(td, ntd);
	}

}morev2a, morev2b, morev2c, morev2d;

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

	uint32_t Check(uint32_t v) {// check oldest first
		if (!nt) return 0;
		register uint32_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if (!((*Rt) & V))return (*Rt);
		return 0;
	}

};
struct MOREVALID {// adds for a valid band 1+2
	BF128 ta[128];
	int nt;
	inline void Init() {  nt = 0; }
	inline void Add(BF128 v) { //add a new more in FIFO 
		if (nt < 128) ta[nt++] = v;
	}	
}morevalidc2,morevalidothers;
struct G17B3HANDLER {
	MINCOUNT smin;
	int known_b3, rknown_b3, active_b3, ib3, nb3,
		active_sub, ndead, wactive0, nmiss, //ncritical,
		irloop, wua, stack;
	uint32_t *uasb3if, nuasb3if, *uasb3of, nuasb3of, andoutf;
	int diagh;
	// ================== entry in the proces
	inline int AddCell_Of(uint32_t cell, int bit) {
		nmiss--;
		known_b3 |= bit;
		return 1;
	}
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
	//===================== process not critical
	void ShrinkUasOfB3();
	void Go_miss1_b3();
	void Go_miss2_b3();
	void Go_miss3_b3();

};


struct GCHK {
	//=========== kwon puzzle filters
	BF128 puzknown,puzknown_perm;
	uint32_t kpfilt[4]; //initial statu 3 
	int kn_ir1, kn_ir2;
	int aigstop, aigstopxy, start_perm, *tpw, *tsortw,
		a_18_seen;
	uint32_t band_order[3];

	int diagbugchk;
	uint64_t debugvalbf;
	//___________________ studied solution 
	char * ze;// given solution grid
	char ze_diag[164];
	char * zp;// first 18 if any
	char zsol[82];// solution grid morphed
	int grid0[81];// same 0 based
	int grid0_diag[81];// diagonal symmetry
	STD_B416 * bands_abc[3],bA,bB;

	//____________structs hosting the  GUA2/3 entries
	struct SG2 {// 27 possible GUA2  
		int col1, col2;// columns  
		int i_27; //index 0 26 of the pair
		int i9, bit9;// locate the minirox index and bit
		int ip1, ip2;// cells
		int dig1, dig2;// digits
		uint32_t pat,digs;// pattern digits
		void Dump() {
			cout<<i_27<< " " << Char27out(pat)<< " ";
			cout <<dig1+1<<dig2+1<<" "<< Char9out(digs) << " ";
			cout << "\tcols " << col1 + 1 << col2 + 1 << " ";
			cout << "\tip 0-26  " << ip1 << " " << ip2 << endl;;

		}
	}tg2[27];
	int GetI27(int bf) {
		for (int i = 0; i < 27; i++)
			if (tg2[i].pat == bf) return i;
		return 0;// should never be
	}
	void Sg2Setup() {//load permanent data
		for (int i = 0; i < 27; i++) {// initial socket 2
			SG2 & w = tg2[i];		w.i_27 = i;
			register int row = i / 9, drow = 9 * row;
			int tpcol[3][2] = { {1,2},{0,2},{0,1} };
			int col = (i % 9), dcol = 3 * (col / 3), rdcol = col % 3, *p = tpcol[rdcol];
			w.i9 = col;
			w.col1 = dcol + p[0];		w.col2 = dcol + p[1];
			w.ip1 = drow + w.col1;		w.ip2 = drow + w.col2;
			w.bit9 = 1 << (i / 3);
			w.pat = (1 << w.ip1) | (1 << w.ip2);
			w.dig1 = grid0[w.ip1+54];	
			w.dig2 = grid0[w.ip2+54];
			w.digs = (1 << w.dig1) | (1 << w.dig2);
		}
	}

	//____ new  ua collector at start
	int StartIs18();
	void UaCollector();
	void Adduab12(int print = 0);
	void FirstUasCollect();
	void SecondUasCollect();
	void UasCollect4box();
	void Guas2Collect();
	void PutUasStartInVector();
	//_____ bands 1+2 direct expansion
	uint64_t v12_v0[64], v12_c[54][64]; //maxi 4096 uas

	struct V12_64 {
		uint64_t v;// current 64 bit vector
		uint32_t ind,// vector ind
		 u_start;//64*ind first ua of the vector
	};
	uint32_t nv12_64_spot[16]; // maxi 64 spots in expand
	V12_64 tv12_64[16][64]; // max 16 spot and 4096 uas

	//___________  first 3/4 then expand
	struct CPT_4 {
		uint64_t t[11]; // box band stack
		CPT_4();// init t
		int GetCount(uint64_t bf) {
			int n = 0;
			for (int i = 0; i < 11; i++)
				if (bf & t[i])n++;
			return n;
		}

	}cpt_4c;
	uint64_t  v12_4_c[54],nt4_to_expand;  
	struct T4_TO_EXPAND {
		uint64_t bf, active, vsort;
		void dump() {
			cout << Char2Xout(bf) << " ";
			cout << Char2Xout(active) << " ";
			cout << (vsort & 31) << " ";
			cout << (vsort >> 32) << endl;;

		}
	}t4_to_expand[5000];
	uint64_t tua4[2048];
	uint32_t ntua4;

	void BuildVectorsForExpand4B12();//64 uas
	void Expand4B12();
	void Do_phase2(T4_TO_EXPAND w);
	void Do_phase2Expand(uint64_t bf, uint64_t ac);

	//_____________ validb12
	uint64_t myb12, myac, myb12add;
	uint32_t  mynclues;//valid status
	uint32_t tadd[50], ntadd;
	struct VB12
	{	BF128 tmore[384], t2[128],tof128[50];
		uint64_t ort2, orof;
		uint32_t ntmore, nt2;
		uint32_t tg2ok[27], ntg2ok,
			tmore27[128],ntmore27;
		uint32_t tclues[15], nclues, bfbf2; //assign compulsory bf2
		uint32_t tof[50], ntof; // outfield

		MINCOUNT smin;
		int Getsmin();
		inline void ApplyBf2();// forcing common cell 2 pairs
		inline void BuildOf();
		inline void BuildOrOf(uint64_t ac);
		inline uint32_t GetAnd() {// called with ntof>0
			register uint32_t andw = tof[0];
			for(uint32_t i=1;i<ntof;i++)andw&= tof[i];
			return andw;
		}
		inline uint32_t GetAndExcept(uint32_t f) {
			register uint32_t andw = BIT_SET_27;
			for (uint32_t i = 0; i < ntof; i++) {
				register uint32_t U = tof[i];
				if(!(f&U))andw &= U;
			}
			return andw;
		}
		inline uint32_t GetMinOf() {
			uint32_t  uamin = tof[0];
			{
				register uint32_t min = _popcnt32(uamin);
				for (uint32_t i = 1; i < ntof; i++) {
					register uint32_t Ru = tof[i], cc = _popcnt32(Ru);
					if (cc < min) { min = cc;	uamin = Ru; }
				}
			}
			return uamin;
		}
		void CleanTmore();

		inline void GetOld(VB12 & vbo, uint64_t bf54);

		void Dumpt2() {
			ort2 = 0;
			for (uint32_t i = 0; i < nt2; i++) {
				cout << Char54out(t2[i].bf.u64[0]) 
					<< " " << t2[i].bf.u32[2] << endl;
				ort2 |= t2[i].bf.u64[0];
			}
			cout << Char54out(ort2) << " or" <<endl;
		}
		void Dumptmore() {
			for (uint32_t i = 0; i < ntmore; i++) {
				cout << Char54out(tmore[i].bf.u64[0]) << " ";
				cout << Char27out(tmore[i].bf.u32[2]) << endl;
			}
		}
		void Dumptof128() {
			for (uint32_t i = 0; i < ntof; i++) {
				cout << Char54out(tof128[i].bf.u64[0]) << " ";
				cout << Char27out(tof128[i].bf.u32[2]) << endl;
			}
		}

	}svb12,svb12add;
	// bands 1+2 valid epansion
	struct VADD {//max 64 c2; 384 cmroe
		uint64_t vc2, vcmore[2];
		inline void Apply(VADD vo, VADD vc) {
			vc2 = vo.vc2 & vc.vc2;
			vcmore[0] = vo.vcmore[0] & vc.vcmore[0];
			vcmore[1] = vo.vcmore[1] & vc.vcmore[1];
		}
	};
	struct VADD_HANDLER {
		VADD vadd0, vaddsteps[10], vaddcell[50];
		uint32_t mapcell[54];
		void SetUpAdd0(VB12& vb12, uint64_t ac);
		inline void Apply(uint64_t stepold, uint64_t icur) {
			vaddsteps[stepold + 1].Apply(vaddsteps[stepold],
				vaddcell[icur]);
		}

	}vaddh;


	void AfterExpandB12(uint64_t bf, uint64_t ac, int ncl);
	
	void ExpandAddB1B2();
	void ExpandAddB1B2Go(int step);
	void InitGoB3(uint64_t bf, uint64_t ac, int ncl);

	void GoB3(  int ncl, VB12 & vbx);
	void BuildExpandB3Vect( uint32_t cl0bf, uint32_t active0, VB12 & vbx);
	void ExpandB3Vect( uint32_t cl0bf=0,
		uint32_t active0=BIT_SET_27);


	uint32_t tclues[40], *tcluesxy;// mini 25+band a
	int nclues_step, nclues,nclf,nmiss;
	uint64_t n_to_clean, n_to_clean2,nwc;
	int  ncluesb3,mincluesb3;


	uint32_t nua_add, nua_3x3y, nxy_filt1,nvalid,n128add;

	//================== clean process
	uint64_t wb12bf, wb12active, myua;
	GINT64  stack_count_step, stack_count, stack_countf;

	uint32_t	ua_of_seen;


	//________ clean and valid
	uint32_t uasb3_1[10000], uasb3_2[2000], 
		nuasb3_1, nuasb3_2,  b3_andout;
	uint32_t nb64_1,nb64_2;

	MORE32 moreuas_b3;
	uint32_t ua_out_seen,clean_valid_done;
	//==================== current band 3 to process
	uint32_t tcluesb3[10],ntcl3,ntcl3_bf2,*tclues3;


	G17B3HANDLER hh0;

	//____ start the process for a triplet banda band b bandc 
	void Start(STD_B416 * bandx, int * tsort,int ip);

	//============ vectors 64 128 bits 
	struct UB2 {// for 2560 uas over 128 in step b1b2
		BF128 vx[20], vcx[20][54];
		void Init() {
			memset(vx, 0, sizeof vx);
			memset(vcx, 255, sizeof vcx);
		}
		inline void Add(uint64_t & ua, uint32_t i) {
			uint32_t ibloc = i >> 7, ir = i & 127;
			vx[ibloc].Set(ir);
			uint32_t cc64;// build cells vectors 
			register uint64_t Rw = ua;
			BF128 * vc = vcx[ibloc];
			while (bitscanforward64(cc64, Rw)) {
				Rw ^= (uint64_t)1 << cc64;// clear bit
				if (cc64 > 26)cc64 -= 5;
				vc[cc64].clearBit(ir);
			}
		}
		void ApplyStep(uint64_t & bf, uint32_t nuas) {
			// put bf in table
			uint32_t tcells[54], ntcells = 0, cc64;
			register uint64_t Rw = bf;
			while (bitscanforward64(cc64, Rw)) {
				Rw ^= (uint64_t)1 << cc64;// clear bit
				if (cc64 > 26)cc64 -= 5;
				tcells[ntcells++] = cc64;
			}
			uint32_t ibloc = 0;
			while (1) {
				BF128 w = vx[ibloc], *vc = vcx[ibloc];
				for (uint32_t i = 0; i < ntcells; i++)
					w &= vc[tcells[i]];
				vx[ibloc] = w;
				if (nuas > 128) { nuas -= 128; ibloc++; }
				else break;
			}
		}
		inline int ApplyXY(uint32_t *tcells, uint32_t ntcells, uint32_t nuas) {
			uint32_t ibloc = 0;
			while (1) {
				BF128 w = vx[ibloc], *vc = vcx[ibloc];
				for (uint32_t i = 0; i < ntcells; i++)
					w &= vc[tcells[i]];
				if (w.isNotEmpty()) return 1;
				if (nuas > 128) { nuas -= 128; ibloc++; }
				else break;
			}
			return 0;
		}
	}ub2;


	//_______ processing potential valid bands 1+2


	void Out17(uint32_t bfb3);

	void Debugifof();
};


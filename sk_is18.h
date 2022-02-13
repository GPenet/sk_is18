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
	uint64_t tua[4096], t2[3000];
	// 4096 size of vectors used in expansion
	uint32_t nua, nt2;
	void SwitchTo54Mode() {
		for (uint32_t i = 0; i < nua; i++) {
			register uint64_t R = tua[i];
			R = (R & BIT_SET_27) | ((R& BIT_SET_B2) >> 5);
			tua[i]=R;// now r54
		}
	}
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
		//if (1 && u27 == 1) {
		//	cout <<Char54out(u54)<< " add nt=" << nt << endl;
		//}
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

	void Get2(BF128 * td, uint32_t&  ntd,
		uint32_t * tc, uint32_t ntc) {
		register uint64_t V = v0;
		uint32_t ir;
		for (uint32_t i = 0; i < ntc; i++)
			V &= vc[tc[i]];
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			BF128 & w = td[ntd++];
			w.bf.u64[0] = tu54[ir] ;
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
	CHUNK3B c2[100],  cmore[100];
	CHUNK1B band3[2], misc[2];// up to 81 band3 up to 128 others
	uint32_t ic2,  icmore, iband3, imisc;
	void Init() {
		ic2=  icmore= iband3= imisc=0;
		c2[0].Init();  cmore[0].Init();
		band3[0].Init(); band3[1].Init();
		misc[0].Init();
	}
	inline int GetC2Count() { return 64 * ic2 +(int) c2[ic2].nt; }
	inline void Add128(BF128 w54, uint32_t cc) {
		if (cc == 2) {// add to c2 switch to index0-26
			if (ic2 == 99 && c2[99].nt > 63) return;
			if (c2[ic2].nt > 63)c2[++ic2].Init();
			c2[ic2].Add(w54.bf.u64[0], w54.bf.u32[2]);
		}

		else {// add to cmore
			if (icmore == 99 && cmore[99].nt > 63) return;
			if (cmore[icmore].nt > 63)
				cmore[++icmore].Init();
			cmore[icmore].Add(w54.bf.u64[0], w54.bf.u32[2]);

		}
	}


	void Addc2(uint64_t ua12, uint32_t ua) {
		uint64_t ua54 = (ua12 & BIT_SET_27) |
			((ua12 & BIT_SET_B2) >> 5);
		if (ic2 == 99 && c2[99].nt > 63) return;
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

	void GetMoreClean(BF128 * td, uint32_t&  ntd,
		uint32_t * tc, uint32_t ntc ) {
		ntd = 0;
		for (uint32_t i = 0; i <= icmore; i++)
			if (ntd < 1024)// table risk of overflow
			cmore[i].Get2(td, ntd, tc, ntc);
	}

	void Get2Clean(BF128 * td, uint32_t&  ntd,
		 uint32_t * tc, uint32_t ntc) {
		ntd = 0;
		for (uint32_t i = 0; i <= ic2; i++)
			if(ntd<1024)
			 c2[i].Get2(td, ntd, tc, ntc);
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

	inline void Add54(uint64_t v) {//add a new more if <128 
		if (nt < 128) {
			uint32_t cc54;// build cells vectors 
			vect.Set(nt);
			register uint64_t Rw = v;
			while (bitscanforward64(cc54, Rw)) {
				Rw ^= (uint64_t)1 << cc54;// clear bit
				vc[cc54].clearBit(nt);
			}
			t[nt++] = v;
		}
	}


	inline int ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		if (!nt) return 0;
		BF128 w = vect;
		for (uint32_t i = 0; i < ntcells; i++)
			w &= vc[tcells[i]];
		return (w.isNotEmpty());
	}
	inline void Extract(uint32_t *tc, uint32_t ntc,
		uint64_t *td, uint64_t & ntd) {
		if (!nt) return;
		BF128 w = vect;
		uint32_t ir;
		for (uint32_t i = 0; i < ntc; i++)	w &= vc[tc[i]];
		{
			register uint64_t R = w.bf.u64[0];
			while (bitscanforward64(ir, R)) {
				R ^= (uint64_t)1 << ir;
				td[ntd++] = t[ir];
			}
		}
		{
			register uint64_t R = w.bf.u64[1];
			while (bitscanforward64(ir, R)) {
				R ^= (uint64_t)1 << ir;
				td[ntd++] = t[ir+64];
			}
		}
	}

};
struct MOREV2 {// 2 more64vect paired
	MORE64VECT mv1, mv2;
	uint32_t sw12;
	inline void Init() { mv1.Init(); mv2.Init(); sw12 = 0; }

	inline void Add54(uint64_t v) {//add a new more if <128 
		if (sw12 && (mv2.nt == 128)) {
			mv1.Init();
			sw12 = 0;// mv1 active
		}
		if ((!sw12) && (mv1.nt == 128)) {
			mv2.Init();
			sw12 = 1;// mv2 active
		}
		if (sw12)mv2.Add54(v);
		else mv1.Add54(v);
	}	
	inline int ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		if (mv1.ApplyXY(tcells, ntcells)) return 1;
		return mv2.ApplyXY(tcells, ntcells);
	}
	inline void Extract(uint32_t *tc, uint32_t ntc,
		uint64_t *td, uint64_t & ntd) {
		mv1.Extract(tc, ntc, td, ntd);
		mv2.Extract(tc, ntc, td, ntd);
	}

	void Status(const char * lib) {
		cout << "status for more uas " << lib << endl;
		cout << mv1.nt << " " << mv2.nt << " " << mv1.nt + mv2.nt << endl;
	}


}morev2a,  morev2c, morev2d;

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
#define BUFVALIDS 8000
struct GCHK {
	//=========== kwon puzzle filters
	BF128 puzknown;
	uint64_t pk54;
	uint32_t kpfilt[4]; //initial statu 3 
	int kn_ir1, kn_ir2;
	int aigstop, aigstopxy, start_perm, *tpw, *tsortw,
		a_18_seen,minb1b2, minb1, minb2;
	uint32_t band_order[3];

	int diagbugchk;
	uint64_t debugvalbf;
	//___________________ studied solution 
	char * ze;// given solution grid
	char zes[164],ze_diag[164];
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
	inline void Adduab12();
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
	V12_64 *tv12_64f[16];
	V12_64 tv12_64[16][64]; // max 16 spot and 4096 uas

	uint32_t nv12_64_spot[16]; // maxi 64 spots in expand

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
	uint64_t bufvalid[BUFVALIDS+1], *pbufvalid,*pendbufvalid;
	int nt4ok,okcheck;// for known
	// bands 1+2 valid epansion
	struct CHUNKVB12
	{
		BF128 tmore[1024 + 64], t2[1024 + 64];
		//+64 margin against overflow
		uint32_t ntmore, nt2;
		void SortTmore() {// sort by size
			BF128 tt[5][400];
			uint32_t ntt[5];
			memset(ntt, 0, sizeof ntt);
			for (uint64_t i = 0; i < ntmore; i++) {
				BF128 w = tmore[i];
				uint32_t cc = _popcnt32(w.bf.u32[2]) - 3;
				if (cc > 4)cc = 4;
				tt[cc][ntt[cc]++] = w;
			}
			ntmore = 0;
			for (int i = 0; i < 5; i++) //3,4;5;6; more
				for (uint32_t j = 0; j < ntt[i]; j++)
					tmore[ntmore++] = tt[i][j];

		}
		void Dumpt2(int ix = -1) {
			cout << "dumpt2 nt2=" << nt2 << endl;
			for (uint32_t i = 0; i < nt2; i++) {
				if (ix >= 0 && (t2[i].bf.u32[2] != ix)) continue;
				cout << Char54out(t2[i].bf.u64[0])
					<< " " << t2[i].bf.u32[2] << endl;
			}
		}
		void Dumptmore() {
			cout << "dumptmore ntmore=" << ntmore << endl;
			for (uint32_t i = 0; i < ntmore; i++) {
				cout << Char54out(tmore[i].bf.u64[0]) << " ";
				cout << Char27out(tmore[i].bf.u32[2]) << endl;
			}
		}
	}chvb12;
	struct CHUNKVB12ADD
	{
		BF128 t[128], t2a[64], tma[128];
		uint32_t nt, nt2, nt2a, ntma;
		void Shrink(CHUNKVB12 & o, uint64_t bf, uint64_t ac) {// apply band 1+2 to old
			nt = 0;
			register uint64_t F = bf;
			for (uint64_t i = 0; i < o.nt2; i++) {
				BF128 w = o.t2[i];
				if (!(w.bf.u64[0] & F)) {
					w.bf.u64[0] &= ac;
					t[nt++] = w;
				}
			}
			nt2 = nt;
			// keep the best first 64 (smallest band3)
			BF128 tt[5][100];
			uint32_t ntt[5];
			memset(ntt, 0, sizeof ntt);
			for (uint64_t i = 0; i < o.ntmore; i++) {
				BF128 w = o.tmore[i];
				if (!(w.bf.u64[0] & F)) {
					w.bf.u64[0] &= ac;
					uint32_t cc = _popcnt32(w.bf.u32[2]) - 3;
					if (cc > 4) cc = 4;
					tt[cc][ntt[cc]++] = w;	//tmore[ntmore++] = w;
				}
			}
			for (int i = 0; i < 5; i++) //3,4;5;6
				for (uint32_t j = 0; j < ntt[i]; j++)
					if (nt < 128) t[nt++] = tt[i][j];
		}
		void Dumpt2(int ix = -1) {
			cout << "dumpt2 nt2=" << nt2 << endl;
			for (uint32_t i = 0; i < nt2; i++) {
				if (ix >= 0 && (t[i].bf.u32[2] != ix)) continue;
				cout << Char54out(t[i].bf.u64[0])
					<< " " << t[i].bf.u32[2] << endl;
			}
		}
		void Dumptmore() {
			cout << "dumptmore ntmore=" << nt - nt2 << endl;
			for (uint32_t i = nt2; i < nt; i++) {
				cout << Char54out(t[i].bf.u64[0]) << " ";
				cout << Char27out(t[i].bf.u32[2]) << endl;
			}
		}
	}  chvb12b;
	struct VB12
	{
		BF128 tof128[50], *myt2;// , *mytmore;
		uint32_t mynt2;// , myntmore;
		//+64 margin against overflow
		uint64_t ort2, orof;
		uint32_t tg2ok[27], ntg2ok,
			tmore27[384], ntmore27;// max is from tmore
		uint32_t tclues[15], nclues, bfbf2; //assign compulsory bf2
		uint32_t tof[50], ntof; // outfield

		MINCOUNT smin;
		int GetsminF(BF128 * t2, uint32_t nt2);
		void CleanTmore(BF128 * tmore, uint32_t ntmore);
		inline void ApplyBf2();// forcing common cell 2 pairs
		inline void BuildOf();
		inline void BuildOrOf(CHUNKVB12& o, uint64_t ac);
		inline uint32_t GetAnd() {// called with ntof>0
			register uint32_t andw = tof[0];
			for (uint32_t i = 1; i < ntof; i++)andw &= tof[i];
			return andw;
		}
		inline uint32_t GetAndExcept(uint32_t f) {
			register uint32_t andw = BIT_SET_27;
			for (uint32_t i = 0; i < ntof; i++) {
				register uint32_t U = tof[i];
				if (!(f&U))andw &= U;
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


		void Dumptof128() {
			for (uint32_t i = 0; i < ntof; i++) {
				cout << Char54out(tof128[i].bf.u64[0]) << " ";
				cout << Char27out(tof128[i].bf.u32[2]) << endl;
			}
		}
		void Dumptmore27() {
			cout << "tmore27" << endl;
			for (uint32_t i = 0; i < ntmore27; i++) {
				cout << Char27out(tmore27[i]) << " " << i << endl;
			}
		}
		void Dumptof() {
			cout << "tof" << endl;
			for (uint32_t i = 0; i < ntof; i++) {
				cout << Char27out(tof[i]) << endl;
			}
		}



	}svb12, svb12b, svb12add;
	struct GADDB3 { //Add during a clean chunk
		BF128 tmore[384], t2[256];
		uint32_t ntmore, nt2;
		inline void Init() { ntmore = nt2 = 0; }
		inline void Add(BF128 w, int cc) {
			if (cc < 3) {
				if (nt2 < 256)t2[nt2++] = w;
			}
			else if (ntmore < 384)
				tmore[ntmore++] = w;
		}

		void Dumpt2() {
			cout << "dumpt2 nt2=" << nt2 << endl;
			for (uint32_t i = 0; i < nt2; i++) {
				cout << Char54out(t2[i].bf.u64[0])
					<< " " << t2[i].bf.u32[2] << endl;
			}
		}
		void Dumptmore() {
			cout << "dumptmore ntmore=" << ntmore << endl;
			for (uint32_t i = 0; i < ntmore; i++) {
				cout << Char54out(tmore[i].bf.u64[0]) << " ";
				cout << Char27out(tmore[i].bf.u32[2]) << endl;
			}
		}

	}gaddb3;
	struct VADD_HANDLER {
		BF128 v0, v2, vm, vsteps[10], vcells[50];
		uint32_t mapcell[54];
		void SetUpAdd0(CHUNKVB12ADD& vb12, uint64_t ac);
		inline void Apply(uint64_t stepold, uint64_t icur) {
			vsteps[stepold + 1] = vsteps[stepold] & vcells[icur];
		}
		void Dump();

	}vaddh;
	struct MOREAND {
		uint64_t tm[500], ntm;
		inline int Check(uint64_t bf) {
			register uint64_t F = bf, i;
			for (i = 0; i < ntm; i++)
				if (!(F&tm[i]))return 1;
			return 0;
		}
		inline uint64_t GetAndUnHit(uint64_t bf) {
			register uint64_t andw = ~0;
			register uint64_t F = bf, i;
			for (i = 0; i < ntm; i++)
				if (!(F&tm[i]))andw &= tm[i];
			return andw;
		}
		void Dump() {
			for (uint64_t i = 0; i < ntm; i++)
				cout << Char54out(tm[i]) << " " << i << " "
				<< _popcnt64(tm[i]) << endl;
		}

	}moreand;
	struct MOREANDB {
		uint64_t tm[100], ntm;
		inline void GetdUnHit(uint64_t bf, uint64_t *tmo, uint64_t ntmo) {
			ntm = 0;
			register uint64_t F = bf, i;
			for (i = 0; i < ntmo; i++)
				if (!(F&tmo[i]))tm[ntm++] = tmo[i];
		}
		void Dump() {
			for (uint64_t i = 0; i < ntm; i++)
				cout << Char54out(tm[i]) << endl;
		}

	}mbisvalid;

	void BuildVectorsForExpand4B12();//64 uas
	void Expand4B12();
	void Do_phase2(T4_TO_EXPAND w);
	void Do_phase2Expand(uint64_t bf, uint64_t ac);
	int IsValidB12();
	void CheckValidBelow(uint64_t bf,uint64_t ac);
	void CleanBufferAndGo(uint64_t andvalid, uint64_t orvalid);
	void CleanMoreUas(uint64_t bf, uint64_t ac, int ncl, MOREANDB & mabo);
	//_____________ validb12
	uint64_t myb12, myb12f, myac, myacf, myb12add;
	uint32_t  mynclues;//valid status
	int 	limb12;
	uint32_t tadd[50], ntadd;


	void ExpandAddB1B2(uint64_t bf);
	void ExpandAddB1B2Go(int step);
	void InitGoB3(uint64_t bf, uint64_t ac54);
	void InitGoB3below(uint64_t bf, uint64_t ac54);

	void GoB3(  int ncl,  VB12 & vbx);
	void BuildExpandB3Vect( uint32_t cl0bf, uint32_t active0,
		 VB12 & vbx);
	void ExpandB3Vect( uint32_t cl0bf=0,
		uint32_t active0=BIT_SET_27);
	inline void NewGua(BF128 w);

	uint32_t tclues[40], *tcluesxy;// mini 25+band a
	int nclues_step, nclues,nclf,nmiss;
	int  ncluesb3,mincluesb3;

	//================== clean process
	uint64_t wb12bf, wb12active, myua;
	uint32_t taddgob3[100], ntaddgob3, clean_valid_done;
	MORE32 moreuas_b3;

	void Out17(uint32_t bfb3);

};


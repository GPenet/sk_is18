struct CBS {// clues band stack
	uint32_t b[3];
	uint32_t s[3];
	inline void Add(uint32_t cell) {
		b[cell / 27]++;
		s[C_stack[cell]]++;
	}
	inline int IsFilt11() {
		if (b[0] > 7 || b[1] > 6)return 1;
		if(s[0] > 7 || s[1] > 7 || s[2] > 7)return 1;
		return 0;
	}
	inline int IsFilt12() {
		if (b[0] != 6)return 1;
		if (s[0] > 6 || s[1] > 6 || s[2] >6)return 1;
		return 0;
	}
}cbs_4,cbs_5;
struct SPB03;
struct TUAS81 {
	BF128 tall[5000],
		told[3000], // later  has band 3
		t2_2[3000];//band 3 2/3 after 5 clues 
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
	void Debug(BF128 *t, uint32_t nt, const char * lib) {
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

}tuas81;
struct TUASB12 {
	uint64_t tua[2560], t2[3000];
	uint32_t nua, nt2;
	// tua is initial set of uas bands 1+2 plus local adds
	//t2 is for guas generation
	struct TUAVECT {// vector 128 for tua
		BF128 v0, vc[54];
		uint64_t t[128];
		void Init() {
			v0.SetAll_0();
			memset(vc, 255, sizeof vc);
		}
		void Dumpv(int lim = 32) {
			cout << Char32out(v0.bf.u32[0]) << " v0" << endl;
			for(int i=0;i<54;i++ )
				cout << Char32out(vc[i].bf.u32[0]) << " cell="<<i << endl;
		}
		void Dump(int lim = 128) {
			for (int i = 0; i < lim; i++)
				if (v0.On(i))
					cout << Char54out(t[i])<<" " << i << endl;
				else return;
		}
	};
	struct TUAVECTH {
		TUAVECT tv128[20];// max start 20*128=2560 uas 
		uint32_t na128, nablocs, nta128[20];
		TUAVECT tvb128[10];// max start 10*128=1280 uas 
		uint32_t nb128, nbblocs, ntb128[20];
		TUAVECT tvc128[10];// max start 10*128=1280 uas 
		uint32_t nc128, ncblocs, ntc128[20];
		void Build_tv128();
		void Build_tvb128(uint32_t * ntt);
		void Build_tvc128(uint32_t* ntt);
		inline void AddA(uint64_t u) {
			if (na128 >= 2560) return;
			register uint32_t bloc = na128 >> 7, ir = na128 - 128 * bloc;
			na128++; nablocs = bloc; nta128[bloc]++;
			tv128[bloc].v0.setBit(ir);
			BF128* myvc = tv128[bloc].vc;
			register uint64_t R = u;
			tv128[bloc].t[ir] = R;
			uint32_t cell;
			while (bitscanforward64(cell, R)) {
				R ^= (uint64_t)1 << cell; //clear bit
				myvc[cell].clearBit(ir);
			}
		}
		inline void AddB(uint64_t u) {
			if (nb128 >= 1280) return;
			register uint32_t bloc = nb128 >> 7,
				ir = nb128 - 128 * bloc;
			nb128++; nbblocs = bloc; ntb128[bloc]++;
			tvb128[bloc].v0.setBit(ir);
			BF128* myvc = tvb128[bloc].vc;
			register uint64_t R = u;
			tvb128[bloc].t[ir] = R;
			register uint32_t cell;
			while (bitscanforward64(cell, R)) {
				R ^= (uint64_t)1 << cell; //clear bit
				myvc[cell].clearBit(ir);
			}
		}
		inline void AddC(uint64_t u) {
			if (nc128 >= 1280) return;
			register uint32_t bloc = nc128 >> 7,
				ir = nc128 - 128 * bloc;
			nc128++; ncblocs = bloc; ntc128[bloc]++;
			tvc128[bloc].v0.setBit(ir);
			BF128* myvc = tvc128[bloc].vc;
			register uint64_t R = u;
			tvc128[bloc].t[ir] = R;
			register uint32_t cell;
			while (bitscanforward64(cell, R)) {
				R ^= (uint64_t)1 << cell; //clear bit
				myvc[cell].clearBit(ir);
			}
		}
		void Status(int mode = 0) {
			cout << "nax status " << na128
				<< " " << nb128 << " " << nc128 << endl;
		}
	}tv128h;
	void SwitchTo54Mode() {
		for (uint32_t i = 0; i < nua; i++) {
			register uint64_t R = tua[i];
			R = (R & BIT_SET_27) | ((R& BIT_SET_B2) >> 5);
			tua[i]=R;// now r54
		}
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
	void Stats() {
		cout << "tuasb12 status end of perm"
			<< " nua=" << nua 	<< endl;
	}
}tuasb12;

//___________ all uas in bands 1+2 and band 3
struct GUA {
	uint64_t tu54[64], killer;
	uint32_t nt, uab3, index, igua2_gua3;
	void Init(uint32_t i, uint32_t u3, uint32_t i23, uint64_t ua0) {
		nt = 1; index = i; igua2_gua3 = i23; uab3 = u3;
		killer = tu54[0]=ua0;
	}
	void ReInit() {
		nt = 1; 	killer = tu54[0] = 0;
	}
	inline void Add( uint64_t ua0) {
		if (nt < 64) {	killer &= ua0;	tu54[nt++] = ua0;	}
	}
	void Dump() {
		cout << Char27out(uab3) << "pat nt="<<nt 
			<<" i9/27="<< igua2_gua3 << endl;
		cout << Char54out(killer) << " killer" << endl;
		if(nt>1)	for (uint32_t i = 0; i < nt; i++)
			cout << Char54out(tu54[i]) << " " << i << endl;
	}
	void DumpShort() {
		cout << Char27out(uab3) << "pat nt=" << nt<< endl;
	}
};

struct G2_256 {//vector 256 chunk level gua2 gua3
	BF128 v[2], vc[2][54], vp[2][27];
	uint32_t bits[256];//bit  i27 or i9
	uint32_t pat27[256];// i27 or i9

	uint32_t nv;
	void Init() {
		nv = 0;
		memset(v, 0, sizeof v);
		memset(vc, 255, sizeof(vc));
		memset(vp, 255, sizeof(vp));
	}

	void Add(BF128 gx) { // gx is a still valid gua2 gua3 
		if (nv >= 256)return; // should never be
		int ix = 0, nx = nv;
		if (nv >= 128) {	ix = 1; nx -= 128;	}
		BF128* vcx = vc[ix];
		v[ix].Set(nx);
		bits[nv] = gx.bf.u32[2];
		register uint64_t V = gx.bf.u64[0];
		register uint32_t x;
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;// clear bit
			vcx[x].clearBit(nx);
		}
		register uint32_t i27 = gx.bf.u32[3];
		pat27[nv]=i27;
		if(i27<27 )vp[ix][i27].clearBit(nx);
		nv++;
	}
	void GetActive0(BF128 vw, uint32_t& gbf,
		uint32_t* tg, uint32_t& ntg) {
		for (register uint32_t i = 0; i < ntg; i++) 
			vw &= vp[0][tg[i]];//apply olds
		register uint64_t V = vw.bf.u64[0];
		register uint32_t  x;
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;
			gbf |= bits[x];
			register uint32_t i27 = pat27[x];
			tg[ntg++] = i27;
			vw &= vp[0][i27];
			V &= vp[0][i27].bf.u64[0];
		}
		V = vw.bf.u64[1];
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;
			gbf |= bits[x+64];
			register uint32_t i27 = pat27[x+64];
			tg[ntg++] = i27;
			V &= vp[0][i27].bf.u64[1];
		}
	}
	void GetActive1(BF128 vw, uint32_t& gbf,
		uint32_t* tg, uint32_t& ntg) {
		for (register uint32_t i = 0; i < ntg; i++)
			vw &= vp[1][tg[i]];//apply olds
		if (vw.isEmpty()) return;
		register uint64_t V = vw.bf.u64[0];
		register uint32_t  x;
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;
			gbf |= bits[x+128];
			register uint32_t i27 = pat27[x+128];
			tg[ntg++] = i27;
			vw &= vp[1][i27];
			V &= vp[1][i27].bf.u64[0];
		}
		V = vw.bf.u64[1];
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;
			gbf |= bits[x + 192];
			register uint32_t i27 = pat27[x + 192];
			tg[ntg++] = i27;
			V &= vp[1][i27].bf.u64[1];
		}
	}
	uint32_t Apply(uint32_t * tc, uint32_t ntc) {
		if (!nv) return 0;
		BF128 vw[2];
		memcpy(vw, v, sizeof v);

		uint32_t g2 = 0;
		for (uint32_t i = 0; i < ntc; i++) {// apply
			uint32_t cell = tc[i];
			vw[0] &= vc[0][cell];
			vw[1] &= vc[1][cell];
		}
		for (uint32_t i = 0; i < 4; i++) {// extract
			register uint64_t V = vw[0].bf.u64[i];
			register uint32_t * b = &bits[64 * i], x;
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				g2 |= b[x];
			}
		}
		return g2;
	}
	void ApplyMore(uint32_t * tc, uint32_t ntc,
		uint32_t * tmore, uint32_t & ntmore, uint32_t asb3) {
		if (!nv) return ;
		register uint32_t Rasb3 = asb3;
		BF128 vw[2];
		memcpy(vw, v, sizeof v);
		for (uint32_t i = 0; i < ntc; i++) {// apply
			uint32_t cell = tc[i];
			vw[0] &= vc[0][cell];
			vw[1] &= vc[1][cell];
		}
		register uint64_t* pV = vw[0].bf.u64;
		for (uint32_t i = 0; i < 4; i++) {// extract
			register uint64_t V = pV[i];
			register uint32_t * b = &bits[64 * i], x;
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint32_t U = b[x];
				if (!(Rasb3&U))	tmore[ntmore++] = U;
			}
		}
	}
	void Dump() {
		for (uint32_t i = 0; i < nv; i++) {
			cout << Char27out(bits[i]) << " ";
			if (i < 128) {
				for (int j = 0; j < 54; j++)
					if (vc[0][j].On(i)) cout << ".";
					else cout << "1";
			}
			else {
				for (int j = 0; j < 54; j++)
					if (vc[1][j].On(i-128)) cout << ".";
					else cout << "1";
			}
			cout<<" "<<i << endl;
		}
	}
	
}g2_256[2], g3_256,gm_256[2];
struct G2_64Add {
	uint64_t u[64];//ua
	uint32_t p[64];// i27 or i9 or pat
	uint32_t nv;
	void Add(BF128 w) {
		if (nv >= 64)return;
		u[nv] = w.bf.u64[0];
		p[nv++]=w.bf.u32[2];
	}
	void Getactive(uint64_t bf, uint32_t& gbf,
		uint32_t* tg, uint32_t& ntg) {
		for (uint32_t i = 0; i < nv; i++) {
			if (u[i] & bf) continue;
			register uint32_t P= p[i],bit=1<<P;
			if (gbf & bit) continue;
			gbf |= bit;
			tg[ntg++] = P;
		}
	}

}g2_64, g3_64, gm_64;
struct G2_256Handler {
	uint32_t n2,  nm;
	inline void GetStart(uint32_t n2e, uint32_t nme) {
		n2 = n2e;  nm = nme;
		BF128* V = vcl[0];
		memcpy(V, g2_256[0].v, sizeof g2_256[0].v);
		memcpy(&V[2], g2_256[1].v, sizeof g2_256[1].v);
		memcpy(&V[4], g3_256.v, sizeof g3_256.v);
		g2_64.nv=0; g3_64.nv=0; gm_64.nv=0;
	}

	inline void Add128(BF128 w, uint32_t cc,int add=0) {
		if (add) {
			if (cc == 2) g2_64.Add(w);
			else if (cc == 3) g2_64.Add(w);
			else gm_64.Add(w);
		}
		if (cc < 4) {
			w.bf.u32[2]=1<<w.bf.u32[2];// switch to bit field
			if (cc == 2) {
				//cout << "add 2 in g2_256 n2=" << n2 << endl;
				if (n2++ < 256)g2_256[0].Add(w);
				else g2_256[1].Add(w);
			}
			else { g3_256.Add(w); }
			return;
		}
		if (nm++ < 256)gm_256[0].Add(w);
		else gm_256[1].Add(w);


	}


	void ApplyMore(uint32_t * tc, uint32_t ntc,
		uint32_t * tmore, uint32_t & ntmore, uint32_t asb3) {
		gm_256[0].ApplyMore( tc, ntc,tmore,  ntmore,asb3);
		if(nm>256)
			gm_256[1].ApplyMore(tc, ntc, tmore, ntmore, asb3);
		if (gm_64.nv) {
			for (uint32_t i = 0; i < gm_64.nv; i++) {
				if (gm_64.u[i] & bf12) continue;
				register uint32_t P = gm_64.p[i];
				if (asb3 & P) continue;
				tmore[ntmore++] = P;
			}
		}
	}

	BF128 vcl[7][6];//6 clues 4 g2 2 g3
	uint64_t bf12;
	uint32_t g2, g3, tg2[27], tg3[9], ntg2 , ntg3 ;
	inline void NewVcl(int ncl, int clue) {
		register int i1 = ncl - 7;
		register BF128* V1 = vcl[i1], * V2 = vcl[i1 + 1];
		V2[0] = V1[0] & g2_256[0].vc[0][clue];
		V2[1] = V1[1] & g2_256[0].vc[1][clue];
		V2[2] = V1[2] & g2_256[1].vc[0][clue];
		V2[3] = V1[3] & g2_256[1].vc[1][clue];
		V2[4] = V1[4] & g3_256.vc[0][clue];
		V2[5] = V1[5] & g3_256.vc[1][clue];
	}
	void GetActive(int nclues, uint64_t bf) {
		BF128* V2 = vcl[nclues - 6];
		bf12 = bf;
		g2 = g3 = ntg2=ntg3=0;
		g2_256[0].GetActive0(V2[0], g2, tg2, ntg2);
		if(V2[1].isNotEmpty())
			g2_256[0].GetActive1(V2[1], g2, tg2, ntg2);
		if (V2[2].isNotEmpty())
			g2_256[1].GetActive0(V2[2], g2, tg2, ntg2);
		if (V2[3].isNotEmpty())
			g2_256[1].GetActive1(V2[3], g2, tg2, ntg2);
		if (V2[4].isNotEmpty())
			g3_256.GetActive0(V2[4], g3, tg3, ntg3);
		if (V2[5].isNotEmpty())
			g3_256.GetActive1(V2[5], g3, tg3, ntg3);
		if(g2_64.nv)g2_64.Getactive(bf, g2, tg2, ntg2);
		if (g3_64.nv)g3_64.Getactive(bf, g3, tg3, ntg3);
	}
	void Dumpvcl(int nclues) {
		cout << "Dump vcl nclues=" << nclues << endl;
		BF128* V2 = vcl[nclues - 6];
		uint64_t* v64 = V2[0].bf.u64;
		for (int i = 0; i < 12; i++) {
			register uint64_t v = v64[i];
			if (v)	cout << Char64out(v) << " i=" << i
				<< " cpt=" << _popcnt64(v) << endl;

		}
	}

}g_256h;
struct CHUNK3B {// storing 64 uas and vectors 3 bands
	uint64_t tu54[64], v0, 
		vc[54], //cell for tu54
		vclean, // v0 after clean group
		nt;
	uint32_t tu27[64];// index 0-26 or pattern
	inline void Init() {
		nt=0,		v0 = vclean=0;
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
	inline void ApplyClean(uint32_t * tc, uint32_t ntc) {
		register uint64_t V = v0;
		for (uint32_t i = 0; i < ntc; i++)
			V &= vc[tc[i]];
		vclean = V;
	}

	void Get2x(BF128 * td, uint32_t&  ntd,
		uint32_t * tc, uint32_t ntc) {
		register uint64_t V = v0;
		uint32_t ir;
		for (uint32_t i = 0; i < ntc; i++)
			V &= vc[tc[i]];
		vclean = V;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			BF128 & w = td[ntd++];
			w.bf.u64[0] = tu54[ir] ;
			w.bf.u32[2] = tu27[ir];
		}
	}

	void GetAddClean(uint64_t * td54, uint32_t * td27, uint32_t&  ntd) {
		ntd = 0;
		register uint64_t V = vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			td27[ntd] = tu27[ir];
			td54[ntd++] = tu54[ir];
		}
	}
	void DebugC2(int all=1) {
		for (uint32_t i = 0; i < nt; i++) {
			cout << Char54out(tu54[i]) << " " << tu27[i] << "\ti=" << i << endl;
		}
		if (!all) return;
		cout << Char64out(v0) << " v0" << endl;
		for(int i=0;i<54;i++) 
			if ((vc[i] & v0) != v0)
				cout << Char64out(vc[i] & v0) << " cell=" <<i<< endl;
	}
	void DebugC2Clean() {
		for (uint32_t i = 0; i < nt; i++) {
			uint64_t bit = (uint64_t)1 << i;
			if (vclean & bit) 
			cout << Char54out(tu54[i]) << " " << tu27[i] << " i=" << i<< endl;
		}
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
		
		{//__ check redundancy
			register uint32_t vn = ~v;
			for (uint64_t i = 0; i < nt; i++)
				if (!(tua[i] & vn))return; // == or subset
		}

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
struct CHUNKS_HANDLER {
#define CSIZE 99
	CHUNK3B c2[CSIZE + 1], c3[CSIZE + 1],
		c4[CSIZE + 1], c5[CSIZE + 1], cmore[CSIZE + 1];
	CHUNK1B band3[2];// up to 81 band3 up to 128 others
	uint32_t ic2, ic3, ic4, ic5, icmore, iband3;
	//___________ valid to process t4 5 6 more
	uint32_t t12[1100], nt12;
	void Init() {
		ic2 = ic3 = ic4 = ic5 = icmore = iband3 = 0;
		c2[0].Init();  c3[0].Init();
		c4[0].Init();  c5[0].Init(); cmore[0].Init();
		band3[0].Init(); band3[1].Init();
	}
	inline int GetC2Count() { return 64 * ic2 + (int)c2[ic2].nt; }

	//inline int  Check2(BF128 w54);
	//inline int  Check3(BF128 w54);
	//inline int  Check4(BF128 w54);
	inline void Add128(BF128 w54, uint32_t cc) {
		switch (cc) {
		case 2:
			//if(Check2(w54 ))return;
			if (ic2 == CSIZE && c2[CSIZE].nt > 63) return;
			if (c2[ic2].nt > 63)c2[++ic2].Init();
			c2[ic2].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		case 3:
			//if (Check3(w54))return;
			if (ic3 == CSIZE && c3[CSIZE].nt > 63) return;
			if (c3[ic3].nt > 63)c3[++ic3].Init();
			c3[ic3].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		case 4:
			//if (Check4(w54))return;
			if (ic4 == CSIZE && c4[CSIZE].nt > 63) return;
			if (c4[ic4].nt > 63)c4[++ic4].Init();
			c4[ic4].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		case 5:
			if (ic5 == CSIZE && c5[CSIZE].nt > 63) return;
			if (c5[ic5].nt > 63)c5[++ic5].Init();
			c5[ic5].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		default:
			if (icmore == CSIZE && cmore[CSIZE].nt > 63) return;
			if (cmore[icmore].nt > 63)cmore[++icmore].Init();
			cmore[icmore].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		}

	}

	void Addc2(uint64_t ua12, uint32_t i27) {
		uint64_t ua54 = (ua12 & BIT_SET_27) |
			((ua12 & BIT_SET_B2) >> 5);
		if (ic2 == CSIZE && c2[CSIZE].nt > 63) return;
		if (c2[ic2].nt > 63)c2[++ic2].Init();
		c2[ic2].Add(ua54, i27);

	}
	void Addband3(uint32_t w) {
		if (iband3 == 1 && band3[1].nt > 63) return;
		if (band3[iband3].nt > 63)band3[++iband3].Init();
		band3[iband3].Add(w);
	}

	inline void ApplyClean(uint32_t * tc, uint32_t ntc) {
		for (uint32_t i = 0; i <= ic2; i++)
			c2[i].ApplyClean(tc, ntc);
		for (uint32_t i = 0; i <= ic3; i++)
			c3[i].ApplyClean(tc, ntc);
		for (uint32_t i = 0; i <= ic4; i++)
			c4[i].ApplyClean(tc, ntc);
		for (uint32_t i = 0; i <= ic5; i++)
			c5[i].ApplyClean(tc, ntc);
		for (uint32_t i = 0; i <= icmore; i++)
			cmore[i].ApplyClean(tc, ntc);
	}

	void DebugAll(int full = 0) {
		cout << "chunkh debug all"
			<< "\nic2/nc2 " << ic2 << " " << c2[ic2].nt
			<< "\nic3/nc3 " << ic3 << " " << c3[ic3].nt
			<< "\nic4/nc4 " << ic4 << " " << c4[ic4].nt
			<< "\nic5/nc5 " << ic5 << " " << c5[ic5].nt
			<< "\nicm/ncm " << icmore << " " << cmore[icmore].nt << endl;
		if (!full) return;
		if (full == 2) {
			cout << "band 3" << endl;
			for (uint32_t i = 0; i <= iband3; i++)
				band3[i].Debug(0);
			cout << "guas 2 cells" << endl;
			for (uint32_t i = 0; i <= ic2; i++)
				c2[i].DebugC2(0);
		}
		cout << "guas 3 cells" << endl;
		for (uint32_t i = 0; i <= ic3; i++)
			c3[i].DebugC2(0);
		cout << "guas more cells" << endl;
		for (uint32_t i = 0; i <= ic4; i++)
			c4[i].DebugMore(0);
		for (uint32_t i = 0; i <= ic5; i++)
			c5[i].DebugMore(0);
		for (uint32_t i = 0; i <= icmore; i++)
			cmore[i].DebugMore(0);

	}
	void Debug4Clean(int full = 0) {
		cout << "chunkh debugclean all"
			<< "\nic2/nc2 " << ic2 << " " << c2[ic2].nt
			<< "\nic3/nc3 " << ic3 << " " << c3[ic3].nt
			<< "\nic4/nc4 " << ic4 << " " << c4[ic4].nt
			<< "\nic5/nc5 " << ic5 << " " << c5[ic5].nt
			<< "\nicm/ncm " << icmore << " " << cmore[icmore].nt << endl;
		cout << Char64out(c2[0].vclean) << " c2" << endl;
		cout << Char64out(c3[0].vclean) << " c3" << endl;
		cout << Char64out(c4[0].vclean) << " c4" << endl;
		cout << Char64out(c5[0].vclean) << " c5" << endl;
		cout << Char64out(cmore[0].vclean) << " cmore" << endl;
		if (!full) return;
		cout << "guas 2 cells" << endl;
		for (uint32_t i = 0; i <= ic2; i++)
			c2[i].DebugC2Clean();
	}

	void C2Status() {
		for (uint32_t i = 0; i <= ic2; i++)
			c2[i].DebugC2(0);
	}
	void C4Status() {
		for (uint32_t i = 0; i <= ic4; i++)
			c4[i].DebugMore(0);
	}

	void C3StatusClean() {
		cout << "status for c3 clean" << endl;
		for (uint32_t i = 0; i <= ic3; i++)
			c3[i].DebugC2Clean();
	}
	void Debug27() {
		cout << "chunk extract status n=" << nt12 << endl;
		for (uint32_t i = 0; i < nt12; i++)
			cout << Char27out(t12[i]) << endl;
	}
	void Stats() {
		cout << " chunks status ";
		cout << " ic2 " << ic2 << ";" << c2[ic2].nt;
		cout << " ic3 " << ic3 << ";" << c3[ic3].nt;
		cout << " ic4 " << ic4 << ";" << c4[ic4].nt;
		cout << " ic5 " << ic5 << ";" << c5[ic5].nt;
		cout << " icmore " << icmore << ";" << cmore[icmore].nt;
		cout << endl;
	}
}chunkh;
//_________________ minimum clues required in band 3 bands 1+2 locked
struct MINCOUNT {
	uint32_t mini_bf1, mini_bf2, mini_bf3,mini2all, mini_triplet,
		critbf, pairs27, mincount;
	inline void SetMincountG2() {// after direct setting minis
		mini2all = mini_bf1;
		mincount = _popcnt32(mini_bf1 | mini_triplet)
			+ _popcnt32(mini_bf3);
		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs
	}
	inline void AddMincountG3(uint32_t tr) {// after direct setting minis
		mini_triplet=tr;// shoul be  >0
		for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
			if (tr&bit)				critbf |= field;
	}
	void Status(const char * lib) {
		cout << lib << "critical Status mincount =" << mincount << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		if (mini_bf1)cout << Char9out(mini_bf1) << "     minis bf1" << endl;
		if (mini_bf2)cout << Char9out(mini_bf2) << "     minis bf2" << endl;
		if (mini_bf3)cout << Char9out(mini_bf3) << "     minis bf3" << endl;
		if (mini_triplet)cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;

	}
};
//================ Bands 

#define myband3 bax[2]
struct STD_B416 {
	char band[28];
	int band0[27], i416, gangster[9], map[27], map81[27], dband;
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
	void PrintUas() {
		cout << " band3 uas" << endl;
		for (uint32_t i = 0; i < nua; i++)
			cout << Char27out(tua[i]) << endl;
	}
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


struct GCHK {
	//=========== kwon puzzle filters
	BF128 puzknown;
	uint64_t pk54;
	int kn_ir1, kn_ir2;
	int aigstop, aigstopxy, start_perm, *tpw, *tsortw,
		a_18_seen,minb1b2, minb1, minb2,iperm;
	uint64_t diagknown7p ,diagtestadd;
	uint64_t debugvalbf;
	//___________________ studied solution 
	char * ze;// given solution grid
	char zes[164],ze_diag[164];
	char * zp;// first 18 if any
	char zsol[82];// solution grid morphed
	int grid0[81];// same 0 based
	int grid0_diag[81];// diagonal symmetry

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

	int nt4ok,okcheck;// for known
	// bands 1+2 valid epansion

	void Expand_03();

	int SetupExpand_46();
	void Expand_46();

	int SetupExpand_7p();
	void Init7p_guas();
	uint64_t ua_ret7p;
	int IsValid7p(SPB03*sn);
	inline int GetNextCell(SPB03* s);
	inline void GetNextUa( SPB03* sn);
	inline void GetNextUaAdd(SPB03* sn);
	inline int GetLastAndUa(SPB03* sn,int d=0);
	void Expand_7_11();
	void Expand_7_12();

	void GoBelow(SPB03* sn);
	void ExpandAddB1B2(SPB03* sn);

	void GoAfterExpand(SPB03* sn, uint32_t nadd=0);
	int IsValidB12(uint32_t ncl);

	//_____________ validb12
	uint64_t myb12, myb12f, 
		 myac, myacf, 
		myb12add, myacadd,anduab12;
	uint32_t  mynclues;//valid status
	int 	limb12;
	uint64_t tuaddb12[50];
	uint32_t tadd[50], ntadd, ntuaddb12;// add b1 b2
	uint32_t anduab3,stopexpandb3;// b3 expand

	// start B3 process 
	uint32_t tw3[1000], ntw3, andout ;
	uint32_t tw3_2[1000], ntw3_2,is1_tw3;

	int BuildB2Table();
	void BuildFinalTable();

	void B3FinalNoExpand();

	void ExpandB3Direct(int ntoass);
	void ExpandB3Vect(int ntoass);
	uint32_t NoExpandB3(uint32_t cl0bf );
	inline void BuildGua(BF128 & w, int cc);
	uint32_t IsValidB3(uint32_t bf);
	uint32_t NotThisPermB3(uint32_t bf);
	uint32_t tclues[40];// mini 25+band a
	int nclgo, nmiss;
	int  ncluesb3,mincluesb3;

	//================== clean process
	uint64_t wb12bf, wb12active, myua;
	uint32_t clean_valid_done,rclean_valid_done;
	uint32_t taddgob3[100], ntaddgob3;

	void Out17(uint32_t bfb3);

};


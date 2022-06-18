typedef union WUBUF_t{
	uint64_t tdophase2[30][4000];
	uint32_t tbuilexpand3[12][500];// 6+6
	uint32_t tg2[15][100];// in fact 10-(19 and more) 
}WUBUF;
WUBUF wubuf;
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

};
struct TUASB12 {
	uint64_t tua[4096], t2[3000];
	uint64_t   ta17[2000], ta18[2000], ta19[4096], ta20[4096], 
		ta21[4096], ta22[4096], tamore[4096];
	// 4096 size of vectors used in expansion
	uint32_t nua, nt2;
	uint32_t  nta17, nta18, nta19, nta20, nta21, nta22, ntamore;
	void SwitchTo54Mode() {
		nta17 = nta18 = nta19 = nta20 = nta21 = nta22 = ntamore = 0;
		for (uint32_t i = 0; i < nua; i++) {
			register uint64_t R = tua[i];
			R = (R & BIT_SET_27) | ((R& BIT_SET_B2) >> 5);
			tua[i]=R;// now r54
		}
	}
	inline void Add(uint64_t u, uint64_t cc);

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
			<< " nua=" << nua << "\tnta17=" << nta17
			<< " nta18=" << nta18 << " nta19=" << nta19
			<< " nta20=" << nta20 << " nta21=" << nta21
			<< " nta22=" << nta22 << " tot=" 
			<< nua+nta17+nta18+nta19+nta20+nta21+nta22	<< endl;
	}
};

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

struct GUASHORT {
	uint64_t killer;
	uint32_t   igua2_gua3, uab3, istart, iend;
}; 
struct GUAHANDLER {
	GUA guabuf[12*64];
	uint32_t nbuf,
		c2_i27[27],c3_i9[9], 
		c2index[27], c3index[9], nc2i, nc3i,
		c4index[200], c4uab3[200],nc4i,
		c5index[200], c5uab3[200], nc5i,
		cmindex[400], cmuab3[400], ncmi,
		g2,g3, nstartx;
	void Init() {
		nc2i = nc3i = nc4i = nc5i = ncmi = nbuf= 0;
	}
	void AddG2(uint32_t u3, uint32_t i27, uint64_t ua0) {
		register uint32_t bit = 1 << i27;
		if (!(g2&bit)) {
			c2index[nc2i++] = nbuf;
			guabuf[nbuf].Init(nbuf, u3, i27, ua0);
			c2_i27[i27] = nbuf++;
			g2 |=bit;
		}
		else {
			if (!ua0)guabuf[c2_i27[i27]].ReInit();
			else guabuf[c2_i27[i27]].Add(ua0);
		}

	}
	void AddG3( uint32_t i9, uint64_t ua0) {
		register uint32_t bit = 1 << i9;
		if (!(g3&bit)) {
			guabuf[nbuf].Init(nbuf, 7 << (3 * i9), i9, ua0);
			c3index[nc3i++] = nbuf;
			c3_i9[i9] = nbuf++;
			g3 |= bit;
		}
		else {
			if (!ua0)guabuf[c3_i9[i9]].ReInit();
			else guabuf[c3_i9[i9]].Add(ua0);
		}

	}
	void AddG4(uint32_t pat, uint64_t ua0) {
		int ipat = -1;
		for (uint32_t i = 0; i < nc4i; i++)		
			if (pat ==c4uab3[i]) {
				ipat = c4index[i];
				break;
			}
		if (ipat < 0) {
			guabuf[nbuf].Init(nbuf,	pat,0, ua0);
			c4uab3[nc4i] = pat;
			c4index[nc4i++] = nbuf++;
		}
		else {
			if (!ua0)guabuf[ipat].ReInit();
			else guabuf[ipat].Add(ua0);
		}
	}
	void AddG5(uint32_t pat, uint64_t ua0) {
		int ipat = -1;
		for (uint32_t i = 0; i < nc5i; i++)
			if (pat == c5uab3[i]) {
				ipat = c5index[i];
				break;
			}
		if (ipat < 0) {
			guabuf[nbuf].Init(nbuf, pat, 0, ua0);
			c5uab3[nc5i] = pat;
			c5index[nc5i++] = nbuf++;
		}
		else {
			if (!ua0)guabuf[ipat].ReInit();
			else guabuf[ipat].Add(ua0);
		}
	}
	void AddGm(uint32_t pat, uint64_t ua0) {
		int ipat = -1;
		for (uint32_t i = 0; i < ncmi; i++)
			if (pat == cmuab3[i]) {
				ipat = cmindex[i];
				break;
			}
		if (ipat < 0) {
			guabuf[nbuf].Init(nbuf,pat, 0, ua0);
			cmuab3[ncmi] = pat;
			cmindex[ncmi++] = nbuf++;
		}
		else {
			if (!ua0)guabuf[ipat].ReInit();
			else guabuf[ipat].Add(ua0);
		}
	}
	void AddGx(uint32_t pat, uint64_t ua0) {
		//cout << Char27out(pat) << " ";
		//cout << Char54out(ua0) << " to add guashadd "<<endl;

		uint64_t cc = __popcnt64(pat);
		switch (cc) {
		case 4: AddG4(pat, ua0); return;
		case 5: AddG5(pat, ua0); return;
		}
		AddGm(pat, ua0);
	}

	void Status(const char * lib) {
		cout << lib << " nbuf=" << nbuf ;
		cout << "\t nc2i=" << nc2i << " nc3i=" << nc3i << " nc4i=" << nc4i
			<< " nc5i=" << nc5i << " ncmi=" << ncmi << endl;

	} 
}guash;
struct GUASHORTHANDLER {
	uint64_t shortbuf[5000],
	vg2, vg3, vgm;// active guas2 and guas3 in first vkill
	GUASHORT guashort[12 * 64];
	uint32_t nbuf,nc2i, nc3i, 
		nshort, nvkill,g2, g3, nstartx,
		g2cell[54],g3cell[54];
	struct VKILL {// storing 64 uas and vectors 3 bands
		uint64_t  v0,
			vc[54], //cell for tu54
			vf, // vclean after a given band 1+2
			nt;
		inline void Init() {
			nt = 0, v0 = 0;
			memset(vc, 255, sizeof vc);
		}
		inline void Add(uint64_t u54) {//add a new killer
			uint64_t bit = (uint64_t)1 << nt++;
			v0 |= bit;
			uint32_t cc54;// build cells vectors
			register  uint64_t Rw = u54;
			while (bitscanforward64(cc54, Rw)) {
				Rw ^= (uint64_t)1 << cc54;// clear bit
				vc[cc54] ^= bit;
			}
		}
		inline void ApplyClean(uint32_t * tc, uint32_t ntc) {
			register uint64_t V = v0;
			for (uint32_t i = 0; i < ntc; i++)
				V &= vc[tc[i]];
			vf = V;
		}
	}vkill[12];

	void BuidAdd(GUASHORTHANDLER & o,  
		uint32_t * tc, uint32_t ntc, //cells to apply
		uint64_t bf, uint64_t ac) {//filter and active
		uint32_t lim2 = o.nc2i,lim3 = o.nc2i + o.nc3i;
		nbuf = nc2i = nc3i = nshort=0;
		for (uint32_t i64 = 0; i64 <= o.nvkill; i64++) {
			o.vkill[i64].ApplyClean(tc, ntc);
			register uint64_t V = o.vkill[i64].vf;
			uint32_t x;// build cells vectors
			while (bitscanforward64(x, V)) {
				uint64_t bit = (uint64_t)1 << x;
				V ^= bit;// clear bit
				GUASHORT wo = o.guashort[x + 64 * i64],
					w=wo;
				w.killer &= ac;
				if (wo.iend == wo.istart){//  aready one ua
					w.iend = w.istart = nshort;
					guashort[nbuf++] = w;
				}
				else {// minimum 2 uas
					w.istart = w.iend = nshort;
					w.killer = ~ 0;
					for (uint32_t j = wo.istart; j < wo.iend; j++) {
						register uint64_t U = o.shortbuf[j];
						if (U & bf) continue;
						U &= ac;
						if (!U) {// empty ua never hit
							w.killer = 0;
							w.iend = w.istart;
							break;
						}
						w.killer &= U;
						shortbuf[w.iend++]=U;
					}
					if (w.killer != ~0) {// a new to store
						if ((w.iend - w.istart) > 1) 	nshort = w.iend;
						else w.iend = w.istart;	// no det if one ua		 
						guashort[nbuf++] = w;
					}
				}
				if (w.killer != ~0 && (!i64) && x<lim3) {
					if (x < lim2) nc2i++;
					else  nc3i++;
				}
			}
		}
		vg2 = maskLSB[nc2i].u64[0];
		vg3 = maskLSB[nc3i].u64[0];
		vg3 <<= nc2i;
		vgm = ~(vg2 | vg3);
	}

	void DoVkill() {
		nvkill = (nbuf + 63) >> 6;
		for (uint32_t i = 0; i < nvkill; i++)vkill[i].Init();
		//cout << "nbuf="<<nbuf << " nvkill=" << nvkill << endl;
		for (uint32_t i = 0; i < nbuf; i++) {
			vkill[i >> 6].Add(guashort[i].killer);
		}
	}


	//___________ extract others for band 3
	uint32_t t12[1100], nt12;
	inline void GetB12CX(uint32_t * tc, uint32_t ntc, uint64_t bf) {
		for (uint32_t i = 1; i <= nvkill; i++)
			vkill[i].ApplyClean(tc, ntc);
		nt12 = 0;
		for (uint32_t i64 = 0; i64 < nvkill; i64++) {
			register uint64_t V = vkill[i64].vf;
			uint32_t x;
			if (!i64)V &= vgm;
			while (bitscanforward64(x, V)) {
				uint64_t bit = (uint64_t)1 << x;
				V ^= bit;// clear bit
				GUASHORT w = guashort[x + 64 * i64];
				register  uint32_t pat = w.uab3;
				if (w.istart == w.iend) // killer is the only ua (not hit )
					t12[nt12++] = pat;
				else {// look for not hit uas
					for (uint32_t i = w.istart; i <  w.iend; i++) {
						register uint64_t U = shortbuf[i];
						if (!(bf & U)) {  // not hit first ok
							t12[nt12++] = pat; break;
						}
					}
				}
			}
		}
	}

	void DumpPack() {
		cout << "dump pack2 nshort=" << nshort 
			<<" nbuf="<< nbuf << endl;
		for (uint32_t i = 0; i < nbuf; i++) {
			GUASHORT & ws = guashort[i];
			cout << Char27out(ws.uab3) << "\ti9/27=" << ws.igua2_gua3 
				<< "\t" << ws.istart << "\t" << ws.iend << "\t ";
			cout << Char54out(ws.killer) << " killer ndet=" << ws.iend - ws.istart << endl;
			for (uint32_t j = ws.istart; j < ws.iend; j++)
				cout << Char54out(shortbuf[j]) << "  " << j << endl;
		}

	}
	void Debug27() {
		cout << "guash2 extract status n=" << nt12 << endl;
		for (uint32_t i = 0; i < nt12; i++)
			cout << Char27out(t12[i]) << endl;
	}

	void Status(const char * lib) {
		cout << lib << " nshort=" << nshort
			<< " nbuf=" << nbuf << " nvkill=" << nvkill;
		cout << "\t nc2i=" << nc2i << " nc3i=" << nc3i  << endl;

	}
}guash2, guash2add;
struct GUASHORTHANDLER3 {
	struct G2G3 {// storing detail for a g2/g3
		uint64_t killer;
		uint32_t uab3, nt;
		void Init(uint32_t u) { killer = ~0; nt = 0; uab3 = u; }
	}tg2g3[36];
	uint64_t  tug2g3[36][30];
	uint64_t shortbuf[5000];
	GUASHORT guashort[12 * 64];
	uint32_t nbuf, nc2i, nc3i,
		nshort, nvkill, nstartx;
	uint32_t ng2_128, ng3_128, nblocg2, nblocg3;
	struct VKILL {// storing 64 uas and vectors 3 bands
		uint64_t  v0,
			vc[54], //cell for tu54
			vf, // vclean after a given band 1+2
			nt;
		inline void Init() {
			nt = 0, v0 = 0;		memset(vc, 255, sizeof vc);	}
		inline void Add(uint64_t u54) {//add a new killer
			uint64_t bit = (uint64_t)1 << nt++;
			v0 |= bit;
			uint32_t cc54;// build cells vectors
			register  uint64_t Rw = u54;
			while (bitscanforward64(cc54, Rw)) {
				Rw ^= (uint64_t)1 << cc54;// clear bit
				vc[cc54] ^= bit;
			}
		}
		inline void ApplyClean(uint32_t * tc, uint32_t ntc) {
			register uint64_t V = v0;
			for (uint32_t i = 0; i < ntc; i++)
				V &= vc[tc[i]];
			vf = V;
		}
		void Dump() {
			cout << "dump vkill bloc" << endl;
			cout << Char64out(v0) << " v0" << endl;
			for (int i = 0; i < 54; i++)
				cout << Char64out(vc[i]) << " cell=" << i << endl;
		}
	}vkill[12];
	void Pack(GUAHANDLER & o);
	void Add128(BF128 w, int cc) {
		if (cc == 2) {
			uint32_t i27 = w.bf.u32[2];
			G2G3 & wg2g3 = tg2g3[i27];
			if (wg2g3.nt < 30) {
				register uint64_t U = w.bf.u64[0];
				wg2g3.killer &= U;
				tug2g3[i27][wg2g3.nt]=U;
			}
		}
		else if (cc == 3) {
			uint32_t i9 = w.bf.u32[2];
			G2G3 & wg2g3 = tg2g3[i9+27];
			if (wg2g3.nt < 30) {
				register uint64_t U = w.bf.u64[0];
				wg2g3.killer &= U;
				tug2g3[i9 + 27][wg2g3.nt] = U;
			}
		}
		else {//add it as new if room for it
			if (nbuf >= (12 * 64)) return;
			register uint32_t nbufr = nbuf;
			GUASHORT & ws = guashort[nbuf++];
			ws.igua2_gua3 = 0; //unused here
			ws.uab3 =  w.bf.u32[2];
			ws.killer = w.bf.u64[0];
			ws.istart = nshort;
			ws.iend = nshort;
			if ((nbuf & 63) == 1) // we have a new 64 bloc
				vkill[nvkill++].Init();
			vkill[nbufr >> 6].Add(guashort[nbufr].killer);
		}
	}
	int BuildAddApply(int last);
	void BuidAdd(GUASHORTHANDLER & o,
		uint32_t * tc, uint32_t ntc, //cells to apply
		uint64_t bf, uint64_t ac) {//filter and active
		uint32_t lim2 = o.nc2i, lim3 = o.nc2i + o.nc3i;
		nbuf = nc2i = nc3i = nshort = 0;
		for (uint32_t i64 = 0; i64 <= o.nvkill; i64++) {
			o.vkill[i64].ApplyClean(tc, ntc);
			register uint64_t V = o.vkill[i64].vf;
			uint32_t x;// build cells vectors
			while (bitscanforward64(x, V)) {
				uint64_t bit = (uint64_t)1 << x;
				V ^= bit;// clear bit
				GUASHORT wo = o.guashort[x + 64 * i64],
					w = wo;
				w.killer &= ac;
				if (wo.iend == wo.istart) {//  aready one ua
					w.iend = w.istart = nshort;
					guashort[nbuf++] = w;
				}
				else {// minimum 2 uas
					w.istart = w.iend = nshort;
					w.killer = ~0;
					for (uint32_t j = wo.istart; j < wo.iend; j++) {
						register uint64_t U = o.shortbuf[j];
						if (U & bf) continue;
						U &= ac;
						if (!U) {// empty ua never hit
							w.killer = 0;
							w.iend = w.istart;
							break;
						}
						w.killer &= U;
						shortbuf[w.iend++] = U;
					}
					if (w.killer != ~0) {// a new to store
						if ((w.iend - w.istart) > 1) 	nshort = w.iend;
						else w.iend = w.istart;	// no det if one ua		 
						guashort[nbuf++] = w;
					}
				}
				if (w.killer != ~0 && (!i64) && x < lim3) {
					if (x < lim2) nc2i++;
					else  nc3i++;
				}
			}
		}
	}
	//___________ extract others for band 3
	uint32_t t12[1100], nt12;
	inline void GetB12CX(uint32_t * tc, uint32_t ntc, uint64_t bf) {
		for (uint32_t i = 1; i <= nvkill; i++)
			vkill[i].ApplyClean(tc, ntc);
		nt12 = 0;
		for (uint32_t i64 = 0; i64 < nvkill; i64++) {
			register uint64_t V = vkill[i64].vf;
			uint32_t x;
			while (bitscanforward64(x, V)) {
				uint64_t bit = (uint64_t)1 << x;
				V ^= bit;// clear bit
				GUASHORT w = guashort[x + 64 * i64];
				register  uint32_t pat = w.uab3;
				if (w.istart == w.iend) // killer is the only ua (not hit )
					t12[nt12++] = pat;
				else {// look for not hit uas
					for (uint32_t i = w.istart; i < w.iend; i++) {
						register uint64_t U = shortbuf[i];
						if (!(bf & U)) {  // not hit first ok
							t12[nt12++] = pat; break;
						}
					}
				}
			}
		}
	}
	void Dumpg2g3() {
		cout << "dump g2g3" << endl;
		for (uint32_t i = 0; i < 36; i++) {
			G2G3 & w = tg2g3[i];
			if (w.killer == (~0))continue;
			cout <<i<<"\t"<< Char27out(w.uab3) << " ";
			cout << Char54out(w.killer) 
				<<" nt="<<w.nt<< endl;

		}

	}
	void Dumpg2g3Det() {
		cout << "dump g2g3 detail" << endl;
		for (uint32_t i = 0; i < 36; i++) {
			G2G3 & w = tg2g3[i];
			uint64_t * t = tug2g3[i];
			if (!w.nt) continue;
			cout << i << "\t" << Char27out(w.uab3) << " ";
			cout << Char54out(w.killer)
				<< " nt=" << w.nt << endl;
			for (uint32_t j = 0; j < w.nt; j ++ )
				cout << Char54out(t[j]) << endl;

		}

	}	
	void DumpPack() {
		cout << "dump pack2 nshort=" << nshort
			<< " nbuf=" << nbuf << endl;
		for (uint32_t i = 0; i < nbuf; i++) {
			GUASHORT & ws = guashort[i];
			cout << Char27out(ws.uab3) << "\ti9/27=" << ws.igua2_gua3
				<< "\t" << ws.istart << "\t" << ws.iend << "\t ";
			cout << Char54out(ws.killer) << " killer ndet=" << ws.iend - ws.istart << endl;
			for (uint32_t j = ws.istart; j < ws.iend; j++)
				cout << Char54out(shortbuf[j]) << "  " << j << endl;
		}

	}
	void Debug27() {
		cout << "guash2 extract status n=" << nt12 << endl;
		for (uint32_t i = 0; i < nt12; i++)
			cout << Char27out(t12[i]) << endl;
	}
	void Status(const char * lib) {
		cout << lib << " nshort=" << nshort
			<< " nbuf=" << nbuf << " nvkill=" << nvkill<< endl;

	}
}guash3, guash3add;
struct G2_256 {//vector 256 chunk level gua2 gua3
	BF128 v[2], vc[2][54],bfadd;
	uint32_t bits[256];// i27 or i9
	uint32_t nv;
	void Init() {
		nv = 0;
		memset(v, 0, sizeof v);
		memset(vc, 255, sizeof(vc));
	}
	void Addx(BF128 &vx, BF128 *vcx, uint32_t* bitx, uint32_t n) {
		vx.Set(n);
		bitx[n] = bfadd.bf.u32[2];
		register uint64_t V = bfadd.bf.u64[0];
		register uint32_t x;
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;// clear bit
			vcx[x].clearBit(n);
		}
	}
	void Add(BF128 gx) { // gx is a still valid gua2 gua3 
		bfadd = gx;
		if (nv > 256)return; // should never be
		if (nv < 128)Addx(v[0], vc[0], bits, nv);
		else Addx(v[1], vc[1], 	&bits[128], nv-128);
		nv++;
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
		for (uint32_t i = 0; i < 4; i++) {// extract
			register uint64_t V = vw[0].bf.u64[i];
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
			cout << endl;
		}
	}
	
}g2_256[2], g3_256,gm_256[2];
struct G2_256Handler {
	uint32_t n2, n3, nm;
	inline void GetStart(uint32_t n2e, uint32_t n3e, uint32_t nme) {
		n2 = n2e; n3 = n3e, nm = nme;
	}
	inline void Add128(BF128 w, uint32_t cc) {
		if (cc < 4) {
			w.bf.u32[2]=1<<w.bf.u32[2];// switch to bit field
			if (cc == 2) {
				//cout << "add 2 in g2_256 n2=" << n2 << endl;
				if (n2++ < 256)g2_256[0].Add(w);
				else g2_256[1].Add(w);
			}
			else { g3_256.Add(w);n3++; }
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
	}


}g_256h;
struct CHUNK3B {// storing 64 uas and vectors 3 bands
	uint64_t tu54[64], v0, 
		vc[54], //cell for tu54
		vclean, // v0 after claean groupo
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

	inline int  Check2(BF128 w54);
	inline int  Check3(BF128 w54);
	inline int  Check4(BF128 w54);
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
		register uint64_t R1 = vect.bf.u64[0];
		for (uint32_t i = 0; i < ntcells; i++)
			if(!(R1 &= vc[tcells[i]].bf.u64[0])) break;
		if(R1) return 1;
		if (nt > 64) {
			register uint64_t R2 = vect.bf.u64[1];
			for (uint32_t i = 0; i < ntcells; i++)
				if (!(R2 &= vc[tcells[i]].bf.u64[1]))
					return 0;
			return 1;
		}
		return 0;
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
				if(ntd<100)td[ntd++] = t[ir];
			}
		}
		{
			register uint64_t R = w.bf.u64[1];
			while (bitscanforward64(ir, R)) {
				R ^= (uint64_t)1 << ir;
				if (ntd < 100)td[ntd++] = t[ir+64];
			}
		}
	}
	void Dump() {
		for (uint32_t i = 0; i < nt; i++)
			cout << Char54out(t[i])<< " "<<i<<"\t"
			<<_popcnt64(t[i])<< endl;
	}

};
struct MORE64VHandler {// vectors for a given chunk
	MORE64VECT mv[10];
	uint32_t nmv,n;
	inline void Init() { mv[0].Init(); nmv = 0; }

	inline void Add54(uint64_t v) {//add a new more if <128 
		if (mv[nmv].nt == 128) {
			if (nmv < 10)mv[++nmv].Init();
			else return;
		}
		mv[nmv].Add54(v);
	}
	inline int ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		for(uint32_t i=0;i<=nmv;i++)
		if (mv[i].ApplyXY(tcells, ntcells)) return 1;
		return 0;
	}
	inline void Extract(uint32_t *tc, uint32_t ntc,
		uint64_t *td, uint64_t & ntd) {
		for (uint32_t i = 0; i <= nmv; i++)
			mv[i].Extract(tc, ntc, td, ntd);
	}
	void Dump(int all = 1) {
		cout << " nmv=" << nmv	<< " mv[nmv].nt=" << mv[nmv].nt << endl;
		if (!all)return;
		for (uint32_t i = 0; i <= nmv; i++) mv[i].Dump();
	}

}m64vh;

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

	int diagbugchk,diagtestadd;
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
	}t4_to_expand[10000];
	uint64_t tua4[4096];// also tuaclean
	uint32_t ntua4;
	uint64_t tuaclean[4096]; 
	uint32_t ntuaclean;
	BF128 vph2[32], vph2c[32][54];

	uint64_t v12_5_v0[64], v12_5_c[64][54]; //maxi 4096 uas
	uint32_t nv12_5,nvph2;

	void AddTua4(uint64_t v);

	uint64_t bufvalid[BUFVALIDS+50], *pbufvalid,*pendbufvalid;
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
	struct GADDB3 { //Add during a clean chunk
		BF128 tmore[384], t2[256], t3[256];
		uint32_t ntmore, nt2, nt3;
		inline void Init() { ntmore = nt2 = nt3 = 0; }
		inline void Add(BF128 w, int cc) {
			if (cc < 3) {
				if (nt2 < 256)t2[nt2++] = w;
			}
			else if (cc == 3) {
				if (nt3 < 256)t3[nt3++] = w;
			}
			else if (ntmore < 384)
				tmore[ntmore++] = w;
		}
		inline void Get3(uint64_t F, uint32_t &guas3) {
			for (uint32_t i = 0; i < nt3; i++)
				if (!(F & t3[i].bf.u64[0]))
					guas3 |= 1 << t3[i].bf.u32[2];
		}
		inline void GetMore(uint64_t F,
			uint32_t *t, uint32_t &nt) {
			for (uint32_t i = 0; i < ntmore; i++)
				if (!(F & tmore[i].bf.u64[0]))
					t[nt++] =tmore[i].bf.u32[2];
		}
		void Dumpt2() {
			cout << "dumpt2 nt2=" << nt2 << endl;
			for (uint32_t i = 0; i < nt2; i++) {
				cout << Char54out(t2[i].bf.u64[0])
					<< " " << t2[i].bf.u32[2] << endl;
			}
		}
		void Dumpt3() {
			cout << "dumpt3 nt3=" << nt3 << endl;
			for (uint32_t i = 0; i < nt3; i++) {
				cout << Char54out(t3[i].bf.u64[0])
					<< " " << t3[i].bf.u32[2] << endl;
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
	struct MOREAND {
		uint64_t tm[700], ntm;
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
		inline uint64_t GetAnd(uint64_t x) {
			for (uint64_t i = 0; i < ntm; i++)
				 x &=tm[i];
			return x;
		}
		inline uint64_t SetAnd() {
			if (ntm < 2) return 1;
			for (uint64_t i = 1; i < ntm; i++)
				tm[0] &= tm[i];
			ntm = 1;
			return tm[0];
		}
		void Dump() {
			for (uint64_t i = 0; i < ntm; i++)
				cout << Char54out(tm[i]) << endl;
		}
		void DumpAnd(uint64_t a) {
			for (uint64_t i = 0; i < ntm; i++)
				cout << Char54out(tm[i]&a) << endl;
		}
	}mbisvalid, mab;

	void BuildVectorsForExpand4B12();//64 uas
	void Expand4B12();
	void Do_phase2(T4_TO_EXPAND w);
	void Do_phase2Expand(uint64_t bf, uint64_t ac);
	void Do_phase2Expand_128(uint64_t bf, uint64_t ac);
	inline int ApplyCountFilterB(uint64_t bf);
	inline void LoadRemainningClues(uint64_t bf);
	int IsValidB12(uint32_t ncl);


	void CleanBufferAndGo(uint64_t andvalid, uint64_t orvalid);
	void CheckValidBelow(uint64_t bf,uint64_t ac);
	void BuildReducedGuasVectors23();

	void CleanMoreUas(uint64_t bf, uint64_t ac, int ncl, MOREANDB & mabo);
	//_____________ validb12
	uint64_t myandvalid, myorvalid, 
		myb12, myb12f, 
		myac_4, myac, myacf,
		myb12add, myacadd;
	uint32_t  mynclues;//valid status
	int 	limb12;
	uint32_t tadd[50], ntadd;// add b1 b2
	uint32_t anduab3,stopexpandb3;// b3 expand


	void ExpandAddB1B2(uint64_t bf);
	void ExpandAddB1B2Go(int step);

	// start B3 process 
	uint32_t tw3[1000], ntw3, andout ;
	uint32_t tw3_2[1000], ntw3_2,is1_tw3;
	int GetB3BaseB();

	int BuildB2Table();
	int BuildFinalTable();
	void InitGoB3(uint64_t bf,uint32_t ncl);

	void B3FinalNoExpand();

	void ExpandB3Direct(int ntoass);
	void ExpandB3Vect(int ntoass);
	uint32_t NoExpandB3(uint32_t cl0bf );
	inline void BuildGua(BF128 & w, int cc);
	uint32_t IsValidB3(uint32_t bf);
	uint32_t tclues[40], *tcluesxy;// mini 25+band a
	int nclues_step,ncluesb12,nclgo, nclf,nmiss;
	int  ncluesb3,mincluesb3;
	uint32_t gguas2, gguas3;// 27;9 bits active guas 2 3

	//================== clean process
	uint64_t wb12bf, wb12active, myua;
	uint32_t clean_valid_done;
	uint32_t taddgob3[100], ntaddgob3;
	MORE32 moreuas_b3;

	void Out17(uint32_t bfb3);

};


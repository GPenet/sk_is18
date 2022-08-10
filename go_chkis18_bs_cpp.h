//#define HAVEKNOWN
#ifdef HAVEKNOWN
int GCHK::Kn0() {
	// put knownn in the right order
	if (strlen(ze) < 163) return -1;// skip blank lines
	char* w = &ze[82], * ww = &zes[82];
	for (int ib = 0; ib < 3; ib++) {
		int* m = bax[ib].map81;
		for (int c = 0; c < 27; c++) {
			w[27 * ib + c] = ww[m[c]];
		}
	}
	puzknown.SetAll_0();
	for (int i = 0; i < 81; i++) {
		char c = w[i];
		if (c < '1' || c>'9') continue;
		puzknown.Set_c(i);// store the pattern in 3X mode
		// this can be a pattern, no check of the digit with the solution
	}
	register uint64_t R = puzknown.bf.u64[0];
	pk54 = (R & BIT_SET_27) | ((R & BIT_SET_B2) >> 5);
	int n1 = (int)_popcnt32(puzknown.bf.u32[0]),
		n2 = (int)_popcnt32(puzknown.bf.u32[1]),
		n3 = (int)_popcnt32(puzknown.bf.u32[2]);
	//for (int i = 0; i < 81; i++) cout << ze[i + 82];
	//cout << " known reshaped  "<<n1<<n2<<n3 << endl;

	if (n3 < mincluesb3) return  0;
	if (n3 == 6 && n2 != 6) return  0;
	if (n3 < n1 || n3 < n2) return  0;
	if (n3 == 7 && n2 < 4) return  0;
	BF128 ws1 = puzknown & band3xBM[3],
		ws2 = puzknown & band3xBM[4],
		ws3 = puzknown & band3xBM[5];
	if (n3 < ws1.Count() || n3 < ws2.Count() || n3 < ws3.Count()) return 0;
	cout << ws1.Count() << " " << ws2.Count() << " " << ws3.Count() << endl;
	return 1;
}
void GCHK::Kn02() {
	for (int ibs = 0; ibs < 3; ibs++) {
		STD_B416& b = bax[ibs];
		cout << b.band ;
	}
	cout << " reshaped  "  << endl;
	for (int i = 0; i < 81; i++) cout << ze[i + 82];
	cout << " known reshaped iperm="<< iperm << endl;
}
#endif

//_______________ start  a band3 perm and UAs GUAs search 
#define AddUaB12UN(X,Y)AddUA64(tuasb12.tua,tuasb12.nua,X|((uint64_t)Y<<59))
void STD_B416::InitBand2_3(char * ze, BANDMINLEX::PERM & p
	, int iband) {
	i416 = p.i416;
	dband = 27 * iband;
	GetUAs();
	strncpy(band, ze, 27);
	for (int i = 0; i < 27; i++) band0[i] = band[i] - '1';
	// create the cell map in output
	for (int i = 0; i < 3; i++) {
		int vr = 9 * p.rows[i], vr0 = 9 * i;
		for (int j = 0; j < 9; j++)
			map[vr0 + j] = vr + p.cols[j];
	}
	MorphUas();// morph all uas
	SetGangster();
}

//_______________ process band ordered
int  GCHK::StartIs18() {
	if (aigstop) return 0;
	p_cpt2g[0]++;
#ifdef HAVEKNOWN
	if(!Kn0()) return 0;
	//return 0;
#endif
	int * zs0 = grid0, *zs0_diag = grid0_diag;
	BANDMINLEX::PERM perm_ret;
	memcpy(zs0, bax[0].band0, sizeof  bax[0].band0);
	memcpy(&zs0[27], bax[1].band0, sizeof  bax[0].band0);
	memcpy(&zs0[54], bax[2].band0, sizeof  bax[0].band0);
	minb1 = t16_min_clues[bax[0].i416] ;
	minb2 =  t16_min_clues[bax[1].i416];
	minb1b2 = minb1 + minb2;
	limb12 = 18 - mincluesb3;
	for (int i = 0; i < 81; i++)zs0_diag[i] = zs0[C_transpose_d[i]];
	for (int i = 0; i < 81; i++)ze_diag[i] = ze[C_transpose_d[i]];
	for (int ibs = 0; ibs < 3; ibs++) {
		bandminlex.Getmin(&zs0_diag[27 * ibs], &perm_ret);
		bax[ibs + 3].InitBand2_3(&ze_diag[27 * ibs], perm_ret, ibs);
	}
#ifdef HAVEKNOWN
	Kn02();
#endif
	//__________________________ start uas search
	UaCollector();
	p_cpt2g[1] = tuasb12.nua;
	p_cpt2g[2] = chunkh.GetC2Count();
	Expand_03();
	return nok;
}
//________________ uas generation and store in chunkh
void GCHK::UaCollector() {
	Sg2Setup();
	FirstUasCollect();// uas 2 digits bands and stacks
	SecondUasCollect();// uas 345 bands 1+2
	UasCollect4box();// 4 box in bands 1+2
	PutUasStartInVector();
	Guas2Collect();
}

inline void GCHK::Adduab12() {
	if (tuasb12.nua > 2550)return;
	for (uint32_t i = 0; i < zh2gxn.nua; i++) {
		uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);
		AddUaB12UN(w, cc);
	}
}

BF128 w_tfua[24][200];// put it out of the stack
void GCHK::FirstUasCollect() {// produce all uas 2/3 digits
	tuas81.ntall = 0;// no ua so far
	zhgxn.SetupFsol(grid0);
	for (int i = 0; i < 36; i++) {
		uint32_t myf = floors_2d[i];
		zhou2[0].GoZ2(myf);
		for (uint32_t j = 0; j < zhgxn.nua; j++) {
			BF128 w = zhgxn.tua[j];
			tuas81.Add(w, myf);
		}
	}
	// insert bands and stacks and apply subsets
	//BF128 wubuf_tfua[24][200]; 
	uint32_t ntsort[27 - 4];
	memset(ntsort, 0, sizeof ntsort); // empty table
	{ //sort all by size
		BF128 * t = tuas81.tall;
		uint32_t nt = tuas81.ntall;
		for (uint32_t iua = 0; iua < nt; iua++) {
			BF128 wt = t[iua];
			uint32_t cc = wt.bf.u16[7] - 4;// index is 0 for 4
			wt.bf.u32[3] = 0; // clear digits and count
			w_tfua[cc][ntsort[cc]++] = wt;
		}
	}
	// here insert all missing uas one band/stack
	for (int i = 0; i < 3; i++) {
		BF128 wt;
		wt.SetAll_0();
		uint32_t *tu = bax[i].tua, nu = bax[i].nua;
		for (uint32_t j = 0; j < nu; j++) {
			register uint32_t u = tu[j] & BIT_SET_27, cc = _popcnt32(u);
			//if (cc < 7)continue;// all 2/3 digits covered
			cc -= 4;
			wt.bf.u32[i] = u;
			w_tfua[cc][ntsort[cc]++] = wt;
		}
	}

	for (int i = 3; i < 6; i++) {
		uint32_t dk = 3 * (i - 3); // stack adjustment
		BF128 wt;
		uint32_t *tu = bax[i].tua, nu = bax[i].nua;
		for (uint32_t j = 0; j < nu; j++) {
			register uint32_t u = tu[j] & BIT_SET_27, cc = _popcnt32(u);
			//if (cc < 7)continue;// all 2/3 digits covered
			wt.SetAll_0();
			for (uint32_t k = 0, bit = 1; k < 27; k++, bit <<= 1) if (u&bit) {
				uint32_t cell = C_transpose_d[k] + dk;
				wt.Set_c(cell);
			}
			cc -= 4;
			w_tfua[cc][ntsort[cc]++] = wt;
		}
	}
	// reload tall and check subsets/redundancy
	tuas81.ntall = 0;
	for (int i = 0; i < 23; i++)if (ntsort[i]) {
		BF128 * tt = w_tfua[i];
		for (uint32_t j = 0; j < ntsort[i]; j++) {
			tuas81.Add2(tt[j], i + 4);
		}
	}
	// split uas band 12 and others
	tuas81.ntold = tuasb12.nua = 0;
	for (uint32_t iua = 0; iua < tuas81.ntall; iua++) {
		BF128 wt = tuas81.tall[iua];
		if (wt.bf.u32[2]) tuas81.told[tuas81.ntold++] = wt;
		else {
			uint64_t w = wt.bf.u64[0], cc = _popcnt64(w);
			w |= cc << 59;
			tuasb12.tua[tuasb12.nua++] = w;
		}
	}
}
void GCHK::SecondUasCollect() {// collect 345 digits in bands 1+2
	zh2gxn.InitKnown(tuasb12.t2, &tuasb12.nt2);
	zh2gxn.SetupFsol(grid0);
	for (int i = 0; i < 84; i++) {
		uint32_t myf = floors_3d[i];
		int ir = zh2_3[0].GoZ3(myf);// cells unsolved count
		if (ir < 6) continue;// minimum for a fresh ua 3 digits
		uint64_t F = zh2_3[0].cells_unsolved.bf.u64;
		tuasb12.GetT2(F);
		zh2_3[0].DoZ3Go();
		if (zh2gxn.nua) Adduab12();
	}
	for (int i = 0; i < 126; i++) {
		int ir = zh2_4[0].GoZ4(floors_4d[i]);
		if (ir < 8) continue;// minimum for a fresh ua 4 digits
		uint64_t F = zh2gxn.unsolved_field;
		tuasb12.GetT2(F);
		zh2_4[0].DoZ4Go();
		if (zh2gxn.nua) Adduab12();
	}

	for (int i = 0; i < 126; i++) {
		int ir = zh2_5[0].GoZ5(0777 ^ floors_4d[i]);
		uint64_t F = zh2gxn.unsolved_field;
		tuasb12.GetT2(F);
		zh2_5[0].DoZ5Go();
		if (zh2gxn.nua) Adduab12();
	}
}
void GCHK::UasCollect4box() {
	//___________________ box 1245
	BF64 wf4;
	wf4.bf.u32[0] = wf4.bf.u32[1] = 077077077;// stack 3 bands 12
	tuasb12.GetT2(wf4.bf.u64);
	//tuasb12.DumpT2();
	zh2b[0].InitB1245(grid0);
	zh2b[0].Do4bGo();
	if (zh2gxn.nua) Adduab12();
	//___________________ box 1346
	wf4.bf.u32[0] = wf4.bf.u32[1] = 0707707707;
	tuasb12.GetT2(wf4.bf.u64);
	zh2b[0].InitB1346(grid0);
	zh2b[0].Do4bGo();
	if (zh2gxn.nua) Adduab12();
	//___________________ box 2356
	wf4.bf.u32[0] = wf4.bf.u32[1] = 0770770770;
	tuasb12.GetT2(wf4.bf.u64);
	zh2b[0].InitB2356(grid0);
	zh2b[0].Do4bGo();
	if (zh2gxn.nua) Adduab12();
}
void GCHK::Guas2Collect() {

	// guas 2 digits and stack g2s already there 
	// extract known G2s in ua64 mode

	uint64_t tg2uas[27][5];
	uint32_t ntg2uas[27];
	memset(ntg2uas, 0, sizeof ntg2uas);
	{
		uint64_t * tu54 = chunkh.c2[0].tu54, nt = chunkh.c2[0].nt;
		uint32_t * tu27 = chunkh.c2[0].tu27;
		for (uint64_t i = 0; i < nt; i++) {
			int i27 = tu27[i];
			uint64_t R = tu54[i], ua2x = (R&BIT_SET_27) | ((R >> 27) << 32);
			tg2uas[i27][ntg2uas[i27]++] = ua2x;
		}
	}
	for (int i27 = 0; i27 < 27; i27++) {// 27 g2
		SG2 w = tg2[i27];
		uint64_t uasi27[200], nuasi27 = ntg2uas[i27];
		memcpy(uasi27, tg2uas[i27], sizeof tg2uas[0]);
		if (nuasi27 == 1 && _popcnt64(uasi27[0]) == 2) continue;// nothing to do

		// __________________find guas 3 digits
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			int fl = floors_3d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			int ir = zh2_3[0].GoZ3G2(fl, w.col1, w.dig1, w.col2, w.dig2);
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				uint64_t U = zh2gxn.tua[0], cc = _popcnt64(U);
				if (cc < 5) continue;
				uasi27[nuasi27++] = U;
				for (uint64_t i = 0; i < nuasi27 - 1; i++)
					if (U == uasi27[i]) { nuasi27--; break; }
				continue;
			}
			uint64_t F = zh2gxn.unsolved_field;
			tuasb12.GetT2(F, uasi27, nuasi27);
			int ir2 = zh2_3[0].DoZ3Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					uasi27[nuasi27++] = U;
				}
		}
		//_________________find uas 4 digits
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			int fl = floors_4d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			int ir = zh2_4[0].GoZ4G2(fl, w.col1, w.dig1, w.col2, w.dig2);
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				uint64_t U = zh2gxn.tua[0], cc = _popcnt64(U);
				if (cc < 7) continue;
				uasi27[nuasi27++] = U;
				for (uint64_t i = 0; i < nuasi27 - 1; i++)
					if (U == uasi27[i]) { nuasi27--; break; }
				continue;
			}
			uint64_t F = zh2gxn.unsolved_field;
			tuasb12.GetT2(F, uasi27, nuasi27);
			int ir2 = zh2_4[0].DoZ4Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					uasi27[nuasi27++] = U;
					}
		}
		//_________________find uas 5 digits
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			int fl = 0777 ^ floors_4d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			int ir = zh2_5[0].GoZ5G2(fl, w.col1, w.dig1, w.col2, w.dig2);
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				uint64_t U = zh2gxn.tua[0];
				if (_popcnt64(U) < 9) continue;
				uasi27[nuasi27++] = U;
				for (uint64_t i = 0; i < nuasi27 - 1; i++)
					if (U == uasi27[i]) { nuasi27--; break; }
				continue;
			}
			uint64_t F = zh2gxn.unsolved_field;
			tuasb12.GetT2(F, uasi27, nuasi27);
			int ir2 = zh2_5[0].DoZ5Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i], cc = _popcnt64(U);
					if (cc > 20) continue;
					uasi27[nuasi27++] = U;
				}
		}		
		// store results in chunkh.c2
		uint64_t *tk = tg2uas[i27], ntk = ntg2uas[i27];
		for (uint64_t i = 0; i < nuasi27; i++) {
			register uint64_t U2x = uasi27[i];
			for (uint64_t j = 0; i < ntk; j++) // redundancy
				if (U2x == tk[j]) { U2x = 0; break; }
			if (U2x) chunkh.Addc2(U2x, i27);
		}
	}
}

inline void GCHK::BuildGua(BF128 & w, int cc) {
	register uint64_t ua12 = w.bf.u64[0];
	uint64_t ua54 = (ua12 & BIT_SET_27) |
		((ua12 & BIT_SET_B2) >> 5);
	w.bf.u64[0] = ua54;// store in 54 cells mode
	if (cc < 4) {
		register uint32_t A = w.bf.u32[2], B, i;
		if (cc == 2) 	w.bf.u32[2] = GetI27(A);
		else {
			for (i = 0, B = 7; i < 9; i++, B <<= 3) 
				if (A&B) { w.bf.u32[2] = i; break; }
		}
	}
}
void  GCHK::PutUasStartInVector() {	// split guas 2/3 and others
	chunkh.Init();
	for (uint32_t iua = 0; iua < tuas81.ntold; iua++) {
		BF128 w = tuas81.told[iua];
		if (w.bf.u64[0]) {
			int cc = _popcnt32(w.bf.u32[2]);
			BuildGua(w, cc);
			chunkh.Add128(w, cc);
		}
		else chunkh.Addband3(w.bf.u32[2]);
	}
}

//________ start expand bands 1+2 3 steps

struct SPB03 {// spots to first 7 clues
	BF128 v;
	uint64_t  possible_cells, all_previous_cells, active_cells;
	CBS cbs;
	uint32_t ncl;
	void Dump(uint64_t x) {
		cout << "spb03 status ncl=" << ncl << " " << x << endl;
		cout << Char54out(all_previous_cells) << " assigned" << endl;
		cout << Char54out(active_cells) << " active" << endl;
		cout << Char54out(possible_cells) << " possible" << endl;
		cout << Char64out(v.bf.u64[0]) << " 64 v" << endl;
	}
}spb_0_15[16];
uint64_t ts_47[30][300];
//____  expand 03
void TUASB12::TUAVECTH::Build_tv128() {
	for (int i = 0; i < 20; i++)tv128[i].Init();
	memset(nta128, 0, sizeof nta128);
	na128 = nablocs = 0;
	for (uint32_t i = 0; i < tuasb12.nua; i++)
		AddA(tuasb12.tua[i]);
}
void GCHK::Expand_03() {
	SPB03 *s,*sn;
	if (aigstop) return;
	//cout << iperm << "  exp03" << endl;
	//if (iperm) return;
	tuasb12.SwitchTo54Mode();// switch tua to 54 mode
	tuasb12.tv128h.Build_tv128();
	zh2b[0].InitBands12(grid0);
	TUASB12::TUAVECT& tuv128 = tuasb12.tv128h.tv128[0];
	uint64_t* twu = tuv128.t;
	s = spb_0_15;	memset(s, 0, sizeof spb_0_15[0]);
	s->active_cells = maskLSB[54].u64[0];
	s->possible_cells = twu[0];
	s->v=tuv128.v0;// initial nothing done

	//_______ start search 3 first clues
next:
	// catch and apply cell in bitfields
	int cell;
	uint64_t p = s->possible_cells;
	if (!p)	if (--s >= spb_0_15)goto next; else return;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	uint64_t ac = s->active_cells;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v &= tuv128.vc[cell];
	if (sn->ncl == 3) {// 3 cellsfirst step
		p_cpt2g[3]++;
		if (SetupExpand_46()) goto next;// dead
#ifdef HAVEKNOWN
		if (!((~pk54) & sn->all_previous_cells)) {
			cout << Char54out(sn->all_previous_cells) << " expected 3" << endl;
			Expand_46();
			aigstop = 1;
		}
		else 	Expand_46();
#else 
		Expand_46();

#endif
		if (aigstop) return;
		goto next;
	}
	// find next ua
	int ir = sn->v.getFirst128();
	if (ir < 0) return;//never
	uint64_t Ru = twu[ir] & ac;
	if (!Ru)goto next;//dead branch unlikely
	sn->possible_cells = Ru;
	s++; // switch to next spot
	goto next;
}

//____ expand 46
//#define TEST11
void TUASB12::TUAVECTH::Build_tvb128(uint32_t* ntt) {
	for (int i = 0; i < 10; i++)tvb128[i].Init();
	memset(ntb128, 0, sizeof ntb128);
	nb128 = nbblocs = 0;
	for (uint32_t i = 0; i <=20; i++) {
		register uint32_t nn = ntt[i];
		if (nn) {
			register uint64_t* t = ts_47[i];
			for (uint32_t j = 0; j<nn; j++) 
				AddB(t[j]);			
		}
	}

}
int GCHK::SetupExpand_46() {
	BF128 tvw[20];
	uint32_t lastbloc = tuasb12.tv128h.nablocs;
	tvw[0] = spb_0_15[3].v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUASB12::TUAVECT & vv= tuasb12.tv128h.tv128[i];
		BF128 v = vv.v0 , * vc=vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[tclues[ic]];
		tvw[i]=v;
	}
	// apply active on uas and sort by size
	uint32_t  ntt[21];
	memset(ntt, 0, sizeof ntt);
	{
		register uint64_t Ac = spb_0_15[3].active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t * t= tuasb12.tv128h.tv128[i].t;
			BF128 V= tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir>=0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 20)cc = 20;
					ts_47[cc][ntt[cc]++] = R;
				}
				else break;
			}
		}
	}
	tuasb12.tv128h.Build_tvb128(ntt);
	return 0;
}
void GCHK::Expand_46() {
	if (aigstop) return;
	SPB03*sl= &spb_0_15[4] ,* s=sl, * sn;
	TUASB12::TUAVECT& tuv128 = tuasb12.tv128h.tvb128[0];
	uint64_t* twu = tuv128.t;
	*s = spb_0_15[3];	// duplicate 3 for new vector
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done

	//_______ start search clues 4-6
next:	// catch and apply cell in bitfields
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)	if (--s >= sl)goto next; else return;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	uint64_t ac = s->active_cells;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v &= tuv128.vc[cell];
	if (sn->ncl == 6) {// 6 cells 
		SetupExpand_7p();
#ifdef HAVEKNOWN
		if (!((~pk54) & sn->all_previous_cells)) {
			cout << Char54out(sn->all_previous_cells) << " expected 6" << endl;
			if (mincluesb3 == 6)Expand_7_12();
			else 		Expand_7_11();
			aigstop = 1;
		}
		else {
			if (mincluesb3 == 6)Expand_7_12();
			else 		Expand_7_11();
		}
#else
		if (mincluesb3 == 6)Expand_7_12();
		else 		Expand_7_11();
#endif
		if (aigstop) return;
		goto next;
	}
	// find next ua
	int ir = sn->v.getFirst128();
	uint64_t Ru;
	if (ir < 0) {//never valid
		if (tuasb12.tv128h.nbblocs) {// more uas to check
			for (uint32_t i = 1; i <= tuasb12.tv128h.nbblocs; i++) {
				TUASB12::TUAVECT& vv = tuasb12.tv128h.tvb128[i];
				BF128 v = vv.v0, * vc = vv.vc;
				for (uint32_t ic = 3; ic < sn->ncl; ic++)
					v &= vc[tclues[ic]];
				if (v.isNotEmpty()) {
					int ir2 = v.getFirst128();
					uint64_t Ru = vv.t[ir2] & ac;
					if (!Ru)goto next;//dead branch
					sn->possible_cells = Ru;
					s++; // switch to next spot
					goto next;
				}
			}

		}
		if (zh2b[1].IsValid(tclues, sn->ncl)) {
			uint32_t i = zh2gxn.nua - 1;
			register uint64_t ua = zh2gxn.tua[i],
				cc = _popcnt64(ua),
				ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			tuasb12.tv128h.AddA(ua54);
			tuasb12.tv128h.AddB(ua54);
			Ru = ua;
		}
		else {
			cout << "bug exp 4-7 lower 7" << endl;
			aigstop = 1; return;
		}

	}
	else  Ru = twu[ir] & ac;
	if (!Ru)goto next;//dead branch 
	sn->possible_cells = Ru;
	s++; // switch to next spot
	goto next;
}

//____ expand 7_12

void TUASB12::TUAVECTH::Build_tvc128(uint32_t* ntt) {
	uint64_t tlow[50], ntlow = 0;
	for (int i = 0; i < 10; i++)tvc128[i].Init();
	memset(ntc128, 0, sizeof ntc128);
	nc128 = ncblocs = 0;
	for (uint32_t i = 0; i <= 20; i++) {
		register uint32_t nn = ntt[i];
		if (nn) {
			register uint64_t* t = ts_47[i];
			for (uint32_t j = 0; j < nn; j++) {
				register uint64_t U = t[j],
					Un = ~U;
				for (uint32_t k = 0; k < j; k++)//is it redundant
					if (t[k] == U) { U = 0; break; }
				if (!U)continue;
				for (uint64_t k = 0; k < ntlow; k++)//has it subset
					if (!(tlow[k] & Un)) { U = 0; break; }
				if (U) {
					if (nn < 8 && ntlow < 50) tlow[ntlow++] = U;
					AddC(U);
				}
			}
		}
	}

}
int GCHK::SetupExpand_7p() {
	BF128 tvw[20];
	uint32_t lastbloc = tuasb12.tv128h.nbblocs;
	tvw[0] = spb_0_15[7].v;// v for 6 clues
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUASB12::TUAVECT& vv = tuasb12.tv128h.tvb128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 3; ic < 6; ic++)
			v &= vc[tclues[ic]];
		tvw[i] = v;
	}
	// apply active on uas and sort by size
	uint32_t  ntt[21];
	memset(ntt, 0, sizeof ntt);
	{
		register uint64_t Ac = spb_0_15[7].active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = tuasb12.tv128h.tvb128[i].t;
			BF128 V = tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir >= 0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					//if (cc < 10)
						//cout << Char54out(R) << " Ri=" << ir
						//<< " cc=" << cc << " bloc " << i << endl;
					if (cc > 20)cc = 20;
					ts_47[cc][ntt[cc]++] = R;
				}
				else break;
			}
		}
	}
	tuasb12.tv128h.Build_tvc128(ntt);
	return 0;
}
void GCHK::Init7p_guas() {
	register uint64_t Ac = spb_0_15[7].active_cells;
	chunkh.ApplyClean(tclues, 6);
	// _______________extract guas2 and vectors
	BF128 t[2000], tw;// temp storage for still active
	uint32_t  tt[15][100],ntt[15], 
		nt = 0, n2 = 0,  nm = 0;
	memset(ntt, 0, sizeof ntt);

	// split guas2   in 2 x 256 vectors ( room for more )
	g2_256[0].Init();  g2_256[1].Init(); g3_256.Init();
	gm_256[0].Init();  gm_256[1].Init();
	for (uint32_t ic2 = 0; ic2 <= chunkh.ic2; ic2++) {
		CHUNK3B& w = chunkh.c2[ic2];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			register uint32_t i27 = w.tu27[ir];
			register uint64_t pat12 = w.tu54[ir] & Ac;
			tw.bf.u64[0] = pat12;
			tw.bf.u32[3] = i27;
			tw.bf.u32[2] = 1 << i27;
			register uint32_t cc = (uint32_t)_popcnt64(pat12);
			// below cc=10 load direct else sort by size
			if (cc < 10) { g2_256[0].Add(tw); n2++; continue; }
			if (cc < 20) cc -= 10;
			else if (nt < 500)cc = 9;
			else continue;
			tt[cc][ntt[cc]++] = nt;
			if (nt < 2000)t[nt++] = tw;
			else break;
		}
	}

	for (uint32_t i = 0; i < 10; i++) {
		uint32_t* tti = tt[i];
		for (uint32_t j = 0; j < ntt[i]; j++) {
			BF128 w = t[tti[j]];
			if (n2++ < 256)g2_256[0].Add(w);
			else g2_256[1].Add(w);
		}
	}
	// ___________________extract guas3
	for (uint32_t ic3 = 0; ic3 <= chunkh.ic3; ic3++) {
		CHUNK3B& w = chunkh.c3[ic3];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			uint32_t i9 = w.tu27[ir];
			register uint64_t pat12 = w.tu54[ir] & Ac;
			tw.bf.u64[0] = pat12;
			tw.bf.u32[3] = i9;
			tw.bf.u32[2] = 1 << i9;
			g3_256.Add(tw);
		}
	}
	tw.bf.u32[3] = 0;
	// _____________________________________________extract guas4
	for (uint32_t ic4 = 0; ic4 <= chunkh.ic4; ic4++) {
		CHUNK3B& w = chunkh.c4[ic4];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & Ac;
			tw.bf.u32[2] = w.tu27[ir];
			gm_256[0].Add(tw);
			nm++;
		}
	}
	// _____________________________________________extract guas5
	for (uint32_t ic5 = 0; ic5 <= chunkh.ic5; ic5++) {
		CHUNK3B& w = chunkh.c5[ic5];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & Ac;
			tw.bf.u32[2] = w.tu27[ir];
			gm_256[0].Add(tw);
			nm++;
		}
	}
	// _____________________________________________extract guasmore
	for (uint32_t icm = 0; icm <= chunkh.icmore; icm++) {
		CHUNK3B& w = chunkh.cmore[icm];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & Ac;
			tw.bf.u32[2] = w.tu27[ir];
			if (nm++ < 256)gm_256[0].Add(tw);
			else gm_256[1].Add(tw);
		}
	}
	g_256h.GetStart(n2, nm);
}

int GCHK::IsValid7p(SPB03* sn) {
	if (zh2b[1].IsValid(tclues, sn->ncl)) {
		anduab12 = ~0;
		register uint64_t ua54;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i], cc = _popcnt64(ua);
			ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			anduab12 &= ua54;
			if (cc < 21) {
				tuasb12.tv128h.AddA(ua54);
				tuasb12.tv128h.AddB(ua54);
				tuasb12.tv128h.AddC(ua54);
			}
		}
		ua_ret7p = ua54;// return last (smaller)
		return 1;
	}
	return 0;
}
inline int GCHK::GetNextCell(SPB03* s ) {
	SPB03* sn;
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)return 1;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v &= tuasb12.tv128h.tvc128[0].vc[cell];
	return 0;
}
inline void GCHK::GetNextUa(SPB03* sn) {
	register uint64_t  V;
	if ((V = sn->v.bf.u64[0])) {// next ua
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p=  tuasb12.tv128h.tvc128[0].t[ir];
	}
	else {// next ua must be here
		V = sn->v.bf.u64[1];
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p= tuasb12.tv128h.tvc128[0].t[ir + 64];
	}
}
inline void GCHK::GetNextUaAdd(SPB03* sn) {
	if (tuasb12.tv128h.nc128<=128) return ;
	// more uas to check
	for (uint32_t i = 1; i <= tuasb12.tv128h.ncblocs; i++) {
		TUASB12::TUAVECT& vv = tuasb12.tv128h.tvc128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 6; ic < sn->ncl; ic++)
			v &= vc[tclues[ic]];
		if (v.isNotEmpty()) {
			int ir2 = v.getFirst128();
			ua_ret7p = vv.t[ir2];
			return ;
		}
	}
}
inline int GCHK::GetLastAndUa(SPB03* sn,int d) {
	int aig = 0;
	register uint64_t  V, And = ~0;
	register uint32_t ir;
	if ((V = sn->v.bf.u64[0])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= tuasb12.tv128h.tvc128[0].t[ir];
		}
	}
	if (!And) return 1 ;
	if ((V = sn->v.bf.u64[1])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= tuasb12.tv128h.tvc128[0].t[ir + 64];
			if (!And) return 1;
		}
	}
	if (!And) return 1;
	if (tuasb12.tv128h.ncblocs) {
		// more uas to check
		for (uint32_t i = 1; i <= tuasb12.tv128h.ncblocs; i++) {
			TUASB12::TUAVECT& vv = tuasb12.tv128h.tvc128[i];
			BF128 v = vv.v0, * vc = vv.vc;
			for (uint32_t ic = 6; ic < sn->ncl; ic++)
				v &= vc[tclues[ic]];
			if (v.isNotEmpty()) {
				aig = 1;
				int ir2;
				while ((ir2 = v.getFirst128()) >= 0) {
					v.clearBit(ir2);
					And &= vv.t[ir2];
				}
			}
		}
	}
	if (aig) ua_ret7p = And;
	return aig;
}

void GCHK::Expand_7_11() {
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	TUASB12::TUAVECT& tuv128 = tuasb12.tv128h.tvc128[0];
	uint64_t* twu = tuv128.t;
	if (tuasb12.tv128h.nc128 < 128)tuasb12.tv128h.nc128 = 128;// force adds outside 
	*s = spb_0_15[7];	// duplicate 6 for new vector
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done
	Init7p_guas();// do init 7 clues

	//_______ start search clues 7_11
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))
		if (--s >= sl)goto next;	else return;
	sn = s + 1;
#ifdef HAVEKNOWN
	if (!((~pk54) & sn->all_previous_cells)) {
		cout << Char54out(sn->all_previous_cells) << " on the path 7_11 ncl=" << sn->ncl << endl;
	}
#endif
	ua_ret7p = 0;
	if (sn->ncl == 11) {// 11 cells 
		if (sn->cbs.IsFilt11()) goto next;
		GetNextUaAdd(sn);// check adds 
		if(ua_ret7p)goto next;
		clean_valid_done = 0;
		p_cpt2g[5]++;
		g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
		GoAfterExpand(sn);
		goto next;
	}
	if (sn->ncl == 10) {// last step 
		if (GetLastAndUa(sn)) // last not empty
			if (!ua_ret7p) goto next;
	}
	else {
		if (sn->v.isNotEmpty())GetNextUa(sn);// first 128
		else GetNextUaAdd(sn);
	}
	if (!ua_ret7p) {// add ua or go below
		if (IsValid7p(sn)) {// got uas to use
			if (sn->ncl == 10)ua_ret7p = anduab12;
		}
		else {// valid below
			g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
			clean_valid_done = 1;
			GoAfterExpand(sn);//try direct
			GoBelow(sn);
			goto next;
		}
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) { // switch to next spot
		g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
		s++;
	}
	goto next;
}
void GCHK::Expand_7_12() {
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	if (aigstop) return;
	TUASB12::TUAVECT& tuv128 = tuasb12.tv128h.tvc128[0];
	uint64_t* twu = tuv128.t;
	if (tuasb12.tv128h.nc128 < 128)tuasb12.tv128h.nc128 = 128;// force adds outside 
	*s = spb_0_15[7];	// duplicate 6 for new vector
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done
	Init7p_guas();// do init 7 clues

	//_______ start search clues 7_12
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))	if (--s >= sl)goto next; else return;
	sn = s + 1;
#ifdef HAVEKNOWN
	if (!((~pk54) & sn->all_previous_cells)) {
		cout << Char54out(sn->all_previous_cells) << " on the path 7_12 ncl=" << sn->ncl << endl;
	}
#endif
	ua_ret7p = 0;
	if (sn->ncl == 12) {// 12 cells 
		if (sn->cbs.IsFilt12()) goto next;
		p_cpt2g[5]++;
		GetNextUaAdd(sn);// check adds
		if (ua_ret7p)goto next;
		g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
		clean_valid_done = 0;
		GoAfterExpand(sn);
		if (aigstop) {
			cout << Char54out(sn->all_previous_cells) << " stop 7-12 tuasb12.tv128h.ncblocs="
				<< tuasb12.tv128h.ncblocs  << endl;
			tuasb12.tv128h.tvc128[1].Dump();
			tuasb12.tv128h.tvc128[1].Dumpv();
			return;
		}
		goto next;
	}
	ua_ret7p = 0;
	if (sn->ncl == 11) {// last step 
		if (GetLastAndUa(sn)) {// first 128
			if (!ua_ret7p) goto next;
		}
	}
	else {
		if (sn->v.isNotEmpty())		GetNextUa(sn);// first 128
		else GetNextUaAdd(sn);
	}
	if (!ua_ret7p) {
#ifdef HAVEKNOWN
		if (!((~pk54) & sn->all_previous_cells)) {
			cout << Char54out(sn->all_previous_cells) << " after on the path 7_12 ncl=" << sn->ncl << endl;
		}
#endif
		if (sn->ncl < 11) {
			if (!IsValid7p(sn)) {// valid below
				g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
				clean_valid_done = 1;
				GoAfterExpand(sn);//try direct
				GoBelow(sn);		goto next;
			}
		}
		else {// 11 clues try direct first
			if (sn->cbs.IsFilt11())goto next;
			g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
			clean_valid_done = 0;
			GoAfterExpand(sn);//try direct
			if (clean_valid_done == 2)  // bach not validb12
				ua_ret7p = anduab12;
			else { GoBelow(sn); goto next; }
		}
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) {
		g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
		s++; // switch to next spot
	}
	goto next;
}

int GCHK::IsValidB12(uint32_t ncl) {
	//mbisvalid.ntm = 0;
	p_cpt2g[29]++;
#ifdef DEBUG5
	int locdiag = 0;
	if (p_cpt2g[5] == DEBUG5) {
		cout << "entry isvalid" << endl;
		locdiag = 1;
	}
#endif
	anduab12 = ~0;
	if (zh2b[1].IsValid(tclues, ncl)) {
		p_cpt2g[30]++;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i],
				cc = _popcnt64(ua),
				ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
#ifdef DEBUG5
			if (locdiag) 		
			cout << Char54out(ua54) << " i=" << i << endl;

#endif
			if(ntuaddb12<50)tuaddb12[ntuaddb12++]=ua54;// used in add mode
			anduab12 &= ua54;
			if (cc < 20) {
				tuasb12.tv128h.AddC(ua54);
				if (cc < 20) {
					tuasb12.tv128h.AddB(ua54);
					if (cc < 19)tuasb12.tv128h.AddA(ua54);
				}
			}
		}
	}
	//nclf = ncl;
	return zh2gxn.nua; 
}

#define stack1_54 07007007007007007
int tstack27[3] = { 07007007 ,070070070,0700700700 };
struct CRITB3 {
	uint32_t minix[4],// triplet bf1 bf2 bf3  
		critbf, pairs27, mincount,
		critstack, stackmin[3], stackf[3], stackb12[3],
		assigned, active,
		ncl, nb3, nmiss;
	inline void Init(int ncb12, uint64_t bf, CBS& cbsx) {
		memset(this, 0, sizeof(*this));
		ncl = ncb12, nb3 = 18 - ncl;
		memcpy(stackb12, cbsx.s, sizeof stackb12);
		//minb3 = minb3e;
	}
	inline int  CanNotAddi27(int i27) {
		int bit27 = 1 << i27;
		if (assigned & bit27)return 1;// critical stack common clue assigned
		int imini = i27 / 3, bitmini = 1 << imini, mask = 7 << (3 * imini);
		if (minix[0] & bitmini) {// was a triplet
			minix[0] ^= bitmini;
			minix[1] |= bitmini;
			pairs27 |= bit27;
			critbf ^= bit27;
			return 0;
		}
		if (minix[1] & bitmini) {// was a one pair
			minix[1] ^= bitmini;
			minix[2] |= bitmini;
			pairs27 |= bit27;
			critbf |= mask;
			return 0;
		}

		if (minix[3] & bitmini) {// can not be
			cout << "critb3 can not add to 3 pairs bug" << endl;
			gchk.aigstop = 1;
			return 1;
		}
		// now true add 0->1 or 2_>3
		int  stack = C_stack[i27], bitstack = 1 << stack;
		if (critstack & bitstack)return 1;// not possible

		if (minix[2] & bitmini) {// was 2 pairs
			minix[2] ^= bitmini;
			minix[3] |= bitmini;
			pairs27 |= bit27;
		}
		else {// first pair in the mini row no triplet
			minix[1] |= bitmini;
			pairs27 |= bit27;
			critbf |= (mask ^ bit27);
		}
		stackf[stack]++;
		nmiss--;
		if (!nmiss) {
			critstack = 7;
			active &= critbf;// no more outfield
		}
		else if (stackf[stack] == nb3) {
			critstack |= bitstack;
			active &= (~tstack27[stack]) | critbf;
		}
		return 0;
	}
	inline int  CanNotAddi9(int i9) {
		//can not have a pair active
		int  stack = i9 % 3, bitstack = 1 << stack;
		if (critstack & bitstack)return 1;// not possible
		int mask = 7 << (3 * i9);
		critbf |= mask;
		minix[0] |= 1 << i9;
		stackf[stack]++;
		nmiss--;
		if (!nmiss) {
			critstack = 7;
			active &= critbf;// no more outfield
		}
		else if (stackf[stack] == nb3) {
			critstack |= bitstack;
			active &= (~tstack27[stack]) | critbf;
		}
		return 0;
	}

	inline int Addone(uint32_t i27) {// back 1 if not possible
		int bit27 = 1 << i27, stack = C_stack[i27];
		int imini = i27 / 3, bitmini = 1 << imini, bitstack = 1 << stack;
		assigned |= bit27;
		if (!(bit27 & critbf)) {// clue added outfield
			if (critstack & bitstack)return 1;// not possible
			nmiss--;
			if (!nmiss) {
				critstack = 7;
				active &= critbf;// no more outfield
			}
			else {
				stackf[stack]++;
				if (stackf[stack] == nb3)critstack |= bitstack;
			}
			return 0;;
		}
		// now add in field within mincount
		if (minix[3] & bitmini) {// 2 clues expected
			critbf ^= bit27;minix[3] ^= bitmini;return 0;
		}
		if ((minix[2] & bitmini) && (pairs27 & bit27)) {
			// 2 clues if not common clue one more clue
			if (critstack & bitstack)return 1;// not possible
			critbf ^= bit27;	minix[2] ^= bitmini;
			nmiss--;
			if (!nmiss) {
				critstack = 7;
				active &= critbf;// no more outfield
			}
			else {
				stackf[stack]++;
				if (stackf[stack] == nb3)critstack |= bitstack;
			}
			return 0;
		}
		register int mask = ~(7 << (3 * imini));// clear minirow
		critbf &= mask;
		return 0;
	}
	inline int AddAssign(uint32_t bfa) {// back 1 if not possible
		if (assigned & bfa)return 1; //should never be
		active &= ~bfa; // minimum is to kill new assign
		register int i27, X = bfa;
		while (bitscanforward(i27, X)) {
			X ^= 1 << i27;// clear bit
			if (Addone(i27))return 1;
		}
		return 0;
	}

	void AssignBf2(int bf) {// all or part of the bf2
		minix[2] &= ~bf;// clean the bf2 
		int imini;
		while (bitscanforward(imini, bf)) {
			int  mask = 7 << (3 * imini), m27 = pairs27 & mask;
			bf ^= 1 << imini;// clean the bf
			critbf &= ~mask;// clean the mini row as crit field
			active &= ~mask;// no more clue in this minirow
			pairs27 ^= m27;// clean the 27 pairs to assign
			assigned |= mask ^ m27;// assign the third cell of the mini row
		}
	}
	inline void AssignCritical() {
		active = critbf;
		if (critstack == 7) {
			if (minix[2])AssignBf2(minix[2]);
			return;
		}
		else if (critstack) {
			for (int i = 0, bit = 1, mask = 07007007; i < 3; i++, bit <<= 1, mask <<= 3)
				if (!(critstack & bit)) active |= mask;
		}
		else active = BIT_SET_27;
		if ((!critstack) || (!minix[2])) return;
		for (int i = 0, bit = 1, mask = 0111; i < 3; i++, bit <<= 1, mask <<= 1) {
			if (critstack & bit) {
				active &= ~(stack1_54 << 3 * i); //  (stack154 << 3*i)active must be in critbf for this; stack
				if (minix[2] & mask)		AssignBf2(minix[2] & mask);

			}
		}
		active |= critbf;// be sure to keep it

	}

	inline void MinVmini(int imini, int vmini) {
		uint32_t shift = 3 * imini, bit9 = 1 << imini, cmin = 1,
			mask = 7 << shift, vminishift = vmini << shift;
		if (vmini == 8) {
			minix[0] |= bit9;// mini triplet
			critbf |= mask;
		}
		else {
			pairs27 |= vminishift;
			uint32_t cc = _popcnt32(vmini);
			minix[cc] |= bit9;
			if (cc > 1)critbf |= mask;
			else critbf |= (mask ^ vminishift);
			if (cc == 3) 	cmin = 2;
		}
		mincount += cmin;
		stackmin[imini % 3] += cmin;
	}
	inline int Stackx() {
		register uint32_t x = nb3;
		if (stackf[0] > x)return 1;
		if (stackf[1] > x)return 1;
		if (stackf[2] > x)return 1;
		if (stackf[0] == x)critstack |= 1;
		if (stackf[1] == x)critstack |= 2;
		if (stackf[2] == x)critstack |= 4;
		return 0;
	}
	inline void SetStackf(uint64_t bf) {
		memcpy(stackf, stackmin, sizeof stackf);
		stackf[0] += stackb12[0];
		stackf[1] += stackb12[1];
		stackf[2] += stackb12[2];
		nmiss = nb3 - mincount;
		if (!nmiss) critstack = 7;
	}
	inline int BadPerm() {
		if(ncl!=18) return 0;// skip if 17
		for (int i = 0; i < 3; i++)	if (stackf[i] > nb3) return 1;
		if (nb3 == 6) return 0; // finished for 12 clues
		int nequal = 0;
		for (int i = 0; i < 3; i++)	if (stackf[i] == nb3) nequal++;
		if (!nequal) return 0;// keep it
		return 0;
		return (gchk.iperm > 2);
	}
	void Status(const char* lib) {
		cout << lib << "critical Status mincount =" << mincount << " nmiss=" << nmiss
			<< " critstack=" << critstack << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char27out(assigned) << " assigned" << endl;
		cout << Char27out(active) << " active" << endl;
		if (minix[1])cout << Char9out(minix[1]) << "     minis bf1" << endl;
		if (minix[2])cout << Char9out(minix[2]) << "     minis bf2" << endl;
		if (minix[3])cout << Char9out(minix[3]) << "     minis bf3" << endl;
		if (minix[0])cout << Char9out(minix[0]) << " mini triplets" << endl << endl;
		cout << "b12 stacks " << stackb12[0] << stackb12[1] << stackb12[2] << endl;
		cout << "minis stacks " << stackmin[0] << stackmin[1] << stackmin[2] << endl;
		cout << "final stacks " << stackf[0] << stackf[1] << stackf[2] << endl;
	}

}scritb3;

void GCHK::GoBelow(SPB03* sn) {
#ifdef HAVEKNOWN
	if (!((~pk54) & sn->all_previous_cells)) {
		cout << Char54out(sn->all_previous_cells)<< " go below" << endl;
		//g2_256[0].Dump(); g2_256[1].Dump();
		//g_256h.Dumpvcl(sn->ncl);
	}
#endif
	rclean_valid_done=  clean_valid_done;
	g2_64.nv = 0; g3_64.nv = 0; gm_64.nv = 0;
	int nclues = sn->ncl;
	{// check limit per band  11 clues 747 567 657 12 clues 666
		myacadd = sn->active_cells;
		if (nclues == 11) {// 11 clues 747 567 657 12 clues 666
			if (sn->cbs.b[1] < 4 || sn->cbs.b[1]>6) return;
			if (sn->cbs.b[0] == 7)myacadd = 0;
			else if (sn->cbs.b[1] == 6)myacadd &= BIT_SET_27;//add in b1
			else myacadd &= ~BIT_SET_27;//add in b2
		}
		else if (nclues == 10) {// add 1 or add 2 depending on limb12
			if (sn->cbs.b[1] > 6)myacadd = 0; // 288 378
			else {// 73 46 64 55
				if (sn->cbs.b[1] > 5)myacadd &= BIT_SET_27;//add in b1
				else if (sn->cbs.b[1] < 4)myacadd &= ~BIT_SET_27;//add in b2
			}
		}
	}
	if (!myacadd) return;
	{//______________  build tadd from active cells
		ntadd = 0;
		int cell;
		register uint64_t V = myacadd;// still valid c2
		while (bitscanforward64(cell, V)) {
			V ^= (uint64_t)1 << cell;
			tadd[ntadd++] = cell;
		}
	}
#ifdef HAVEKNOWN
	if (!((~pk54) & sn->all_previous_cells)) {
		cout  << " ntadd="<< ntadd <<"  g2_64.nv =" << g2_64.nv << endl;
	}
#endif
	ExpandAddB1B2(sn);
}
void GCHK::ExpandAddB1B2(SPB03* sne ) {// add up to n cells
	mynclues = sne->ncl;
	uint32_t lim = 18 - mincluesb3;
	SPB03* sl = sne, * s = sl, * sn;
	s->possible_cells = 0;// here index
next:	//____________ here start the search ad after valid
	if (s->possible_cells >= ntadd)
		if (--s >= sl)goto next; else return;
	{
		register int cell = tadd[s->possible_cells++];
		register uint64_t bit = (uint64_t)1 << cell;
		tclues[s->ncl] = cell;
		s->active_cells ^= bit;
		sn = s + 1; *sn = *s; sn->ncl++;
		sn->cbs.Add(cell);
		sn->all_previous_cells |= bit;
		sn->v &= tuasb12.tv128h.tvc128[0].vc[cell];
		g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
#ifdef HAVEKNOWN
		if (!((~pk54) & sn->all_previous_cells)) {
			cout << Char54out(sn->all_previous_cells) 
				<< " on the path below ncl=" << sn->ncl << "  g2_64.nv =" << g2_64.nv << endl;
			//g_256h.Dumpvcl(sn->ncl);
		}
#endif
		if (sn->ncl == 12 && sn->cbs.IsFilt12()) goto next;
		else if (sn->ncl == 11 && sn->cbs.IsFilt11()) goto next;
		if (!rclean_valid_done) {// must check fresh uas b3 (redundancy)
			clean_valid_done = 0;
			if (ntuaddb12) {
				register uint64_t bf = sn->all_previous_cells;
				for (uint32_t i = 0; i < ntuaddb12; i++) {
					if (!(tuaddb12[i] & bf)) {// skip a not valid b12
						if (sn->ncl >= lim) 	goto next;
						s++; // switch to next spot (next cell)
						goto next;
					}
				}
			}
		}
		GoAfterExpand(sn,ntuaddb12);//keep old uas add
		if (aigstop) return;
		if (sn->ncl >= lim) 	goto next;
		s++; // switch to next spot (next cell)
		goto next;
	}
}

//_______________________
void GCHK::GoAfterExpand(SPB03* sn, uint32_t nadd){
	if (aigstop) return;
#ifdef HAVEKNOWN
	if (pk54 == sn->all_previous_cells) {
		cout << Char54out(sn->all_previous_cells)<< " expected GoAfterExpand nadd="<<nadd 
			<< " g2_64.nv=" << g2_64.nv << " g3_64.nv=" << g3_64.nv << endl;
		if (g2_64.nv)g2_64.Dump();
	}
#endif

	ntuaddb12 = nadd;// used only if add mode
	nclgo = sn->ncl;
	ncluesb3 = 18 - nclgo;
	myb12f = sn->all_previous_cells;
	g_256h.GetActive(nclgo, myb12f);
	if ((int)g_256h.ntg2 > 2 * ncluesb3)return;
	register uint64_t Bf =sn->all_previous_cells;
	register uint32_t  vmini, Mg2 = g_256h.g2,	Mg3 = g_256h.g3;
	scritb3.Init(nclgo, sn->all_previous_cells, sn->cbs);
	for (int imini = 0; imini < 9; imini++, Mg2 >>= 3, Mg3 >>= 1) {
		if (!(vmini = Mg2 & 7))	if (Mg3 & 1)vmini = 8;
		if (vmini)scritb3.MinVmini(imini, vmini);
	}
	if ((int)scritb3.mincount > ncluesb3) return;
	scritb3.SetStackf(Bf);// stackf; nmisss; critstack if nmiss=0
	if (scritb3.Stackx())return;
	scritb3.AssignCritical();
#ifdef HAVEKNOWN
	if (pk54 == sn->all_previous_cells) {
		scritb3.Status("aaaa");
	}
#endif

	ntw3 = 0; // temporary before  split in/out
	g_256h.ApplyMore(&tclues[6], nclgo - 6,
		tw3, ntw3, scritb3.assigned);
	register uint32_t  critbf = scritb3.critbf,
		nmiss = scritb3.nmiss, n1 = 0, n2 = 0;
	if (!nmiss) {// first out is "dead"
		for (uint32_t i = 0; i < ntw3; i++)
			if (!(tw3[i] & critbf)) return;
	}
	else if (nmiss == 1) {//use "and" as out
		register uint32_t andout = scritb3.active & (~critbf);
		for (uint32_t i = 0; i < ntw3; i++) {
			register uint32_t U = tw3[i];
			if (U & critbf) tw3[n1++] = U;
			else { andout &= U; n2 = 1; }
		}
		if (n2) {
			if (!andout)return;// no possibility to add
			tw3[n1++] = andout; // dummy ua hitting all "out"
		}
		ntw3 = n1;
	}
	ntaddgob3 = 0;
#ifdef HAVEKNOWN
	if (pk54 == sn->all_previous_cells) cout << "call buildb2" << endl;
#endif
	if (BuildB2Table()) return;
	int nass = _popcnt32(scritb3.assigned);
	if (nass == scritb3.nb3) {// already ok to go
		if (ntw3_2)return;
		if (scritb3.BadPerm()) return;
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
		}
		if (!IsValidB3(scritb3.assigned)) Out17(scritb3.assigned);
		return;
	}
#ifdef HAVEKNOWN
	if (pk54 == myb12f) cout << "call BuildFinalTable" << endl;
#endif
	BuildFinalTable();
#ifdef HAVEKNOWN
	//if (pk54 == sn->all_previous_cells) aigstop = 1;
#endif
}


//_____ process a chunk of potential valid bands 1+2
uint32_t tbuilexpand3[12][500];// 6+6


int GCHK::BuildB2Table() {// extract and reorder still possible guas
	uint32_t is1 = 0, ntt[6], i27;
	memset(ntt, 0, sizeof ntt);
	register uint32_t AC = scritb3.active,
		F = scritb3.assigned;
	{
		register uint32_t V = scritb3.pairs27;
		while (bitscanforward64(i27, V)) {
			V ^= 1 << i27;
			tbuilexpand3[2][ntt[2]++] = tg2[i27].pat;
		}
		V = scritb3.minix[0];// triplets
		if (V) {// add guas3 if any (usually 0)
			for (int i = 0; i < 9; i++)if (V & (1 << i))
				tbuilexpand3[3][ntt[3]++] = 7 << (3 * i);
		}

	}

	register uint32_t * t2 = tbuilexpand3[2];
	for (uint32_t i = 0; i < ntw3; i++) {
		register uint32_t U = tw3[i];
		if (!(U&F)) {
			U &= AC;
			int ccu = _popcnt32(U);
			if (!ccu) return 1;
			if (ccu == 1) { is1 |= U; continue; }
			{// clear redundancy
				register uint32_t nu = ~U, x = 0;
				for (uint32_t j = 0; j < ntt[2]; j++)
					if (!(nu&t2[j])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			tbuilexpand3[ccu][ntt[ccu]++] = U;
		}
	}
	// add small band 3 residual uas
	for (uint32_t i = 0; i < myband3.nua; i++) {
		register uint32_t U = myband3.tua[i];
		if (!(U&F)) {
			U &= AC;
			int ccu = _popcnt32(U);
			if (!ccu) return 1;
			if (ccu == 1) { is1 |= U; continue; }
			{// clear redundancy
				register uint32_t nu = ~U, x = 0;
				for (uint32_t j = 0; j < ntt[2]; j++)
					if (!(nu&t2[j])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			tbuilexpand3[ccu][ntt[ccu]++] = U;
		}
	}


	for (uint32_t i = 0; i < ntaddgob3; i++) {
		register uint32_t U = taddgob3[i];
		if (!(U&F)) {
			U &= AC;
			int ccu = _popcnt32(U);
			if (!ccu) return 1;
			if (ccu == 1) { is1 |= U; continue; }
			{// clear redundancy
				register uint32_t nu = ~U, x = 0;
				for (uint32_t j = 0; j < ntt[2]; j++)
					if (!(nu&t2[j])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			tbuilexpand3[ccu][ntt[ccu]++] = U;
		}
	}
	ntw3_2 = is1_tw3=0;
	if (is1) {// now singles to assign
		if (scritb3.AddAssign(is1)) 	return 2;		
		AC = scritb3.active;	F = scritb3.assigned;
		int cc = _popcnt32(F);		if (cc > ncluesb3)return 1;
		// after added clues, revise the others
		uint32_t  ntt2[6];//tt2[i][500]is wubuf.tbuilexpand3[i+5][500]
		memset(ntt2, 0, sizeof ntt2);
		for (int i = 2; i < 6; i++)
			for (uint32_t j = 0; j < ntt[i]; j++) {
				register uint32_t U = tbuilexpand3[i][j];
				if (!(U&F)) {
					U &= AC;
					register int cc = _popcnt32(U);
					if (!cc) return 3; // dead branch
					if (cc > 5) cc = 5;
					tbuilexpand3[cc + 6][ntt2[cc]++] = U;
				}
			}
		if (ntt2[1]) {// send back the bit field status
			for (uint32_t j = 0; j < ntt2[1]; j++) 
				is1_tw3 |= tbuilexpand3[7][j];			
		}
		for (int i = 1; i < 6; i++)
			for (uint32_t j = 0; j < ntt2[i]; j++) {
				// clear redundancy
				register uint32_t U = tbuilexpand3[i + 6][j],
					nu = ~U, x = 0;
				for (uint32_t i = 0; i < ntw3_2; i++)
					if (!(nu&tw3_2[i])) { x = 1; break; }
				if (x) continue;
				tw3_2[ntw3_2++] = U;
			}
	}
	else {// apply directly stored uas
		for (int i = 2; i < 6; i++)
			for (uint32_t j = 0; j < ntt[i]; j++) {
				// clear redundancy
				register uint32_t U = tbuilexpand3[i][j],
					nu = ~U, x = 0;
				for (uint32_t k = 0; k < ntw3_2; k++)
					if (!(nu & tw3_2[k])) { x = 1; break; }
				if (x) continue;
				tw3_2[ntw3_2++] = U;
			}
	}

	return 0;
}
void GCHK::BuildFinalTable() {

	if (is1_tw3) {// still some singles to assign
		if (scritb3.AddAssign(is1_tw3)) {
			//cout << Char54out(is1_tw3) << "is1_tw3 in excess" << endl;
			return;
		}// then pack and filter the rest
		register uint32_t AC = scritb3.active, F = scritb3.assigned,
			n=0;
		for (uint32_t i = 0; i < ntw3_2; i++) {
			register uint32_t U = tw3_2[i];
			if (!(U & F)) {
				if (!(U &= AC))	return;
				tw3_2[n++] = U;
			}
		}
		ntw3_2=n;

	}
	int nass = _popcnt32(scritb3.assigned),ntoass= scritb3.nb3-nass;
#ifdef HAVEKNOWN
	if (pk54 == myb12f) {
		cout << "ntoass ="<<ntoass << endl;
	}
#endif
	 ntaddgob3 = 0;
	 if( (!ntoass) && ntw3_2)return;
	if (!ntoass) {// enough clues assigned in b3
		if (ntw3_2) return;
		B3FinalNoExpand();//must be a valid b12 and valid "all" 
		return;
	}
	if (!ntw3_2) {//could be a valid 17 or not valid below 17
		if (!clean_valid_done) {// must be a valid b12 anyway to go
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
		}
		if (!IsValidB3(scritb3.assigned)) {// a valid 17
			Out17(scritb3.assigned);
			return;
		}
		if (stopexpandb3) return;
		// now ntw3_2>0
	}
	if (ntoass == 1) {
		//add can be any of active scritb3 hitting all uas
		register uint32_t wadd = scritb3.active;
		for (uint32_t i = 0; i < ntw3_2; i++)wadd &= tw3_2[i];
		if (!wadd) return;
		if (!clean_valid_done) {// must be a valid b12 anyway to go
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
		}
		register int cell,
			Ass = scritb3.assigned;
		while (bitscanforward(cell, wadd)) {
			CRITB3 critb3 = scritb3;
			register uint32_t bit = 1 << cell;
			wadd ^= bit; //clear bit
			if (critb3.Addone(cell)) 	continue;
			if (critb3.BadPerm()) 	continue;
			if (IsValidB3(Ass | bit)) {
				if (stopexpandb3) return;
				wadd &= anduab3;
			}
			else 	Out17(Ass | bit);
		}
		return;
	}

	if (ntoass == 2) {
		uint32_t cell1, first = tw3_2 [0];// must exist checked upstream
		while (bitscanforward(cell1, first)) {
			CRITB3 critb3 = scritb3;
			register uint32_t bit = 1 << cell1;
			first ^= bit; //clear bit
			if (critb3.Addone(cell1)) 	continue;
			//register uint32_t bf1 = scritb3.assigned | bit,cell;
			//add can be any of active scritb3 hitting all uas
			//register uint32_t wadd = scritb3.active;
			register uint32_t wadd = critb3.active;
			register uint32_t bf1 = critb3.assigned ,cell;
			for (uint32_t i = 1; i < ntw3_2; i++) {
				register uint32_t U= tw3_2[i];
				if(!(U&bf1))	wadd &= U;
			}
			if (!wadd) continue;
			if (!clean_valid_done) {// ust be a valid b12 anyway to go
				clean_valid_done = 1;
				if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
			}
			while (bitscanforward(cell, wadd)) {
				CRITB3 critb3_2=critb3;
				register uint32_t bit = 1 << cell;
				wadd ^= bit; //clear bit
				if (critb3_2.Addone(cell)) 	continue;
				if (critb3_2.BadPerm()) 	continue;
				if (IsValidB3(bf1 | bit)) {
					if (stopexpandb3) return;
					wadd &= anduab3;
				}
				else 	Out17(bf1 | bit);
			}
		}
		return;
	}
	// now 3 and more clues to add
#ifdef HAVEKNOWN
	if (pk54 == myb12f) {
		cout << "ntw3_2 =" << ntw3_2 << endl;
	}
#endif
	if (ntw3_2 < 10) {
		ExpandB3Direct(ntoass);
		return;
	}	
	b3direct.Init();
	for (uint32_t i = 0; i < ntw3_2; i++) {
		b3direct.Add(tw3_2[i]);
	}
	ExpandB3Vect(ntoass);
}

void GCHK::B3FinalNoExpand() {
	if (scritb3.BadPerm()) return;
	if (!clean_valid_done) {
		clean_valid_done = 1;
		if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
	}
	if(IsValidB3(scritb3.assigned)) return;
	//_________________________   valid 18 
	Out17(scritb3.assigned);
}


//============ clues in band 3 (no more clue in bands 1+2)

#define LIMADD 17



void GCHK::ExpandB3Direct(int ntoass) {
	uint64_t limspot = ntoass - 1, limm1 = limspot - 1;
	struct SPB {
		CRITB3 critb3;
		uint32_t  possible_cells, indtw3;// , x;
	}spb[12];
	register SPB *s, *sn;
	register uint64_t ispot;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->critb3 = scritb3;
	s->indtw3 = 0;// initial first ua
	s->possible_cells = tw3_2[0];
next:
	ispot = s - spb;
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = s->possible_cells;
		if (!p)goto  back;
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		s->possible_cells ^= bit;
		s->critb3.active ^= bit;
		sn = s + 1; *sn = *s;
		if (sn->critb3.Addone(cell )) 	goto next;		
	}
	if (ispot < limm1) {// here max 16 clues never valid b3
		//if (locdiag) cout << " not last step s->indtw3="<< s->indtw3 << endl;
		register uint32_t F = sn->critb3.assigned;
		for (uint32_t i = s->indtw3 + 1; i < ntw3_2; i++) {
			register uint32_t U = tw3_2[i];
			if (!(U&F)) {
				U &= s->critb3.active;
				if (!U)goto next;//dead branch
				sn->possible_cells = U;
				sn->indtw3 = i;
				s++; // switch to next spot
				goto next;
			}
		}
		// no ua available must check( in 18 mode can not be valid)
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) {clean_valid_done = 2; return; }
		}
		if (IsValidB3(sn->critb3.assigned)) {
			if (stopexpandb3) return;
			s->possible_cells &= anduab3;// same spot hitting  new ua
			sn->possible_cells = anduab3;
			s = sn;
			goto next;
		}
		cout << "bug seen valid 16 or less [5] " << endl;
		aigstop = 1;
		return;
	}
	// this is the last step must hit all pending uas
	{ // find next cells hitting all uas
		int aig = 1;
		register uint32_t andw = sn->critb3.active;
		register uint32_t F = sn->critb3.assigned;
		{// and still valid   uas
			for (uint32_t i = s->indtw3 + 1; i < ntw3_2; i++) {
				register uint32_t U = tw3_2[i];
				if (!(U&F)) { // not hit ua
					if (!(andw &= U))goto next;//dead branch
					aig = 0;
				}
			}
		}
		// no more ua or "and all uas" not empty
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
		}
		if (aig) {// no ua could be a 17 valid
			if (IsValidB3(F)) {
				if (stopexpandb3) return;
				andw &= anduab3;
			}
			else { Out17(F);	goto next; }// this is a valid 17
		}
		register int cell;
		while (bitscanforward(cell, andw)) {
			CRITB3 critb3 = sn->critb3;
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
			if (critb3.Addone(cell)) 	continue;
			if (critb3.BadPerm()) 	continue;
			if (IsValidB3(F | bit)) {
				if (stopexpandb3) return;
				andw &= anduab3;
			}
			else 	Out17(F | bit);
		}
		goto next;
	}
	goto next;// never here
back:
	if (--s >= spb)goto next;
}
void GCHK::ExpandB3Vect(int ntoass) {
#ifdef HAVEKNOWN
	if (pk54 == myb12f) {		
		cout << "entry ExpandB3Vect"  << endl;
		//b3direct.Debug(1);
	}
#endif
	uint64_t limspot = (uint64_t)ntoass - 1, limm1 = limspot - 1;
	struct SPB {
		CRITB3 critb3;
		uint64_t vect;
		uint32_t  possible_cells, x;
	}spb[12];
	register SPB *s, *sn;
	register uint64_t ispot;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->critb3 = scritb3;
	s->vect = b3direct.v0;// initial vectsno ua hit here
	s->possible_cells = b3direct.tua[0];
	ntw3_2 = 0;// now add of this expand
next:
	ispot = s - spb;
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = s->possible_cells;
		if (!p)goto back;
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		s->possible_cells ^= bit;
		s->critb3.active ^= bit;
		sn = s + 1; *sn = *s;
		if (sn->critb3.Addone(cell)) 	goto next;
		
		sn->vect &= b3direct.vc[cell];
	}
	if (ispot < limm1) {// here max 16 clues never valid b3
		//if (locdiag) cout << " not last step s->indtw3="<< s->indtw3 << endl;
		register uint64_t V = sn->vect;
		register int ir;
		if (V) {// most often
			bitscanforward64(ir, V);//ir ua to load
			register uint32_t U = b3direct.tua[ir] & s->critb3.active;
			if (!U)goto next;//dead branch
			sn->possible_cells = U;
			s++; // switch to next spot
			goto next;
		}
		if (ntw3_2) {// see if one fresh ua can be used
			register uint32_t Bf = sn->critb3.assigned;
			for (uint32_t ia = 0; ia < ntw3_2; ia++) {
				register uint32_t U = tw3_2[ia];
				if (!(U & Bf)) {
					U &= s->critb3.active;
					if (!U)goto next;//dead branch
					sn->possible_cells = U;
					s++; // switch to next spot
					goto next;
				}
			}
		}
		// no ua available must check( in 18 mode can not be valid)
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
		}
		if (IsValidB3(sn->critb3.assigned)) {
			if (stopexpandb3) return;
			s->possible_cells &= anduab3;// same spot hitting  new ua
			sn->possible_cells = anduab3;
			s = sn;
			goto next;
		}
		cout << "bug seen valid 16 or less [5] " << endl;
		aigstop = 1;
		return;
	}
	// this is the last step must hit all pending uas
	{ // find next cells hitting all uas
		int aig = 1;
		register uint32_t andw = sn->critb3.active;
		register uint32_t F = sn->critb3.assigned;
		if (ntw3_2) {// likely quicker 0
			for (uint32_t ia = 0; ia < ntw3_2; ia++) {
				register uint32_t U = tw3_2[ia];
				if (!(U & F)) {
					if (!(andw &= U)) goto next;
					aig = 0;
				}
			}
		}
		{// and still valid   uas
			register int ir;
			register uint64_t vect = sn->vect;
			while (bitscanforward64(ir, vect)) {//ir ua to load
				vect ^= (uint64_t)1 << ir;
				if (!(andw &= b3direct.tua[ir])) goto next;
				aig = 0;
			}
		}

		// no more ua or "and all uas" not empty
		if (!clean_valid_done) {
			clean_valid_done = 1;
			//cout << " checkvalid " << endl;
			if (IsValidB12(nclgo))		{ clean_valid_done = 2; return; }
		}
		if (aig) {// no ua could be a 17 valid
			if (IsValidB3(F)) {
				if (stopexpandb3) return;
				andw &= anduab3;
			}
			else { Out17(F);	goto next; }// this is a valid 17
		}
		register int cell;
		while (bitscanforward(cell, andw)) {
			CRITB3 critb3 = sn->critb3;
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
			if (critb3.Addone(cell)) 	continue;
			if (critb3.BadPerm()) 	continue;
			if (IsValidB3(F | bit)) {
				if (stopexpandb3) return;
				andw &= anduab3;
			}
			else 	Out17(F | bit);
		}
		goto next;
	}
	goto next;// never here
back:
	if (--s >= spb)goto next;
}


uint32_t GCHK::IsValidB3(uint32_t bf) {
	if (zhou[0].CallCheckB3(tclues, nclgo, bf)) {
		anduab3 = BIT_SET_27;
		stopexpandb3 = 0;
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			if (!cc) {
				cout << " bug no b3 IsValidB3 [8]" << p_cpt2g[8]<< endl;
				aigstop = 1;	return 1;
			}
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			register uint32_t ua = w.bf.u32[2];
			tw3_2[ntw3_2++] = ua;
			taddgob3[ntaddgob3++] = ua;// more for a given b1+2
			anduab3 &= ua;
			BuildGua(w, cc);
			g_256h.Add128(w, cc,1);//clean level new if possible
			if (cc < 4) {
				if ((cc == 2 && scritb3.CanNotAddi27(w.bf.u32[2])) ||
					(cc == 3 && scritb3.CanNotAddi9(w.bf.u32[2]))) {
					stopexpandb3 = 1;
				}
			}
			if (cc0 > 20 || cc > 6)continue;
			chunkh.Add128(w, cc);
			if (cc < 3)p_cpt2g[35]++; else p_cpt2g[36]++;
		}
		return 1;
	}
	else return 0;
}

/*
uint32_t t2clues_band1[9] = {
0400000001,
0200000004,
0100000002,
 040000400,
 010000100,
 020000200,
  02000020,
  04000040,
  01000010,
};
*/

//=========brute force specific to this
int ZHOU::CallCheckB3(uint32_t * t, int n, uint32_t bf, int nogo) {// 17 search mode
	zh_g.diag = nogo;
	memcpy(this, zhoustart, sizeof zhoustart);
	misc.SetAll_0();
	BF128 dca[9];
	memset(dca, 0, sizeof dca);
	int digitsbf = 0;
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell], digit = zh_g2.grid0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		digitsbf |= 1 << digit;
		Assign(digit, cell, xcell);
		dca[digit].Set(xcell);
	}
	{
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = gchk.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	//__________end assign last lot start solver
	zhgxn.nua = 0;
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	int ir = Full17Update();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0);
	return zhgxn.nua;
}

int ZHOU::PartialInitSearch17(uint32_t * t, int n) {
	zh_g2.digitsbf = 0;
	memset(zh_g2.Digit_cell_Assigned, 0, sizeof zh_g2.Digit_cell_Assigned);
	memcpy(this, zhoustart, sizeof zhoustart);
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell], digit = zh_g2.grid0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		zh_g2.digitsbf |= 1 << digit;
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		zh_g2.Digit_cell_Assigned[digit].Set(xcell);
	}
	//cout << "ok init" << endl;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | zh_g2.Digit_cell_Assigned[i];
	return 0;
}

int ZHOU::Apply17SingleOrEmptyCellsB12() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint64_t R1, R2, R3, R4;
	{
		register uint64_t * P = FD[0][0].bf.u64, M = *P;
		R1 = M;
		P += 4; M = *P;	                            R2 = R1 & M;  R1 |= M;
		P += 4; M = *P;              R3 = R2 & M;   R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 = R3 & M;  R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P;	R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
	}
	if (cells_unsolved.bf.u64[0] & (~R1)) 	return 1; // empty cells
	R1 &= ~R2; // now true singles
	R1 &= cells_unsolved.bf.u64[0]; // these are new singles
	if (R1) {// singles to apply
		uint32_t xcell;
		while (bitscanforward64(xcell, R1)) {
			R1 ^= (uint64_t)1 << xcell;
			uint32_t cell = From_128_To_81[xcell];
			for (int idig = 0; idig < 9; idig++) {
				if (FD[idig][0].On(xcell)) {
					Assign(idig, cell, xcell);
					goto nextr1;
				}
			}
			return 1; // conflict with previous assign within this lot
		nextr1:;
		}
		zh_g.single_applied = 1;
		return 0;
	}
	// no single store apply  pair in priority ??
	R2 &= ~R3; // now true singles
	if (!R2) {
		R3 &= ~R4;
		if (R3) R2 = R3;
		else R2 = R4;
	}
	bitscanforward64(zh_g2.xcell_to_guess, R2);
	return 0;
}
int ZHOU::Apply17SingleOrEmptyCellsB3() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint32_t R1, R2, R3, R4;
	{
		register uint32_t * P = &FD[0][0].bf.u32[2], M = *P;
		R1 = M;
		P += 8; M = *P;	                            R2 = R1 & M;  R1 |= M;
		P += 8; M = *P;              R3 = R2 & M;   R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 = R3 & M;  R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P;	R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
	}
	if (cells_unsolved.bf.u32[2] & (~R1)) 	return 1; // empty cells
	R1 &= ~R2; // now true singles
	R1 &= cells_unsolved.bf.u32[2]; // these are new singles
	if (R1) {// singles to apply in band 3
		uint32_t xcell;
		while (bitscanforward(xcell, R1)) {
			R1 ^= (uint64_t)1 << xcell;
			uint32_t cell = xcell + 54;
			xcell += 64;
			for (int idig = 0; idig < 9; idig++) {
				if (FD[idig][0].On(xcell)) {
					Assign(idig, cell, xcell);
					goto nextr1;
				}
			}
			return 1; // conflict with previous assign within this lot
		nextr1:;
		}
		zh_g.single_applied = 1;
		return 0;
	}
	// no single store apply  pair in priority ??
	R2 &= ~R3; // now true pairs
	if (!R2) {
		R3 &= ~R4; // now true tripletss
		if (R3) R2 = R3;
		else R2 = R4;
	}
	bitscanforward(zh_g2.xcell_to_guess, R2);
	zh_g2.xcell_to_guess += 64;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		//if(diag &&cells_unsolved.bf.u32[2]==0)ImageCandidats();
		if (!Unsolved_Count()) return 2;
		if (cells_unsolved.bf.u32[2]) {// fill B3 first
			if (Apply17SingleOrEmptyCellsB3())	return 0; //  empty cell or conflict singles in cells
		}
		else {
			if ((!ISFALSEON))return 0;
			if (Apply17SingleOrEmptyCellsB12())	return 0;
		}
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Compute17Next(int index) {
	int ir = Full17Update();
	//if (zh_g.diag) {
		//cout << "index=" << index << endl;
		//ImageCandidats();
	//}
	if (!ir) return;// locked 
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128  wua;
			int * sol = gchk.grid0;
			wua.SetAll_0();;
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			if (wua.isNotEmpty()) {
				//if (zh_g.diag) {
					//cout << "zhgxn.nua=" << zhgxn.nua << endl;
					//wua.Print3(" ");
				//}
				int cc = _popcnt32(wua.bf.u32[2]);
				if ((!zhgxn.nua) || cc < 3)
					zhgxn.tua[zhgxn.nua++] = wua;
				if (cc < 3 || zhgxn.nua>5)	zh_g.go_back = 1;
			}
		}
		return;
	}
	Guess17(index);// continue the process
}
void ZHOU::Guess17(int index) {
	if (zh_g.go_back) return;
	int xcell = zh_g2.xcell_to_guess,
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell];
	// true first if possible
	if (FD[digit][0].On(xcell)) {
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->Compute17Next(index + 1);
		if (zh_g.go_back) return;
	}
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (ISFALSEON   && zhgxn.nua) continue;
		if (FD[idig][0].On(xcell)) {
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			if (cell >= 54)mynext->ISFALSEON++;
			mynext->SetaCom(idig, cell, xcell);
			mynext->Compute17Next(index + 1);
			if (zh_g.go_back) return;
		}
	}
}



void GCHK::Out17(uint32_t bfb3) {
	// mapping of the output 
	uint32_t map[81];
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 27; j++)
		map[27 * i + j] = bax[i].map81[j];
	char ws[82];

	uint32_t tcf[40], ntcf = 0;
	for (int i = 0; i < nclgo ; i++) {
		tcf[ntcf++] = map[tclues[i]];
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		tcf[ntcf++] = map[54 + i] ;
	strcpy(ws, empty_puzzle);
	for (uint32_t i= 0; i < ntcf; i++) {
		int cell = tcf[i];
		ws[cell] = ze[cell] ;
	}
	cout << ws << " one sol  size " <<  ntcf  << endl;
	if (modefirst) { nok = n18seen=aigstop = 1; return; }
	n18seen++;
	for (int i = 0; i < nok; i++) if (!strcmp(ws,t18found[i])) return;
	if (nok < 100) strcpy(t18found[nok++],ws);
	fout1 << ws << " nok=" << nok << endl;
}





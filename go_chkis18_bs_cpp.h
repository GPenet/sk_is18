#define DEBUG0 1
//#define DEBUG3 39
//#define DEBUG4 10
//#define DEBUG5 26113
//#define DEBUG5 1596
//#define DEBUG6 13 cas interessant 15 rien interessant
//#define DEBUG3 36
//#define DEBUG26 4905500
//#define DEBUG14 5
#define DEBUGB3N
#define DEBUGV4
#define DEBUGPH2

inline int  CHUNKS_HANDLER::Check2(BF128 w54) {
#ifdef DEBUGPH2
	if (gchk.aigstop)return 1;
	register uint32_t i27 = w54.bf.u32[2];
	register uint64_t u = w54.bf.u64[0];
	//cout <<Char54out(u)<< " try add i27=" << i27 << endl;
	for (uint32_t i = 0; i <= ic2; i++) {
		CHUNK3B & wc2 = c2[i];
		for (uint32_t j = 0; j < wc2.nt; j++) {
			CHUNK3B & wc2 = c2[i];
			if (wc2.tu27[j] != i27) continue;
			if (wc2.tu54[j] != u) continue;
			gchk.aigstop = 1;
			cout << " bugstop chexk redundant add c2 i="<<i<<" j="<<j
				<< "ic2 =" <<ic2 
				<< " [3] " << p_cpt2g[3] << " [4] " << p_cpt2g[4]
				<< " [5] " << p_cpt2g[5]<< "i 27="<<i27
				//<< " [4] " << p_cpt2g[4]
				<< endl;
			cout << Char9out(i27) << " ";
			cout << Char54out(u) << " exists " << endl;
			return 1;
		}
	}

#endif
	return 0;
}
inline int  CHUNKS_HANDLER::Check3(BF128 w54) {
#ifdef DEBUGPH2
	if (gchk.aigstop)return 1;
	register uint32_t i9 = w54.bf.u32[2];
	register uint64_t u = w54.bf.u64[0];
	//cout << Char54out(u) << " try add i9=" << i9 
	//	<< " [3] " << p_cpt2g[3] << " [4] " << p_cpt2g[4]
	//	<< " [5] " << p_cpt2g[5]		<< endl;
	for (uint32_t i = 0; i <= ic3; i++) {
		CHUNK3B & wc3 = c3[i];
		for (uint32_t j = 0; j < wc3.nt; j++) {
			CHUNK3B & wc3 = c3[i];
			if (wc3.tu27[j] != i9) continue;
			if (wc3.tu54[j] != u) continue;
			gchk.aigstop = 1;
			cout << " bugstop chexk redundant add c3 i=" << i << " j=" << j << "ic3 =" << ic3
				<< " [3] " << p_cpt2g[3] << " [4] " << p_cpt2g[4]
				<< " [5] " << p_cpt2g[5]
				//<< " [4] " << p_cpt2g[4]
				<< endl;
			cout << Char9out(i9) << " ";
			cout << Char54out(u) << " exists " << endl;
			return 1;
		}
	}

#endif
	return 0;
}
inline int  CHUNKS_HANDLER::Check4(BF128 w54) {
	if (gchk.aigstop)return 1;
	register uint32_t u3 = w54.bf.u32[2];
	register uint64_t u = w54.bf.u64[0];
	//cout << Char54out(u)<< " ";
	//cout << Char54out(u3) << " try add4" 
		//<< " [3] " << p_cpt2g[3] << " [4] " << p_cpt2g[4]
		//<< " [5] " << p_cpt2g[5]		<< endl;
	for (uint32_t i = 0; i <= ic4; i++) {
		CHUNK3B & wc3 = c4[i];
		for (uint32_t j = 0; j < wc3.nt; j++) {
			CHUNK3B & wc4 = c4[i];
			if (wc3.tu27[j] != u3) continue;
			if (wc3.tu54[j] != u) continue;
			gchk.aigstop = 1;
			cout << " bugstop chexk redundant add c4 "
				<< " [3] " << p_cpt2g[3] << " [4] " << p_cpt2g[4]
				<< " [5] " << p_cpt2g[5]
				//<< " [4] " << p_cpt2g[4]
				<< endl;
			return 1;
		}
	}
#ifdef DEBUGPH2

#endif
	return 0;
}

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
#ifdef DEBUG0 
	if (p_cpt2g[0] != DEBUG0) return 0;
#endif

#ifdef HAVEKNOWN
	// put knownn in the right order
	if (strlen(ze) < 163) return -1;// skip blank lines
	char * w = &ze[82],* ww=&zes[82];
	memcpy(w, &ww[27 * band_order[0]], 27);
	memcpy(&w[27], &ww[27 * band_order[1]], 27);
	memcpy(&w[54], &ww[27 * band_order[2]], 27);
	puzknown.SetAll_0();
	for (int i = 0; i < 81; i++) {
		char c = w[i];
		if (c<'1' || c>'9') continue;
		puzknown.Set_c(i);// store the pattern in 3X mode
		// this can be a pattern, no check of the digit with the solution
	}
	register uint64_t R = puzknown.bf.u64[0];
	pk54= (R & BIT_SET_27) |	((R & BIT_SET_B2) >> 5);
	if ((int)_popcnt32(puzknown.bf.u32[2] ) < mincluesb3) return  0;
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
#ifdef TEST_ON
	for (int i = 0; i < 81; i++) cout << zs0[i] + 1;
	cout << " reshaped  " << band_order[0] << band_order[1]
		<< band_order[2] << endl;
	for (int ibs = 0; ibs < 3; ibs++) {
		STD_B416 & b = bax[ibs];
		cout << b.band << "\t" << b.i416 << " " << t416n6[b.i416] << endl;
	}
	cout << " minb1=" << minb1 << " minb2=" << minb2 << " mincluesb3=" << mincluesb3 << endl;
#ifdef HAVEKNOWN
	for (int i = 0; i < 81; i++) cout << ze[i + 82];
	cout << " known reshaped" << endl;
#endif
#endif

	//__________________________ start uas search
	UaCollector();
	p_cpt2g[1] = tuasb12.nua;
	p_cpt2g[2] = chunkh.GetC2Count();
	//chunkh.DebugAll(1);
#ifdef TEST_ON
	tuasb12.Stats();
#endif
	BuildVectorsForExpand4B12();
	Expand4B12();
#ifdef TEST_ON
	tuasb12.Stats();
	chunkh.Stats();
	//chunkh.DebugAll(1);
#endif
	return a_18_seen;
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
	for (uint32_t i = 0; i < zh2gxn.nua; i++) {
		uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);
		AddUaB12UN(w, cc);
	}
}
BF128 wubuf_tfua[24][200];// put it out of the stack
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
			wubuf_tfua[cc][ntsort[cc]++] = wt;
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
			wubuf_tfua[cc][ntsort[cc]++] = wt;
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
			wubuf_tfua[cc][ntsort[cc]++] = wt;
		}
	}
	cout << "call FirstUasCollect reload" << endl;
	//for (int ii = 0; ii < 23; ii++) {
		//if (ii >= 20) break;
		//cout << ii << " " << ntsort[ii] << endl;
	//}

	// reload tall and check subsets/redundancy
	tuas81.ntall = 0;
	for (int i = 0; i < 23; i++)if (ntsort[i]) {
		BF128 * tt = wubuf_tfua[i];
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

	// check final status for c2
	cout << " last ic2 " << chunkh.ic2 << " last index "
		<< chunkh.c2[chunkh.ic2].nt
		<< " " << 64 * chunkh.ic2 + chunkh.c2[chunkh.ic2].nt << endl;
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

//________________________________ start expand 3/4 + no more uas
#define nphase1 5
GCHK::CPT_4::CPT_4() {
	t[0] = units3xBM[18].u64[0];	t[1] = units3xBM[19].u64[0];
	t[2] = units3xBM[20].u64[0];	t[3] = units3xBM[21].u64[0];
	t[4] = units3xBM[22].u64[0];	t[5] = units3xBM[23].u64[0];// 6 box
	t[6] = band3xBM[0].u64[0];		t[7] = band3xBM[1].u64[0];// 2 bands
	t[8] = band3xBM[3].u64[0];		t[9] = band3xBM[4].u64[0];
	t[10] = band3xBM[5].u64[0];		 //3 stack

}
void GCHK::BuildVectorsForExpand4B12() {
	//morev2a.Init();
	//morev2b.Init(); morev2c.Init();	morev2d.Init();
	// switch here tua to 54 mode and init add
	tuasb12.SwitchTo54Mode(); 

	for (int i = 0; i < 54; i++)memset(v12_4_c, 255, sizeof v12_4_c);
	for (uint32_t i = 0; i < 64; i++) {
		register uint64_t R = tuasb12.tua[i];
		uint32_t cell;
		while (bitscanforward64(cell, R)) {
			uint64_t bit = (uint64_t)1 << i;
			R ^= (uint64_t)1 << cell; //clear bit
			v12_4_c[cell] ^= bit;
		}
	}
}
void  GCHK::Expand4B12() {
	if (aigstop) return;

	zh2b[0].InitBands12(grid0);
	uint64_t *twu = tuasb12.tua,
		limspot = nphase1-1;// expand 5 clues
	nt4_to_expand = 0;
	// _______________ expand to have all  3/4 cells to expand more
	struct SPB {// spots to find starts to extract uas 3 digits
		uint64_t  possible_cells, all_previous_cells, active_cells,
			v;
	}spb[6], *s, *sn;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->active_cells = maskLSB[54].u64[0];
	s->possible_cells = twu[0] ;
	s->v = ~0;// initial nothing done
	//____________ here start the search nclues
next:
	uint64_t ispot = s - spb;
	// catch and apply cell in bitfields
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)goto back;
	{
		bitscanforward64(cell, p);
		register uint64_t bit = (uint64_t)1 << cell;
		s->possible_cells ^= bit;
		tclues[ispot] = cell;
		s->active_cells ^= bit;
		uint64_t ac = s->active_cells;
		sn = s + 1; *sn = *s; // (copy the stack count)
		sn->all_previous_cells = s->all_previous_cells | bit;
		if (ispot == limspot) {// 4 cells to use later
			T4_TO_EXPAND & w = t4_to_expand[nt4_to_expand++];
			BF64 vsort; vsort.bf.u64 = 0;
			w.bf = sn->all_previous_cells;
			w.active = ac;
			vsort.bf.u32[1] = cpt_4c.GetCount(w.bf);
			vsort.bf.u32[0] = 27 - (uint32_t)_popcnt64(ac);
			w.vsort = vsort.bf.u64;
			goto next;
		}
		// find next ua
		sn->v &= v12_4_c[cell];
		if (!sn->v) {// no more uas should never be
			aigstop = 1;
			cout << "no more uas in expand 5" << endl;		
			return;
		}
		uint32_t ir;
		bitscanforward64(ir, sn->v);//ir ua to load
		uint64_t Ru = twu[ir] & ac;
		if (!Ru)goto next;//dead branch unlikely
		sn->possible_cells = Ru;
		s++; // switch to next spot
		goto next;
	}
	// going back, for a non empty index, count it back
back:
	if (--s >= spb)goto next;
	// sort output
	for (uint64_t i = 0; i < nt4_to_expand - 1; i++) {
		for (uint64_t j = i + 1; j < nt4_to_expand; j++) {
			if (t4_to_expand[i].vsort < t4_to_expand[j].vsort) {
				GCHK::T4_TO_EXPAND temp = t4_to_expand[i];
				t4_to_expand[i] = t4_to_expand[j];
				t4_to_expand[j] = temp;
			}
		}
	}
#ifdef HAVEKNOWN
	nt4ok = -1;
	{
		register uint64_t K = pk54, Kn = ~K;
		for (uint64_t i = 0; i < nt4_to_expand; i++) {
			T4_TO_EXPAND w=t4_to_expand[i];
			register uint64_t B = w.bf, A = w.active;
			if (B & Kn)continue;// not partial known
			if (K & (~(B | A))) continue;
			if (nt4ok < 0) nt4ok = (int)i;
			cout << Char54out(B) << " bf" << endl;
			cout <<Char54out(A) << " ac" <<endl;
		}
	}
	cout << "expected i=" << nt4ok << endl;
#endif
#ifdef TEST_ON
	cout << "nt4_to_expand=" << nt4_to_expand << endl;
#endif
	for (uint64_t i = 0; i < nt4_to_expand; i++) {
#ifdef HAVEKNOWN
		if(nt4ok>=0)		okcheck = (i == nt4ok);
		else okcheck = 0;
		if (okcheck) cout <<i  << " start expand the expected ckunk" << endl;
#endif
#ifdef DEBUG3
		if (p_cpt2g[3] > DEBUG3) return;
#endif
		Do_phase2(t4_to_expand[i]);
		if (aigstop) 			return;
	}
}

#define PH2_128
//______  process a "no more uas bands 1+2
void GCHK::Do_phase2(T4_TO_EXPAND w) {
	p_cpt2g[3]++;
	register uint64_t bf54 = w.bf, ac54 = w.active;
	myac_4 = ac54;


	{	//____ load the first clues
		register uint64_t R = bf54;
		nclues_step = 0;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			tclues[nclues_step++] = cell;
		}
		//tcluesxy = &tclues[nclues_step];
	}
	uint32_t   n = 0;
	{// build still valid uas sorted by size 

		//uint64_t wubuf.tdophase2[30][4000];
		uint64_t  ntt[30];
		memset(ntt, 0, sizeof ntt);
		{
			for (uint32_t i = 0; i < tuasb12.nua; i++) {
				register uint64_t R = tuasb12.tua[i];
				if (R&bf54) continue;
				R &= ac54 ;
				if (!R)return; //dead
				register uint64_t cc = _popcnt64(R);
				n++;
				if (cc > 20)cc = 20;
				wubuf.tdophase2[cc][ntt[cc]++] = R;
			}
			for (uint32_t i = 0; i < tuasb12.nta17; i++) {
				register uint64_t R = tuasb12.ta17[i];
				if (R&bf54) continue;
				R &= ac54;
				if (!R)return; //dead
				register uint64_t cc = _popcnt64(R);
				n++;
				wubuf.tdophase2[cc][ntt[cc]++] = R;
			}
			if (n < 5000)for (uint32_t i = 0; i < tuasb12.nta18; i++) {
				register uint64_t R = tuasb12.ta18[i];
				if (R&bf54) continue;
				R &= ac54;
				if (!R)return; //dead
				register uint64_t cc = _popcnt64(R);
				n++;
				wubuf.tdophase2[cc][ntt[cc]++] = R;
			}
			if (n < 5000)for (uint32_t i = 0; i < tuasb12.nta19; i++) {
				register uint64_t R = tuasb12.ta19[i];
				if (R&bf54) continue;
				R &= ac54;
				if (!R)return; //dead
				n++;
				register uint64_t cc = _popcnt64(R);
				n++;
				wubuf.tdophase2[cc][ntt[cc]++] = R;
			}
			if (n < 5000)for (uint32_t i = 0; i < tuasb12.nta20; i++) {
				register uint64_t R = tuasb12.ta20[i];
				if (R&bf54) continue;
				R &= ac54;
				if (!R)return; //dead
				register uint64_t cc = _popcnt64(R);
				n++;
				wubuf.tdophase2[cc][ntt[cc]++] = R;
			}
			if (n < 5000)for (uint32_t i = 0; i < tuasb12.nta21; i++) {
				register uint64_t R = tuasb12.ta21[i];
				if (R&bf54) continue;
				R &= ac54;
				if (!R)return; //dead
				register uint64_t cc = _popcnt64(R);
				n++;
				wubuf.tdophase2[cc][ntt[cc]++] = R;
			}
			if (n < 5000)for (uint32_t i = 0; i < tuasb12.nta22; i++) {
				register uint64_t R = tuasb12.ta22[i];
				if (R&bf54) continue;
				R &= ac54;
				if (!R)return; //dead
				register uint64_t cc = _popcnt64(R);
				n++;
				wubuf.tdophase2[cc][ntt[cc]++] = R;
			}

		}

		{
			ntua4 = 0;
			uint32_t lim = 0;
			for (int i = 0; i < 30; i++)if (ntt[i]) {
				register uint64_t * tj = wubuf.tdophase2[i];
				if (i < 14) lim = (ntua4 < 64) ? ntua4 : 64;
				for (uint64_t j = 0; j < ntt[i]; j++) {
					register uint64_t R = tj[j], Rn = ~R;
					//___ clear redundancy
					for (uint32_t k = 0; k < j; k++)
						if (R == tj[k])goto nextj;
					//___clear if subset
					for (uint32_t k = 0; k < lim; k++)
						if (!(tua4[k] & Rn)) goto nextj;
					if (ntua4 < 4000)tua4[ntua4++] = R;
				nextj:;
				}
			}
		}
	}

	//______________________________________________end pick up
	//cout << "ntua4=" << ntua4 << endl;
	{
		register int wnua = ntua4;
		memset(vph2, 0, sizeof vph2);
		for (int i = 0; i < 32; i++) {//max 32*128=4096
			if (wnua >= 128) { vph2[i].SetAll_1();  wnua -= 128; }
			else {//0 127 last
				vph2[i] = maskLSB[wnua];
				break;
			}
		}
	}
	nvph2 = (ntua4 + 127) >> 7; // number of vectors size 64
	memset(vph2c, 255, sizeof vph2c);// each vector 8 bytes
	for (uint32_t i = 0; i < ntua4; i++) {
		register uint64_t R54 = tua4[i];
		//cout << Char54out(R54) << " " << i << endl;
		register uint32_t cell, bloc = i >> 7, ir = i - 128 * bloc;
		while (bitscanforward64(cell, R54)) {
			R54 ^= (uint64_t)1 << cell; //clear bit
			vph2c[bloc][cell].clearBit(ir);
		}
	}
#ifdef DEBUG3
	if (p_cpt2g[3] == DEBUG3) {
		cout << "start [3]=" << p_cpt2g[3] << " ntua4=" << ntua4 << endl;

	}
#endif
	if (ntua4 < 128) ntua4 = 128; // push adds out of firstbloc
	Do_phase2Expand_128(bf54, ac54);
	
}
void  GCHK::Do_phase2Expand_128(uint64_t bf, uint64_t ace) {// mode 54 not 2x
#ifdef HAVEKNOWN
	if (okcheck) {
		okcheck = 1;
		if ((bf&pk54) == bf) {
			cout << Char54out(bf) << " possible expand p_cpt2g[3]=" << p_cpt2g[3] << endl;
			cout << Char54out(ace) << " active here " << _popcnt64(ace) << endl;
			uint64_t w = (ace &pk54) | bf;
			cout << Char54out(w) << " ace & sol";
			if (w == pk54) { cout << " will be ok"; okcheck = 2; }
			cout << endl;
		}
	}
#endif
	pendbufvalid = &bufvalid[BUFVALIDS];
	pbufvalid = bufvalid;
	uint64_t limspot = 17 - nphase1 - mincluesb3 - 1, // 11 or 12 clues
		limclean = 5 - nphase1;// set to more if low {min b1+b2}
	if (minb1b2 < 9)limclean = 6 - nphase1;
	uint64_t andvalid = ~0, orvalid = 0;
	// _______________ expand to have all  minimal < n clues
	struct SPB {// spots to find starts to extract uas 3 digits
		BF128 v128;
		uint64_t filler, possible_cells, all_previous_cells, active_cells;
	}spb[13];
	register SPB  *s, *sn;
	register uint64_t ispot;
	register int cell;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->all_previous_cells = bf;
	s->active_cells = ace;
	s->possible_cells = tua4[0];
	s->v128 = vph2[0];
	//____________ here start the search nclues
next:
	ispot = s - spb;
	{// catch and apply cell in bitfields
		{// next cell
			register uint64_t p = s->possible_cells;
			if (!p)goto back;
			bitscanforward64(cell, p);
			register uint64_t bit = (uint64_t)1 << cell;
			s->possible_cells ^= bit;
			s->active_cells ^= bit;
			sn = s + 1; *sn = *s;
			sn->all_previous_cells |= bit;
			sn->v128 &=  vph2c[0][cell]; 
		}

		register uint64_t Ru,V;
		if ((V=sn->v128.bf.u64[0])) {// next ua
			register uint32_t ir;
			bitscanforward64(ir, V);//relative index first active
			Ru = tua4[ir] & s->active_cells;
		}
		else if ((V = sn->v128.bf.u64[1])) {// next ua
			register uint32_t ir;
			bitscanforward64(ir, V);//relative index first active
			Ru = tua4[ir+64] & s->active_cells;
		}
		else {
			// possible valid band 1+2 to store
			{//__________ check count per band
				register uint64_t n12 = nphase1 + 1 + ispot,
					n1 = _popcnt64(sn->all_previous_cells & BIT_SET_27);
				if (n1 > 7 || (n12 - n1) > 8)goto next;//can be 288
			}
			if (pbufvalid >= pendbufvalid) {// if limite do
				CleanBufferAndGo(andvalid, orvalid);
				andvalid = ~0; orvalid = 0;
			}
			andvalid &= sn->all_previous_cells;
			orvalid |= sn->all_previous_cells;
			*pbufvalid++ = sn->all_previous_cells;
			if (ispot != limspot + 1)*pbufvalid++ = s->active_cells;
			goto next;
		}

		//_________________ now a new ua to process
		if (!Ru)goto next;
		//____if last special process
		if (ispot == limspot) {// 4 for 11 ; 5 for 12
			// clear potential redundancy assumed nphase1=5
			{ 
				register uint64_t n12 = nphase1 + 1 + ispot,
					n1 = _popcnt64(sn->all_previous_cells & BIT_SET_27),
					n2 = 10 - n1; // 10 or 11 already there
				if (limspot == 5) {//12 clues only 666
					n2++;
					if (n1 != 6)Ru &= BIT_SET_27;// can only be band1
					if (n2 != 6)Ru &= ~BIT_SET_27;// can only be band2
				}
				else {// 11 clues can  be  747 567 657
					if (n1 > 7 || n2 > 6)goto next;
					if (n1 == 7)Ru &= ~BIT_SET_27;// only  band2
					if (n2 == 6)Ru &= BIT_SET_27;// only  band1
				}
			}
			if (!Ru)goto next;
			// now all cells hitting all uas can be added
			{
				register uint64_t Bf = sn->all_previous_cells;
				while (bitscanforward64(cell, Ru)) {
					register uint64_t bit = (uint64_t)1 << cell
						, Bff = Bf | bit;
					Ru ^= bit;
					BF128 vv= sn->v128 & vph2c[0][cell];
					if (vv.bf.u64[0] || vv.bf.u64[1]) continue;
					andvalid &= Bff;
					orvalid |= Bff;
					*pbufvalid++ = Bff;
				}
			}
			if (pbufvalid >= pendbufvalid) {// if limite do
				CleanBufferAndGo(andvalid, orvalid);
				andvalid = ~0; orvalid = 0;
			}
			goto next;
		}
		else {
			if (ispot == limclean) {// clean 2+3/4 clues or more
				if (pbufvalid > bufvalid) {
					CleanBufferAndGo(andvalid, orvalid);
					andvalid = ~0; orvalid = 0;
				}
			}

			sn->possible_cells = Ru;
			s++; // switch to next spot
			goto next;

		} 
	}
	goto next;
back:
	if (--s >= spb)goto next;
	if (pbufvalid > bufvalid) // if limite do
		CleanBufferAndGo(andvalid, orvalid);
}

// ________using guas2 guas3 after expand


void GCHK::BuildReducedGuasVectors23() {
	chunkh.ApplyClean(tclues, nclues_step);
	// _______________extract guas2 and vectors
	BF128 t[500], tw;// temp storage for still active
	//use wubuf.tg2[25][100];	
	uint32_t  ntt[25],nt=0,n2=0,n3=0,nm=0; 
	memset(ntt, 0, sizeof ntt);
	// split guas2   in 2 x 256 vectors ( room for more )
	g2_256[0].Init();  g2_256[1].Init(); g3_256.Init();
	gm_256[0].Init();  gm_256[1].Init();

	for (uint32_t ic2 = 0; ic2 <= chunkh.ic2; ic2++) {
		CHUNK3B &w = chunkh.c2[ic2];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			register uint32_t i27 = w.tu27[ir]; 
			register uint64_t pat12= w.tu54[ir]& myac_4;
			tw.bf.u64[0] = pat12;
			tw.bf.u32[3] = i27;
			tw.bf.u32[2] =1<< i27;
			register uint32_t cc= (uint32_t)_popcnt64(pat12);
			// below cc=10 load direct else sort by size
			if (cc < 10) { g2_256[0].Add(tw); n2++; continue;}
			if (cc < 20) cc -= 10; else cc = 9;
			wubuf.tg2[cc][ntt[cc]++]=nt;
			t[nt++] = tw;
		}
	}
	for (uint32_t i = 0; i < 10; i++) {
		uint32_t * tti = wubuf.tg2[i];
		for (uint32_t j = 0; j < ntt[i]; j++) {
			BF128 w = t[tti[j]];
			if(n2++<256)g2_256[0].Add(w);
			else g2_256[1].Add(w);
		}
	}
	// ________________________________________extract guas3
	for (uint32_t ic3 = 0; ic3 <= chunkh.ic3; ic3++) {
		CHUNK3B &w = chunkh.c3[ic3];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			uint32_t i9 = w.tu27[ir];
			guash.AddG3(i9, w.tu54[ir]);
			register uint64_t pat12 = w.tu54[ir] & myac_4;
			tw.bf.u64[0] = pat12;
			tw.bf.u32[3] = i9;
			tw.bf.u32[2] = 1 << i9;
			g3_256.Add(tw);
			n3++;
		}
	}
	tw.bf.u32[3] = 0;
	// _____________________________________________extract guas4
	for (uint32_t ic4 = 0; ic4 <= chunkh.ic4; ic4++) {
		CHUNK3B &w = chunkh.c4[ic4];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & myac_4;
			tw.bf.u32[2] = w.tu27[ir];
			gm_256[0].Add(tw);
			nm++;
		}
	}
	// _____________________________________________extract guas5
	for (uint32_t ic5 = 0; ic5 <= chunkh.ic5; ic5++) {
		CHUNK3B &w = chunkh.c5[ic5];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & myac_4;
			tw.bf.u32[2] = w.tu27[ir];
			gm_256[0].Add(tw);
			nm++;
		}
	}
	// _____________________________________________extract guasmore
	for (uint32_t icm = 0; icm <= chunkh.icmore; icm++) {
		CHUNK3B &w = chunkh.cmore[icm];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & myac_4;
			tw.bf.u32[2] = w.tu27[ir];
			if (nm++ < 256)gm_256[0].Add(tw);
			else gm_256[1].Add(tw);
		}
	}
	g_256h.GetStart(n2, n3, nm);
	if(n2>340)
	cout << Char54out(myac_4)<< "BuildReducedGuasVectors23 "
		<<"n2=" << n2 << " n3=" << n3	<< " nm=" << nm << endl;

}

#define stack1_54 07007007007007007
int tstack27[3] = { 07007007 ,070070070,0700700700 };
struct CRITB3 {
	uint32_t minix[4],// triplet bf1 bf2 bf3  
		critbf, pairs27,  mincount,
		critstack,stackmin[3]	, stackf[3],
		assigned,active,
		ncl,nb3,nmiss;
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

		if (minix[3] & bitmini) {// can noo be
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
			critbf |= (mask^bit27);
		}
		stackf[stack]++;
		nmiss--;
		if (!nmiss) {
			critstack = 7;
			active &= critbf;// no more outfield
		}
		else if (stackf[stack] == nb3) {
			critstack |= bitstack;
			active &= (~tstack27[stack]) || critbf;
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
			active &= (~tstack27[stack]) || critbf;
		}
		return 0;
	}

	inline int Addone(uint32_t i27,int diag=0) {// back 1 if not possible
		int bit27 = 1 << i27, stack = C_stack[i27];
		int imini = i27 / 3, bitmini = 1 << imini, bitstack = 1 << stack;
		assigned |= bit27;
		if (!(bit27 & critbf)) {// clue added outfield
			if (diag) cout << "added outfield" << endl;
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
			if (diag) cout << "2 clues expected" << endl;
			critbf ^= bit27;
			minix[3] ^= bitmini;
			return 0;
		}
		if ((minix[2] & bitmini) && (pairs27 & bit27)) {
			if (diag) {
				cout << Char27out(bit27) << " bit27 2 pairs not common" << endl;
				cout << Char27out(pairs27) << " pairs 27" << endl;
			}
			// 2 clues if not common clue one more clue
			if (critstack & bitstack)return 1;// not possible
			critbf ^= bit27;
			minix[2] ^= bitmini;
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
		if (diag) cout << "one clue in field" << endl;

		register int mask = ~(7 << (3 * imini));// clear minirow
		critbf &= mask;
		//active &= mask;// likely not		
		return 0;
	}

	inline int AddAssign(uint32_t bfa) {// back 1 if not possible
		if (assigned & bfa)return 1; //should never be
		active &= ~bfa; // minimum is to kill new assign
		register int i27, X = bfa;
		//cleanmini= scritb3.minix[0] | scritb3.minix[1];// one pair or triplet
		while (bitscanforward(i27, X)) {
			if (Addone(i27))return 1;
		}
		return 0;
	}

	inline void Init(int ncb12) {
		memset(this, 0, sizeof(*this));
		ncl = ncb12, nb3 = 18 - ncl;
	}
	void AssignBf2(int bf) {// all or part of the bf2
		minix[2] &= ~bf;// clean the bf2 
		int imini;
		while (bitscanforward(imini, bf)) {
			int  mask = 7 << (3 * imini),m27=pairs27&mask;
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
				if(! (critstack & bit)) active |= mask;
		}
		else active = BIT_SET_27;
		if ((!critstack)|| (!minix[2])) return;
		for (int i = 0, bit = 1, mask = 0111; i < 3; i++, bit <<= 1, mask <<= 1) {
			if (critstack & bit) {
				active &= ~(stack1_54 << 3 * i); //  (stack154 << 3*i)active must be in critbf for this; stack
				if(minix[2] & mask)		AssignBf2(minix[2] & mask);

			}
		}
		active |= critbf;// be sure to keep it

	}

	inline void MinVmini(int imini,int vmini) {
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
			else critbf |= (mask^vminishift);
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
		stackf[0] += (uint32_t)_popcnt64(bf&stack1_54);
		stackf[1] += (uint32_t)_popcnt64(bf&(stack1_54 << 3));
		stackf[2] += (uint32_t)_popcnt64(bf&(stack1_54 << 6));
		nmiss = nb3 - mincount;
		if (!nmiss) critstack = 7;
	}
	void Status(const char * lib) {
		cout << lib << "critical Status mincount =" << mincount <<" nmiss="<<nmiss
			<<" critstack="<< critstack << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char27out(assigned) << " assigned" << endl;
		cout << Char27out(active) << " active" << endl;
		if (minix[1])cout << Char9out(minix[1]) << "     minis bf1" << endl;
		if (minix[2])cout << Char9out(minix[2]) << "     minis bf2" << endl;
		if (minix[3])cout << Char9out(minix[3]) << "     minis bf3" << endl;
		if (minix[0])cout << Char9out(minix[0]) << " mini triplets" << endl << endl;
		cout << "minis stacks " << stackmin[0] << stackmin[1] << stackmin[2] << endl;
		cout << "final stacks " << stackf[0] << stackf[1] << stackf [2] << endl;
	}

}scritb3;


inline void TUASB12::Add(uint64_t u, uint64_t cc) {
	if (cc < 18) {
#ifdef DEBUGPH2
		if (cc == 16) {
			for (uint32_t i = 0; i < nta17; i++) {
				if (ta17[i]  != u) continue;
				gchk.aigstop = 1;
				cout << Char54out(u) << " exists redundant 16 stop "
					<< " [3] " << p_cpt2g[3] << " [4] " << p_cpt2g[4]
					<< " [5] " << p_cpt2g[5] << endl;
				return;
			}
		}
#endif

		p_cpt2g[31]++;
		if (cc < 14) {	if (nua < 4096)	tua[nua++] = u;	}
		else if (nta17 < 2000)ta17[nta17++] = u;
		return;
	}
	switch (cc) {
	case 18:if (nta18 < 2000)ta18[nta18++] = u;
		p_cpt2g[32]++; return;
	case 19:if (nta19 < 4096) ta19[nta19++] = u;
		p_cpt2g[32]++; return;
	case 20:if (nta20 < 4096) ta20[nta20++] = u;
		p_cpt2g[33]++; return;
	case 21:if (nta21 < 4096) ta21[nta21++] = u;
		p_cpt2g[33]++; return;
	case 22:if (nta22 < 4096) ta22[nta22++] = u;
		p_cpt2g[33]++; return;
	}
}
inline void GCHK::LoadRemainningClues(uint64_t bf) {
	ncluesb12 = 0;
	register uint64_t R = bf ^ myandvalid;
	uint32_t cell;
	while (bitscanforward64(cell, R)) {
		R ^= (uint64_t)1 << cell; //clear bit
		tcluesxy[ncluesb12++] = cell;
	}
}
void GCHK::AddTua4(uint64_t v) {
	if (ntua4 >= 4096)return;
	register uint32_t cell, bloc = ntua4 >> 7,
		ir = ntua4 - 128 * bloc;
	tua4[ntua4++] = v;
#ifdef DEBUG3
	if (p_cpt2g[3] == DEBUG3) {
		cout << Char54out(v) << "add [3]=" << p_cpt2g[3] << " ntua4=" << ntua4
			<< " bloc=" << bloc << " ir=" << ir << endl;

	}
#endif
	vph2[bloc].setBit(ir);
	if (!ir) nvph2++; // a new 64 bits bloc
	while (bitscanforward64(cell, v)) {
		v ^= (uint64_t)1 << cell;; //clear bit
		vph2c[bloc][cell].clearBit(ir);
	}

}
int GCHK::IsValidB12(uint32_t ncl) {
	p_cpt2g[29]++;
	//mbisvalid.ntm = 0;

	if (zh2b[1].IsValid(tclues, ncl)) {
		p_cpt2g[30]++;
		int locdiag = 0;
#ifdef DEBUG5
		if (p_cpt2g[5] == DEBUG5) {
			cout << "not valid b12 zh2gxn.nua= " << zh2gxn.nua << endl;
			locdiag = 1;
		}
#endif
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i],
				cc = _popcnt64(ua),
				ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			if (locdiag) {
				cout << Char54out(ua54)<< " i="<<i << endl;
			}
			if (cc == 12 || cc==16)		cout << Char54out(ua54) << " add"<<cc<<" [3]="
				<< p_cpt2g[3] << " [4]=" << p_cpt2g[4] << " [5]=" << p_cpt2g[5] << endl;
			tuasb12.Add(ua54, cc);  // global level
			m64vh.Add54(ua54); // 5 clues  level
			m64vh.n++;
			AddTua4(ua54);// chunk level
			mbisvalid.tm[mbisvalid.ntm++] = ua54;
			//if (p_cpt2g[30]<10)		cout << Char54out(ua54) << " added validb12 [4]="
			//	<< p_cpt2g[4] <<" [5] " << p_cpt2g[5] << endl;
			//if (moreand.ntm < 500)moreand.tm[moreand.ntm++] = ua54;
		}
	}
	nclf = ncl;
	return (int)mbisvalid.ntm;
}

//_____ process a chunk of potential valid bands 1+2

void GCHK::CleanBufferAndGo(uint64_t andvalid, uint64_t orvalid) {
	if (aigstop) return;
	p_cpt2g[4]++;
	uint64_t *myend = pbufvalid, *pw = bufvalid; pbufvalid = bufvalid;	
	myandvalid = andvalid;	myorvalid = orvalid;
	uint64_t nn8 = (myend - bufvalid),	mincb2 = _popcnt64(andvalid >> 27);
#ifdef HAVEKNOWN
	if ((andvalid&pk54) == andvalid) 
		cout << Char54out(andvalid) << " p_cpt2g[4]=" << p_cpt2g[4] << endl;
#endif
#ifdef DEBUG4
	if (p_cpt2g[4] > DEBUG4) {	aigstop = 1; return;}
	if (p_cpt2g[4] <= DEBUG4)
		cout << Char54out(andvalid) << " andv  p_cpt2g[4] " << p_cpt2g[4] << " size" << nn8 << endl;
#endif
	{//____ load the common clues
		register uint64_t R = myandvalid;
		nclues_step = 0;
		uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			tclues[nclues_step++] = cell;
		}
		tcluesxy = &tclues[nclues_step];
	}
	{//___________ check other uas >64 and "more"
		m64vh.Init();
		register uint32_t * tcl = tclues, ntcl = nclues_step, x;
		// re visit first 128 for adds (forcedout )
		for (uint32_t i = 1; i < nvph2; i++) {
			BF128 vw = vph2[i];
			for (uint32_t j = 0; j < ntcl; j++) vw &= vph2c[i][tcl[j]];
			register uint64_t V = vw.bf.u64[0];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x; //clear bit
				register uint64_t U = tua4[x + 128 * i];
				m64vh.Add54(U);
			}
			V = vw.bf.u64[1];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x; //clear bit
				register uint64_t U = tua4[x + 128 * i+64];
				m64vh.Add54(U);
			}		
		}

		m64vh.n = m64vh.nmv * 128 + m64vh.mv[m64vh.nmv].nt;
		if (m64vh.n > p_cpt2g[49])p_cpt2g[49] = m64vh.n;
	}
	moreand.ntm = 0;
	BuildReducedGuasVectors23();
#ifdef DEBUG4
	if (p_cpt2g[4] == DEBUG4) {
		//guash3.Dumpg2g3();		
		//guash3.Dumpg2g3Det();
	}
	//else {		aigstop = 1; return;	}
#endif
	gaddb3.Init();
	ntuaclean = 0;// checking uas bands 1+2 added in this cycle

	//___________________________________ potential valid bands +2
	while (pw < myend) {
		if (aigstop) return;
		p_cpt2g[5]++;
		myb12f = myb12 = *pw++;
		int nclues = (int)_popcnt64(myb12);
		if (nclues < limb12)myac = *pw++;	else myac = 0;
#ifdef HAVEKNOWN
		if ((myb12&pk54) == myb12) {
			cout << Char54out(myb12) << " p_cpt2g[5]=" << p_cpt2g[5]
				<< " ncl=" << ncl << " limb12=" << limb12 << endl;
			cout << Char54out(myac) << " still active " << _popcnt64(myac) << endl;
		}
#endif

#ifdef DEBUG5
		int locdiag = 0;
		if (p_cpt2g[5] > DEBUG5) return;
		if (p_cpt2g[5] == DEBUG5) {
			cout << endl << endl << Char54out(myb12) << " myb12  (4 ) " << p_cpt2g[4]
				<< " [5] " << p_cpt2g[5] << endl;
			cout << "tclues xy";
			for (int i = 0; i < ncluesb12; i++)cout << tcluesxy[i] << " ";
			cout << " m64vh.nmv=" << m64vh.nmv << endl;
			//m64vh.Dump();
			locdiag = 1;
		}
#endif
		LoadRemainningClues(myb12);
		gaddb3.Init();
		ntaddgob3 = 0;
		mbisvalid.ntm = 0;
		clean_valid_done = 0;
		//_______________ draft of the new process
		if (nclues == limb12 ) {//(! (p_cpt2g[5] &7)) {
			if (m64vh.n && m64vh.ApplyXY(tcluesxy, ncluesb12)) 	continue;
			p_cpt2g[6]++;
			InitGoB3(myb12,nclues);
			continue;
		}
		else {
			mab.ntm = 0;
			if (m64vh.n)	m64vh.Extract(tcluesxy, ncluesb12, mab.tm, mab.ntm);
#ifdef DEBUG5
			if (locdiag) {
				cout << " mab.ntm= " << mab.ntm << endl;
				for (uint64_t i = 0; i < mab.ntm; i++) {
					register uint64_t R = mab.tm[i];
					cout << Char54out(R) << "  " << _popcnt64(R) << endl;
				}
			}
#endif
			return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			if (!mab.ntm) {// check b12 now
				if (!IsValidB12(nclues)) {// add one cell hitting all uas
					clean_valid_done = 1;
					InitGoB3(myb12,nclues);// try direct
					CheckValidBelow(myb12, myac);
					continue;
				}
			}
			p_cpt2g[7]++;	// not limit can add clues in bands 1+2
			// now at least one ua b12 still valid
			int nfree = limb12 - nclues;
#ifdef DEBUG5
			if (p_cpt2g[5] == DEBUG5) {
				cout << " mab.ntm=" << mab.ntm << " nfree=" << nfree << endl;
			}
#endif
			if (nfree == 1) {//11 clues 747 567 657 12 clues 666
				uint64_t a = mab.GetAnd(myac);
				if (!a) { return; }// not one clue
				//__ apply limits depending on the number of clues in bands 1+2
				uint64_t n1 = _popcnt64(myb12 & BIT_SET_27),
					n2 = nclues - n1; // 10 or 11 already there
				if (n2 < 3 || n2>6)return; // none of the 4 possibles
				if (n2 == 6)a &= BIT_SET_27;//can only add in band 1
				else if (n2 == 3)a &= ~BIT_SET_27;//must add in band 2
				else if (limb12 == 12 && n2 == 5)a &= ~BIT_SET_27;//must add in band 2
				if (!a) { return; }// no expected valid
				int x, ncluer = nclues++;
				while (bitscanforward64(x, a)) {
					tclues[ncluer] = x;
					uint64_t bit = (uint64_t)1 << x;
					a ^= bit; //clear bit
					myb12f = myb12 | bit;
					InitGoB3(myb12f,nclues);
				}
				continue;
			}
			p_cpt2g[17]++;
			p_cpt2g[50 + nfree]++;
			CleanMoreUas(myb12, myac, nclues, mab);
		}  
	}
}


int GCHK::BuildB2Table() {// extract and reorder still possible guas
	uint32_t is1 = 0, ntt[6], i27;
	memset(ntt, 0, sizeof ntt);
	register uint32_t AC = scritb3.active,
		F = scritb3.assigned;
	{
		register uint32_t V = scritb3.pairs27;
		while (bitscanforward64(i27, V)) {
			V ^= 1 << i27;
			wubuf.tbuilexpand3[2][ntt[2]++] = tg2[i27].pat;
		}
		V = scritb3.minix[0];// triplets
		if (V) {// add guas3 if any (usually 0)
			for (int i = 0; i < 9; i++)if (V & (1 << i))
				wubuf.tbuilexpand3[3][ntt[3]++] = 7 << (3 * i);
		}

	}

	register uint32_t * t2 = wubuf.tbuilexpand3[2];
	for (uint32_t i = 0; i < ntw3; i++) {
		register uint32_t U = tw3[i];
		if (!(U&F)) {
			U &= AC;
			int ccu = _popcnt32(U);
			if (!ccu) return 1;
			if (ccu == 1) { is1 |= U; continue; }
			{// clear redundancy
				register uint32_t nu = ~U, x = 0;
				for (uint32_t i = 0; i < ntt[2]; i++)
					if (!(nu&t2[i])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			wubuf.tbuilexpand3[ccu][ntt[ccu]++] = U;
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
				for (uint32_t i = 0; i < ntt[2]; i++)
					if (!(nu&t2[i])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			wubuf.tbuilexpand3[ccu][ntt[ccu]++] = U;
			//cout << Char27out(U) << "b3 added" << endl;
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
				for (uint32_t i = 0; i < ntt[2]; i++)
					if (!(nu&t2[i])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			wubuf.tbuilexpand3[ccu][ntt[ccu]++] = U;
		}
	}

	ntw3_2 = is1_tw3=0;
	if (is1) {// now singles to assign
		if (scritb3.AddAssign(is1)) {
			//cout << Char54out(is1) << "is1 in excess" << endl;
			return 2;
		}
		AC = scritb3.active;	F = scritb3.assigned;
		int cc = _popcnt32(F);		if (cc > ncluesb3)return 1;
		// after added clues, revise the others
		uint32_t  ntt2[6];//tt2[i][500]is wubuf.tbuilexpand3[i+5][500]
		memset(ntt2, 0, sizeof ntt2);
		for (int i = 2; i < 6; i++)
			for (uint32_t j = 0; j < ntt[i]; j++) {
				register uint32_t U = wubuf.tbuilexpand3[i][j];
				if (!(U&F)) {
					U &= AC;
					register int cc = _popcnt32(U);
					if (!cc) return 3; // dead branch
					if (cc > 5) cc = 5;
					wubuf.tbuilexpand3[cc + 6][ntt2[cc]++] = U;
				}
			}
		if (ntt2[1]) {// send back the bit field status
			for (uint32_t j = 0; j < ntt2[1]; j++) 
				is1_tw3 |= wubuf.tbuilexpand3[7][j];
			
		}
		for (int i = 1; i < 6; i++)
			for (uint32_t j = 0; j < ntt2[i]; j++) {
				// clear redundancy
				register uint32_t U = wubuf.tbuilexpand3[i + 6][j],
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
				register uint32_t U = wubuf.tbuilexpand3[i][j],
					nu = ~U, x = 0;
				for (uint32_t i = 0; i < ntw3_2; i++)
					if (!(nu&tw3_2[i])) { x = 1; break; }
				if (x) continue;
				tw3_2[ntw3_2++] = U;
			}
	}

	return 0;
}
int GCHK::BuildFinalTable() {
#ifdef DEBUG5
	int locdiag = 0;
	if (p_cpt2g[5] == DEBUG5) {
		cout  << " BuildFinalTable() is1_tw3="<< is1_tw3 << endl;
		locdiag = 1;
	}
#endif
	if (is1_tw3) {// still some singles to assign
		if (scritb3.AddAssign(is1_tw3)) {
			cout << Char54out(is1_tw3) << "is1_tw3 in excess" << endl;
			return 1;
		}// then pack and filter the rest
		register uint32_t AC = scritb3.active, F = scritb3.assigned,
			n=0;
		for (uint32_t i = 0; i < ntw3_2; i++) {
			register uint32_t U = tw3_2[i];
			if (!(U & F)) {
				if (!(U &= AC))	return 1;
				tw3_2[n++] = U;
			}
		}
		ntw3_2=n;

	}
	int nass = _popcnt32(scritb3.assigned),ntoass= scritb3.nb3-nass;
	 ntaddgob3 = 0;
#ifdef DEBUG5
	 if (locdiag) {
		 cout << " ntoass=" << ntoass<< " ntw3_2="<< ntw3_2 << endl;
	 }
#endif
	 if( (!ntoass) && ntw3_2)return 1;
	if (!ntoass) {// enough clues assigned in b3
		if (ntw3_2) return 1;
		p_cpt2g[10 ]++;
		B3FinalNoExpand();//must be a valid b12 and valid "all" 
		return 0;
	}
	if(ntoass<6)p_cpt2g[10 + ntoass]++; else p_cpt2g[16]++;
	if (!ntw3_2) {//could be a valid 17 or not valid below 17
		if (!clean_valid_done) {// must be a valid b12 anyway to go
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) return 0;
		}
		if (!IsValidB3(scritb3.assigned)) {// a valid 17
			Out17(scritb3.assigned);
			return 0;
		}
		if (stopexpandb3) return 1;
		// now ntw3_2>0
	}
	if (ntoass == 1) {
		//cout << "do ntoass=1" << endl;
		//add can be any of active scritb3 hitting all uas
		register uint32_t wadd = scritb3.active;
		for (uint32_t i = 0; i < ntw3_2; i++)wadd &= tw3_2[i];
		if (!wadd) return 0;
		if (!clean_valid_done) {// must be a valid b12 anyway to go
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) return 0;
		}
		register int cell,
			Ass = scritb3.assigned;
		while (bitscanforward(cell, wadd)) {
			register uint32_t bit = 1 << cell;
			wadd ^= bit; //clear bit
			if (IsValidB3(Ass | bit)) {
				if (stopexpandb3) return 1;
				wadd &= anduab3;
			}
			else 	Out17(Ass | bit);
		}
		return 0;
	}

	if (ntoass == 2) {
		//cout << "do ntoass=2" << endl;
		uint32_t cell1, first = tw3_2 [0];// must exist checked upstream
		while (bitscanforward(cell1, first)) {
			register uint32_t bit = 1 << cell1;
			first ^= bit; //clear bit
			register uint32_t bf1 = scritb3.assigned | bit,cell;
			//add can be any of active scritb3 hitting all uas
			register uint32_t wadd = scritb3.active;
			for (uint32_t i = 1; i < ntw3_2; i++) {
				register uint32_t U= tw3_2[i];
				if(!(U&bf1))	wadd &= U;
			}
			if (!wadd) continue;
			if (!clean_valid_done) {// ust be a valid b12 anyway to go
				clean_valid_done = 1;
				if (IsValidB12(nclgo)) return 0;
			}
			while (bitscanforward(cell, wadd)) {
				register uint32_t bit = 1 << cell;
				wadd ^= bit; //clear bit
				if (IsValidB3(bf1 | bit)) {
					if (stopexpandb3) return 1;
					wadd &= anduab3;
				}
				else 	Out17(bf1 | bit);
			}
		}
		return 0;
	}
	// now 3 and more clues to add
	//cout << " more than 2 clues to add" << endl;
	if (ntw3_2 < 10) {
		ExpandB3Direct(ntoass);
		return 0;
	}
	//cout << " not direct" << endl;
	b3direct.Init();
	for (uint32_t i = 0; i < ntw3_2; i++) {
		//cout << Char27out(tw3_2[i]) << endl;
		b3direct.Add(tw3_2[i]);
	}
	//scritb3.Status(" not direct ");
	ExpandB3Vect(ntoass);
	return 0;
}

void GCHK::B3FinalNoExpand() {
	if (!clean_valid_done) {
		clean_valid_done = 1;
		if (IsValidB12(nclgo)) return ;
	}
	if(IsValidB3(scritb3.assigned)) return;
	//_________________________   valid 18 
	Out17(scritb3.assigned);
}


void GCHK::InitGoB3(uint64_t bf54,  uint32_t ncl) {// bands1+2 locked
#ifdef HAVEKNOWN
	if (okcheck == 2) {
		if ((bf54 & pk54) == bf54) {
			cout << Char54out(bf54) << " p_cpt2g[6]=" << p_cpt2g[6] << endl;
		}
	}
#endif
#ifdef DEBUG5
	int locdiag = 0;
	if (p_cpt2g[5] == DEBUG5) {
		cout << Char54out(bf54) << " bf init gob3   nclues=" << ncl<< endl;
		locdiag = 1;
	}
#endif
	nclgo = ncl;
	uint32_t 	nclb12 = ncl - nclues_step;	
	ncluesb3 = 18 - ncl;
	{
		uint32_t mg2 = g2_256[0].Apply(tcluesxy, nclb12),
			mg3 = g3_256.Apply(tcluesxy, nclb12);
		register uint64_t Bf = bf54;
		register uint32_t  vmini, Mg2 = mg2, Mg3 = mg3;
		scritb3.Init(ncl);
		for (int imini = 0; imini < 9; imini++, Mg2 >>= 3, Mg3 >>= 1) {
			if (!(vmini = Mg2 & 7))	if (Mg3 & 1)vmini = 8;
			if (vmini)scritb3.MinVmini(imini, vmini);
		}
		if (ncl == limb12)
			if ((int)scritb3.mincount > ncluesb3) return ;
		scritb3.SetStackf(Bf);// stackf; nmisss; critstack if nmiss=0
		if (scritb3.Stackx())return ;
		if (g2_256[1].nv) {
			uint32_t mg2_2 = g2_256[1].Apply(tcluesxy, nclb12) | mg2;
			if (mg2 != mg2_2) { // redo it, not a common situation
				Mg2 = mg2_2; Mg3 = mg3;
				scritb3.Init(ncl);
				for (int imini = 0; imini < 9; imini++, Mg2 >>= 3, Mg3 >>= 1) {
					if (!(vmini = Mg2 & 7))		if (Mg3 & 1)vmini = 8;
					if (vmini)scritb3.MinVmini(imini, vmini);
				}
				if (ncl == limb12)
					if ((int)scritb3.mincount > ncluesb3) return ;
				scritb3.SetStackf(Bf);// stackf; nmisss; critstack if nmiss=0
				if (scritb3.Stackx())return ;
			}
		}
	}
	p_cpt2g[8]++;

	myb12f = bf54;
	scritb3.AssignCritical();
	ntw3 = 0; // temporary before  split in/out
	g_256h.ApplyMore(tcluesxy, ncluesb12,
		tw3, ntw3, scritb3.assigned);

	{
		register uint32_t  critbf = scritb3.critbf,
			nmiss = scritb3.nmiss, n1 = 0, n2 = 0;
		if (!nmiss) {// first out is "dead"
			for (uint32_t i = 0; i < ntw3; i++)
				if (!(tw3[i] & critbf)) return ;
		}
		else if (nmiss == 1) {//use "and" as out
			register uint32_t andout = scritb3.active & (~critbf);
			for (uint32_t i = 0; i < ntw3; i++) {
				register uint32_t U = tw3[i];
				if (U & critbf) tw3[n1++] = U;
				else { andout &= U; n2 = 1; }
			}
			if (n2) {
				if (!andout)return ;// no possibility to add
				tw3[n1++] = andout; // dummy ua hitting all "out"
			}
			ntw3 = n1;
		}
	}
	if(BuildB2Table()) return;

	int nass= _popcnt32(scritb3.assigned);
#ifdef DEBUG5
	if (locdiag) {
		cout  << " nass=" << nass << endl;
		scritb3.Status(" after  BuildB2Table ");
	}
#endif
	if (nass == scritb3.nb3) {// already ok to go
		if (ntw3_2)return;
		cout << "cas saturé, vérifier de suite validité ==>>a faire" << endl;
		return;
	}
	//for (uint32_t i = 0; i < ntw3_2; i++)
	//	cout << Char27out(tw3_2[i]) << endl;
	if (is1_tw3) {
		cout << Char27out(is1_tw3) << " singles to assign " << endl;
		scritb3.Status(" before singles to assign  ");

	}
	BuildFinalTable();
}

void GCHK::CleanMoreUas(uint64_t bf, uint64_t ac, int ncl, MOREANDB & mabo) {
	p_cpt2g[64]++;
	uint64_t ua = mabo.tm[0];
	MOREANDB mabn;
#ifdef DEBUG5
	if (p_cpt2g[5] == DEBUG5) {
		cout << Char54out(bf) << "entry CleanMoreUas mab.ntm=" << mabo.ntm << " ncl=" << ncl << endl;
		cout << Char54out(ac) << "ac"  << endl;
		cout << Char54out(ua) << "ua" << endl;
	}
#endif	// apply ua and loop if more uas or call add process
	uint32_t cell;
	while (bitscanforward64(cell, ua)) {
		uint64_t bit = (uint64_t)1 << cell;
		ua ^= bit;
		tclues[ncl] = cell;
		int nclues = ncl + 1;
		ac ^= bit;// not active downstream
		uint64_t bf2 = bf | bit;

// myb12f nclues nluesb12 needed to call ApplyCountFilterB

		{// build MOREAND for the next cycle
			mabn.ntm = 0;
			for (uint64_t i = 1; i < mabo.ntm; i++) {
				register uint64_t R = mabo.tm[i];
				if(!(R&bit))		mabn.tm[mabn.ntm++] =R;
			}
			// add 		mbisvalid.ntm = 0;
			for (uint64_t i = 0; i < mbisvalid.ntm; i++) {
				register uint64_t R = mbisvalid.tm[i];
				if (!(R & bit))		mabn.tm[mabn.ntm++] = R;
			}

		}
		int nfree = limb12 - nclues;
#ifdef DEBUG5
		if (p_cpt2g[5] == DEBUG5) {
			cout << Char54out(bf2) << " bf2 mabn.ntm=" << mabn.ntm << " nclues="<<nclues << endl;
			//mabn.Dump();
		}
#endif
		if (!nfree) {
			if (mabn.ntm)continue;
			else {
				clean_valid_done = 0;			myacf = ac;
				mynclues = nclues;		ncluesb3 = 18 - mynclues;
				InitGoB3(bf2,nclues);
			}
		}
		else {
			if (!mabn.ntm) {
				if (!IsValidB12(nclues)) {// add one cell hitting all uas
					clean_valid_done = 1;
					InitGoB3(bf2,nclues);// try direct
					CheckValidBelow(bf2, myac);
				}
				else for (uint64_t i = 0; i < mbisvalid.ntm; i++) {
					register uint64_t R = mbisvalid.tm[i];
					if (!(R & bit))		mabn.tm[mabn.ntm++] = R;
				}

				nclues = ncl + 1;//restore nclues

			}
			if (nfree == 1) { // last step must hit all uas
				uint64_t a = mabn.GetAnd(ac);
				if (!a) { continue; }// not one clue

				//__ apply limits depending on the number of clues in bands 1+2
				uint64_t n1 = _popcnt64(bf2 & BIT_SET_27),
					n2 = nclues - n1; // 10 or 11 already there
#ifdef DEBUG5
				//if (p_cpt2g[5] == DEBUG5) {
					//cout << Char54out(bf2) << "CleanMoreUas n1= " << n1 << " n2=" << n2 << endl;
					//cout << Char54out(a) << " a mabn.ntm="<< mabn.ntm << endl;
					//mabn.Dump();
				//}
#endif
				int lim = (limb12 == 12) ? 6 : 7;
				if (n1 > lim || n2 > lim) continue;// no final ok
				if (n2 < 3 || n2>6)continue; // none of the 4 possibles
				if (n2 == 6)a &= BIT_SET_27;//can only add in band 1
				else if (n2 == 3)a &= ~BIT_SET_27;//must add in band 2
				else if (limb12 == 12 && n2 == 5)a &= ~BIT_SET_27;//must add in band 2
#ifdef DEBUG5
				if (p_cpt2g[5] == DEBUG5) {
					cout << Char54out(bf2) << "CleanMoreUas n1= " << n1 << " n2=" << n2 
						<< " nclues="<<nclues << endl;
					cout << Char54out(a) << " a2 mabn.ntm=" << mabn.ntm << " [5] "<< p_cpt2g[5] << endl;
					mabn.Dump();
				}
#endif
				if (!a) { continue; }// no expected valid

				int x, ncluer = nclues++;
				while (bitscanforward64(x, a)) {
					tclues[ncluer] = x;
					uint64_t bit = (uint64_t)1 << x;
					a ^= bit; //clear bit
					myb12f = bf2 | bit;
					InitGoB3(myb12f,nclues);
				}
				continue;
			}
			else CleanMoreUas(bf | bit, ac, nclues, mabn);
		}


	}
}
void GCHK::CheckValidBelow(uint64_t bf, uint64_t ac) {
	p_cpt2g[65]++;
#ifdef HAVEKNOWN
	int locdiag = 0;
	if (okcheck == 2) {
		if ((bf&pk54) == bf) {
			cout << Char54out(bf) << " start below p_cpt2g[7]=" << p_cpt2g[7] << endl;
			locdiag = 1;
		}
	}
#endif	
	int nclues = (int)_popcnt64(bf);
	{// check limit per band  11 clues 747 567 657 12 clues 666
		uint64_t n12 = _popcnt64(bf),
			n1 = _popcnt64(bf& BIT_SET_27),
			n2 = n12 - n1, lim = 8;
		//if (n1 > 7 || n2 > 8) return;// should not be here
		myacadd = myac_4 & ~bf;
		if (n12 ==11) {// 11 clues 747 567 657 12 clues 666
			if (n2 < 4 || n2>6) return;
			if (n1 == 7)myacadd = 0;
			else if (n2 == 6)myacadd &= BIT_SET_27;//add in b1
			else myacadd &= ~BIT_SET_27;//add in b2
		}
		else if (n12 == 10) {// add 1 or add 2 depending on limb12
			if (n2 > 6)myacadd = 0; // 288 378 
			else {// 73 46 64 55
				if (n2>5)myacadd &= BIT_SET_27;//add in b1
				else if (n2 <4)myacadd &= ~BIT_SET_27;//add in b2
			}
		}
		p_cpt2g[66]++;
#ifdef DEBUG5
		if (p_cpt2g[5] == DEBUG5) {
			cout << "CheckValidBelow call isvalid nclues="<<nclues 
				<< " n1=" << n1 << " n2=" << n2 << endl;
			cout << "tclues ";
			for (int i = 0; i < nclues; i++)cout << tclues[i] << " ";
			cout << endl;
			cout << Char54out(bf) << " bf" << endl;
			cout << Char54out(ac) << " ac" << endl;
			cout << Char54out(myacadd) << " myacadd" << endl;
		}
#endif
		if ((!clean_valid_done) &&IsValidB12(nclues)) {// add one cell hitting all uas
			int nfree = limb12 - nclues;
			if (!nfree) return;
			if (nfree == 1)if (!mbisvalid.SetAnd()) return;
			//if(myacadd)
			CleanMoreUas(bf, ac, nclues, mbisvalid);
			return;
		}
	}
	clean_valid_done = 1;
	p_cpt2g[67]++;

	//InitGoB3(bf);// try first direct with more clues in b3
	if (!myacadd) return;
	p_cpt2g[46]++;

	{//______________  build tadd from active cells
		ntadd = 0;
		int cell;
		register uint64_t V = myacadd;// still valid c2
		while (bitscanforward64(cell, V)) {
			V ^= (uint64_t)1 << cell;
			tadd[ntadd++] = cell;
		}
	}
	//cout << "\t\t" << p_cpt2g[46] << " " << p_cpt2g[5]
	//	<<" "<<ntadd<<" "<<nclues << endl;
	ExpandAddB1B2(bf);
}
void GCHK::ExpandAddB1B2Go(int step) {
	if (aigstop) return;
	p_cpt2g[68]++;
	myb12f = myb12add;
	{// clean redundancy
		uint64_t n12 = _popcnt64(myb12f), n2 = _popcnt64(myb12f >> 27),
			n1 = n12 - n2, n3 = 18 - n12;
		if (n3 < n1 || n3 < n2) return;
	}
	uint32_t *tgo = &tclues[mynclues];
	int nclues = mynclues + step;
	ncluesb12 = nclues - nclues_step;
#ifdef DEBUG5
	if (p_cpt2g[5] == DEBUG5) {
		cout << Char54out(myb12f) << "ExpandAddB1B2Go p_cpt2g[68]=" << p_cpt2g[68]
			<< " step=" << step << " clue=" << tgo[step] << endl;
	}
#endif

#ifdef HAVEKNOWN
	int locdiag = 0;
	if (myb12f == pk54) {
		cout << Char54out(myb12f) << " expected bands 1+2 in add"
			<< " p_cpt2g[26]=" << p_cpt2g[26] << endl;
		locdiag = 1;
	}
#endif	
	clean_valid_done = 0;
	InitGoB3(myb12f,nclues);
}
void GCHK::ExpandAddB1B2(uint64_t bf) {// add up to n cells
	mynclues = (uint32_t)_popcnt64(bf);
	uint32_t *tgo=&tclues[mynclues],	lim=17- mincluesb3- mynclues;

	struct SPB {
		uint64_t  all_previous_cells, icur;
	}spb[20], * s, * sn;
	s = spb;
	memset(s, 0, sizeof spb[0]);//	s->icur=0
	s->all_previous_cells = bf;// mode 54 previous cells
	//____________ here start the search nclues
next:
	if (s->icur >= ntadd) goto back;
	{
		uint64_t ispot = s - spb;
		register int cell = tadd[s->icur++];
		register uint64_t bit = (uint64_t)1 << cell;
		tgo[ispot] = cell;
		sn = s + 1; *sn = *s;
		sn->all_previous_cells = s->all_previous_cells | bit;
		myb12add = sn->all_previous_cells;
		ExpandAddB1B2Go((int)ispot + 1);// call the process for this 
		if (ispot >= lim) 	goto next;
		s++; // switch to next spot (icur+1)
		goto next;
	}
back:
	if (--s >= spb)goto next;

}



//============ clues in band 3 (no more clue in bands 1+2)
/*
int GCHK::VB12::GetsminF(uint32_t v27) {// mincoutn still active guas2
	register uint64_t F = gchk.myb12f;
	memset(&smin, 0, sizeof smin);
	ntg2ok = 0;
	//  guas2 added in this chunk
	for (uint32_t i = 0; i < gchk.gaddb3.nt2; i++) 
		if (!(F & gchk.gaddb3.t2[i].bf.u64[0]))
			v27 |= 1 << gchk.gaddb3.t2[i].bf.u32[2];
	for (register uint32_t i27 = 0,bit=1; i27 < 27; i27++,bit <<=1) {
		if (!(v27&bit)) continue;
		tg2ok[ntg2ok++] = i27;
		v27 |= bit;
		GCHK::SG2 w2 = gchk.tg2[i27];
		register uint32_t	bit9 = w2.bit9;
		smin.mini_bf3 |= smin.mini_bf2&bit9;
		smin.mini_bf2 |= smin.mini_bf1&bit9;
		smin.mini_bf1 |= bit9;
		smin.critbf |= w2.pat;
		smin.pairs27 |= 1 << i27;
	}
	smin.SetMincountG2();
	return smin.mincount;
}
inline void GCHK::VB12::ApplyBf2() {
	nclues = bfbf2 = 0;
	if (!smin.mini_bf2) return;
	for (int imini = 0, bit = 1; imini < 9; imini++, bit <<= 1)
		if (smin.mini_bf2&bit) {
			uint32_t Mask = 7 << (3 * imini), cell;
			smin.critbf &= (~Mask); // clear the minirow
			uint32_t bfcell = Mask & (~smin.pairs27);//bit common cell 
			bfbf2 |= bfcell;
			bitscanforward(cell, bfbf2);
			tclues[nclues++] = cell;
		}
}
inline void GCHK::VB12::BuildOf() {
	ntof = 0;
	register uint32_t F = smin.critbf | bfbf2; // assigned or critbf
	for (uint32_t i = 0; i <ntmore27; i++) {
		register uint32_t U =tmore27[i];
		if (!(U&F))tof[ntof++] = U;
		if (ntof > 20) return;//enough to check higher min
	}
	// same with band 3
	for (uint32_t i = 0; i < bax[2].nua; i++) {
		register uint32_t U = bax[2].tua[i] & BIT_SET_27;
		if (!(U&F))tof[ntof++] = U;
		if (ntof > 20) return; 
	}
}
*/


#define LIMADD 17


uint32_t GCHK::NoExpandB3(uint32_t cl0bf) {
	b3direct.nt = 0;
	if (!clean_valid_done) {
		clean_valid_done = 1;
		if (IsValidB12(nclgo)) return 0;
	}
	if (zhou[0].CallCheckB3(tclues, nclf, cl0bf)) {
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			register uint32_t ua = w.bf.u32[2];
			taddgob3[ntaddgob3++] = ua;
			b3direct.Add(ua);
			BuildGua(w, cc);
			gaddb3.Add(w, cc);//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			if (cc0 > 20 || cc > 8)continue;
			if (cc < 7 && cc0 < LIMADD) {
				if (cc < 3)p_cpt2g[35]++; else p_cpt2g[36]++;
				chunkh.Add128(w, cc);
			}
		}
		return 1;
	}
	//_________________________   valid 18 
	Out17(cl0bf);
	return 0;
}



void GCHK::ExpandB3Direct(int ntoass) {
	int locdiag = 0;

	p_cpt2g[38]++;
#ifdef DEBUG5
	if (p_cpt2g[5] == DEBUG5) {
		cout << "do ExpandB3Direct(ntoass) p_cpt2g[38]= " << p_cpt2g[38] << endl;
		locdiag = 1;
	}
#endif
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
#ifdef DEBUG5
	if (locdiag) {
		cout << "expand direct ntoass=" << ntoass << " lims1/2 " << limspot << " " << limm1
			<< " ntw3_2=" << ntw3_2 << endl;
		scritb3.Status("start expand direct");
	}
#endif
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
		if (sn->critb3.Addone(cell )) {
			if (locdiag) cout << " killed add one" << endl;
			goto next;
		}
#ifdef DEBUG5
		if (locdiag) {
			cout  << " " << cell << " " << ispot << " cell spot\t";
			cout << Char27out(sn->critb3.active) << " act after assign" << endl;
		}
#endif
		//sn->critb3.assigned |= bit;
		//__________________check for next uaor no more ua
		//sn->vect &= b3direct.vc[cell];
	}
	//if (locdiag) 	sn->critb3.Status("after assign");
	
	if (ispot < limm1) {// here max 16 clues never valid b3
		//if (locdiag) cout << " not last step s->indtw3="<< s->indtw3 << endl;
		register uint32_t F = sn->critb3.assigned;
		for (uint32_t i = s->indtw3 + 1; i < ntw3_2; i++) {
			register uint32_t U = tw3_2[i];
			//if (locdiag)cout << Char27out(U) << " next ua" << endl;
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
		p_cpt2g[19]++;
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
		}
		p_cpt2g[20]++;
		if (IsValidB3(sn->critb3.assigned)) {
			if (stopexpandb3) return;
			s->possible_cells &= anduab3;// same spot hitting  new ua
			sn->possible_cells = anduab3;
			s = sn;
			goto next;
		}
		cout << "bug seen valid 16 or less [5] "
			<< p_cpt2g[5] << " [19] " << p_cpt2g[19] << endl;
		aigstop = 1;
		return;
	}
	// this is the last step must hit all pending uas
	{ // find next cells hitting all uas
		int aig = 1;
		if (locdiag) {
			cout << Char27out(sn->critb3.assigned)
				<< " last step indtw3=" << s->indtw3 << endl;
			cout << Char27out(sn->critb3.active) << " act" << endl;
		}
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
		if (locdiag)cout << Char27out(andw) << "and clean_valid_done="<< clean_valid_done << endl;
		// no more ua or "and all uas" not empty
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (locdiag)cout << " checkvalid " << endl;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
			if (locdiag)cout << " back checkvalid " << endl;
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
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
			if (IsValidB3(F | bit)) {
				if (stopexpandb3) return;
				if (locdiag)cout << Char27out(anduab3) << " new ua anduab3" << endl;
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
	p_cpt2g[39]++;
	int locdiag = 0;
	//if (p_cpt2g[39] ==1) locdiag = 1;
	//if (p_cpt2g[5] ==393) locdiag = 1;
	if (locdiag)cout << "do ExpandB3Vect(ntoass) p_cpt2g[39]= " << p_cpt2g[39] << endl;
	uint64_t limspot = ntoass - 1, limm1 = limspot - 1;
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
	if (locdiag) {
		cout << " ntoass=" << ntoass << " lims1/2 " << limspot << " " << limm1 << endl;
		scritb3.Status("start expand vect");
	}
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
		if (sn->critb3.Addone(cell)) {
			if (locdiag) cout << " killed add one" << endl;
			goto next;
		}
		sn->vect &= b3direct.vc[cell];
		if (locdiag) {
			cout << Char27out(sn->critb3.assigned) << " " << cell << " " << ispot << " cell spot\t";
			cout << Char27out(sn->critb3.active) << " act after assign" << endl;
		}
		//sn->critb3.assigned |= bit;
		//__________________check for next uaor no more ua
		//sn->vect &= b3direct.vc[cell];
	}
	//if (locdiag) 	sn->critb3.Status("after assign");

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
		p_cpt2g[19]++;
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
		}
		p_cpt2g[20]++;
		if (IsValidB3(sn->critb3.assigned)) {
			if (stopexpandb3) return;
			s->possible_cells &= anduab3;// same spot hitting  new ua
			sn->possible_cells = anduab3;
			s = sn;
			goto next;
		}
		cout << "bug seen valid 16 or less [5] "
			<< p_cpt2g[5] << " [19] " << p_cpt2g[19] << endl;
		aigstop = 1;
		return;
	}
	// this is the last step must hit all pending uas
	{ // find next cells hitting all uas
		int aig = 1;
		if (locdiag) {
			cout << Char27out(sn->critb3.assigned)
				<< " last step "  << endl;
			cout << Char27out(sn->critb3.active) << " act stacks " 
				<< sn->critb3.stackf[0] << sn->critb3.stackf[1] 
				<< sn->critb3.stackf[2] << endl;
		}
		register uint32_t andw = sn->critb3.active;
		register uint32_t F = sn->critb3.assigned;
		if (ntw3_2) {// likely quicker 0
			for (uint32_t ia = 0; ia < ntw3_2; ia++) {
				register uint32_t U = tw3_2[ia];
				if (!(U & F)) {
					if (locdiag)cout << Char27out(U) << " from ia= " << ia << endl;
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
				if (locdiag)cout << Char27out(b3direct.tua[ir]) << " solde " << endl;
				if (!(andw &= b3direct.tua[ir])) goto next;
				aig = 0;
			}
		}
		if (locdiag)cout << Char27out(andw) << "and" << endl;
		// no more ua or "and all uas" not empty
		if (!clean_valid_done) {
			clean_valid_done = 1;
			//cout << " checkvalid " << endl;
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
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
			if (IsValidB3(F | bit)) {
				if (stopexpandb3) return;
				if (locdiag)cout << Char27out(anduab3) << " new ua anduab3" << endl;
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
	int locdiag = 0;
	//if (p_cpt2g[38] < 20) locdiag = 1;
	if (zhou[0].CallCheckB3(tclues, nclf, bf)) {
		anduab3 = BIT_SET_27;
		stopexpandb3 = 0;
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			if (locdiag) {
				cout << Char27out(w.bf.u32[2]) << " ";
				cout << Char54out(w.bf.u64[0]) << " added b3" << endl;
			}
			int cc = _popcnt32(w.bf.u32[2]);
			if (!cc) {
				cout << " bug no ua b3 IsValidB3 p_cpt2g[5]" << p_cpt2g[5]<< endl;
				aigstop = 1;
				return 1;
			}
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			register uint32_t ua = w.bf.u32[2];
			tw3_2[ntw3_2++] = ua;
			taddgob3[ntaddgob3++] = ua;// more for a given b1+2
			anduab3 &= ua;
			BuildGua(w, cc);
			guash3.Add128(w, cc);//clean level new if possible
			g_256h.Add128(w, cc);//clean level new if possible
			if (cc < 4) {
				if ((cc == 2 && scritb3.CanNotAddi27(w.bf.u32[2])) ||
					(cc == 3 && scritb3.CanNotAddi9(w.bf.u32[2]))) {
					stopexpandb3 = 1;
					/*
						cout << " forced stopexpand b3 cc="<<cc <<" cc0 ="<< cc0 << endl;
					if (cc == 2) {
						cout << Char27out(tg2[w.bf.u32[2]].pat) << " ";
						cout << Char54out(w.bf.u64[0]) << " forced stop" << endl;
					}
					if (cc == 3) {
						cout << Char27out(7<<(3*w.bf.u32[2])) << " ";
						cout << Char54out(w.bf.u64[0]) << " forced stop" << endl;
					}
				*/
				}
			}
			if (cc0 > 20 || cc > 6)continue;
			chunkh.Add128(w, cc);
			if (cc < 3)p_cpt2g[35]++; else p_cpt2g[36]++;
			/*
			if (cc == 2) {
				cout << Char27out(tg2[w.bf.u32[2]].pat) << " ";
				cout << Char54out(w.bf.u64[0]) << " added validb3"
					<< " [4] " << p_cpt2g[4]
					<< " [5] " << p_cpt2g[5] << " [35] " << p_cpt2g[35] << endl;
			}
			*/
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
	//zh_g2.isfalse_on = -1;
	//ImageCandidats();
	int ir = FullUpdate();
	if (nogo==1) { ImageCandidats(); return 0; }// test
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
	R2 &= ~R3; // now true singles
	if (!R2) {
		R3 &= ~R4; // now true singles
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
	if (zh_g.diag) {
		cout << "index=" << index << endl;
		ImageCandidats();
	}
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
				if (zh_g.diag) {
					cout << "zhgxn.nua=" << zhgxn.nua << endl;
					wua.Print3(" ");
				}
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
	// mapping of the output on band_order[ib] 
	uint32_t map[81];
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 27; j++)
		map[27 * i + j] = 27 * band_order[i] + j;
	char ws[82];

	uint32_t tcf[40], ntcf = 0;
	for (int i = 0; i < nclf ; i++) {
		tcf[ntcf++] = map[tclues[i]];
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		tcf[ntcf++] = map[54 + i] ;
	strcpy(ws, empty_puzzle);
	for (uint32_t i= 0; i < ntcf; i++) {
		int cell = tcf[i];
		ws[cell] = ze[cell] ;
	}
	cout << ws << " one sol  size " <<  ntcf 
		<< " [3]=" << p_cpt2g[3] << " [4]=" << p_cpt2g[4] << " [5]=" << p_cpt2g[5] << endl;
 

}




#define DEBUG0 1
//#define DEBUG3 4
//#define DEBUG4 167
//#define DEBUG5 11829
//#define DEBUG6 46
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
	cout << " minb1=" << minb1 << " minb2=" << minb2 << endl;
#ifdef HAVEKNOWN
	for (int i = 0; i < 81; i++) cout << ze[i + 82];
	cout << " known reshaped" << endl;
#endif
#endif

	//__________________________ start uas search
	UaCollector();
	p_cpt2g[1] = tuasb12.nua;
	p_cpt2g[2] = chunkh.GetC2Count();
	BuildVectorsForExpand4B12();
	Expand4B12();
#ifdef TEST_ON
	tuasb12.Stats();
	chunkh.Stats();
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
		Do_phase2(t4_to_expand[i]);
		if (aigstop) 			return;
	}
}


//______  process a "no more uas bands 1+2
void GCHK::Do_phase2(T4_TO_EXPAND w) {
	p_cpt2g[3]++;
	{	register uint64_t bf54 = w.bf, ac54 = w.active;
	myac_4 = ac54;
#ifdef DEBUG3
	if (p_cpt2g[3] > DEBUG3) return;
	//if (p_cpt2g[3] >= DEBUG3 -3) 
	if (p_cpt2g[3] == DEBUG3) {
		//cout << Char54out(bf54) << " p_cpt2g[3]= " << p_cpt2g[3] << endl;
		//cout << Char54out(ac54) << " ac54 " << endl;
		//for (uint32_t i = 0; i < tuasb12.nua; i++) {
			//register uint64_t R = tuasb12.tua[i];
			//register uint64_t cc = _popcnt64(R);
			//if (cc == 12)cout << Char54out(R)  << endl;
		//}
	}

#endif

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
	//morev2a.Init(); morev2b.Init(); // now in vector
	//morev2c.Init(); morev2d.Init(); // now in vector
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
#ifdef DEBUG3
		if (p_cpt2g[3] == DEBUG3) {
			//cout << " tua4 end pick up sort" << endl;
			//ntua4 = 0;
			//uint32_t lim = 0;
			//for (int i = 0; i < 30; i++)if (ntt[i]) {
				//register uint64_t * tj = wubuf.tdophase2[i];
				//for (uint64_t j = 0; j < ntt[i]; j++) {
					//register uint64_t R = tj[j];
					//cout << Char54out(R) << " i=" << i << " n=" << ntua4++ << " "
						//<< _popcnt64(R) << endl;
				//}
			//}
		}

#endif	
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
#ifdef DEBUG3
	if (p_cpt2g[3] == DEBUG3) {
		//cout <<  " tua4 end pick up ntua4="<<ntua4  << endl;
		//for (uint32_t i = 0; i < ntua4; i++) {
			//register uint64_t R = tua4[i];
			//register uint64_t cc = _popcnt64(R);
			//cout << Char54out(R) <<" " << cc << " " <<i<< endl;
		//}
	}

#endif

	{
		register int wnua = ntua4;
		memset(v12_5_v0, 0, sizeof v12_5_v0);
		for (int i = 0; i < 64; i++) {//max 64*64=4096
			if (wnua >= 64) {	v12_5_v0[i] = ~0; wnua -= 64;	}
			else {//0 63 last
				v12_5_v0[i] = ((uint64_t)1 << wnua) - 1;// 0 if wnua=0
				break;
			}
		}
	}
	nv12_5 = (ntua4 + 63) >> 6; // number of vectors size 64
	memset(v12_5_c, 255, sizeof v12_5_c) ;// each vector 8 bytes
	for (uint32_t i = 0; i < ntua4; i++) {
		register uint64_t R54 = tua4[i];
		register uint32_t cell, bloc = i >> 6, ir = i - 64 * bloc;
		uint64_t bit = (uint64_t)1 << ir;
		while (bitscanforward64(cell, R54)) {
			R54 ^= (uint64_t)1 << cell; //clear bit
			v12_5_c[bloc][cell] ^= bit;
		}
	}
#ifdef DEBUG3
	if (p_cpt2g[3] == DEBUG3) {
		//cout << " vectors start  64_128 nv12_5=" << nv12_5 << endl;
		//for (uint32_t i = 1; i < nv12_5; i++) {
			//cout << Char64out(v12_5_v0[i])<<" vect "<< endl;
			//uint64_t * vc = v12_5_c[i];
			//for (int j = 0; j < 54; j++) if (vc[j] != ~0) {
				//cout << Char64out(vc[j] ) << " j= "<<j << endl;
			//}
			//register uint64_t R = tua4[i];
			//register uint64_t cc = _popcnt64(R);
			//cout << Char54out(R) << " " << cc << " " << i << endl;
			//break;
		//}
	}

#endif
	Do_phase2Expand(bf54, ac54);
	}
}
void  GCHK::Do_phase2Expand(uint64_t bf, uint64_t ace) {// mode 54 not 2x
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
	uint64_t limspot = 17 - nphase1 - mincluesb3-1, // 11 or 12 clues
		limclean = 5 - nphase1;// set to more if low {min b1+b2}
	if (minb1b2 < 9)limclean = 6 - nphase1;
	uint64_t andvalid = ~0, orvalid = 0;
	// _______________ expand to have all  minimal < n clues
	struct SPB {// spots to find starts to extract uas 3 digits
		uint64_t V64, possible_cells, all_previous_cells, active_cells;
	}spb[13];
	register SPB  *s, *sn;
	register uint64_t ispot;
	register int cell;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->all_previous_cells = bf;
	s->active_cells = ace;
	s->possible_cells = tua4[0];
	s->V64 = v12_5_v0[0];

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
			sn->V64 &= v12_5_c[0][cell];
		}
		if (sn->V64) {// next ua
			register uint32_t ir;
			bitscanforward64(ir, sn->V64);//relative index first active
			register uint64_t Ru = tua4[ir] & s->active_cells;
			if (!Ru)goto next;
			//____if last special process
			if (ispot == limspot) {// 4 for 11 ; 5 for 12
				{// "and"" remaining uas must all be hit
					register uint64_t V = sn->V64;
					V ^= (uint64_t)1 << ir; //clear first ua
					while (bitscanforward64(ir, V)) {
						V ^= (uint64_t)1 << ir;
						if (!(Ru &= tua4[ir])) goto next;;
					}
				}
				//if (!Ru)goto next;
				// clear potential redundancy
				{//__________ check count per band assumed nphase1=5
					register uint64_t n12 = nphase1 + 1 + ispot,  
						n1 = _popcnt64(sn->all_previous_cells & BIT_SET_27),
						n2=10-n1; // 10 or 11 already there
					if (limspot == 5) {//12 clues only 666
						n2++;
						if (n1 != 6)Ru &= BIT_SET_27;// can only be band1
						if (n2 != 6)Ru &= ~BIT_SET_27;// can only be band2
					}
					else {// 11 clues can  be  747 567 657
						if (n1 > 7 || n2 > 6)goto next; 
						if (n1 ==7)Ru &= ~BIT_SET_27;// only  band2
						if (n2 == 6)Ru &= BIT_SET_27;// only  band1
					}
				}
				if (!Ru)goto next;
				// now all cells can be added
				{
					register uint64_t Bf= sn->all_previous_cells;
					while (bitscanforward64(cell, Ru)) {
						register uint64_t bit = (uint64_t)1 << cell
							,Bff=Bf|bit;
						Ru ^= bit;
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
	}
	{// possible valid band 1+2 to store
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
		if (ispot != limspot+1)*pbufvalid++ = s->active_cells;
	}
	goto next;
back:
	if (--s >= spb)goto next;
	if (pbufvalid > bufvalid) // if limite do
		CleanBufferAndGo(andvalid, orvalid);
}

void GCHK::CleanBufferAndGo(uint64_t andvalid, uint64_t orvalid) {
	if (aigstop) return;
	p_cpt2g[4]++;
	uint64_t *myend = pbufvalid, *pw = bufvalid;
	limb12 = 18 - mincluesb3;
	pbufvalid = bufvalid;
	myandvalid = andvalid;
	myorvalid = orvalid;
	//	if (1) return; // check min phase 1 gen valids bands 1+2
	uint64_t nn8 = (myend - bufvalid),
		mincb2 = _popcnt64(andvalid >> 27);

#ifdef DEBUG3
#ifndef DEBUG4
	if (p_cpt2g[3] == DEBUG3)
		cout << Char54out(andvalid) << " CleanBufferAndGo p_cpt2g[4] debug " << p_cpt2g[4] << " size" << nn8 << endl;
#endif
#endif
#ifdef HAVEKNOWN
	if ((andvalid&pk54) == andvalid) {
		cout << Char54out(andvalid) << " p_cpt2g[4]=" << p_cpt2g[4] << endl;
	}
#endif
#ifdef DEBUG4
	if (p_cpt2g[4] > DEBUG4)return;
	if (p_cpt2g[4] == DEBUG4)
		cout << Char54out(andvalid) << " andv  p_cpt2g[4] " << p_cpt2g[4] << " size" << nn8 << endl;
#endif

	{//____ load the common clues
		register uint64_t R = andvalid;
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
		register uint32_t * tcl = tclues, ntcl = nclues_step;
		for (uint32_t i = 1; i < nv12_5; i++) {
			register uint64_t V = v12_5_v0[i],* Vc = v12_5_c[i];
			uint32_t x;
			//cout << Char64out(V) << " v12_5_v0[i] " << endl;
			for (uint32_t j = 0; j < ntcl; j++) {
				//cout << Char64out(Vc[tcl[j]]) << " vc cl="<< tcl[j] << endl;
				V &= Vc[tcl[j]];
			}
			while (bitscanforward64(x,V)) {
				V ^= (uint64_t)1 << x; //clear bit
				register uint64_t U = tua4[x + 64 * i];// not if free in b12  &orvalid;
				//if (!U) return;
				m64vh.Add54(U);
			}
		}
		uint32_t mm = m64vh.nmv * 128 + m64vh.mv[m64vh.nmv].nt;
		if (mm > p_cpt2g[49])p_cpt2g[49] = mm;
	}
	moreand.ntm = 0;
	BuildReducedGuas();
#ifdef DEBUG4
	if (p_cpt2g[4] == DEBUG4) {
		//guash.DumpPack();
		guash2.Pack(guash);
		guash2.DoVkill();
		guash2.DumpPack();

	}
#endif
	gaddb3.Init();
	ntuaclean = 0;// checking uas bands 1+2 added in this cycle

	//___________________________________ potential valid bands +2
	while (pw < myend) {
		if (aigstop) return;
		p_cpt2g[5]++;
		myb12 = *pw++;
		int ncl = (int)_popcnt64(myb12);
		if (ncl < limb12)myac = *pw++;		else myac = 0;
#ifdef HAVEKNOWN
		int locdiag = 0;
		if ((myb12&pk54) == myb12) {
			cout << Char54out(myb12) << " p_cpt2g[5]=" << p_cpt2g[5]
				<< " ncl=" << ncl << " limb12=" << limb12 << endl;
			cout << Char54out(myac) << " still active " << _popcnt64(myac) << endl;
			locdiag = 1;
		}
#endif
		{//____ load the remaining  clues
			ncluesb12 = 0;
			register uint64_t R = myb12 ^ andvalid;
			uint32_t cell;
			while (bitscanforward64(cell, R)) {
				R ^= (uint64_t)1 << cell; //clear bit
				tcluesxy[ncluesb12++] = cell;
			}
			nclues = nclues_step + ncluesb12;
		}
#ifdef DEBUG5
		if (p_cpt2g[5] > DEBUG5) return;
		if (p_cpt2g[5] == DEBUG5) {
			cout<<endl<<endl << Char54out(myb12) << " myb12  (4 ) " << p_cpt2g[4] 
				<<" [5] "<< p_cpt2g[5] << endl;
			cout << "tclues xy";
			for (int i = 0; i < ncluesb12; i++)cout << tcluesxy[i] << " ";
			cout << " m64vh.nmv=" << m64vh.nmv << endl;
			cout << " ncl=" << ncl << " limb12=" << limb12 << endl;
		}
#endif
		clean_valid_done = 0;
		if (ncl == limb12) {
			if (m64vh.ApplyXY(tcluesxy, ncluesb12)) continue;
			InitGoB3(myb12);
		}
		else CleanBelow();
	}
}
//_____ process a chunk of potential valid bands 1+2

void GCHK::BuildReducedGuas() {
	chunkh.ApplyClean(tclues, nclues_step);
	guash.Init();
	// _______________________________________extract guas2
	guash.g2 = 0;
	for (uint32_t ic2 = 0; ic2 <= chunkh.ic2; ic2++) {
		CHUNK3B &w = chunkh.c2[ic2];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			uint32_t i27 = w.tu27[ir];
			guash.AddG2(tg2[i27].pat, i27, w.tu54[ir]);
		}
	}
	// ________________________________________extract guas3
	guash.g3 = 0;
	for (uint32_t ic3 = 0; ic3 <= chunkh.ic3; ic3++) {
		CHUNK3B &w = chunkh.c3[ic3];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			uint32_t i9 = w.tu27[ir];
			guash.AddG3(	i9, w.tu54[ir]);
		}
	}
	// _____________________________________________extract guas4
	for (uint32_t ic4 = 0; ic4 <= chunkh.ic4; ic4++) {
		CHUNK3B &w = chunkh.c4[ic4];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			guash.AddG4(w.tu27[ir], w.tu54[ir]);
		}
	}
	guash.nstartx = guash.nbuf;
	// _____________________________________________extract guas5
	for (uint32_t ic5 = 0; ic5 <= chunkh.ic5; ic5++) {
		CHUNK3B &w = chunkh.c5[ic5];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			guash.AddG5(w.tu27[ir], w.tu54[ir]);
		}
	}
	// _____________________________________________extract guasmore
	for (uint32_t icm = 0; icm <= chunkh.icmore; icm++) {
		CHUNK3B &w = chunkh.cmore[icm];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			guash.AddGm(w.tu27[ir], w.tu54[ir]);
		}
	}
	// pack the results 
	guash2.Pack(guash);
	guash2.DoVkill();
	if (p_cpt2g[4] == 1) {
		//cout << "analyse cas 1 changement mode stockage " << endl;
		//guash.DumpPack();
		//guash2.DumpPack();
	}

}
inline void TUASB12::Add(uint64_t u, uint64_t cc) {
	if (cc < 18) {
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


int GCHK::IsValidB12() {
	p_cpt2g[15]++;
	mbisvalid.ntm = 0;
	if (zh2b[1].IsValid(tclues, nclues)) {
		p_cpt2g[30]++;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i],
				cc = _popcnt64(ua),
				ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			//if (ntuaclean < 4096) tuaclean[ntuaclean++] = ua54;
			tuasb12.Add(ua54, cc);  // global level
			m64vh.Add54(ua54); // 5 clues  level
			AddTua4(ua54);// chunk level
			mbisvalid.tm[mbisvalid.ntm++] = ua54;
			if(cc==12)		cout << Char54out(ua54)	<< " add12 p_cpt2g[4]=" 
				<< p_cpt2g[4] 	 <<  endl;
			//if (moreand.ntm < 500)moreand.tm[moreand.ntm++] = ua54;
		}
	}
	return (int)mbisvalid.ntm;
}
void GCHK::InitGoB3(uint64_t bf54 ) {// bands1+2 locked
	p_cpt2g[6]++;
	return;
#ifdef HAVEKNOWN
	int locdiag = 0;
	if (okcheck == 2) {
		if ((bf54&pk54) == bf54) {
			cout << Char54out(bf54) << " p_cpt2g[6]=" << p_cpt2g[6] << endl;
			locdiag = 1;
		}
	}
	if (okcheck == 3) {
		cout << Char54out(bf54) << " diag5 en coursp_cpt2g[6]=" << p_cpt2g[6] << endl;
		locdiag = 1;
	}
#endif

	myb12f = bf54;
	// get the active GUAs2 in vf then active guas2
	uint32_t * tgo = &tclues[nclues_step], ntgo = nclues - nclues_step;
	// 4.3 s to 13.8 s
	register uint64_t R = guash2.vkill[0].v0;
	{
		gguas2 = 0;
		register uint64_t bit64 = 1;
		for (uint32_t i = 0; i < ntgo; i++)
			R &= guash2.vkill[0].vc[tgo[i]];
		register uint32_t i = 0, bit27 = 1;
		for (; i < guash2.nc2i; i++,bit64<<=1) {
			if (R&bit64) {// active gua2 
				register GUASHORT & w = guash2.guashort[i];
				register uint32_t  bit27 = 1 << w.igua2_gua3;
				if (w.istart == w.iend)gguas2 |= bit27;
				else {// more than one ua first not hit ok
					for (uint32_t j = w.istart; j < w.iend; j++)
						if (!(bf54 & guash2.shortbuf[j]))
							{	gguas2 |= bit27; break; }
				}
			}
		}
	}
	//guash2.ApplyClean0(tgo, ntgo, myb12f);
	//gguas2 = guash2.g2;
	if (1) return;


#ifdef DEBUG5 
	int locdiag = 0;
	if (p_cpt2g[5] == DEBUG5) {
		locdiag = 1;
		cout << Char54out(myb12) << " clean p_cpt2g[6] debug " << p_cpt2g[6] << endl;
		cout << Char27out(gguas2) << " gguas2" << endl;
	}
#endif

	int nmin = svb12.GetsminF(gguas2);
	ncluesb3 = 18 - nclues;
	nmiss = ncluesb3 - nmin;
	//____________ passing the first filter nclues band 3 see Guas3
#ifdef HAVEKNOWN
	if (locdiag) {
		cout << Char54out(bf54) << " nmin=" << nmin << "nmiss=" << nmiss << endl;
		svb12.smin.Status(" guas2 only");
	}
#endif
	if (nmiss < 0) return;

	// ___ check here the limit in stack here 11 or 12 clues
	uint32_t min1 = svb12.smin.mini2all, min2 = svb12.smin.mini_bf3,
		mins1= _popcnt32(min1 & 0111)+ _popcnt32(min2 & 0111),
		mins2 = _popcnt32(min1 & 0222) + _popcnt32(min2 & 0222),
		mins3 = _popcnt32(min1 & 0444) + _popcnt32(min2 & 0444);
	uint64_t ns1 = _popcnt64(band3xBM[3].u64[0] & bf54),
		ns2 = _popcnt64(band3xBM[4].u64[0] & bf54),
		ns3=nclues-ns1-ns2,
		nt1=mins1+ns1,nt2=mins2+ns2,nt3=mins3+ns3;
	if (nt1 > 7 || nt2 > 7 || nt3 > 7)return;
	if (nclues == 12 &&	(nt1 > 6 ||nt2 > 6 || nt3 > 6))
		return;


	guash2.GetG3(myb12f);
	gguas3 = guash2.g3;

	gaddb3.Get3(myb12f, gguas3);// add guas3 from this add loop
	gguas3 &= ~svb12.smin.mini2all;// clean redundancy
	if (gguas3) {
		register int n3= _popcnt32(gguas3);
		if (nmiss < n3) return;
		nt1 += _popcnt32(gguas3 & 0111);
		nt2 += _popcnt32(gguas3 & 0222);
		nt3 += _popcnt32(gguas3 & 0444);
		if (nt1 > 7 || nt2 > 7 || nt3 > 7)return;
		if (nclues == 12 && (nt1 > 6 || nt2 > 6 || nt3 > 6))return;
		svb12.smin.mincount += n3;
		nmin += n3;
		nmiss -= n3;
		svb12.smin.AddMincountG3(gguas3);
	}
	if (++p_cpt2g[40] < 10) {
		cout << Char54out(bf54) << " stack debug nclues=" << nclues
			<< " miss " << nmiss << " " << p_cpt2g[6] << endl;
		svb12.smin.Status("after triplets");
		cout << "bands12 " << ns1 << ns2 << ns3 << endl;
		cout << "min " << mins1 << mins2 << mins3 << endl;
	}
	return;
#ifdef HAVEKNOWN
	if (locdiag && gguas3) 
		cout  << "after guas3 final  nmin=" << nmin << "nmiss=" << nmiss << endl;
#endif

	guash2.GetB12CX(tgo, ntgo, myb12f);

	gaddb3.GetMore(myb12f, guash2.t12, guash2.nt12);
	svb12.AttachMore(guash2.t12, guash2.nt12);


	// add more  from this add loop
#ifdef DEBUG5
//#define DEBUG14 5
	if (locdiag) {
		cout << Char54out(bf54) << " nmin=" << nmin << "nmiss=" << nmiss << endl;
		svb12.smin.Status(" call gob3 ");
	}
#endif
#ifdef DEBUG14
	if(DEBUG14 ==p_cpt2g[14]) {
		cout << Char54out(bf54) << " nmin=" << nmin << "nmiss=" << nmiss << endl;
		svb12.smin.Status(" call gob3 ");
	}
#endif
	GoB3(nclues, svb12);
}

void GCHK::CleanBelow() {
	if (1) return;//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	p_cpt2g[7]++;	// not limit can add clues in bands 1+2
	mab.ntm = 0;
	m64vh.Extract(tcluesxy, ncluesb12, mab.tm, mab.ntm);
	int nfree = limb12 - nclues;
#ifdef DEBUG5
	if (p_cpt2g[5] == DEBUG5) {
		cout << " mab.ntm=" << mab.ntm << " nfree="<<nfree << endl;
	}
#endif
	if (nfree == 1 && mab.ntm) {//11 clues 747 567 657 12 clues 666
		uint64_t a = mab.GetAnd(myac);
		if (!a) { p_cpt2g[51]++; return; }// not one clue
		//__ apply limits depending on the number of clues in bands 1+2
		uint64_t n1 = _popcnt64(myb12 & BIT_SET_27),
			n2 = nclues - n1; // 10 or 11 already there
		if (n2 < 3 || n2>6)return; // none of the 4 possibles
		if (n2 == 6)a &= BIT_SET_27;//can only add in band 1
		else if (n2 == 3)a &= ~BIT_SET_27;//must add in band 2
		else if (limb12 == 12 && n2 == 5)a &= ~BIT_SET_27;//must add in band 2
		if (!a) { p_cpt2g[51]++; return; }// no expected valid
		int x, ncluer = nclues++;
		while (bitscanforward64(x, a)) {
			tclues[ncluer] = x;
			uint64_t bit = (uint64_t)1 << x;
			a ^= bit; //clear bit
			InitGoB3(myb12 | bit);
		}
		return;
	}
	p_cpt2g[17]++;
	if (mab.ntm) 	CleanMoreUas(myb12, myac, nclues, mab);
	else CheckValidBelow(myb12, myac);
}
void GCHK::CleanMoreUas(uint64_t bf, uint64_t ac, int ncl, MOREANDB & mabo) {
	p_cpt2g[64]++;
	uint64_t ua = mabo.tm[0];
	MOREANDB mabn;
	// apply ua and loop if more uas or call add process
	uint32_t cell;
	while (bitscanforward64(cell, ua)) {
		uint64_t bit = (uint64_t)1 << cell;
		ua ^= bit;
		tclues[ncl] = cell;
		nclues = ncl + 1;
		ac ^= bit;// not active downstream
		{// build MOREAND for the next cycle
			mabn.ntm = 0;
			for (uint64_t i = 1; i < mabo.ntm; i++) {
				register uint64_t R = mabo.tm[i];
				if(!(R&bit))		mabn.tm[mabn.ntm++] =R;
			}
		}
		int nfree = limb12 - nclues;
		if (!nfree) {
			if (mabn.ntm)continue;
			else {
				clean_valid_done = 0;			myacf = ac;
				mynclues = nclues;		ncluesb3 = 18 - mynclues;
				InitGoB3(bf | bit);
			}
		}
		else if (nfree == 1 && mabn.ntm) { // last step must hit all uas
			uint64_t a = mabn.GetAnd(ac);
			if (!a) { p_cpt2g[51]++; return; }// not one clue
			//__ apply limits depending on the number of clues in bands 1+2
			uint64_t n1 = _popcnt64((bf | bit) & BIT_SET_27),
				n2 = nclues - n1; // 10 or 11 already there
			if (n2 < 3 || n2>6)return; // none of the 4 possibles
			if (n2 == 6)a &= BIT_SET_27;//can only add in band 1
			else if (n2 == 3)a &= ~BIT_SET_27;//must add in band 2
			else if (limb12 == 12 && n2 == 5)a &= ~BIT_SET_27;//must add in band 2
			if (!a) { p_cpt2g[51]++; return; }// no expected valid
			int x, ncluer = nclues++;
			while (bitscanforward64(x, a)) {
				tclues[ncluer] = x;
				uint64_t bit = (uint64_t)1 << x;
				a ^= bit; //clear bit
				InitGoB3(myb12 | bit);
			}
			return;
		}
		else {
			if (mabn.ntm)CleanMoreUas(bf | bit, ac, nclues, mabn);
			else  CheckValidBelow(bf | bit, ac);
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
			cout << "CheckValidBelow call isvalid nclues="<<nclues << endl;
			cout << "tclues ";
			for (int i = 0; i < nclues; i++)cout << tclues[i] << " ";
			cout << endl;
			cout << Char54out(bf) << " bf" << endl;
			cout << Char54out(ac) << " ac" << endl;
		}
#endif
		if (IsValidB12()) {// add one cell hitting all uas
			int nfree = limb12 - nclues;
			if (!nfree) return;
			if (nfree == 1)if (!mbisvalid.SetAnd()) return;
			if(myacadd)CleanMoreUas(bf, ac, nclues, mbisvalid);
			return;
		}
	}
	clean_valid_done = 1;
	p_cpt2g[67]++;

	InitGoB3(bf);// try first direct with more clues in b3
	if (!myacadd) return;
	p_cpt2g[46]++;
	p_cpt2g[47]+= _popcnt64(myacadd );

	uint32_t * tgo = &tclues[nclues_step], ntgo = nclues - nclues_step;
	guash2add.BuidAdd(guash2, tgo, ntgo, myb12f, myacadd);
	guash2add.DoVkill();
#ifdef DEBUG5
	if (p_cpt2g[5] == DEBUG5) {
		cout << "isvalid go b3" << endl;
		guash2.Status("guash2");
		guash2add.Status("guash2add");
		guash2add.DumpPack();
		//zh2b[1].ImageCandidats();
		//return;
	}
#endif
	{//______________  build tadd from active cells
		ntadd = 0;
		int cell;
		register uint64_t V = myacadd;// still valid c2
		while (bitscanforward64(cell, V)) {
			V ^= (uint64_t)1 << cell;
			tadd[ntadd++] = cell;
		}
	}
	ExpandAddB1B2(bf);
}
void GCHK::ExpandAddB1B2Go(int step) {
	if (1) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if (aigstop) return;
	p_cpt2g[68]++;
	myb12f = myb12add;
	{// clean redundancy
		uint64_t n12 = _popcnt64(myb12f), n2 = _popcnt64(myb12f >> 27),
			n1 = n12 - n2, n3 = 18 - n12;
		if (n3 < n1 || n3 < n2) return;
	}
	uint32_t *tgo = &tclues[mynclues];
	nclues = mynclues + step;
	// get the active GUAs2 in vf then active guas2
	guash2add.ApplyClean0(tgo, step, myb12f);
	gguas2 = guash2add.g2;

#ifdef DEBUG5
	if (p_cpt2g[5] == DEBUG5) {
		cout << Char54out(myb12f) << "ExpandAddB1B2Go p_cpt2g[68]=" << p_cpt2g[68]
			<< " step=" << step << " clue=" << tgo[step] << endl;
		if ( p_cpt2g[68] == 1) {
			cout << Char27out(gguas2) << "  gguas2 limb12=" << limb12 << endl;
		}
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
	int nmin = svb12.GetsminF(gguas2);
	ncluesb3 = 18 - nclues;
	nmiss = ncluesb3 - nmin;
#ifdef HAVEKNOWN
	if (locdiag) {
		svb12.smin.Status(" after guas2 ");
		cout << Char27out(gguas2) << " gguas2 from chunkadd nmiss="<<nmiss << endl;
	}
#endif		
	if (nmiss < 0) return;
	p_cpt2g[27]++;

	//____________ passing the first filter nclues band 3 see Guas3
	guash2add.GetG3(myb12f);
	gguas3 = guash2add.g3;
	gaddb3.Get3(myb12f,gguas3);// add guas3 from this add loop
#ifdef DEBUG5
	if (p_cpt2g[5] == DEBUG5) {
		cout << Char27out(gguas3) << "  gguas3" << endl;
		svb12.smin.Status("before gguas3");
		guash2add.GetG3(myb12f,1);

	}

#endif	
	gguas3 &= ~svb12.smin.mini2all;// clean redundancy


#ifdef HAVEKNOWN
	if (locdiag)
		cout << Char9out(gguas3) << "n3=" << _popcnt32(gguas3) << " nmiss=" << nmiss << endl;
#endif	
	if (gguas3) {
		register int n3 = _popcnt32(gguas3);
		if (nmiss < n3) return;
		svb12.smin.mincount += n3;
		nmin += n3;
		nmiss -= n3;
		svb12.smin.AddMincountG3(gguas3);
	}
	//chunkhadd.GetB12CX();
	// add more  from this add loop
	//gaddb3.GetMore(myb12f, chunkhadd.t12, chunkhadd.nt12);
	//svb12.AttachMore(chunkhadd.t12, chunkhadd.nt12);


	guash2add.GetB12CX(tgo, step, myb12f);
	gaddb3.GetMore(myb12f, guash2add.t12, guash2add.nt12);
	svb12.AttachMore(guash2add.t12, guash2add.nt12);
	p_cpt2g[28]++;
	//if (guash2add.nt12 >= 100) {
		//cout << Char54out(myb12f) << guash2add.nt12 << "uas [28]" << p_cpt2g[28] << endl;
		//guash2add.Debug27();
	//}
#ifdef HAVEKNOWN
	if (locdiag) {
		svb12.smin.Status(" after guas3 ");
		cout<< "[28]"<< p_cpt2g[28] << "after guas3 " << " nmin=" << nmin << " nmiss=" << nmiss
			<< " nt2=" << chunkhadd.nt12 << endl;
		//chunkhadd.Debug27();
	}
#endif	
	GoB3(nclues, svb12);
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



void GCHK::GoB3(  int ncl, VB12 & vbx) {
	p_cpt2g[8]++;
	int isdirect = 0;
	ntaddgob3 = 0;
	nclf = ncl;
#ifdef HAVEKNOWN
	int locdiag = 0;
	if ((myb12f&pk54) == myb12f) {
		cout << Char54out(myb12f) << "GoB3 p_cpt2g[8]=" << p_cpt2g[8] << endl;
		locdiag = 1;
	}
#endif
	// _direct process unless compulsory outfield clues
	if (nmiss > 2)isdirect = 1;
	else {
		vbx.bfbf2 = 0;
		if(!nmiss)	vbx.ApplyBf2();// assign critical 2 pairs
		vbx.BuildOf();// outfield after 2 pairs
#ifdef HAVEKNOWN
		if (locdiag) 
			cout << Char27out(vbx.bfbf2) << "bfb2 ntof="<<vbx.ntof << endl;
#endif	


		if (!nmiss) {
			p_cpt2g[10]++;
			if (vbx.ntof) return;
			BuildExpandB3Vect(vbx.bfbf2, vbx.smin.critbf, vbx); 
		}
		else if (nmiss == 1) {
#ifdef DEBUG5
			if (p_cpt2g[5] == DEBUG5) {
				cout <<  " miss 1 vbx.ntof=" << vbx.ntof << endl;
			}
#endif
			if (! vbx.ntof) isdirect = 1;
			else {
				p_cpt2g[11]++;
				if (vbx.smin.mini_bf2) {
					vbx.ApplyBf2();// assign critical 2 pairs
					vbx.BuildOf();// outfield after 2 pairs
				}
				uint32_t ua = vbx.GetAnd();
#ifdef HAVEKNOWN
				if (locdiag) {
					vbx.Dumptof();
					cout << Char27out(vbx.bfbf2) << " nmiss=1 bfb2 ntof=" << vbx.ntof << endl;
					cout << Char27out(ua) << "ua get and"  << endl;
				}

#endif					
				if (!ua)return;// to many clues				
				uint32_t  c;// use ua as first clue and go std
#ifdef DEBUG5
				if (p_cpt2g[5] == DEBUG5) {
					cout << Char27out(ua) << "ua get and go b3 miss 1" << endl;
				}
#endif
				while (bitscanforward(c, ua)) {
					uint32_t bit = 1 << c, bfbf2 = vbx.bfbf2 | bit;
					ua ^= bit;
					BuildExpandB3Vect(bfbf2, vbx.smin.critbf,vbx);
					if (clean_valid_done == 2)return;
				}
			}
		}
		else {// now nmiss=2
			if (vbx.ntof < 3 || vbx.GetAnd()) isdirect = 1;
			else {
				if (vbx.smin.mini_bf2) {
					vbx.ApplyBf2();// assign critical 2 pairs
					vbx.BuildOf();// outfield after 2 pairs
				}
				uint32_t ua1 = vbx.GetMinOf();// start with smallest
				uint32_t  c;// use ua as first clue and go std
				while (bitscanforward(c, ua1)) {
					uint32_t bit = 1 << c, bfbf2 = vbx.bfbf2 | bit;
					ua1 ^= bit;
					uint32_t ua2 = vbx.GetAndExcept(bfbf2);
					while (bitscanforward(c, ua2)) {
						uint32_t bit2 = 1 << c, bfbf2_2 = bfbf2 | bit2;
						ua2 ^= bit2;
						BuildExpandB3Vect(bfbf2_2, vbx.smin.critbf, vbx);
						if (clean_valid_done == 2)return;
					}
				}
			}
		}
	}
	if (isdirect) {// go direct
		p_cpt2g[13]++;
#ifdef HAVEKNOWN 
		if (locdiag ) {
			cout << Char54out(myb12) << "  go b3 is direct p_cpt2g[13]  "
				<< p_cpt2g[13] << endl;
		}
#endif
#ifdef DEBUG5
		if (p_cpt2g[5] == DEBUG5) {
			cout << " go isdirect"   << endl;
		}
#endif
		BuildExpandB3Vect(0, BIT_SET_27, vbx);
		if (b3direct.nt)ExpandB3Vect( );
	}
}

#define LIMADD 17
void  GCHK::BuildExpandB3Vect( uint32_t cl0bf, uint32_t active0,
	 VB12 & vbx) {// when assigned before
	b3direct.Init();
	uint32_t is1 = 0;
	{// try to improve the process avoiding redundancy
		//uint32_t wubuf.tbuilexpand3[5][500];
		uint32_t  ntt[5];
		memset(ntt, 0, sizeof ntt);
		register uint32_t AC = active0,F=cl0bf; 
		for (uint32_t i = 0; i < vbx.ntg2ok; i++) {
			register uint32_t U = tg2[vbx.tg2ok[i]].pat;
			if (!(U&F))  wubuf.tbuilexpand3[2][ntt[2]++] = U;
		}
		if (gguas3) {// add guas3 if any
			for (int i = 0; i < 9; i++)if (gguas3 & (1 << i))
				wubuf.tbuilexpand3[3][ntt[3]++] = 7 << (3 * i);
		}
		register uint32_t * t2 = wubuf.tbuilexpand3[2], nt2 = ntt[2];
		for (uint32_t i = 0; i < vbx.ntmore27; i++) {
			register uint32_t U = vbx.tmore27[i];
			if (!(U&F)) {
				U &= AC;
				int ccu = _popcnt32(U);
				if (!ccu) return;
				if (ccu == 1) {		is1 |= U;continue;			}
				{// clear redundancy
					register uint32_t nu = ~U, x = 0;
					for (uint32_t i = 0; i < nt2; i++)
						if (!(nu&t2[i])) { x = 1; break; }
					if (x) continue;
				}
				if (ccu > 4) ccu = 4;
				if (ntt[ccu] < 450) wubuf.tbuilexpand3[ccu][ntt[ccu]++] = U;
			}
		}
		// add small band 3 residual uas
		for (uint32_t i = 0; i < myband3.nua; i++) {
			register uint32_t U = myband3.tua[i];
			if (!(U&F)) {
				U &= AC;
				int ccu = _popcnt32(U);
				if (!ccu) return;
				if (ccu == 1) { is1 |= U; continue; }
				if (ccu > 4) ccu = 4;
				wubuf.tbuilexpand3[ccu][ntt[ccu]++] = U;
				//cout << Char27out(U) << "b3 added" << endl;
			}
		}
		for (uint32_t i = 0; i < ntaddgob3; i++) {
			register uint32_t U = taddgob3[i];
			if (!(U&F)) {
				U &= AC;
				int ccu = _popcnt32(U);
				if (!ccu) return;
				if (ccu == 1) { is1 |= U; continue; }
				if (ccu > 4) ccu = 4;
				wubuf.tbuilexpand3[ccu][ntt[ccu]++] = U;
			}
		}
#ifdef HAVEKNOWN 
		if (locdiag) {
			cout << Char27out(is1) << " is1   build ntt ";
			for (int i = 0; i < 5; i++)cout << " " << ntt[i];
			cout << " ntmore27=" << vbx.ntmore27 << endl;
			cout << Char27out(is1) << " to assign" << endl;
		}
#endif
#ifdef DEBUG14
		if (DEBUG14 == p_cpt2g[14]) {
			cout << Char27out(is1) << " is1   build ntt "<<endl;
			for (int i = 0; i < 5; i++)
				for (uint32_t j = 0; j < ntt[i]; j++)
					cout << Char27out(wubuf.tbuilexpand3[i][j]) 
					<< endl;

		}
#endif

		if (is1) {// now singles to assign
			AC &= ~is1;// minimum is to clear is1
			F |= is1;// add to assign
			MINCOUNT wmc = vbx.smin;
			int i27;
			register int X = is1;
			while (bitscanforward(i27, X)) {
				int bit27 = 1 << i27;
				X ^= bit27; //clear bit
				int i = i27 / 3, bit = 1 << i;
				if (!(bit & (wmc.mini_bf2 | wmc.mini_bf3))) {// clear mini row
					register uint32_t mask = (7 << (3 * i));
					if (!(AC & mask)) {// one clue added!!
						if (!nmiss) return;// no room for it
					}
					AC &= ~mask;
				}
			}	
		}
		cl0bf = F;// be sure to be after new assign
		active0 = AC;
		int cc = _popcnt32(F);
		if (cc > ncluesb3)return;
		if (!is1) {// apply directly stored uas
			for (int i = 2; i < 5; i++)
				for (uint32_t j = 0; j < ntt[i]; j++)  
					b3direct.Add(wubuf.tbuilexpand3[i][j]);
		}
		else {// after added clues, revise the others
			uint32_t  ntt2[5];//tt2[i][500]is wubuf.tbuilexpand3[i+5][500]
			memset(ntt2, 0, sizeof ntt2);
			for (int i = 2; i < 5; i++)
				for (uint32_t j = 0; j < ntt[i]; j++) {
					register uint32_t U = wubuf.tbuilexpand3[i][j];
					if (!(U&F)) {
						U &= AC;
						register int cc = _popcnt32(U);
						if (cc > 4) cc = 4;
						wubuf.tbuilexpand3[cc+5][ntt2[cc]++] = U;
					}
				}
			if (ntt2[0]) return; // dead branch
			for (int i = 1; i < 5; i++)
				for (uint32_t j = 0; j < ntt2[i]; j++)
					b3direct.Add(wubuf.tbuilexpand3[i+5][j]);
		}
		if (cc == ncluesb3) {
			if (b3direct.nt) return;
			NoExpandB3(F);
			return;
		}

#ifdef HAVEKNOWN 
		if (locdiag) {
			cout << " final ntt =0 from uas b3direct.nt="<< b3direct.nt << endl;
			cout << Char27out(cl0bf) << " cl0bf" << endl;
			cout << Char27out(active0) << " active0" << endl;
			b3direct.Debug();
		}
#endif
#ifdef DEBUG5
		if (p_cpt2g[5] == DEBUG5) {
			cout << " final ntt =0 from uas b3direct.nt=" << b3direct.nt << endl;
			cout << Char27out(cl0bf) << " cl0bf" << endl;
			cout << Char27out(active0) << " active0" << endl;
			b3direct.Debug(0);
		}
#endif

	}
 	if((!b3direct.nt)&& (!NoExpandB3(cl0bf))) return;
	if (!b3direct.nt) {
		cout <<"bug back noexpand3 nt=0" << " [5]=" << p_cpt2g[5] << endl;
		aigstop = 1;		return;	}
	ExpandB3Vect( cl0bf, active0);
	b3direct.nt = 0;
}
uint32_t GCHK::NoExpandB3(uint32_t cl0bf) {
	b3direct.nt = 0;
	if (!clean_valid_done) {
		clean_valid_done = 1;
		if (IsValidB12()) return 0;
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

void  GCHK::ExpandB3Vect(uint32_t cl0bf, uint32_t active0) {
	p_cpt2g[14]++;
	//if (p_cpt2g[14] > 10) return;
	uint64_t limspot = 17 - nclf - _popcnt32(cl0bf),
		limm1 = limspot - 1;
	ntaddb3 = 0;

	if (limspot<2) {
		if(limspot)ExpandB3Vect2(cl0bf, active0);
		else ExpandB3Vect1(cl0bf, active0);
		return;
	}
	if (b3direct.nt > p_cpt2g[44])p_cpt2g[44] = b3direct.nt;
	if (b3direct.nt >= 64) {
		p_cpt2g[43]++; return;
	}
#ifdef DEBUG5
	int locdiag = 0;
	if (p_cpt2g[5] == DEBUG5) {
		cout << "start p_cpt2g[14]=" << p_cpt2g[14] << " limspot=" << limspot << endl;		
		cout << " ExpandB3Vect b3direct.nt=" << b3direct.nt << endl;
		b3direct.Debug(0);
		locdiag = 1;
	}
#endif
	//cout << "start p_cpt2g[14]=" << p_cpt2g[14] << " limspot=" << limspot << endl;
	//b3direct.Debug(0);

	//uint32_t tadd[500], ntadd = 0;
	// uas basis in b3direct and chunkh.band3[2]
	struct SPB {
		uint64_t vect;
		uint32_t  possible_cells, all_previous_cells, active_cells, x;
	}spb[12];
	register SPB *s, *sn;
	register uint64_t ispot;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->all_previous_cells = cl0bf;// if former assign
	s->active_cells = active0;// default bitsets27
	s->vect = b3direct.v0;// initial vectsno ua hit here
	s->possible_cells = b3direct.tua[0] & active0;// must be if no ua2

	//____________ here start the search nclues
next:
	ispot = s - spb;
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = s->possible_cells;
		if (!p)goto back;
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		s->possible_cells ^= bit;
		s->active_cells ^= bit;
		sn = s + 1; *sn = *s;
		sn->all_previous_cells  |= bit;
		//__________________check for next uaor no more ua
		sn->vect &= b3direct.vc[cell];
	}
	if (ispot < limm1) {// here max 16 clues never valid b3

		{
			p_cpt2g[52]++;
#ifdef DEBUG5
			if (locdiag) {
				cout << Char27out(sn->all_previous_cells) << " ispot=" << ispot<<" "<< p_cpt2g[52] << endl;
				//cout << Char64out(sn->vect) << "sn" << endl;

			}

#endif
			register uint64_t V = sn->vect;
			register int ir;
			if (V) {// most often
				bitscanforward64(ir, V);//ir ua to load
				register uint32_t U = b3direct.tua[ir] & s->active_cells;
				if (!U)goto next;//dead branch
				p_cpt2g[53]++;
				sn->possible_cells = U;
				s++; // switch to next spot
				goto next;
			}
			if (ntaddb3) {// see if one fresh ua can be used
				register uint32_t Bf = sn->all_previous_cells;
				for (uint32_t ia = 0; ia < ntaddb3; ia++) {
					register uint32_t U = taddb3[ia];
					if (!(U & Bf)) {
						U &= s->active_cells;
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
				if (IsValidB12()) { clean_valid_done = 2; return; }
			}
			p_cpt2g[20]++;
			if (IsValidB3(sn->all_previous_cells)) {
				s->possible_cells &= anduab3;// same spot hitting  new ua
				sn->possible_cells = anduab3;
				s = sn;
				goto next;
			}
			cout <<"bug seen valid 16 or less [5] "
				<< p_cpt2g[5]<< " [19] " << p_cpt2g[19] <<endl;
			aigstop = 1;
			return;
		}
	}
	// this is the last step must hit all pending uas
	//cout << Char27out(sn->all_previous_cells) << "last" << endl;
	//cout << Char64out(sn->vect) << "sn" << endl;
#ifdef DEBUG5
	if (locdiag) {
		cout << Char27out(sn->all_previous_cells) << "last ispot=" << ispot << " " << p_cpt2g[52] << endl;
		//cout << Char64out(sn->vect) << "sn" << endl;

	}

#endif
	{ // find next cells hitting all uas
		int aig = 1;
		register uint32_t andw = sn->active_cells;
		register uint32_t bf = sn->all_previous_cells;
		if (ntaddb3) {// likely quicker 0
			for (uint32_t ia = 0; ia < ntaddb3; ia++) {
				register uint32_t U = taddb3[ia];
				if (!(U & bf)) {
					aig = 0;
					if (!(andw &= U)) goto next;
				}
			}
		}
		{// and still valid   uas
			register int ir;
			register uint64_t vect = sn->vect;
			while (bitscanforward64(ir, vect)) {//ir ua to load
				vect ^= (uint64_t)1 << ir;
				aig = 0;
				if(!(andw &= b3direct.tua[ir])) goto next;
			}
		}
		// no more ua or "and all uas" not empty
		if (!clean_valid_done) {
			clean_valid_done = 1;
			//cout << " checkvalid " << endl;
			if (IsValidB12()) { clean_valid_done = 2; return; }
			//cout << "valid " << endl;
		}

		//cout << Char27out(andw) << "and" << endl;

		if (aig) {// no ua could be a 17 valid
			if (IsValidB3(bf)) andw &= anduab3;
			else { Out17(bf);	goto next; }// this is a valid 17
		}
		register int cell;
		while (bitscanforward(cell, andw)) {
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
			if (IsValidB3(bf | bit))andw &= anduab3;
			else 	Out17(bf | bit);
		}
		goto next;

	}

	goto next;// never here
	// going back, for a non empty index, count it back
back:
	if (--s >= spb)goto next;
}
void  GCHK::ExpandB3Vect1(uint32_t cl0bf, uint32_t active0) {
	// only one clue free must hit all pending uas
	//cout << Char27out(cl0bf) << " cl0bf limspot=0" << endl;
	//b3direct.Debug(0);
	for (uint32_t i = 0; i < b3direct.nt; i++)// min one ua here
		active0 &= b3direct.tua[i];
	//cout << Char27out(active0) << " active0 final" << endl;
	if (!active0) return; // nothing to add
	if (!clean_valid_done) {
		clean_valid_done = 1;
		if (IsValidB12()) { clean_valid_done = 2; return; }
		//cout << "valid " << endl;
	}
	register int cell;
	while (bitscanforward(cell, active0)) {
		register uint32_t bit = 1 << cell;
		active0 ^= bit; //clear bit
		if (IsValidB3(cl0bf | bit))active0 &= anduab3;
		else 	Out17(cl0bf | bit);
	}		
}
void  GCHK::ExpandB3Vect2(uint32_t cl0bf, uint32_t active0) {
	// 2 clues free one ua then  must hit all pending uas
#ifdef DEBUG5
	int locdiag = 0;
	if (p_cpt2g[5] == DEBUG5) {
		cout << Char27out(cl0bf) << " cl0bf limspot=1 2 clues" << endl;
		b3direct.Debug(0);
		locdiag = 1;
	}
#endif		
	uint32_t cell1, first= b3direct.tua[0];// must exist checked upstream
	//cout << Char27out(first) << " first"  << endl;
	while (bitscanforward(cell1, first)) {
		register uint32_t bit = 1 << cell1;
		first ^= bit; //clear bit
		uint32_t bf1 = cl0bf | bit;
		uint64_t vect= b3direct.vc[cell1]& b3direct.v0;
#ifdef DEBUG5
		if (locdiag) {
			cout << Char27out(bf1) << "last ntaddb3=" << ntaddb3  << endl;
			cout << Char64out(b3direct.v0) << "v0" << endl;
			cout << Char64out(b3direct.vc[cell1]) << "vcx" << endl;

		}

#endif		
		//cout << Char27out(bf1) << " bf1 cell="<<cell1 << endl;
		register uint32_t andw = active0;
		if ((!vect)&& (!ntaddb3)) {// could be a 17 valid
			if (!clean_valid_done) {
				clean_valid_done = 1;
				if (IsValidB12()) { clean_valid_done = 2; return; }
				//cout << "valid no vect" << endl;
			}
			if (IsValidB3(bf1)) andw &= anduab3;
			else {Out17(bf1);	continue;	}// this is a valid 17
		}
		else {// and valid other uas
			register int ir;
			while(	bitscanforward64(ir, vect)){//ir ua to load
				vect ^= (uint64_t)1 << ir;
				andw &= b3direct.tua[ir] ;
			}
			if (ntaddb3) {// see if one fresh ua can be used
				register uint32_t Bf = cl0bf;
				for (uint32_t ia = 0; ia < ntaddb3; ia++) {
					register uint32_t U = taddb3[ia];
					if (!(U & Bf)) 	andw &= U;
				}
			}
		}
		if (!andw) continue; // nothing to add
		if (!clean_valid_done) { //if not done earleir
			clean_valid_done = 1;
			if (IsValidB12()) { clean_valid_done = 2; return; }
			//cout << "valid with vect" << endl;
		}
		register int cell;
#ifdef DEBUG5
		if (locdiag)			cout << Char27out(andw) << " andw" << endl;
#endif	
		while (bitscanforward(cell, andw)) {
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
#ifdef DEBUG5
			if (locdiag)			cout << Char27out(bf1 | bit) << " check" << endl;
#endif	
			if (IsValidB3(bf1 | bit))andw &= anduab3;
			else 	Out17(bf1 | bit);
		}
	}
}
uint32_t GCHK::IsValidB3(uint32_t bf) {
	if (zhou[0].CallCheckB3(tclues, nclf, bf)) {
		anduab3 = BIT_SET_27;
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			register uint32_t ua = w.bf.u32[2];
			taddb3[ntaddb3++] = ua;
			taddgob3[ntaddgob3++] = ua;
			anduab3 &= ua;
			if (cc0 > 20 || cc > 8)continue;
			BuildGua(w, cc);
			if (cc < 7 && cc0 < LIMADD) {
				if (cc < 3)p_cpt2g[35]++; else p_cpt2g[36]++;
				chunkh.Add128(w, cc);
				//if (cc==3 ) {
				//	cout << Char9out(w.bf.u32[2]) << " ";
				//	cout << Char54out(w.bf.u64[0]) << " " 
				//		<< " [5] " << p_cpt2g[5]<< " [35] "<<p_cpt2g[35] <<endl;
				//}
			}
			gaddb3.Add(w, cc);
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




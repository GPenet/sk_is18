
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
	if ((int)__popcnt(puzknown.bf.u32[2] ) < mincluesb3) return  0;
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
	if (1) {
		for (int i = 0; i < 81; i++) cout << zs0[i] + 1;
		cout << " reshaped  "<< band_order[0]<< band_order[1]
		<< band_order[2] << endl;
		for (int ibs = 0; ibs < 3; ibs++) {
			STD_B416 & b = bax[ibs];
			cout << b.band << "\t" << b.i416 << " " << t416n6[b.i416] << endl;
		}
		cout << " minb1=" << minb1 << " minb2=" << minb2 << endl;
	}
#ifdef HAVEKNOWN
	if (1)
		for (int i = 0; i < 81; i++) cout << ze[i + 82];
	cout << " known reshaped" << endl;
	//return 0;
#endif
	//__________________________ start uas search

	UaCollector();
	p_cpt2g[1] = tuasb12.nua;
	p_cpt2g[2] = chunkh.C2Count();
	//chunkh.DebugAll();
	BuildVectorsForExpand4B12();
	Expand4B12();
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
	BF128 tsort[27 - 4][200]; // max is 27  min is 4
	uint32_t ntsort[27 - 4];
	memset(ntsort, 0, sizeof ntsort); // empty table
	{ //sort all by size
		BF128 * t = tuas81.tall;
		uint32_t nt = tuas81.ntall;
		for (uint32_t iua = 0; iua < nt; iua++) {
			BF128 wt = t[iua];
			uint32_t cc = wt.bf.u16[7] - 4;// index is 0 for 4
			wt.bf.u32[3] = 0; // clear digits and count
			tsort[cc][ntsort[cc]++] = wt;
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
			tsort[cc][ntsort[cc]++] = wt;
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
			tsort[cc][ntsort[cc]++] = wt;
		}
	}
	// reload tall and check subsets/redundancy
	tuas81.ntall = 0;
	for (int i = 0; i < 24; i++)if (ntsort[i]) {
		BF128 * tt = tsort[i];
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
			if (U2x) 	chunkh.Addc2(U2x, i27);
		}
	}

	// check final status for c2
	cout << " last ic2 " << chunkh.ic2 << " last index "
		<< chunkh.c2[chunkh.ic2].nt
		<< " " << 64 * chunkh.ic2 + chunkh.c2[chunkh.ic2].nt << endl;
}
inline void GCHK::NewGua(BF128 w) {
	register uint64_t ua12 = w.bf.u64[0];
	uint64_t ua54 = (ua12 & BIT_SET_27) |
		((ua12 & BIT_SET_B2) >> 5);
	int cc = _popcnt32(w.bf.u32[2]);
	uint64_t cc0 = _popcnt64(w.bf.u64[0]);
	w.bf.u64[0] = ua54;// store in 54 cells mode
	if (cc < 3) 	w.bf.u32[2] = GetI27(w.bf.u32[2]);
	chunkh.Add128(w, cc);
}
void  GCHK::PutUasStartInVector() {	// split guas 2/3 and others
	chunkh.Init();
	for (uint32_t iua = 0; iua < tuas81.ntold; iua++) {
		BF128 wt = tuas81.told[iua];
		if (wt.bf.u64[0])		NewGua(wt);
		else chunkh.Addband3(wt.bf.u32[2]);
	}
}

//________________________________ start expand 3/4 + no more uas
#define nphase1 4

GCHK::CPT_4::CPT_4() {
	t[0] = units3xBM[18].u64[0];	t[1] = units3xBM[19].u64[0];
	t[2] = units3xBM[20].u64[0];	t[3] = units3xBM[21].u64[0];
	t[4] = units3xBM[22].u64[0];	t[5] = units3xBM[23].u64[0];// 6 box
	t[6] = band3xBM[0].u64[0];		t[7] = band3xBM[1].u64[0];// 2 bands
	t[8] = band3xBM[3].u64[0];		t[9] = band3xBM[4].u64[0];
	t[10] = band3xBM[5].u64[0];		 //3 stack

}
void GCHK::BuildVectorsForExpand4B12() {
	morev2a.Init();
	morev2c.Init();	morev2d.Init();
	// switch here tua to 54 mode
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
		limspot = nphase1-1;// expand 3/4 clues
	nt4_to_expand = 0;
	// _______________ expand to have all  3/4 cells to expand more
	struct SPB {// spots to find starts to extract uas 3 digits
		uint64_t  possible_cells, all_previous_cells, active_cells,
			v;
	}spb[5], *s, *sn;
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
		cout << "no more uas" << endl;		goto next;
	}
	uint32_t ir;
	bitscanforward64(ir, sn->v);//ir ua to load
	uint64_t Ru = twu[ir] & ac;
	if (!Ru)goto next;//dead branch unlikely
	sn->possible_cells = Ru;
	s++; // switch to next spot
	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s >= spb)goto next;
	//cout << "count expand4b12 n="  << nt4_to_expand << endl;
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

//	if(nt4ok<0) 	return;
//	else nt4_to_expand = nt4ok+1;
#endif
	cout << "nt4_to_expand=" << nt4_to_expand << endl;
	for (uint64_t i = 0; i < nt4_to_expand; i++) {
#ifdef HAVEKNOWN
		if(nt4ok>=0)		okcheck = (i == nt4ok);
		else okcheck = 0;
		if (okcheck) cout <<i  << " start expand the expected ckunk" << endl;
#endif
		Do_phase2(t4_to_expand[i]);
		//if(p_cpt2g[0] == 3) cout << " back expand aigstop="<<aigstop << " for i=" << i << endl;
		if (aigstop) 			return;
	
	}

}
//______  process a "no more uas bands 1+2


void GCHK::Do_phase2(T4_TO_EXPAND w) {
	p_cpt2g[3]++;
	{	register uint64_t bf54 = w.bf, ac54 = w.active;

	{	//____ load the first clues
		register uint64_t R = bf54;
		nclues_step = 0;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			tclues[nclues_step++] = cell;
		}
		tcluesxy = &tclues[nclues_step];
	}
	morev2a.Init(); // now in vector
	{// build still valid uas sorted by size 
		uint64_t tt[30][500], ntt[30];
		memset(ntt, 0, sizeof ntt);
		{
			for (uint32_t i = 0; i < tuasb12.nua; i++) {
				register uint64_t R = tuasb12.tua[i];
				if (R&bf54) continue;
				R &= ac54;
				if (!R)return; //dead
				register uint64_t cc = _popcnt64(R);
				tt[cc][ntt[cc]++] = R;
			}
		}
		{
			ntua4 = 0;
			for (int i = 0; i < 30; i++)if (ntt[i]) {
				for (uint64_t j = 0; j < ntt[i]; j++) {
					register uint64_t R = tt[i][j], Rn = ~R;
					for (uint32_t k = 0; k < ntua4; k++)
						if (!(tua4[k] & Rn)) goto nextj;
					tua4[ntua4++] = R;
				nextj:;
				}
			}
		}
	}
	//______________________________________________end pick up
	{
		register int wnua = ntua4;
		memset(v12_v0, 0, sizeof v12_v0);
		for (int i = 0; i < 64; i++) {
			if (wnua >= 64) {
				v12_v0[i] = ~0; wnua -= 64;
			}
			else {//0 63 last
				v12_v0[i] = ((uint64_t)1 << wnua) - 1;// 0 if wnua=0
				break;
			}
		}
	}
	for (int i = 0; i < 54; i++)memcpy(v12_c[i], v12_v0, sizeof v12_c[0]);
	for (uint32_t i = 0; i < ntua4; i++) {
		register uint64_t R54 = tua4[i];
		register uint32_t cell, bloc = i >> 6, ir = i - 64 * bloc;
		while (bitscanforward64(cell, R54)) {
			uint64_t bit = (uint64_t)1 << ir;
			R54 ^= (uint64_t)1 << cell; //clear bit
			v12_c[cell][bloc] ^= bit;
		}
	}
	uint32_t n = 0;
	for (uint32_t i = 0; i < 64; i++) {
		register uint64_t A = v12_v0[i];
		if (!A)break;
		V12_64 &w = tv12_64[0][i];
		w.ind = i;
		w.u_start = i << 6;// 64*i
		w.v = A;
		n++;
	}
	nv12_64_spot[0] = n;// nb chunks at start
	tv12_64f[0] = &tv12_64[0][n];
	Do_phase2Expand(bf54, ac54);
	}
}
void  GCHK::Do_phase2Expand(uint64_t bf, uint64_t ace) {// mode 54 not 2x
#ifdef HAVEKNOWN
	if (okcheck) {
		okcheck = 1;
		if ((bf&pk54) == bf) {
			cout << Char54out(bf) << " possible expand p_cpt2g[3]=" << p_cpt2g[3] << endl;
			uint64_t w = (ace &pk54) | bf;
			cout << Char54out(w) << " ace & sol";
			if (w == pk54) { cout << " will be ok"; okcheck = 2; }
			cout << endl;
		}
	}

#endif
	pendbufvalid = &bufvalid[BUFVALIDS];
	pbufvalid = bufvalid;
	uint64_t limspot = 17 - nphase1 - mincluesb3,
		limclean=1;// set to more if low {min b1+b2}
	//if(minb1<5 && minb2<5)limclean = 2;
	if (minb1b2 < 9)limclean = 2;
	//if (minb1 < 4 && minb2 < 4)limclean = 3;
	if (minb1b2 <7)limclean = 3;
	uint64_t andvalid = ~0,orvalid=0;
	// _______________ expand to have all  minimal < n clues
	struct SPB {// spots to find starts to extract uas 3 digits
		uint64_t  possible_cells, all_previous_cells, active_cells;
	}spb[13];
	register SPB  *s, *sn;
	register uint64_t ispot;
	register int cell;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->all_previous_cells = bf;
	s->active_cells = ace;
	s->possible_cells = tua4[0];
	//____________ here start the search nclues
next:
	ispot = s - spb;
	{// catch and apply cell in bitfields
		register uint64_t p = s->possible_cells;
		if (!p)goto back;
		bitscanforward64(cell, p);
		register uint64_t bit = (uint64_t)1 << cell;
		s->possible_cells ^= bit;
		s->active_cells ^= bit;
		sn = s + 1; *sn = *s;
		sn->all_previous_cells = s->all_previous_cells | bit;
	}
	{// find next ua
		register V12_64 *tv1 = tv12_64[ispot], *tv2 = tv12_64[ispot+1],
			*tv1f = tv12_64f[ispot], *tv2w = tv2;
		register uint64_t *vc = v12_c[cell];
		while (tv1 < tv1f) {// 64 bits uas vectors active
			register uint64_t V = tv1->v & vc[tv1->ind];
			if (V) { tv2w->v = V; tv2w++->ind = tv1->ind; }
			tv1++;
		}
		if (tv2 != tv2w) {// still at least one ua not hit
			tv12_64f[ispot + 1] = tv2w;// new last active ua vector
			register uint32_t ir;
			bitscanforward64(ir, tv2->v);//relative index first active
			register uint64_t Ru = tua4[ir + tv2->u_start] & s->active_cells;
			if (!Ru)goto next;
			if (ispot >= limspot)goto next;
			// here dynamic limit depending on minb1/mminb2
			if (ispot == limclean) {// clean 2+3/4 clues or more
				if (pbufvalid > bufvalid) {
					CleanBufferAndGo(andvalid,orvalid);
					andvalid = ~0; orvalid = 0;
				}
			}
			sn->possible_cells = Ru;
			s++; // switch to next spot
			goto next;
		}
	}
	{// possible valid band 1+2 to store
		if (pbufvalid >= pendbufvalid) {// if limite do
			CleanBufferAndGo(andvalid,orvalid);
			andvalid = ~ 0; orvalid = 0;
		}
		andvalid &= sn->all_previous_cells;
		orvalid |= sn->all_previous_cells;
		*pbufvalid++ = sn->all_previous_cells;
		if (ispot != limspot)*pbufvalid++ = s->active_cells;
	}
	goto next;
back:
	if (--s >= spb)goto next;
	if (pbufvalid > bufvalid) // if limite do
		CleanBufferAndGo(andvalid,orvalid);
}
void GCHK::CleanBufferAndGo(uint64_t andvalid, uint64_t orvalid) {
	if (aigstop) return;
	p_cpt2g[4]++;
	uint64_t *myend = pbufvalid, *pw = bufvalid;
	limb12 = 18 - mincluesb3;
	pbufvalid = bufvalid;
	uint64_t nn8 = (myend - bufvalid),
		mincb2=_popcnt64(andvalid>>27);
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
	{// catch more still valid after andvalid
		moreand.ntm = 0;
		morev2a.Extract(tclues, nclues_step, moreand.tm, moreand.ntm);
		morev2c.Extract(tclues, nclues_step, moreand.tm, moreand.ntm);
		morev2d.Extract(tclues, nclues_step, moreand.tm, moreand.ntm);
	}
	chunkh.Get2Clean(chvb12.t2, chvb12.nt2, tclues, nclues_step);
	chunkh.GetMoreClean(chvb12.tmore, chvb12.ntmore, tclues, nclues_step);
	//chvb12.SortTmore();
	gaddb3.Init();
#ifdef HAVEKNOWN
	if (okcheck==2) {
		if ((andvalid&pk54) == andvalid) {
			cout << Char54out(andvalid) << " p_cpt2g[4]=" << p_cpt2g[4] << endl;
		}
	}

#endif
	while (pw < myend) {
		p_cpt2g[5]++;
		myb12 = *pw++;
		int ncl = (int)_popcnt64(myb12);
		if (ncl < limb12)myac = *pw++;		else myac = 0;
		if ((mincluesb3 + _popcnt64(myb12 >> 32)) > 12) continue;
#ifdef HAVEKNOWN
		int locdiag = 0;
		if (okcheck == 2) {
			if ((myb12&pk54) == myb12) {
				cout << Char54out(myb12) << " p_cpt2g[5]=" << p_cpt2g[5] 
					<<" ncl=" << ncl<<" limb12="<<limb12<< endl;
				locdiag = 1;
				//moreand.Dump();
			}
		}			
#endif
		p_cpt2g[16]++;
		// don't apply moreand and if free cells available in bands 1+2
		if (ncl == limb12) 
			if (moreand.Check(myb12))continue;// check "more" not common
#ifdef HAVEKNOWN
		if (locdiag)		cout << "ok moreand" << endl;
#endif
		{//____ load the remaining  clues
			register uint64_t R = myb12 ^ andvalid;
			nclues = nclues_step;
			uint32_t cell;
			while (bitscanforward64(cell, R)) {
				R ^= (uint64_t)1 << cell; //clear bit
				tclues[nclues++] = cell;
			}
		}
#ifdef COLOIN
		if (IsValidB12()) continue;
		{
			if (a_18_seen++ > 5) continue;;
			char zs[82]; zs[81] = 0;		memcpy(zs, ze, 81);
			for (uint64_t i = 0, bit = 1; i < 54; i++, bit <<= 1)
				if (!(myb12&bit))zs[i] = '.';
			cout << zs << endl;			
			continue;
		}
#endif
		clean_valid_done = 0;
		if (nclues == limb12) 	InitGoB3(myb12, myac);	
		else {// apply moreand mis uas if anysing
			MOREANDB mab;
			mab.GetdUnHit(myb12, moreand.tm, moreand.ntm);
			if (mab.ntm) {
#ifdef HAVEKNOWN
				if (locdiag) {
					cout << "have unhit ua'(s)" << endl;
					//mab.Dump();
				}
#endif
				CleanMoreUas(myb12, myac, nclues, mab);
				continue;
			}
			else CheckValidBelow(myb12, myac);
		}
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
			moreand.tm[moreand.ntm++] = ua54;
			mbisvalid.tm[mbisvalid.ntm++] = ua54;
			if (cc < 20) {
				p_cpt2g[31]++;
				morev2a.Add54(ua54);  //12-20
				if(tuasb12.nua<4096)tuasb12.tua[tuasb12.nua++] = ua54;
				//cout << Char54out(ua54) << " add ";
				//cout << Char54out(myb12) << " nclues=" << nclues << endl;
			}
			else if (cc < 21)morev2c.Add54(ua54); // 19 20
			else if (cc < 23) morev2d.Add54(ua54);// 21-22
		}
	}
	return (int)mbisvalid.ntm;
}
void GCHK::InitGoB3(uint64_t bf54, uint64_t ac54) {// bands1+2 locked
	p_cpt2g[6]++;
#ifdef HAVEKNOWN
	int locdiag = 0;
	if (okcheck == 2) {
		//cout << Char54out(bf54) << endl;
		if ((bf54&pk54) == bf54) {
			cout << Char54out(bf54) << " p_cpt2g[6]=" << p_cpt2g[6] << endl;
			locdiag = 1;
		}
	}
#endif
	clean_valid_done = 0;
	myb12f = bf54;
	int nmin = svb12.GetsminF(chvb12.t2, chvb12.nt2);
	ncluesb3 = 18 - limb12;
	nmiss = ncluesb3 - nmin;
#ifdef HAVEKNOWN
	if (locdiag) 
		cout << Char54out(bf54) << " nmin=" << nmin  << "nmiss=" << nmiss<< endl;
#endif
	if (nmiss >= 0) {
		svb12.CleanTmore(chvb12.tmore, chvb12.ntmore);
		GoB3(limb12,  svb12);
	}
}
void GCHK::InitGoB3below(uint64_t bf54, uint64_t ac54) {// bands1+2 locked
	p_cpt2g[6]++;
#ifdef HAVEKNOWN
	int locdiag = 0;
	if (okcheck == 2) {
		//cout << Char54out(bf54) << endl;
		if ((bf54&pk54) == bf54) {
			cout << Char54out(bf54) << " below  p_cpt2g[6]=" << p_cpt2g[6] << endl;
			locdiag = 1;
		}
		//else return;
	}
#endif
	myb12f = bf54;
	int nmin = svb12b.GetsminF(chvb12b.t, chvb12b.nt2);
	ncluesb3 = 18 - nclues;
	nmiss = ncluesb3 - nmin;
#ifdef HAVEKNOWN
	if (locdiag)
		cout << Char54out(bf54) << " nmin=" << nmin << "nmiss=" << nmiss << endl;
#endif
	if (1) return;
	if (nmiss >= 0) {
		BF128 * tm = &chvb12b.t[chvb12b.nt2];
		uint32_t	ntm = chvb12b.nt- chvb12b.nt2;
		svb12.CleanTmore(tm, ntm);
		GoB3(nclues, svb12b);
	}
}
void GCHK::CleanMoreUas(uint64_t bf, uint64_t ac, int ncl, MOREANDB & mabo) {
	p_cpt2g[18]++;
#ifdef HAVEKNOWN
	int locdiag=0;
	if (okcheck == 2) {
		if ((bf&pk54) == bf) {
			cout << Char54out(bf) << " p_cpt2g[18]=" << p_cpt2g[18]
				<< " miss=" << (limb12 - ncl) << endl;
			locdiag = 1;
		}
	}
#endif
	uint64_t ua = mabo.tm[0];
	MOREANDB mabn;
	if (ncl == limb12 - 1) // last step must hit all uas
		for (uint64_t i = 1; i < mabo.ntm; i++)ua &= mabo.tm[i];
	ua &= ac;// only active cells used
#ifdef HAVEKNOWN
	if (locdiag)cout << Char54out(ua) << "ua to add" << endl;
#endif
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
#ifdef HAVEKNOWN
		//if (locdiag) 
		//	cout << Char54out(bf | bit) << " mabn.ntm="<< mabn.ntm << " nclues="<<nclues<< endl;
#endif

		if (nclues == limb12) {
			if (mabn.ntm)continue;
			else {
				myb12f = bf | bit;
				myacf = ac;
				mynclues = nclues;
#ifdef HAVEKNOWN
				if (locdiag) {
					locdiag = 1;
					if (myb12f==pk54) {
						cout << Char54out(myb12f) << " expected band 12 after more uas"  << endl;
						locdiag = 2;
					}					
				}
#endif
				ncluesb3 = 18 - mynclues;
				chvb12b.Shrink(chvb12, myb12f, myac);
				int nmin = svb12add.GetsminF(chvb12b.t, chvb12b.nt2);
				nmiss = ncluesb3 - nmin;
				if (nmin > ncluesb3) continue;
				svb12add.CleanTmore(&chvb12b.t[chvb12b.nt2], 
					chvb12b.nt- chvb12b.nt2);

#ifdef HAVEKNOWN
				if (locdiag==2) {
					//chvb12b.Dumpt2();
					//chvb12b.Dumptmore();
					svb12add.smin.Status(" ");
					cout << "nmin=" << nmin << endl;
					svb12add.Dumptmore27();
				}
#endif
				GoB3(mynclues , svb12add);
			}
		}
		else {
			continue;
			if (mabn.ntm)CleanMoreUas(bf | bit, ac, nclues, mabn);
			else  CheckValidBelow(bf | bit, ac);
		}
	}
}
void GCHK::CheckValidBelow(uint64_t bf, uint64_t ac) {
	p_cpt2g[7]++;
#ifdef HAVEKNOWN
	int locdiag = 0;
	if (okcheck == 2) {
		if ((bf&pk54) == bf) {
			cout << Char54out(bf) << " start below p_cpt2g[7]=" << p_cpt2g[7] << endl;
			locdiag = 1;
		}
	}

#endif	
	if (IsValidB12()) {// add one cell hitting all uas
#ifdef HAVEKNOWN
		if (locdiag)  
			cout  << " not valid go to clean more uas"  << endl;
#endif	
		if (nclues == limb12) return;
		CleanMoreUas( bf,  ac, nclues, mbisvalid);
		return;
	}
	clean_valid_done = 1;
#ifdef HAVEKNOWN
	if (locdiag) cout << "  valid enter the add process" << endl;
#endif	
	p_cpt2g[17]++;
	chvb12b.Shrink(chvb12, bf, ac);
	InitGoB3below(bf, ac);// try first direct with more clurs in b3
#ifdef HAVEKNOWN
	if (locdiag) cout << "  back from direct" << endl;
#endif	
	{//______________  build tadd from active cells
		ntadd = 0;
		memset(vaddh.mapcell, 255, sizeof vaddh.mapcell);
		int cell;
		register uint64_t V = ac;// still valid c2
		while (bitscanforward64(cell, V)) {
			V ^= (uint64_t)1 << cell;
			vaddh.mapcell[cell] = ntadd;
			tadd[ntadd++] = cell;
		}
	}
	vaddh.SetUpAdd0(chvb12b, ac);
	ExpandAddB1B2(bf);
}
void GCHK::VADD_HANDLER::Dump() {
	cout << Char64out((v0&v2).bf.u64[0]) << " c2" << endl;
	BF128 w=v0&vm;
	cout << Char64out(w.bf.u64[0]) << " cmore1" << endl;
	cout << Char64out(w.bf.u64[1]) << " cmore2" << endl;

}
void GCHK::VADD_HANDLER::SetUpAdd0(CHUNKVB12ADD& vbo, uint64_t ac) {// build vector for still valid in vbo
	memset(vcells, 255, sizeof vcells);
	uint64_t nt = vbo.nt;
	v0 = maskLSB[nt];
	for (uint32_t i = 0; i < nt; i++) {
		uint64_t w = vbo.t[i].bf.u64[0] & ac;
		uint32_t cell;
		while (bitscanforward64(cell, w)) {
			w ^= (uint64_t)1 << cell;
			vcells[mapcell[cell]].Clear (i);
		}
	}
	v2= maskLSB[vbo.nt2];
	vm = v0 - v2;
}

// en attente de révision

void GCHK::ExpandAddB1B2Go(int step) {
	p_cpt2g[26]++;
	myb12f = myb12add;
	BF128 va = vaddh.vsteps[step];// call the process for this
	ncluesb3 = 18 - mynclues - step;
	chvb12b.nt2a = chvb12b.ntma = 0;
	uint32_t i64;
	{ // get active c2
		register uint64_t V = (va&vaddh.v2).bf.u64[0];// still valid c2
		while (bitscanforward64(i64, V)) {
			V ^= (uint64_t)1 << i64;
			chvb12b.t2a[chvb12b.nt2a++] = chvb12b.t[i64];
		}
	}
	int nmin = svb12add.GetsminF(chvb12b.t2a, chvb12b.nt2a);
	nmiss = ncluesb3 - nmin;
	if (nmin > ncluesb3) return;
	p_cpt2g[27]++;
	BF128 wm = va & vaddh.vm;
	register uint64_t V = wm.bf.u64[0];
	while (bitscanforward64(i64, V)) {
		V ^= (uint64_t)1 << i64;
		chvb12b.tma[chvb12b.ntma++] = chvb12b.t[i64 ];
	}
	V = wm.bf.u64[1];
	while (bitscanforward64(i64, V)) {
		V ^= (uint64_t)1 << i64;
		chvb12b.tma[chvb12b.ntma++] = chvb12b.t[i64+64];
	}
	svb12add.CleanTmore(chvb12b.tma, chvb12b.ntma);
	GoB3(mynclues + step, svb12add);
}
void GCHK::ExpandAddB1B2(uint64_t bf) {// add up to n cells
	mynclues = (uint32_t)_popcnt64(bf);
	uint32_t *tgo=&tclues[mynclues],	lim=17- mincluesb3- mynclues;
	vaddh.vsteps[0] = vaddh.v0;
	struct SPB {// spots to find starts to extract uas 3 digits
		uint64_t  all_previous_cells, icur;
	}spb[20], * s, * sn;
	s = spb;
	memset(s, 0, sizeof spb[0]);//	s->icur=0
	s->all_previous_cells = bf;// mode 54 previous cells
	//____________ here start the search nclues
next:
	if (s->icur >= ntadd) goto back;
	uint64_t ispot = s - spb;
	register int cell = tadd[s->icur++];
	register uint64_t bit = (uint64_t)1 << cell;
	tgo[ispot] = cell;
	sn = s + 1; *sn = *s; 
	vaddh.Apply(ispot, s->icur-1);
	sn->all_previous_cells = s->all_previous_cells | bit;
	myb12add = sn->all_previous_cells;
	ExpandAddB1B2Go((int)ispot + 1);// call the process for this 
	if (ispot >= lim) 	goto next;
	s++; // switch to next spot (icur+1)
	goto next;
back:
	if (--s >= spb)goto next;

}



//============ clues in band 3 (no more clue in bands 1+2)

int GCHK::VB12::GetsminF(	BF128 * t2,	uint32_t nt2) {// apply filter myb12
	myt2 = t2;
	mynt2 = nt2;
	register uint64_t F = gchk.myb12f;
	memset(&smin, 0, sizeof smin);
	ntg2ok = 0;
	int v27 = 0;
	for (uint32_t i = 0; i < nt2; i++) {
		if (F & t2[i].bf.u64[0]) continue;
		register int i27 = t2[i].bf.u32[2],
			bit = 1 << i27;
		if (!(v27&bit)) {
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
	}
// same with guas2 added in this chunk
	for (uint32_t i = 0; i < gchk.gaddb3.nt2; i++) {
		if (F & gchk.gaddb3.t2[i].bf.u64[0]) continue;
		register int i27 = gchk.gaddb3.t2[i].bf.u32[2],
			bit = 1 << i27;
		if (!(v27&bit)) {
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
	}
	smin.SetMincount();
	return smin.mincount;

}
void GCHK::VB12::CleanTmore(BF128 * tmore, uint32_t ntmore) {// clean subsets redundancy and sort
	//mytmore = tmore;
	//myntmore = ntmore;
	uint32_t tw[400],tt[6][100], ntt[6];
	int ntw = 0;
	memset(ntt, 0, sizeof ntt);
	for (uint32_t i = 0; i < ntg2ok; i++)
		tw[ntw++] = gchk.tg2[tg2ok[i]].pat;
	register uint64_t F = gchk.myb12f;

	for (uint32_t i = 0; i < gchk.gaddb3.ntmore; i++) {
		if (F & gchk.gaddb3.tmore[i].bf.u64[0]) continue;
		register uint32_t R = gchk.gaddb3.tmore[i].bf.u32[2],
			cc = _popcnt32(R)-3;// 3 is the minimum in tmore
		if(cc>5 )cc = 5;// should not be more then 8 clues in b3
		tt[cc][ntt[cc]++] = R;
	}

// same with guasmore added in the chunk
	for (uint32_t i = 0; i < ntmore; i++) {
		if (F & tmore[i].bf.u64[0]) continue;
		register uint32_t R = tmore[i].bf.u32[2],
			cc = _popcnt32(R) - 3;// 3 is the minimum in tmore
		if (cc > 5)cc = 5;// should not be more then 8 clues in b3
		tt[cc][ntt[cc]++] = R;
	}


	ntmore27 = 0;// store still valid
	for (uint32_t i = 0; i < 6; i++) {
		for (uint32_t j = 0; j < ntt[i]; j++) {
			register uint32_t R = tt[i][j],
				Rn = ~R;
			for (int k = 0; k < ntw; k++) 
				if (!(tw[k] & Rn)) goto nextj;
			tw[ntw++] = R;
			tmore27[ntmore27++] = R;
		
		nextj:;
		}
	}
	 
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
inline void GCHK::VB12::BuildOrOf(CHUNKVB12& o, uint64_t ac) {
	ntof = 0; orof = 0;
	register uint32_t F = smin.critbf; // assigned or critbf
	for (uint32_t i = 0; i < o.ntmore; i++) {
		BF128 w = o.tmore[i];
		w.bf.u64[0] &= ac;
		register uint32_t U = w.bf.u32[2];
		if (!(U&F)) {
			tof128[ntof++] = w;
			orof |= w.bf.u64[0];
		}
		if (ntof > 20) return;//should never be
	}
	// same with band 3
	BF128 w; w.SetAll_0();
	for (uint32_t i = 0; i < bax[2].nua; i++) {
		register uint32_t U = bax[2].tua[i] & BIT_SET_27;
		w.bf.u32[2] = U;
		if (!(U&F))tof128[ntof++] = w;
		if (ntof > 20) return;//should never be
	}
}



void GCHK::GoB3(  int ncl, VB12 & vbx) {
	p_cpt2g[8]++;
	int isdirect = 0;
	ntaddgob3 = 0;
	nclf = ncl;
#ifdef HAVEKNOWN
	int locdiag = 0;
	if (okcheck == 2) {
		//cout << Char54out(bf54) << endl;
		if ((myb12f&pk54) == myb12f) {
			cout << Char54out(myb12f) << " p_cpt2g[8]=" << p_cpt2g[8] << endl;
			vbx.smin.Status("");
			vbx.Dumptmore27();
			locdiag = 1;
		}
	}
#endif
	// _direct process unless compulsory outfield clues
	if (nmiss > 2)isdirect = 1;
	else {
		vbx.bfbf2 = 0;
		if(!nmiss)	vbx.ApplyBf2();// assign critical 2 pairs
		vbx.BuildOf();// outfield after 2 pairs
#ifdef HAVEKNOWN
		if (locdiag) {
			cout << Char27out(vbx.bfbf2) << "bfb2" << endl;
			vbx.Dumptof();
		}

#endif		
		if (!nmiss) {
			p_cpt2g[10]++;
			if (vbx.ntof) return;
			BuildExpandB3Vect(vbx.bfbf2, vbx.smin.critbf, vbx);
		}
		else if (nmiss == 1) {
			if (! vbx.ntof) isdirect = 1;
			else {
				p_cpt2g[11]++;
				vbx.ApplyBf2();// assign critical 2 pairs
				uint32_t ua = vbx.GetAnd();
				if (!ua)return;// to many clues				
				uint32_t  c;// use ua as first clue and go std
				while (bitscanforward(c, ua)) {
					uint32_t bit = 1 << c, bfbf2 = vbx.bfbf2 | bit;
					ua ^= bit;
					BuildExpandB3Vect(bfbf2, vbx.smin.critbf,vbx);
				}
			}
		}
		else {// now nmiss=2
			if (vbx.ntof < 2 || vbx.GetAnd()) isdirect = 1;
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
				}
			}
		}
	}
	if (isdirect) {// go direct
		p_cpt2g[13]++;
		b3direct.Init();
		for (uint32_t i = 0; i < vbx.ntg2ok; i++)
			b3direct.Add(tg2[vbx.tg2ok[i]].pat);
		for (uint32_t i = 0; i < vbx.ntmore27; i++)
			b3direct.Add(vbx.tmore27[i]);
		ExpandB3Vect( );
	}
}
void  GCHK::BuildExpandB3Vect( uint32_t cl0bf, uint32_t active0,
	 VB12 & vbx) {// when assigned before
	int cc= _popcnt32(cl0bf);
	if (cc > ncluesb3)return;
	if (cc == ncluesb3) {// nothing to expand
		p_cpt2g[19]++;
		if (!clean_valid_done) {
			clean_valid_done = 1;
			//cout << " check b12"	<< endl;
			if (IsValidB12()) return;
			//cout << " check b12 valid" << endl;
		}
		if (zhou[0].CallCheckB3(tclues, nclf, cl0bf)) {
			BF128 w = zhgxn.tua[0];// take the first
			int cc = _popcnt32(w.bf.u32[2]);
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			if (cc <= 6 && cc0 < 17)NewGua(w);
			return;
		}
		//_________________________   valid 18 
		cout << " valid 18 to process BuildExpandB3Vect" << endl;

		return;
	}
	b3direct.Init();
	for (uint32_t i = 0; i < vbx.ntg2ok; i++) {
		register uint32_t U = tg2[vbx.tg2ok[i]].pat;
		if(!(U&cl0bf))	b3direct.Add(U);
	}
	for (uint32_t i = 0; i <vbx.ntmore27; i++) {
		register uint32_t U = vbx.tmore27[i];
		if (!(U&cl0bf))	b3direct.Add(U&active0);
	}
	for (uint32_t i = 0; i < ntaddgob3; i++) {
		b3direct.Add(taddgob3[i]);
	}
	ExpandB3Vect( cl0bf, active0);

}
void  GCHK::ExpandB3Vect(uint32_t cl0bf, uint32_t active0) {
	p_cpt2g[14]++;
	int limspot = 17 - nclf - _popcnt32(cl0bf);
	uint32_t tadd[500], ntadd = 0;
	// uas basis in b3direct and chunkh.band3[2]
	struct SPB {
		uint64_t vects[3];
		uint32_t  possible_cells, all_previous_cells, active_cells, x;
	}spb[12];
	register SPB *s, *sn;
	register uint64_t ispot;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->all_previous_cells = cl0bf;// if former assign
	s->active_cells = active0;// default bitsets27
	s->vects[0] = b3direct.v0;// initial vectsno ua hit here
	s->vects[1] = chunkh.band3[0].v0;// initial vects
	s->vects[2] = chunkh.band3[1].v0;// initial vects
	if (cl0bf) {// apply cl0bf to band3
		register uint32_t w = cl0bf, c;
		while (bitscanforward(c, w)) {
			w ^= 1 << c;
			s->vects[1] &= chunkh.band3[0].vc[c];
			s->vects[2] &= chunkh.band3[1].vc[c];
		}
	}
	s->possible_cells = b3direct.tua[0];
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
		sn->all_previous_cells = s->all_previous_cells | bit;
		//__________________check for next uaor no more ua
		sn->vects[0] &= b3direct.vc[cell];
		sn->vects[1] &= chunkh.band3[0].vc[cell];
		sn->vects[2] &= chunkh.band3[1].vc[cell];
	}
	{ // find next ua
		register int ir = -1;
		{
			register uint32_t U;
			if (sn->vects[0]) {// most often
				bitscanforward64(ir, sn->vects[0]);//ir ua to load
				U = b3direct.tua[ir];	}
			else if (ntadd) {// 
				for (uint32_t ia = 0; ia < ntadd; ia++) {
					if (!(tadd[ia] & sn->all_previous_cells)) {
						ir = 1; U = tadd[ia]; break;		}
				}
			}
			if (ir < 0) {
				if (sn->vects[1]) {// allways if less than 64 uas for the band
					bitscanforward64(ir, sn->vects[1]);//ir ua to load
					U = chunkh.band3[0].tua[ir];		}
				else if (sn->vects[2]) {
					bitscanforward64(ir, sn->vects[2]);//ir ua to load
					U = chunkh.band3[1].tua[ir];	}
			}
			if (ir >= 0) {// still at least one ua
				if (ispot >= limspot) goto next;// to many clues
				U &= s->active_cells;
				if (!U)goto next;//dead branch  
				sn->possible_cells = U;
				s++; // switch to next spot
				goto next;
			}
		}
	}
	//___________  possible nclues do final check
	p_cpt2g[19]++;
	if (!clean_valid_done) {
		clean_valid_done = 1;
		if (IsValidB12()) return;
	}
	p_cpt2g[20]++;

	if (zhou[0].CallCheckB3(tclues, nclf, sn->all_previous_cells)) {
		register uint32_t  anduan = BIT_SET_27;
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			register uint32_t ua = w.bf.u32[2];
			tadd[ntadd++] =ua;
			taddgob3[ntaddgob3++] = ua;

			anduan &= ua;
			if (cc0 > 20|| cc>8)continue;
			register uint64_t ua12 = w.bf.u64[0];
			uint64_t ua54 = (ua12 & BIT_SET_27) |
				((ua12 & BIT_SET_B2) >> 5);
			w.bf.u64[0] = ua54;// store in 54 cells mode
			if (cc < 3) 	w.bf.u32[2] = GetI27(ua);
			gaddb3.Add(w, cc);
			if (cc < 7 && cc0 < 17) {
				if(cc<3)p_cpt2g[35]++;else p_cpt2g[36]++;
				chunkh.Add128(w, cc);
			}
		}
		s->possible_cells &= anduan;// same spot hitting  new ua		
		if (ispot < limspot) {// new ua for next spot
			sn->possible_cells = anduan;
			s = sn;
		}
		goto next;
	}
	//_________________________   valid 18 
	Out17(sn->all_previous_cells);
	/*		return; // stop at first	*/
	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s >= spb)goto next;
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
	cout << ws << " one sol "  << endl;
#ifdef TEST_ON
#endif
	//if(zp)strcpy(zp, ws);
	//a_18_seen = 1;
	//aigstop = 1; 

}



//#define DEBUGA 4842934
//#define DEBUG3 80
//#define DEBUGB 107551
//#define DEBUGC 111462
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
	if ((int)__popcnt(puzknown.bf.u32[2] ) < mincluesb3) return  0;
#endif		
	int * zs0 = grid0, *zs0_diag = grid0_diag;
	BANDMINLEX::PERM perm_ret;
	memcpy(zs0, bax[0].band0, sizeof  bax[0].band0);
	memcpy(&zs0[27], bax[1].band0, sizeof  bax[0].band0);
	memcpy(&zs0[54], bax[2].band0, sizeof  bax[0].band0);


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
	BuildVectorsForExpand4B12();
	Expand4B12();
	return a_18_seen;

}


//________________ uas generation and store in chunkh
void GCHK::UaCollector() {

	Sg2Setup();
	//for (int i = 0; i < 27; i++)tg2[i].Dump();

	FirstUasCollect();// uas 2 digits bands and stacks
	SecondUasCollect();// uas 345 bands 1+2
	UasCollect4box();// 4 box in bands 1+2
	PutUasStartInVector();
	Guas2Collect();
}
void GCHK::Adduab12(int print) {
	for (uint32_t i = 0; i < zh2gxn.nua; i++) {
		uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);

		if (AddUaB12UN(w, cc))
			if (print)cout << Char2Xout(w) << " added to table" << endl;;
	}
}
void GCHK::FirstUasCollect() {
	//____________________________ produce all uas 2/3 digits
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
	//tuas81.Debug2();
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
	//cout << "end FirstUasCollect() tuas81.ntold=" << tuas81.ntold
	//	<< " tuasb12.nua=" << tuasb12.nua << endl;


}
void GCHK::SecondUasCollect() {
	//cout << "final all 2 digits n=" << tuasb12.nua << endl;
	//__________ start collect 345 digits in bands 1+2
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
	//cout << "final all 3 digits n=" << tuasb12.nua << endl;
	for (int i = 0; i < 126; i++) {
		int ir = zh2_4[0].GoZ4(floors_4d[i]);
		if (ir < 8) continue;// minimum for a fresh ua 4 digits
		uint64_t F = zh2gxn.unsolved_field;
		tuasb12.GetT2(F);
		//cout << Char9out(floors_4d[i]) << "study ua  i=" << i << endl << endl;
		zh2_4[0].DoZ4Go();
		if (zh2gxn.nua) Adduab12();
	}
	//cout << "final  all 4 digits n=" << tuasb12.nua << endl;

	for (int i = 0; i < 126; i++) {
		int ir = zh2_5[0].GoZ5(0777 ^ floors_4d[i]);
		uint64_t F = zh2gxn.unsolved_field;
		tuasb12.GetT2(F);
		//if (i == 0) {
			//cout << Char9out(0777 ^ floors_4d[i]) << "study ua  i=" << i 
			//	<< " nt2=" << tuasb12.nt2 << " nua=" << tuasb12.nua << endl;
			//tuasb12.DumpT2();
			//zh2_5[0].ImageCandidats();
		//}
		//cout << Char9out(0777 ^ floors_4d[i]) << "study ua  i=" << i 
		//	<< " nt2=" << tuasb12.nt2 << " nua=" << tuasb12.nua << endl;
		//tuasb12.DumpT2();
		//zh2_5[0].ImageCandidats();
		zh2_5[0].DoZ5Go();
		if (zh2gxn.nua) Adduab12();
	}
	//cout << "final  all 5 digits n=" << tuasb12.nua << endl;

}
void GCHK::UasCollect4box() {
	//___________________ box 1245
	//cout << "stack 3 empty" << endl;
	BF64 wf4;
	wf4.bf.u32[0] = wf4.bf.u32[1] = 077077077;// stack 3 bands 12
	tuasb12.GetT2(wf4.bf.u64);
	//tuasb12.DumpT2();
	zh2b[0].InitB1245(grid0);
	//zh2b[0].ImageCandidats();
	zh2b[0].Do4bGo();
	if (zh2gxn.nua) Adduab12();
	//cout << "nt2 end=" << tuasb12.nt2 << "tuasb12.nua=" << tuasb12.nua << endl;


	//___________________ box 1346
	//cout << "stack 2 empty" << endl;
	wf4.bf.u32[0] = wf4.bf.u32[1] = 0707707707;
	tuasb12.GetT2(wf4.bf.u64);
	//tuasb12.DumpT2();
	zh2b[0].InitB1346(grid0);
	//zh2b[0].ImageCandidats();
	zh2b[0].Do4bGo();
	if (zh2gxn.nua) Adduab12();
	//cout << "nt2 end=" << tuasb12.nt2 << "tuasb12.nua=" << tuasb12.nua << endl;

	//___________________ box 2356
	//cout << "stack 1 empty" << endl;
	wf4.bf.u32[0] = wf4.bf.u32[1] = 0770770770;
	tuasb12.GetT2(wf4.bf.u64);
	//tuasb12.DumpT2();
	zh2b[0].InitB2356(grid0);
	//zh2b[0].ImageCandidats();
	zh2b[0].Do4bGo();
	if (zh2gxn.nua) Adduab12();
	cout << "nt2 end=" << tuasb12.nt2 << "tuasb12.nua=" << tuasb12.nua << endl;

}
void GCHK::Guas2Collect() {
	//cout << "entry guas2 collect" << endl;
	//chunkh.C2Status();

	// guas 2 digits and stack g2s already there 
	//for (int i = 0; i < 27; i++)  tg2[i].Dump();
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
			//cout << Char2Xout(ua2x) << " i27=" << i27 << endl;
			tg2uas[i27][ntg2uas[i27]++] = ua2x;
		}
	}
	//return;
	for (int i27 = 0; i27 < 27; i27++) {// 27 g2
		SG2 w = tg2[i27];
		uint64_t uasi27[200], nuasi27 = ntg2uas[i27];
		memcpy(uasi27, tg2uas[i27], sizeof tg2uas[0]);
		if (nuasi27 == 1 && _popcnt64(uasi27[0]) == 2) continue;// nothing to do
		// __________________find guas 3 digits

		//cout << "\t\t\ti27=" << i27 << endl;
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
				//cout << Char2Xout(U) << "added " << _popcnt64(U) << endl;
				continue;
			}
			//cout << Char9out(fl) << " digit processed i27=" << i27 << endl;
			uint64_t F = zh2gxn.unsolved_field;
			tuasb12.GetT2(F, uasi27, nuasi27);
			//tuasb12.DumpT2();
			int ir2 = zh2_3[0].DoZ3Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					uasi27[nuasi27++] = U;
					//cout << Char2Xout(U) << "added "<< _popcnt64(U) << endl;
				}
		}
		//_________________find uas 4 digits
		//cout << "\t\t\ti27=" << i27 <<"  go4"<< endl;

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
			//cout << Char9out(fl) << "\t\t digit processed i27=" << i27 << endl;
			uint64_t F = zh2gxn.unsolved_field;
			tuasb12.GetT2(F, uasi27, nuasi27);
			//tuasb12.DumpT2();
			int ir2 = zh2_4[0].DoZ4Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					uasi27[nuasi27++] = U;
					//cout << Char2Xout(U) << "added " <<  _popcnt64(U) << endl;
				}
		}

		
		//_________________find uas 5 digits
		//cout << "\t\t\ti27=" << i27 << "  go5" << endl;
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
			//cout << Char9out(fl) << "\t\t digit processed i27=" << i27 << endl;
			uint64_t F = zh2gxn.unsolved_field;
			tuasb12.GetT2(F, uasi27, nuasi27);
			//tuasb12.DumpT2();
			int ir2 = zh2_5[0].DoZ5Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i], cc = _popcnt64(U);
					if (cc > 20) continue;
					uasi27[nuasi27++] = U;
					//cout << Char2Xout(U) << "added " << cc << endl;
				}
		}
		
		// store results in chunkh.c2
		uint64_t *tk = tg2uas[i27], ntk = ntg2uas[i27];
		for (uint64_t i = 0; i < nuasi27; i++) {
			register uint64_t U2x = uasi27[i];
			for (uint64_t j = 0; i < ntk; j++) {// redundancy
				if (U2x == tk[j]) { U2x = 0; break; }
			}
			if (U2x) 	chunkh.Addc2(U2x, i27);
		}
	}

	//chunkh.C2Status();// check final status for c2
	cout << " last ic2 " << chunkh.ic2 << " last index "
		<< chunkh.c2[chunkh.ic2].nt
		<< " " << 64 * chunkh.ic2 + chunkh.c2[chunkh.ic2].nt << endl;
}
void CHUNKS_HANDLER::Add128(BF128 w) {
	register uint64_t ua12 = w.bf.u64[0];
	uint32_t ua = w.bf.u32[2],
		cc = _popcnt32(ua);
	uint64_t ua54 = (ua12 & BIT_SET_27) |
		((ua12 & BIT_SET_B2) >> 5);
	if (cc == 2) {// add to c2 switch to index0-26
		ua = gchk.GetI27(ua);
		if (ic2 == 39 && c2[39].nt > 63) return;
		if (c2[ic2].nt > 63)c2[++ic2].Init();
		c2[ic2].Add(ua54, ua);
	}

	else {// add to cmore
		if (icmore == 19 && cmore[19].nt > 63) return;
		if (cmore[icmore].nt > 63)
			cmore[++icmore].Init();
		cmore[icmore].Add(ua54, ua);

	}
}
void  GCHK::PutUasStartInVector() {
	chunkh.Init();
	// split guas 2/3 and others
	for (uint32_t iua = 0; iua < tuas81.ntold; iua++) {
		BF128 wt = tuas81.told[iua];
		if (wt.bf.u64[0])		chunkh.Add128(wt);
		else chunkh.Addband3(wt.bf.u32[2]);
	}
	//chunkh.band3[0].Debug();
	//chunkh.band3[1].Debug();
}

//__ start expand 3/4 + no more uas
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
	morev2a.Init();	morev2b.Init();
	morev2c.Init();	morev2d.Init();
	for (int i = 0; i < 54; i++)memset(v12_4_c, 255, sizeof v12_4_c);
	for (uint32_t i = 0; i < 64; i++) {
		register uint64_t R = tuasb12.tua[i] & BIT_SET_2X;
		uint32_t xcell;
		while (bitscanforward64(xcell, R)) {
			uint32_t cell = From_128_To_81[xcell];
			uint64_t bit = (uint64_t)1 << i;
			R ^= (uint64_t)1 << xcell; //clear bit
			v12_4_c[cell] ^= bit;
		}
	}
}
void  GCHK::Expand4B12() {

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
	s->active_cells = BIT_SET_2X;
	s->possible_cells = twu[0] & BIT_SET_2X;
	s->v = ~0;// initial nothing done
	//____________ here start the search nclues
next:
	uint64_t ispot = s - spb;
	// catch and apply cell in bitfields
	register int xcell;
	uint64_t p = s->possible_cells;
	if (!p)goto back;
	bitscanforward64(xcell, p);
	register uint64_t bit = (uint64_t)1 << xcell;
	s->possible_cells ^= bit;
	int cell = From_128_To_81[xcell];
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
		register uint64_t K = puzknown.bf.u64[0], Kn = ~K;
		for (uint64_t i = 0; i < nt4_to_expand; i++) {
			T4_TO_EXPAND w=t4_to_expand[i];
			register uint64_t B = w.bf, A = w.active;
			if (B & Kn)continue;// not partial known
			if (K & (~(B | A))) continue;
			if (nt4ok < 0) nt4ok = (int)i;
			cout << Char2Xout(B) << " bf" << endl;
			cout <<Char2Xout(A) << " ac" <<endl;
		}
	}
	cout << "expected i=" << nt4ok << endl;
	if(nt4ok<0) 	return;
	else nt4_to_expand = nt4ok+1;
#endif
	cout << "nt4_to_expand=" << nt4_to_expand << endl;
	for (uint64_t i = 0; i < nt4_to_expand; i++) {
#ifdef HAVEKNOWN
		okcheck = (i == nt4ok);
#endif
		Do_phase2(t4_to_expand[i]);
	}

}
//______  process a "no more uas bands 1+2
inline void GCHK::VB12::GetOld(VB12 & vbo, uint64_t bf54) {
	// t2 and tmore reduced by fresh cells in bf54
	cout << Char64out(bf54) << " bf54 en 64" << endl;
	//vbo.Dumpt2();
	//cout << "end vbo t2" << endl;
	nt2 = ntmore = 0;
	for (uint32_t i = 0; i < vbo.nt2; i++) {
		BF128 w = vbo.t2[i];
		if (!(w.bf.u64[0] & bf54))	t2[nt2++] = w;
	}
	for (int i = 0; i < morevalidc2.nt; i++) {
		BF128 w = morevalidc2.ta[i];
		if (!(w.bf.u64[0] & bf54))	t2[nt2++] = w;
	}
	for (uint32_t i = 0; i < vbo.ntmore; i++) {
		BF128 w = vbo.tmore[i];
		if (!(w.bf.u64[0] & bf54))tmore[ntmore++] = w;
	}
	for (int i = 0; i < morevalidothers.nt; i++) {
		BF128 w = morevalidothers.ta[i];
		if (!(w.bf.u64[0] & bf54))tmore[ntmore++] = w;
	}
}


void GCHK::Do_phase2(T4_TO_EXPAND w) {
	p_cpt2g[3]++;
	{	register uint64_t bf54 = (w.bf & BIT_SET_27) | ((w.bf & BIT_SET_B2) >> 5),
		ac54 = (w.active & BIT_SET_27) | ((w.active & BIT_SET_B2) >> 5);

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
	{// build still valid uas sorted by size 
		uint64_t tt[30][500], ntt[30];
		memset(ntt, 0, sizeof ntt);
		{
			register uint64_t F = w.bf, A = w.active;

			for (uint32_t i = 0; i < tuasb12.nua; i++) {
				register uint64_t R = tuasb12.tua[i];
				if (R&F) continue;
				R &= A;
				if (!R)return; //dead
				R = (R & BIT_SET_27) | ((R& BIT_SET_B2) >> 5);// now r54
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
	if (0) {
		cout << Char2Xout(w.bf) << " p_cpt2g[3] " << p_cpt2g[3] << " ntua4=" << ntua4
			<< " ntuas=" << tuasb12.nua
			<< " nc2=" << chunkh.GetC2Count()
			<< " p_cpt2g[4] " << p_cpt2g[4] << endl;
	}
	//______________________________________________end pick up
#ifdef DEBUGA
	if (p_cpt2g[3] > DEBUG3)return;
	if (p_cpt2g[3] > DEBUG3) {
		cout << Char2Xout(w.bf)
			<< " p_cpt2g[3] " << p_cpt2g[3]
			<< " ntua4=" << ntua4
			<< " ntuas=" << tuasb12.nua
			<< " nc2=" << chunkh.GetC2Count()
			<< " p_cpt2g[4] " << p_cpt2g[4] << endl;
	}
#endif
	//return;
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
	pendbufvalid = &bufvalid[BUFVALIDS];
	pbufvalid = bufvalid;
	uint64_t limspot = 17 - nphase1 - mincluesb3;
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
			if (ispot == 1) {// clean 2+3/4 clues
				if (pbufvalid > bufvalid) {
					CleanBufferAndGo(andvalid,orvalid);
					andvalid = ~0; orvalid = 0;
				}
			}
			sn->possible_cells = Ru;
			s++; // switch to next spot
			//if (ispot < 2)cout << ispot << endl;
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
	p_cpt2g[8]++;
	uint64_t *myend = pbufvalid, *pw = bufvalid;
	limb12 = 18 - mincluesb3;
	pbufvalid = bufvalid;
	uint64_t nn8 = (myend - bufvalid),
		mincb2=_popcnt64(andvalid>>27);

	//____ load the common clues
	{
		register uint64_t R = andvalid;
		nclues_step = 0;
		//nclues_step = (int)_popcnt64(w.bf);
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
		morev2b.Extract(tclues, nclues_step, moreand.tm, moreand.ntm);
		morev2c.Extract(tclues, nclues_step, moreand.tm, moreand.ntm);
		morev2d.Extract(tclues, nclues_step, moreand.tm, moreand.ntm);
	}
#ifdef DEBUGA
	if (p_cpt2g[3] > DEBUG3) {
		cout << "CleanBufferAndGo p_cpt2g[8]= " << p_cpt2g[8] << " nn8=" << nn8 << endl;
		cout << Char54out(andvalid) << " andvalid" << endl;
		cout << Char54out(orvalid) << " orvalid" << endl;
		cout << "extract ntm=" << moreand.ntm << " ncstep=" << nclues_step << endl;
		moreand.Dump();
		//return;
	}
#endif
	while (pw < myend) {
		p_cpt2g[4]++;
		myb12 = *pw++;
		int ncl = (int)_popcnt64(myb12);
		if (ncl < limb12)myac = *pw++;
		else myac = 0;
		if ((mincluesb3 + _popcnt64(myb12 >> 32)) > 12) continue;
#ifdef HAVEKNOWN
		//if (okcheck)cout << Char2Xout(bf) << " p_cpt2g[4]=" << p_cpt2g[4] << endl;
		if (okcheck&& p_cpt2g[4] == 4842934) {
			cout << Char2Xout(puzknown.bf.u64[0]) << " known" << endl;
			cout << Char2Xout(bf) << " p_cpt2g[4]=" << p_cpt2g[4] << endl;
			cout << Char2Xout(ac) << " ac=" << endl;
	}
#endif

		p_cpt2g[16]++;
		if (moreand.Check(myb12))continue;// check remaining "add" 
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
		morevalidc2.Init(); morevalidothers.Init();
		if (nclues == limb12) {// 	InitGoB3(bf54, ac54, nclues);
			continue;		}
		CheckValidBelow();
		/*	12 936 966		4  entry clean	*/
	}
}
void GCHK::CheckValidBelow() {
	register uint64_t uar = IsValidB12();// loop on add if not valid
	clean_valid_done = 1;
	if (uar) {// add one cell hitting all uas
		int ncluemore = nclues++;
		uar &= myac;// no redundancy
		uint32_t cell;
		while (bitscanforward64(cell, uar)) {
			clean_valid_done = 0;
			uar ^= (uint64_t)1 << cell; //clear bit
			tclues[ncluemore++] = cell;
			if (nclues == limb12) {	// InitGoB3(bf54, ac54, nclues);
				continue;
			}
			CheckValidBelow();// check more missing uas
		}
	}
	// continue in add mode 
}
uint64_t GCHK::IsValidB12() {
	p_cpt2g[5]++;
	if (zh2b[1].IsValid(tclues, nclues)) {
		register uint64_t uaand = ~0; 
			for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i],
				cc = _popcnt64(ua),
				ua54= (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			uaand &= ua54;
			moreand.tm[moreand.ntm++] = ua54;
			if (cc < 19) {
				morev2a.Add54(ua54);  //12-20
				//cout << Char54out(ua54) << " add ";
				//cout << Char54out(myb12) << " nclues=" << nclues << endl;
			}
			else if (cc < 21)morev2c.Add54(ua54); // 19 20
			else if (cc < 23) morev2d.Add54(ua54);// 21-22
		}
		return uaand;
	}
	return 0;
}
void GCHK::AfterExpandB12(uint64_t bf54, uint64_t ac54) {// no more ua
	//return;
	//if (p_cpt2g[42] > 500) return;

	p_cpt2g[5]++;
	InitGoB3(bf54, ac54, nclues);// valid puzzle  band3 bands1+2 locked
	if (nclues >= 18 - mincluesb3) return;// 2 is min clues in b3
	p_cpt2g[7]++;

	//_______________ room for add ib bands 1+2 
	if (p_cpt2g[22] < svb12.nt2)p_cpt2g[22] = svb12.nt2;
	if (p_cpt2g[23] < svb12.ntmore)p_cpt2g[23] = svb12.ntmore;
	myb12 = bf54;
	myac = ac54;
	mynclues = nclues;
	//______________  build tadd from ac54
	{
		ntadd = 0;
		memset(vaddh.mapcell, 255, sizeof vaddh.mapcell);
		int cell;
		register uint64_t V = ac54;// still valid c2
		while (bitscanforward64(cell, V)) {
			V ^= (uint64_t)1 << cell;
			vaddh.mapcell[cell] = ntadd;
			tadd[ntadd++] = cell;
		}
	}
	vaddh.SetUpAdd0(svb12, ac54);
	ExpandAddB1B2();


}


void GCHK::VADD_HANDLER::SetUpAdd0(VB12& vbo, uint64_t ac) {
	memset(vaddcell, 255, sizeof vaddcell);
	uint64_t nt2 = vbo.nt2; if (nt2 > 64) nt2 = 64;// safety should never be
	vadd0.vc2 = maskLSB[nt2].u64[0];
	for (uint64_t i = 0, bit = 1; i < nt2; i++, bit <<= 1) {
		uint64_t w = vbo.t2[i].bf.u64[0]&ac;
		//cout << Char54out(w) << " c2 "<< vbo.t2[i].bf.u32[2] << endl;
		uint32_t cell,mcell;
		while (bitscanforward64(cell, w)) {
			w ^= (uint64_t)1 << cell;
			mcell= mapcell[cell];
			vaddcell[mcell].vc2 ^= bit;
		}
	}
	uint32_t ntmore = vbo.ntmore; if (ntmore > 128) ntmore = 128;
	uint32_t imore = ntmore >> 6;
	if (imore) {
		vadd0.vcmore[0] = ~0;
		vadd0.vcmore[1] = maskLSB[ntmore - 64 ].u64[0];
	}
	else {
		vadd0.vcmore[1] = 0;
		vadd0.vcmore[0] = maskLSB[ntmore ].u64[0];
	}
	for (uint32_t i = 0; i < ntmore; i++) {
		uint64_t w = vbo.tmore[i].bf.u64[0]&ac;
		//cout << Char54out(w) << "cmore " << i<< endl;
		uint32_t cell, mcell,
			bloc=i>>6,ir=i-64*bloc;
		uint64_t bit = (uint64_t)1 << ir;
		while (bitscanforward64(cell, w)) {
			w ^= (uint64_t)1 << cell;
			mcell = mapcell[cell];
			vaddcell[mcell].vcmore[bloc] ^= bit;
		}
	}	
}

void GCHK::ExpandAddB1B2Go(int step) {
	p_cpt2g[26]++;

#ifdef DEBUG
	//if (p_cpt2g[26] > 1725)return;
	//if (p_cpt2g[26] > 1747)return; //c6 ok sol
	if (p_cpt2g[26] > 1748)return; //c7 pas ok sol
	if (p_cpt2g[26] >=1748)cout << " p_cpt2g[26] " << p_cpt2g[26] << endl;

#endif
	//if (p_cpt2g[26] < 50) {
		//cout << Char54out(myb12add) << " stepo= " << step << " " << _popcnt64(myb12add) << endl;
	//}
	VADD& va = vaddh.vaddsteps[step];// call the process for this 
	ncluesb3 = 18 - mynclues - step;
	svb12add.nt2 = svb12add.ntmore = 0;
	uint32_t i64;
	{ // get active c2
		register uint64_t V = va.vc2;// still valid c2
		while (bitscanforward64(i64, V)) {
			V ^= (uint64_t)1 << i64;
			svb12add.t2[svb12add.nt2++] = svb12.t2[i64];
		}
	}
	// insert more c2
	{
		BF128* w128 = morevalidc2.ta;
		int n128 = morevalidc2.nt;
		register uint64_t F = myb12add;
		for (int i = 0; i < n128; i++) {
			BF128 w = w128[i];
			register uint64_t V = w.bf.u64[0];
			if(!(V&F))svb12add.t2[svb12add.nt2++] = w;
		}
	}

	int nmin = svb12add.Getsmin();
	nmiss = ncluesb3 - nmin;
#ifdef DEBUG
	if (p_cpt2g[26] == DEBUG)
	{
		cout << Char54out(myb12add) << " stepo= " << step << " " << _popcnt64(myb12add)
			<< "  p_cpt2g[26]=" << p_cpt2g[26] << endl;
		cout << " nmin=" << nmin << "ncluesb3=" << ncluesb3<< endl;
	}


#endif

	if (nmin > ncluesb3) return;
	p_cpt2g[27]++;
	for (uint32_t i = 0; i <2; i++) {
		register uint64_t V = va.vcmore[i],i0=64*i;
		while (bitscanforward64(i64, V)) {
			V ^= (uint64_t)1 << i64;
			svb12add.tmore[svb12add.ntmore++] = svb12.tmore[i64+i0];
		}
	}
	// insert more others
	{
		BF128* w128 = morevalidothers.ta;
		int n128 = morevalidothers.nt;
		register uint64_t F = myb12add;
		for (int i = 0; i < n128; i++) {
			BF128 w = w128[i];
			register uint64_t V = w.bf.u64[0];
			if (!(V & F))svb12add.tmore[svb12add.ntmore++] = w;
		}
	}
#ifdef DEBUG
	if (p_cpt2g[26] == DEBUG) 	{
		cout 	<< "ntmore=" << svb12add.ntmore	<< " nmin="<<nmin  << endl;
	}
#endif
	GoB3(mynclues + step, svb12add);
}


void GCHK::ExpandAddB1B2() {
	uint32_t *tgo=&tclues[mynclues],	lim=17- mincluesb3- mynclues;
	vaddh.vaddsteps[0] = vaddh.vadd0;
	struct SPB {// spots to find starts to extract uas 3 digits
		uint64_t  all_previous_cells, icur;
	}spb[20], * s, * sn;
	s = spb;
	memset(s, 0, sizeof spb[0]);//	s->icur=0
	s->all_previous_cells = myb12;// mode 54
#ifdef DEBUG
	//if (p_cpt2g[7] >= 131)
	//	cout << "bbbb ntadd=" << ntadd << " lim=" << lim<< endl;
	
#endif

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
#ifdef DEBUG
	//if (p_cpt2g[26] == DEBUG)
	//	cout << "bbbb back" << endl;
	
#endif
}



//============ clues in band 3 (no more clue in bands 1+2)

int GCHK::VB12::Getsmin() {
	memset(&smin, 0, sizeof smin);
	ntg2ok = 0;
	int v27 = 0;
	for (uint32_t i = 0; i < nt2; i++) {
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
	smin.SetMincount();
	return smin.mincount;

}
void GCHK::VB12::CleanTmore() {// clean subsets redundancy and sort
	uint32_t tw[400], ntw = 0,
		tt[6][100],ntt[20];
	memset(ntt, 0, sizeof ntt);
	for (uint32_t i = 0; i < ntg2ok; i++)
		tw[ntw++] = gchk.tg2[tg2ok[i]].pat;
	for (uint32_t i = 0; i < ntmore; i++) {
		register uint32_t R = tmore[i].bf.u32[2],
			cc = _popcnt32(R)-3;// 3 is the minimum in tmore
		if(cc>5 )continue;// should not be more then 8 clues in b3
		tt[cc][ntt[cc]++] = R;
	}
	ntmore27 = 0;// store still valid
	for (uint32_t i = 0; i < 6; i++) {
		for (uint32_t j = 0; j < ntt[i]; j++) {
			register uint32_t R = tt[i][j],
				Rn = ~R;
			for (uint32_t k = 0; k < ntw; k++) 
				if (!(tw[k] & Rn)) goto nextj;
			tw[ntw++] = R;
			tmore27[ntmore27++] = R;
		
		nextj:;
		}
	}
//	for (uint32_t i = 0; i < ntmore27; i++)  
//		cout << Char27out(tmore27[i]) << " " << i << endl;
	 
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
	for (uint32_t i = 0; i < ntmore; i++) {
		register uint32_t U = tmore[i].bf.u32[2];
		if (!(U&F))tof[ntof++] = U;
		if (ntof > 20) return;//should never be
	}
	// same with band 3
	for (uint32_t i = 0; i < bax[2].nua; i++) {
		register uint32_t U = bax[2].tua[i] & BIT_SET_27;
		if (!(U&F))tof[ntof++] = U;
		if (ntof > 20) return;//should never be
	}
}
inline void GCHK::VB12::BuildOrOf(uint64_t ac) {
	ntof = 0; orof = 0;
	register uint32_t F = smin.critbf; // assigned or critbf
	for (uint32_t i = 0; i < ntmore; i++) {
		BF128 w = tmore[i];
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


void GCHK::InitGoB3(uint64_t bf54, uint64_t ac54, int ncl) {// bands1+2 locked
#ifdef DEBUGA
	if (p_cpt2g[4] == DEBUGA) {
		cout << Char54out(bf54) << " bf54 ncl=" << ncl << endl;
		cout << Char54out(ac54) << " ac54 "	 << endl;
	}
#endif	
	chunkh.Get2(svb12.t2, svb12.nt2, bf54, ac54, tclues, ncl);
	chunkh.GetMore(svb12.tmore, svb12.ntmore, bf54, ac54, tclues, ncl);
	int nmin = svb12.Getsmin();
	ncluesb3 = 18 - ncl;
	nmiss = ncluesb3 - nmin;
	if (nmiss >= 0)GoB3(ncl, svb12);
}



void GCHK::GoB3(  int ncl,  VB12 & vbx) {
#ifdef DEBUGA
	if (p_cpt2g[4] == DEBUGA) {
		cout << "entry gob3 nmiss=" << nmiss << " ncl=" << ncl << endl;
		cout << vbx.nt2 << " " << vbx.ntmore << endl;
	}
#endif
	int isdirect = 0;
	vbx.CleanTmore();
	nclf = ncl;

	// _direct process unless compulsory outfield clues
	if (nmiss > 2)isdirect = 1;
	else {
		vbx.ApplyBf2();// assign critical 2 pairs
		vbx.BuildOf();// outfield after 2 pairs
		if (!nmiss) {
			p_cpt2g[10]++;
			if (vbx.ntof) return;
			BuildExpandB3Vect(vbx.bfbf2, vbx.smin.critbf,vbx);
		}
		else if (nmiss == 1) {
			if (! vbx.ntof) isdirect = 1;
			else {
				p_cpt2g[11]++;
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
#ifdef DEBUGA
			if (p_cpt2g[4] == DEBUGA) {
				cout << " nmiss2 vbx.ntof=" << vbx.ntof << endl;
				vbx.smin.Status(" ");
				vbx.Dumpt2();
				vbx.Dumptmore();
				vbx.Dumptof128();
			}
#endif
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
					BuildExpandB3Vect(bfbf2_2, vbx.smin.critbf,vbx);
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
		return;
	}
}
void  GCHK::BuildExpandB3Vect( uint32_t cl0bf, uint32_t active0, VB12 & vbx) {// when assigned before
	int cc= _popcnt32(cl0bf);
	if (cc > ncluesb3)return;
	if (cc == ncluesb3) {// nothing to expand
		p_cpt2g[19]++;
		if (zhou[0].CallCheckB3(tclues, nclf, cl0bf)) {
			BF128 w = zhgxn.tua[0];// take the first
			int cc = _popcnt32(w.bf.u32[2]);
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			if (cc <= 8 && cc0 < 20)chunkh.Add128(w);
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
	for (uint32_t i = 0; i < vbx.ntmore; i++) {
		register uint32_t U = vbx.tmore[i].bf.u32[2];
		if (!(U&cl0bf))	b3direct.Add(U&active0);
	}
	ExpandB3Vect( cl0bf, active0);

}
void  GCHK::ExpandB3Vect( uint32_t cl0bf ,uint32_t active0 ) {
	p_cpt2g[14]++;
#ifdef DEBUGA
	if (p_cpt2g[4] == DEBUGA) {
		cout << "call expandb3 p_cpt2g[14]" << p_cpt2g[14]
			<< " cl0bf=" << cl0bf << " active0=" << active0
			<< "\tb3direct.nt=" << b3direct.nt << " nclf=" << nclf<< endl;
		b3direct.Debug(0);
		//return;
	}
#endif
	int limspot = 17 - nclf-_popcnt32(cl0bf);


	uint32_t tadd[500], ntadd = 0;
	// uas basis in b3direct and chunkh.band3[2]
	struct SPB {
		uint64_t vects[3];
		uint32_t  possible_cells, all_previous_cells, active_cells, x;
	}spb[12], *s, *sn;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->all_previous_cells = cl0bf;// if former assign
	s->active_cells = active0;// default bitsets27
	s->vects[0] = b3direct.v0;// initial vectsno ua hit here
	s->vects[1] = chunkh.band3[0].v0;// initial vects
	s->vects[2] = chunkh.band3[1].v0;// initial vects
	if (cl0bf) {// apply cl0bf to band3
		register uint32_t w = cl0bf,c;
		while (bitscanforward(c, w)) {
			w ^= 1 << c;
			s->vects[1] &= chunkh.band3[0].vc[c];
			s->vects[2] &= chunkh.band3[1].vc[c];
		}
	}
	s->possible_cells = b3direct.tua[0];
	//____________ here start the search nclues
next:
	uint64_t ispot = s - spb;
	// catch and apply cell in bitfields
	register int cell;
	uint32_t p = s->possible_cells;
	if (!p)goto back;
	bitscanforward(cell, p);
	register uint32_t bit = 1 << cell;
	s->possible_cells ^= bit;
	s->active_cells ^= bit;
	uint32_t ac = s->active_cells;
	sn = s + 1; *sn = *s;
	sn->all_previous_cells = s->all_previous_cells | bit;
	//__________________check for next uaor no more ua
	sn->vects[0] &= b3direct.vc[cell];
	sn->vects[1] &= chunkh.band3[0].vc[cell];
	sn->vects[2] &= chunkh.band3[1].vc[cell];
	int ir = -1;
	{
		register uint32_t U;
		if (sn->vects[0]) {// most often
			bitscanforward64(ir, sn->vects[0]);//ir ua to load
			U = b3direct.tua[ir];
		}
		else if (ntadd) {// 
			for (uint32_t ia = 0; ia < ntadd; ia++) {
				if (!(tadd[ia] & sn->all_previous_cells)) {
					ir = 1; U = tadd[ia]; break;
				}
			}
		}
		if (ir < 0) {
			if (sn->vects[1]) {// allways if less than 64 uas for the band
				bitscanforward64(ir, sn->vects[1]);//ir ua to load
				U = chunkh.band3[0].tua[ir];
			}
			//else  {
				//if(bitscanforward64(ir, sn->vects[2]))//ir ua to load
				//U = chunkh.band3[1].tua[ir];
			//}
			else if (sn->vects[2]) {
				bitscanforward64(ir, sn->vects[2]);//ir ua to load
				U = chunkh.band3[1].tua[ir];
			}
		}
		if (ir>=0) {// still at least one ua
			if (ispot >= limspot) goto next;// to many clues
			U &= ac;
			if (!U)goto next;//dead branch  
			sn->possible_cells = U;
			s++; // switch to next spot
			goto next;
		}
	}
	//___________  possible nclues do final check
	p_cpt2g[19]++;
#ifdef DEBUGA
	if (p_cpt2g[4] == DEBUGA) {
		cout << Char27out(sn->all_previous_cells) << " call check"   << endl;
	}
#endif	
	if (zhou[0].CallCheckB3(tclues, nclf, sn->all_previous_cells)) {
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			if(ntadd<100)tadd[ntadd++] = w.bf.u32[2];
			if (cc <= 8 && cc0 < 20) {
				chunkh.Add128(w);
				register uint64_t ua12 = w.bf.u64[0];
				uint64_t ua54 = (ua12 & BIT_SET_27) |
					((ua12 & BIT_SET_B2) >> 5);
				w.bf.u64[0]=ua54;// store in 54 cells mode
				if (cc < 3) {
					w.bf.u32[2] = GetI27(w.bf.u32[2]); 
					morevalidc2.Add(w);
				}
				else morevalidothers.Add(w);
			}
		}
		uint32_t ua = zhgxn.tua[0].bf.u32[2];
#ifdef DEBUGA
		if (p_cpt2g[4] == DEBUGA) {
			cout << Char27out(ua) << " new ua ispot/limspot "<<ispot <<" "<<limspot << endl;
			if (!ua) {
				cout << "bug expandb3 ua null p_cpt2g[4]=" << p_cpt2g[4] << endl;
				zhou[0].CallCheckB3(tclues, nclf, sn->all_previous_cells, 2);
				zhgxn.tua[0].Print3(" back ua");
			}
		}
#endif	
		if (ispot < limspot) {// new ua for next spot
			sn->possible_cells = ua;
			s = sn;
		}
		else {// continue same spot hitting the new ua
			s->possible_cells &= ua;
		}
		goto next;
	}
	//_________________________   valid 18 
	//zhou[0].CallCheckB3(tclues, nclf, sn->all_previous_cells, 1);
	Out17(sn->all_previous_cells);
#ifdef DEBUG
	cout << "call expandb3 p_cpt2g[14]" << p_cpt2g[14]
		<< " cl0bf=" << cl0bf 		<< "\tb3direct.nt=" << b3direct.nt 
		<<"  p_cpt2g[4]="<< p_cpt2g[4] << b3direct.nt << " nclf=" << nclf << endl;
#endif
	/*
		return; // stop at first
	*/
	goto next;


	// going back, for a non empty index, count it back
back:
	if (--s >= spb)goto next;


}
void  GCHK::ExpandB3VectV2(uint32_t cl0bf, uint32_t active0) {
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
	if (zhou[0].CallCheckB3(tclues, nclf, sn->all_previous_cells)) {
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			//if (ntadd < 100)
				tadd[ntadd++] = w.bf.u32[2];
			if (cc <= 8 && cc0 < 20) {
				chunkh.Add128(w);
				register uint64_t ua12 = w.bf.u64[0];
				uint64_t ua54 = (ua12 & BIT_SET_27) |
					((ua12 & BIT_SET_B2) >> 5);
				w.bf.u64[0] = ua54;// store in 54 cells mode
				if (cc < 3) {
					w.bf.u32[2] = GetI27(w.bf.u32[2]);
					morevalidc2.Add(w);
				}
				else morevalidothers.Add(w);
			}
		}
		uint32_t ua = zhgxn.tua[0].bf.u32[2];
		if (ispot < limspot) {// new ua for next spot
			sn->possible_cells = ua;
			s = sn;
		}
		else s->possible_cells &= ua;// same spot hitting  new ua			
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
	cout << ws << " one solentry mode"  << endl;
#ifdef TEST_ON
#endif
	//if(zp)strcpy(zp, ws);
	//a_18_seen = 1;
	//aigstop = 1; 

}





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
	mincluesb3 = t16_min_clues[bax[2].i416];
	gchk.UaCollector();
	p_cpt2g[1] = tuasb12.nua;
	p_cpt2g[2] = chunkh.C2Count();
	gchk.BuildVectorsForExpand4B12();
	gchk.Expand4B12();
	morev2c.Status(" b 19-20  ");
	morev2d.Status(" d 21-22  ");
	cout << "count phase2expand=" << p_cpt2g[41]
		<< " " << p_cpt2g[4] << "\n10*av chunks=" << (10 * p_cpt2g[43]) / p_cpt2g[41]
		<< "\nto check=" << p_cpt2g[19] << "\nvalid=" << p_cpt2g[6]
		<< "\nb3 expand=" << p_cpt2g[14] << endl;


	return a_18_seen;

}


//________________ uas generation and store in chunkh
void GCHK::UaCollector() {
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
		cout << " reshaped" << endl;


		for (int ibs = 0; ibs < 3; ibs++) {
			STD_B416 & b = bax[ibs];
			cout << b.band << "\t" << b.i416 << " " << t416n6[b.i416] << endl;
		}
	}
#ifdef HAVEKNOWN
	if (1)
		for (int i = 0; i < 81; i++) cout << ze[i + 82];
	cout << " known reshaped" << endl;
#endif
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
	cout << "final all 2 digits n=" << tuasb12.nua << endl;
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
	cout << "final all 3 digits n=" << tuasb12.nua << endl;
	for (int i = 0; i < 126; i++) {
		int ir = zh2_4[0].GoZ4(floors_4d[i]);
		if (ir < 8) continue;// minimum for a fresh ua 4 digits
		uint64_t F = zh2gxn.unsolved_field;
		tuasb12.GetT2(F);
		//cout << Char9out(floors_4d[i]) << "study ua  i=" << i << endl << endl;
		zh2_4[0].DoZ4Go();
		if (zh2gxn.nua) Adduab12();
	}
	cout << "final  all 4 digits n=" << tuasb12.nua << endl;

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
	cout << "final  all 5 digits n=" << tuasb12.nua << endl;

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

}

//__ start expand 3/4 + no more uas

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
	morev2c.Init(); morev2d.Init();

	//for (uint32_t i = 0; i < 20; i++) {
		//register uint64_t R = tuasb12.tua[i] & BIT_SET_2X;
		//cout << Char2Xout(R) << " i=" << i << endl;

	//}
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
	cout << "entry Expand4B12()" << endl;
	zh2b[0].InitBands12(grid0);
	uint64_t *twu = tuasb12.tua,
		limspot = 2;// expand 4 clues
	nt4_to_expand = 0;
	// _______________ expand to have all  4 cells to expand more
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
	p_cpt2g[3]+= nt4_to_expand;
	for (uint64_t i = 0; i < nt4_to_expand; i++) {
		cout << i << "\t";
		//t4_to_expand[i].dump();
		Do_phase2(t4_to_expand[i]);
	}

}
void GCHK::Do_phase2(T4_TO_EXPAND w) {
	//____ load the first clues
	{
		register uint64_t R = w.bf;
		nclues_step = 0;
		//nclues_step = (int)_popcnt64(w.bf);
		uint32_t xcell;
		while (bitscanforward64(xcell, R)) {
			uint32_t cell = From_128_To_81[xcell];
			R ^= (uint64_t)1 << xcell; //clear bit
			tclues[nclues_step++] = cell;
		}
		tcluesxy = &tclues[nclues_step];
	}
	ntua4 = 0;
	{// build still valid uas
		register uint64_t F = w.bf, A = w.active;
		for (uint32_t i = 0; i < 128; i++) {
			register uint64_t R = tuasb12.tua[i];
			if (R&F) continue;
			R &= A;
			if (!R)return; //dead
			register uint64_t cc = _popcnt64(R);
			AddUA64(tua4, ntua4, R | (cc << 59));
		}
		/*
		if (0) {
			// insert morev2a more v2b if any (small uas )
			uint64_t tt4[1000];
			uint32_t ntt4 = 0;
			morev2a.Add_to(tt4, ntt4);
			morev2b.Add_to(tt4, ntt4);
			for (uint32_t i = 0; i < ntt4; i++) {
				register uint64_t R = tt4[i];
				//cout << Char2Xout(R) << " inserted" << endl;
				if (R&F) continue;
				R &= A;
				if (!R)return; //dead
				tua4[ntua4++] = R;
			}
			// add tua4 to tuasb12.tua
			for (uint32_t i = 0; i < ntt4; i++) {
				register uint64_t R = tt4[i];
				//cout << Char2Xout(R) << " inserted" << endl;
				tuasb12.tua[tuasb12.nua++] = R;
			}
			morev2a.Init();	morev2b.Init();

		}
		*/
		for (uint32_t i = 128; i < tuasb12.nua; i++) {
			register uint64_t R = tuasb12.tua[i];
			if (R&F) continue;
			R &= A;
			if (!R)return; //dead
			tua4[ntua4++] = R;
		}

	}
	//______________________________________________end collect
	cout << Char2Xout(w.bf) << " ntua4=" << ntua4
		<< " ntuas=" << tuasb12.nua << endl;
	//return;
	int wnua = ntua4;
	memset(v12_v0, 0, sizeof v12_v0);
	for (int i = 0; i < 64; i++) {
		if (wnua > 64) {
			v12_v0[i] = ~0; wnua -= 64;
		}
		else {
			v12_v0[i] = maskLSB[wnua].u64[0];
			break;
		}
	}
	for (int i = 0; i < 54; i++)memcpy(v12_c[i], v12_v0, sizeof v12_c[0]);
	for (uint32_t i = 0; i < ntua4; i++) {
		register uint64_t R = tua4[i] & BIT_SET_2X;
		uint32_t xcell, bloc = i >> 6, ir = i - 64 * bloc;
		while (bitscanforward64(xcell, R)) {
			uint32_t cell = From_128_To_81[xcell];
			uint64_t bit = (uint64_t)1 << ir;
			R ^= (uint64_t)1 << xcell; //clear bit
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
	Do_phase2Expand(w.bf);
}
void  GCHK::Do_phase2Expand(uint64_t bf) {
	int mincluesb3 = 6,//<<<<<<<<<<<<<<<<<<<< must be dynamic and linked ot band3 i416
		limspot;
	// _______________ expand to have all  minimal < n clues
	struct SPB {// spots to find starts to extract uas 3 digits
		uint64_t  possible_cells, all_previous_cells, active_cells;
	}spb[20], *s, *sn;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->all_previous_cells = bf;
	s->active_cells = BIT_SET_2X^bf;
	s->possible_cells = tua4[0] & BIT_SET_2X;
	if (1) {// debugging of start
		limspot = 7;// 7+3=10
	}
#define nphase1 3
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
		tcluesxy[ispot] = cell;
		s->active_cells ^= bit;
		uint64_t ac = s->active_cells;
		sn = s + 1; *sn = *s; // (copy the stack count)
		sn->all_previous_cells = s->all_previous_cells | bit;
		// find next ua
		V12_64 *tv1 = tv12_64[ispot], *tv2 = tv12_64[ispot + 1];
		uint32_t nd = 0;
		uint64_t *vc = v12_c[cell];
		for (uint32_t i = 0; i < nv12_64_spot[ispot]; i++) {
			V12_64 w = tv1[i];
			w.v &= vc[w.ind];
			if (w.v) tv2[nd++] = w;
		}
		if (nd) {// still at least one ua not hit
			nv12_64_spot[ispot + 1] = nd;
			V12_64 w = tv2[0];// here is next ua
			uint32_t ir;
			bitscanforward64(ir, w.v);//relative index
			uint64_t Ru = tua4[ir + w.u_start] & ac;
			if (!Ru)goto next;// what does it mean??
			if (ispot >= limspot) {
				p_cpt2g[41]++;
				p_cpt2g[43] += nd;
				goto next;
			}
			sn->possible_cells = Ru;
			s++; // switch to next spot
			goto next;

		}
		// no more uas below count
		//cout << Char2Xout(sn->all_previous_cells) << " " << p_cpt2g[42]<<endl;
		AfterExpandB12(sn->all_previous_cells, s->active_cells,
			(uint32_t)ispot + 1 + nclues_step);

		goto next;
		// going back, for a non empty index, count it back
	back:
		if (--s >= spb)goto next;

		//morev2a.Status(" a 12-16");
		//morev2b.Status(" b 17-18  ");

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

void GCHK::AfterExpandB12(uint64_t bf, uint64_t ac, int ncl) {// no more ua
	//return;
	//if (p_cpt2g[42] > 500) return;
	p_cpt2g[4]++;
	if (morev2a.ApplyXY(tclues, ncl)) return;
	if (morev2b.ApplyXY(tclues, ncl)) return;
	if (morev2c.ApplyXY(tclues, ncl)) return;
	if (morev2d.ApplyXY(tclues, ncl)) return;
	p_cpt2g[5]++;
	nclues = ncl;
	//for (uint32_t i = 0; i < ncl; i++) cout << tclues[i] << " ";
	if (zh2b[1].IsValid(tclues, ncl)) {
		//cout << Char2Xout(bf) << " " << ncl << " " << p_cpt2g[42] << " ";
		//cout << "not valid 2 bands " << zh2gxn.nua << " uas to use "<< p_cpt2g[42] << endl;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i],
				cc = _popcnt64(ua);
			if (cc < 19)morev2a.Add(ua);  //12-20

			//if(cc<17)morev2a.Add(ua);  //12-16
			//else if (cc < 19)morev2b.Add(ua); // 17-18
			else if (cc < 21)morev2c.Add(ua); // 19 20
			else if (cc < 23) morev2d.Add(ua);// 21-22
			// if >22 use it only to solve locally
		}
		// process fresh uas  in a special loop till valid
		return;
	}
	morevalidc2.Init(); morevalidothers.Init();
	p_cpt2g[6]++;
	uint64_t bf54 = (bf & BIT_SET_27) | ((bf & BIT_SET_B2) >> 5),
		ac54 = (ac & BIT_SET_27) | ((ac & BIT_SET_B2) >> 5);
	InitGoB3(bf54, ac54, ncl );// valid puzzle  band3 bands1+2 locked
	if (ncl >= 18 - mincluesb3) return;// 2 is min clues in b3
	p_cpt2g[7]++;
	//_______________ room for add ib bands 1+2 
	if (p_cpt2g[22] < svb12.nt2)p_cpt2g[22] = svb12.nt2;
	if (p_cpt2g[23] < svb12.ntmore)p_cpt2g[23] = svb12.ntmore;
	myb12 = bf54;
	myac = ac54;
	mynclues = ncl;
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
	if (p_cpt2g[7] < 1000 || p_cpt2g[7] > 1100) return;
	cout << Char54out(bf54) << " bf54 ncl=" << ncl << endl;
	cout << Char54out(ac54) << " ac54 ntadd="	<< ntadd << endl;
	cout << "count=\tnt2=" << svb12.nt2 << "\tnmore=" << svb12.ntmore << "\tmin=" << svb12.smin.mincount << endl;
	vaddh.SetUpAdd0(svb12,ac54);
	//cout << Char64out(vaddh.vadd0.vc2)<< " vc2" << endl;
	//for (uint32_t i = 0; i < ntadd; i++)
		//cout << Char64out(vaddh.vaddcell[i].vc2) << " icell=" << i << endl;
	//cout << Char64out(vaddh.vadd0.vcmore[0]) << " vcmore" << endl;
	//for (uint32_t i = 0; i < ntadd; i++)
		//cout << Char64out(vaddh.vaddcell[i].vcmore[0]) << " icell=" << i << endl;
	ExpandAddB1B2();
}

void GCHK::VADD_HANDLER::SetUpAdd0(VB12& vbo, uint64_t ac) {
	//memset(&vadd0, 0, sizeof vadd0);
	memset(vaddcell, 255, sizeof vaddcell);
	// setup vadd0
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
	if (p_cpt2g[26] & 31)return;
	//cout << "entry ExpandAddB1B2Go" << endl;
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
	//cout << Char54out(myb12add) << " stepo= " << step << " " << _popcnt64(myb12add) 	
	//	<<" nmin=" <<nmin<< endl;
	if (nmin > ncluesb3) return;
	//cout << Char27out((uint32_t)vaddh.vaddsteps[step].vc2) << " still ok c2" << endl;
	//cout << Char64out(vaddh.vaddsteps[step].vcmore[0]) << " still ok more" << endl;
	//for(int i=0;i<=step;i++)
	//	cout << Char27out((uint32_t)vaddh.vaddsteps[i].vc2) << " ok i=" << endl;


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
	//cout << "ExpandAddB1B2Go svb12add.ntmore=" << svb12add.ntmore << endl;
	//GoB3(mynclues + step, svb12add);
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

	//____________ here start the search nclues
next:
	if (s->icur > ntadd) goto back;
	uint64_t ispot = s - spb;
	register int cell = tadd[s->icur++];
	register uint64_t bit = (uint64_t)1 << cell;
	tgo[ispot] = cell;
	sn = s + 1; *sn = *s; 
	vaddh.Apply(ispot, s->icur-1);
	sn->all_previous_cells = s->all_previous_cells | bit;
	myb12add = sn->all_previous_cells;
	if (ispot <2) {
		ExpandAddB1B2Go((int)ispot + 1);// call the process for this 
	}
	if (ispot >= lim) 	goto next;
	s++; // switch to next spot (icur+1)
	goto next;
back:
	if (--s >= spb)goto next;
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
	chunkh.Get2(svb12.t2, svb12.nt2, bf54, ac54, tclues, ncl);
	chunkh.GetMore(svb12.tmore, svb12.ntmore, bf54, ac54, tclues, ncl);
	int nmin = svb12.Getsmin();
	ncluesb3 = 18 - ncl;
	nmiss = ncluesb3 - nmin;
	if (nmiss >= 0)GoB3(ncl, svb12);
}
void GCHK::GoB3(  int ncl,  VB12 & vbx) {

	int isdirect = 0;

	// _direct process unless compulsory outfield clues
	if (nmiss > 2)isdirect = 1;
	else {
		return;
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
		for (uint32_t i = 0; i < vbx.ntmore; i++)
			b3direct.Add(vbx.tmore[i].bf.u32[2]);
		ExpandB3Vect(17 - ncl);
		return;
	}
}
void  GCHK::BuildExpandB3Vect( uint32_t cl0bf, uint32_t active0, VB12 & vbx) {// when assigned before
	int cc= _popcnt32(cl0bf);
	if (cc > ncluesb3)return;
	if (cc == ncluesb3) {// nothing to expand
		p_cpt2g[19]++;
		if (zhou[0].CallCheckB3(tclues, nclues, cl0bf)) {
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
	int nn = ncluesb3 -1 - cc;
	b3direct.Init();
	for (uint32_t i = 0; i < vbx.ntg2ok; i++) {
		register uint32_t U = tg2[vbx.tg2ok[i]].pat;
		if(!(U&cl0bf))	b3direct.Add(U);
	}
	for (uint32_t i = 0; i < vbx.ntmore; i++) {
		register uint32_t U = vbx.tmore[i].bf.u32[2];
		if (!(U&cl0bf))	b3direct.Add(U&active0);
	}
	ExpandB3Vect(nn, cl0bf, active0);

}
void  GCHK::ExpandB3Vect(int limspot , uint32_t cl0bf ,uint32_t active0 ) {
	p_cpt2g[14]++;
	uint32_t tadd[100], ntadd = 0;
	// uas basis in b3direct and chunkh.band3[2]
	struct SPB {
		uint64_t vects[3];
		uint32_t  possible_cells, all_previous_cells, active_cells, x;
	}spb[10], *s, *sn;
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
			else {
				bitscanforward64(ir, sn->vects[2]);//ir ua to load
				U = chunkh.band3[1].tua[ir];
			}
		}
		if (ir) {// still at least onet ua
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
	if (zhou[0].CallCheckB3(tclues, nclues, sn->all_previous_cells)) {
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			if (cc < 3 && cc0 < 20) {
				cout << Char2Xout(w.bf.u64[0]) << "\t";
				cout << Char27out(w.bf.u32[2]) << " new ua "
					<< p_cpt2g[6] << " " << ++p_cpt2g[20] << endl;
			}
			tadd[ntadd++] = w.bf.u32[2];
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
		uint32_t ua = zhgxn.tua[0].bf.u32[0];
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
	zhou[0].CallCheckB3(tclues, nclues, sn->all_previous_cells, 1);
	cout << " isvalid18" << endl;
	Out17(sn->all_previous_cells);
	/*
		return; // stop at first
	*/
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

on passe en état final plus de changement en bandes 1+2

			orw = 0;
			for (uint32_t i = 0; i < svb12.ntmore; i++) {
				register uint32_t R3 = svb12.tmore[i].bf.u32[2];
				if (!(R3&infield)){
					cout << Char54out(svb12.tmore[i].bf.u64[0]) << " ";
					cout << Char27out(svb12.tmore[i].bf.u32[2]) << endl;
					orw |= svb12.tmore[i].bf.u64[0];
				}
			}
			if(orw)		cout << Char54out(orw) << " or" << endl;
*/

//=========brute force specific to this
int ZHOU::CallCheckB3(uint32_t * t, int n, uint32_t bf, int nogo) {// 17 search mode
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
	if (nogo) { ImageCandidats(); return 0; }// test
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
int ZHOU::CallMultipleB3x(ZHOU & o, uint32_t bf, int diagx) {
	*this = o;
	BF128 dca[9];
	int digitsbf = zh_g2.digitsbf;
	memcpy(dca, zh_g2.Digit_cell_Assigned, sizeof dca);
	{
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = gchk.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 0;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	zhgxn.nua = 0;
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	//__________end assign last lot start solver
	zh_g.go_back = 0;	// modevalid is set to  1
	//zh_g2.isfalse_on = -1;
	int ir = Full17Update();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0);
	return zhgxn.nua;
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
	//cout << "index=" << index << endl;
	//ImageCandidats();
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
				//cout << "zhgxn.nua=" << zhgxn.nua << endl;
				//wua.Print3(" ");
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





//= old code to kill


//__________________ end of phase 1 process the file and clean
inline int Check128uas(uint32_t ncl, uint32_t *tcl,	BF128  bv, BF128 * cellsv) {
	for (uint32_t i = 0; i < ncl; i++) {//all cells common to the step
		bv &= cellsv[tcl[i]];
	}
	return bv.isNotEmpty();
}

int GCHK::Is_B12_Not_Unique() {
	//myua = zh2b[0].Valid_XY(tcluesxy, nclues);
	if (myua) {//not unique do the best with the UA
		NewUaB12();		return 1;
	}
	return 0;
}


//___________ potential valid bands 1+2  after 128 uas

void GCHK::CleanAll() {

}
int GCHK::Clean_Valid() {
	clean_valid_done = 1;
	myua = 0;// zh2b[0].ValidXY(tclues, nclues + nclues_step);
	if (myua) { NewUaB12();	return 1; }
	return 0;
}


void G17B3HANDLER::Critical2pairs() {// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (smin.mini_bf2) {// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern
		for (int ist = 0; ist < 3; ist++) {
			int shrink = TblShrinkMask[smin.mini_bf2 & (0111 << ist)];
			if (shrink) {// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
			}
		}
		smin.mini_bf2 = 0;
	}
}

void G17B3HANDLER::CriticalLoop() {// after optional assignment
	if (gchk.aigstopxy)return;
	//if (gchk.ua_out_seen)return;
	while (1) {// first shrink uas in field
		irloop = 0;
		uint32_t * tn = &uasb3if[nuasb3if], n = 0;
		register uint32_t Ra = active_b3,
			Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuasb3if; iua++) {
			register int Ru = uasb3if[iua];
			if (Ru & known_b3) continue;// already hit, forget it
			Ru &= active_b3;
#ifdef DBB
			if (diagh) 				cout << Char2Xout(Ru) << " n=" <<n<< endl;

#endif
			if (!Ru) return;// dead branch
			if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
				CriticalAssignCell(Ru);
				Ra = active_b3; //can be  modified
				irloop = 1;// should loop for new singles
			}
			else tn[n++] = Ru;
		}
		uasb3if = tn;		nuasb3if = n;
		if (!n) irloop = 0;// no need to loop again
		if (!irloop) break;
	}
	if (!active_b3) {// must be here expected number of clues
		if (nuasb3if) return; //can not be valid
		gchk.FinalCheckB3(known_b3);
		return; // branch closed
	}
	if (_popcnt64(known_b3) >=gchk.ncluesb3) return; // branch closed
	int wua = uasb3if[0] & active_b3, cell;
	while (bitscanforward(cell, wua)) {
		register int bit = 1 << cell;
		wua ^= bit;// clear bit
		// clean the bit in active_b3, this is now a dead cell downstream
		active_b3 ^= bit;
		G17B3HANDLER hn = *this;
		hn.CriticalAssignCell(bit);
		hn.CriticalLoop();
		if (gchk.aigstop ||gchk.aigstopxy|| gchk.ua_out_seen)return;
	}
}

void GCHK::ExpandB3() {
#ifdef DEBUGKNOWN
	cout << "entry expandb3 for nlues="<<ncluesb3 << endl;
	//for (uint32_t i = 0; i < nuasb3_1; i++)		cout << Char27out(uasb3_1[i])<<endl;
#endif

	struct SPB3 {// spots to find band 3 minimum valid solutions
		// ====================== constant after initialization
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3;
	}spb3[15], *s3, *sn3;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = uasb3_1[0] & s3->active_cells;
	int tcells[15];
	//_______________ assign first forced in stack locked


	//____________ here start the search nclues
next:
	uint64_t ispot = s3 - spb3;
	{// catch and apply cell in bitfields
		register uint32_t cell, p = s3->possible_cells;
		if (!p)goto back;
		bitscanforward(cell, p);
		register int bit = 1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1; *sn3 = *s3; // (copy the stack count)
		sn3->all_previous_cells = filter;
		sn3->active_cells = s3->active_cells = ac;

		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuab3 + 1; i < nuasb3_1; i++) {
			if (uasb3_1[i] & filter)continue;
			if (ispot >= ncluesb3-1) 	goto next;//passing the limit
			sn3->iuab3 = i;
			register uint32_t Ru = uasb3_1[i] & sn3->active_cells;
			if (!Ru)goto next;
			if (ispot == ncluesb3-2) {// last must hit all remaining uas
				for (uint32_t i2 = i + 1; i2 < nuasb3_1; i2++) {
					register uint32_t Ra = uasb3_1[i2];
					if (Ra & filter)continue;
					Ru &= Ra;
					if (!Ru)goto next;
				}
			}
			sn3->possible_cells = Ru;
			s3 = sn3; // switch to next spot
			goto next;
		}
	}
	// this is a possible nclues do final check

	//p_cpt2g[19]++;
	if (!clean_valid_done)
		if (Clean_Valid()) {
			moreuas_b3.Add(0);//lock the call 
			aigstopxy = 1;		return;
		}
//	uint32_t ir = zhou[0].CallCheckB3(tclues, nclues + nclues_step, sn3->all_previous_cells);

	if (zhou[0].CallCheckB3(tclues, nclues + nclues_step, sn3->all_previous_cells)) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		if (nuasb3_1 < 300)uasb3_1[nuasb3_1++] = ua;
		NewUaB3();
		if (ispot < (ncluesb3 - 1)) {// new ua for next spot
			sn3->possible_cells = ua;
			s3 = sn3; 
		}
		else {// continue same spot hitting the new ua
			s3->possible_cells &= ua;
		}
	}
	else {
		Out17(sn3->all_previous_cells);
		return; // stop at first  
	}
	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s3 >= spb3)goto next;
}

uint32_t G17B3HANDLER::IsMultiple(int bf) {
	p_cpt2g[29] ++;
	if (bf == rknown_b3) return 0;
	if (_popcnt32(bf) > 25) return 0;
	uint32_t ua = 0;
	rknown_b3 = bf;
	GCHK & bab = gchk;
	// check first if all tuab3 is hit
	ua = gchk.moreuas_b3.Check(bf);
	if(ua)	return ua;//  send back the ua

	uint32_t ir = zhou[0].CallCheckB3(gchk.tclues, gchk.nclues + gchk.nclues_step, bf);
	//int ir = zhou[1].CallMultipleB3(zhou[0], bf, 0);
	if (ir) {// setup the ua found 
		ua = zh_g2.cells_assigned.bf.u32[2];
		gchk.moreuas_b3.Add(ua);
		gchk.NewUaB3();
	}

	return ua;
}

//================= critical process
void G17B3HANDLER::CriticalAssignCell(int Ru) {// assign a cell within the critical cells
	// Ru is usually a regidster containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_b3 |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & smin.mini_bf3) {// the cell is in a minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		smin.mini_bf3 ^= bit; // now only a pairto hit
		smin.mini_bf1 |= bit;
	}
	else {// either one pair or a triplet in the minirow
		active_b3 &= (~Mask); // kill the minirow as active
		smin.mini_bf1 &= ~bit;
		smin.mini_triplet &= ~bit;
	}
}


void G17B3HANDLER::Go_Critical() {// critical situation all clues in pairs tripl:ets
	active_b3 = smin.critbf;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}


//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow() {
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1 << ndead, mask = 7 << (3 * ndead);
	for (int i = ndead; i < 9; i++, bit <<= 1, mask <<= 3) {
		stack = i % 3;
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (bit & smin.mini_bf1) {// it was a gua2 pair assign both
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf1 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & smin.mini_bf2)// it was 2 gua2 pair assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_bf2 ^= bit;
				hn.SubMini(M, mask);
			}
		else if (bit & smin.mini_bf3) {// it was 3 gua2 pair assign 3 out of 3
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf3 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & smin.mini_triplet)// it was a gua3 triplet assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_triplet ^= bit;
				hn.SubMini(M, mask);
			}
		else {// second add in the mini row one residual cell take it
			G17B3HANDLER hn = *this;
			hn.SubMini(M, mask);
		}
	}
}
void G17B3HANDLER::SubMini(int M, int mask) {
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added
	active_b3 &= ~mask;
	active_sub ^= M;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}
void G17B3HANDLER::Go_Subcritical() {// nmiss to select in the critical field
	active_b3 = active_sub = smin.critbf;
	ndead = 0;
	Go_SubcriticalMiniRow();// find the first miss
}

//======================================================================= not critical sequence
void G17B3HANDLER::ShrinkUasOfB3() {
	if (known_b3) {// shrink the out field table
		p_cpt2g[53] ++;
		uint32_t * tn = &uasb3of[nuasb3of], n = 0;
		andoutf = BIT_SET_27;
		register uint32_t Ra = wactive0,
			Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuasb3of; iua++) {
			register int Ru = uasb3of[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			andoutf &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
		uasb3of = tn;
		nuasb3of = n;
	}
}

void G17B3HANDLER::Go_miss1_b3() {// not called if more than 1 needed
	wua = wactive0;
	if (known_b3 &&nuasb3of)ShrinkUasOfB3();// if known from up stream
	if (smin.mincount) {
		if (!nuasb3of) {// subcritical in hn if solved
			int uabr = IsMultiple(known_b3 | active_b3);
			if (uabr) {// on ua outfield seen
				p_cpt2g[56] ++;
				andoutf = uabr;
				nuasb3of = 1;
			}
			else {// confirmed subcritical possible
				p_cpt2g[57] ++;
				G17B3HANDLER hn = *this;
				hn.Go_Subcritical();
			}
		}
	}
	//if (1) return;
	if (nuasb3of) wua &= andoutf;
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		p_cpt2g[58] ++;
		while (bitscanforward(res, wua)) {
			p_cpt2g[59] ++;
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;
			if (hn.AddCell_Of(res, bit))hn.Go_Critical();
		}
	}
}
void G17B3HANDLER::Go_miss2_b3() {
	ShrinkUasOfB3();// if known from upstream
	if (smin.mincount) {
		if (!nuasb3of) {// subcritical in hn if solved
			int uabr = IsMultiple(known_b3 | active_b3);
			if (uabr) {// on ua outfield seen
				uasb3of[0] = uabr;
				nuasb3of = 1;
			}
			else {// confirmed subcritical possible
				G17B3HANDLER hn = *this;
				hn.Go_Subcritical();
			}
		}
	}
	wua = wactive0;
	if (nuasb3of)wua &= uasb3of[0];// use first ua
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;
			if (hn.AddCell_Of(res, bit))hn.Go_miss1_b3();
		}
	}
}
void G17B3HANDLER::Go_miss3_b3() {// always first direct entry
	ShrinkUasOfB3();// if known from upstream
	if (smin.mincount) {
		if (!nuasb3of) {// subcritical in hn if solved
			int uabr = IsMultiple(known_b3 | active_b3);
			if (uabr) {// on ua outfield seen
				uasb3of[0] = uabr;
				nuasb3of = 1;
			}
			else {// confirmed subcritical possible
				G17B3HANDLER hn = *this;
				hn.Go_Subcritical();
			}
		}
	}
	wua = wactive0;
	if (nuasb3of)wua &= uasb3of[0];// use first ua
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;
			if (hn.AddCell_Of(res, bit))hn.Go_miss2_b3();
		}
	}
}

//_________ final called by all branches 
void GCHK::FinalCheckB3(uint32_t bfb3) {
#ifdef DBB
	if (diagbugchk) {
		cout<<Char27out(bfb3) << "call final diagbug active "  << endl;
		
	}

#endif	
	p_cpt2g[14]++;
	if (aigstopxy) return;
	if ((int)_popcnt32(bfb3) >ncluesb3) 		return;	
	if (!clean_valid_done) {
		if (Clean_Valid()) {
			moreuas_b3.Add(0);//lock the call 
			aigstopxy = 1;		return;
		}
	}
	if (moreuas_b3.Check(bfb3))return;
	p_cpt2g[15]++;
	register uint32_t ir = zhou[0].CallCheckB3(tclues, nclues + nclues_step, bfb3);
#ifdef DBB
	if (diagbugchk) 	cout << "call diagbug active  ir="<<ir<<" p_cpt2g[14]="<< p_cpt2g[14] << endl;
#endif
	if (ir) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		if (ua && (!(ua & hh0.smin.critbf)))	ua_out_seen = ua;
		NewUaB3();
		return;
	}

	Out17(bfb3);
}
void GCHK::Out17(uint32_t bfb3) {
	// mapping of the output on band_order[ib] 
	uint32_t map[81];
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 27; j++)
		map[27 * i + j] = 27 * band_order[i] + j;
	char ws[82];

	uint32_t tcf[40], ntcf = 0;
	for (int i = 0; i < nclues ; i++) {
		tcf[ntcf++] = map[tclues[i]];
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		tcf[ntcf++] = map[54 + i] ;
	strcpy(ws, empty_puzzle);
	for (uint32_t i= 0; i < ntcf; i++) {
		int cell = tcf[i];
		ws[cell] = ze[cell] ;
	}
	cout << ws << " one sol in entry mode p_cpt2g[58]=" << p_cpt2g[58]
		<< "   p_cpt2g[18]=" << p_cpt2g[18] << endl;
	cout << Char2Xout(wb12bf) << " b12  ip=" << start_perm << endl;
#ifdef TEST_ON
#endif
	//if(zp)strcpy(zp, ws);
	//a_18_seen = 1;
	//aigstop = 1; 

}
void GCHK::NewUaB12() {
}


void GCHK::NewUaB3() {// new ua from final check zh_g2.cells_assigned
	BF128 ua128 = zh_g2.cells_assigned;
	register uint64_t ua12 = ua128.bf.u64[0];
	register uint32_t ua = ua128.bf.u32[2],
		cc = _popcnt32(ua),
		cc0 = (uint32_t)_popcnt64(ua12);
	if (!(ua&hh0.smin.critbf))ua_out_seen = ua;
#ifdef DBB
	if (diagbugchk) {
		cout << Char27out(ua) << "uab3 to add of="<< ua_out_seen << endl;
		cout << Char2Xout(ua12) << " ua12 cc=" << cc0 << endl;
		cout << tguas.nvg2 << " " << tguas.nvg3 << " " << g4t_start.ntua4 << endl;
	}

#endif	
	moreuas_b3.Add(ua);
	if (cc0 > GUALIMSIZE+1) return;
	if (cc > 3) {
		if ((cc0 + cc) > GUALIMSIZE + 2) return;
	}
	uint64_t ua54 = (ua12 & BIT_SET_27) | ((ua12 & BIT_SET_B2) >> 5);

	// find the digits pattern from the current band 3
	int * cur_b3 = &gchk.grid0[54], wdigs = 0, c27;
	{
		register uint32_t wua = ua;
		while (bitscanforward(c27, wua)) {// look for  possible cells
			wua ^= 1 << c27;// clear bit
			wdigs |= 1 << cur_b3[c27];
		}
	}



}


void GCHK::Debugifof() {
#ifdef TEST_ON
	cout << "in field out field status) "<< nuasb3_1<<" "<< nuasb3_2 << endl;
	for (uint32_t i = 0; i < 15; i++)
		cout << Char27out(uasb3_1[i] )<< " if " << i << endl;
	for (uint32_t i = 0; i < nuasb3_2; i++)
		cout << Char27out(uasb3_2[i]) << " of " << i << endl;
#endif

}

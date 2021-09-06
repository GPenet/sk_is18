
//_______________ start  a band3 perm and UAs GUAs search 
void BINDEXN::Copy(STD_B1_2 & b) {
	ntvb=b.nvalidb;
	for (uint32_t i = 0; i < ntvb; i++) tvb[i] = b.my_validb[i];
	nt2 = b.nbi2;
	memcpy(t2, b.my_bi2, nt2 * sizeof t2[0]);
}
void BINDEXN::Copy_no7clues(STD_B1_2 & b) {
	nt2 = ntvb = 0;
	VALIDB *vb = tvb;
	for (uint32_t i = 0; i < b.nbi2; i++) { // all index 2
		BI2_32 wi = b.my_bi2[i], win = wi;
		win.istart =ntvb;
		for (uint32_t j = wi.istart; j < wi.iend; j++) {
			VALIDB wj = b.my_validb[j];
			if (_popcnt32(wj.bf) < 7)vb[ntvb++] = wj;
		}
		if (ntvb != win.istart) {// group with <7 clues
			win.iend = ntvb;
			t2[nt2++] = win;
		}
	}		
}
void GCHK::Copy_Check_7clues_needed() {
	int minb1 = t16_min_clues[myband1.i416],
		minb2 = t16_min_clues[myband2.i416];
	if(minb2<5)bin_b1.Copy(myband1);
	else bin_b1.Copy_no7clues(myband1);
	if (minb1 < 5)bin_b2.Copy(myband2);
	else bin_b2.Copy_no7clues(myband2);
}
void GCHK::Start(STD_B416 * bandx, int * tsort, int ip) {
	if (aigstop)return;
	int tperm18[3][3] = { {0,1,2},{0,2,1} , { 1,2,0 } };// keep  increasing order
#ifdef TEST_ON
	cout << "start ip=" << ip << endl;
#endif
	start_perm = ip; //needed later
	int * tp = tperm18[ip];
	if (aigstop) return;
	memset(gchk.kpfilt, 0, sizeof gchk.kpfilt);// debugging code for known
	{   //___________________________________ build the new puzzle and mode
		tpw = tp;// to print later in the right order
		tsortw = tsort;
		for (int ib = 0; ib < 3; ib++) {
			int iperm = tp[ib], iband = tsort[iperm] & 7;
			band_order[ib] = iband;
			bands_abc[ib] = &bandx[iband];
			memcpy(&zsol[27 * ib], bands_abc[ib]->band, 27);
#ifdef DEBUGKNOWN
			// in known chekc mode, morph the known
			puzknown_perm.bf.u32[ib] = puzknown.bf.u32[iband];
#endif
		}

#ifdef DEBUGKNOWN
		puzknown_perm.Print3("known of the perm");
		uint64_t cc = _popcnt64(puzknown_perm.bf.u64[0]);
		uint32_t	cc1= _popcnt32(puzknown_perm.bf.u32[0]),
				cc2 = _popcnt32(puzknown_perm.bf.u32[0]);
		if (cc > 12|| cc1>7 || cc2>7 ||	(ip && cc>11)
			|| (cc==12 && cc1!=6)) {
			cout << "known not seenhere" << endl;
			return;
		}

		
#endif
		// _____________ init bands 1+2 status 

		// set up ordered solution
		//______________ setup brute force initial status
		zsol[81] = 0;
		zh_g.modevalid = 1;
		zh_g2.grid0 = grid0;
		zh_g2.zsol = zh_g2.stdfirstsol;
		//zh_g5.Init(zsol);// reinit brute force
		for (int i = 0; i < 81; i++)grid0[i] = zsol[i] - '1';
		memcpy(genb12.zsol, zsol, sizeof zsol);
		memcpy(genb12.grid0, grid0, sizeof grid0);
		memcpy(zh2b_g.puz0, grid0, sizeof zh2b_g.puz0);
		//___________ setup bands 1 2
		memset(&myband1, 0, sizeof myband1);
		memcpy(&myband1, bands_abc[0], sizeof bax[0]);

		memset(&myband2, 0, sizeof myband2);
		memcpy(&myband2, bands_abc[1], sizeof bax[1]);

		memset(&myband3, 0, sizeof myband3);
		memcpy(&myband3, bands_abc[2], sizeof bax[2]);
		myband3.EndInitBand3();

		memcpy(zh2b_g.puzc, myband1.band, 27);
		memcpy(&zh2b_g.puzc[27], myband2.band, 27);
		zh2b_g.puzc[54] = 0;
		//____________ setup working index for external loop
		Copy_Check_7clues_needed();

#ifdef DEBUGKNOWN
		// check if the valids are there
		int aig = 1;
		for (uint32_t i1 = 0; i1 < bin_b1.ntvb; i1++) {
			uint32_t w1 = vab1w[i1].bf;
			if (w1 == puzknown_perm.bf.u32[0]) {
				aig = 0;
				cout << "seen valid b1 i1=" << i1 << endl;
			}
		}
		if (!aig) {
			for (uint32_t i2 = 0; i2 < bin_b2.ntvb; i2++) {
				uint32_t w2 = vab2w[i2].bf;
				if (w2 == puzknown_perm.bf.u32[1]) {
					aig = 0;
					cout << "seen valid b2 i2=" << i2 << endl;
				}
			}
		}
		else { // band 1 failed print index
			for (uint32_t ii = 0; ii < bin_b2.nt2; ii++) {
				BI2_32 wi = bin_b2.t2[ii];
				uint32_t bf = wi.bf;
				if ((bf & puzknown_perm.bf.u32[0]) == bf) {
					cout << "this is the expected 2 clues index" << endl;
					cout << Char27out(wi.bf) << " start" << wi.istart << " endt" << wi.iend << endl;
					for (uint32_t i = wi.istart; i < wi.iend; i++)
						cout << Char27out(vab1w[i].bf) << endl;
					return;
				}
			}
			cout << "failed to see the expected  valid b1  ntvb=" << bin_b1.ntvb  << endl;
			return;
		}
		if (aig) {
			cout << "failed to see the expected  valid b1 B2 " 
				<< " bin_b1.ntvb=" << bin_b1.ntvb << " bin_b2.ntvb=" << bin_b2.ntvb << endl;

			return;
		}

#endif

		genb12.ValidInitGang();// set up gang bandc and gang bands ab
		genb12.BuildGang9x3();// setup ordered gang digits band C
		zh2b[0].Init_std_bands();
		Go1_Collect_Uas();// next task in the process 
	}

}
void G4TABLE::DebugFout(  ) {
	for (uint32_t i = 0; i < ntua4; i++) {
		fout1 << gchk.start_perm << " "  ;
		fout1 << Char27out(tua4[i].bf.u32[2]) << "\t";
		fout1 << Char2Xout(tua4[i].bf.u64[0]) << " ";
		fout1 << _popcnt64(tua4[i].bf.u64[0]) << endl;
	}
}

void GCHK::Go1_Collect_Uas() {// catch uas and guas 
	//=========================== collect UAs  old process 
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	if (genuasb12.Initgen()) return;
	Go1_Build_V12();

	// _____ GUAs 
	memset(&gvcells, 255, sizeof gvcells); // all gua cells vectors to no hit
	memset(&gvs_start, 0, sizeof gvs_start);
	memset(&gvs_b1, 0, sizeof gvs_start);
	memset(&gvs_b2, 0, sizeof gvs_start);
	g4t_start.Init();
	tguas.InitStart();
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more

	// process dopied from 17 search
	genb12.SecondSockets2Setup();// collect GUA2s 
	genb12.SecondSockets3Setup();// collect GUA3s 

	Go1_Stacks_Uas();// collect stack uas in band 3 in band 3
	//tguas.DebugStart(1);	g4t_start.Debug(1);
#ifdef TEST_ON
	if (TEST_ON) {// debug special
		char z[82]; z[81 ]= 0; 
		memcpy(z, myband1.band, 27);
		memcpy(&z[27], myband2.band, 27);
		memcpy(&z[54], myband3.band, 27);
		z[81] = 0;
		cout << z << " processed" << endl;
		//g4t_start.Debug(0);
	}

#endif
	Go2_Ext_Loop();
}

void STD_B3::InitStack(int i16, int * z0, BANDMINLEX::PERM & p, int iband) {
	i416 = i16;
	int dstack = 3 * iband;
	GetUAs();
	// create the cell map in output
	for (int i = 0; i < 3; i++) {
		int vr = 9 * p.rows[i], vr0 = 9 * i;
		for (int j = 0; j < 9; j++)
			map[vr0 + j] = vr + p.cols[j];
	}
	// morph all uas
	for (uint32_t i = 0; i < nua; i++) {
		register int uao = tua[i] & BIT_SET_27;
		BF128  ua; ua.SetAll_0();
		for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)
			if (bit & uao) {// set the bit in the band form
				int cell = C_transpose_d[map[i]];// +dstack;
				ua.Set_c(cell + dstack);
			}
		int cc3 = _popcnt32(ua.bf.u32[2]);
		uint64_t cc12 = _popcnt64(ua.bf.u64[0]);
		if (cc3 < 2) {// could be missing in the collection of uas
			if (cc12 < 12) continue; // must be there
			genuasb12.AddUACheck(ua.bf.u64[0] | (cc12 << 59));
			continue;
		}
		if (cc3 > 3 ) {
			g4t_start.tua4[g4t_start.ntua4++] = ua;
			continue;
		}
		if (cc12 <= 6) continue; // must be there

		if (cc3 == 2) {// be sure to have gua 2 clues in table
			uint32_t my_i27 = myband3.GetI27(ua.bf.u32[2]);
			tguas.tgua_start[my_i27].Adduacheck(ua.bf.u64[0] | (cc12 << 59));
		}
		if (cc3 == 3) {
			uint32_t my_i9 = gchk.GetI9(ua.bf.u32[2])+27;
			tguas.tgua_start[my_i9].Adduacheck(ua.bf.u64[0] | (cc12 << 59));
		}
	}
}

void GCHK::Go1_Stacks_Uas() {// band 3 uas used as gangsters via  C_transpose_d[81]
	STD_B3 wbs;
	//g4t_start.ntua4 = 0;
	int zt[81];
	for (int i = 0; i < 81; i++) {
		zt[i] = genb12.grid0[C_transpose_d[i]];
	}
	BANDMINLEX::PERM perm_ret;
	//___ stack 1
	bandminlex.Getmin(zt, &perm_ret);
	int ib1 = perm_ret.i416;
	wbs.InitStack(ib1, zt, perm_ret, 0);
	//___ stack 2
	bandminlex.Getmin(&zt[27], &perm_ret);
	int ib2 = perm_ret.i416;
	wbs.InitStack(ib2, &zt[27], perm_ret, 1);
	//___ stack 3
	bandminlex.Getmin(&zt[54], &perm_ret);
	int ib3 = perm_ret.i416;
	wbs.InitStack(ib3, &zt[54], perm_ret, 2);

}

void GCHK::Go1_Build_V12() {// Vector for the first 128 
	register uint64_t * tua = genuasb12.tua;
	register uint32_t nua = genuasb12.nua;
#ifdef TEST_ON
	cout << "nua=" << nua << endl;
#endif	
	ntua12_r = 0;
	uint32_t ntua_128 =nua;
	if (nua > 128) {// store others in ntua12_r
		ntua12_r = nua - 128;
		ntua_128 = 128;
		memcpy(tua12_r, &tua[128], ntua12_r * sizeof tua12_r[0]);
	}

	v128uas = maskLSB[ntua_128];// Uas vector
	memset(vc128, 255, sizeof vc128);// all bits to 1
	uint32_t cc64;// build cells vectors A
	for (uint32_t i = 0; i < ntua_128; i++) {
		register uint64_t Rw = tua[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc128[From_128_To_81[cc64]].clearBit(i);
		}
	}
	
}
//________________  external loop 
void BINDEXN::Copy(BINDEXN & b) {
	ntvb = b.ntvb;
	// for (uint32_t i = 0; i < ntvb; i++) tvb[i] = b.tvb[i];
	memcpy(tvb, b.tvb, ntvb * sizeof tvb[0]);
	nt2 = b.nt2;
	memcpy(t2, b.t2, (nt2+1) * sizeof t2[0]);
}

uint64_t GCHK::FindSockets(uint64_t active, uint64_t lim) {
	nextl1 = 0, nextl2 = 0;
	uint64_t *t = genuasb12.tua;
	uint32_t n = genuasb12.nua;
	for (uint32_t i = 0; i < n; i++) {
		register uint64_t U = t[i] & active,
			U1 = U & BIT_SET_27, U2 = U & BIT_SET_B2;
		if ((!U1) || (!U2))continue;
		uint64_t n = _popcnt64(U& BIT_SET_27);
		if ((lim==2 && n < lim)|| n == lim) {
			uint32_t i1, aig = 1;

			for (i1 = 0; i1 < nextl1; i1++)
				if (U1 == extl1[i1].bfx) { aig = 0; break; }
			if (aig) {
				if (nextl1 < 6) { // open a new window (limit 10
					i1 = nextl1++; aig = 0;
					extl1[i1].Init((uint32_t)U1, 1);
				}
			}
			if (!aig)if (extl1[i1].ntbfy < 20)
				extl1[i1].tbfy[extl1[i1].ntbfy++] = (uint32_t)(U2 >> 32);

		}
		n = _popcnt64(U& BIT_SET_B2);
		if ((lim == 2 && n < lim) || n == lim) {
			register uint32_t U2a = (uint32_t)(U2 >> 32);
			//if (nt2sk < 5)t2sk[nt2sk++] = (uint32_t)(U2 >> 32);
			uint32_t i2, aig = 1;
			for (i2 = 0; i2 < nextl2; i2++)
				if (U2a == extl2[i2].bfx) { aig = 0; break; }
			if (aig) {
				if (nextl2 < 6) { // open a new window (limit 10
					i2 = nextl2++; aig = 0;
					extl2[i2].Init(U2a, 2);
				}
			}
			if (!aig)if (extl2[i2].ntbfy < 20)
				extl2[i2].tbfy[extl2[i2].ntbfy++] = (uint32_t)U1;
		}
	}
#ifdef DEBUGEXL
	//cout << "FindSockets lim=" << lim << "\tnt1=" << nextl1 << "\tnt2=" << nextl2 << endl;
#endif
	if (nextl1 | nextl2) {
		//for (uint32_t i = 0; i < nextl1; i++) extl1[i].Debug();
		//for (uint32_t i = 0; i < nextl2; i++) extl2[i].Debug();
		return lim-1;
	}
	else	return 0;
}
void GetCpt_6clues(uint64_t & n_no, uint64_t & nt, VALIDB *vb, uint32_t nvb, uint32_t bf) {
	n_no = nt = 0;
	register uint32_t F = bf;
	for (uint32_t i = 0; i < nvb; i++) {
		VALIDB &w = vb[i];
		if (w.nval == 4) {
			nt++;
			if (!(bf&w.bf))n_no++;
		}
	}
}
void GetCpt_6clues(uint64_t & n_no, uint64_t & nt, VALIDB *vb, uint32_t nvb, uint32_t *tbf, uint32_t ntbf) {
	n_no = nt = 0;
	for (uint32_t i = 0; i < nvb; i++) {
		VALIDB &w = vb[i];
		if (w.nval == 4) {
			nt++;
			uint32_t aig = 1;
			register uint32_t U = w.bf;
			for (uint32_t j = 0; j < ntbf; j++) {
				if ((tbf[j] & U)) { aig = 0; break; }
			}
			n_no += aig;
		}
	}
}
void GCHK::ExtractMin(uint64_t active, BINDEXN & bin1, BINDEXN & bin2) {
	uint64_t limb1 = _popcnt64(active & BIT_SET_27) - 3,
		limb2 = _popcnt64(active & BIT_SET_B2) - 3;
	for (uint32_t i = 0; i < nextl1; i++) {
		uint32_t orw = 0,*t= extl1[i].tbfy;
		for (uint32_t j = 0; j < extl1[i].ntbfy; j++) {
			orw |= t[j];
		}
		if (_popcnt32( orw)>(uint32_t)limb1) continue;//min 3 clues will end as 100% ratio
		GetCpt_6clues(n_nob1, ntotb1, bin1.tvb, bin1.ntvb, extl1[i].bfx);

		if (extl1[i].ntbfy == 1)
			GetCpt_6clues(n_nob2, ntotb2, bin2.tvb, bin2.ntvb, extl1[i].tbfy[0]);
		else GetCpt_6clues(n_nob2, ntotb2, bin2.tvb, bin2.ntvb, extl1[i].tbfy, extl1[i].ntbfy);
		uint64_t ratio1 = (ntotb1 - n_nob1)*ntotb2 + (ntotb2 - n_nob2)* n_nob1,
			ratio = (uint64_t)100 * ratio1 / ntotb1 / ntotb2;

		if (ratio < minratio) {
			minratio = ratio;
			extl1[i].ratio= (uint32_t )ratio;
			extlw = extl1[i];
			extlw.noxyes = (ntotb2 - n_nob2)* n_nob1;
#ifdef DEBUGEXL
			//cout << Char27out(extl1[i].bfx) << " t1";
			//cout << "\t" << n_nob1 << "\t" << ntotb1 << "\t" << n_nob2 << "\t" << ntotb2
			//	<< "\t =" << ratio1 << "\t" << ntotb1 * ntotb2 << "\t" << ratio << endl;
#endif
		}

	}
	for (uint32_t i = 0; i < nextl2; i++) {
		uint32_t orw = 0, *t = extl2[i].tbfy;
		for (uint32_t j = 0; j < extl2[i].ntbfy; j++) {
			orw |= t[j];
		}
		if (_popcnt32(orw ) > (uint32_t)limb2) continue;// will end as 100% ratio
		GetCpt_6clues(n_nob2, ntotb2, bin2.tvb, bin2.ntvb, extl2[i].bfx);

		if (extl2[i].ntbfy == 1)
			GetCpt_6clues(n_nob1, ntotb1, bin1.tvb, bin1.ntvb, extl2[i].tbfy[0]);
		else GetCpt_6clues(n_nob1, ntotb1, bin1.tvb, bin1.ntvb, extl2[i].tbfy, extl2[i].ntbfy);
		uint64_t ratio1 = (ntotb2 - n_nob2)*ntotb1 + (ntotb1 - n_nob1) * n_nob2,
			ratio = (uint64_t)100 * ratio1 / ntotb1 / ntotb2;
		if (ratio < minratio) {
			minratio = ratio;
			extl2[i].ratio =(uint32_t )ratio;
			extlw = extl2[i];
			extlw.noxyes = (ntotb1 - n_nob1) * n_nob2;
#ifdef DEBUGEXL
			//cout << Char27out(extl2[i].bfx) << " t2";
			//cout << "\t" << n_nob1 << "\t" << ntotb1 << "\t" << n_nob2 << "\t" << ntotb2
			//	<< "\t =" << ratio1 << "\t" << ntotb1 * ntotb2 << "\t" << ratio << endl;
#endif
		}
	}
}

void GCHK::ExtSplitY(BINDEXN & binw,uint32_t *tbf, uint32_t ntbf, 
	uint32_t & activer) {// source binw exit yes in binw
	uint32_t lim_source = binw.nt2;// table no is also source
	binw.nt2 = binw.ntvb = 0;
	VALIDB *vb = binw.tvb;
	activer = 0; // or status in new vb1 (no) start null

	for (uint32_t i = 0; i < lim_source; i++) { // all index 2
		BI2_32 wi = binw.t2[i], win = wi;
		win.istart = binw.ntvb;
		for (uint32_t j = wi.istart; j < wi.iend; j++) {
			VALIDB wj = binw.tvb[j];
			register uint32_t Bf = wj.bf;

			for (uint32_t k = 0; k < ntbf; k++) {
				if (Bf & tbf[k]) {// a new "yes" 
					vb[binw.ntvb++] = wj;
					activer |= wj.bf;
					break;// one hit is enough 
				}
			}
		}
		if (binw.ntvb != win.istart) {// new group in "yes"
			win.iend = binw.ntvb;
			binw.t2[binw.nt2++] = win;
		}
	}
}
void  GCHK::ExtSplitX(BINDEXN & bin1no, BINDEXN & bin1yes,
	uint32_t bf, uint32_t & activer) {
	uint32_t lim_source = bin1no.nt2;// table no is also source
	// initial status empty for table yes and no
	bin1no.nt2 = bin1no.ntvb = 0;
	bin1yes.nt2 = bin1yes.ntvb = 0;
	activer = 0; // or status in new vb1 (no) start null
	VALIDB *vb1 = bin1no.tvb, *vb2 = bin1yes.tvb;

	register uint32_t F = bf; // UA to hit to be yes
	for (uint32_t i = 0; i < lim_source; i++) { // all index 2
		BI2_32 wi = bin1no.t2[i], win = wi, win1;
		win.istart = bin1yes.ntvb;
		if (F&wi.bf) {// all the group is "yes"
			uint32_t n = wi.iend - wi.istart;
			win.iend = bin1yes.ntvb + n;
			memcpy(&vb2[bin1yes.ntvb], &vb1[wi.istart], n * sizeof  vb2[0]);
			bin1yes.t2[bin1yes.nt2++] = win;
			bin1yes.ntvb += n;
			continue;
		}
		// must check each validb of the group
		win1 = wi;
		win1.istart = bin1no.ntvb;
		for (uint32_t j = wi.istart; j < wi.iend; j++) {
			VALIDB wj = vb1[j];
			if (F&wj.bf) 	vb2[bin1yes.ntvb++] = wj;
			else { vb1[bin1no.ntvb++] = wj; activer |= wj.bf; }
		}
		if (bin1no.ntvb != win1.istart) {// new group in "no"
			win1.iend = bin1no.ntvb;
			bin1no.t2[bin1no.nt2++] = win1;
		}
		if (bin1yes.ntvb != win.istart) {// new group in "yes"
			win.iend = bin1yes.ntvb;
			bin1yes.t2[bin1yes.nt2++] = win;
		}
	}
	
}


int CheckBf(BINDEXN & binw, uint32_t bfw) {
	register VALIDB * tvb = binw.tvb;
	for (uint32_t i = 0; i < binw.ntvb; i++) {
		if (tvb[i].bf == bfw)return i;
	}
	return -1;
}

/*external loop control
 the external loop cut the process using small UAs bands 1+2
 if 'yes' is the part hiiting a ua in a band
 yes1 hit in band 1 ; yes 2 hit in band2
 the process can then be cut in 2 chunks
 yes1 * all2		all1-yes1 * yes2 
 this is of interest if the count is significantly lower than 
		all1 * all 2
 and if all1*all2 is big enough

 the process can be split several times with different uas

 here, the first chunk yes * all is split again
 the main loop continues to split    all1-yes1 * yes2 
 Note the first band can be band 1 or band 2

 Control of the main loop

 EXLRATIO is the minimal wanted reduction of the count
 EXLNLOOP1 the maximum number of steps in the main loop
 EXLNLOOP2 the maximum number of steps in the first chunk split
 EXLBREAK the count minimum to loop
 */
#define EXLRATIO 80
#define EXLNLOOP1 4
#define EXLNLOOP2 3
#define EXLBREAK 4000000

void GCHK::Go2_Ext_Loop() {	//_____________ outer loop 
#ifdef DEFPHASE
	if (DEFPHASE == -1)return;
#endif

	loopb1 = 0;
	uint32_t activerb1, activerb2;
	uint64_t activeloop = BIT_SET_2X;
	activerb1= activerb2 = BIT_SET_27;

	//_________________ external loop
	while (++loopb1 << EXLNLOOP1) {
#ifdef DEBUGEXL
		if (aigstop)cout << "seen aigstop=1" << endl;
		cout << Char2Xout(activeloop) << "========= loop " << loopb1 << endl;
#endif
		if (aigstop)return;;
		minratio = extlr.ratio=1000;
		uint64_t ir = FindSockets(activeloop,2);
		if (ir) 	ExtractMin(activeloop, bin_b1, bin_b2); 	
		if (!ir || minratio > EXLRATIO) {
			ir = FindSockets(activeloop,3);
			if (ir)ExtractMin(activeloop, bin_b1, bin_b2);
			if (!ir || minratio > EXLRATIO) {
				ir = FindSockets(activeloop,4);
				if (ir)ExtractMin(activeloop, bin_b1, bin_b2);
			}
		}
#ifdef DEBUGEXL
		//else 	cout << "finf2 minratio= " << minratio << endl;
		cout << "final selection minratio=" << minratio << endl;
#endif
		if (minratio > EXLRATIO) break;
		else {
			extlr = extlw;
#ifdef DEBUGEXL
			cout << "entry count=" << extlw.noxyes << endl;
			extlw.Debug();
			cout << "shrink mode=" << extlw.mode << endl;
#endif
			if (extlw.mode == 1) {// this is a band1 X band2 Y
				ExtSplitX(bin_b1, bin_b1yes, extlw.bfx, activerb1);

#ifdef DEBUGKNOWN
#ifdef TEST_ON
				cout << "kpfilt[0]=" << kpfilt[0] << endl;
#endif

				if (!kpfilt[0]) {// forget if the first is there
					if (extlw.bfx & puzknown_perm.bf.u32[0]) {//hit
						kpfilt[0] = 1;// this must be in this ip
#ifdef TEST_ON
						cout << " loopb1 mode 1 to use for known=" << loopb1 << endl;
						cout << Char27out(extlw.bfx) << " extlw.bfx b1" << endl;
#endif
						if (loopb1 == 1)
							Go2b_Ext_Loop(BIT_SET_2X , 1);
						else  Go3(bin_b1yes, bin_b2);

					}
				}
#else

				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X | activerb1, 1);
				else  Go3(bin_b1yes, bin_b2);

#endif
				if (aigstop)return;
				ExtSplitY(bin_b2,			extlr.tbfy, extlr.ntbfy, activerb2);
			}
			else {// this is a band2 X band1 Y
				ExtSplitX(bin_b2, bin_b2yes,	extlw.bfx, activerb2);
#ifdef DEBUGKNOWN
				if (!kpfilt[0]) {// forget if the first is there
					if (extlw.bfx & puzknown_perm.bf.u32[1]) {//hit
						kpfilt[0] = 1;// this must be in this ip
#ifdef TEST_ON
						cout << " loopb1 mode 2 to use for known=" << loopb1 << endl;
						cout << "\t\t" << Char27out(extlw.bfx) << " extlw.bfx b2" << endl;
#endif
						if (loopb1 == 1)
							Go2b_Ext_Loop(BIT_SET_2X, 2);
						else  Go3(bin_b1, bin_b2yes);

					}
				}
#else

				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X, 2);
				else  Go3(bin_b1, bin_b2yes);

#endif
				ExtSplitY(bin_b1,	extlr.tbfy, extlr.ntbfy,  activerb1);
			}
		}
		activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
		if (extlr.noxyes < EXLBREAK) break;
	}
	Go3(bin_b1, bin_b2);// last call 

}
void GCHK::Go2b_Ext_Loop(uint64_t activeloop, uint32_t mode2) {
	//___________init the working areas bin2_b1, bin2_b2
#ifdef DEBUGEXL
	cout << Char2Xout(activeloop) << "entry Go2b_Ext_Loop mode " << mode2 << endl;
#endif
	if (mode2 == 1) {// b1 yes b2 all
		bin2_b1.Copy(bin_b1yes);
		bin2_b2.Copy(bin_b2);
	}
	else {// b1 all b2 yes
		bin2_b1.Copy(bin_b1);
		bin2_b2.Copy(bin_b2yes);
#ifdef DEBUGKNOWN
		if (CheckBf(bin_b1, puzknown_perm.bf.u32[0]) < 0) {
			cout << "entry first not seen b1 valid to find" << endl;
			aigstop = 1;
			return;
		}
		if (CheckBf(bin2_b1, puzknown_perm.bf.u32[0]) < 0) {
			cout << "entry first not seen copie b1 valid to find" << endl;
			aigstop = 1;
			return;
		}
		if (CheckBf(bin_b2yes, puzknown_perm.bf.u32[1]) < 0) {
			cout << "entry first not seen b2 valid to find" << endl;
			aigstop = 1;
			return;
		}
		if (CheckBf(bin2_b2, puzknown_perm.bf.u32[1]) < 0) {
			cout << "entry first not seen copie b2 valid to find" << endl;
			aigstop = 1;
			return;
		}
#endif

	}
	uint32_t loopb2 = 0;
	uint32_t activerb1, activerb2;
	while (++loopb2 <= EXLNLOOP2) {
#ifdef DEBUGKNOWN
		//int CheckBf(BINDEXN & binw, uint32_t bfw) 
		if (CheckBf(bin2_b1, puzknown_perm.bf.u32[0]) < 0) {
			cout << "not seen b1 valid to find" << endl;
			aigstop = 1;
			return;
		}
		if (CheckBf(bin2_b2, puzknown_perm.bf.u32[1]) < 0) {
			cout << "not seen b2 valid to find" << endl;
			aigstop = 1;
			return;
		}
#endif

#ifdef DEBUGEXL
		if (aigstop)cout << "seen aigstop=1" << endl;
		cout << Char2Xout(activeloop) << "==== loop first " << loopb2 << endl;
#endif
		if (aigstop) return;
		minratio = 1000;
		uint64_t ir = FindSockets(activeloop, 2);
		if (ir)  ExtractMin(activeloop, bin2_b1, bin2_b2); 
		if (!ir || minratio > EXLRATIO) {
			ir = FindSockets(activeloop, 3);
			if (ir)ExtractMin(activeloop, bin2_b1, bin2_b2);
			if (!ir || minratio > EXLRATIO) {
				ir = FindSockets(activeloop, 4);
				if (ir)ExtractMin(activeloop, bin2_b1, bin2_b2);
			}
		}
#ifdef DEBUGEXL
		cout << "final selection minratio=" << minratio << endl;
#endif
		if (minratio >= EXLRATIO) break;
		else {
#ifdef DEBUGEXL
			extlw.Debug();
#endif
			if (extlw.mode == 1) {// this is a band1 X band2 Y
				ExtSplitX(bin2_b1, bin2_b1yes, extlw.bfx, activerb1);

#ifdef DEBUGKNOWN
				if (!kpfilt[1]) {// forget if the first is there
					if (extlw.bfx&puzknown_perm.bf.u32[0]) {//hit
						kpfilt[1] = 1;// this must be in this ip
#ifdef TEST_ON
						cout << " loopb1-2 mode 1 to use for known=" << loopb2 << endl;
						cout << "\t\t" << Char27out(extlw.bfx) << " extlw.bfx b1" << endl;
#endif
						Go3(bin2_b1yes, bin2_b2);
				}
			}
#else
				Go3(bin2_b1yes, bin2_b2);
#endif
				if (aigstop) return;

				ExtSplitY(bin2_b2, extlw.tbfy, extlw.ntbfy, activerb2);
			}
			else {// this is a band2 X band1 Y
				ExtSplitX(bin2_b2, bin2_b2yes, extlw.bfx, activerb2);

#ifdef DEBUGKNOWN
				if (!kpfilt[1]) {// forget if the first is there
					if (extlw.bfx&puzknown_perm.bf.u32[1]) {//hit
						kpfilt[1] = 1;// this must be in this ip
#ifdef TEST_ON
						cout << " loopb1-2 mode 2 to use for known=" << loopb2 << endl;
						cout << "\t\t" << Char27out(extlw.bfx) << " extlw.bfx b2" << endl;
#endif
						Go3(bin2_b1, bin2_b2yes);
					}
				}
#else
				Go3(bin2_b1, bin2_b2yes);
#endif
				ExtSplitY(bin2_b1, extlw.tbfy, extlw.ntbfy, activerb1);
			}
		}
		activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
		if (extlw.noxyes < EXLBREAK) break;
	}
	Go3(bin2_b1, bin2_b2);// last call

}

//__________________ process a subset of valid {band1;band2}

/* the subset is made of 2 pieces of the expansion	
	usually the smaller piece is in band 1 

	here 2 loops, outer band1 inner band 2
	for each loop, the process in cut in pieces 
	  depending on the expansion index (2 common cells)
      valid bands are split by size
	  a pair of 2 cells index gives a step

	Uas bands 12 use is the following
		start : lock 128 first uas build a vector for each valid band.
				remaining uas do nothing (UARs)
				the first 128 can be immediately built for each valid band
		band 1 Shrink the UARs table and reshape it (no subset for a ua)
		       use to reshape the Or status of the band 1 piece
		band1 shrink again, but don't clear subsets supersets (too expensive)
		main loop : filter "and of valid band1 valid band2" for the first 128
		Clean all 
		   having the Or status of the clean entries, skrink again the UARs table
		   then for each entry check the residual table 



	 UAs and GUAs tables are reduced to reach a "step size"
	 common cells and possible cells of the step are identified to optimize the process

  
*/

ZS128 * GCHK::G3_Split(int mode, int kill7,	BINDEXN & binw,
	uint32_t ibi2,	INDEXB & indb, ZS128 * p128) {
	//___ sort  the lot per size of valid bands
	uint32_t shift = (mode == 1) ? 0 : 32;
	uint32_t buffer[100000],
		nv[5] = { 0,0,0,0,0 }, // count per size
		// see define ZSTx
		stv[5] = { 0,200,1500,7000,20000 },// start in buffer per size
		*psv[5]; // pointers to starts
	for (int i = 0; i < 5; i++) psv[i] = &buffer[stv[i]];
	BI2_32 biw = binw.t2[ibi2];
	indb.bf =(uint64_t) biw.bf << shift;
	VALIDB * tvb = binw.tvb;
	uint32_t id = biw.istart, iend = biw.iend;
	for (uint32_t iv = id; iv < iend; iv++) {
		VALIDB & wv = tvb[iv];
		uint32_t nc = wv.nval;
		nc--;
		psv[nc][nv[nc]++] = iv;
	}
	if (kill7)nv[4] = 0;
	indb.ntotvb = 0;
	for (int isize = 4; isize >= 0; isize--) {
		if (nv[isize]) {
			indb.ntotvb += nv[isize];
			indb.ncluesmin = isize;
		}
	}
	//_________ build the final tables in 64 mode in indb
	uint64_t	wand = BIT_SET_2X, wor = 0,bf64= biw.bf;
	bf64 <<= shift;
	uint32_t sumpvb = 0;
	BF128 * tvc = vc128;
	if (mode == 2)tvc += 27;
	//BF128 wvect = v128uas;
	//for (int i = 0; i < 2; i++)wvect &= tvc[biw.tval[i]];

	for (int isize = 0; isize < 5; isize++) {//	3 to 7
		uint32_t *psvw = psv[isize];
		INDEXB::ITEM &itemw = indb.titem[isize];
		itemw.ntvb = nv[isize];
		sumpvb += itemw.ntvb;
		itemw.sum_vb = sumpvb;
		itemw.tvb = p128;
		for (uint32_t j = 0; j < nv[isize]; j++) {
			VALIDB vb32 = tvb[psvw[j]];// source valid band
			uint64_t bf = (uint64_t)vb32.bf << shift; // mode 2X
			wand &= bf;
			wor |= bf;
			ZS128 & zb128 = *p128++;
			zb128.bf = bf;
			//zb128.v = wvect;
			//for (uint32_t itv = 0; itv < vb32.nval; itv++)
				//zb128.v &= tvc[vb32.tval[itv]];
		}
	}
	indb.and_g = wand;
	indb.or_g = wor;
	//indb.vand = SetUpVect128(wvect, wand ^bf64);
	//indb.vor = SetUpVect128(wvect, wor ^bf64);
	return  p128;
}

void GCHK::Go3(BINDEXN & bin1, BINDEXN & bin2) {
	p_cpt2g[1]++;
#ifdef TEST_ON
	cout << "go3\t" << bin1.ntvb << "\t " << bin2.ntvb << "\t p_cpt2g[1]="<< p_cpt2g[1]
		<< " aigstop="<<aigstop<< endl;
	debugtest = 0;
#ifdef TEST_ON_KNOWN
	debugtest = 1;
	//if (p_cpt2g[1] == 2)debugtest = 1;
	//if (p_cpt2g[1] > 2)return;
#endif
#endif
#ifdef DEFPHASE
	if (DEFPHASE == -2)return;
#endif
	if (aigstop) return;
	if ((!bin1.ntvb) || (!bin2.ntvb)) return;
	int no7_b2 = t16_min_clues[myband1.i416] > 4;

#ifdef DEBUGKNOWN
	// fb2 must be in the lot
	kn_ir1 = CheckBf(bin1, puzknown_perm.bf.u32[0]);
#ifdef TEST_ON
	cout << "kn_ir1=" << kn_ir1 << endl;
#endif
	if (kn_ir1 < 0)return;
	kn_ir2 = CheckBf(bin2, puzknown_perm.bf.u32[1]);
#ifdef TEST_ON
	cout << "kn_ir2=" << kn_ir2 << endl;
#endif
	if (kn_ir2 < 0)return;
	cout << "got the right Go3 set" << endl;
	if (DEBUGKNOWN) {
		cout << Char27out(bin1.tvb[kn_ir1].bf) << " bfb1 ir1=" << kn_ir1 << endl;
		cout << "\t\t" << Char27out(bin2.tvb[kn_ir2].bf) << " bfb2 ir2=" << kn_ir2 << endl;
	}
#endif


	//______________________________________ loop b1 main loop

#ifdef TEST_ON
	cout << "start loopb1" << endl;
#endif
	for (uint32_t ib1 = 0; ib1 < bin1.nt2; ib1++) {
		if (aigstop) return;
		p_cpt2g[2]++;
#ifdef DEBUGKNOWN
		BI2_32 wi2 = bin1.t2[ib1];
		if ((uint32_t)kn_ir1 >= wi2.istart && (uint32_t)kn_ir1 <= wi2.iend) {
#ifdef TEST_ON
			cout << Char27out(wi2.bf) << " bf bi2 b1 kn_ir1=" << kn_ir1
				<< " this is the expected bi2 band 1 " << wi2.istart << " " << wi2.iend << endl;
#endif
		}
		else continue;
#endif
		G3_Split(1, 0, bin1, ib1, indb1, zs128b1) ;
#ifdef DEBUGKNOWN
		if (DEBUGKNOWN) {
			cout << "global check ntvb=" << indb1.ntotvb << " expected " << wi2.iend - wi2.istart << endl;
			uint32_t cc1 = _popcnt32(puzknown_perm.bf.u32[0]);
			for (uint32_t i = wi2.istart; i < wi2.iend; i++) {
				register uint32_t R = bin1.tvb[i].bf;
				if (R == puzknown_perm.bf.u32[0])
					cout << Char27out(R) << "i=" << i << " check bin1 ok" << endl;
			}
			cc1 -= 3;
			ZS128 * tb = indb1.titem[cc1].tvb;
			for (uint32_t i = 0; i < indb1.titem[cc1].ntvb; i++) {
				register uint32_t R = (uint32_t)tb[i].bf;
				if (R == puzknown_perm.bf.u32[0])
					cout << Char27out(R) << "i=" << i << " check item b1 ok " << endl;
			}
		}
#endif

		fb1 = indb1.and_g;// including more common cells
		acb1 = indb1.or_g;// all these cells can be active in the step
#ifdef TEST_ON_KNOWN
		cout << Char2Xout(indb1.bf) << "\t b1    p_cpt2g[2]=" << p_cpt2g[2] << endl;
		if (debugtest) {// mans tes
			if (p_cpt2g[2] == 3) {
				//cout << "this is the right band1 step" << endl;
				//cout << "this is the  band1 in bug" << endl;
				debugtest = 2;
			}
			//if (p_cpt2g[2] > 3) return;
		}

#endif
		Apply_B1_V();
		p_cpt2g[50] += ntusb1;
		tguas.ApplyLoopB1();


#ifdef DEBUGL1L2
		cout << " bf ib1=" << ib1;		index_xy_b1.Debug();
		cout << "\tnuas=" << ntusb1;		cout << endl;
#endif

		//______________________________________ loop b2

		for (uint32_t ib2 = 0; ib2 < bin2.nt2; ib2++) {
			if (aigstop) return;
			p_cpt2g[3]++;
#ifdef DEBUGKNOWN
			BI2_32 wi2_2 = bin2.t2[ib2];
			if ((uint32_t)kn_ir2 >= wi2_2.istart && (uint32_t)kn_ir2 <= wi2_2.iend) {
#ifdef TEST_ON
				cout << "\t\t" << Char27out(wi2_2.bf) << " bf bi2 b2 kn_ir2=" << kn_ir2
					<< " this is the expected bi2 band 2 "
					<< wi2_2.istart << " " << wi2_2.iend << endl;
#endif
			}
			else continue;
#endif
			G3_Split(2, no7_b2, bin2, ib2, indb2, zs128b2);
#ifdef DEBUGKNOWN
			if (DEBUGKNOWN) {
				cout << "global check ntvb=" << indb2.ntotvb << " expected " << wi2_2.iend - wi2_2.istart << endl;
				uint32_t cc2 = _popcnt32(puzknown_perm.bf.u32[1]);
				for (uint32_t i = wi2_2.istart; i < wi2_2.iend; i++) {
					register uint32_t R = bin2.tvb[i].bf;
					if (R == puzknown_perm.bf.u32[1])
						cout << Char27out(R) << "i=" << i << " check bin2 ok" << endl;
				}
				cc2 -= 3;
				ZS128 * tb = indb2.titem[cc2].tvb;
				for (uint32_t i = 0; i < indb2.titem[cc2].ntvb; i++) {
					register uint32_t R = (uint32_t)(tb[i].bf >> 32);
					if (R == puzknown_perm.bf.u32[1])
						cout << Char27out(R) << "i=" << i << " check item b2 ok " << endl;
				}
			}
#endif				
			fb2 = indb2.and_g;// including more common cells
			acb2 = indb2.or_g;// all these cells can be active in the step
			fb12 = fb1 | fb2;
			acb12 = acb1 | acb2;
			if (Apply_B2_V())continue;// check dead branch set vector
#ifdef TEST_ON_KNOWN
			if (debugtest > 1) {
				cout << Char2Xout(indb2.bf) << "\t b2    p_cpt2g[3]=" << p_cpt2g[3] << endl;
/*
				if (p_cpt2g[3] == 44) {
					//cout << "this is the right band2 step" << endl;
					cout << "this is the  band2 in bug" << endl;
					debugtest = 3;
				}
				if (p_cpt2g[3] > 44) {
					cout << " stop debug on 45" << endl;
					debugtest = 0;
					//aigstop = 1; return;
				}
*/
			}
#endif
#ifdef DEBUGL1L2
			cout << Char2Xout(fb12) << "\t bf ib2=" << ib2 << " p_cpt2g[3]=" << p_cpt2g[3] << endl;;
			g4t_b2.Debug();
#endif			
			p_cpt2g[5]++;
#ifdef TEST_ON
			p_cpt2g[51] += ntusb2;
			p_cpt2g[52] += ntusb1;
#endif
			Do128uas();
		}
	}
#ifdef TEST_ON
	//cout << "back go3" << endl;
#endif

}


//______________main loop 128 first uas
void GCHK::Do128uas() {//apply first 128 filter
#ifdef TEST_ON_KNOWN
	if (debugtest ==3) {
		cout << "step to debug start the main loop nuas= " << ntusb2 << " cpt5=" << p_cpt2g[5] << endl;
		cout << "index b1"; indb1.Debug(); cout << endl;
		cout << "index b2"; indb2.Debug(); cout << endl;
		cout << Char2Xout(fb12) << "fb12 and " << endl;;
		cout << Char2Xout(acb12) << "acb12 or " << endl;;
	}
#endif


#ifdef DEBUGKNOWN
	uint32_t cc1 = _popcnt32(puzknown_perm.bf.u32[0]),
		cc2 = _popcnt32(puzknown_perm.bf.u32[1]);
	cout << "do128 expected step  ncluesb1="<<cc1 	<<" ncluesb2=" << cc2 << endl;
	cc1 -= 3;// kill unused band1
	for (uint32_t i = 0; i < 5; i++) 	if (i != cc1) indb1.titem[i].ntvb = 0;
	cc2 -= 3;// kill band2 over cc2
	INDEXB::ITEM *t = indb2.titem;
	uint32_t s = t[cc2].sum_vb;
	for (uint32_t i = cc2 + 1; i < 5; i++) {
		t[i].ntvb = 0;
		t[i].sum_vb = s;
	}
	if (DEBUGKNOWN) {
		indb1.Debug();
		indb2.Debug();
	}
#endif

	// reset more uas tables
	moreuas_12_13.Init();	moreuas_14.Init();
	moreuas_15.Init();	moreuas_AB_small.Init();
	moreuas_AB.Init();	moreuas_AB_big.Init();
	zb = indb2.titem[0].tvb;
	if (start_perm) {// clues count band 3 strictly higher
		nza = indb1.titem[0].ntvb;
		if (nza) { // 3 b1 all b2
			za = indb1.titem[0].tvb;
			nzb = indb2.ntotvb;
			Do128Chunk();
		}
		nza = indb1.titem[1].ntvb;
		if (nza) { // 4 b1   b2<7
			za = indb1.titem[1].tvb;
			nzb = indb2.titem[3].sum_vb;
			Do128Chunk();
		}
		nza = indb1.titem[2].ntvb;
		if (nza) { // 5 b1   b2<7
			za = indb1.titem[2].tvb;
			nzb = indb2.titem[3].sum_vb;
			Do128Chunk();
		}
		nza = indb1.titem[3].ntvb;
		if (nza) { // 6 b1   b2<6
			za = indb1.titem[3].tvb;
			nzb = indb2.titem[2].sum_vb;
			Do128Chunk();
		}
		nza = indb1.titem[4].ntvb;
		if (nza) { // 7 b1   
			za = indb1.titem[4].tvb;
			if (start_perm == 2) {//3 b2
				nzb = indb2.titem[0].ntvb;
				Do128Chunk();
			}
			else {//	3;4 b2
				nzb = indb2.titem[1].sum_vb;
				Do128Chunk();
			}
		}
	}
	else { //higher or equal for band 3 add 747 477 666  (666 worst case)
		nza = indb1.titem[1].sum_vb;
		if (nza) { // 3 4 b1 all b2
			za = indb1.titem[0].tvb;
			nzb = indb2.ntotvb;
			Do128Chunk();
		}
		nza = indb1.titem[2].ntvb;
		if (nza) { // 5 b1   b2<7
			za = indb1.titem[2].tvb;
			nzb = indb2.titem[3].sum_vb;
			Do128Chunk();
		}
		nza = indb1.titem[3].ntvb;
		if (nza) { // 6 b1   b2<7
			za = indb1.titem[3].tvb;
			nzb = indb2.titem[3].sum_vb;
			Do128Chunk();
		}
		nza = indb1.titem[4].ntvb;
		if (nza) { // 7 b1   
			za = indb1.titem[4].tvb;
			nzb = indb2.titem[1].sum_vb;
			Do128Chunk();
		}
	}
}


void GCHK::Do128Chunk() {
	if (!nzb) return;
	if (aigstop) return;
	p_cpt2g[54]+= (nza * nzb);
	n_to_clean = 0;
#ifdef DEBUGKNOWN
	if (DEBUGKNOWN) {
		cout << "chunk " << nza << " x " << nzb << endl;
		for (uint32_t i = 0; i < nza; i++) {
			register uint32_t R = (uint32_t)za[i].bf;
			if (R == puzknown_perm.bf.u32[0])
				cout << Char2Xout(za[i].bf) << " i=" << i << " za ok" << endl;
		}
		for (uint32_t i = 0; i < nzb; i++) {
			register uint32_t R = (uint32_t)(zb[i].bf >> 32);
			if (R == puzknown_perm.bf.u32[1])
				cout << Char2Xout(zb[i].bf) << " i=" << i << " zb ok" << endl;
		}
	}
#endif
	if ((nza * nzb) < 5000) {
		if (nza < nzb)Do128Go(za, zb, nza, nzb);
		else Do128Go(zb, za, nzb, nza);
		if (n_to_clean) CleanAll();
		return;
	}
	// cut in chunks max Xchunk Ychunk
	ZS128 * z1 = za, *z2 = zb;
	uint32_t n1 = nza, n2 = nzb;
	uint32_t  ideb2 = 0, iend2 = YCHUNK64;
	if (iend2 > n2)iend2 = n2;
	while (ideb2 < n2) { //Y chunk
		uint32_t ny = iend2 - ideb2;
		uint32_t ideb1 = 0, iend1 = XCHUNK64;
		if (iend1 > n1)iend1 = n1;

		while (ideb1 < n1) {// X chunk
			uint32_t nx = iend1 - ideb1;
			if (nx < ny)Do128Go(&z1[ideb1], &z2[ideb2],  nx, ny);
			else Do128Go( &z2[ideb2], &z1[ideb1], ny, nx);
			ideb1 = iend1; iend1 += XCHUNK64;
			if (iend1 > n1)iend1 = n1;
		}
		ideb2 = iend2; iend2 += YCHUNK64;
		if (iend2 > n2)iend2 = n2;
	}
	if (n_to_clean) CleanAll();
}

void GCHK::Do128Go(ZS128 * a, ZS128 * b, uint32_t na, uint32_t nb) {
	if (aigstop)return;
#ifdef DEBUGKNOWN
	uint32_t aig = 2,i1,i2;
	register uint64_t R1 = puzknown_perm.bf.u32[0], R2 = puzknown_perm.bf.u32[1];
	R2 <<= 32;
	for (uint32_t i = 0; i < na; i++) {
		register uint64_t R =a[i].bf;
		if (R == R1 || R==R2) {
			//cout << Char2Xout(R) << " i=" << i << " a ok" << endl;
			aig = 1; i1 = i; break;
		}
	}
	if (aig == 2) return;
	for (uint32_t i = 0; i < nb; i++) {
		register uint64_t R = b[i].bf;
		if (R == R1 || R == R2) {
			//cout << Char2Xout(R) << " i=" << i << " b ok" << endl;
			aig = 0; i2 = i; break;
		}
	}
	if (aig)return;
	cout << "chunk go " << na << " x " << nb << endl;
	if (DEBUGKNOWN) {
		if ((a[i1].v& b[i2].v).isNotEmpty()) {
			cout << "bug vectors" << endl;
			a[i1].v.Print("a vect");
			b[i2].v.Print("b vect");
			(a[i1].v& b[i2].v).Print("a&b vect");
			aigstop = 1;
			return;
		}
	}
#endif
	//check a matrix band 1 band2 for potential 2 bands valid 11 clues
	register ZS128 * Ra = &a[na - 1];
	register uint64_t * Rs = &to_clean[n_to_clean];
	for (; Ra >= a; Ra--) {
		register ZS128 * Rb = &b[nb - 1];
		BF128 va = (*Ra).v;
		register uint64_t bfa = (*Ra).bf;
		for (; Rb >= b; Rb--)
			if ((Rb->v&va).isEmpty()) 		*Rs++ = bfa | Rb->bf;
	}
	n_to_clean = Rs - to_clean;
	if (n_to_clean > 10000)CleanAll();
}


void BuildBaseAndCellsVector(uint32_t nuas, BF128 & bv, BF128 * cellsv, uint64_t * tu) {
	bv = maskLSB[nuas];// Uas vector
	memset(cellsv, 255, sizeof gchk.vc64_192);// all bits to 1
	uint32_t cc64;// build cells vectors A
	for (uint32_t i = 0; i < nuas; i++) {
		register uint64_t Rw = tu[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			cellsv[From_128_To_81[cc64]].clearBit(i);
		}
	}
}

void AddUaToVector(uint64_t ua12, BF128 * cellsv, uint32_t iloc) {
	register uint64_t Rw = ua12 & BIT_SET_2X;
	uint32_t cc64;// build cells vectors A
	while (bitscanforward64(cc64, Rw)) {// look for  possible cells
		Rw ^= (uint64_t)1 << cc64;// clear bit
		cellsv[From_128_To_81[cc64]].clearBit(iloc);
	}
}

inline void SetUpBaseVector(uint32_t ncl, uint32_t *tcl,
	BF128 & bv, BF128 * cellsv, BF128 & bvf) {
	bvf = bv;// start with valid uas
	for (uint32_t i = 0; i < ncl; i++) {//all cells common to the step
		bvf &= cellsv[tcl[i]];
	}
}
inline int SetUpStepV(uint32_t * tc, uint32_t ntc, BF128 & vb, BF128 & vd, BF128 * tvc) {
	vd = vb;
	if (vb.isEmpty())return 1;
	for (uint32_t i = 0; i < ntc; i++)vd &= tvc[tc[i]];
	return vd.isEmpty();
}


void GCHK::Apply_B1_V() {
	ntusb1 = ntusb1_128 = 0;
	register uint64_t * tua = genuasb12.tua;
	register uint32_t nua = genuasb12.nua;
	register uint64_t filter = fb1, Ra = acb1 | BIT_SET_B2;
	g4t_b1.Shrink(g4t_start, filter);
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint64_t Ru = tua[iua];
		if (Ru&filter) continue;
		Ru &= Ra;
		if (ntusb1_128 < 128) {
			tusb1_128[ntusb1_128++] = Ru;
		}
		else tusb1[ntusb1++] = Ru;
	}
	// Vector for the first 128
	v128uas = maskLSB[ntusb1_128];// Uas vector
	memset(vc128, 255, sizeof vc128);// all bits to 1
	uint32_t cc64;// build cells vectors A
	for (uint32_t i = 0; i < ntusb1_128; i++) {
		register uint64_t Rw = tusb1_128[i];
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc128[From_128_To_81[cc64]].clearBit(i);
		}
	}
	// Apply vect to b1 step bi2
	ZS128 * t128 = indb1.titem[0].tvb;
	for (uint32_t idetb1 = 0; idetb1 < indb1.ntotvb; idetb1++) {
		ZS128 & w128 = t128[idetb1];
		BF128 wvect = v128uas;// only uas not hit by fb1
		register uint64_t W = w128.bf^fb1;
		while (bitscanforward64(cc64, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc64;// clear bit
			wvect &=vc128[From_128_To_81[cc64]];
		}
		w128.v = wvect;
	}
	if (debugtest==2) {
		cout << "debug end Apply_B1_V ntusb1="<<ntusb1 << endl;
		//for (int i = 0; i < 128; i++) {
			//cout << Char2Xout(tusb1_128[i]) << " i128=" << i << endl;
		//}
	}
#ifdef DEBUGKNOWN
	cout << "debug end Apply_B1_V ntusb1="<<ntusb1 << endl;
	//for (int i = 0; i < 128; i++) {
		//cout << Char2Xout(tusb1_128[i]) << " i128=" << i << endl;
	//}
#endif


}
int GCHK::Apply_B2_V() {
	register uint64_t Ra = acb12, filter = fb2;
	register uint64_t * tua = tusb1;
	g4t_b2.Shrink(g4t_b1, filter);
	ntusb2 = 0;
	for (uint32_t iua = 0; iua < ntusb1; iua++) {
		register uint64_t Ru = tua[iua];
		if (!(Ru&filter)) {
			Ru &= Ra;
			if (!Ru) { return 1; }
			tusb2[ntusb2++] = Ru;
		}
	}
	uint32_t cc64;// build cells vectors A
	BF128 wvectc = v128uas;// setup common cells 
	register uint64_t F = fb2;
	while (bitscanforward64(cc64, F)) {// look for  possible cells
		F ^= (uint64_t)1 << cc64;// clear bit
		wvectc &= vc128[From_128_To_81[cc64]];
	}
	// Apply vect to b1 step bi2
	ZS128 * t128 = indb2.titem[0].tvb;
	for (uint32_t idetb2 = 0; idetb2 < indb2.ntotvb; idetb2++) {
		ZS128 & w128 = t128[idetb2];
		BF128 wvect = wvectc;// only uas not hit by fb1
		register uint64_t W = w128.bf^fb2;
		while (bitscanforward64(cc64, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc64;// clear bit
			wvect &= vc128[From_128_To_81[cc64]];
		}
		w128.v = wvect;
	}

	//_________ setup the brute force start
	nclues_step = 0;
	uint64_t w = fb12;
	uint32_t xcell;
	while (bitscanforward64(xcell, w)) {
		w ^= (uint64_t)1 << xcell;
		tclues[nclues_step++] = From_128_To_81[xcell];
	}
	tcluesxy = &tclues[nclues_step];
	zh2b_i1.ValidXY_Step(tclues, nclues_step);

	return 0;
}

void TGUAS::ApplyLoopB1() {
	register uint64_t Bf = gchk.fb1;
	register uint64_t Ac = gchk.acb1 | BIT_SET_B2;

	for (uint32_t i = 0; i < 36; i++) {
		GUA & wg = tgua_start[i],&wgd= tgua_b1[i];
		wgd.Init(wg.i36,0);
		for (uint32_t j = 0; j < wg.nua; j++) {// apply new subsets
			register uint64_t Ua = wg.tua[j];
			if (Ua & Bf) continue;
			Ua &= Ac;// could be empty
			uint64_t cc = _popcnt64(Ua);
			wgd.Adduacheck(Ua | (cc << 59)); // no redundancy
		}
	}

}
void TGUAS::ApplyLoopB2() {
	// relay table per size
	uint32_t ntt[38]; // count for tt
	BF128 tt[38][1000];// guas 54 ua12 plus code 0_161 for the gangster
	memset(ntt, 0, sizeof ntt);
	register uint64_t Bf = gchk.fb12c;
	register uint64_t Ac = gchk.acb12c;

	for (uint32_t i36 = 0; i36 < 36; i36++) {
		GUA & wg = tgua_b1[i36];
		guaw.Init(wg.i36, 0);
		for (uint32_t j = 0; j < wg.nua; j++) {// apply new subsets
			register uint64_t Ua = wg.tua[j];
			if (Ua & Bf) continue;
			Ua &= Ac;// could be empty
			if (!Ua) {
				guaw.tua[0] = 0;
				guaw.nua = 1;
				break;
			}
			guaw.tua[guaw.nua++]=Ua;
		}
		if (guaw.nua) {
			// split per size in 54 + index 0-26 or 0-7
			BF128 w;
			w.bf.u64[1] = wg.i36;
			for (uint32_t j = 0; j < guaw.nua; j++) {
				register uint64_t Ua = guaw.tua[j],
					cc = _popcnt64(Ua);
				if (cc > 18) continue; // should not be
				if (i36 > 26) {// this is a g3
					cc += 19;
				}
				w.bf.u64[0] = (Ua & BIT_SET_27) | ((Ua & BIT_SET_B2) >> 5);
				tt[cc][ntt[cc]++] = w;
			}
		}
	}
	//_____________ create vectors
	nvg2 = 0;
	for (uint32_t i = 0; i < 19; i++) {// all guas2
		uint32_t n = ntt[i];
		BF128 * tw = tt[i];
		for (uint32_t j = 0; j < n; j++) {
			BF128 w = tw[j];
			AddVG2(w.bf.u64[0], w.bf.u32[2]);
			if (nvg2 >= 512) break;
		}
		if (nvg2 >= 512) break;
	}
	nvg3 = 0;
	for (uint32_t i = 19; i < 38; i++) {// all guas2
		uint32_t n = ntt[i];
		BF128 * tw = tt[i];
		for (uint32_t j = 0; j < n; j++) {
			BF128 w = tw[j];
			AddVG3(w.bf.u64[0], w.bf.u32[2]);
			if (nvg3 >= 256) break;
		}
		if (nvg3 >= 256) break;
	}

}

//__________________ end of phase 1 process the file and clean
inline int Check128uas(uint32_t ncl, uint32_t *tcl,	BF128  bv, BF128 * cellsv) {
	for (uint32_t i = 0; i < ncl; i++) {//all cells common to the step
		bv &= cellsv[tcl[i]];
	}
	return bv.isNotEmpty();
}

int GCHK::Is_B12_Not_Unique() {
	myua = zh2b[0].Valid_XY(tcluesxy, nclues);
	if (myua) {//not unique do the best with the UA
		uint64_t cc64 = _popcnt64(myua&BIT_SET_2X);
		p_cpt2g[9]++;
		if (cc64 < 12) {// this should never be check for a bug
			cerr << endl << endl << Char2Xout(myua) << " ua < 12 to add   clean" << endl;
			cout << endl << endl << Char2Xout(myua) << " ua < 12 to add   clean" << endl;
			aigstop = 1;
			return 1;
		}

		if (cc64 < 18) {
			if(cc64 < 14)moreuas_12_13.Add(myua);
			else if (cc64 == 14)moreuas_14.Add(myua);
			else if (cc64 == 15)moreuas_15.Add(myua);
			else moreuas_AB_small.Add(myua);
			register uint64_t ua_add = myua | (cc64 << 59);
			genuasb12.AddUACheck(ua_add);
			tusb1[ntusb1++] = myua;
		}
		else if (cc64 < 21)			moreuas_AB.Add(myua);
		else moreuas_AB_big.Add(myua);
		return 1;
	}
	return 0;
}


int GCHK::Clean_1() {// check filter other uas 12
	register uint64_t And = BIT_SET_2X, Or = 0;
	for (uint64_t i = 0; i < n_to_clean; i++) {// setup and / or for this set
		register uint64_t bf = to_clean[i];
		And &= bf; Or |= bf;
	}
	{ //  collect still valid uas and check dead branch (skip 0_63 checked)
		ntusr = 0;
		for (uint32_t iua = 64; iua < ntusb2; iua++) {
			register uint64_t Ru = tusb2[iua];
			if (!(Ru&And)) {
				Ru &= Or;
				if (!Ru) 	return 1;
				else tusr[ntusr++] = Ru;
			}
		}
	}
	while (1) {// loop if > 128 uas in tusb2
		if (!ntusr) {	n_to_clean2 = n_to_clean;	return 0;	}
		n_to_clean2 = 0;
		if (ntusr <= 10 || n_to_clean < 10) return Clean_1_all();
		if (ntusr <= 64) return Clean_1_64(And);
		if (ntusr <= 128) return Clean_1_128(And, ntusr);
		// more than 128, do the first 128
		if (Clean_1_128(And, 128)) return 1; // at least one ua not hit
		n_to_clean = n_to_clean2; // new status for clean (>0)
		And = BIT_SET_2X; Or = 0;// redo and/or
		for (uint64_t i = 0; i < n_to_clean; i++) {
			register uint64_t bf = to_clean[i];
			And &= bf; Or |= bf;
		}
		for (uint32_t i = 128; i < ntusr; i++)
			tusr[i - 128] = tusr[i] & Or;
		ntusr -= 128;
	}
}
int GCHK::Clean_1_all() {// <= 10 remaining uas
	for (uint64_t i = 0; i < n_to_clean; i++) {// direct check
		register uint64_t Bf = to_clean[i];
		int aig = 1;
		for (uint32_t iua = 0; iua < ntusr; iua++)
			if (!(Bf & tusr[iua])) { aig = 0; break; }
		if (aig)  to_clean[n_to_clean2++] = Bf;
	}
	return (n_to_clean2 == 0);
}
int GCHK::Clean_1_64(uint64_t bfc) {// <= 64 uas
	uint64_t ve = maskLSB[ntusr].u64[0], vc[54];
	memset(vc, 255, sizeof vc);// all bits to 1
	uint32_t cc64;// build cells vectors A
	uint64_t biti = 1;
	for (uint32_t i = 0; i < ntusr; i++, biti <<= 1) {
		register uint64_t Rw = tusr[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc[From_128_To_81[cc64]] ^= biti;
		}
	}
	{// apply common cells to ve
		register uint64_t Bf = bfc;
		while (bitscanforward64(cc64, Bf)) {
			Bf ^= (uint64_t)1 << cc64;
			ve &= vc[From_128_To_81[cc64]];
		}
	}
	for (uint64_t i = 0; i < n_to_clean; i++) {
		register uint64_t Bf = to_clean[i] ^ bfc,
			v = ve;// initial vector
		while (bitscanforward64(cc64, Bf)) {
			Bf ^= (uint64_t)1 << cc64;
			v &= vc[From_128_To_81[cc64]];
		}
		if (!v) to_clean[n_to_clean2++] = to_clean[i];
	}
	return (n_to_clean2 == 0);
}
int GCHK::Clean_1_128(uint64_t bfc, uint32_t nu) {
	BF128 ve, vc[54];
	ve = maskLSB[nu];// Uas vector
	memset(vc, 255, sizeof vc);// all bits to 1
	uint32_t cc64;
	for (uint32_t i = 0; i < nu; i++) {
		register uint64_t Rw = tusr[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc[From_128_To_81[cc64]].clearBit(i);
		}
	}

	{// apply common cells to ve
		register uint64_t Bf = bfc;
		while (bitscanforward64(cc64, Bf)) {
			Bf ^= (uint64_t)1 << cc64;
			ve &= vc[From_128_To_81[cc64]];
		}
	}
	for (uint64_t i = 0; i < n_to_clean; i++) {
		register uint64_t Bf = to_clean[i] ^ bfc;
		BF128 	v = ve;// initial vector
		while (bitscanforward64(cc64, Bf)) {
			Bf ^= (uint64_t)1 << cc64;
			v &= vc[From_128_To_81[cc64]];
		}
		if (v.isEmpty()) to_clean[n_to_clean2++] = to_clean[i];
	}
	return (n_to_clean2 == 0);
}

int GCHK::CleanAllFifo() {
	n_to_clean2 = 0;
	for (uint64_t i = 0; i < n_to_clean; i++) {
		register uint64_t bf = to_clean[i];
		if (moreuas_12_13.Check(bf))continue;
		if (moreuas_14.Check(bf))continue;
		if (moreuas_15.Check(bf))continue;
		if (moreuas_AB_small.Check(bf))continue;
		if (moreuas_AB.Check(bf)) continue;
		if (moreuas_AB_big.Check(bf)) continue;
		to_clean[n_to_clean2++] = bf;
	}
	return (n_to_clean2 == 0);
}


void GCHK::Clean_2() {
	register uint64_t And = BIT_SET_2X, Or = 0; // re do and or
	for (uint64_t i = 0; i < n_to_clean2; i++) {// setup and / or for this set
		register uint64_t bf = to_clean[i];
		And &= bf; Or |= bf;
	}
	fb12c = And;	acb12c = Or;
	tguas.ApplyLoopB2();// create vectors sockets 2 3
	g4t_clean.Shrink(g4t_b2, fb12c);
	{	//_________ setup the brute force start
		nclues_step = 0;
		register uint32_t xcell;
		while (bitscanforward64(xcell, And)) {
			And ^= (uint64_t)1 << xcell;
			uint32_t cell = From_128_To_81[xcell];
			tclues[nclues_step++] = cell;
		}
		tcluesxy = &tclues[nclues_step];
	}

	clean_valid_done = 0;



#ifdef TEST_ON_KNOWN
	if (debugtest == 4) {

		cout << " clean2 status before loop nvg2=" << tguas.nvg2 << "  nvg3=" << tguas.nvg3
			<< " g4t_clean.ntua4=" << g4t_clean.ntua4 <<" n_to_clean2=" << n_to_clean2 << endl;
		cout << Char2Xout(fb12c) << " and" << endl;
		cout << Char2Xout(Or) << " or" << endl;
		//for (uint64_t i = 0; i < n_to_clean2; i++)
		//	cout << Char2Xout(to_clean[i]) << endl;
	}
#endif


	// chek the minimum clues in guas2 guas3
	n_to_clean2 = 0;
	moreuasy.Init();

#ifdef DEFPHASE
	if (DEFPHASE == -4) return; 
#endif

	for (uint64_t iclean = 0; iclean < n_to_clean; iclean++) {// now
		if (aigstop) {
			cout << "aigstop=1 wbf=" << wb12bf << endl;
			return;
		}
		wb12bf = to_clean[iclean];
		if (moreuasy.Check(wb12bf))continue;
		Clean_3();
	}
}

void GCHK::Clean_3() {// clean 2 for a given band 1+2
	ncluesb3 = 18 - (uint32_t)_popcnt64(wb12bf);
#ifdef TEST_ON
	int locdiag = 0;
	//if (wb12bf == debugvalbf)locdiag = 1;
	if (locdiag) cout << Char2Xout(wb12bf) << "start clean2 on expected bf ncluesb3 ="
		<< ncluesb3 << endl;
#endif

	nclues = 0;
	{
		register uint64_t w = wb12bf ^ fb12c;
		register uint32_t xcell, cell;
		while (bitscanforward64(xcell, w)) {
			w ^= (uint64_t)1 << xcell;
			cell = From_128_To_81[xcell];
			tcluesxy[nclues++] = cell;
		}
	}

	aigstopxy = 0;
	g2ok = g3ok = g2moreok = nuasb3_1 = 0;
	memset(&smin, 0, sizeof smin);
	Clean_3G2();						//______ apply gua2s
	sminr = smin;
	smin.SetMincount();
#ifdef TEST_ON 
	if (locdiag)smin.Status(" status g2");	
#endif
	if (smin.mincount > ncluesb3)return;
	p_cpt2g[58]++;
	smin = sminr;
	Clean_3G3();						// add gua3s
	sminr = smin;
	smin.SetMincount();
#ifdef TEST_ON 
	if (locdiag) smin.Status(" status g3");	
#endif				
	if (smin.mincount > ncluesb3)return;
	p_cpt2g[59]++;
	nmiss = ncluesb3 - smin.mincount; 

	if (!clean_valid_done) {
		clean_valid_done = 1;
		myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
		if (myua) {
			NewUaB12();		aigstopxy = 1; return;
		}
	}

	Clean_3BuildIfOf();					//_Build  in/out field 

#ifdef DEBUGKNOWN
#ifdef TEST_ON
	cout << "nmiss=" << nmiss << " nif=" << nuasb3_1 << " nof=" << nuasb3_2 << endl;
	smin.Status("final status for call b3");
	Debugifof();
#endif
#endif

#ifdef TEST_ON 
	if (locdiag) {
		cout << p_cpt2g[58] << "nmiss=" << nmiss << " nif=" << nuasb3_1 << " nof=" << nuasb3_2 << endl;
		smin.Status("final status for call b3");
	}
#endif

#ifdef DEFPHASE
	if (DEFPHASE == -5) continue;
#endif

						//_______________________ start band 3 process 
	moreuas_b3.Init();
	ua_out_seen=0;
	hh0.diagh = 0;
	uint32_t min = 100, uamin;

	if (!nmiss) {
		if (nuasb3_2) return;
		zhou[0].PartialInitSearch17(tclues, nclues + nclues_step);
		hh0.GoMiss0();
		p_cpt2g[10]++;
		p_cpt2g[57]++;
		return;
	}
	if (nmiss == 1) {
		register uint32_t andout = BIT_SET_27;
		for (uint32_t i = 0; i < nuasb3_2; i++) {
			register uint32_t	Ru = uasb3_2[i];
			andout &= Ru;
			if (!andout) return;
		}
		andmiss1 = andout;
		zhou[0].PartialInitSearch17(tclues, nclues + nclues_step);
		p_cpt2g[11]++;
		p_cpt2g[57]++;
		hh0.GoMiss1();
		return;
	}
	if (nmiss == 2) {//skip if not 2 out needed
		/*  bug to see here
		if (!nuasb3_2) {// subcritical in hn if solved
			zhou[0].PartialInitSearch17(tclues, nclues + nclues_step);
			if (zhou[1].CallMultipleB3(zhou[0], smin.critbf, 0)) {
				// setup the ua found
				uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
				moreuas_b3.Add(ua);
				NewUaB3();
				uasb3_2[nuasb3_2++] = ua;
			}
		}
		*/
		if (nuasb3_2>1) {
			register uint32_t andout = BIT_SET_27;
			for (uint32_t i = 0; i < nuasb3_2; i++) {
				register uint32_t	Ru = uasb3_2[i];
				andout &= Ru;
				uint32_t cc = _popcnt32(Ru);
				if (cc < min) { min = cc; uamin = Ru; }
			}
			if (!andout) {
				nmiss2ok = 0;
				register uint32_t Ru = uamin, bit, cc, andx = 0, n = 0;
				while (bitscanforward(cc, Ru)) {// look for  possible cells
					bit = 1 << cc;
					Ru ^= bit;// clear bit
					andx = BIT_SET_27 ^ smin.critbf;
					andx &= ~bit;// be sure to have it cleared
					for (uint32_t i = 1; i < nuasb3_2; i++) {
						register uint32_t Ru2 = uasb3_2[i];
						if (!(Ru2 & bit)) {
							andx &= Ru2; n = 1;
						}
					}
					if (andx) {
						uint32_t * t = miss2ok[nmiss2ok++];
						t[0] = cc;		t[1] = andx;	t[2] = n;
					}
				}
				if (nmiss2ok) {// if sminplus is ncluess must have an out field 2 clues valid -
					zhou[0].PartialInitSearch17(tclues, nclues + nclues_step);
					p_cpt2g[12]++;
					p_cpt2g[57]++;
					hh0.GoMiss2(uamin);
				}
				return;
			}
		}
	}
	//  go direct expand nmss > 2 or no out field
	for (uint32_t i = 0; i < nuasb3_2; i++) {
		uasb3_1[nuasb3_1++] = uasb3_2[i];
	}
	//  sort all on size
	uint32_t ntt[27]; // count for tt
	uint32_t tt[27][100];// guas 54 ua12 plus code 0_161 for the gangster
	memset(ntt, 0, sizeof ntt);
	for (uint32_t i = 0; i < nuasb3_1; i++) {
		register uint32_t Ua = uasb3_1[i] & BIT_SET_27,
			cc = _popcnt32(Ua);
		tt[cc][ntt[cc]++] = Ua;
	}
	nuasb3_1 = 0;
	for (uint32_t cc = 0; cc < 27; cc++) if (ntt[cc])
		for (uint32_t j = 0; j < ntt[cc]; j++)
			uasb3_1[nuasb3_1++] = tt[cc][j];
#ifdef TEST_ON
	if (locdiag) Debugifof();

#endif
	zhou[0].PartialInitSearch17(tclues, nclues + nclues_step);

	p_cpt2g[13]++;
	p_cpt2g[57]++;
	ExpandB3();
}

void GCHK::Clean_3G2() {
	nb64_1 = (tguas.nvg2 + 63) >> 6;
	for (uint32_t iv = 0; iv < nb64_1; iv++) {
		TVG64 &vv = tvg64g2[iv];
		register uint64_t V = vv.v;
		for (int j = 0; j < nclues; j++)
			V &= vv.cells[tcluesxy[j]];
		register uint32_t cc64;
		while (bitscanforward64(cc64, V)) {
			V ^= (uint64_t)1 << cc64;// clear bit
			uint32_t i27 = vv.ti162[cc64], bit27 = 1 << i27;
			if (g2ok&bit27) continue;
			g2ok |= bit27;
			register uint32_t bit = 1 << (i27 / 3);
			smin.mini_bf3 |= smin.mini_bf2&bit;
			smin.mini_bf2 |= smin.mini_bf1&bit;
			smin.mini_bf1 |= bit;
			smin.critbf |= myband3.ua27_bf[i27];
			smin.pairs27 |= myband3.ua2pair27[i27];
			uasb3_1[nuasb3_1++] = guapats2[i27];
		}
	}
}
void GCHK::Clean_3G3() {
	nb64_2 = (tguas.nvg3 + 63) >> 6;
	for (uint32_t iv = 0; iv < nb64_2; iv++) {
		TVG64 &vv = tvg64g3[iv];
		register uint64_t V = vv.v;
		for (int j = 0; j < nclues; j++)
			V &= vv.cells[tcluesxy[j]];
		register uint32_t cc64;
		while (bitscanforward64(cc64, V)) {
			V ^= (uint64_t)1 << cc64;// clear bit
			uint32_t i9 = vv.ti162[cc64] - 27, bit9 = 1 << i9;
			if (g3ok&bit9) continue;
			g3ok |= bit9;
			smin.mini_triplet |= bit9;
			uasb3_1[nuasb3_1++] = guapats3[i9];
		}
	}

}
void GCHK::Clean_3BuildIfOf() {
	//_________Build now tables in/out field 
	nuasb3_2 = 0;
	{
		register int  Rfilt = smin.critbf;
		register uint64_t Rhit = wb12bf;
		for (uint32_t i = 0; i < g4t_clean.ntua4; i++) {
			register uint32_t	Ru = g4t_clean.tua4[i].bf.u32[2];
			if (g4t_clean.tua4[i].bf.u64[0] & Rhit) continue;
			if (Ru & Rfilt)	uasb3_1[gchk.nuasb3_1++] = Ru;
			else uasb3_2[nuasb3_2++] = Ru;
		}
		for (uint32_t i = 0; i < myband3.nua; i++) {
			register uint32_t Ru = myband3.tua[i] & BIT_SET_27;
			if (Ru & Rfilt)	uasb3_1[nuasb3_1++] = Ru;
			else uasb3_2[nuasb3_2++] = Ru;
		}
	}

}


void GCHK::CleanAll() {
	if (!n_to_clean) return;
	p_cpt2g[55] += n_to_clean;
	p_cpt2g[6]++;

#ifdef DEBUGKNOWN
	cout << " entry clean " << n_to_clean << " p_cpt2g[6]=" << p_cpt2g[6] << endl;
	if (kpfilt[2]) {
		n_to_clean = 0;
		return;
	}
	for (uint64_t i = 0; i < n_to_clean; i++) {
		uint64_t bf = to_clean[i];
		if (bf == puzknown_perm.bf.u64[0]) {
			cout << Char2Xout(bf)
				<< "\t\tclean all seen bf i=" << i << " forced to one clean p_cpt2g[10]=" << p_cpt2g[10] 
				<< "value="<< bf << endl;
			cout << Char2Xout(fb1|fb2) << " step applied" << endl;
			to_clean[0] = bf;
			n_to_clean = 1;
			kpfilt[2] = 1;
			break;
		}
	}
	if (!kpfilt[2]) {// skip if not seen
		n_to_clean = 0;
		return;
	}

#endif
	debugvalbf = 225195374533435462;

#ifdef TEST_ON_KNOWN
	if (debugtest >= 3) {
		debugtest = 3;
		for (uint64_t i = 0; i < n_to_clean; i++) {
			uint64_t bf = to_clean[i];
			if (bf == debugvalbf) {
				cout << "seen expected bf in clean" << endl;
				debugtest = 4;
					//n_to_clean = i + 1;
			}
		}

	}
#endif

	if (Clean_1()) {n_to_clean = 0; return;	}// filter all uars
	n_to_clean = n_to_clean2; // remaining valid
	if (CleanAllFifo()) return;// filter Fifo tables
	n_to_clean = n_to_clean2; // remaining valid
	p_cpt2g[56] += n_to_clean;
#ifdef TEST_ON_KNOWN
	if (debugtest == 4) {
		for (uint64_t i = 0; i < n_to_clean; i++) {
			uint64_t bf = to_clean[i];
			if (bf == debugvalbf) {
				cout << "seen expected bf in clean before clean2" << endl;
				//n_to_clean = i + 1;
			}
		}
}
#endif	
#ifdef DEFPHASE
	if (DEFPHASE == -3) {n_to_clean = 0;return;	}
#endif
	Clean_2();
	n_to_clean = 0;
}

void G17B3HANDLER::GoMiss0() {
	smin = gchk.smin;
	uasb3if = gchk.uasb3_1;
	nuasb3if = gchk.nuasb3_1;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}
void G17B3HANDLER::GoMiss1() {
	nmiss = 1;
	smin = gchk.smin;
	uasb3if = gchk.uasb3_1;
	nuasb3if = gchk.nuasb3_1;
	nuasb3of =gchk.nuasb3_2;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wua = gchk.andmiss1;
#ifdef DEBUGKNOWN
	cout <<Char27out(wua)<< " call do-miss1 wua nuasb3of="<< nuasb3of
		<<" nuasb3if=" << nuasb3if << endl;
#endif
	Do_miss1();
}

void G17B3HANDLER::Do_miss1() {
#ifdef DEBUGKNOWN
#endif
	if (!nuasb3of) {// subcritical in hn if solved
		G17B3HANDLER hn = *this;
		hn.CriticalLoop();// try direct in field
		if (gchk.aigstop || gchk.aigstopxy)return;
		if (gchk.ua_out_seen) {
			wua &= gchk.ua_out_seen;
			gchk.ua_out_seen = 0;
		}
	}
	uint32_t res;
	Critical2pairs();// assign 2 pairs in minirow to common cell
#ifdef DEBUGKNOWN
	wua &= gchk.puzknown_perm.bf.u32[2];
	cout << Char27out(wua) << "  do-miss1 wua" << endl;
	cout << Char27out(known_b3) << " known" << endl;
#endif
	while (bitscanforward(res, wua)) {
		int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
		G17B3HANDLER hn = *this;
		hn.known_b3 |= bit;
		hn.CriticalLoop();
		if (gchk.aigstop|| gchk.aigstopxy)return;
		if (gchk.ua_out_seen) {
			wua &= gchk.ua_out_seen;
			gchk.ua_out_seen = 0;
		}

	}
}

inline void G17B3HANDLER::AddCell_Miss2(uint32_t * t) {//uint32_t cell, int bit) {
	nuasb3of = t[2];
	wua = t[1];
	nmiss--;
	known_b3 |= 1 << t[0];
	uint32_t res;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	while (bitscanforward(res, wua)) {
		int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
		G17B3HANDLER hn = *this;
		hn.known_b3 |= bit;
		hn.CriticalLoop();
		if (gchk.aigstop || gchk.aigstopxy)return;
		if (gchk.ua_out_seen) {
			wua &= gchk.ua_out_seen;
			gchk.ua_out_seen = 0;
		}
	}
}
void G17B3HANDLER::GoMiss2( uint32_t uamin) {
	//if (p_cpt2g[12] < 10) diagh = 1;
	//if (diagh)cout << Char27out(uamin) << " min at entry miss2" << endl;
	//else return;
	smin = gchk.smin;
	if (diagh) smin.Status("go miss2");
	if (p_cpt2g[12] > 2)	return;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wactive0 = BIT_SET_27 ^ active_b3;//  active out field
	nmiss = 2;
	uasb3if = gchk.uasb3_1;
	nuasb3if = gchk.nuasb3_1;
	uasb3of = gchk.uasb3_2;
	nuasb3of = gchk.nuasb3_2;
	for (uint32_t i = 0; i < gchk.nmiss2ok; i++) {// now call 
		G17B3HANDLER hn = *this;
		hn.AddCell_Miss2(gchk.miss2ok[i]);
		if (gchk.aigstop || gchk.aigstopxy)return;
	}

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
	if (gchk.ua_out_seen)return;
	while (1) {// first shrink uas in field
		irloop = 0;
		uint32_t * tn = &uasb3if[nuasb3if], n = 0;
		register uint32_t Ra = active_b3,
			Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuasb3if; iua++) {
			register int Ru = uasb3if[iua];
			if (Ru & known_b3) continue;// already hit, forget it
			Ru &= active_b3;
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

	p_cpt2g[19]++;
	/*
	if (!clean_valid_done) {
		clean_valid_done = 1;
		myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
		if (myua) {
			NewUaB12();		aigstopxy = 1; return;
		}
	}
	p_cpt2g[62]++;
	*/
	if (zhou[1].CallMultipleB3(zhou[0], sn3->all_previous_cells, 0)) {
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

//================ part 2  band 3 processing using guas2/3
int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diag) {
	 zh_g2.isfalse_on = 0;

	*this = o;
	BF128 dca[9];
	int digitsbf = zh_g2.digitsbf;
	memcpy(dca, zh_g2.Digit_cell_Assigned, sizeof dca);
	{
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = genb12.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 0;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	if (_popcnt32(digitsbf < 8)) {
		if (diag) {
			cout << "not eight digits" << endl;
			ImageCandidats();
		}
		//return 1;// can not be one solution but must now send back a ua
	}
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	if (diag) {
		ImageCandidats();
		for (int i = 0; i < 81; i++)cout << zh_g2.grid0[i]+1;
		cout << endl;
	}
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; // modevalid is set to  1
	int ir = Full17Update();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0, diag);

	return zh_g.nsol;
}
int ZHOU::Apply17SingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched 
	BF128 Map;
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0]; 	R1 |= FD[1][0];
	BF128 R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	Map = FD[3][0];	BF128 R4 = R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = FD[4][0]; BF128 R5 = R4 & Map ; R4 |= R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map=  FD[5][0];  BF128 R6 = R5 & Map; 
	R5 |= R4 & Map; R4 |= R3 & Map;R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map = FD[6][0];  BF128 R7 = R6	& Map; R6 |= R5 & Map;
	R5 |= R4 & Map; R4 |= R3 & Map;R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map = FD[7][0];  R7 |= R6 & Map; R6 |= R5 & Map; 
	R5 |= R4 & Map; R4 |= R3 & Map;R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map = FD[8][0]; R7 |= R6 & Map; R6 |= R5 & Map; 
	R5 |= R4 & Map; R4 |= R3 & Map;R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	if ((cells_unsolved - R1).isNotEmpty()) 	return 1; // empty cells
	R1 -= R2; // now true singles
	R1 &= cells_unsolved; // these are new singles
	if (R1.isEmpty()) {// no single store pairs and more
		if (cells_unsolved.bf.u32[2]) {// use only b3
			if (R7.bf.u32[2])zh_g2.cells_for_guess = R7;
			else if (R6.bf.u32[2])zh_g2.cells_for_guess = R6;
			else if (R5.bf.u32[2])zh_g2.cells_for_guess = R5;
			else if (R4.bf.u32[2])zh_g2.cells_for_guess = R4;
			else if (R3.bf.u32[2])zh_g2.cells_for_guess = R3;
			else zh_g2.cells_for_guess = R2;
		}
		else {
			if (R7.bf.u64[0])zh_g2.cells_for_guess = R7;
			else if (R6.bf.u64[0])zh_g2.cells_for_guess = R6;
			else if (R5.bf.u64[0])zh_g2.cells_for_guess = R5;
			else if (R4.bf.u64[0])zh_g2.cells_for_guess = R4;
			else if (R3.bf.u64[0])zh_g2.cells_for_guess = R3;
			else zh_g2.cells_for_guess = R2;
		}
		return 0;
	}
	int tcells[80], ntcells = R1.Table3X27(tcells);
	for (int i = 0; i < ntcells; i++) {
		int cell = tcells[i];
		for (int idig = 0; idig < 9; idig++) {
			if (FD[idig][0].On_c(cell)) {
				Assign(idig, cell, C_To128[cell]);
				goto nextr1;
			}
		}
		return 1; // conflict with previous assign within this lot
	nextr1:;
	}
	zh_g.single_applied = 1;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (Apply17SingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Guess17(int index, int diag) {
	if (zh_g.go_back) return;
	if (index > 20) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< debugging temporary code
	if (diag) ImageCandidats();
	uint32_t xcell,cell, digit;
	if (cells_unsolved.bf.u32[2]) {// fill in priority band 3
		uint32_t w3 = zh_g2.cells_for_guess.bf.u32[2];
		bitscanforward(cell, w3);
		cell += 54;
		xcell = cell + 10;
		if (diag) {
			cout << Char27out(w3) << " w3 go for " << zh_g2.grid0[cell] + 1 << cellsFixedData[cell].pt << endl;
		}
		//if (!zh_g2.isfalse_on) {// skip if band 3 full
			//if (_popcnt32(cells_unsolved.bf.u32[2]) < 2) {
				//digit = zh_g2.grid0[cell];
				//FD[digit][0].Clear(xcell);// force false
			//}
		//}
	}
	// stop if all true 
	else {// fill band 12 highest first 
		if (!zh_g2.isfalse_on)	return; // already checked
		uint64_t w12 = zh_g2.cells_for_guess.bf.u64[0];
		bitscanforward64(xcell, w12);
		cell = From_128_To_81[xcell];
		if (diag) cout << zh_g2.grid0[cell] + 1 << cellsFixedData[cell].pt << " in band 12 index=" << index << endl;
	}
	digit = zh_g2.grid0[cell];
	// true first if possible
	if (FD[digit][0].On(xcell)) {
		if (diag) cout << digit+1 << cellsFixedData[cell].pt  << " true index=" << index << endl;
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->Compute17Next(index + 1, diag);
		if (zh_g.go_back) return;

	}
	// if first step try first false
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig][0].On(xcell)) {
			if (diag) cout << idig + 1 << cellsFixedData[cell].pt << " false "  << endl;
			if(xcell>=64)zh_g2.isfalse_on = 1;
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			mynext->SetaCom(idig, cell, xcell);
			mynext->Compute17Next(index + 1, diag);
			if (zh_g.go_back) return;
		}
	}
}
void ZHOU::Compute17Next(int index, int diag) {
	int ir = Full17Update();
	if (!ir) return;// locked 
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128 & wua = zh_g2.cells_assigned;
			int * sol = genb12.grid0;
			wua.SetAll_0();;
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			if (wua.isNotEmpty()) {// ignore true solution
				zh_g.nsol++;
				zh_g.go_back = 1;// closed anyway
			}
		}
		return;
	}
	Guess17(index, diag);// continue the process
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


	int ir = zhou[1].CallMultipleB3(zhou[0], bf, 0);
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
	// now adjust the stack count
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else {	// leave sub critical mode and enter the critical mode
		Critical2pairs();// assign 2 pairs in minirow to common cell
		CriticalLoop();
	}
}
void G17B3HANDLER::Go_Subcritical() {// nmiss to select in the critical field
	active_b3 = active_sub = smin.critbf;
	// check first if a global solution  is still possible
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
	if (aigstopxy) return;
	if (_popcnt32(bfb3) >ncluesb3) 		return;	
	if (moreuas_b3.Check(bfb3))return;
	p_cpt2g[18]++;
	/*
	if (!clean_valid_done) {
		clean_valid_done = 1;
		myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
		if (myua) {
			NewUaB12();		aigstopxy = 1; return;
		}
	}
	p_cpt2g[61]++;
	*/
	register uint32_t ir = zhou[1].CallMultipleB3(zhou[0], bfb3, 0);
	if (ir) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		if (ua && (!(ua & hh0.smin.critbf)))	ua_out_seen = ua;
		//DebugOut(zh_g2.cells_assigned);
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
	for (int i = 0; i < nclues + nclues_step; i++) {
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
	if(zp)strcpy(zp, ws);
	a_18_seen = 1;
	aigstop = 1; 

}
void GCHK::NewUaB12() {
	uint64_t cc64 = _popcnt64(myua&BIT_SET_2X);
	if (cc64 < 12) {// this should never be  => bug
		aigstop = 1;
		cerr << "ua < 12 to add clean" << endl;
		cout << endl << endl << Char2Xout(myua) << " ua < 12 to add   clean" << endl;
		return;
	}
	p_cpt2g[9]++;
	moreuasy.Add(myua);

	if (cc64 < 16) {
		if (cc64 < 14)moreuas_12_13.Add(myua);
		else if (cc64 == 14)moreuas_14.Add(myua);
		else if (cc64 == 15)moreuas_15.Add(myua);		
	}
	else if (cc64 < 18) moreuas_AB_small.Add(myua);
	else if (cc64 < 21)	moreuas_AB.Add(myua);
	else moreuas_AB_big.Add(myua);
	if (cc64 < 20) {
		register uint64_t ua_add = myua | (cc64 << 59);
		register uint32_t nua = genuasb12.nua;
		if (cc64 < 18 || nua < 800 || (cc64 == 18 && nua < 1200)) {
			genuasb12.AddUACheck(ua_add);
			//cout << genuasb12.nua << " nua genb12" << endl;
			tusb1[ntusb1++] = myua;
			p_cpt2g[31]++;
		}
	}


}


void GCHK::NewUaB3() {// new ua from final check zh_g2.cells_assigned
	BF128 ua128 = zh_g2.cells_assigned;
	register uint64_t ua12 = ua128.bf.u64[0];
	register uint32_t ua = ua128.bf.u32[2],
		cc = _popcnt32(ua),
		cc0 = (uint32_t)_popcnt64(ua12);
	//if (!cc) {// bands 1+2 not valid
		//myua = ua12;		NewUaB12();		aigstopxy = 1;
		//return;
	//}
	moreuas_b3.Add(ua);
	if (cc >3 ) {
		if ((cc0+cc) <= 15 || (cc == 4 && cc0 < 16)) {
			p_cpt2g[8]++;
			g4t_start.Add(ua128);
			g4t_b1.Add(ua128);
			g4t_b2.Add(ua128);
			g4t_clean.Add(ua128);
		}
		return;
	}
	if ( cc0 > 17) return;

	uint64_t ua54= (ua12 & BIT_SET_27) | ((ua12 & BIT_SET_B2) >> 5);

	if (cc == 2) {// one of the 27 GUA2s add to the table
		p_cpt2g[22]++;
		uint32_t i27 = myband3.GetI27(ua);
		tguas.tgua_start[i27].Add(ua12);
		tguas.tgua_b1[i27].Add(ua12);
		if (tguas.nvg2 < 512) {
			uint32_t ibloc = tguas.nvg2 >> 6, ir = tguas.nvg2 - 64 * ibloc;
			tvg64g2[ibloc].SetVect54(ua54, ir, i27);
			tguas.nvg2++;
		}
		return;
	}
	if (cc == 3) {// one of the 3 GUA3s add to the table
		p_cpt2g[23]++;
		uint32_t i9 = GetI9(ua),i36=i9+27;
		tguas.tgua_start[i36].Add(ua12);
		tguas.tgua_b1[i36].Add(ua12);
		if (tguas.nvg3 < 256) {
			uint32_t ibloc = tguas.nvg3 >> 6, ir = tguas.nvg3 - 64 * ibloc;
			tvg64g3[ibloc].SetVect54(ua54, ir, i36);
			tguas.nvg3++;
		}
		return;
	}


}
void GCHK::DebugOut(BF128 w) {
	uint32_t ua = w.bf.u32[2];
	uint64_t cc12 = _popcnt64(w.bf.u64[0]);
	uint64_t cc = _popcnt32(ua);
	int * sol = genb12.grid0;	
	if (cc12 < 18 && cc < 10 ) {
		cout << _popcnt32(ua) << " " << Char27out(ua) << " _ " << Char2Xout(zh_g2.cells_assigned.bf.u64[0])
			<< " " << cc12 << " " << p_cpt2g[61];
		if (cc == 2)cout << " " << tguas.nvg2;
		if (cc == 3)cout << " " << tguas.nvg3;
		cout << endl;
	}
	/*
	if (cc == 2 && cc12==17) {
		uint32_t x =(uint32_t) p_cpt2g[61];
		if (x == 7 || x == 322 || x == 495 || x == 614) {
			cout << "p_cpt2g[2]=" << p_cpt2g[2] << "p_cpt2g[3]=" << p_cpt2g[3]
				<< "p_cpt2g[5]=" << p_cpt2g[5]
				<< "p_cpt2g[6]=" << p_cpt2g[6] << "  p_cpt2g[58]=" << p_cpt2g[58] << endl;
		}
	}
	*/


}


void GCHK::Debugifof() {
#ifdef TEST_ON
	cout << "in field out field status)" << endl;
	for (uint32_t i = 0; i < nuasb3_1; i++)
		cout << Char27out(uasb3_1[i] )<< " if " << i << endl;
	for (uint32_t i = 0; i < nuasb3_2; i++)
		cout << Char27out(uasb3_2[i]) << " of " << i << endl;
#endif

}

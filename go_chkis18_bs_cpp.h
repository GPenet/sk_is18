#define EXLBREAK 4000000

//_______________ start  a band3 perm and UAs GUAs search 
void BINDEXN::Copy(STD_B1_2 & b) {
	ntvb=b.nvalidb;
	memcpy(tvb, b.my_validb, ntvb * sizeof tvb[0]);
	nt2 = b.nbi2;
	memcpy(t2, b.my_bi2, nt2 * sizeof t2[0]);
}
void BINDEXN::Copy_no7clues(STD_B1_2 & b) {
	nt2 = ntvb = 0;
	VALIDB *vb = tvb;
	for (uint32_t i = 0; i < b.nbi2; i++) { // all index 2
		BI2 wi = b.my_bi2[i], win = wi;
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
	int tperm18[3][3] = { {0,1,2},{0,2,1} , { 1,2,0 } };// keep  increasing order
	//if (ip) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	cout << "start ip=" << ip << endl;
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
		//cout << zsol << " traité order "
		//	<< band_order[0] << band_order[1] << band_order[2] << endl;
#ifdef DEBUGKNOWN
		puzknown_perm.Print3("known of the perm");
		uint64_t cc = _popcnt64(puzknown_perm.bf.u64[0]);
		uint32_t	cc1= _popcnt32(puzknown_perm.bf.u32[0]),
				cc2 = _popcnt32(puzknown_perm.bf.u32[0]);
		if (cc > 12|| cc1>7 || cc2>7 ||	(ip && cc>11)) {
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

		//____________ setup working index for external loop
		Copy_Check_7clues_needed();
#ifdef DEBUGINIT
		cout << "bin_b1.ntvb=" << bin_b1.ntvb << " bin_b2.ntvb=" << bin_b2.ntvb
		<< endl;
		uint64_t tcpt[2][5][2];// initial count b1 b2 3;4;5;6;7 + sum
		memset(tcpt, 0, sizeof tcpt);

		for (uint32_t i1 = 0; i1 < bin_b1.ntvb; i1++)
			tcpt[0][vab1w[i1].nval-1][0]++;
		tcpt[0][0][1] = tcpt[0][0][0];
		for (uint32_t i = 1; i < 5; i++)
			tcpt[0][i][1] = tcpt[0][i][0]+ tcpt[0][i-1][1];

		for (uint32_t i2 = 0; i2 < bin_b2.ntvb; i2++)
			tcpt[1][vab2w[i2].nval - 1][0]++;
		tcpt[1][0][1] = tcpt[1][0][0];
		for (uint32_t i = 1; i < 5; i++)
			tcpt[1][i][1] = tcpt[1][i][0] + tcpt[1][i - 1][1];

		cout << 31 + 2 * start_perm << "b1 count 3-6 base;cum  ";
		for (uint32_t i = 0; i < 5; i++) cout << tcpt[0][i][0] << ";" << tcpt[0][i][1] << "  ";
		cout << endl;
		cout << "b2 count 3-6 base;cum  ";
		for (uint32_t i = 0; i < 5; i++) cout << tcpt[1][i][0] << ";" << tcpt[1][i][1] << "  ";
		cout << endl;
		//G3_SplitBi2(1, 0, bin1, ib1, index_xy_b1, vab64b1);
		uint64_t & pp = p_cpt2g[31 + 2 * start_perm];
		if (start_perm) {// clues count band 3 strictly higher
			pp +=
				tcpt[0][0][0] * tcpt[1][4][1] //33 34 35 36 37
				+ tcpt[0][1][0] * tcpt[1][3][1] //43 44 45 46 (not 477)
				+ tcpt[0][2][0] * tcpt[1][3][1] //53 54 55 56
				+ tcpt[0][3][0] * tcpt[1][2][1]; //63 64 65 (not 666)

			if (start_perm)	pp+= tcpt[0][4][0] * tcpt[1][0][0]; //73  (not 747)
			else pp += tcpt[0][4][0] * tcpt[1][1][1]; //73   747

		}
		else {
			pp +=
				tcpt[0][0][0] * tcpt[1][4][1] //33 34 35 36 37
				+ tcpt[0][1][0] * tcpt[1][4][1] //43 44 45 46  477
				+ tcpt[0][2][0] * tcpt[1][3][1] //53 54 55 56
				+ tcpt[0][3][0] * tcpt[1][3][1] //63 64 65   666 
				+ tcpt[0][4][0] * tcpt[1][1][1]; //73   747
		}


		//return;
#endif

		//zh2b_g.GetBands(); ?????
		genb12.ValidInitGang();// set up gang bandc and gang bands ab
		genb12.BuildGang9x3();// setup ordered gang digits band C
		zh2b[0].Init_std_bands();
		if (0) {
			zh2b[0].ImageCandidats();
			return;
		}
		Go1_Collect_Uas();// next task in the process 
	}
}
void GCHK::Go1_Collect_Uas() {// catch uas and guas 
	//=========================== collect UAs  old process 
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	if (genuasb12.Initgen()) return;
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	tuguan.Init();
	genb12.SecondSockets2Setup();// collect GUA2s 
	tuguan.ng2 = tuguan.nguan;
	genb12.SecondSockets3Setup();// collect GUA3s 

#ifdef DO2X2
	GoM10GUas2x2();// find guas 2x2
#endif
#ifdef DEBUGEINIT
	cout << "\tua bands1+2   \t" << genuasb12.nua << endl;;
	//cout << "table uas" << endl;
	//uint64_t *t = genuasb12.tua;
	//uint32_t n = genuasb12.nua;
	//for (uint32_t i = 0; i < n; i++) cout << Char2Xout(t[i]) << " " << i << endl;
	cout << "\nnguas socket2  \t" << genb12.ntua2
		<< "\tnguas socket3  \t" << genb12.ntua3
		<< "  sockets  2;3\t" << tuguan.ng2 << "\t" << tuguan.nguan << endl;
	//if (1) return;
#endif
	myband3.Build_tsock();
	Go2_Ext_Loop();
}

#ifdef DO2X2
void G17B::GoM10GUas2x2() {// collect guas 2 sockets 2 digits
	//cout << "entry GoM10GUas2x2()" << endl;
	//for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
	//	cout << genb12.bands3[ib3].band << " " << ib3 << endl;
	socks2x2.SetAll_0();
	nts2x2 = 0;
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)genb12.bands3[ib3].issocket2x2.SetAll_0();
	for (int ib1 = 0, i81_1 = 0; ib1 < 3; ib1++) {
		for (int i1 = 0; i1 < 27; i1++, i81_1++) {// first pair of column and digits
			GEN_BANDES_12::SGUA2 & w1 = genb12.tsgua2[i81_1];
			int ib2 = (ib1 + 1) % 3, i81_2 = 27 * ib2;
			for (int i2 = 0; i2 < 27; i2++, i81_2++) {// secaond pair 
				GEN_BANDES_12::SGUA2 & w2 = genb12.tsgua2[i81_2];
				if (w1.digs != w2.digs)continue;
				int aig = 0;// one potentiel 2sockets 2 digits
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {// bands 3??
					BF128 wbf = genb12.bands3[ib3].guas.isguasocket2;
					int is_mini1 = wbf.On_c(i81_1), is_mini2 = wbf.On_c(i81_2);
					if (!(is_mini1 && is_mini2)) {
						if (is_mini1 || is_mini2) continue;
						if (!genb12.bands3[ib3].Is2Rows(i81_1, i81_2)) continue;
					}
					else {
						genb12.bands3[ib3].pat2x2[i81_1] =
							genb12.bands3[ib3].guas.ua_pair[i81_1] |
							genb12.bands3[ib3].guas.ua_pair[i81_2];
					}
					aig = 1;// minimum one active b3 
					genb12.bands3[ib3].issocket2x2.Set_c(i81_1);
					//cout << Char27out(genb12.bands3[ib3].pat2x2[i81_1])<<  "  pat " << ib3 <<" i81_1="<< i81_1 << endl;
				}
				if (aig) {
					//cout << " try sock4 for " << i81_1 << endl;
					int gangcols[9];// revised gangster
					memcpy(gangcols, genb12.gangb12, sizeof genb12.gangb12);
					gangcols[w1.col1] ^= w1.digs;
					gangcols[w1.col2] ^= w1.digs;
					gangcols[w2.col1] ^= w1.digs;
					gangcols[w2.col2] ^= w1.digs;
					zh2b_g.InitGangster(genb12.gangb12, gangcols);
					zh2b5_g.sizef5 = GUALIMSIZE;
					zh2b5_g.modevalid = 0;
					uint64_t solved_cells = zh2b5_g.FindUAsInit(w1.digs, 2);
					if (solved_cells) {// no digit solved true
						zh2b5_g.CollectUas5();// collect uas for this set of floors
						//cout << "zh2b5_g.CollectUas5() zh2b5_g.nuaf5="<< zh2b5_g.nuaf5 << endl;
						if (zh2b5_g.nuaf5) {
							//cout << Char2Xout(zh2b5_g.tuaf5[0].bf.u64) << " first gua" << endl;
							socks2x2.Set_c(i81_1);
							tuas2x2[nts2x2] = zh2b5_g.tuaf5[0].bf.u64;
							ts2x2[nts2x2++] = i81_1;
						}
					}
				}
			}
		}
	}
	if (0) {
		cout << "GoM10GUas2x2() nts2x2=" << nts2x2 << endl;
		socks2x2.Print3("socks2x2");
	}
}

#endif
void STD_B3::Build_tsock() {
	guas.isguasocketc2.Convert3X27to81(guas.isguasocket2);
	guas.isguasocketc3.Convert3X27to81(guas.isguasocket3);
	guas.isguasocketc2_46.Convert3X27to81(guas.isguasocket2_46);
	isguasocketc246 = guas.isguasocketc2 | guas.isguasocketc2_46;
	ntsock = 0;
	memset(index2, 255, sizeof index2); // all to unused
	memset(index246, 255, sizeof index246); // all to unused
	for (int i = 0; i < 81; i++) if (guas.isguasocketc2.On(i)) {
		index2[i] = index246[i] = ntsock;
		tsock[ntsock++] = guas.ua_pair[i];
	}
	memset(index3, 255, sizeof index3); // all to unused
	for (int i = 0; i < 81; i++) if (guas.isguasocketc3.On(i)) {
		index3[i] = ntsock;
		tsock[ntsock++] = guas.ua_triplet[i];
	}
	memset(index4_6, 255, sizeof index4_6); // all to unused
	for (int i = 0; i < 81; i++) if (guas.isguasocketc2_46.On(i)) {
		index4_6[i] = index246[i] = ntsock;
		tsock[ntsock++] = guas.ua_pair[i];
		if (ntsock == MAXSOCKB3) break;
	}
}

//________________  external loop 
void BINDEXN::Copy(BINDEXN & b) {
	ntvb = b.ntvb;
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
	cout << "FindSockets lim=" << lim << "\tnt1=" << nextl1 << "\tnt2=" << nextl2 << endl;
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
			cout << Char27out(extl1[i].bfx) << " t1";
			cout << "\t" << n_nob1 << "\t" << ntotb1 << "\t" << n_nob2 << "\t" << ntotb2
				<< "\t =" << ratio1 << "\t" << ntotb1 * ntotb2 << "\t" << ratio << endl;
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
			cout << Char27out(extl2[i].bfx) << " t2";
			cout << "\t" << n_nob1 << "\t" << ntotb1 << "\t" << n_nob2 << "\t" << ntotb2
				<< "\t =" << ratio1 << "\t" << ntotb1 * ntotb2 << "\t" << ratio << endl;
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
		BI2 wi = binw.t2[i], win = wi;
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
		BI2 wi = bin1no.t2[i], win = wi, win1;
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


void GCHK::Go2_Ext_Loop() {	//_____________ outer loop 
	// Copy valid bands in working areas
	//______________________________

	loopb1 = 0;
	uint32_t activerb1, activerb2;
	uint64_t activeloop = BIT_SET_2X;
	activerb1= activerb1 = BIT_SET_27;

	//_________________ external loop
	while (++loopb1 < 7) {
#ifdef DEBUGEXL
		if (aigstop)cout << "seen aigstop=1" << endl;
		cout << Char2Xout(activeloop) << "========= loop " << loopb1 << endl;
#endif
		if (aigstop)break;
		minratio = extlr.ratio=1000;
		uint64_t ir = FindSockets(activeloop,2);
		if (ir) 	ExtractMin(activeloop, bin_b1, bin_b2); 	
		if (!ir || minratio > 95) {
			ir = FindSockets(activeloop,3);
			if (ir)ExtractMin(activeloop, bin_b1, bin_b2);
			if (!ir || minratio > 95) {
				ir = FindSockets(activeloop,4);
				if (ir)ExtractMin(activeloop, bin_b1, bin_b2);
			}
#ifdef DEBUGEXL
			else 	cout << "finf3 minratio= " << minratio << endl;
#endif
		}
#ifdef DEBUGEXL
		else 	cout << "finf2 minratio= " << minratio << endl;
		cout << "final selection minratio=" << minratio << endl;
#endif
		if (minratio >= 100) break;
		else {
			extlr = extlw;
#ifdef DEBUGEXL
			cout << "entry count=" << extlw.noxyes << endl;
			extlw.Debug();
			cout << "shrink mode=" << extlw.mode << endl;
#endif
			if (extlw.mode == 1) {// this is a band1 X band2 Y
				ExtSplitX(bin_b1, bin_b1yes, extlw.bfx, activerb1);
#ifdef DEBUGEXL
				cout << "back X " << bin_b1yes.ntvb << " " << bin_b1.ntvb << endl;
				cout << Char27out(activerb1) << " back active b1" << endl;
#endif

#ifdef DEBUGKNOWN
				if (!kpfilt[0]) {// forget if the first is there
					if (extlw.bfx & puzknown_perm.bf.u32[0]) {//hit
						kpfilt[0] = 1;// this must be in this ip
						cout << " loopb1 mode 1 to use for known=" << loopb1 << endl;
						cout << Char27out(extlw.bfx) << " extlw.bfx b1" << endl;

						if (loopb1 == 1)
							Go2b_Ext_Loop(BIT_SET_2X | activerb1, 1);
						else  Go3(bin_b1yes, bin_b2);

					}
				}
#else

				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X | activerb1, 1);
				else  Go3(bin_b1yes, bin_b2);

#endif

				ExtSplitY(bin_b2,			extlw.tbfy, extlw.ntbfy, activerb2);
#ifdef DEBUGEXL
				cout << "back Y " << bin_b2.nt2 << " " << bin_b2.ntvb << endl;
				cout << Char27out(activerb2) << " back active b2" << endl;
#endif
			}
			else {// this is a band2 X band1 Y
				ExtSplitX(bin_b2, bin_b2yes,	extlw.bfx, activerb2);
#ifdef DEBUGEXL
				cout << "back X " << bin_b2yes.ntvb << " " << bin_b2.ntvb << endl;
				cout << Char27out(activerb2) << " back active b2" << endl;
#endif
#ifdef DEBUGKNOWN
				if (!kpfilt[0]) {// forget if the first is there
					if (extlw.bfx&puzknown_perm.bf.u32[1]) {//hit
						kpfilt[0] = 1;// this must be in this ip
						cout << " loopb1 mode 2 to use for known=" << loopb1 << endl;
						cout<<"\t\t" << Char27out(extlw.bfx) << " extlw.bfx b2" << endl;
						if (loopb1 == 1)
							Go2b_Ext_Loop(((uint64_t)activerb2 << 32) | BIT_SET_27, 2);
						else  Go3(bin_b1, bin_b2yes);

					}
				}
#else

				if (loopb1 == 1)
					Go2b_Ext_Loop(((uint64_t)activerb2 << 32) | BIT_SET_27, 2);
				else  Go3(bin_b1, bin_b2yes);

#endif


				ExtSplitY(bin_b1,	extlw.tbfy, extlw.ntbfy,  activerb1);
#ifdef DEBUGEXL
				cout <<"end loopb1 "<<loopb1<< " back Y " << bin_b1.nt2 << " " << bin_b1.ntvb << endl;
				cout << Char27out(activerb1) << " back active b1" << endl;
#endif
			}
		}
		activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
#ifdef DEBUGEXL
		cout << "count=" << extlr.noxyes << endl;
#endif
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
	}
	uint32_t loopb2 = 0;
	uint32_t activerb1, activerb2;
	while (++loopb2 < 6) {
#ifdef DEBUGEXL
		if (aigstop)cout << "seen aigstop=1" << endl;
		cout << Char2Xout(activeloop) << "==== loop first " << loopb2 << endl;
#endif
		if (aigstop)break;
		minratio = 1000;
		uint64_t ir = FindSockets(activeloop, 2);
		if (ir)  ExtractMin(activeloop, bin2_b1, bin2_b2); 
		if (!ir || minratio > 95) {
			ir = FindSockets(activeloop, 3);
			if (ir)ExtractMin(activeloop, bin2_b1, bin2_b2);
			if (!ir || minratio > 95) {
				ir = FindSockets(activeloop, 4);
				if (ir)ExtractMin(activeloop, bin2_b1, bin2_b2);
			}
#ifdef DEBUGEXL
			else 	cout << "finf3 minratio= " << minratio << endl;
#endif
		}
#ifdef DEBUGEXL
		else 	cout << "finf2 minratio= " << minratio << endl;
		cout << "final selection minratio=" << minratio << endl;
#endif
		if (minratio >= 100) break;
		else {
#ifdef DEBUGEXL
			extlw.Debug();
			cout << "shrink mode=" << extlw.mode << endl;
#endif
			if (extlw.mode == 1) {// this is a band1 X band2 Y
				ExtSplitX(bin2_b1, bin2_b1yes, extlw.bfx, activerb1);
#ifdef DEBUGEXL
				cout << "back X " << bin2_b1yes.ntvb << " " << bin2_b1.ntvb << endl;
				cout << Char27out(activerb1) << " back active b1" << endl;
#endif

#ifdef DEBUGKNOWN
				if (!kpfilt[1]) {// forget if the first is there
					if (extlw.bfx&puzknown_perm.bf.u32[0]) {//hit
						kpfilt[1] = 1;// this must be in this ip
						cout << " loopb1-2 mode 1 to use for known=" << loopb2 << endl;
						cout << "\t\t" << Char27out(extlw.bfx) << " extlw.bfx b1" << endl;
						Go3(bin2_b1yes, bin2_b2);
				}
			}
#else
				Go3(bin2_b1yes, bin2_b2);
#endif


				//Go3(bin2_b1yes, bin2_b2);
				ExtSplitY(bin2_b2, extlw.tbfy, extlw.ntbfy, activerb2);
#ifdef DEBUGEXL
cout << "back Y " << bin2_b2.nt2 << " " << bin2_b2.ntvb << endl;
				cout << Char27out(activerb2) << " back active b2" << endl;
#endif
			}
			else {// this is a band2 X band1 Y
				ExtSplitX(bin2_b2, bin2_b2yes, extlw.bfx, activerb2);
#ifdef DEBUGEXL
				cout << "back X " << bin2_b2yes.ntvb << " " << bin2_b2.ntvb << endl;
				cout << Char27out(activerb2) << " back active b2" << endl;
#endif

#ifdef DEBUGKNOWN
				if (!kpfilt[1]) {// forget if the first is there
					if (extlw.bfx&puzknown_perm.bf.u32[1]) {//hit
						kpfilt[1] = 1;// this must be in this ip
						cout << " loopb1-2 mode 2 to use for known=" << loopb2 << endl;
						cout << "\t\t" << Char27out(extlw.bfx) << " extlw.bfx b2" << endl;
						Go3(bin2_b1, bin2_b2yes);
					}
				}
#else
				Go3(bin2_b1, bin2_b2yes);
#endif
				ExtSplitY(bin2_b1, extlw.tbfy, extlw.ntbfy, activerb1);
#ifdef DEBUGEXL
				cout << "back Y " << bin2_b1.nt2 << " " << bin2_b1.ntvb << endl;
				cout << Char27out(activerb1) << " back active b1" << endl;
#endif
			}
		}
		activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
#ifdef DEBUGEXL
		cout << "first count=" << extlw.noxyes << endl;
#endif
		if (extlw.noxyes < EXLBREAK) break;
	}
	Go3(bin2_b1, bin2_b2);// last call

}

//__________________ process a subset of valid {band1;band2}

/* the sub lot is made of 2 pieces of the expansion	
	usually the smaller piece is in band 1 

	here 2 loops, outer band1 inner band 2
	for each loop, the process in cut in sub lots  
	  depending on the expansion index (2 common cells)
      valid bands are split by size
	  a pair of 2 cells index gives a step

	 UAs and GUAs tables are reduced to reach a "step size"
	 common cells and possible cells of the step are identified to optimize the process


struct INDEX_XY {
	uint64_t bf,and_g,or_g;// 2 cells common to the lot
	struct ITEM {// one of the tables per lot
		VALIDB64 * tvb;
		uint32_t ntvb;
	}titem[5];

  
*/


int GCHK::G3_SplitBi2( int mode,int kill7,BINDEXN & binw ,uint32_t ibi2,
INDEX_XY & ixyw,VALIDB64 * pvb ){// build tables per size 
	//___ sort  the lot per size of valid bands
	uint32_t shift = (mode == 1) ? 0 : 32;
	uint32_t buffer[100000],
		nv[5] = { 0,0,0,0,0 }, // count per size
		// see define ZSTx
		stv[5] = { 0,200,1000,5000,10000 },// start in buffer per size
		*psv[5]; // pointers to starts
	for (int i = 0; i < 5; i++) psv[i] = &buffer[stv[i]];
	BI2 biw = binw.t2[ibi2];
	ixyw.bf = biw.bf << shift;
	VALIDB * tvb = binw.tvb;
	uint32_t id = biw.istart, iend = biw.iend;
	for (uint32_t iv = id; iv < iend; iv++) {
		VALIDB & wv = tvb[iv];
		uint32_t nc = wv.nval;
		nc--;
		psv[nc][nv[nc]++] = iv;
	}
	if (kill7)nv[4] = 0;
	ixyw.ntotvb = 0;
	for (int isize = 4; isize >= 0; isize--) {
		if (nv[isize]) {
			ixyw.ntotvb += nv[isize]; 
			ixyw.ncluesmin = isize;	}
	}
	if (!ixyw.ntotvb)return 1;// contains only 7clues killed

	//_________ build the final tables in 64 mode in ixyw
	uint64_t	wand= BIT_SET_2X,wor=0;
	VALIDB64 * pvb64=pvb;
	uint32_t sumpvb = 0;
	for (int isize = 0; isize < 5; isize++) {//	3 to 7
		//cout << "isize " << isize << endl;
		uint32_t *psvw = psv[isize];
		INDEX_XY::ITEM &itemw = ixyw.titem[isize];
		itemw.ntvb = nv[isize];
		sumpvb += itemw.ntvb;
		itemw.sum_vb = sumpvb;
		itemw.tvb= pvb64;
		for (uint32_t j = 0; j < nv[isize]; j++) {
			//cout << " j " << j << endl;
			VALIDB vb32 = tvb[psvw[j]];// source valid band
			VALIDB64 & vb64 = *pvb64++;
			uint64_t bf = (uint64_t)vb32.bf << shift; // mode 2X
			wand &= bf;
			wor |= bf;
			vb64.bf = bf;
			vb64.nval = vb32.nval+2;
			memcpy(vb64.tval, biw.tval, sizeof biw.tval);
			memcpy(&vb64.tval[2], vb32.tval, sizeof vb32.tval);
			if (mode == 2)
				for (uint32_t itv = 0; itv < vb64.nval; itv++)
					vb64.tval[itv] += 27;
		}
	}
	ixyw.and_g = wand;
	ixyw.or_g = wor;
	return 0;
}

int CheckBf(BINDEXN & binw, uint32_t bfw) {
	register VALIDB * tvb= binw.tvb;
	for (uint32_t i = 0; i < binw.ntvb; i++) {
		if (tvb[i].bf == bfw)return i;
	}
	return -1;
}

void GCHK::Go3(BINDEXN & bin1, BINDEXN & bin2) {
	p_cpt2g[1]++;
#ifdef DEBUGL1L2
	cout << "go3\t" << bin1.ntvb << "\t " << bin2.ntvb << endl;
#endif
	if (aigstop) return;
	if ((!bin1.ntvb) || (!bin2.ntvb)) return;
#ifdef DEBUGKNOWN
	// fb2 must be in the lot
	kn_ir1 = CheckBf(bin1, puzknown_perm.bf.u32[0]);
	if(kn_ir1<0)return;
	kn_ir2 = CheckBf(bin2, puzknown_perm.bf.u32[1]);
	if (kn_ir2<0)return;
	cout << "got the right Go3 set" << endl;
	cout << Char27out(bin1.tvb[kn_ir1].bf) << " bfb1 ir1=" << kn_ir1 << endl;
	cout <<"\t\t"<< Char27out(bin2.tvb[kn_ir2].bf) << " bfb2 ir2=" << kn_ir2 << endl;

#endif

	//__________ loop on B1
	for (uint32_t ib1 = 0; ib1 < bin1.nt2; ib1++) {
#ifdef DEBUGKNOWN
		BI2 wi2 = bin1.t2[ib1];
		if ((uint32_t)kn_ir1 >= wi2.istart && (uint32_t)kn_ir1 <= wi2.iend) {
			cout << Char27out(wi2.bf) << " bf bi2 b1 kn_ir1=" << kn_ir1 
				<<" this is the expected bi2 band 1"<< endl;
		}
		else continue;
#endif

		if (aigstop) return;
		if(G3_SplitBi2(1, 0,bin1, ib1, index_xy_b1, vab64b1))
			continue;// X should never be killed
		int no7_b2 = index_xy_b1.ncluesmin > 4;

		fb1 = index_xy_b1.and_g;// including more common cells
		acb1= index_xy_b1.or_g;// all these cells can be active in the step
		Apply_Band1_Step();// reduce uas
		tuguan.DoStepb1();// reduce guas
		// reset more uas tables 				   
		moreuas_AB_small.Init();
		moreuas_AB.Init();
		moreuas_AB_big.Init();
		p_cpt2g[2]++;
#ifdef DEBUGL1L2
		cout << " bf ib1=" << ib1;
		index_xy_b1.Debug();
		cout << "\tnuas=" << ntusb1;
		cout << endl;
#endif						   //___________ loop on B2
		for (uint32_t ib2 = 0; ib2 < bin2.nt2; ib2++) {
#ifdef DEBUGKNOWN
			BI2 wi2_2 = bin2.t2[ib2];
			if ((uint32_t)kn_ir2 >= wi2_2.istart && (uint32_t)kn_ir2 <= wi2_2.iend) {
				cout <<"\t\t"<< Char27out(wi2_2.bf) << " bf bi2 b2 kn_ir2=" << kn_ir2
					<< " this is the expected bi2 band 2" << endl;
			}
			else continue;
#endif			
			if (aigstop) return;
			p_cpt2g[3]++;			
			if(G3_SplitBi2(2,no7_b2, bin2, ib2, index_xy_b2, vab64b2))
				continue;
			// check if 7 clues needed in band 1
			ntvb1go = (uint32_t)index_xy_b1.ntotvb;
			if (index_xy_b2.ncluesmin > 4)
				ntvb1go = index_xy_b1.titem[3].sum_vb;//limit 6

			fb2 = index_xy_b2.and_g;// including more common cells
			acb2 = index_xy_b2.or_g;// all these cells can be active in the step
			fb12 = fb1 | fb2;
					   
			if(Apply_Band2_Step())continue;// reduce uas
			tuguan.DoStepb2();// reduce guas
#ifdef DEBUGL1L2
			cout<<Char2Xout(fb12) << "\t bf ib2=" << ib2;
			index_xy_b2.Debug();
			cout << "\tnuas=" << ntusb2;
			cout << endl;
#endif	
			p_cpt2g[5+((ntusb2-1)>>6)]++;
			p_cpt2g[10]++;
			//_________ setup the brute force start
			nclues_step = 0;
			uint64_t w = fb12;
			uint32_t xcell;
			while (bitscanforward64(xcell, w)) {
				w ^= (uint64_t)1 << xcell;
				tclues[nclues_step++]= From_128_To_81[xcell];
			}
			tcluesxy = &tclues[nclues_step];
			zh2b_i1.ValidXY_Step(tclues, nclues_step);

			// and process the step depending on the uas size
			moreuas_AB_small.Init(); // added to tua_2b2
			n_to_clean = 0;// count "to clean" to 0
#ifdef DEBUGKNOWN
			cout << "start the main loop nuas= " << ntusb2<< " start_perm="<< start_perm << endl;
			cout << "index b1"; index_xy_b1.Debug(); cout << endl;
			cout << "index b2"; index_xy_b2.Debug(); cout << endl;
			//for (uint32_t i = 0; i < ntusb2; i++)
			//	cout << Char2Xout(tusb2[i]) << "  " << i << endl;
			//return;
#endif
#ifdef DEBUGONE
			if ( loopb1==2) {
				cout << "start the main loop nuas= " << ntusb2 << " cpt10=" << p_cpt2g[10] << endl;
				cout << "index b1"; index_xy_b1.Debug(); cout << endl;
				cout << "index b2"; index_xy_b2.Debug(); cout << endl;
				//for (uint32_t i = 0; i < ntusb2; i++)
					//cout << Char2Xout(tusb2[i]) << "  " << i << endl;
			}
			else if (loopb1   >2){
				aigstop = 1;
				return;
			}
#endif


			if (ntusb2 <= 64) Do64uas();
			else if (ntusb2 <= 128)  Do128uas();
			else if (ntusb2 <= 256)  Do256uas();
			else				 Do256uas();
			CleanAll();
		}
	}
}

void GCHK::Apply_Band1_Step() {// shrink the uas table	 no dead branch
	register uint64_t filter = fb1,
		Ra = acb1 | BIT_SET_B2;
	register uint64_t * tua = genuasb12.tua;
	register uint32_t nua = genuasb12.nua;
	ntusb1 = 0;
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint64_t Ru = tua[iua];
		if (Ru&filter) continue;
		Ru &= Ra;
		Ru |= ((uint64_t)_popcnt64(Ru) << 59);
		AddUA64(tusb1, ntusb1, Ru);
	}
}
int GCHK::Apply_Band2_Step() {// prepare the main loop
	acb12 = acb1 | acb2;
	register uint64_t Ra = acb12, filter = fb2;
	register uint64_t * tua = tusb1;
	uint32_t nua = ntusb1;
	ntusb2 = 0;
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint64_t Ru = tua[iua];
		if (!(Ru&filter)) {
			Ru &= Ra;
			if (!Ru) { p_cpt2g[4]++;	return 1; }
			Ru |= (uint64_t)_popcnt64(Ru) << 59;
			AddUA64(tusb2, ntusb2, Ru);
		}
	}
	return 0;
}
void TU_GUAN::DoStepb1() {// reduce tables after B1 
	pguabufb1 = guabufb1;
	nguanb1 = 0;
	register uint64_t F = gchk.fb1, Ra = gchk.acb1 | BIT_SET_B2;
	// go through the first table to extract still valid
	for (uint32_t ig = 0; ig < nguan; ig++) {
		GUAN & g = tguan[ig];
		if (g.killer&F)continue;
		GUAN gn = g;// catch permanent data
		gn.tua = pguabufb1;
		gn.nua = 0;
		register uint64_t killer = BIT_SET_2X;
		for (uint32_t iua = 0; iua < g.nua; iua++) {
			register uint64_t ua = g.tua[iua];
			if (ua & F)continue;
			ua &= Ra;// active cells only
			killer &= ua;	// if empty the GUA2/GUA3 is valid
			register uint64_t cc = _popcnt64(ua);
			ua |= (cc << 59);
			AddUA64(pguabufb1, gn.nua, ua);
		}
		if (!gn.nua)continue;
		gn.killer = killer;
		tguanb1[nguanb1++] = gn;
		pguabufb1 += gn.nua;
	}
}
void TU_GUAN::DoStepb2() {// guas filter => already killed/forced plus sub table
	p_cpt2g[34]++;
	ntsockets2_2 = ntsockets3_2 = ng2 = 0;
	pguabufr = guabufr;
	//vstep.SetAll_0(); vstep_ok.SetAll_0();
	register uint64_t F = gchk.fb12, Ra = gchk.acb12;
	{   // go through the stepb1 table to extract still valid
		for (uint32_t ig = 0; ig < nguanb1; ig++) {
			GUAN & g = tguanb1[ig];
			GUAN gn = g;// catch permanent data
			if (g.killer&F)continue;
			gn.tua = pguabufr;
			gn.nua = 0;
			register uint64_t killer = BIT_SET_2X;
			for (uint32_t iua = 0; iua < g.nua; iua++) {
				register uint64_t ua = g.tua[iua];
				if (ua & F)continue;
				ua &= Ra;// active cells only
				if (!ua) {// one empty forces socket to true
					if (g.ncol == 2)  tsockets2[ntsockets2_2++] = g.i81;
					else tsockets3[ntsockets3_2++] = g.i81;
					goto nextig;
				}
				if (gn.nua < 15) {
					register uint64_t cc = _popcnt64(ua);
					killer &= ua;
					//pguabufr[gn.nua++] = ua;
					ua |= (cc << 59);
					AddUA64(pguabufr, gn.nua, ua);
				}
			}
			if (!gn.nua)continue;
			gn.killer = killer;
			g2[ng2++] = gn;
			pguabufr += gn.nua;
		nextig:;
		}
	}
	p_cpt2g[35] += pguabufr - guabufr;
#ifdef DODEBUG
	if (g17b.diagbug) {
		cout << "exit DoStepb2()" << endl;
		Debug3();
	}
#endif

}


//______________main loop 64
void GCHK::Do64uas() {//<=64 uas in the step
	v64uas = maskLSB[ntusb2].u64[0];// Uas vector
	memset(vc64, 255, sizeof vc64);// all bits to 1
	uint32_t cc64;// build cells vectors A
	uint64_t biti = 1;
	for (uint32_t i = 0; i < ntusb2; i++, biti <<= 1) {
		register uint64_t Rw = tusb2[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc64[From_128_To_81[cc64]] ^= biti;
		}
	}

	// _____ apply cells vectors to  band 1 step 3_6 or 3_7

	for (uint32_t iv = 0; iv < ntvb1go; iv++) {
		VALIDB64 vb = vab64b1[iv];
		register uint64_t bf = vb.bf;// must hit all "empty"
		ZS64 & w= zs64b1[iv];
		w.bf = bf;
		register uint64_t V = v64uas;// apply clues 2_5/6
		for (uint64_t j = 0; j < vb.nval; j++)
			V &= vc64[vb.tval[j]];
		w.v = V;
	}
	// _____ apply cells vectors to  band 2 step no filter empty
	for (uint32_t iv = 0; iv < index_xy_b2.ntotvb; iv++) {
		VALIDB64 vb = vab64b2[iv];
		register uint64_t bf = vb.bf;
		ZS64 & w = zs64b2[iv];
		w.bf = bf;
		register uint64_t V = v64uas;
		for (uint64_t j = 0; j < vb.nval; j++)
			V &= vc64[vb.tval[j]];
		w.v = V;
	}
	if (1) {// stats to delete later 
		uint64_t & pp = p_cpt2g[30 + 2 * start_perm];
		if (start_perm) {// clues count band 3 strictly higher
			pp =
				index_xy_b1.titem[0].ntvb* index_xy_b2.ntotvb//33 34 35 36 37
				+ index_xy_b1.titem[1].ntvb* index_xy_b2.titem[3].sum_vb//43 44 45 46
				+ index_xy_b1.titem[3].ntvb* index_xy_b2.titem[2].sum_vb;//63 64 65
			if (start_perm == 2)
				pp += index_xy_b1.titem[4].ntvb* index_xy_b2.titem[0].ntvb;//73
			else pp += index_xy_b1.titem[4].ntvb* index_xy_b2.titem[1].sum_vb;//73  74
		}
		else {
			pp =
				index_xy_b1.titem[0].ntvb* index_xy_b2.ntotvb//33 34 35 36 37
				+ index_xy_b1.titem[1].ntvb* index_xy_b2.titem[3].sum_vb//43 44 45 47
				+ index_xy_b1.titem[2].ntvb* index_xy_b2.titem[3].sum_vb//53 54 55 56
				+ index_xy_b1.titem[3].ntvb* index_xy_b2.titem[3].sum_vb//63 64 65 66
				+ index_xy_b1.titem[4].ntvb* index_xy_b2.titem[1].sum_vb;//73  74
		}

	}
	// start chunks one size band2 all bands1 to use

	if (start_perm) {// clues count band 3 strictly higher
		if (index_xy_b1.titem[0].ntvb) // valid bands size 3 in band 1 all band 2
			DoChunk64(zs64b1, zs64b2, index_xy_b1.titem[0].ntvb, index_xy_b2.ntotvb);
		if (index_xy_b1.titem[1].ntvb&& index_xy_b2.titem[3].sum_vb) // size 4 in band 1  band 2 <7
			DoChunk64(&zs64b1[index_xy_b1.titem[0].sum_vb], zs64b2,
				index_xy_b1.titem[1].ntvb, index_xy_b2.titem[3].sum_vb);
		if (index_xy_b1.titem[2].ntvb && index_xy_b2.titem[3].sum_vb) //  size 5 in band 1  band 2 <7
			DoChunk64(&zs64b1[index_xy_b1.titem[1].sum_vb], zs64b2,
				index_xy_b1.titem[2].ntvb, index_xy_b2.titem[3].sum_vb);
		if (index_xy_b1.titem[3].ntvb&&index_xy_b2.titem[2].sum_vb)// size 6 in band 1  band 2 < 6
			DoChunk64(&zs64b1[index_xy_b1.titem[2].sum_vb], zs64b2,
				index_xy_b1.titem[3].ntvb, index_xy_b2.titem[2].sum_vb);
		if (start_perm == 2) {
			if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[0].ntvb) //  size 7 band 1  3 band  
				DoChunk64(&zs64b1[index_xy_b1.titem[3].sum_vb], zs64b2,
					index_xy_b1.titem[4].ntvb, index_xy_b2.titem[0].ntvb);
		}
		else {
			if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[1].sum_vb) //  size 7 band 1  3/4 band2  
				DoChunk64(&zs64b1[index_xy_b1.titem[3].sum_vb], zs64b2,
					index_xy_b1.titem[4].ntvb, index_xy_b2.titem[1].sum_vb);
		}

	}

	else { //higher or equal for band 3 add 747 477 666  (666 worst case)
		if (index_xy_b1.titem[1].sum_vb) // size 3/4 band 1 all band 2
			DoChunk64(zs64b1, zs64b2, index_xy_b1.titem[1].sum_vb, index_xy_b2.ntotvb);
		if (index_xy_b1.titem[2].ntvb && index_xy_b2.titem[3].sum_vb) //  size 5 in band 1  band 2 <7
			DoChunk64(&zs64b1[index_xy_b1.titem[1].sum_vb], zs64b2,
				index_xy_b1.titem[2].ntvb, index_xy_b2.titem[3].sum_vb);
		if (index_xy_b1.titem[3].ntvb&&index_xy_b2.titem[3].sum_vb)// size 6 in band 1  band 2 < 7
			DoChunk64(&zs64b1[index_xy_b1.titem[2].sum_vb], zs64b2,
				index_xy_b1.titem[3].ntvb, index_xy_b2.titem[3].sum_vb);
		if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[1].sum_vb) //  size 7 band 1  3/4 band2  
			DoChunk64(&zs64b1[index_xy_b1.titem[3].sum_vb], zs64b2,
				index_xy_b1.titem[4].ntvb, index_xy_b2.titem[1].sum_vb);
	}



}
void  GCHK::DoChunk64(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb) {
	if (aigstop) return;
#ifdef DEBUGKNOWN
	cout << "chunk " << na << " x " << nb << " = " << na * nb << endl;
#endif
	ZS64 * z1, *z2;
	uint64_t n1, n2;
	if (na < nb) { z1 = a; z2 = b; n1 = na; n2 = nb; }
	else { z1 = b; z2 = a; n1 = nb; n2 = na; }
	if ((nb * na) < 5000) Do64uas_11(z1, z2, n1, n2);
	else {// cut in chunks max Xchunk Ychunk
		uint64_t  ideb2 = 0, iend2 = YCHUNK64;
		if (iend2 > n2)iend2 = n2;
		while (ideb2 < n2) { //Y chunk
			uint64_t ny = iend2 - ideb2;
			uint64_t ideb1 = 0, iend1 = XCHUNK64;
			if (iend1 > n1)iend1 = n1;

			while (ideb1 < n1) {// X chunk
				uint64_t nx = iend1 - ideb1;
				Do64uas_11(&z1[ideb1], &z2[ideb2], nx, ny);
				ideb1 = iend1; iend1 += XCHUNK64;
				if (iend1 > n1)iend1 = n1;
			}
			ideb2 = iend2; iend2 += YCHUNK64;
			if (iend2 > n2)iend2 = n2;
		}
	}
}
inline void GCHK::Do64uas_11(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb) {
	//check a matrix band 1 band2 for potential 2 bands valid 11 clues
#ifdef DODEBUG
	p_cpt2g[50] ++;
#endif
	register ZS64 * Ra = &a[na - 1];
	register uint64_t * Rs = &to_clean[n_to_clean];
	for (; Ra >= a; Ra--) {
		register ZS64 * Rb = &b[nb - 1];
		register uint64_t va = (*Ra).v, bfa = (*Ra).bf | fb12;
		for (; Rb >= b; Rb--)
			if (!(Rb->v&va))		*Rs++ = bfa | Rb->bf;
	}
	n_to_clean = Rs - to_clean;
	if (n_to_clean > 10000)CleanAll();
}

//______________main loop 128
void GCHK::Do128uas() {//<=128 uas in the step
	ntua_128 = ntusb2;
	if (ntusb2 > 128)ntua_128 = 128;
	v128uas = maskLSB[ntua_128];// Uas vector
	memset(vc128, 255, sizeof vc128);// all bits to 1
	uint32_t cc64;// build cells vectors A
	for (uint32_t i = 0; i < ntua_128; i++) {
		register uint64_t Rw = tusb2[i] & BIT_SET_2X;
		//if (locdiag>1) cout << Char2Xout(Rw) << " ua i="<<i << endl;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc128[From_128_To_81[cc64]].clearBit(i);
		}
	}
	// _____ apply cells vectors to  band 1 step

	for (uint32_t iv = 0; iv < ntvb1go; iv++) {
		VALIDB64 vb = vab64b1[iv];
		register uint64_t bf = vb.bf;// must hit all "empty"
		ZS128 & w = zs128b1[iv];
		w.bf = bf;
		BF128 V = gchk.v128uas;// apply clues 2_5/6
		for (uint64_t j = 0; j < vb.nval; j++)
			V &= gchk.vc128[vb.tval[j]];
		w.v = V;
	}

	// _____ apply cells vectors to  band 2 step no filter empty

	for (uint32_t iv = 0; iv < index_xy_b2.ntotvb; iv++) {
		VALIDB64 vb = vab64b2[iv];
		register uint64_t bf = vb.bf;// must hit all "empty"
		ZS128 & w = zs128b2[iv];
		w.bf = bf;
		BF128 V = gchk.v128uas;
		for (uint64_t j = 0; j < vb.nval; j++)
			V &= gchk.vc128[vb.tval[j]];
		w.v = V;
	}

	//if (p_cpt2g[10] == 1147) {
		//cout << "1147 analysis  6 5   35;340" << endl;
		//ZS128 * w128 = &zs128b1[index_xy_b1.titem[2].sum_vb];
		//for (uint32_t i1 = 0; i1 < index_xy_b1.titem[3].ntvb; i1++) {
			//ZS28[i1].dump();
		//}
		//for (uint32_t i2 = 0; i2 < index_xy_b2.titem[2].sum_vb; i2++) {
		//	zs128b2[i2].dump();
		//}
	//}

	if (1) {// stats to delete later 
		uint64_t & pp = p_cpt2g[30 + 2 * start_perm];
		if (start_perm) {// clues count band 3 strictly higher
			pp =
				index_xy_b1.titem[0].ntvb* index_xy_b2.ntotvb//33 34 35 36 37
				+ index_xy_b1.titem[1].ntvb* index_xy_b2.titem[3].sum_vb//43 44 45 46
				+ index_xy_b1.titem[3].ntvb* index_xy_b2.titem[2].sum_vb;//63 64 65
			if (start_perm == 2)
				pp += index_xy_b1.titem[4].ntvb* index_xy_b2.titem[0].ntvb;//73
			else pp += index_xy_b1.titem[4].ntvb* index_xy_b2.titem[1].sum_vb;//73  74
		}
		else {
			pp =
				index_xy_b1.titem[0].ntvb* index_xy_b2.ntotvb//33 34 35 36 37
				+ index_xy_b1.titem[1].ntvb* index_xy_b2.titem[3].sum_vb//43 44 45 47
				+ index_xy_b1.titem[2].ntvb* index_xy_b2.titem[3].sum_vb//53 54 55 56
				+ index_xy_b1.titem[3].ntvb* index_xy_b2.titem[3].sum_vb//63 64 65 66
				+ index_xy_b1.titem[4].ntvb* index_xy_b2.titem[1].sum_vb;//73  74
		}

	}
	// start chunks one size band2 all bands1 to use

	if (start_perm) {// clues count band 3 strictly higher
		if (index_xy_b1.titem[0].ntvb) // size 3  band 1 all band 2
			DoChunk128(zs128b1, zs128b2, index_xy_b1.titem[0].ntvb, index_xy_b2.ntotvb);

		if (index_xy_b1.titem[1].ntvb&& index_xy_b2.titem[3].sum_vb) // size 4 in band 1  band 2 <7
			DoChunk128(&zs128b1[index_xy_b1.titem[0].sum_vb], zs128b2,
				index_xy_b1.titem[1].ntvb, index_xy_b2.titem[3].sum_vb);

		if (index_xy_b1.titem[2].ntvb && index_xy_b2.titem[3].sum_vb) //  size 5 in band 1  band 2 <7
			DoChunk128(&zs128b1[index_xy_b1.titem[1].sum_vb], zs128b2,
				index_xy_b1.titem[2].ntvb, index_xy_b2.titem[3].sum_vb);

		if (index_xy_b1.titem[3].ntvb&&index_xy_b2.titem[2].sum_vb)// size 6 in band 1  band 2 < 6
			DoChunk128(&zs128b1[index_xy_b1.titem[2].sum_vb], zs128b2,
				index_xy_b1.titem[3].ntvb, index_xy_b2.titem[2].sum_vb);
		if (start_perm == 2) {
			if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[0].ntvb) //  size 7 band 1  3 band  
				DoChunk64(&zs64b1[index_xy_b1.titem[3].sum_vb], zs64b2,
					index_xy_b1.titem[4].ntvb, index_xy_b2.titem[0].ntvb);

		}
		else {
			if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[1].sum_vb) //  size 7 band 1  3/4 band2  
				DoChunk128(&zs128b1[index_xy_b1.titem[3].sum_vb], zs128b2,
					index_xy_b1.titem[4].ntvb, index_xy_b2.titem[1].sum_vb);
		}
	}

	else { //higher or equal for band 3 add 747 477 666  (666 worst case)
		if (index_xy_b1.titem[1].sum_vb) // size 3/4 band 1 all band 2
			DoChunk128(zs128b1, zs128b2, index_xy_b1.titem[1].sum_vb, index_xy_b2.ntotvb);

		if (index_xy_b1.titem[2].ntvb && index_xy_b2.titem[3].sum_vb) //  size 5 in band 1  band 2 <7
			DoChunk128(&zs128b1[index_xy_b1.titem[1].sum_vb], zs128b2,
				index_xy_b1.titem[2].ntvb, index_xy_b2.titem[3].sum_vb);

		if (index_xy_b1.titem[3].ntvb&&index_xy_b2.titem[3].sum_vb)// size 6 in band 1  band 2 < 7
			DoChunk128(&zs128b1[index_xy_b1.titem[2].sum_vb], zs128b2,
				index_xy_b1.titem[3].ntvb, index_xy_b2.titem[3].sum_vb);

		if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[1].sum_vb) //  size 7 band 1  3/4 band2  
			DoChunk128(&zs128b1[index_xy_b1.titem[3].sum_vb], zs128b2,
				index_xy_b1.titem[4].ntvb, index_xy_b2.titem[1].sum_vb);
	}



}
void  GCHK::DoChunk128(ZS128 * a, ZS128 * b, uint64_t na, uint64_t nb) {
	if (aigstop) return;
#ifdef DEBUGKNOWN
	cout << "chunk " << na << " x " << nb << " = " << na * nb << endl;
#endif
	ZS128 * z1, *z2;
	uint64_t n1, n2;
	if (na < nb) { z1 = a; z2 = b; n1 = na; n2 = nb; }
	else { z1 = b; z2 = a; n1 = nb; n2 = na; }
	if ((nb * na) < 5000)		Do128uas_11(z1, z2, n1, n2);
	else {// cut in chunks max Xchunk Ychunk
		uint64_t  ideb2 = 0, iend2 = YCHUNK128;
		if (iend2 > n2)iend2 = n2;
		while (ideb2 < n2) { //Y chunk
			uint64_t ny = iend2 - ideb2;
			uint64_t ideb1 = 0, iend1 = XCHUNK128;
			if (iend1 > n1)iend1 = n1;

			while (ideb1 < n1) {// X chunk
				uint64_t nx = iend1 - ideb1;
				Do128uas_11(&z1[ideb1], &z2[ideb2], nx, ny);
				ideb1 = iend1; iend1 += XCHUNK128;
				if (iend1 > n1)iend1 = n1;
			}
			ideb2 = iend2; iend2 += YCHUNK128;
			if (iend2 > n2)iend2 = n2;
		}
	}
}
inline void GCHK::Do128uas_11(ZS128 * a, ZS128 * b, uint64_t na, uint64_t nb) {
	//check a matrix band 1 band2 for potential 2 bands valid 11 clues
	register ZS128 * Ra = &a[na - 1];
	register uint64_t * Rs = &to_clean[n_to_clean];
	for (; Ra >= a; Ra--) {
		register ZS128 * Rb = &b[nb - 1];
		ZS128 za = *Ra;
		BF128 va = za.v;
		register uint64_t  bfa = za.bf | fb12;
		for (; Rb >= b; Rb--) 	if ((Rb->v&va).isEmpty())
			//Rs++;
			*Rs++ = bfa | Rb->bf;
	}
	n_to_clean = Rs - to_clean;
	if (n_to_clean > 10000)CleanAll();
}

//______________main loop 256
void GCHK::Do256uas() {// > 128  <=256 uas in the step
	ntua_128 = 128;
	ntua_256 = ntusb2 - 128;
	if (ntusb2 > 256)ntua_256 = 128;
	v128uas = maskLSB[128];// Uas vector
	v256uas = maskLSB[ntua_256];// Uas vector
	memset(vc128, 255, sizeof vc128);// all bits to 1
	memset(vc256, 255, sizeof vc256);// all bits to 1
	// ___________________build first vector


	uint32_t cc64;// build cells vectors A
	for (uint32_t i = 0; i < 128; i++) {
		register uint64_t Rw = tusb2[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc128[From_128_To_81[cc64]].clearBit(i);
		}
	}
	// ___________________build second vector
	for (uint32_t i = 128, j = 0; j < ntua_256; i++, j++) {
		register uint64_t Rw = tusb2[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc256[From_128_To_81[cc64]].clearBit(j);
		}
	}

	// _____ apply cells vectors to  band 1 step
	uint32_t nbfb1 = (uint32_t)_popcnt64(fb1);
	//nzs1_6 = nzs1_5 = 0;
	for (uint32_t iv = 0; iv < ntvb1go; iv++) {
		VALIDB64 vb = vab64b1[iv];
		register uint64_t bf = vb.bf;// must hit all "empty"
		ZS256 & w = zs256b1[iv];
		w.bf = bf;
		BF128 V = v128uas, V2 = v256uas;// apply clues 2_5/6
		for (uint64_t j = 0; j < vb.nval; j++) {
			//for (uint64_t j = nbfb1; j < vb.nval; j++) {
			uint32_t cell = vb.tval[j];
			V &= vc128[cell];
			V2 &= vc256[cell];
		}
		w.v = V;
		w.v2 = V2;
		//if (vb.nval == 5) zs256_1_5[nzs1_5++] = w;
		//else zs256_1_6[nzs1_6++] = w;
	}

	// _____ apply cells vectors to  band 2 step no filter empty

	for (uint32_t iv = 0; iv < index_xy_b2.ntotvb; iv++) {
		VALIDB64 vb = vab64b2[iv];
		register uint64_t bf = vb.bf;// must hit all "empty"

		ZS256  & w = zs256b2[iv];
		w.bf = bf;
		BF128 V = v128uas, V2 = v256uas;
		for (uint64_t j = 0; j < vb.nval; j++) {
			//for (uint64_t j = 2; j < vb.nval; j++) {
			uint32_t cell = vb.tval[j];
			V &= vc128[cell];
			V2 &= vc256[cell];
		}
		w.v = V;
		w.v2 = V2;
	}
	if (1) {// stats to delete later 
		uint64_t & pp = p_cpt2g[30 + 2 * start_perm];
		if (start_perm) {// clues count band 3 strictly higher
			pp =
				index_xy_b1.titem[0].ntvb* index_xy_b2.ntotvb//33 34 35 36 37
				+ index_xy_b1.titem[1].ntvb* index_xy_b2.titem[3].sum_vb//43 44 45 46
				+ index_xy_b1.titem[3].ntvb* index_xy_b2.titem[2].sum_vb;//63 64 65
			if (start_perm==2)
				pp+= index_xy_b1.titem[4].ntvb* index_xy_b2.titem[0].ntvb;//73
			else pp+= index_xy_b1.titem[4].ntvb* index_xy_b2.titem[1].sum_vb;//73  74
		}
		else {
			pp =
				index_xy_b1.titem[0].ntvb* index_xy_b2.ntotvb//33 34 35 36 37
				+ index_xy_b1.titem[1].ntvb* index_xy_b2.titem[3].sum_vb//43 44 45 47
				+ index_xy_b1.titem[2].ntvb* index_xy_b2.titem[3].sum_vb//53 54 55 56
				+ index_xy_b1.titem[3].ntvb* index_xy_b2.titem[3].sum_vb//63 64 65 66
				+ index_xy_b1.titem[4].ntvb* index_xy_b2.titem[1].sum_vb;//73  74
		}

	}
	// start chunks one size band2 all bands1 to use

	if (start_perm) {// clues count band 3 strictly higher
		if (index_xy_b1.titem[0].ntvb) // size 3  band 1 all band 2
			DoChunk256(zs256b1, zs256b2, index_xy_b1.titem[0].ntvb, index_xy_b2.ntotvb);

		if (index_xy_b1.titem[1].ntvb&& index_xy_b2.titem[3].sum_vb) // size 4 in band 1  band 2 <7
			DoChunk256(&zs256b1[index_xy_b1.titem[0].sum_vb], zs256b2,
				index_xy_b1.titem[1].ntvb, index_xy_b2.titem[3].sum_vb);

		if (index_xy_b1.titem[2].ntvb && index_xy_b2.titem[3].sum_vb) //  size 5 in band 1  band 2 <7
			DoChunk256(&zs256b1[index_xy_b1.titem[1].sum_vb], zs256b2,
				index_xy_b1.titem[2].ntvb, index_xy_b2.titem[3].sum_vb);

		if (index_xy_b1.titem[3].ntvb&&index_xy_b2.titem[2].sum_vb)// size 6 in band 1  band 2 < 6
			DoChunk256(&zs256b1[index_xy_b1.titem[2].sum_vb], zs256b2,
				index_xy_b1.titem[3].ntvb, index_xy_b2.titem[2].sum_vb);
		if (start_perm == 2) {
			if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[0].ntvb) //  size 7 band 1  3 band  
				DoChunk256(&zs256b1[index_xy_b1.titem[3].sum_vb], zs256b2,
					index_xy_b1.titem[4].ntvb, index_xy_b2.titem[0].ntvb);
		}
		else {
			if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[1].sum_vb) //  size 7 band 1  3/4 band2  
				DoChunk256(&zs256b1[index_xy_b1.titem[3].sum_vb], zs256b2,
					index_xy_b1.titem[4].ntvb, index_xy_b2.titem[1].sum_vb);
		}
	}

	else { //higher or equal for band 3 add 747 477 666  (666 worst case)
		if (index_xy_b1.titem[1].sum_vb) // size 3/4 band 1 all band 2
			DoChunk256(zs256b1, zs256b2, index_xy_b1.titem[1].sum_vb, index_xy_b2.ntotvb);

		if (index_xy_b1.titem[2].ntvb && index_xy_b2.titem[3].sum_vb) //  size 5 in band 1  band 2 <7
			DoChunk256(&zs256b1[index_xy_b1.titem[1].sum_vb], zs256b2,
				index_xy_b1.titem[2].ntvb, index_xy_b2.titem[3].sum_vb);

		if (index_xy_b1.titem[3].ntvb&&index_xy_b2.titem[3].sum_vb)// size 6 in band 1  band 2 < 7
			DoChunk256(&zs256b1[index_xy_b1.titem[2].sum_vb], zs256b2,
				index_xy_b1.titem[3].ntvb, index_xy_b2.titem[3].sum_vb);

		if (index_xy_b1.titem[4].ntvb && index_xy_b2.titem[1].sum_vb) //  size 7 band 1  3/4 band2  
			DoChunk256(&zs256b1[index_xy_b1.titem[3].sum_vb], zs256b2,
				index_xy_b1.titem[4].ntvb, index_xy_b2.titem[1].sum_vb);
	}

}
void  GCHK::DoChunk256(ZS256 * a, ZS256 * b, uint64_t na, uint64_t nb) {
	if (aigstop) return;
#ifdef DEBUGKNOWN
	cout << "chunk " << na << " x " << nb << " = " << na * nb << endl;
#endif
	ZS256 * z1, *z2;
	uint64_t n1, n2;
	if (na < nb) { z1 = a; z2 = b; n1 = na; n2 = nb; }
	else { z1 = b; z2 = a; n1 = nb; n2 = na; }
	if ((nb * na) < 5000)
		//if (nb < 300 && na < 300)
		Do256uas_11(z1, z2, n1, n2);
	else {// cut in chunks max Xchunk Ychunk
		uint64_t  ideb2 = 0, iend2 = YCHUNK256;
		if (iend2 > n2)iend2 = n2;
		while (ideb2 < n2) { //Y chunk
			uint64_t ny = iend2 - ideb2;
			uint64_t ideb1 = 0, iend1 = XCHUNK256;
			if (iend1 > n1)iend1 = n1;

			while (ideb1 < n1) {// X chunk
				uint64_t nx = iend1 - ideb1;
				Do256uas_11(&z1[ideb1], &z2[ideb2], nx, ny);
				ideb1 = iend1; iend1 += XCHUNK256;
				if (iend1 > n1)iend1 = n1;
			}
			ideb2 = iend2; iend2 += YCHUNK256;
			if (iend2 > n2)iend2 = n2;
		}
	}
}
inline void GCHK::Do256uas_11(ZS256 * a, ZS256 * b, uint64_t na, uint64_t nb) {
	//check a matrix band 1 band2 for potential 2 bands valid 11 clues
	register ZS256 * Ra = &a[na - 1];
	register uint64_t * Rs = &to_clean[n_to_clean];
	for (; Ra >= a; Ra--) {
		register ZS256 * Rb = &b[nb - 1];
		ZS256 za = *Ra;
		BF128 va = za.v, va2 = za.v2;
		register uint64_t  bfa = za.bf | fb12;
		for (; Rb >= b; Rb--) {
			if ((Rb->v&va).isEmpty() && (Rb->v2&va2).isEmpty())
				*Rs++ = bfa | Rb->bf;
		}
	}
	n_to_clean = Rs - to_clean;
	if (n_to_clean > 10000)CleanAll();
}

//__________________ end of phase 1 process the file and clean
inline void TU_GUAN::SetValidV23() {//__ build valid tables 2 3
	register uint64_t F = gchk.wb12bf;

	//___set vector for used sockets2
	ntsockets2 = ntsockets2_2;
	ntsockets3 = ntsockets3_2;
	for (uint32_t ig = 0; ig < ng2; ig++) {
		GUAN gu = g2[ig];
		if (gu.killer & F) continue;
		//if (gu.nua < 2) 			tsockets2[ntsockets2++]=gu.i81;
		//else 
		for (uint32_t iua = 0; iua < gu.nua; iua++) {
			if (!(gu.tua[iua] & F)) {
				if (gu.ncol == 2)	tsockets2[ntsockets2++] = gu.i81;
				else tsockets3[ntsockets3++] = gu.i81;
				break; // first not hit is enough
			}
		}
	}
}
int GCHK::Is_B12_Not_Unique() {
	myua = zh2b[0].Valid_XY(tcluesxy, nclues);
	if (myua) {//not unique do the best with the UA
		uint64_t cc64 = _popcnt64(myua&BIT_SET_2X);
		if (cc64 < 12) {// this should never be check for a bug
			cerr << endl << endl << Char2Xout(myua) << " ua < 12 to add   clean" << endl;

#ifdef DEBUGONE
			cout << endl << endl << Char2Xout(myua) << " ua < 12 to add   clean" << endl;
			cout << Char2Xout(wb12bf) << " b12 at call ntusb2=" << ntusb2 << " stepp_cpt2g[10] " << p_cpt2g[10] << endl;
			for (int i = 0; i < nclues_step; i++) cout << tclues[i] << " ";
			cout << "\t";
			for (int i = 0; i < nclues; i++) cout << tcluesxy[i] << " ";
			cout << endl;
			zh2b_i.ImageCandidats();
			zh2b_i1.ImageCandidats();
			zh2b[0].ImageCandidats();
			cout << "table uas" << endl;
			uint64_t *t = genuasb12.tua;
			uint32_t n = genuasb12.nua;
			for (uint32_t i = 0; i < n; i++) {
				uint64_t cc = _popcnt64(t[i] & BIT_SET_2X);
				if (cc > 12)break;
				cout << Char2Xout(t[i]) << " " << i << " " << cc << endl;

			}
			for (uint32_t i = 0; i < ntusb1; i++) {
				uint64_t cc = _popcnt64(tusb1[i] & BIT_SET_2X);
				if (cc > 12)break;
				cout << Char2Xout(tusb1[i]) << " b1 i=" << i << " " << cc << endl;

			}
			for (uint32_t i = 0; i < ntusb2; i++) {
				uint64_t cc = _popcnt64(tusb2[i] & BIT_SET_2X);
				if (cc > 12)break;
				cout << Char2Xout(tusb2[i]) << " b2 i=" << i << " " << cc << endl;

			}
#endif
			aigstop = 1;
			return 1;
		}
		if (cc64 < 18) {
			moreuas_AB_small.Add(myua);
			register uint64_t ua_add = myua | (cc64 << 59);
			genuasb12.AddUACheck(ua_add);
			tusb1[ntusb1++] = myua;
			p_cpt2g[31]++;
			//if (diag_12_2020)cout << Char2Xout(myua) << " added ua"<< endl;
		}
		else if (cc64 < 21)			moreuas_AB.Add(myua);
		else moreuas_AB_big.Add(myua);
		return 1;
	}
	return 0;
}

void GCHK::CleanAll() {
	if (!n_to_clean) return;
#ifdef DEBUGKNOWN
	if (kpfilt[2]) {
		n_to_clean = 0;
		return;
	}
	//cout << " entry clean " << n_to_clean << endl;
	for (uint64_t i = 0; i < n_to_clean; i++) {
		uint64_t bf = to_clean[i];
		if (bf == puzknown_perm.bf.u64[0]) {
			cout << "\t\tclean all seen bf i="<<i<< " forced to one clean p_cpt2g[10]="<< p_cpt2g[10] <<endl;
			to_clean[0]=bf;
			n_to_clean = 1;
			kpfilt[2] = 1;
			break;
		}
	}
	if (!kpfilt[2]) {// skip if not seen
		n_to_clean = 0;
		return;
	}
		//cout << "exit DoStepb2()" << endl;
		//tuguan.Debug3();

#endif
	uint64_t nw = n_to_clean;
	n_to_clean = 0;
	for (uint64_t i = 0; i < nw; i++) {
		if (aigstop) return;
		uint64_t bf = to_clean[i];
		if (moreuas_AB_small.Check(bf))continue;
		if (moreuas_AB.Check(bf)) continue;
		if (moreuas_AB_big.Check(bf)) continue;
		wb12bf = bf;
		myband3.nclues = 18 - (uint32_t) _popcnt64(bf);

		nclues = 0;
		uint64_t w = wb12bf ^ fb12;
		uint32_t xcell;
		//cout << Char2Xout(wb12bf) << "wb12bf" << endl;
		//cout << Char2Xout(fb12) << "fb12" << endl;
		//cout << Char2Xout(w) << "w" << endl;

		while (bitscanforward64(xcell, w)) {
			w ^= (uint64_t)1 << xcell;
			tcluesxy[ nclues++]= From_128_To_81[xcell];
		}
		tuguan.SetValidV23();
		if (myband3.Clean_valid_bands3A()) {
#ifdef DEBUGKNOWN
			cout << "myband3.nclues=" << myband3.nclues;
			myband3.smin.Status(" status after Clean_valid_bands3A ");
#endif

			if (Is_B12_Not_Unique()) continue;
			else 	zhou[0].PartialInitSearch17(tclues, nclues+ nclues_step);
#ifdef DEBUGKNOWN
			cout << "check valid ok callClean_valid_bands3B " << endl;
#endif

			myband3.Clean_valid_bands3B();
		}

	}
}
//________ apply guas in band3

int STD_B3::Clean_valid_bands3A() {
	memset(&smin, 0, sizeof smin);
	nt2 = nt3 = nt46 = 0;
	for (uint32_t i2 = 0; i2 < tuguan.ntsockets2; i2++) {
		register uint32_t i81 = tuguan.tsockets2[i2];
		if (isguasocketc246.On(i81)) {
			t246[nt246++] = index246[i81];
			if (guas.isguasocketc2.On(i81)) {
				t2[nt2++] = index2[i81];
				register uint32_t	bit = guas.ua2bit[i81];
				smin.mini_bf3 |= smin.mini_bf2&bit;
				smin.mini_bf2 |= smin.mini_bf1&bit;
				smin.mini_bf1 |= bit;
				smin.critbf |= guas.ua_pair[i81];
				smin.pairs27 |= guas.ua2pair27[i81];
			}
			else t46[nt46++] = i81;
		}
	}
	for (uint32_t i3 = 0; i3 < tuguan.ntsockets3; i3++) {
		register uint32_t i81 = tuguan.tsockets3[i3];
		if (guas.isguasocketc3.On(i81)) {
			t3[nt3++] = index3[i81];
			smin.mini_triplet |= guas.ua3bit[i81];
		}
	}
	smin.SetMincount();

	if (smin.mincount > nclues)return 0; 

	register uint32_t  nout = 0, andout = BIT_SET_27;
	wactive0 = BIT_SET_27;

	if (smin.mincount > (nclues-2) ){// quick check of more
			// build free stacks pattern for outfield uas
		register int  Rfilt = smin.critbf;
		if (smin.mincount == nclues) {// first outfield is enough
			for (uint32_t i2 = 0; i2 < nt46; i2++)
				if (!(guas.ua_pair[t46[i2]] & Rfilt))
					return 0; //not ok
			for (uint32_t i = 0; i < nua; i++)
				if (!(tua[i] & Rfilt)) return 0;; //not ok
			//______________ out field guas2x2
#ifdef DO2X2
			for (int i = 0; i < g17b.nts2x2_clean; i++) {
				register int i81_1 = g17b.ts2x2_clean[i];
				if (issocket2x2.Off_c(i81_1)) continue;
				if (!(pat2x2[i81_1] & Rfilt))return 0; //not ok
			}

#endif
		}
		else {//smin.mincount = nclues-1 need min 2 out
			noutmiss1 = 0;
			register uint32_t andout = wactive0 ^ smin.critbf;
			for (uint32_t i2 = 0; i2 < nt46; i2++) {
				register uint32_t	Ru = guas.ua_pair[t46[i2]];
				if (!(Ru & Rfilt)) {
					andout &= Ru;	noutmiss1 = 1;	if (!andout) return 0; //not ok
				}
			}
			for (uint32_t i = 0; i < nua; i++) {
				register uint32_t Ru = tua[i];
				if (!(Ru & Rfilt)) {
					andout &= Ru;	noutmiss1 = 1; if (!andout) return 0; //not ok
				}
			}
#ifdef DO2X2
			for (int i = 0; i < g17b.nts2x2_clean; i++) {
				register int i81_1 = g17b.ts2x2_clean[i];
				if (issocket2x2.Off_c(i81_1)) continue;
				register uint32_t Ru = pat2x2[i81_1];
				if (!(Ru & Rfilt)) {
					andout &= Ru;	noutmiss1 = 1; if (!andout) return 0; //not ok
				}
			}
#endif
			andmiss1 = andout;
		}
	}
	else if (smin.mincount == (nclues-2)) {// need sminplus to do the right process
		register int  Rfilt = smin.critbf;
		uint32_t tout[100], nout = 0;
		register uint32_t andout = wactive0;
#ifdef DO2X2
		for (int i = 0; i < g17b.nts2x2_clean; i++) {// UA size 4 is smallest
			register int i81_1 = g17b.ts2x2_clean[i];
			if (issocket2x2.Off_c(i81_1)) continue;
			register uint32_t Ru = pat2x2[i81_1];
			if (!(Ru & Rfilt)) {
				andout &= Ru;	tout[nout++] = Ru;
			}
		}
#endif
		for (uint32_t i2 = 0; i2 < nt46; i2++) {
			register uint32_t	Ru = guas.ua_pair[t46[i2]];
			if (!(Ru & Rfilt)) {
				andout &= Ru;	tout[nout++] = Ru;
			}
		}
		for (uint32_t i = 0; i < nua; i++) {
			register uint32_t Ru = tua[i];
			if (!(Ru & Rfilt)) {
				andout &= Ru;	tout[nout++] = Ru;
			}
		}
		if (nout) {
			smin.minplus++;
			if (!andout) {// if sminplus is ncluess must have an out field 2 clues valid -
				smin.minplus++;
				register uint32_t Ru = tout[0], bit, cc, andx = 0;
				while (bitscanforward(cc, Ru)) {// look for  possible cells
					bit = 1 << cc;
					Ru ^= bit;// clear bit
					andx = BIT_SET_27;
					for (uint32_t i = 1; i < nout; i++) {
						register uint32_t Ru2 = tout[i];
						if (!(Ru2 & bit)) 	andx &= Ru2;
					}
					if (andx)break;
				}
				if (!andx) return 0;//no 2 clues killing outfield
			}
		}

	}
	return 1;
}
void STD_B3::Clean_valid_bands3B() {
	gchk.moreuas_b3_small.Init();
	gchk.moreuas_b3.Init();
	memcpy(&genb12.grid0[54], band0, 4 * 27);
	//____________________ call the relevant band 3 process

	int nmiss = nclues - smin.mincount;
#ifdef DEBUGKNOWN
	cout << "nmiss=" << nmiss << endl;
	cout << "nt2=" << nt2 << " nt3=" << nt3 << " nt46=" << nt46 << " nua=" << nua 
		<< " tot" <<nt2 + nt3 + nt46 + nua << endl;
#endif
	if ((nmiss == 2 && smin.minplus < nclues) || nmiss > 2) 
		ExpandB3();	
	else {
		gchk.smin = smin;
		G17B3HANDLER hh0;
		hh0.diagh = 0;
		if (!nmiss) { // all is infield (checked upstream)
			gchk.nuasb3_1 = 0;
			for (uint32_t i2 = 0; i2 < nt2; i2++)	gchk.uasb3_1[gchk.nuasb3_1++] = tsock[t2[i2]];
			for (uint32_t i3 = 0; i3 < nt3; i3++)	gchk.uasb3_1[gchk.nuasb3_1++] = tsock[t3[i3]];
			for (uint32_t i46 = 0; i46 < nt46; i46++)gchk.uasb3_1[gchk.nuasb3_1++] = guas.ua_pair[t46[i46]];
#ifdef DO2X2
			//______________ out field guas2x2
			for (int i = 0; i < g17b.nts2x2_clean; i++) {
				register int i81_1 = g17b.ts2x2_clean[i];
				if (issocket2x2.On_c(i81_1))
					g17b.uasb3_1[g17b.nuasb3_1++] = pat2x2[i81_1];
			}
#endif
			for (uint32_t i = 0; i < nua; i++)	gchk.uasb3_1[gchk.nuasb3_1++] = tua[i];
			hh0.GoMiss0((*this));
		}
		else if (nmiss == 1) {//checked 0 ua or andout not empty
			gchk.nuasb3_1 = gchk.nuasb3_2 = 0;
			register int  Rfilt = smin.critbf;
			for (uint32_t i2 = 0; i2 < nt2; i2++)	gchk.uasb3_1[gchk.nuasb3_1++] = tsock[t2[i2]];
			for (uint32_t i3 = 0; i3 < nt3; i3++)	gchk.uasb3_1[gchk.nuasb3_1++] = tsock[t3[i3]];
			for (uint32_t i46 = 0; i46 < nt46; i46++) {
				register uint32_t	Ru = guas.ua_pair[t46[i46]];
				if (Ru & Rfilt)	gchk.uasb3_1[gchk.nuasb3_1++] = Ru;
			}
#ifdef DO2X2
			for (int i = 0; i < g17b.nts2x2_clean; i++) {
				register int i81_1 = g17b.ts2x2_clean[i];
				if (issocket2x2.On_c(i81_1)) {
					register uint32_t Ru = pat2x2[i81_1];
					if (!(Ru & Rfilt))	g17b.uasb3_1[g17b.nuasb3_1++] = Ru;
				}
			}
#endif
			for (uint32_t i = 0; i < nua; i++) {
				register uint32_t Ru = tua[i];
				if (Ru & Rfilt)	gchk.uasb3_1[gchk.nuasb3_1++] = Ru;
			}
			// outfield is reduced to andmiss1
			hh0.GoMiss1((*this));

		}
		else {// now miss2 any out field possibility 
			hh0.GoMiss2Init((*this));
			gchk.nuasb3_1 = gchk.nuasb3_2 = 0;
			uint32_t min = 100, uamin;
			if (!hh0.known_b3) {
				register int  Rfilt = smin.critbf;
				for (uint32_t i2 = 0; i2 < nt2; i2++)gchk.uasb3_1[gchk.nuasb3_1++] = tsock[t2[i2]];
				for (uint32_t i3 = 0; i3 < nt3; i3++)gchk.uasb3_1[gchk.nuasb3_1++] = tsock[t3[i3]];
				for (uint32_t i46 = 0; i46 < nt46; i46++) {
					register uint32_t	Ru = guas.ua_pair[t46[i46]];
					if (Ru & Rfilt)
						gchk.uasb3_1[gchk.nuasb3_1++] = Ru;
					else {
						Ru &= hh0.wactive0;
						uint32_t cc = _popcnt32(Ru);
						if (cc < min) { min = cc; uamin = Ru; }
						gchk.uasb3_2[gchk.nuasb3_2++] = Ru;
					}
				}
#ifdef DO2X2
				for (int i = 0; i < g17b.nts2x2_clean; i++) {
					register int i81_1 = g17b.ts2x2_clean[i];
					if (issocket2x2.On_c(i81_1)) {
						register uint32_t Ru = pat2x2[i81_1];
						if (Ru & Rfilt)
							g17b.uasb3_1[g17b.nuasb3_1++] = Ru;
						else {
							Ru &= hh0.wactive0;
							uint32_t cc = _popcnt32(Ru);
							if (cc < min) { min = cc; uamin = Ru; }
							g17b.uasb3_2[g17b.nuasb3_2++] = Ru;
						}
					}
				}

#endif
				for (uint32_t i = 0; i < nua; i++) {
					register uint32_t Ru = tua[i];
					if (Ru & Rfilt)
						gchk.uasb3_1[gchk.nuasb3_1++] = Ru;
					else {
						Ru &= hh0.wactive0;
						uint32_t cc = _popcnt32(Ru);
						if (cc < min) { min = cc; uamin = Ru; }
						gchk.uasb3_2[gchk.nuasb3_2++] = Ru;
					}
				}

			}
			else {// same ignoring uas hit by known
				register int  Rfilt = hh0.smin.critbf,
					Rk = hh0.known_b3;
				for (uint32_t i2 = 0; i2 < nt2; i2++) {
					register uint32_t	Ru = tsock[t2[i2]];
					if (!(Ru&Rk))gchk.uasb3_1[gchk.nuasb3_1++] = Ru;
				}

				for (uint32_t i3 = 0; i3 < nt3; i3++) {
					register uint32_t	Ru = tsock[t3[i3]];
					if (!(Ru&Rk))gchk.uasb3_1[gchk.nuasb3_1++] = Ru;
				}

				for (uint32_t i46 = 0; i46 < nt46; i46++) {
					register uint32_t	Ru = guas.ua_pair[t46[i46]];
					if (Ru & Rk) continue;
					if (Ru & Rfilt)
						gchk.uasb3_1[gchk.nuasb3_1++] = Ru;
					else {
						Ru &= hh0.wactive0;
						uint32_t cc = _popcnt32(Ru);
						if (cc < min) { min = cc; uamin = Ru; }
						gchk.uasb3_2[gchk.nuasb3_2++] = Ru;
					}
				}

#ifdef DO2X2
				for (int i = 0; i < g17b.nts2x2_clean; i++) {
					register int i81_1 = g17b.ts2x2_clean[i];
					if (issocket2x2.On_c(i81_1)) {
						register uint32_t Ru = pat2x2[i81_1];
						if (Ru & Rk) continue;
						if (Ru & Rfilt)
							g17b.uasb3_1[g17b.nuasb3_1++] = Ru;
						else {
							Ru &= hh0.wactive0;
							uint32_t cc = _popcnt32(Ru);
							if (cc < min) { min = cc; uamin = Ru; }
							g17b.uasb3_2[g17b.nuasb3_2++] = Ru;
						}
					}
				}
#endif
				for (uint32_t i = 0; i < nua; i++) {
					register uint32_t Ru = tua[i];
					if (Ru & Rk) continue;
					if (Ru & Rfilt)
						gchk.uasb3_1[gchk.nuasb3_1++] = Ru;
					else {
						Ru &= hh0.wactive0;
						uint32_t cc = _popcnt32(Ru);
						if (cc < min) { min = cc; uamin = Ru; }
						gchk.uasb3_2[gchk.nuasb3_2++] = Ru;
					}
				}
			}
			hh0.GoMiss2((*this), uamin);

		}
	}

}


void G17B3HANDLER::GoMiss0(STD_B3 & b3) {
	smin = b3.smin;
	uasb3if = gchk.uasb3_1;
	nuasb3if = gchk.nuasb3_1;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}

void G17B3HANDLER::GoMiss1(STD_B3 & b3) {
	nmiss = 1;
	smin = b3.smin;
	uasb3if = gchk.uasb3_1;
	nuasb3if = gchk.nuasb3_1;
	nuasb3of = b3.noutmiss1;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wua = b3.andmiss1;
#ifdef DEBUGKNOWN
	cout <<Char27out(wua)<< " call do-miss1 wua nuasb3of="<< nuasb3of
		<<" nuasb3if=" << nuasb3if << endl;
#endif
	Do_miss1();
}
inline void G17B3HANDLER::AddCell_Miss2(uint32_t * t) {//uint32_t cell, int bit) {
	nuasb3of = t[2];
	wua = t[1];
	nmiss--;
	known_b3 |= 1 << t[0];
	Do_miss1();
}
void G17B3HANDLER::Do_miss1() {
#ifdef DEBUGKNOWN
#endif
	if (!nuasb3of) {// subcritical in hn if solved
		int uabr = IsMultiple(active_b3);
#ifdef DEBUGKNOWN
		cout << Char27out(active_b3) << " entry do-miss1 active_b3 no of "  << endl;
		cout<<Char27out(uabr) << "  do-miss1uabr"   << endl;
#endif
		if (uabr) wua = uabr;// one ua outfield seen			
		else {// confirmed subcritical possible
			G17B3HANDLER hn = *this;
			hn.Go_Subcritical();
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
		if (gchk.aigstop)return;
	}

}
void G17B3HANDLER::GoMiss2Init(STD_B3 & b3) {
	smin = b3.smin;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wactive0 = b3.wactive0 & (BIT_SET_27 ^ active_b3);//  active cells out field
}

void G17B3HANDLER::GoMiss2(STD_B3 & b3, uint32_t uamin) {
	nmiss = 2;
	uasb3if = gchk.uasb3_1;
	nuasb3if = gchk.nuasb3_1;
	uasb3of = gchk.uasb3_2;
	nuasb3of = gchk.nuasb3_2;
	if (!nuasb3of) {// subcritical in hn if solved
		wua = wactive0;
		int uabr = IsMultiple(active_b3 | known_b3);
		if (uabr) // one ua outfield seen
			wua = uabr;
		else {// confirmed subcritical possible
			G17B3HANDLER hn = *this;
			hn.Go_Subcritical();
		}
	}
	else wua = uamin;
	// cells added must produce cells hitting all remaining uas
	uint32_t res, tcellsok[27][3], ntcellsok = 0;
	while (bitscanforward(res, wua)) {
		uint32_t   nout = 0;
		register uint32_t  bit = 1 << res;
		wua ^= bit;
		wactive0 ^= bit;
		register uint32_t andx = wactive0, s = C_stack[res];
		for (uint32_t i = 0; i < nuasb3of; i++) {
			register uint32_t ua = uasb3of[i];
			if (!(ua&bit)) {
				nout = 1;	andx &= ua;	if (!andx) break;
			}
		}
		if (andx || (!nuasb3of)) {
			tcellsok[ntcellsok][0] = res;
			tcellsok[ntcellsok][1] = andx;
			tcellsok[ntcellsok++][2] = nout;

		}
	}

	for (uint32_t i = 0; i < ntcellsok; i++) {// now call 
		G17B3HANDLER hn = *this;
		hn.AddCell_Miss2(tcellsok[i]);
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
	while (1) {// first shrink uas in field
#ifdef DEBUGKNOWN
		//cout <<Char27out(known_b3) << " <<<<<<<<active critical loop start " << endl;
#endif
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
		uasb3if = tn;
		nuasb3if = n;
		if (!n) irloop = 0;// no need to loop again
		if (!irloop) break;
	}
	if (!active_b3) {// must be here expected number of clues
		if (nuasb3if) return; //can not be valid
		gchk.FinalCheckB3(known_b3);
		return; // branch closed
	}

	int wua = uasb3if[0] & active_b3, cell;
	while (bitscanforward(cell, wua)) {
		register int bit = 1 << cell;
		wua ^= bit;// clear bit
		// clean the bit in active_b3, this is now a dead cell downstream
		active_b3 ^= bit;

		G17B3HANDLER hn = *this;
		hn.CriticalAssignCell(bit);
		hn.CriticalLoop();
		if (gchk.aigstop)return;
	}
}

void STD_B3::ExpandB3() {
	//cout <<"entry expandb3"<<endl;
	uint32_t tuaw[300], nuaw = 0;
	for (uint32_t i2 = 0; i2 < nt2; i2++)	tuaw[nuaw++] = tsock[t2[i2]];
	for (uint32_t i3 = 0; i3 < nt3; i3++)	tuaw[nuaw++] = tsock[t3[i3]];
	for (uint32_t i46 = 0; i46 < nt46; i46++)tuaw[nuaw++] = guas.ua_pair[t46[i46]];
#ifdef DO2X2
	for (int i = 0; i < g17b.nts2x2_clean; i++) {
		register int i81_1 = g17b.ts2x2_clean[i];
		if (issocket2x2.On_c(i81_1)) 	tuaw[nuaw++] = pat2x2[i81_1];
	}
#endif
	for (uint32_t i = 0; i < nua; i++)		tuaw[nuaw++] = tua[i];

	//for (uint32_t i = 0; i < nuaw; i++)		cout << Char27out(tuaw[i])<<endl;

	struct SPB3 {// spots to find band 3 minimum valid solutions
		// ====================== constant after initialization
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3;
	}spb3[15], *s3, *sn3;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tuaw[0] & s3->active_cells;
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
		for (uint32_t i = s3->iuab3 + 1; i < nuaw; i++) {
			if (tuaw[i] & filter)continue;
			if (ispot >= nclues-1) 	goto next;//passing the limit
			sn3->iuab3 = i;
			register uint32_t Ru = tuaw[i] & sn3->active_cells;
			if (!Ru)goto next;
			if (ispot == nclues-2) {// last must hit all remaining uas
				for (uint32_t i2 = i + 1; i2 < nuaw; i2++) {
					if (tuaw[i2] & filter)continue;
					Ru &= tuaw[i2];
					if (!Ru)goto next;
				}
			}
			sn3->possible_cells = Ru;
			s3 = sn3; // switch to next spot
			goto next;
		}
	}
	//cout << Char27out(sn3->all_previous_cells) << "sol to check" << endl;
	// this is a possible nclues do final check
	if (zhou[1].CallMultipleB3(zhou[0], sn3->all_previous_cells, 0)) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		//cout << Char27out(ua) << "ua to add" << endl;
		if (nuaw < 300)tuaw[nuaw++] = ua;
		s3->possible_cells &= ua;
		//zhou[1].CallMultipleB3(zhou[0], sn3->all_previous_cells, 1);
	}
	else {
		gchk.Out17(sn3->all_previous_cells);
		//gchk.aigstop = 1;
		return; // stop at first  
	}
	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s3 >= spb3)goto next;
}

//================ part 2  band 3 processing using guas2/3

int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diag) {
	*this = o;
	if (diag) cout << Char27out(bf) << " call multipleb3" << endl;
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
		return 1;// can not be one solution
	}
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	if (diag) ImageCandidats();
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	int ir = Full17Update();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0, diag);

	return zh_g.nsol;
}
int ZHOU::Apply17SingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched 
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0]; 	R1 |= FD[1][0];
	BF128 R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	BF128 R4 = R3 & FD[3][0];
	R3 |= R2 & FD[3][0]; R2 |= R1 & FD[3][0]; R1 |= FD[3][0];
	BF128 R5 = R4 & FD[4][0]; R4 |= R3 & FD[4][0];
	R3 |= R2 & FD[4][0]; R2 |= R1 & FD[4][0]; R1 |= FD[4][0];
	R5 |= R4 & FD[5][0]; R4 |= R3 & FD[5][0];
	R3 |= R2 & FD[5][0]; R2 |= R1 & FD[5][0]; R1 |= FD[5][0];
	R5 |= R4 & FD[5][6]; R4 |= R3 & FD[6][0];
	R3 |= R2 & FD[6][0]; R2 |= R1 & FD[6][0]; R1 |= FD[6][0];
	R5 |= R4 & FD[7][0]; R4 |= R3 & FD[7][0];
	R3 |= R2 & FD[7][0]; R2 |= R1 & FD[7][0]; R1 |= FD[7][0];
	R5 |= R4 & FD[8][0]; R4 |= R3 & FD[8][0];
	R3 |= R2 & FD[8][0]; R2 |= R1 & FD[8][0]; R1 |= FD[8][0];
	if ((cells_unsolved - R1).isNotEmpty()) 	return 1; // empty cells
	R1 -= R2; // now true singles
	R1 &= cells_unsolved; // these are new singles
	if (R1.isEmpty()) {// no single store pairs and more
		zh_g.pairs = R2 - R3;
		zh_g2.triplets = R3 - R4;
		zh_g2.quads = R4 - R5;
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
	BF128 w = zh_g.pairs;
	if (w.isEmpty())w = zh_g2.triplets;
	if (w.isEmpty())w = zh_g2.quads;
	if (w.isEmpty())w = cells_unsolved;
	// here the target is to have ua band 3 as small as possible

	if (w.bf.u32[2]) w.bf.u32[0] = w.bf.u32[1] = 0;// select band 3 in priority
	int xcell = w.getFirst128(),
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell];

	// true first if possible
	if (FD[digit][0].On(xcell)) {
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
				if (diag) ImageCandidats();
				zh_g.nsol++;
				zh_g.go_back = 1;// closed anyway
			}
		}
		return;
	}
	Guess17(index, diag);// continue the process
}


uint32_t G17B3HANDLER::IsMultiple(int bf) {
	p_cpt2g[54] ++;
	if (bf == rknown_b3) return 0;
	if (_popcnt32(bf) > 25) return 0;
	uint32_t ua = 0;
	rknown_b3 = bf;
	GCHK & bab = gchk;
	// check first if all tuab3 is hit
	int ir = zhou[1].CallMultipleB3(zhou[0], bf, 0);
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
	if (_popcnt32(bfb3) > myband3.nclues) {
		cout << Char27out(bfb3) << " too many clues" << endl;
		return;
	}
	if (moreuas_b3_small.Check(bfb3))return;
	if (moreuas_b3.Check(bfb3))return;
	register uint32_t ir = zhou[1].CallMultipleB3(zhou[0], bfb3, 0);
	if (ir) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		if (_popcnt32(ua) < 10)moreuas_b3_small.Add(ua);
		else moreuas_b3.Add(ua);
		//cout << Char27out(ua) << " uab3 added" << endl;
		return;
	}
	Out17(bfb3);
}
void GCHK::Out17(uint32_t bfb3) {
	cout << Char32out(bfb3) << "\t\tone sol to print final check " << endl;
	char ws[82];
	strcpy(ws, empty_puzzle);
	for (int i = 0; i < nclues; i++) {
		int cell = tclues[i];
		ws[cell] = genb12.grid0[cell] + '1';
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		ws[54 + i] = genb12.grid0[54 + i] + '1';
	a_18_seen = 1;
	aigstop = 1; 

}

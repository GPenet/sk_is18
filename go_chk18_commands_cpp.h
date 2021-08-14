
//========================================
const char * libs_c17_00_cpt2g[70] = {
	"0  ",//0
	"1  entries go3 ext loop ",//1
	"2  loop b1 ",//2
	"3  loop b2",//3
	"4  active b2",
	"5  steps ",
	"6  entry n_to_clean",
	"7  pass filters",
	"8  call bf 12",
	"9  addua",
	"10 call 3b",//10
	"11 not critical 1",//11
	"12 not critical 2",//12
	"13 not critical >2",//13
	"14 fout1 count total",//14
	"15 fout1 expand",//15
	"16 critical + sub critical",//16
	"17 add 1 from active",//17
	"18 final checkb3",
	"19 ua to add upstream",
	"20 entry critical",
	"21 chehk critical",
	"22","23",
	"24 bb crit ",//24
	"25 bb miss12 ",//25
	"26 bb miss more ",//26
	"27 ",//27
	"28 XY brute force",//28
	"29 call multipleb3",//29
	"30 ip0 count1 stats ",	"31 ip0 count2",	
	"32 ip1 count1",	"33 ip1 count2",	
	"34 ip2 count1 stepb2",	"35 ip2 count2",	
	"36 stepb2 ",	"37 stepb2 ngua",	"38",	"39",
	"40 adduasb3 gua2",	
	"41 adduasb3 gua3",	
	"42 adduasb3 gua4", 
	"43",	"44",	"45",	"46",	"47",	"48",	
	"49",
	"50 somme ntusb1",	
	"51 somme ntusb2",	
	"52 somme ntusb1 steps",	
	"53 call 128 chunk",	
	"54 count matrix",
	"55 count entry clean",
	"56 count clean1",	"57",	"58",	"59",
	"60",	"61",	"62",	"63",	"64",	"65",	"66",	"67",	"68",	"69",
};

//=========================entry file of solution grids to search
void STD_B416::InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
	, int iband) {
	i416 = i16;
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

//29      9       567     7398    51516   237762  803574  =1100826
#define MAXEXP7 1200000
BI2_32 bi2x[3][250];
VALIDB vabx[3][MAXEXP7];
//VALIDB1 vab1_1[MAXEXP7];
void STD_B416::ExpandOneBand() {
	struct SPOT_E32 {// spots to find band12 valid solutions n clues
		SPOT_E32 * sp;
		uint32_t  all_previous_cells, active_cells;
		uint32_t * start_possibles, n_possibles, ipos, ispot;
		uint32_t * tua;
		uint32_t stack[3], bands[2], missing_clues, nua;
		inline void Copy(SPOT_E32 * old) {
			*this = *old;
			start_possibles += n_possibles;
			ispot++;
			missing_clues--;
			ipos = 0;
			tua += nua;
		}
		inline void AddPossibles(uint64_t v) {
			uint32_t cc;
			while (bitscanforward64(cc, v)) {// look for  possible cells
				register uint64_t bit2 = (uint64_t)1 << cc;
				v ^= (uint64_t)1 << cc;// clear bit
				start_possibles[n_possibles++] = From_128_To_81[cc];
			}
		}
		inline int GetUa(uint64_t v) {
			n_possibles = 0;
			AddPossibles(v);
			return n_possibles;
		}
		inline int GetLast() {
			n_possibles = 0;
			uint32_t andx = BIT_SET_27;
			for (uint32_t iua = 0; iua < nua; iua++) {
				andx &= tua[iua];
				if (!andx)return 0;
			}
			uint32_t cc;
			while (bitscanforward(cc, andx)) {// look for  possible cells
				andx ^= 1 << cc;// clear bit
				start_possibles[n_possibles++] = cc;
			}
			return n_possibles;
		}

	};
	uint32_t tclues[12], bufp[1000], lastspot = 6, notfirst = 0;
	uint32_t bufua12[30000];
	SPOT_E32 spt[9], *s, *sn;
	nbi2 = nvalidb = 0;
	my_bi2[0].istart = my_bi2[0].iend = 0;
	s = spt;
	memset(s, 0, sizeof spt[0]);// init the stack status ...
	s->missing_clues = 7;
	s->active_cells = BIT_SET_27;// all cells active
	s->start_possibles = bufp;
	s->tua = bufua12;// copy init table to the buffer
	s->nua = nua;
	memcpy(s->tua, tua, 4 * nua);
	s->GetUa(s->tua[0] & BIT_SET_27);

next:
	{// catch and apply cell in bitfields
		uint32_t iw = s->ipos++;
		if (iw >= s->n_possibles)goto back;
		register uint32_t cell = s->start_possibles[iw];
		tclues[s->ispot] = cell;
		register uint32_t bit = 1 << cell;
		register  uint32_t filter = s->all_previous_cells | bit,
			ac = s->active_cells ^ bit;
		sn = s + 1; sn->Copy(s); // prepare next spot
		sn->all_previous_cells = filter;
		sn->active_cells = s->active_cells = ac;
		{//level>0 shrink the ua table in new
			register uint32_t nua1 = s->nua, iua;
			sn->nua = 0;
			register uint32_t Ra = sn->active_cells;
			for (iua = 0; iua < nua1; iua++) {
				register uint32_t Ru = s->tua[iua];
				if (Ru&filter)continue;
				Ru &= Ra;
				if (!Ru) goto next;// dead branch
				register uint32_t cc = _popcnt32(Ru);
				Ru |= (cc << 27);
				AddUA32(sn->tua, sn->nua, Ru);
			}
		}
		if (s->ispot == 1) {// open a new index2
			if (notfirst) {// save previous if active
				BI2_32 & pr = my_bi2[nbi2];
				if (pr.istart != pr.iend) {
					nbi2++;
					BI2_32 & pn = my_bi2[nbi2];
					pn.istart = pn.iend = pr.iend;
				}
			}
			notfirst = 1;
			BI2_32 & pn = my_bi2[nbi2];// init the ne status
			pn.bf = sn->all_previous_cells;
			//pn.active = sn->active_cells;
			memcpy(pn.tval, tclues, sizeof pn.tval);
		}
		if (!sn->nua)goto no_more_uas;
		else if (sn->missing_clues == 1) { if (!sn->GetLast())goto next; }
		else sn->GetUa(sn->tua[0] & BIT_SET_27);
		if (s->ispot < lastspot)s++;

		goto next;
	}
no_more_uas: 	
	{	
		BI2_32 & pi = my_bi2[nbi2];
		register uint32_t R0 = sn->all_previous_cells;// ^pi.bf;
		my_validb[pi.iend++].Enter(R0, &tclues[2]);
		if (s->ispot < 6) {//if below 7 loop  for redundant clues
			int tc[32], nt = 0;
			uint32_t tbit[32];
			{	uint32_t register cell;
				register uint32_t ac = s->active_cells&BIT_SET_27;
				while (bitscanforward(cell, ac)) {// put active cells in table
					uint32_t bit = 1 << cell;
					ac ^= bit;	tc[nt] = cell;	tbit[nt++] = bit;
				}
			}
			if (s->ispot == 5) {//6 clues
				for (int i7 = 0; i7 < nt; i7++) {
					tclues[6] = tc[i7];
					my_validb[pi.iend++].Enter(tbit[i7] | R0, &tclues[2]);
				}
			}
			else if (s->ispot == 4) { //5 clues
				for (int i6 = 0; i6 < nt; i6++) {
					register uint32_t R6 = tbit[i6] | R0;
					tclues[5] = tc[i6];
					my_validb[pi.iend++].Enter(R6, &tclues[2]);
					for (int i7 = i6 + 1; i7 < nt; i7++) {
						tclues[6] = tc[i7];
						my_validb[pi.iend++].Enter(tbit[i7] | R6, &tclues[2]);
					}
				}
			}
			else if (s->ispot == 3) { // valid 4 clues
				for (int i5 = 0; i5 < nt; i5++) {
					register uint32_t R5 = tbit[i5] | R0;
					tclues[4] = tc[i5];
					my_validb[pi.iend++].Enter(R5, &tclues[2]);
					for (int i6 = i5 + 1; i6 < nt; i6++) {
						register uint32_t R6 = tbit[i6] | R5;
						tclues[5] = tc[i6];
						my_validb[pi.iend++].Enter( R6, &tclues[2]);
						for (int i7 = i6 + 1; i7 < nt; i7++) {
							tclues[6] = tc[i7];
							my_validb[pi.iend++].Enter(tbit[i7] | R6, &tclues[2]);
						}
					}
				}
			}
			else if (s->ispot == 2) { // valid 3 clues
				for (int i4 = 0; i4 < nt; i4++) {
					register uint32_t R4 = tbit[i4] | R0;
					tclues[3] = tc[i4];
					my_validb[pi.iend++].Enter(R4, &tclues[2]);
					for (int i5 = i4 + 1; i5 < nt; i5++) {
						register uint32_t R5 = tbit[i5] | R4;
						tclues[4] = tc[i5];
						my_validb[pi.iend++].Enter(R5, &tclues[2]);
						for (int i6 = i5 + 1; i6 < nt; i6++) {
							register uint32_t R6 = tbit[i6] | R5;
							tclues[5] = tc[i6];
							my_validb[pi.iend++].Enter( R6, &tclues[2]);
							for (int i7 = i6 + 1; i7 < nt; i7++) {
								tclues[6] = tc[i7];
								my_validb[pi.iend++].Enter(tbit[i7] | R6, &tclues[2]);
							}
						}
					}
				}
			}
			// band 29 not considered no valid 2 clue
		}
	}
	goto next;
back:
	if (--s >= spt)goto next;
	// save the last index if
	BI2_32 & pr = my_bi2[nbi2];
	if (pr.istart != pr.iend) 	nbi2++;
	nvalidb=pr.iend;
}


int  Is18(char * ze, char * zp) {
	gchk.ze = ze;
	gchk.zp = zp;
	gchk.GuapatsInit();
#ifdef DEBUGKNOWN
	if (strlen(ze) < 163) return -1;// skip blank lines
	char * w = &ze[82];
	gchk.puzknown.SetAll_0();
	gchk.puzknown_perm= gchk.puzknown;
	for (int i = 0; i < 81; i++) {
		char c = w[i];
		if(c<'1' || c>'9') continue;
		gchk.puzknown.Set_c(i);// store the pattern in 3X mode
		// this can be a pattern, no check of the digit with the solution
	}
#endif
	gchk.a_18_seen = 0;
	//ze[81] = 0; //kill outfield
	memset(p_cpt2g, 0, sizeof p_cpt2g);
	zh_g.modevalid = 1;
	zh_g2.grid0 = gchk.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	bin_b1.Attach(bi2_b1w, vab1w);
	bin_b2.Attach(bi2_b2w, vab2w);
	bin_b1yes.Attach(bi2_b1yes , vab1yes);
	bin_b2yes.Attach(bi2_b2yes , vab2yes);
	bin2_b1.Attach(bi2_b1w2, vab1w2);
	bin2_b2.Attach(bi2_b2w2, vab2w2);
	bin2_b1yes.Attach(bi2_b1yes2, vab1yes2);
	bin2_b2yes.Attach(bi2_b2yes2, vab2yes2);
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	int * zs0 = gchk.grid0;
	gchk.aigstop = 0;
	//gchk.DiagPuz();
	long tdeb = GetTimeMillis();
	// ====catch entry uas and rank
	for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
	BANDMINLEX::PERM perm_ret;
	bandminlex.Getmin(zs0, &perm_ret);
	int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
	if (ib1 == 28)return -1;
	bax[0].InitBand2_3(ib1, ze, perm_ret, 0);
	bandminlex.Getmin(&zs0[27], &perm_ret);
	int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
	if (ib2 == 28)return -1;
	bax[1].InitBand2_3(ib2, &ze[27], perm_ret, 1);
	bandminlex.Getmin(&zs0[54], &perm_ret);
	int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
	if (ib3 == 28)return -1;
	bax[2].InitBand2_3(ib3, &ze[54], perm_ret, 2);
	//_________________________________ expand bands 
	for (int ib = 0; ib < 3; ib++) {
		bax[ib].InitExpand(bi2x[ib], vabx[ib]);
		bax[ib].ExpandOneBand();
#ifdef DEBUGONE
		cout << ib1 << "\tib1 mins=" << t16_min_clues[ib1] << endl;
		cout << ib2 << "\tib2 mins=" << t16_min_clues[bax[1].i416] << endl;
		cout << ib3 << "\tib3 mins=" << t16_min_clues[bax[2].i416] << endl;
		cout << "end expand i=" << ib << "\tnbi2=" << bax[ib].nbi2
			<< "\tnvalid=" << bax[ib].nvalidb << " nlast ndex=" << bax[ib].my_bi2[bax[ib].nbi2 - 1].iend << endl;
#endif

#ifdef DEBUGKNOWN
#ifdef DEBUGONE
		uint32_t R = gchk.puzknown.bf.u32[ib];
		if (_popcnt32(R) > 7) continue;
		cout << Char27out(R) << " searched valid band" << endl;
		uint32_t aig = 1;
		STD_B416 &b= bax[ib];
		for (uint32_t i = 0; i < b.nbi2; i++) { // all index 2
			BI2 wi = b.my_bi2[i];
			if ((R&wi.bf) == wi.bf) {
				cout <<Char27out(wi.bf )<< "valid must be index " << i << endl;
				for (uint32_t j = wi.istart; j < wi.iend; j++) {
					VALIDB wj = b.my_validb[j];
					if ((wj.bf& (~R))) continue;
					if (!wj.bf) {
						cout << "bug empty bf located in j="<< j << endl;
						return 0;
					}
					cout << Char27out(wj.bf& (~R))<<"\t" ;

					cout << Char27out(wj.bf)<<" j=" << j << endl;
					if (R == wj.bf) {
						cout << "expected seen in in j=" << j << endl;
						aig = 0;
						break;
					}
				}
				if (aig) {
					cout << "missing valid" << endl;
					return 0;
				}
			}
		}

#endif
#endif
	}

	int tsort[3];
	{// sort entry increasing order of min clues
		for (int i = 0; i < 3; i++)
			tsort[i] = i + (t416n6[bax[i].i416] << 8);
		for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
			if (tsort[i] > tsort[j]) {
				int temp = tsort[i];
				tsort[i] = tsort[j];
				tsort[j] = temp;
			}
	}

	// collect uas for the process 
	//zh_g5.Init(ze);
	//gchk.uacollector.Collect1();
	//genb12.InitialSockets2Setup();
	//genb12.InitialSockets3Setup();
	{// now start all perms in bands a b c mode 
		for (int ip = 0; ip < 3; ip++) {// 3 bands perms
			gchk.Start(bax, tsort, ip);
			//break;
		}
	}

	cout << "print final stats" << endl;
	for (int i = 0; i < 70; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
	return gchk.a_18_seen;
}

//=== buffers to store valid bands and vectors
// max found in test all bands 243 steps 17063 valid bands

/* processing abc increasing order tperm6
une fois
666
trois fois
667 558/9 44A/B 33C/D 577 477 388
6 fois 
567/8 459/A 468/9 478 34B/C 35A/B 369/A 378/9
*/


/* UAs collector phase 1 tous en mode 3 
36 x 2 digits
84 x 3 digits
126 x 4 digits  

deux bandes vides, on les a tous
une bande vide on peut utiliser zhb1b2 
pas de bande vide ?? on peut démarrer sur ce mode

une mini row on force inversé sur général 
*/

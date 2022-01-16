//========================================
const char * libs_c17_00_cpt2g[70] = {
	"0  ",//0
	"1  uas b12 start ",//1
	"2  uas g2 start",//2
	"3  steps first 3/4 ",//3
	"4  entry clean",
	"5  pass filter ",
	"6  valid b 12 ",
	"7  look for more clues b1+2 ",
	"8  ",
	"9  all band 3 calls",
	"10 not critical 0",//10
	"11 not critical 1",//11
	"12 not critical 2",//12
	"13 b3 direct",//13

	"14 b3 expand",//14
	"15 ",//15
	"16 ",//16
	"17 ",//17
	"18  ",
	"19  final checkb3",
	"20  count add",
	"21 ",
	"22  max c2",
	"23  max cmore",
	"24 ",//24
	"25 entry add b1 b2 ",//25
	"26 add one cell b1 b2 ",//26
	"27 add one more",//27
	"28 ",//28
	"29 ",//29

	"30  new ua b12",	
	"31 add to uasb12",	
	"32 add to uasb12 +1",	"33 add to uasb12 +2",	
	"34 add to uasb12 ++",	
	
	"35 ",	"36  ",	"37 ",	"38",	"39",

	"40 ",	
	"41 limspot",	
	"42 ", 
	"43 limspot nd",	
	"44","45","46","47","48",	
	"49","50 ","51 ","52 ",	"53 ",	"54 ",	"55 ",
	"56 count clean1",	
	"57 count clean2 actif",	
	"58 clean2 > gua2",	
	"59 clean2 > gua3",
	"60 clean2 > get multiple",	
	"61 final check b3",	
	"62 final check b3 expand",	"63",	"64",	"65",	"66",	"67",	"68",	"69",
};

//=========================entry file of solution grids to search

//29      9       567     7398    51516   237762  803574  =1100826
#define MAXEXP7 1200000

void BandReShape(int * s, int * d, BANDMINLEX::PERM p) {
	int * pc = p.cols;
	for (int irow = 0; irow < 3; irow++) {
		int drow = 9 * irow;
		for (int icol = 0; icol < 9; icol++)
			d[drow + icol] = p.map[s[drow + pc[icol]]];
	}
	int temp[9];// now reorder 
	if (d[0] > d[9]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[9], sizeof temp);
		memcpy(&d[9], temp, sizeof temp);
	}
	if (d[0] > d[18]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
	if (d[9] > d[18]) {
		memcpy(temp, &d[9], sizeof temp);
		memcpy(&d[9], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
}
int  Is18Blue(char * ze, char * zp) {
	gchk.ze = ze;	gchk.zp = zp;
	char ze_diag[82]; ze_diag[81] = 0;
	memset(p_cpt2g, 0, sizeof p_cpt2g);
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	int * zs0 = gchk.grid0,*zs0_diag= gchk.grid0_diag;
	// ====catch entry uas and rank
	for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
	for (int i = 0; i < 81; i++)zs0_diag[i] = zs0[C_transpose_d[i]];
	for (int i = 0; i < 81; i++)ze_diag[i] = ze[C_transpose_d[i]];
	BANDMINLEX::PERM perm_ret;
	for (int ibs = 0; ibs < 3; ibs++) {
		bandminlex.Getmin(&zs0[27*ibs], &perm_ret);
		bax[ibs].InitBand2_3( &ze[27 * ibs], perm_ret, ibs);
		bandminlex.Getmin(&zs0_diag[27 * ibs], &perm_ret);
		bax[ibs+3].InitBand2_3(&ze_diag[27 * ibs], perm_ret, ibs);
	}
	int i416min = 500;
	for (int ibs = 0; ibs < 6; ibs++)
		if (bax[ibs].i416 < i416min)i416min = bax[ibs].i416;
	STD_B416 bw[3],myb1;
	int nfirst = 0, gend[81],gw[81];
	for (int ibs0 = 0; ibs0 < 6; ibs0++) {
		STD_B416 & b = bax[ibs0];
		int ibsw;
		if (b.i416 == i416min) {
			nfirst++;
			//cout << "check ibs=" << ibs0 << endl;
			for (int ibs = 0; ibs < 3; ibs++)
				if (ibs0 > 2) bw[ibs] = bax[ibs+3];
				else  bw[ibs] = bax[ibs];
			if (ibs0 > 2) ibsw = ibs0 - 3;
			else ibsw = ibs0;
			myb1 = bw[ibsw];
			int * ss = myb1.band0, dd[27];
			bandminlex.Getmin(myb1.band0, &perm_ret); // redo perm
			for (int i = 0; i < 27; i++) {
				int iss = 9 * perm_ret.rows[i / 9] + perm_ret.cols[i % 9];
				dd[i] = perm_ret.map[ss[iss]];
			}
			int iband2 = (ibsw + 1) % 3, iband3 = (ibsw + 2) % 3
				, dd2[27], dd3[27];
			BandReShape(bw[iband2].band0, dd2, perm_ret);
			BandReShape(bw[iband3].band0, dd3, perm_ret);
			memcpy(gw, dd, sizeof dd);
			if (dd2[0] == 1) {// must have '2' in r4c1
				memcpy(&gw[27], dd2, sizeof dd);
				memcpy(&gw[54], dd3, sizeof dd);
			}
			else {
				memcpy(&gw[27], dd3, sizeof dd);
				memcpy(&gw[54], dd2, sizeof dd);
			}
			//for (int i = 0; i < 81; i++) cout << gw[i] + 1;
			//cout << endl;
			if (nfirst > 1) {
				for (int i = 0; i < 81; i++){
					if (gend[i] < gw[i]) break;
					if (gend[i] > gw[i]) {
						nfirst = 1; break;
					}
				}
			}
			if (nfirst == 1)memcpy(gend, gw, sizeof gend);
		}

	}
	for (int i = 0; i < 81; i++) fout1 << gend[i] + 1;
	fout1 << endl;
	//for (int i = 0; i < 81; i++) cout << gend[i] + 1;
	//cout << endl;

	return 0;
}

int  Is18(char * ze, char * zp) {
	gchk.ze = ze;
	gchk.zp = zp;
	gchk.diagbugchk = 0;
	char ze_diag[164]; ze_diag[163] = 0;
	gchk.a_18_seen = 0;
	memset(p_cpt2g, 0, sizeof p_cpt2g);
	zh_g.modevalid = 1;
	zh_g2.grid0 = gchk.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	int * zs0 = gchk.grid0, *zs0_diag = gchk.grid0_diag;
	gchk.aigstop = 0;
	long tdeb = GetTimeMillis();
	// ====catch entry uas and rank
	for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
	for (int i = 0; i < 81; i++)zs0_diag[i] = zs0[C_transpose_d[i]];
	for (int i = 0; i < 81; i++)ze_diag[i] = ze[C_transpose_d[i]];
	BANDMINLEX::PERM perm_ret;
	for (int ibs = 0; ibs < 3; ibs++) {
		bandminlex.Getmin(&zs0[27 * ibs], &perm_ret);
		bax[ibs].InitBand2_3(&ze[27 * ibs], perm_ret, ibs);
		bandminlex.Getmin(&zs0_diag[27 * ibs], &perm_ret);
		bax[ibs + 3].InitBand2_3(&ze_diag[27 * ibs], perm_ret, ibs);
	}
	/*
	if (0) {
		cout << ze << "contrôle acquisition" << endl;
		int imin = 7, min = 1000000;
		for (int ibs = 0; ibs < 6; ibs++) {
			STD_B416 & b = bax[ibs];
			cout << ibs << "\t" << b.band << "\t" << b.i416 << " " << t416n6[b.i416]
				<< " index " << t416_to_n6[b.i416] << endl;
			if (t416_to_n6[b.i416] < min) {
				min = t416_to_n6[b.i416]; imin = ibs;
			}
		}
		cout << "min expected i=" << imin << endl;
		if (imin > 2)cout << "should use digonal symmetry" << endl;
	}
	*/


	int tsort[3];
	{// sort entry increasing order of min clues
		for (int i = 0; i < 3; i++)
			tsort[i] = i + (t416n6[bax[i].i416] << 8);
			//tsort[i] = i + (t16_min_clues[bax[i].i416] << 8);
		for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
			if (tsort[i] > tsort[j]) {
				int temp = tsort[i];
				tsort[i] = tsort[j];
				tsort[j] = temp;
			}
	}
	// put bands in the right order
	bax[3] = bax[tsort[0] & 7];
	bax[4] = bax[tsort[1] & 7];
	bax[5] = bax[tsort[2] & 7];
	bax[0] = bax[3]; bax[1] = bax[4]; bax[2] = bax[5];

	gchk.band_order[0] = tsort[0] & 7;
	gchk.band_order[1] = tsort[1] & 7;
	gchk.band_order[2] = tsort[2] & 7;
	gchk.mincluesb3 = 6;//t16_min_clues[bax[2].i416];

	// put knownn in the right order

#ifdef HAVEKNOWN
	if (strlen(ze) < 163) return -1;// skip blank lines
	char * w = &ze[82],ww[82];
	memcpy(ww, w, 81);
	memcpy(w, &ww[27 * (tsort[0] & 7)], 27);
	memcpy(&w[27], &ww[27 * (tsort[1] & 7)], 27);
	memcpy(&w[54], &ww[27 * (tsort[2] & 7)], 27);
	gchk.puzknown.SetAll_0();
	gchk.puzknown_perm = gchk.puzknown;
	for (int i = 0; i < 81; i++) {
		char c = w[i];
		if (c<'1' || c>'9') continue;
		gchk.puzknown.Set_c(i);// store the pattern in 3X mode
		// this can be a pattern, no check of the digit with the solution
	}
#endif
	// reshape known

	int irs= gchk.StartIs18();
#ifdef TEST_ON
	cout << "print final stats" << endl;
	for (int i = 0; i < 70; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
#endif

	return irs;
}


/*
void  GCHK::ExpandB12() {
	int mincluesb3 = 6;//<<<<<<<<<<<<<<<<<<<< must be dynamic and linked ot band3 i416
	zh2b[0].InitBands12(grid0);
	uint64_t *twu = tuasb12.tua,
		limspot = 18 - mincluesb3 - 1;// 5 clues done
	// _______________ expand to have all  minimal < n clues
	struct SPB {// spots to find starts to extract uas 3 digits
		uint64_t  possible_cells, all_previous_cells, active_cells;
	}spb[20], *s, *sn;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->active_cells = BIT_SET_2X;
	s->possible_cells = twu[0] & BIT_SET_2X;
	if (1) {// debugging of start
		limspot = 10;
	}

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
	uint64_t Ru = twu[ir + w.u_start] & ac;
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
p_cpt2g[42]++;
//cout << Char2Xout(sn->all_previous_cells) << " " << p_cpt2g[42]<<endl;
AfterExpandB12(sn->all_previous_cells, s->active_cells, (uint32_t)ispot + 1);

goto next;
// going back, for a non empty index, count it back
back:
if (--s >= spb)goto next;
cout << "count expandb12=" << p_cpt2g[41]
<< " " << p_cpt2g[42] << " 10*av chunks=" << (10 * p_cpt2g[43]) / p_cpt2g[41]
<< "\nto check=" << p_cpt2g[44] << "\n valid=" << p_cpt2g[45]
<< "\n direct expand b3=" << p_cpt2g[46] << endl;

morev2a.Status(" a 12-16");
morev2b.Status(" b 17-18  ");
morev2c.Status(" b 19-20  ");
morev2d.Status(" d 21-22  ");

}
*/

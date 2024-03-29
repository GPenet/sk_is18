//========================================
const char * libs_c17_00_cpt2g[70] = {
	"0  compte entr�e",//0
	"1  uas b12 start ",//1
	"2  uas g2 start",//2
	"3  steps first 3 ",//3
	"4  steps first 6",
	"5  raw b12 11/12 ",
	"6  check mincount",
	"7  after mincount ",
	"8  after minplus ",
	"9  init gob3 after mincount",
	"10 toass 0 1",//10
	"11       2 3",//11
	"12       4 5",//12
	"13       6 7",//13
	"14       8 9",//14
	"15      10 11",//15
	"16 entry go below",//16
	"17 count 11 go after",//17
	"18 count 12 go after  ",//18
	"19  checkb3 below 18",
	"20 entry go after ",
	"21 ",
	"22  max c2",
	"23  max cmore",
	"24 ",//24
	"25 entry add b1 b2 ",//25
	"26  redo mincount",//26
	"27 add after g2",//27
	"28 add after g3",//28
	"29 call checkb12",//29

	"30  new ua b12",	
	"31 add to <18",	
	"32 add to <20",	
	"33 add to <23",	
	"34  666",	
	
	"35 add c2 ",	
	"36 ad other ",	
	"37 filter one of 3 perms",	
	"38 expand b3 direct",	
	"39 expand b3 vector",

	"40 pass the 666 filter",	
	"41 test vecteurs g2 g3",	
	"42 debug buildapply", 
	"43 compte >=64",	
	"44 max uas b3 direct",
	"45 must recheck min",
	"46 after 67 go",
	"47 46 somme of add",
	"48 apply mv64vh",	
	"49 max mv64h",

	"50 ctl loop next",
	"51 7cl",
	"52 8cl",	
	"53 9cl",	"54 10l",	"55 11cl",
	"56 count clean1",	
	"57 count clean2 actif",	
	"58 clean2 > gua2",	
	"59 clean2 > gua3",
	"60 clean2 > get multiple",	
	"61 final check b3",	
	"62 final check b3 expand",	
	"63 max uas b3",	
	"64 clean more uas",	 
	"65 entry below check valid",	
	"66 do check valid below",	
	"67 ok check valid below",	
	"68 expand go",	"69",
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

int  Is18(char* ze, char* zp, int first) {
	gchk.ze = ze;
	gchk.zp = zp;
	gchk.aigstop = 0;
	char ze_diag[164]; ze_diag[163] = 0;
	memcpy(gchk.zes, ze, 164);
	gchk.nok = gchk.n18seen = 0;
	gchk.modefirst = first;
	memset(p_cpt2g, 0, sizeof p_cpt2g);
	zh_g.modevalid = 1;
	zh_g2.grid0 = gchk.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	int* zs0 = gchk.grid0, * zs0_diag = gchk.grid0_diag;
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

	int tsort[3], tsort2[3];
	{// sort entry increasing order of min clues
		for (int i = 0; i < 3; i++) {
			tsort[i] = i + (t416n6[bax[i].i416] << 8);
			tsort2[i] = i + 3 + (t416n6[bax[i + 3].i416] << 8);
		}
		for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++) {
			if (tsort[i] > tsort[j]) {
				int temp = tsort[i];	tsort[i] = tsort[j]; tsort[j] = temp;
			}
			if (tsort2[i] > tsort2[j]) {
				int temp = tsort2[i];	tsort2[i] = tsort2[j]; tsort2[j] = temp;
			}
		}
	}
	uint64_t mi1 = t416n6[bax[tsort[0] & 7].i416]; mi1 *= t416n6[bax[tsort[1] & 7].i416];
	uint64_t mi2 = t416n6[bax[tsort2[0] & 7].i416]; mi2 *= t416n6[bax[tsort2[1] & 7].i416];
	// setup map81
	for (int ibs = 0; ibs < 3; ibs++) {
		STD_B416& b = bax[ibs], & s = bax[ibs + 3];
		for (int c = 0; c < 27; c++) {
			b.map81[c] = c + 27 * ibs;
			s.map81[c] = C_transpose_d[c + 27 * ibs];
		}
	}

	// put bands in the right order  012
	for (int i = 0; i < 3; i++) {
		if (mi1 < mi2) {
			baxs[i] = bax[tsort[i] & 7];
			baxs[3 + i] = bax[tsort2[i] & 7];
		}
		else {
			baxs[i] = bax[tsort2[i] & 7];
			baxs[3 + i] = bax[tsort[i] & 7];
		}

	}
	int tppp[6][3] = { {3,4,5},{4,5,3},{5,3,4} ,
		{1,2,0},{2,0,1},{0,1,2} };
	gchk.mincluesb3 = 7;// first 7 and more in b3
	int irs = 0;
	for (gchk.iperm = 0; gchk.iperm < 6; gchk.iperm++) {// 6 triplets
		int* ptp = tppp[gchk.iperm];
		for (int i = 0; i < 3; i++)	bax[i] = baxs[ptp[i]];
		if (gchk.iperm == 5)gchk.mincluesb3 = 6;
		irs = gchk.StartIs18();
		if (first & irs) return 1;
	}
	if (gchk.nok != gchk.n18seen) cout << "redundancy to clear" << endl;
	return irs;

}

//========================================
const char * libs_c17_00_cpt2g[70] = {
	"0  compte entrée",//0
	"1  uas b12 start ",//1
	"2  uas g2 start",//2
	"3  steps first 3/4 ",//3
	"4  entry clean",
	"5  raw b12 ",
	"6  try b12 full",
	"7  try b12  not full ",
	"8  gob3 ",
	"9  all band 3 calls",
	"10 not critical 0",//10
	"11 not critical 1",//11
	"12 not critical 2",//12
	"13 b3 direct",//13

	"14 b3 expand",//14
	"15 isvalid b12",//15
	"16 clean not 6 b2 forbidden",//16
	"17 add after valid",//17
	"18 clean more uas ",
	"19  final checkb3",
	"20  call checkb3",
	"21 ",
	"22  max c2",
	"23  max cmore",
	"24 ",//24
	"25 entry add b1 b2 ",//25
	"26 cells added b1 b2 ",//26
	"27 add after g2",//27
	"28 add after g3",//28
	"29 ",//29

	"30  new ua b12",	
	"31 add to <18",	
	"32 add to <20",	
	"33 add to <23",	
	"34  ",	
	
	"35 add c2 ",	
	"36 ad other ",	
	"37 filter one of 3 perms",	"38",	"39",

	"40 ",	
	"41 limspot",	
	"42 ", 
	"43 d",	
	"44","45","46","47","48",	
	"49",
	"50 below -1",
	"51 below -1 kill",
	"52 below -1 assign",	
	"53 ",	"54 ",	"55 ",
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
	gchk.diagbugchk = gchk.aigstop = 0;
	char ze_diag[164]; ze_diag[163] = 0;
	memcpy(gchk.zes, ze, 164);
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
	// put bands in the right order  012
	baxs[0] = bax[tsort[0] & 7];
	baxs[1] = bax[tsort[1] & 7];
	baxs[2] = bax[tsort[2] & 7];
	bax[0] = baxs[0]; bax[1] = baxs[1]; bax[2] = baxs[2];
	int irs=0;
	gchk.band_order[0] = tsort[0] & 7;
	gchk.band_order[1] = tsort[1] & 7;
	gchk.band_order[2] = tsort[2] & 7;
	gchk.mincluesb3 = 6;
	irs= gchk.StartIs18();
#ifndef COLOIN

	// second perm  exchange bands 2/3  021
	gchk.mincluesb3 = 7;
	gchk.band_order[1] = tsort[2] & 7; // old band 3
	gchk.band_order[2] = tsort[1] & 7;
	bax[0] = baxs[0]; bax[1] = baxs[2]; bax[2] = baxs[1];
	irs += gchk.StartIs18();
	// third perm    120
	gchk.band_order[0] = tsort[1] & 7;// old band 2
	gchk.band_order[2] = tsort[0] & 7; // old band 1
	bax[0] = baxs[1]; baxs[1] = baxs[2]; bax[2] = baxs[0];
	irs += gchk.StartIs18();
#endif
#ifdef TEST_ON
	cout << "print final stats" << endl;
	for (int i = 0; i < 70; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
#endif

	return irs;
}

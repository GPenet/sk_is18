

// standard first band (or unique band)

void STD_B416::SetGangster() {
	memset(gangster, 0, sizeof gangster);
	for (int ir = 0, i = 0; ir < 3; ir++)
		for (int ic = 0; ic < 9; ic++, i++)
			gangster[ic] |= 1 << band0[i];
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = band0[cell];
			fd_sols[1][dig] |= Zhoucol << i; // add candidates in the column
			fd_sols[0][dig] |= 1 << cell;
		}
	}


}

void STD_B416::MorphUas() {
	// morph all uas
	for (uint32_t i = 0; i < nua; i++) {
		register int uao = tua[i]&BIT_SET_27, ua = 0;
		register uint32_t cc;
		while (bitscanforward(cc, uao)) {
			uao ^= 1 << cc;
			ua |= 1 << map[cc];
		}
		tua[i] = ua | _popcnt32(ua) << 27;// store with count
	}
}


void STD_B416::PrintStatus() {
	cout << "band status i=" << i416 << "\tstart=" << dband<<endl<<"map ";
	for (int i = 0; i < 27; i++)cout << map[i] << " ";
	cout <<endl;
	cout << band << endl<<"gangster status"<<endl;;
	cout << "UAs table" << endl;
	for (uint32_t i = 0; i < nua; i++)
		cout << Char27out(tua[i]) << endl;
}
void STD_B1_2::FillMiniDigsMiniPairs(STD_B1_2 & bb) {
	if (0) {
		cout << "gangster other band" << oct << endl;
		for (int i = 0; i < 9; i++)cout << bb.gangster[i] << " ";
		cout << endl << "my gangster  " << endl;
		for (int i = 0; i < 9; i++)cout << gangster[i] << " ";
		cout << dec << endl;
		cout << band << " my band" << endl;
	}

	nvpairs = 0;
	for (int i = 0, j = 0; i < 9; i++, j += 3) {
		int a = (1 << band0[j]), b = (1 << band0[j + 1]), c = (1 << band0[j + 2]);
		mini_digs[i] = a | b | c;
		mini_pairs[j] = b | c;// missing a  relative columns 2,3
		mini_pairs[j+1] = a | c;// missing b
		mini_pairs[j+2] = a | b;// missing c 
		int jcol = j % 9;// start col for the mini row
		int * gg = bb.gangster;
		if ((gg[jcol + 1] & c) && (gg[jcol + 2] & b))
			tv_pairs[nvpairs++] = j;
		if ((gg[jcol ] & c) && (gg[jcol + 2] & a))
			tv_pairs[nvpairs++] = j+1;
		if ((gg[jcol] & b) && (gg[jcol + 1] & a))
			tv_pairs[nvpairs++] = j + 2;
	}
	if (0) {
		cout << " valid pairs table ";
		for (int i = 0; i < nvpairs; i++)
			cout << tv_pairs[i] << " ";
		cout << endl;
		cout << " valid pairs pairs " << oct;
		for (int i = 0; i < nvpairs; i++)
			cout << mini_pairs[tv_pairs[i]] << " ";
		cout << dec << endl;
	}
}
int STD_B1_2::ReviseG_triplet(int imini, int ip, STD_B1_2 * bb) {
	int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
	int dcell=3*imini,dcol= dcell % 9,
		*cols = &bb->gangster[dcol],
		*myp= tp3f[ip];
	int digit[3], digit_bit[3], pdigit[3], pdigit_bit[3], digp = 0;
	for (int i = 0; i < 3; i++) {// collect digits
		digit[i] = band0[dcell + i];
		int bit = 1 << digit[i];
		digit_bit[i] = bit;
		digp |= bit;
	}
	for (int i = 0; i < 3; i++) {// collect digits perm
		pdigit[i] = digit[myp[i]];
		pdigit_bit[i] = 1 << pdigit[i];
	}
	int *g12 = &zh2b_g.gangster[dcol];
	if (!(pdigit_bit[0] & g12[0]) ||
		!(pdigit_bit[1] & g12[1]) ||
		!(pdigit_bit[2] & g12[2])) return 0;
	for (int ic = 0; ic < 3; ic++)
		bb->revised_g[dcol + ic] ^= (pdigit_bit[ic] | digit_bit[ic]);
	return digp;
}

uint32_t  STD_B1_2::GetMiniData(int index, uint32_t & bcells, STD_B1_2 *bb) {
	//index is cell 0-26 in the band assumed free in the mini-row
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	uint32_t tcells[3] = { 6,5,3 };// corresponding pairs
	int imini= index / 3,	dmini = 3*imini,dcol=dmini%9,
		 perm = index % 3,	*pcol = tpcol[perm];
	bcells = tcells[perm] << (3 * imini);
	uint32_t digs= mini_pairs[index];
	//bb->InitRevisedg();// must be done by the caller
	bb->revised_g[dcol + pcol[0]] ^= digs;
	bb->revised_g[dcol + pcol[1]] ^= digs;
	return digs;
}


void STD_B1_2::PrintShortStatus() {
	cout  << band << "\t\tband main data i0-415="<<i416 << endl;
	cout << "nua    \t" << nua << endl;

}

struct Ntable{
	int A;
	int C;
	int G;
	int U;
public:
	Ntable(){
		A=0;
		C=0;
		G=0;
		U=0;
	}
};

struct Ctable{
	int AU;
	int GC;
	int GU;
public:
	Ctable(){
		AU=0;
		GC=0;
		GU=0;
	}
};


inline void addNtable(Ntable &N, char n){
	switch (n)
	{
	case 'A':
		N.A++;
		break;
	case 'C':
		N.C++;
		break;
	case 'G':
		N.G++;
		break;
	case 'U':
		N.U++;
		break;
	default:
		cerr << "Unexpected nucleotide"  << n << endl;
	}
}
inline void subtNtable(Ntable &N, char n){
	switch (n)
	{
	case 'A':
		N.A--;
		break;
	case 'C':
		N.C--;
		break;
	case 'G':
		N.G--;
		break;
	case 'U':
		N.U--;
		break;
	default:
		cerr << "Unexpected nucleotide"  << n << endl;
	}
}

inline void addCtable(Ctable &C, char n1, char n2){
	if((n1 == 'A' && n2 == 'U') || (n1 == 'U' && n2 == 'A')){
		C.AU++;
	}
	else if((n1 == 'G' && n2 == 'C') || (n1 == 'C' && n2 == 'G')){
		C.GC++;
	}
	else if((n1 == 'G' && n2 == 'U') || (n1 == 'U' && n2 == 'G')){
		C.GU++;
	}
}

void showNtable(Ntable N){
	cout << "A=" << N.A << endl;
	cout << "C=" << N.C << endl;
	cout << "G=" << N.G << endl;
	cout << "U=" << N.U << endl;
}
void showCtable(Ctable C){
	cout << "AU=" << C.AU << endl;
	cout << "GC=" << C.GC << endl;
	cout << "GU=" << C.GU << endl;
}

float calcPseudoEnergy(const Ntable &N, const Ctable &C){
	float energy = 0;

	energy  = N.A * N.U * -1;
	energy += N.G * N.C * -3.12;
	energy += N.G * N.U * -1;

//	cout << "chk1:" << energy << endl;

	energy -= C.AU * -1;
//	cout << "chk2-1:" << energy << endl;
	energy -= C.GC * -3.12;
//	cout << "chk2-2:" << energy << ":" << C.GC << endl;
	energy -= C.GU * -1;
//	cout << "chk2-3:" << energy << endl;

	return energy;
}

string rev_fold_step1(const char *aaseq, const int aalen,
			codon &codon_table, const string &exc_codons){
	int nuc_len = aalen * 3 + 1;

	string optseq_r;
	optseq_r.resize(nuc_len, 'N');
	optseq_r[0] = ' '; // 1-based

	Ntable Ntab; //= {0,0,0,0};
	Ctable Ctab; // = {0,0,0};
	//Ntab.A = 0; Ntab.C = 0; Ntab.G = 0; Ntab.U = 0;
	//Ctab.AU = 0;Ctab.GC = 0;Ctab.GU = 0;

	//最初のコドンはランダムに選ぶ。
	InitRand();
	vector<string> codons1 = codon_table.getCodons(aaseq[0], exc_codons);
	shuffleStr(&codons1, codons1.size());

	optseq_r[1] = codons1[0][0];	addNtable(Ntab, optseq_r[1]);
	optseq_r[2] = codons1[0][1];	addNtable(Ntab, optseq_r[2]);
	optseq_r[3] = codons1[0][2];	addNtable(Ntab, optseq_r[3]);
	addCtable(Ctab, optseq_r[1], optseq_r[2]);
	addCtable(Ctab, optseq_r[2], optseq_r[3]);
	addCtable(Ctab, optseq_r[1], optseq_r[3]);

	//2個目以降のコドン
	for(int i = 1; i < aalen; i++){
		//cout << "ENTER" << endl;
		float maxP = -INF;
		string maxPcodon = "";
		Ntable maxN;
		Ctable maxC;

		vector<string> cand_codons = codon_table.getCodons(aaseq[i], exc_codons);
		for(unsigned int l = 0; l < cand_codons.size(); l++){
			string codon = cand_codons[l];
			Ntable N = Ntab;
			Ctable C = Ctab;
			//			Ntable N = {0,0,0,0};
			//			Ctable C = {0,0,0};
//			N = Ntab;
//			C = Ctab;

			for(int j = 0; j <=2; j++){
				int nuc_pos = i * 3 + j + 1;
				optseq_r[nuc_pos] = codon[j];
				addNtable(N, codon[j]);
				for(int k = MAX2(nuc_pos - 3, 1); k < nuc_pos - 1; k++){
					addCtable(C, optseq_r[k], codon[j]);
				}
			}

			float P = calcPseudoEnergy(N, C);
			//cout << Pe << endl;
			if(P > maxP){
				maxP = P;
				maxPcodon = codon;
				maxN = N;
				maxC = C;
			}
		}
		Ntab = maxN;
		Ctab = maxC;
		optseq_r[i*3+1] = maxPcodon[0];
		optseq_r[i*3+2] = maxPcodon[1];
		optseq_r[i*3+3] = maxPcodon[2];

	}

	//showNtable(Ntab);
	//showCtable(Ctab);

	cout << "step1:" << optseq_r << endl;

	//cout << "ok" << endl;
	return optseq_r;
}

Ntable countNtable(string &seq, int F){
	Ntable N;
	N.A = 0; N.C = 0; N.G = 0; N.U = 0;
	for(unsigned int i = F; i < seq.size(); i++){
		switch (seq[i])
		{
		case 'A':
			N.A++;
			break;
		case 'C':
			N.C++;
			break;
		case 'G':
			N.G++;
			break;
		case 'U':
			N.U++;
			break;
		default:
			cerr << "Unexpected nucleotide"  << seq[i] << endl;
		}
	}
	return N;
}
Ctable countCtable(string &seq, int F){
	Ctable C;
	C.AU = 0; C.GC = 0; C.GU = 0;

    for(unsigned int i = F; i < seq.size() - 1; i++){
        for(unsigned int j = i+1; j <= MIN2(i+3, seq.size()-1); j++){
        	//cout << i << ";" << j << endl;
        	addCtable(C, seq[i], seq[j]);
        }
    }

    return C;

}

void rev_fold_step2(string *optseq_r, const char *aaseq, const int aalen,
		codon &codon_table, const string &exc_codons){

	Ntable Ntab = countNtable(*optseq_r, 1);
	showNtable(Ntab);
	Ctable Ctab = countCtable((*optseq_r), 1);
	showCtable(Ctab);

	float max_energy_prev = -INF;
	float max_energy = calcPseudoEnergy(Ntab, Ctab);
	cout << "step2: " << max_energy << endl;
	string max_codon = "   ";
	Ntable max_N;
	Ctable max_C;
	int max_i = 0;
	string max_codon_from = "   ";

	int MAX_CYCLE = aalen * 3;
	int cycle = 0;
	while(cycle <= MAX_CYCLE){
		for(int i = 0; i < aalen; i++){
			//cout << aalen << endl;
			vector<string> codons = codon_table.getCodons(aaseq[i], exc_codons);
			if(codons.size() == 1) continue; // コドンが一つしか無いところは変異しない。

			int codon_fm = i * 3 + 1;
			int codon_to = i * 3 + 3;
			int local_fm = MAX2(1,codon_fm - 3);
			int local_to = MIN2((*optseq_r).size()-1, codon_to + 3);

//			cout << "TEST fm=" << local_fm << ":" << local_to << endl;

			// 変異前情報の整理
			string local_region_org;
			for(int j = local_fm; j <= local_to; j++){
				local_region_org.push_back((*optseq_r)[j]);
			}

			string codon_org;
			codon_org.push_back((*optseq_r)[codon_fm]);
			codon_org.push_back((*optseq_r)[codon_fm+1]);
			codon_org.push_back((*optseq_r)[codon_fm+2]);

//					if(i == 1){
//						cout << codon_org << endl;
//					}

			Ctable C_local_org = countCtable(local_region_org, 0);

			for(unsigned int m = 0; m < codons.size(); m++){
				string codon = codons[m];
				// 変異後情報の整理
				string local_region = local_region_org;
				// コドン部分の上書き
				if(i == 0){ // 最初のコドンを置換する
					local_region[0] = codon[0];
					local_region[1] = codon[1];
					local_region[2] = codon[2];
				}
				else{
					local_region[3] = codon[0];
					local_region[4] = codon[1];
					local_region[5] = codon[2];
				}

				Ctable C;
				Ctable C_local = countCtable(local_region, 0);
				C.AU = Ctab.AU + C_local.AU - C_local_org.AU;
				C.GC = Ctab.GC + C_local.GC - C_local_org.GC;
				C.GU = Ctab.GU + C_local.GU - C_local_org.GU;

				Ntable N = Ntab;
				for(int p = 0; p < 3; p++){
					char n_org = codon_org[p];
					char n= codon[p];
					if(n_org != n){
						addNtable(N, n);
						subtNtable(N, n_org);
					}
				}

				float en = calcPseudoEnergy(N, C);
				//cout << i << ":" << en << "->" << max_energy << endl;
				if(en >= max_energy){
					max_energy = en;
					max_codon = codon;
					max_N = N;
					max_C = C;
					max_i = i;
					max_codon_from = codon_org;
				}
			}
		}
		cycle++;

		// 配列のアップデート
		int codon_fm = max_i * 3 + 1;
		int codon_to = max_i * 3 + 3;
		int j = 0;
		for(int i = codon_fm; i <= codon_to; i++){
			//cout << "ok:" << (*optseq_r)[i]  << "->" << max_codon[j] << endl;
			(*optseq_r)[i] = max_codon[j++];
		}
		Ntab = max_N;
		Ctab = max_C;
		cout << "pos = " << max_i << "," << max_codon_from << "->" << max_codon << endl;
		cout << (*optseq_r) << "\t" << max_energy << endl;

//		showNtable(Ntab);
// 		showCtable(Ctab);

		if(max_energy == max_energy_prev) break;
		max_energy_prev = max_energy;
	}

}

/*
void rev_fold_step2_bk(string *optseq_r, const char *aaseq, const int aalen,
		codon &codon_table, const string &exc_codons, const int *ofm, const int *oto, const int n_interval){

	Ntable Ntab = countNtable(*optseq_r, 1);
	showNtable(Ntab);
	Ctable Ctab = countCtable((*optseq_r), 1);
	showCtable(Ctab);

	float max_energy_prev = -INF;
	float max_energy = calcPseudoEnergy(Ntab, Ctab);
	cout << "step2: " << max_energy << endl;
	string max_codon = "   ";
	Ntable max_N;
	Ctable max_C;
	int max_i = 0;
	string max_codon_from = "   ";

	int MAX_CYCLE = aalen * 3;
	int cycle = 0;
	while(cycle <= MAX_CYCLE){
		for(int i = 0; i < aalen; i++){
			for(int l = 0; l < n_interval; l++){
				if(ofm[l]/3 <= i+1 && i+1 <= oto[l]/3){ //インターバルに入っているので、最適化する。
					//cout << aalen << endl;
					vector<string> codons = codon_table.getCodons(aaseq[i], exc_codons);
					if(codons.size() == 1) break; // コドンが一つしか無いところは変異しない。
					int codon_fm = i * 3 + 1;
					int codon_to = i * 3 + 3;
					int local_fm = MAX2(1,codon_fm - 3);
					int local_to = MIN2((*optseq_r).size()-1, codon_to + 3);

					// 変異前情報の整理
					string local_region_org;
					for(int j = local_fm; j <= local_to; j++){
			            local_region_org.push_back((*optseq_r)[j]);
			        }

					string codon_org;
					codon_org.push_back((*optseq_r)[codon_fm]);
					codon_org.push_back((*optseq_r)[codon_fm+1]);
					codon_org.push_back((*optseq_r)[codon_fm+2]);

//					if(i == 1){
//						cout << codon_org << endl;
//					}

					Ctable C_local_org = countCtable(local_region_org, 0);

					for(unsigned int m = 0; m < codons.size(); m++){
						string codon = codons[m];
						// 変異後情報の整理
						string local_region = local_region_org;
						local_region[3] = codon[0];
						local_region[4] = codon[1];
						local_region[5] = codon[2];

						Ctable C;
						Ctable C_local = countCtable(local_region, 0);
						C.AU = Ctab.AU + C_local.AU - C_local_org.AU;
						C.GC = Ctab.GC + C_local.GC - C_local_org.GC;
						C.GU = Ctab.GU + C_local.GU - C_local_org.GU;

						Ntable N = Ntab;
						for(int p = 0; p < 3; p++){
							char n_org = codon_org[p];
							char n= codon[p];
							if(n_org != n){
								addNtable(N, n);
								subtNtable(N, n_org);
							}
						}

						float en = calcPseudoEnergy(N, C);
						//cout << i << ":" << en << "->" << max_energy << endl;
						if(en >= max_energy){
							max_energy = en;
							max_codon = codon;
							max_N = N;
							max_C = C;
							max_i = i;
							max_codon_from = codon_org;
						}
					}
				}
			}
		}
		cycle++;

		// 配列のアップデート
		int codon_fm = max_i * 3 + 1;
		int codon_to = max_i * 3 + 3;
		int j = 0;
		for(int i = codon_fm; i <= codon_to; i++){
			//cout << "ok:" << (*optseq_r)[i]  << "->" << max_codon[j] << endl;
			(*optseq_r)[i] = max_codon[j++];
		}
		Ntab = max_N;
		Ctab = max_C;
		cout << "pos = " << max_i << "," << max_codon_from << "->" << max_codon << endl;
		cout << (*optseq_r) << "\t" << max_energy << endl;

//		showNtable(Ntab);
// 		showCtable(Ctab);

		if(max_energy == max_energy_prev) break;
		max_energy_prev = max_energy;
	}

}
*/

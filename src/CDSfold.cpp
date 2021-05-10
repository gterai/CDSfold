/*
 * CDSfold.cpp

 *
 *  Created on: Sep 2, 2014
 *      Author: terai
 */
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))
#define TURN 3
#include <cstdio>
#include <iostream>
#include <sstream>
#include <time.h>
#include <unistd.h>
#include <string>

extern "C" {
#include  "utils.h"
#include  "fold_vars.h"
#include  "fold.h"
#include  "part_func.h"
#include  "params.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "ctype.h"
#include "limits.h"
}

#include "codon.hpp"
#include "fasta.hpp"
#include "CDSfold.hpp"
#include "CDSfold_rev.hpp"
#include "AASeqConverter.hpp"
//#include <algorithm>
//#include <sys/time.h>
//#include <sys/resource.h>

int *indx;
int BP_pair[5][5] =
/* _  A  C  G  U  */
{ { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 5 }, { 0, 0, 0, 1, 0 }, { 0, 0, 2, 0, 3 }, {
		0, 6, 0, 4, 0 } };
int rtype[7] = { 0, 2, 1, 4, 3, 6, 5 };

//#define MAXLOOP 20
#define noGUclosure  0
int test;

using namespace std;
int main(int argc, char *argv[]) {
	//printf("%d\n%ld", INT_MAX, LONG_MAX);
	int W = 0; // -w
	string exc = ""; // -e
	int m_disp = 0; // -M
	int rand_tb_flg = 0; //-R
	int rev_flg = 0; // -r
	int part_opt_flg = 0; // -f and -t
	int opt_fm = 0; // -f
	int opt_to = 0; // -t
	// get options
	{
		int opt;
		while((opt=getopt(argc,argv,"w:e:f:t:rMR"))!=-1){
			switch(opt){
			case 'w':
				W = atoi(optarg);
				break;
			case 'e':
				exc = string(optarg);
				break;
			case 'M':
				m_disp = 1;
				break;
			case 'R':
				rand_tb_flg = 1;
				break;
			case 'r':
				rev_flg = 1;
				break;
			case 'f':
				opt_fm = atoi(optarg);
				part_opt_flg = 1;
				break;
			case 't':
				opt_to = atoi(optarg);
				part_opt_flg = 1;
				break;

			}
		}
	}
	//exit(0);


	// -R オプションに関するチェック
	if(rand_tb_flg){
		if(W!=0 || exc !="" || m_disp || rev_flg || part_opt_flg){
			cerr << "The -R option must not be used together with other options." << endl;
			return 0;
		}
	}


	//int energy = E(5, 1, 1, 2, "ATGCATGC");
	//int energy = E_Hairpin(5, 1, 1, 2, "ATGCATGC");

	map<char, int> n2i = make_n2i();
	char i2n[20];
	make_i2n(i2n);

	int i2r[20], ii2r[100];
	make_i2r(i2r);
	make_ii2r(ii2r);

	AASeqConverter conv;

	codon codon_table;
	//codon_table.Table();

	const char dummy_str[10] = "XXXXXXXXX";
	//int NCflg = 1;
//	int TB_CHK_flg = 0;
//	int preHPN_flg = 1;
	int TEST = 1;
	int DEPflg = 1;
	int NCflg = 0;
	//	char *NucDef = "*AUGGGUCUUCCAGUGUCAUUACGAGCUGACACCAUUCGAGAUUUAUUACUUGGUGUCAGCUCGAUAAUGACCUGGAAGACCCUUGCUCUUGUGUUAGCUGUGAUCAAUCUCAAGAAUCUGCCACUAGUGUGGCACCCGGGGGAUCCUCAUUUCCCCCGGGGGAAGGCGCUGGUGACGCAUACGGGCAAACCCACUCAUCCGGUGUUUGUCCCGUAUGCGAUCACCAGUCGCACUCCGAUUCUUGAGACUGAUUACAACUUUCACAAGAGCAAUUCCACGUAUUUUAGCGAUUUGGAUAUU";
	const char NucDef[] = "*AUGGAGGGGAUUGUCACGGGAGAUCGGCUUGCUUGCGUGGCGCUUCAUGGAAGCUCUUUGCUCCAUGAAGCGUCCGUAAGCAAGUAUACCGAUAUCCCGGGCAUUCUCCUCCAAUACAUCGAUGAAUUUCCCCUCACUGAUAUUGCCGCGCACGCGCCACGCGAGGCGUGGCAAAGCCUGUGCGAACAGGCGAUCUGUAUCGUCCAUCAUAUUAGCGACCGGGGCAUCCUCAAUGAGGAUGUUAAAACCCGGUCGCUGACGAUACAGAUCAACAGUGAGGGGAUGUUCAAGAUGUUUAUG";
	//string tmp_def = "*AUGGCCCCCAUACAGCAGAAGGCACUAAUCAACUGCGAUAUGGGGGAAGCUUACGGGAACUGGGCCUGCGGCCCAGAUCUCGAGCUCCUCCCCAUGAUCGACAUCGCCAACGUGGCGUGUGGAUUUCAUGGGGGGGAUCCAUUAAUAAUGAUGGAAACGGUGCGCAACUGUAAAGCGCACAAUGUGCGCAUAGGGGCGCACCCUGGCCUCCCGGACCUGCAGGGGUUCGGGAGGCGGGAGAUGAAACUCUCCCCUGAAGAGCUCACCGCCAUGACUAUUUAUCAGGUGGGAGCUCUUCAG";
	//char *NucDef = "*AUGAGUCUGGCGUGCAUGGCCAAGUAG";
	//char *NucDef = "*AUGUCUCUCGCGUGCAUGGCCAAGUGA";
	//char *NucDef = "*AUGUCUUUAGCCUGUAUGGCUAAAUAA";
	//cout << "optind is " << optind << endl;
	//const char *NucDef = tmp_def.c_str();
	fasta all_aaseq(argv[optind]); // get all sequences

	cout << "W = " << W << endl;
	cout << "e = " << exc << endl;
	do {
		char *aaseq = all_aaseq.getSeq();
		int aalen = all_aaseq.getSeqLen();

		if(aalen <= 2){
			cerr << "The amino acid sequence is too short.\n";
			exit(1);
		}

		int n_inter = 0; //今の実装では、n_inter=1 or 2となる。
		int ofm[100];
		int oto[100];
		if(part_opt_flg){
			// 部分最適化が指定された。
			if(opt_fm == 0 ||opt_to == 0){
				cerr << "The -f and -t option must be used together." << endl;
				exit(1);
			}
			if(opt_fm < 1){
				cerr << "The -f value must be 1 or more." << endl;
				exit(1);
			}
			if(opt_to < 1){
				cerr << "The -t value must be 1 or more." << endl;
				exit(1);
			}
			if(opt_to < opt_fm){
				cerr << "The -f value must be smaller than -t value." << endl;
				exit(1);
			}
			if(opt_to > aalen){
				opt_to = aalen;
			}

			// 部分逆最適化情報の作成
			if(rev_flg){ // 指定された領域の構造除去
				ofm[0] = (opt_fm-1) * 3 + 1;
				oto[0] = opt_to * 3;
				n_inter = 1;
			}
			else{ // 指定された領域の構造安定化
				int l = 0;
				if(opt_fm != 1){
					ofm[l] = 1;
					oto[l++] = (opt_fm-1) * 3;
				}
				if(opt_to != aalen){
					ofm[l] = opt_to * 3 + 1;
					oto[l++] = aalen * 3;
				}
				n_inter = l;
			}
			// 最適化領域のチェック
			//for(int I = 0; I < n_inter; I++){
			//	cout << ofm[I] << "-" << oto[I] << endl;
			//}
			//exit(0);
		}

		int nuclen = aalen * 3;
		int w_tmp;

		if(W == 0){
			w_tmp = nuclen;
		}
		else if(W < 10){
			cerr << "W must be more than 10" << "(you used " << W << ")" <<endl;
			exit(1);
		}
		else if(W > nuclen){
			w_tmp = nuclen;
		}
		else{
			w_tmp = W;
		}

		//		w_tmp = 50;// test!
//		vector<vector<vector<string> > >  substr = conv.getBases(string(aaseq),8, exc);
		vector<vector<vector<string> > >  substr = conv.getOriginalBases(string(aaseq), exc);
		vector<vector<int> > Dep1;
		vector<vector<int> > Dep2;

		Dep1 = conv.countNeighborTwoBase(string(aaseq), exc);
		Dep2 = conv.countEveryOtherTwoBase(string(aaseq), exc);


		//		cout << ptotal_Mb_base << endl;

//		pid_t pid2 = getpid();
//		stringstream ss2;
//		ss2 << "/proc/" << pid2 << "/status";
//		int m2 = getMemoryUsage(ss2.str());
//		cout << "Memory(VmRSS): "  << float(m2)/1024 << " Mb" << endl;
		//cout << "Estimate: "  << ptotal_Mb_base << " Mb" << endl;
		//exit(0);

		//map<string, int> predefHPN_E;
		map<string, int> predefHPN_E = conv.getBaseEnergy();
		//		vector<vector<vector<vector<pair<int, string> > > > > predefHPN = conv.calcQueryOriginalBaseEnergy(string(aaseq), "");
		vector<vector<vector<vector<pair<int, string> > > > > predefHPN;
		//vector<vector<vector<vector<pair<int, string> > > > > predefHPN = conv.calcQueryOriginalBaseEnergy(string(aaseq), "");
/*
		for(unsigned int i = 1; i < predefHPN.size(); i++){
 			cout << ">>position " << i << endl;
			for(unsigned int l = 0; l < predefHPN[i].size(); l++){
				for(unsigned int li = 0; li < predefHPN[i][l].size(); li++){
					for(unsigned int rj = 0; rj < predefHPN[i][l][li].size(); rj++){
						if(predefHPN[i][l][li][rj].first != INF)
							continue;
						cout << predefHPN[i][l][li][rj].second << ":" << predefHPN[i][l][li][rj].first << endl;
					}
				}
			}
 		}
 		cout << INF << endl;
		exit(0);
*/

		stack sector[500];

		//vector<int> NucConst = createNucConstraint(NucDef, nuclen, n2i);

		vector<int> NucConst;
		if(NCflg){
			NucConst = createNucConstraint(NucDef, nuclen, n2i);
		}
		//		for(int i = 1; i <= nuclen; i++)
		//		printf("%d %d\n", i, NucConst[i]);

		//createNucConstraint

		cout << aaseq << endl;
//		cout << aalen << endl;

		vector<vector<int> > pos2nuc = getPossibleNucleotide(aaseq, aalen, codon_table, n2i, exc);
//		vector<vector<int> > pos2nuc = getPossibleNucleotide(aaseq, aalen, codon_table, n2i, 'R');
//		showPos2Nuc(pos2nuc, i2n);
//		exit(0);
		indx = new int[nuclen + 1];

		set_ij_indx(indx, nuclen, w_tmp);
		//set_ij_indx(indx, nuclen);

		string optseq;
		optseq.resize(nuclen+1, 'N');
		optseq[0] = ' ';

		string optseq_org;
		optseq_org.resize(nuclen+1, 'N');
		optseq_org[0] = ' ';


		//	  int ***C, ***Mbl, ***Mbr, ***Mbb, ***M, ***F, ***Fbr, ***tFbr;
		int ***C, ***M, ***F;
		int ***F2;
		int ***DMl,***DMl1,***DMl2;
		int *chkC, *chkM;
		bond *base_pair;

			//		int n_inter = 1;


		paramT *P = NULL;
		P = scale_parameters();
		update_fold_params();

//		rev_flg = 0;
//		if(rev_flg && num_interval == 0){
		if(rev_flg && !part_opt_flg){
			// reverse mode
			string optseq_rev = rev_fold_step1(aaseq, aalen, codon_table, exc);
			//			rev_fold_step2(&optseq_rev, aaseq, aalen, codon_table, exc, ofm, oto, 1);
			rev_fold_step2(&optseq_rev, aaseq, aalen, codon_table, exc);
			fixed_fold(optseq_rev, indx, w_tmp, predefHPN_E, BP_pair, P, aaseq, codon_table);
			free(P);
			break; //returnすると、実行時間が表示されなくなるためbreakすること。
		}





		//		allocate_arrays(nuclen, indx, pos2nuc, pos2nuc, &C, &M, &F);
		allocate_arrays(nuclen, indx, w_tmp, pos2nuc, &C, &M, &F, &DMl, &DMl1, &DMl2, &chkC, &chkM, &base_pair);
		if(rand_tb_flg){
			allocate_F2(nuclen, indx, w_tmp, pos2nuc, &F2);
		}
		//float ptotal_Mb = ptotal_Mb_alloc + ptotal_Mb_base;


//		pid_t pid1 = getpid();
//		stringstream ss1;
//		ss1 << "/proc/" << pid1 << "/status";
//		int m1 = getMemoryUsage(ss1.str());
//		cout << "Memory(VmRSS): "  << float(m1)/1024 << " Mb" << endl;
//		exit(0);

		// main routine
		for (int l = 2; l <= 4; l++) {
			for (int i = 1; i <= nuclen - l + 1; i++) {
//				test = 1;
				int j = i + l - 1;
				//int ij = indx[j] + i;
				int ij = getIndx(i, j, w_tmp, indx);

				chkC[ij] = INF;
				chkM[ij] = INF;

				for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
					int L_nuc = pos2nuc[i][L];
					if(NCflg == 1 && i2r[L_nuc] != NucConst[i]){	continue;}
					for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {
						int R_nuc = pos2nuc[j][R];
						if(NCflg == 1 && i2r[R_nuc] != NucConst[j]){	continue;}
						//L-R pair must be filtered
//						if(j-1==1){
//							cout << i << ":" << j << " " << L_nuc << "-" << R_nuc << " " << Dep1[ii2r[L_nuc*10+R_nuc]][i] <<endl;
//						}
						if(DEPflg && j-i == 1 && i <= nuclen - 1 && Dep1[ii2r[L_nuc*10+R_nuc]][i] == 0){continue;} // nuclen - 1はいらないのでは？
						if(DEPflg && j-i == 2 && i <= nuclen - 2 && Dep2[ii2r[L_nuc*10+R_nuc]][i] == 0){continue;}

						C[ij][L][R] = INF;
						M[ij][L][R] = INF;
						if(rand_tb_flg)
							F2[ij][L][R] = 0;

					}
				}
			}
		}

		//		cout << "TEST" << M[13][0][0] << endl;
		// main routine
		for (int l = 5; l <= nuclen; l++) {
			if(l > w_tmp) break;
			cout << "process:" << l << endl;

			//	  for(int l = 5; l <= 5; l++){
			for (int i = 1; i <= nuclen - l + 1; i++) {
				int j = i + l - 1;

				int opt_flg_ij = 1;
				if(part_opt_flg){
					for(int I = 0; I < n_inter; I++){
						if((ofm[I] <= i && oto[I] >= i) ||
								(ofm[I] <= j && oto[I] >= j)){
							opt_flg_ij = 0;
							break;
						}
					}
				}


				for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
					int L_nuc = pos2nuc[i][L];
//					cout << NCflg << endl;
					if(NCflg == 1 && i2r[L_nuc] != NucConst[i]){	continue;}
//					cout << "ok" << endl;

					for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {

						int R_nuc = pos2nuc[j][R];

						if(NCflg == 1 && i2r[R_nuc] != NucConst[j]){	continue;}

						//int ij = indx[j] + i;
						int ij = getIndx(i,j,w_tmp,indx);

						C[ij][L][R] = INF;
						M[ij][L][R] = INF;
						//						cout << i << " " << j << ":" << M[ij][L][R] << endl;

						int type = BP_pair[i2r[L_nuc]][i2r[R_nuc]];


						if (type && opt_flg_ij) {
							// hairpin
							if((l == 5 || l ==6 || l == 8) && TEST){
								for(unsigned int s = 0; s < substr[i][l].size(); s++){
									string hpn = substr[i][l][s];
									int hL_nuc  = n2i[hpn[0]];
									int hL2_nuc = n2i[hpn[1]];
									int hR2_nuc = n2i[hpn[l-2]];
									int hR_nuc  = n2i[hpn[l-1]];
									if(hL_nuc != i2r[L_nuc]) continue;
									if(hR_nuc != i2r[R_nuc]) continue;

									if(NCflg == 1){
										string s1 = string(NucDef).substr(i, l);
										if(hpn != s1) continue;
									}

//									cout << hpn << endl;
//
									if(DEPflg && L_nuc > 4 && Dep1[ii2r[L_nuc*10+hL2_nuc]][i] == 0){continue;}   // Dependencyをチェックした上でsubstringを求めているので, hpnの内部についてはチェックする必要はない。
									if(DEPflg && R_nuc > 4 && Dep1[ii2r[hR2_nuc*10+R_nuc]][j-1] == 0){continue;} // ただし、L_nuc、R_nucがVWXYのときだけは、一つ内側との依存関係をチェックする必要がある。
																											      // その逆に、一つ内側がVWXYのときはチェックの必要はない。既にチェックされているので。
									if(predefHPN_E.count(hpn) > 0){
										C[ij][L][R] = MIN2(predefHPN_E[hpn], C[ij][L][R]);

									}
									else{
//										int energy = HairpinE(j - i - 1, type,
//												i2r[hL2_nuc], i2r[hR2_nuc],
//												dummy_str);
										int energy = E_hairpin(j - i - 1, type,
												i2r[hL2_nuc], i2r[hR2_nuc],
												dummy_str, P);
										C[ij][L][R] = MIN2(energy, C[ij][L][R]);
									}
								}
								//exit(0);
							}
							else{
								for (unsigned int L2 = 0;
										L2 < pos2nuc[i + 1].size(); L2++) {
									int L2_nuc = pos2nuc[i + 1][L2];
									if(NCflg == 1 && i2r[L2_nuc] != NucConst[i+1]){	continue;}
									//if(chkDep2){continue:}
									for (unsigned int R2 = 0;
											R2 < pos2nuc[j - 1].size(); R2++) {
										int R2_nuc = pos2nuc[j - 1][R2];
										if(NCflg == 1 && i2r[R2_nuc] != NucConst[j-1]){	continue;}

										if(DEPflg && Dep1[ii2r[L_nuc*10+L2_nuc]][i] == 0){continue;}
										if(DEPflg && Dep1[ii2r[R2_nuc*10+R_nuc]][j-1] == 0){continue;}

										int energy;
										//cout << j-i-1 << ":" << type << ":" << i2r[L2_nuc] << ":" << i2r[R2_nuc] << ":" << dummy_str << endl;
										//										energy = HairpinE(j - i - 1, type,
										//												i2r[L2_nuc], i2r[R2_nuc],
										//												dummy_str);
										energy = E_hairpin(j - i - 1, type,
												i2r[L2_nuc], i2r[R2_nuc],
												dummy_str, P);
										//cout << "HairpinE(" << j-i-1 << "," << type << "," << i2r[L2_nuc] << "," << i2r[R2_nuc] << ")" << " at " << i << "," << j << ":" << energy << endl;
										//cout << i << " " << j  << " " << energy << ":" << i2n[L_nuc] << "-" << i2n[R_nuc] << "<-" << i2n[L2_nuc] << "-" << i2n[R2_nuc] << endl;
										C[ij][L][R] = MIN2(energy, C[ij][L][R]);

										// check predefined hairpin energy
										//if((l == 5 || l == 6 || l == 8) && preHPN_flg == 1){
										//  if(predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].second != ""){
										//		if(NCflg == 1){
										//			string s1 = string(NucDef).substr(i, l);
										//			if(predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].second == s1){
										//				//C[ij][L][R] = MIN2(C[ij][L][R], predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].first);
										//				C[ij][L][R] = predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].first; // Note that predefined hairpin is forced when it is found
										//			}
										//			}
										//		else{
										//			//一つ内側の塩基とのDependencyをチェックする。
										//			string s1 = predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].second;
										//			int preL2_nuc = n2i[s1[1]];
										//			int preR2_nuc = n2i[s1[s1.size()-2]];
//										//			cout << s1 << endl;
										//			if(DEPflg && Dep1[ii2r[L_nuc*10+preL2_nuc]][i] == 0){continue;}
										//			if(DEPflg && Dep1[ii2r[preR2_nuc*10+R_nuc]][j-1] == 0){continue;}
										//			C[ij][L][R] = predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].first; // Note that predefined hairpin is forced when it is found
										//		}
//										//	exit(0);
										//	}
										//}
									}

								}
							}

							// interior loop
							//cout << i+1 << " " <<  MIN2(j-2-TURN,i+MAXLOOP+1) << endl;
							for (int p = i + 1;
									p <= MIN2(j-2-TURN, i+MAXLOOP+1); p++) { // loop for position q, p
								int minq = j - i + p - MAXLOOP - 2;
								if (minq < p + 1 + TURN)
									minq = p + 1 + TURN;
								for (int q = minq; q < j; q++) {

									int pq = getIndx(p,q,w_tmp, indx);

									for (unsigned int Lp = 0;
											Lp < pos2nuc[p].size(); Lp++) {
										int Lp_nuc = pos2nuc[p][Lp];
										if(NCflg == 1 && i2r[Lp_nuc] != NucConst[p]){	continue;}

										if(DEPflg && p == i + 1 && Dep1[ii2r[L_nuc*10+Lp_nuc]][i] == 0){ continue;}
										if(DEPflg && p == i + 2 && Dep2[ii2r[L_nuc*10+Lp_nuc]][i] == 0){ continue;}


										for (unsigned int Rq = 0;
												Rq < pos2nuc[q].size(); Rq++) { // nucleotide for p, q
											int Rq_nuc = pos2nuc[q][Rq];
											if(NCflg == 1 && i2r[Rq_nuc] != NucConst[q]){	continue;}

											if(DEPflg && q == j - 1 && Dep1[ii2r[Rq_nuc*10+R_nuc]][q] == 0){ continue;}
											if(DEPflg && q == j - 2 && Dep2[ii2r[Rq_nuc*10+R_nuc]][q] == 0){ continue;}

											int type_2 =
													BP_pair[i2r[Lp_nuc]][i2r[Rq_nuc]];

											if (type_2 == 0)
												continue;
											type_2 = rtype[type_2];


//											if (noGUclosure)
//												if ((type_2 == 3)
//														|| (type_2 == 4))
//													if ((p > i + 1)
//															|| (q < j - 1))
//														continue; /* continue unless stack *//* no_close is removed. It is related with BONUS */

											//											if(i==8&&j==19){
//												cout << "test:" << p << "-" << q << endl;
//											}

											// for each intloops
											for (unsigned int L2 = 0;
													L2 < pos2nuc[i + 1].size();
													L2++) { // nucleotide for i+1,j-1
												int L2_nuc = pos2nuc[i + 1][L2];
												if(NCflg == 1 && i2r[L2_nuc] != NucConst[i+1]){	continue;}

												if(DEPflg && Dep1[ii2r[L_nuc*10+L2_nuc]][i] == 0){ continue;}


												for (unsigned int R2 = 0;
														R2
																< pos2nuc[j - 1].size();
														R2++) {
													int R2_nuc =
															pos2nuc[j - 1][R2];
													if(NCflg == 1 && i2r[R2_nuc] != NucConst[j-1]){	continue;}

													if(DEPflg && Dep1[ii2r[R2_nuc*10+R_nuc]][j-1] == 0){ continue;}

													for (unsigned int Lp2 = 0;
															Lp2
																	< pos2nuc[p
																			- 1].size();
															Lp2++) { // nucleotide for p-1,q+1
														int Lp2_nuc = pos2nuc[p
																- 1][Lp2];
														if(NCflg == 1 && i2r[Lp2_nuc] != NucConst[p-1]){ continue;}

														if(DEPflg && Dep1[ii2r[Lp2_nuc*10+Lp_nuc]][p-1] == 0){ continue;}
														if(p == i + 2 && L2_nuc != Lp2_nuc){ continue; } // check when a single nucleotide between i and p, this sentence confirm the dependency between Li_nuc and Lp2_nuc
														if(DEPflg && i + 3 == p && Dep1[ii2r[L2_nuc*10+Lp2_nuc]][i+1] == 0){ continue;} // check dependency between i+1, p-1 (i,X,X,p)

														for (unsigned int Rq2 =
																0;
																Rq2
																		< pos2nuc[q
																				+ 1].size();
																Rq2++) {
															int Rq2_nuc =
																	pos2nuc[q
																			+ 1][Rq2];
															if(q == j - 2 && R2_nuc != Rq2_nuc){ continue; } // check when a single nucleotide between q and j,this sentence confirm the dependency between Rj_nuc and Rq2_nuc

															if(NCflg == 1 && i2r[Rq2_nuc] != NucConst[q+1]){	continue;}

															if(DEPflg && Dep1[ii2r[Rq_nuc*10+Rq2_nuc]][q] == 0){ continue;}
															if(DEPflg && q + 3 == j && Dep1[ii2r[Rq2_nuc*10+R2_nuc]][q+1] == 0){ continue;} // check dependency between q+1, j-1 (q,X,X,j)

															int int_energy =
																	E_intloop(
																			p
																			- i
																			- 1,
																			j
																			- q
																			- 1,
																			type,
																			type_2,
																			i2r[L2_nuc],
																			i2r[R2_nuc],
																			i2r[Lp2_nuc],
																			i2r[Rq2_nuc],
																			P);
																	//LoopEnergy(p- i- 1,j- q- 1,type,type_2,i2r[L2_nuc],i2r[R2_nuc],i2r[Lp2_nuc],i2r[Rq2_nuc]);

															//int energy =
															//		int_energy
															//		+ C[indx[q]
															//			+ p][Lp][Rq];

															int energy =
																	int_energy
																	+ C[pq][Lp][Rq];
															C[ij][L][R] =
																	MIN2(energy,
																			C[ij][L][R]);

														}

													}
												}
											}
										}
									}
								} /* end q-loop */
							} /* end p-loop */

							// multi-loop
							for (unsigned int Li1 = 0;
									Li1 < pos2nuc[i + 1].size(); Li1++) {
								int Li1_nuc = pos2nuc[i+1][Li1];
								if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){	continue;}

								if(DEPflg && Dep1[ii2r[L_nuc*10+Li1_nuc]][i] == 0){ continue;}

								for (unsigned int Rj1 = 0;
										Rj1 < pos2nuc[j - 1].size(); Rj1++) {
									int Rj1_nuc = pos2nuc[j-1][Rj1];
									if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){	continue;}

									if(DEPflg && Dep1[ii2r[Rj1_nuc*10+R_nuc]][j-1] == 0){ continue;}
									//if(DEPflg && j-i == 2 && i <= nuclen - 2 && Dep2[ii2r[L_nuc*10+R_nuc]][i] == 0){continue;}
									if(DEPflg && (j-1)-(i+1) == 2 && Dep2[ii2r[Li1_nuc*10+Rj1_nuc]][i+1] == 0){continue;} // 2014/10/8 i-jが近いときは、MLclosingする必要はないのでは。少なくとも3つのステムが含まれなければならない。それには、５＋５＋２（ヘアピン2個分＋2塩基）の長さが必要。

									int energy = DMl2[i+1][Li1][Rj1]; // 長さが2個短いときの、複合マルチループ。i'=i+1を選ぶと、j'=(i+1)+(l-2)-1=i+l-2=j-1(because:j=i+l-1)
									int tt = rtype[type];

									energy += P->MLintern[tt];
									if(tt > 2)
										energy += P->TerminalAU;

									energy += P->MLclosing;
									//cout << "TEST:" << i << " " << j << " " << energy << endl;
									C[ij][L][R] =
											MIN2(energy,
													C[ij][L][R]);

//									if(C[ij][L][R] == -1130 && ij == 10091){
//										exit(0);
//									}


								}
							}


//							cout << "ok" << endl;
						}

						else C[ij][L][R] = INF;


						// fill M
						// create M[ij] from C[ij]
						if(type){
					        int energy_M = C[ij][L][R];
					        if(type > 2)
					          energy_M += P->TerminalAU;

					        energy_M += P->MLintern[type];
					        M[ij][L][R] = energy_M;
						}

						// create M[ij] from M[i+1][j]
						for (unsigned int Li1 = 0;
								Li1 < pos2nuc[i + 1].size(); Li1++) {
							int Li1_nuc = pos2nuc[i + 1][Li1];
							if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]){	continue;}
							if(DEPflg && Dep1[ii2r[L_nuc*10+Li1_nuc]][i] == 0){ continue;}

							//int energy_M = M[indx[j]+i+1][Li1][R]+P->MLbase;
							int energy_M = M[getIndx(i+1, j, w_tmp, indx)][Li1][R]+P->MLbase;
					        M[ij][L][R] = MIN2(energy_M, M[ij][L][R]);
						}

						// create M[ij] from M[i][j-1]
						for (unsigned int Rj1 = 0;
								Rj1 < pos2nuc[j - 1].size(); Rj1++) {
							int Rj1_nuc = pos2nuc[j - 1][Rj1];
							if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]){	continue;}
							if(DEPflg && Dep1[ii2r[Rj1_nuc*10+R_nuc]][j-1] == 0){ continue;}

							//int energy_M = M[indx[j-1]+i][L][Rj1]+P->MLbase;
							int energy_M = M[getIndx(i,j-1, w_tmp,indx)][L][Rj1]+P->MLbase;
					        M[ij][L][R] = MIN2(energy_M, M[ij][L][R]);
						}


						/* modular decomposition -------------------------------*/
						for (int k = i + 2 + TURN; k <= j - TURN - 1; k++) { // Is this correct?
							//cout << k << endl;
							for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size();
									Rk1++) {
								int Rk1_nuc = pos2nuc[k-1][Rk1];
								if(NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]){	continue;}
								//if(DEPflg && k == i + 2 && Dep1[ii2r[L_nuc*10+Rk1_nuc]][k-1] == 0){ continue;} // dependency between i and k - 1(=i+1)
								//if(DEPflg && k == i + 3 && Dep2[ii2r[L_nuc*10+Rk1_nuc]][k-1] == 0){ continue;} // dependency between i and k - 1(=i+2)

								for (unsigned int Lk = 0; Lk < pos2nuc[k].size();
										Lk++) {
									int Lk_nuc = pos2nuc[k][Lk];
									if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){	continue;}
									if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k - 1 and k
									//if(DEPflg && (k-1) - i + 1 == 2 && Dep2[ii2r[Rk1_nuc*10+L_nuc]][k-1] == 0){ continue;} // dependency between i and k - 1

									//cout << i << " " << k-1 << ":" << M[indx[k-1]+i][L][Rk1] << "," << k << " " << j << ":" << M[indx[j]+k][Lk][R] << endl;
									//int energy_M =  M[indx[k-1]+i][L][Rk1]+M[indx[j]+k][Lk][R];
									int energy_M =  M[getIndx(i,k-1,w_tmp,indx)][L][Rk1]+M[getIndx(k,j,w_tmp,indx)][Lk][R];
									DMl[i][L][R] = MIN2(energy_M, DMl[i][L][R]);
							        M[ij][L][R] = MIN2(energy_M, M[ij][L][R]);

								}
							}
						}


//						if(i == 3 && j == 7)
						//cout << i << " " << j << ":" << C[ij][L][R] << " " << L << "-" << R << endl;
						//vwxyがあるので、ここを複数回訪れることがある。
						//なので、MIN2を取っておく。
						//if(i2r[L_nuc] == NucConst[i] && i2r[R_nuc] == NucConst[j]){
							chkC[ij] = MIN2(chkC[ij], C[ij][L][R]);
							chkM[ij] = MIN2(chkM[ij], M[ij][L][R]);
						//} このループは多分意味がない。
					}
				}
			}
			// rotate DMl arrays
			int ***FF;
			FF = DMl2; DMl2 = DMl1; DMl1 = DMl; DMl =FF;
			for(int j = 1; j <= nuclen; j++){
				for(unsigned int L = 0; L < 4; L++){
					fill(DMl[j][L], DMl[j][L]+4,INF);
				}
			}
		}



		// Fill F matrix
		// Initialize F[1]
		for (unsigned int L = 0; L < pos2nuc[1].size(); L++) {
			for (unsigned int R = 0; R < pos2nuc[1].size(); R++) {
				F[1][L][R] = 0;
			}
		}

		for (unsigned int L1 = 0; L1 < pos2nuc[1].size(); L1++) {
			int L1_nuc = pos2nuc[1][L1];
			if(NCflg == 1 && i2r[L1_nuc] != NucConst[1]){	continue;}

			for (int j = 2; j <= nuclen; j++) {

//				int opt_flg_1 = 1;
//				for(int I = 0; I < n_inter; I++){
//					if(ofm[I] <= 1 && oto[I] >= 1){
//						opt_flg_1 = 0;
//						break;
//					}
//				}
				//				int opt_flg_j = 1;
				//				for(int I = 0; I < n_inter; I++){
				//					if(ofm[I] <= j && oto[I] >= j){
				//						opt_flg_j = 0;
				//						break;
				//					}
				//				}

				for (unsigned int Rj = 0; Rj < pos2nuc[j].size(); Rj++) {
					int Rj_nuc = pos2nuc[j][Rj];
					if(NCflg == 1 && i2r[Rj_nuc] != NucConst[j]){	continue;}

					if(DEPflg && j == 2 && Dep1[ii2r[L1_nuc*10+Rj_nuc]][1] == 0){ continue;}
					if(DEPflg && j == 3 && Dep2[ii2r[L1_nuc*10+Rj_nuc]][1] == 0){ continue;}


					F[j][L1][Rj] = INF;


					int type_L1Rj = BP_pair[i2r[L1_nuc]][i2r[Rj_nuc]];
					if (type_L1Rj) {
//						if(opt_flg_1 && opt_flg_j){
							int au_penalty = 0;
							if (type_L1Rj > 2)
								au_penalty = P->TerminalAU;
							if(j <= w_tmp)
								F[j][L1][Rj] = MIN2(F[j][L1][Rj], C[getIndx(1,j,w_tmp,indx)][L1][Rj] + au_penalty); // recc 1
								//F[j][L1][Rj] = MIN2(F[j][L1][Rj], C[indx[j] + 1][L1][Rj] + au_penalty); // recc 1
//						}
					}

					// create F[j] from F[j-1]
					for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size();
							Rj1++) {
						int Rj1_nuc = pos2nuc[j-1][Rj1];
						if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){	continue;}
						if(DEPflg && Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;}

						F[j][L1][Rj] = MIN2(F[j][L1][Rj], F[j - 1][L1][Rj1]); // recc 2
					}

					// create F[j] from F[k-1] and C[k][j]
					//for (int k = 2; k <= j - TURN - 1; k++) { // Is this correct?
					for (int k = MAX2(2, j-w_tmp+1); k <= j - TURN - 1; k++) { // Is this correct?

//						int opt_flg_k = 1;
//						for(int I = 0; I < n_inter; I++){
//							if(ofm[I] <= k && oto[I] >= k){
//								opt_flg_k = 0;
//								break;
//							}
//						}


						for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size();
								Rk1++) {
							int Rk1_nuc = pos2nuc[k-1][Rk1];
							if(NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]){	continue;}
							if(DEPflg && k == 3 && Dep1[ii2r[L1_nuc*10+Rk1_nuc]][1] == 0){ continue;} // dependency between 1(i) and 2(k-1)
							if(DEPflg && k == 4 && Dep2[ii2r[L1_nuc*10+Rk1_nuc]][1] == 0){ continue;} // dependency between 1(i) and 3(k-1)

							for (unsigned int Lk = 0; Lk < pos2nuc[k].size();
									Lk++) {
								int Lk_nuc = pos2nuc[k][Lk];
								if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){	continue;}

								if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k-1 and k

								int type_LkRj =
										BP_pair[i2r[Lk_nuc]][i2r[Rj_nuc]];

								int au_penalty = 0;
								if (type_LkRj > 2)
									au_penalty = P->TerminalAU;
								//int kj = indx[j] + k;
								int kj = getIndx(k,j,w_tmp,indx);

								int energy = F[k - 1][L1][Rk1] + C[kj][Lk][Rj]
										+ au_penalty; // recc 4

								F[j][L1][Rj] = MIN2(F[j][L1][Rj], energy);
							}

						}
					}

					//cout << j << ":" << F[j][L1][Rj] << " " << i2n[L1_nuc] << "-" << i2n[Rj_nuc] << endl;

					//test
					if (j == nuclen) {
						cout << i2n[L1_nuc] << "-" << i2n[Rj_nuc] << ":"
								<< F[j][L1][Rj] << endl;
					}
				}
			}
		}

		int minL, minR, MFE;
		MFE = INF;
		for(unsigned int L = 0; L < pos2nuc[1].size(); L++){
			int L_nuc = pos2nuc[1][L];
			if(NCflg == 1 && i2r[L_nuc] != NucConst[1]){	continue;}
			for(unsigned int R = 0; R < pos2nuc[nuclen].size(); R++){
				int R_nuc = pos2nuc[nuclen][R];
				if(NCflg == 1 && i2r[R_nuc] != NucConst[nuclen]){	continue;}

				if(F[nuclen][L][R] < MFE){
					MFE  = F[nuclen][L][R];
					minL = L;
					minR = R;
				}

			}
		}

		if(MFE == INF){
			printf("Mininum free energy is not defined.\n");
			exit(1);
		}


		if(rand_tb_flg){
			// Fill F2 matrix
			for (int l = 5; l <= nuclen; l++) {
				if(l > w_tmp) break;
			cout << "process F2:" << l << endl;

				for (int i = 1; i <= nuclen - l + 1; i++) {
					int j = i + l - 1;

					for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
						int L_nuc = pos2nuc[i][L];

						for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {
							int R_nuc = pos2nuc[j][R];
							int ij = getIndx(i,j,w_tmp,indx);

							F2[ij][L][R] = 0;

							int type = BP_pair[i2r[L_nuc]][i2r[R_nuc]];

							// from i, j-1 -> i, j
							for (unsigned int R1 = 0; R1 < pos2nuc[j-1].size(); R1++) {
								int R1_nuc = pos2nuc[j-1][R1];
								if(DEPflg && Dep1[ii2r[R1_nuc*10+R_nuc]][j-1] == 0){continue;}
								int ij1 = getIndx(i,j-1,w_tmp,indx);
								F2[ij][L][R] = MIN2(F2[ij][L][R], F2[ij1][L][R1]);
							}
							// from i-1, j -> i, j
							for (unsigned int L1 = 0; L1 < pos2nuc[i+1].size(); L1++) {
								int L1_nuc = pos2nuc[i+1][L1];
								if(DEPflg && Dep1[ii2r[L_nuc*10+L1_nuc]][i] == 0){continue;}
								int i1j = getIndx(i+1,j,w_tmp,indx);
								F2[ij][L][R] = MIN2(F2[ij][L][R], F2[i1j][L1][R]);
							}

							// from C
							int au_penalty = 0;
							if (type > 2)
								au_penalty = P->TerminalAU;
							if(j - i + 1 <= w_tmp){
								F2[ij][L][R] = MIN2(F2[ij][L][R], C[ij][L][R] + au_penalty);
								//cout << "test:" << F2[ij][L][R] << endl;
							}

							// Bifucation
							/* modular decomposition -------------------------------*/
							for (int k = i + 2 + TURN; k <= j - TURN - 1; k++) { // Is this correct?
								//cout << k << endl;
//			    				if((k - 1) - i + 1 > w_tmp ||
//			    					j - k + 1 > w_tmp)
//			    						continue;

								for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size();
										Rk1++) {
									int Rk1_nuc = pos2nuc[k-1][Rk1];

									for (unsigned int Lk = 0; Lk < pos2nuc[k].size();
											Lk++) {
										int Lk_nuc = pos2nuc[k][Lk];
										if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k - 1 and k

										int energy =  F2[getIndx(i,k-1,w_tmp,indx)][L][Rk1]+F2[getIndx(k,j,w_tmp,indx)][Lk][R];
										F2[ij][L][R] = MIN2(F2[ij][L][R], energy);

									}
								}
							}

//							cout << i << "," << j << "," << L << "," << R << "," << F2[ij][L][R] << endl;

						}
					}
				}
			}
		}

//		string optseq;
//		optseq.resize(nuclen+1, 'N');
//		optseq[0] = ' ';
//		backtrackR(&optseq, &*sector, &*base_pair, C, M, F,
//					indx, minL, minR, P, NucConst, pos2nuc, NCflg, i2r, nuclen, w_tmp, BP_pair, i2n, rtype, ii2r, Dep1, Dep2, DEPflg, predefHPN, predefHPN_E, substr, n2i, NucDef);


		if(rand_tb_flg){
			backtrack2(&optseq, &*sector, &*base_pair, C, M, F2,
					indx, minL, minR, P, NucConst, pos2nuc, NCflg, i2r, nuclen, w_tmp, BP_pair, i2n, rtype, ii2r, Dep1, Dep2, DEPflg, predefHPN, predefHPN_E, substr, n2i, NucDef);
		}
		else{
			backtrack(&optseq, &*sector, &*base_pair, C, M, F,
					indx, minL, minR, P, NucConst, pos2nuc, NCflg, i2r, nuclen, w_tmp, BP_pair, i2n, rtype, ii2r, Dep1, Dep2, DEPflg, predefHPN, predefHPN_E, substr, n2i, NucDef);
		}

		//塩基Nの修正
		for(int i = 1; i <= nuclen; i++){
			if(optseq[i] == 'N'){
				for(unsigned int R = 0; R < pos2nuc[i].size(); R++){ // check denendency with the previous and next nucleotide
					int R_nuc = pos2nuc[i][R];
					if(NCflg == 1 && i2r[R_nuc] != NucConst[i]){	continue;}

					if(i != 1 && optseq[i-1] != 'N'){ // check consistensy with the previous nucleotide
						int R_prev_nuc = n2i[optseq[i-1]];
						if(DEPflg && Dep1[ii2r[R_prev_nuc*10+R_nuc]][i-1] == 0){ continue;}
					}
					if(i != nuclen && optseq[i+1] != 'N'){ // check consistensy with the next nucleotide
						int R_next_nuc = n2i[optseq[i+1]];
						if(DEPflg && Dep1[ii2r[R_nuc*10+R_next_nuc]][i] == 0){ continue;}
					}
					if(i < nuclen - 1 && optseq[i+2] != 'N'){ // check consistensy with the next nucleotide
						int R_next_nuc = n2i[optseq[i+2]];
						if(DEPflg && Dep2[ii2r[R_nuc*10+R_next_nuc]][i] == 0){ continue;}
					}

					optseq[i] = i2n[R_nuc];
					break;
				}
				//cout << i << ":" << optseq[i] << endl;
			}
		}
		//塩基V,W,X,Yの修正
		for(int i = 1; i <= nuclen; i++){

			if(optseq[i] == 'V' || optseq[i] == 'W'){
				optseq[i] = 'U';
			}
			else if(optseq[i] == 'X' || optseq[i] == 'Y'){
				optseq[i] = 'G';
			}
		}


		//2次構造情報の表示
		string optstr;
		optstr.resize(nuclen+1, '.');
		optstr[0] = ' ';
		for(int i = 1; i <= base_pair[0].i; i++){
			optstr[base_pair[i].i] = '(';
			optstr[base_pair[i].j] = ')';
		}

		//show original amino acids
		for(int i = 0; i < aalen; i++){
			cout << aaseq[i] << "  ";
		}
		cout << endl;
		//check amino acids of desinged DNA
		int j = 0;
		for(unsigned int i = 1; i < optseq.size(); i = i+3){
			char aa = codon_table.c2a(n2i[optseq[i]], n2i[optseq[i+1]], n2i[optseq[i+2]]);
			cout << aa << "  ";
			if(aaseq[j] != aa){
				cerr << j+1 << "-th amino acid differs:" << aaseq[j] << ":" << aa << endl;
			}
			j++;
		}
		cout << endl;
		//最適塩基と2次構造の表示
		string optseq_disp = optseq;
		optseq_disp.erase(0, 1);
		optstr.erase(0, 1);
		cout << optseq_disp << endl;
		cout << optstr << endl;
		cout << "MFE:" << float(MFE)/100 << " kcal/mol" << endl;


		if(part_opt_flg == 1){
			//部分アミノ酸配列の作成
			for(int I = 0; I< n_inter; I++){
				int aa_fm = (ofm[I]-1)/3 + 1; // 1-based
				int aa_to = oto[I]/3; // 1-based
				int part_aalen = aa_to - aa_fm + 1;

				char part_aaseq[part_aalen];
				part_aaseq[part_aalen] = '\0';
				int j = 0;
				for(int i = aa_fm; i <= aa_to; i++){
					part_aaseq[j++] = aaseq[i-1]; // convert to 0-based
				}

				cout << aa_fm << ":" << aa_to << endl;
				cout << part_aalen << endl;
				cout << part_aaseq << endl;


				string part_optseq = rev_fold_step1(part_aaseq, part_aalen, codon_table, exc);
				rev_fold_step2(&part_optseq, part_aaseq, part_aalen, codon_table, exc);
				// combine optseq_rev and optseq
				j = 1;
				for(int i = ofm[I]; i <= oto[I]; i++){
					optseq[i] = part_optseq[j++];
				}
			}
			fixed_fold(optseq, indx, w_tmp, predefHPN_E, BP_pair, P, aaseq, codon_table);
			//fixed_fold(optseq, indx, w_tmp, predefHPN_E, BP_pair, P, aaseq, codon_table);
		}

		if(m_disp){
			// get process ID

			pid_t pid = getpid();
			stringstream ss;
			ss << "/proc/" << pid << "/status";
			int m = getMemoryUsage(ss.str());
			if(m == -1){
				cerr << "Cannot get memory usage information from " << ss.str() << endl;
			}
			else{
				cout << "Memory(VmRSS): "  << float(m)/1024 << " Mb" << endl;
			}
		}

		free_arrays(nuclen, indx, w_tmp, pos2nuc, &C, &M, &F, &DMl, &DMl1, &DMl2, &chkC, &chkM, &base_pair);
		if(rand_tb_flg)
			free_F2(nuclen, indx, w_tmp, pos2nuc, &F2);

		free(P);

	} while (all_aaseq.next());

	clock_t end = clock();
	float sec = (double)end/CLOCKS_PER_SEC;
	float min = sec/60;
	cout << "Runing time: " <<  min << " minutes" << endl;

//	struct rusage r;
//	if (getrusage(RUSAGE_SELF, &r) != 0) {
//			/*Failure*/
//	}
//	printf("Memory usage: %ld Mb\n", r.ru_maxrss/1024);


	//	printf("Memory usage: %ld Mb\n", r.ru_maxrss/1024);

	return 0;

}

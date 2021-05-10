#include <map>
#include <algorithm> // only to use fill
#include <fstream>
//#include <iostream>
//#include <stdlib.h>
//#include <codon.hpp>
using namespace std;

inline int E_hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P);
inline int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P);

typedef struct stack {
	int i;
	int j;
	int Li;
	int Rj;
	int ml;
} stack;

typedef struct bond {
	int i;
	int j;
} bond;


inline int TermAU(int const &type, paramT * const &P);

int getMatrixSize(int len, int w){
	int size = 0;
	for(int i = 1; i <= w; i++){
		size += len-(i-1); // マトリクスの斜めの要素数を足し合わせるイメージ
							// i=1のときは、対角要素（len個）を足し合わせる。
	}

	cout << "The size of matrix is " << size << endl;
	return size;
}

inline int getIndx(int const &i, int const &j, int const &w, int *const &indx){
	return indx[j] + i - MAX2(0, j - w); // j-wは使わない要素の数。
										  // wが指定されない(=length)と、j列にはj個分(1<i<j)の要素が用意される。
										  // wが指定されると、j列にはで使う要素はw個、使わない要素はj-w個となる。
}

void clear_sec_bp(stack *s, bond *b, int len){
	for(int i = 0; i < 500; i++){
		s[i].i = -INF;
		s[i].j = -INF;
		s[i].Li = -INF;
		s[i].Rj = -INF;
		s[i].ml = -INF;
	}
	for(int i = 0; i < len/2; i++){
		b[i].i = -INF;
		b[i].j = -INF;
	}

}


void allocate_arrays(int len, int *indx, int w, vector <vector<int> > &pos2nuc, int ****c, int ****m, int ****f, int ****dml, int ****dml1, int ****dml2, int **chkc, int **chkm, bond **b)
{
	int size = getMatrixSize(len, w);
//	int n_elm = 0;
	int total_bytes = 0;
	//int test_bytes = 0;

	//	*c   = new int**[len*(len+1)/2+1];
	//	*m   = new int**[len*(len+1)/2+1];
	*c   = new int**[size+1];
	*m   = new int**[size+1];
//	*f2   = new int**[size+1];
	for(int i = 1; i <= len; i++){
		for(int j = i; j <= MIN2(len, i+w-1); j++){
			//cout << i << " " << j << endl;
			//int ij = indx[j]+i;
			int ij = getIndx(i,j,w,indx);
			(*c)[ij]   = new int*[pos2nuc[i].size()];
			(*m)[ij]   = new int*[pos2nuc[i].size()];
//			(*f2)[ij]   = new int*[pos2nuc[i].size()];
			//total_bytes += sizeof(int*) * (pos2nuc[i].size()+1) * 2 * 2; // *2 is empirical constant
			total_bytes += sizeof(int*) * (pos2nuc[i].size()+4) * 2; // +4 is an empirical value
			for(unsigned int L = 0; L < pos2nuc[i].size(); L++){
				(*c)[ij][L]   = new int[pos2nuc[j].size()];
				(*m)[ij][L]   = new int[pos2nuc[j].size()];
//				(*f2)[ij][L]   = new int[pos2nuc[j].size()];
				//total_bytes += sizeof(int) * (pos2nuc[j].size()+1) * 2 * 2; // *2 is empirical constant
//				n_elm += (pos2nuc[j].size()+1);
				total_bytes += sizeof(int) * (pos2nuc[j].size()+4)*2; // +4 is an empirical value
			}
		}
	}

	total_bytes *= 1.2; // 1.2 is an empirical value
//	cout << "sizeof(int): " << sizeof(int) << endl;
//	cout << "sizeof(int*): " << sizeof(int***) << endl;
//	cout << "N_ELM: "<< n_elm << endl;
//	cout << "test_bytes: "<< (float)test_bytes/(1024*1024) << endl;
//	return test_bytes;

//	int total_bytes = n_elm * sizeof(int) * 2; // *2 is m and c
//	float total_Mb = (float)total_bytes/(1024*1024);
//	cout << "estimated memory size: " << total_Mb << " Mb" << endl;

	*f   = new int**[len+1];
	*dml  = new int**[len+1];
	*dml1  = new int**[len+1];
	*dml2  = new int**[len+1];
	for(int j = 1; j <= len; j++){
		(*f)[j]   = new int*[pos2nuc[j].size()];
		for(unsigned int L = 0; L < pos2nuc[1].size(); L++){ // The first position
			(*f)[j][L]   = new int[pos2nuc[j].size()];
		}

		(*dml)[j]   = new int*[4]; // always secure 4 elements, because the maximum number of nucleotides is 4
		(*dml1)[j]   = new int*[4]; // always secure 4 elements, because the maximum number of nucleotides is 4
		(*dml2)[j]   = new int*[4]; // always secure 4 elements, because the maximum number of nucleotides is 4
		for(unsigned int L = 0; L < 4; L++){
			(*dml)[j][L]   = new int[4]; // always secure 4 elements, because the maximum number of nucleotides is 4
			(*dml1)[j][L]   = new int[4]; // always secure 4 elements, because the maximum number of nucleotides is 4
			(*dml2)[j][L]   = new int[4]; // always secure 4 elements, because the maximum number of nucleotides is 4
			fill((*dml)[j][L], (*dml)[j][L]+4,INF);
			fill((*dml1)[j][L], (*dml1)[j][L]+4,INF);
			fill((*dml2)[j][L], (*dml2)[j][L]+4,INF);
		}
	}

	//	*chkc   = new int[len*(len+1)/2+1];
	//	*chkm   = new int[len*(len+1)/2+1];

	//	fill(*chkc, *chkc+len*(len+1)/2, INF);
	//	fill(*chkm, *chkm+len*(len+1)/2, INF);

	*chkc   = new int[size+1];
	*chkm   = new int[size+1];

	fill(*chkc, *chkc+size, INF);
	fill(*chkm, *chkm+size, INF);

	*b      = new bond[len/2];

	//	return (float)total_bytes/(1024*1024) * 1.25; // 1.2 is empirical constant
	//return (float)total_bytes/(1024*1024);

}

void allocate_F2(int len, int *indx, int w, vector <vector<int> > &pos2nuc, int ****f2)
{
	int size = getMatrixSize(len, w);
	*f2   = new int**[size+1];
	for(int i = 1; i <= len; i++){
		for(int j = i; j <= MIN2(len, i+w-1); j++){
			int ij = getIndx(i,j,w,indx);
			(*f2)[ij]   = new int*[pos2nuc[i].size()];
			for(unsigned int L = 0; L < pos2nuc[i].size(); L++){
				(*f2)[ij][L]   = new int[pos2nuc[j].size()];
			}
		}
	}
}


void free_arrays(int len, int *indx, int w, vector <vector<int> > &pos2nuc, int ****c, int ****m, int ****f, int ****dml, int ****dml1, int ****dml2, int **chkc, int **chkm, bond **b)
{
	for(int i = 1; i <= len; i++){
		for(int j = i; j <= MIN2(len, i+w-1); j++){
			//int ij = indx[j]+i;
			int ij = getIndx(i,j,w,indx);
			for(unsigned int L = 0; L < pos2nuc[i].size(); L++){
				delete [] (*c)[ij][L];
				delete [] (*m)[ij][L];
			}
			delete [] (*c)[ij];
			delete [] (*m)[ij];
		}
	}
	delete [] *c;
	delete [] *m;

	for(int j = 1; j <= len; j++){
		for(unsigned int L = 0; L < pos2nuc[1].size(); L++){ // The first position
			delete [] (*f)[j][L];
		}
		delete [] (*f)[j];

		for(unsigned int L = 0; L < 4; L++){
			delete [] (*dml)[j][L];
			delete [] (*dml1)[j][L];
			delete [] (*dml2)[j][L];
		}
		delete [] (*dml)[j];
		delete [] (*dml1)[j];
		delete [] (*dml2)[j];
	}
	delete [] *f;
	delete [] *dml;
	delete [] *dml1;
	delete [] *dml2;

	delete [] *chkc;
	delete [] *chkm;

	delete [] *b;

}

void free_F2(int len, int *indx, int w, vector <vector<int> > &pos2nuc, int ****f2)
{
	for(int i = 1; i <= len; i++){
		for(int j = i; j <= MIN2(len, i+w-1); j++){
			//int ij = indx[j]+i;
			int ij = getIndx(i,j,w,indx);
			for(unsigned int L = 0; L < pos2nuc[i].size(); L++){
				delete [] (*f2)[ij][L];
			}
			delete [] (*f2)[ij];
		}
	}
	delete [] *f2;

}


void set_ij_indx(int *a, int length)
{
	for (int n = 1; n <= length; n++){
		a[n] = (n*(n-1))/2;
	}
}


void set_ij_indx(int *a, int length, int w)
{
	if(w <= 0){
		cerr << "Invalid w:" << w << endl;
		exit(1);
	}
	w = MIN2(length, w);
	int cum = 0;
	for (int n = 1; n <= length; n++){
		a[n] = cum;
		//cout << n << ":" << a[n] << endl;
		if(n < w){
			cum += n;
		}
		else{
			cum += w;
		}
	}
}

void set_arrays(int **a, int length)
{
	*a = new int[length+1]; //

	for (int n = 1; n <= length; n++){
		(*a)[n] = (n*(n-1))/2;
		//cout << *a[n] << endl;
	}
}

void make_i2r(int *n){
	n[1] = 1;
	n[2] = 2;
	n[3] = 3;
	n[4] = 4;
	n[5] = 4; // 2nd position of L is converted to U
	n[6] = 4;
	n[7] = 3; // 2nd position of R is converted to G
	n[8] = 3;

}

void make_ii2r(int *n){
	int s = 1;
	for(int i1 = 1; i1 <= 8; i1++){
		for(int i2 = 1; i2 <= 8; i2++){
			n[i1*10+i2] = s++;
		}
	}

}

map<char, int> make_n2i(){
	map<char, int> m;
	m['A'] = 1;
	m['C'] = 2;
	m['G'] = 3;
	m['U'] = 4;
	m['V'] = 5; // 2nd position of L, before A/G
	m['W'] = 6; // 2nd position of L, before U/C
	m['X'] = 7; // 2nd position of R, before A/G
	m['Y'] = 8; // 2nd position of R, before U/C

	return m;
}

//void view_n2i(map<char, int> n2i, char const &c){
int view_n2i(map<char, int> n2i, char c){
	//cout << c << endl;
	return n2i[c];
}

void make_i2n(char *n){
	n[0] = ' ';
	n[1] = 'A';
	n[2] = 'C';
	n[3] = 'G';
	n[4] = 'U';
	n[5] = 'V';
	n[6] = 'W';
	n[7] = 'X';
	n[8] = 'Y';
}

void showPos2Nuc(vector<vector<int> > &v){

	  for(unsigned int i = 1; i <= v.size(); i++){
		  cout << i;
		  for(unsigned int j = 0; j < v[i].size(); j++){
		  	  cout << "\t" << v[i][j];
		  }
		  cout << endl;
	  }

}

void showPos2Nuc(vector<vector<int> > &v, char i2n[]){

		  for(unsigned int i = 1; i < v.size(); i++){
	//	  for(unsigned int i = 1; i <= v.size(); i++){ // This is segmentaion fault
		  cout << i;
		  for(unsigned int j = 0; j < v[i].size(); j++){
		  	  cout << "\t" << i2n[v[i][j]];
		  }
		  cout << endl;
	  }

}

vector<vector<int> > getPossibleNucleotide(char *aaseq, int aalen, codon &codon_table, map<char,int> &n2i, string excludedCodons){

 vector<vector<int> > v;

 int nuclen = aalen*3;
 v.resize(nuclen+1);

 for(int i = 0; i < aalen; i++){
	  int nucpos = i*3+1;
	  char aa = aaseq[i];
	  //cout << "test: " << aa << endl;
	  vector<string> codons = codon_table.getCodons(aa, excludedCodons);
	  for(int k = 0; k < 3; k++){ // for each codon position
		  nucpos += k;

		  if(aa == 'L' && k == 1){
			  v[nucpos].push_back(n2i['V']);
			  v[nucpos].push_back(n2i['W']);
		  }
		  else if(aa == 'R' && k == 1){
			  v[nucpos].push_back(n2i['X']);
			  v[nucpos].push_back(n2i['Y']);
		  }
		  else{
			  bool flg_A,flg_C,flg_G,flg_U;
			  flg_A = flg_C = flg_G = flg_U = 0;

			  for(unsigned int j = 0; j < codons.size(); j++){ // for each codon corresp. to each aa
				  char nuc = codons[j][k];
				  if(nuc == 'A' && flg_A == 0){
					  v[nucpos].push_back(n2i[nuc]);
					  flg_A = 1;
				  }
				  else if(nuc == 'C' && flg_C == 0){
					  v[nucpos].push_back(n2i[nuc]);
					  flg_C = 1;
				  }
				  else if(nuc == 'G' && flg_G == 0){
					  v[nucpos].push_back(n2i[nuc]);
					  flg_G = 1;
				  }
				  else if(nuc == 'U' && flg_U == 0){
					  v[nucpos].push_back(n2i[nuc]);
					  flg_U = 1;
				  }
			  }

		  }
		  nucpos -=k;

	  }
 }
  return v;
}

void showChkMatrix(int *&m, int *&indx, int len, int w){
	printf("REF:");
	for(int j = 1; j <= len; j++){
		printf("\t%d", j);
	}
	printf("\n");

	for(int i = 1; i <= len; i++){
		printf("REF:");
		printf("%d", i);
		for(int j = 1; j <= len; j++){
			if(i>=j){
				printf("\t-");
			}
			else{
				//printf("\t%d", m[indx[j]+i]);
				printf("\t%d", m[getIndx(i,j,w,indx)]);
			}
		}
		printf("\n");
	}

}

void showFixedMatrix(const int *m, int *indx, const int len, const int w){
	printf("REF:");
	for(int j = 1; j <= len; j++){
		printf("\t%d", j);
	}
	printf("\n");

	for(int i = 1; i <= len; i++){
		printf("REF:");
		printf("%d", i);
		for(int j = 1; j <= len; j++){
			if(i>=j){
				printf("\t-");
			}
			else{
				//printf("\t%d", m[indx[j]+i]);
				printf("\t%d", m[getIndx(i,j,w,indx)]);
			}
		}
		printf("\n");
	}

}


vector<int> createNucConstraint(const char* s, int &len, map<char, int> &n2i){
	vector<int> v;
	v.resize(len+1);
	for(int i = 1; i <= len; i++){
		v[i] = n2i[s[i]];
	}
	return v;
}

inline void InitRand()
{
    srand((unsigned int)time(NULL));
}

void shuffleStr(vector<string> (*ary),int size)
{
    for(int i=0;i<size;i++)
    {
        int j = rand()%size;
        string t = (*ary)[i];
        (*ary)[i] = (*ary)[j];
        (*ary)[j] = t;
    }
}
void shuffle(int ary[],int size)
{
    for(int i=0;i<size;i++)
    {
        int j = rand()%size;
        int t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }
}
vector <pair<int, int> > shufflePair(vector<pair<int, int> > ary, int size)
{
    for(int i=0;i<size;i++)
    {
        int j = rand()%size;
        //cout << rand() << ":" << j << endl;
        pair<int, int> t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }

    return ary;
}



void backtrack(string *optseq, stack *sector, bond *base_pair, int ***const &c, int ***const &m, int*** const &f,
			int *const indx, const int &initL, const int &initR, paramT *const&P, const vector<int> &NucConst,
			const vector<vector <int> > &pos2nuc, const int &NCflg, int *const &i2r, int const &length, int const &w,
			int const (&BP_pair)[5][5], char * const &i2n, int * const &rtype, int *const &ii2r,
			vector<vector<int> > &Dep1, vector<vector<int> > &Dep2, int &DEPflg,
			vector<vector<vector<vector<pair<int, string> > > > > &predefH, map<string, int> &predefE, vector<vector<vector<string> > > &substr, map<char, int> &n2i, const char* nucdef){

	int s = 0;
	int b = 0;
	sector[++s].i = 1;
	sector[s].j = length;
	sector[s].Li = initL;
	sector[s].Rj = initR;
	sector[s].ml = 0;


	OUTLOOP:
	while (s>0) {
	    int fij, fi, ij, cij, traced, traced_Lk, i1, j1, k, p , q;
//	    int canonical = 1;     /* (i,j) closes a canonical structure */

	    //ここでoptseqに値を反映させるべき？
	    int i  = sector[s].i;
	    int j  = sector[s].j;
	    int Li  = sector[s].Li;
	    int Rj  = sector[s].Rj;
	    int ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
	                              	  occur in the M- (1) or in the F-array (0) */
	    int Li_nuc = pos2nuc[i][Li];
	    int Rj_nuc = pos2nuc[j][Rj];

	    if(i + 1 == j && Dep1[ii2r[Li_nuc*10+Rj_nuc]][i] == 0){ continue;}
	    if(i + 2 == j && Dep2[ii2r[Li_nuc*10+Rj_nuc]][i] == 0){ continue;}

	    int type_LiRj = BP_pair[i2r[Li_nuc]][i2r[Rj_nuc]];

	    (*optseq)[i] = i2n[Li_nuc];
	    (*optseq)[j] = i2n[Rj_nuc];

	    if (ml==2) {
	      base_pair[++b].i = i;
	      base_pair[b].j   = j;
	      goto repeat1;
	    }

	    //	    if (j < i+TURN+1) continue; /* no more pairs in this interval */
	    if (j == i) break;

	    //	    fij = (ml == 1)? m[indx[j]+i][Li][Rj] : f[j][Li][Rj];
	    fij = (ml == 1)? m[getIndx(i,j,w,indx)][Li][Rj] : f[j][Li][Rj];
//	    if(TB_CHK_flg == 1)
	    cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << fij << ")" << ":" << *optseq << endl;

	    for(unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size(); Rj1++){
	    	int Rj1_nuc = pos2nuc[j-1][Rj1];
	    	if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){continue;}

		    if(Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;}

		    //	    	fi  = (ml == 1)? m[indx[j-1]+i][Li][Rj1] + P->MLbase: f[j-1][Li][Rj1];
		    fi  = (ml == 1)? m[getIndx(i,j-1,w,indx)][Li][Rj1] + P->MLbase: f[j-1][Li][Rj1];

	        if (fij == fi) {  /* 3' end is unpaired */
	          sector[++s].i = i;
	          sector[s].j   = j-1;
	          sector[s].Li   = Li;
	          sector[s].Rj   = Rj1;
	          sector[s].ml  = ml;
	          //continue;
	          goto OUTLOOP;
	        }

	    }

	    if (ml == 0) { /* backtrack in f */

	    	if(i != 1){
    			cerr << "Traceback failure: i must be 1 during bachtrack in f" << endl;
	    	}

	    	//f[j]とC[1][j]が一致している時の処理。Vieenaでは、次for文に統合されている。
	    	if(type_LiRj && j <= w){
	    		// note i == 1
	    		//int en_c = TermAU(type_LiRj, P) +  c[indx[j]+i][Li][Rj];
	    		int en_c = TermAU(type_LiRj, P) +  c[getIndx(i,j,w,indx)][Li][Rj];
                int en_f = f[j][Li][Rj];
              	if(en_c ==  en_f){
              		k = i;
              		traced = j;
               		traced_Lk = Li;
               		goto LABEL1;
              	}
	    	}

	    	//for(k=j-TURN-1,traced=0; k>=1; k--){
	    	for(k=j-TURN-1,traced=0; k>=MAX2(2,j-w+1); k--){

	    		for(unsigned int Rk1 = 0; Rk1 < pos2nuc[k-1].size(); Rk1++){
	    	    	int Rk1_nuc = pos2nuc[k-1][Rk1];
	    	    	if(NCflg == 1 && i2r[Rk1_nuc] != NucConst[k-1]){continue;}
	    		    if(DEPflg && k == 3 && Dep1[ii2r[Li_nuc*10+Rk1_nuc]][i] == 0){ continue;} // dependency between 1(i) and 2(k-1)
	    		    if(DEPflg && k == 4 && Dep2[ii2r[Li_nuc*10+Rk1_nuc]][i] == 0){ continue;} // dependency between 1(i) and 3(k-1)

	    	    	for(unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++){
		    			int Lk_nuc = pos2nuc[k][Lk];
		    	    	if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){continue;}
		    		    if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k-1 and k

		    	    	int type_LkRj = BP_pair[i2r[Lk_nuc]][i2r[Rj_nuc]];
	                    if(type_LkRj){
	                    	//	int en_c = TermAU(type_LkRj, P) +  c[indx[j]+k][Lk][Rj];
	                    	int en_c = TermAU(type_LkRj, P) +  c[getIndx(k,j,w,indx)][Lk][Rj];
	                    	int en_f = f[k-1][Li][Rk1];
//	                    	int test = en_c + en_f;
//	                    	cout << fij << "=" << test << "(" << en_c << "+" << en_f << ")" << i << ":" << k << ":" << j << endl;
	                    	if(fij ==  en_c + en_f){
	                    		traced = j;
	                    		traced_Lk = Lk;
	                    		/* push back the remaining f portion */
	                    		sector[++s].i = i;
	                    		sector[s].j   = k-1;
	                    		sector[s].Li = Li;
	                    		sector[s].Rj = Rk1;
	                    		sector[s].ml  = 0;

	                    		goto LABEL1;
	                    	}
	                    }
		    		}
	    		}
	    	}
	    	LABEL1:

	    	if (!traced){
	    		fprintf(stderr, "backtrack failed in f\n");
	    		fprintf(stderr, "cannot trace f[%d][%d][%d] Lnuc=%c Rnuc=%c \n", j, Li, Rj, i2n[Li_nuc], i2n[Rj_nuc]);
	    		exit(0);
	    	}

	    	/* trace back the base pair found */
	    	// [1]
	    	i=k;            // iを更新
	    	j=traced;       // この代入は多分必要ない。jに関しては不変だから
	    	Li = traced_Lk; // Liを更新
	    	//Rjは不変
	    	base_pair[++b].i = i;
	    	base_pair[b].j   = j;
	    	goto repeat1;
	    }
	    else { /* trace back in fML array */

		    for(unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size(); Li1++){
		    	int Li1_nuc = pos2nuc[i+1][Li1];
		    	if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){continue;}
    		    if(DEPflg && Dep1[ii2r[Li_nuc*10+Li1_nuc]][i] == 0){ continue;} // dependency between k-1 and k

    		    //	  	if (m[indx[j]+i+1][Li1][Rj]+P->MLbase == fij) { /* 5' end is unpaired */
    		    if (m[getIndx(i+1,j,w,indx)][Li1][Rj]+P->MLbase == fij) { /* 5' end is unpaired */
	    			sector[++s].i = i+1;
	    			sector[s].j   = j;
	    			sector[s].Li  = Li1;
	    			sector[s].Rj  = Rj;
	    			sector[s].ml  = ml;
	    			goto OUTLOOP;
		    	}

    		    //	ij  = indx[j]+i;
    		    ij  = getIndx(i,j,w,indx);

		    	if(fij == c[ij][Li][Rj] + TermAU(type_LiRj, P) + P->MLintern[type_LiRj]){
		    		base_pair[++b].i = i;
		    		base_pair[b].j   = j;
		    		goto repeat1;
		    	}

		    }

		    //		    for(k = i + 1 + TURN; k <= j - 2 - TURN; k++){
		    for(k = i + 2 + TURN; k <= j - 1 - TURN; k++){

	    		for(unsigned int Rk1 = 0; Rk1 < pos2nuc[k-1].size(); Rk1++){
	    	    	int Rk1_nuc = pos2nuc[k-1][Rk1];
	    	    	if(NCflg == 1 && i2r[Rk1_nuc] != NucConst[k-1]){continue;}
	    	    	// check dependency is not needed because i+2<k,k+2<J
	    	    	for(unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++){
		    			int Lk_nuc = pos2nuc[k][Lk];
		    	    	if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){continue;}
		    		    if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k-1 and k

		    		    //if(fij == (m[indx[k-1]+i][Li][Rk1]+m[indx[j]+k][Lk][Rj])){
		    		    if(fij == (m[getIndx(i,k-1,w,indx)][Li][Rk1]+m[getIndx(k,j,w,indx)][Lk][Rj])){
		    	    		sector[++s].i = i;
		    	    		sector[s].j   = k-1;
		    	    		sector[s].Li  = Li;
		    	    		sector[s].Rj  = Rk1;
		    	    		sector[s].ml  = ml;
		    	    		sector[++s].i = k;
		    	    		sector[s].j   = j;
		    	    		sector[s].Li  = Lk;
		    	    		sector[s].Rj  = Rj;
		    	    		sector[s].ml  = ml;
		    	    		goto OUTLOOP;
		    	    	}
		    		}
	    		}

		    }
    		if (k>j-1-TURN){
    			fprintf(stderr, "backtrack failed in fML\n");
    			exit(1);
    		}
	    }

	    repeat1: // いちいちスタックに積まずに、ここで部分的なトレースバックをしてしまう。
	    //continue;
	    /*----- begin of "repeat:" -----*/
	    if(j - i + 1 > w){
	    	cerr << "backtrack failed at " << i << "," << j << " : the length must at most << w << endl";
	    }
	    //	    ij = indx[j]+i; // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
	    ij = getIndx(i,j,w,indx); // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
	    Li_nuc = pos2nuc[i][Li]; // Liは更新されている。
	    Rj_nuc = pos2nuc[j][Rj]; //Rj_nucは更新されていない場合もある[1]。
	    type_LiRj = BP_pair[i2r[Li_nuc]][i2r[Rj_nuc]];
	    cij = c[ij][Li][Rj];
		(*optseq)[i] = i2n[Li_nuc]; //塩基対部分を記録
		(*optseq)[j] = i2n[Rj_nuc];

//		if(TB_CHK_flg == 1)
//			cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << cij << ")" << ":" << *optseq << endl;


	    // predefinedなヘアピンのトレースバック
//		if ((j-i+1 == 5 ||j-i+1 == 6 ||j-i+1 == 8) &&
 //   			cij == predefH[i][j-i+1][i2r[Li_nuc]][i2r[Rj_nuc]].first){
//			string s1 = predefH[i][j-i+1][i2r[Li_nuc]][i2r[Rj_nuc]].second;
//			for(unsigned int k = 0; k < s1.size(); k++){
//				(*optseq)[i+k] = s1[k]; //塩基を記録
//			}
//			cout << "Predefined Hairpin " << s1 << " at " << i << ":" << j << endl;
//			goto OUTLOOP;
 //   	}
		if (j-i+1 == 5 ||j-i+1 == 6 ||j-i+1 == 8){
			int l = j-i+1;
			for(unsigned int s = 0; s < substr[i][l].size(); s++){
				string hpn = substr[i][l][s];
				int hL_nuc  = n2i[hpn[0]];
				int hL2_nuc = n2i[hpn[1]];
				int hR2_nuc = n2i[hpn[l-2]];
				int hR_nuc  = n2i[hpn[l-1]];
				if(hL_nuc != i2r[Li_nuc]) continue;
				if(hR_nuc != i2r[Rj_nuc]) continue;

				if(NCflg == 1){
					string s1 = string(nucdef).substr(i, l);
					if(hpn != s1) continue;
				}


				if(DEPflg && Li_nuc > 4 && Dep1[ii2r[Li_nuc*10+hL2_nuc]][i] == 0){continue;} // Dependency is already checked.
				if(DEPflg && Rj_nuc > 4 && Dep1[ii2r[hR2_nuc*10+Rj_nuc]][j-1] == 0){continue;}// ただし、Li_nuc、Rj_nucがVWXYのときだけは、一つ内側との依存関係をチェックする必要がある。

				// predefinedなヘアピンとの比較
				if(predefE.count(hpn) > 0){
					if(c[ij][Li][Rj] == predefE[hpn]){
						cout << "Predefined Hairpin at " << i << "," << j << endl;
						for(unsigned int k = 0; k < hpn.size(); k++){
							(*optseq)[i+k] = hpn[k]; //塩基を記録
						}
						goto OUTLOOP;
					}

				}
				else{	// 普通のヘアピンとの比較
					//					int energy = HairpinE(j - i - 1, type_LiRj,
					//							i2r[hL2_nuc], i2r[hR2_nuc],
					//							"NNNNNNNNN");
					int energy = E_hairpin(j - i - 1, type_LiRj,
												i2r[hL2_nuc], i2r[hR2_nuc],
												"NNNNNNNNN", P);

					if(c[ij][Li][Rj] == energy){
						for(unsigned int k = 0; k < hpn.size(); k++){
							(*optseq)[i+k] = hpn[k]; //塩基を記録
						}
						goto OUTLOOP;
					}
				}
			}

		}
		else{
			// 普通のヘアピンのトレースバック
			for(unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size(); Li1++){
				int Li1_nuc = pos2nuc[i+1][Li1];
				if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){continue;}
				if(DEPflg && Dep1[ii2r[Li_nuc*10+Li1_nuc]][i] == 0){ continue;} // dependency between i and i+1

				for(unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size(); Rj1++){
					int Rj1_nuc = pos2nuc[j-1][Rj1];
					if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){continue;}
					if(DEPflg && Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;} // dependency between j-1 and j

					//if (cij == HairpinE(j-i-1, type_LiRj, i2r[Li1_nuc], i2r[Rj1_nuc], "NNNNNNNNN")){
					if (cij == E_hairpin(j-i-1, type_LiRj, i2r[Li1_nuc], i2r[Rj1_nuc], "NNNNNNNNN", P)){
						(*optseq)[i+1] = i2n[Li1_nuc]; //塩基対の内側のミスマッチ塩基を記録
						(*optseq)[j-1] = i2n[Rj1_nuc];
						goto OUTLOOP;
					}
					else{
						continue;
					}
				}
			}
		}
	    // Hairpinに該当がなければ、Internal loopのトレースバック。もっとも手強い。
	    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
		    for (unsigned int Lp = 0; Lp < pos2nuc[p].size() ; Lp++) {
		    	int Lp_nuc = pos2nuc[p][Lp];
		    	if(NCflg == 1 && i2r[Lp_nuc] != NucConst[p]){continue;}
		    	if(DEPflg && p == i + 1 && Dep1[ii2r[Li_nuc*10+Lp_nuc]][i] == 0){ continue;} // dependency between i and q
		    	if(DEPflg && p == i + 2 && Dep2[ii2r[Li_nuc*10+Lp_nuc]][i] == 0){ continue;} // dependency between i and q

		    	int minq = j-i+p-MAXLOOP-2;
		    	if (minq<p+1+TURN) minq = p+1+TURN;
		    	for (q = j-1; q >= minq; q--) {
				    for (unsigned int Rq = 0; Rq < pos2nuc[q].size() ; Rq++) {
				    	int Rq_nuc = pos2nuc[q][Rq];
				    	if(NCflg == 1 && i2r[Rq_nuc] != NucConst[q]){continue;}
				    	if(DEPflg && q == j - 1 && Dep1[ii2r[Rq_nuc*10+Rj_nuc]][q] == 0){ continue;} // dependency between q and j
				    	if(DEPflg && q == j - 2 && Dep2[ii2r[Rq_nuc*10+Rj_nuc]][q] == 0){ continue;} // dependency between q and j

				    	int type_LpRq = BP_pair[i2r[Lp_nuc]][i2r[Rq_nuc]];
				    	if (type_LpRq==0) continue;
				    	type_LpRq = rtype[type_LpRq];

					    for (unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size() ; Li1++) {
					    	int Li1_nuc = pos2nuc[i+1][Li1];
					    	if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){continue;}
					    	if(DEPflg && Dep1[ii2r[Li_nuc*10+Li1_nuc]][i] == 0){ continue;} // dependency between i and i+1
					    	if(i+1 == p && Li1_nuc != Lp_nuc){ continue; } // i,pの時は、i+1の塩基とpの塩基は一致していないといけない。(1)

				    		for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size() ; Rj1++) {
						    	int Rj1_nuc = pos2nuc[j-1][Rj1];
						    	if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){continue;}
						    	if(DEPflg && Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;} // dependency between j-1 and j
					    		if(q == j - 1 && Rj1_nuc != Rq_nuc){ continue;} // q,jの時は、qの塩基とj-1の塩基は一致していないといけない。(2)

						    	for (unsigned int Lp1 = 0; Lp1 < pos2nuc[p-1].size() ; Lp1++) {
							    	int Lp1_nuc = pos2nuc[p-1][Lp1];
							    	if(NCflg == 1 && i2r[Lp1_nuc] != NucConst[p-1]){continue;}

							    	if(DEPflg && Dep1[ii2r[Lp1_nuc*10+Lp_nuc]][p-1] == 0){ continue;} // dependency between p-1 and p
							    	if(DEPflg && i == p-2 && Dep1[ii2r[Li_nuc*10+Lp1_nuc]][i] == 0){ continue; }  // i,X,p: dependency between i and p-1
							    	if(DEPflg && i == p-3 && Dep1[ii2r[Li1_nuc*10+Lp1_nuc]][i+1] == 0){ continue; }  // i,X,X,p: dependency between i+1 and p-1

							    	if(i == p-1 && Li_nuc != Lp1_nuc){ continue; }   // i,pの時は、iの塩基とp-1の塩基は一致していないといけない。(1)の逆
						    		if(i == p-2 && Li1_nuc != Lp1_nuc){ continue; }  // i,X,pの時は、i+1の塩基とp-1の塩基(X)は一致していないといけない。

							    	for (unsigned int Rq1 = 0; Rq1 < pos2nuc[q+1].size() ; Rq1++) {
								    	int Rq1_nuc = pos2nuc[q+1][Rq1];
								    	if(NCflg == 1 && i2r[Rq1_nuc] != NucConst[q+1]){continue;}
								    	if(DEPflg && Dep1[ii2r[Rq_nuc*10+Rq1_nuc]][q] == 0){ continue;} // dependency between q and q+1

								    	if(DEPflg && j == q+2 && Dep1[ii2r[Rq1_nuc*10+Rj_nuc]][q+1] == 0){ continue; }   // q,X,j: dependency between j and q-1
								    	if(DEPflg && j == q+3 && Dep1[ii2r[Rq1_nuc*10+Rj1_nuc]][q+1] == 0){ continue; }  // q,X,X,j: dependency between j+1 and q-1


								    	if(q+1 == j && Rq1_nuc != Rj_nuc){ continue;} // q,jの時は、q+1の塩基とjの塩基は一致していないといけない。(2)の逆
								    	if(q+2 == j && Rq1_nuc != Rj1_nuc){ continue;} // q,X,jの時は、q+1の塩基とj-1の塩基は一致していないといけない。


								    	//int energy = LoopEnergy(p-i-1, j-q-1, type_LiRj, type_LpRq,
								    	// 			i2r[Li1_nuc], i2r[Rj1_nuc], i2r[Lp1_nuc], i2r[Rq1_nuc]);
								    	int energy = E_intloop(p-i-1, j-q-1, type_LiRj, type_LpRq,
								    			i2r[Li1_nuc], i2r[Rj1_nuc], i2r[Lp1_nuc], i2r[Rq1_nuc], P);

								    	//	int energy_new = energy+c[indx[q]+p][Lp][Rq];
								    	int energy_new = energy+c[getIndx(p,q,w,indx)][Lp][Rq];
								    	traced = (cij == energy_new);
								    	if (traced) {
								    		base_pair[++b].i = p;
								    		base_pair[b].j   = q;

								    		(*optseq)[p] = i2n[Lp_nuc];
								    		(*optseq)[q] = i2n[Rq_nuc];

								    		(*optseq)[i+1] = i2n[Li1_nuc];
								    		(*optseq)[p-1] = i2n[Lp1_nuc];
								    		(*optseq)[j-1] = i2n[Rj1_nuc];
								    		(*optseq)[q+1] = i2n[Rq1_nuc];

								    		i = p, j = q; // i,jの更新
								    		Li = Lp;
								    		Rj = Rq;

								    		goto repeat1;
								    	}

								    }
							    }
					    	}
					    }
				    }
		    	}
		    }

	    }
	    /* end of repeat: --------------------------------------------------*/

	    /* (i.j) must close a multi-loop */

	    int rtype_LiRj = rtype[type_LiRj];
	    i1 = i+1; j1 = j-1;

	    sector[s+1].ml  = sector[s+2].ml = 1;

	    int en = cij - TermAU(rtype_LiRj, P) - P->MLintern[rtype_LiRj] - P->MLclosing;
	    //	    for(k = i+2+TURN; k < j-2-TURN; k++){
	    int Li1_save, Rk1_save, Lk_save, Rj1_save;
	    Li1_save = Rk1_save = Lk_save = Rj1_save = -1;
	    for(k = i+3+TURN; k < j-1-TURN; k++){
		    for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k-1].size(); Rk1++) {
		    	int Rk1_nuc = pos2nuc[k-1][Rk1];
		    	if(NCflg == 1 && i2r[Rk1_nuc] != NucConst[k-1]){continue;}
		    	// i,k,jでの塩基の矛盾はチェックする必要がない。なぜなら、i+4<k, k+4<jだから。

		    	for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
			    	int Lk_nuc = pos2nuc[k][Lk];
			    	if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){continue;}
			    	if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k-1 and k

			    	for (unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size(); Li1++) {
				    	int Li1_nuc = pos2nuc[i+1][Li1];
				    	if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){continue;}
				    	if(DEPflg && Dep1[ii2r[Li_nuc*10+Li1_nuc]][i] == 0){ continue;} // dependency between i and i+1

				    	for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size(); Rj1++) {
					    	int Rj1_nuc = pos2nuc[j-1][Rj1];
					    	if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){continue;}
					    	if(DEPflg && Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;} // dependency between j-1 and j

					    	//マルチループを閉じるところと、bifucationを同時に探している。
					    	//if(en == m[indx[k-1]+i+1][Li1][Rk1] + m[indx[j-1]+k][Lk][Rj1]){
					    	if(en == m[getIndx(i+1,k-1,w,indx)][Li1][Rk1] + m[getIndx(k,j-1,w,indx)][Lk][Rj1]){
					    		Li1_save = Li1;
					    		Rk1_save = Rk1;
					    		Lk_save = Lk;
					    		Rj1_save = Rj1;
					    		goto LABEL2;
					    	}
					    }
				    }
			    }
		    }
	    }
	    LABEL2:

	    if (k<=j-2-TURN) { /* found the decomposition successfully*/
	    	sector[++s].i = i1;
	    	sector[s].j   = k-1;
	    	sector[s].Li  = Li1_save;
	    	sector[s].Rj  = Rk1_save;

	    	sector[++s].i = k;
	    	sector[s].j   = j1;
	    	sector[s].Li  = Lk_save;
	    	sector[s].Rj  = Rj1_save;

	    } else {
	    	fprintf(stderr, "backtracking failed in repeat %d %d\n", i , j);
//	    	exit(1);
	    }

	}

	base_pair[0].i = b;    /* save the total number of base pairs */
//	cout << base_pair[0].i << endl;

}

void backtrack2(string *optseq, stack *sector, bond *base_pair, int ***const &c, int ***const &m, int*** const &f2,
			int *const indx, const int &initL, const int &initR, paramT *const&P, const vector<int> &NucConst,
			const vector<vector <int> > &pos2nuc, const int &NCflg, int *const &i2r, int const &length, int const &w,
			int const (&BP_pair)[5][5], char * const &i2n, int * const &rtype, int *const &ii2r,
			vector<vector<int> > &Dep1, vector<vector<int> > &Dep2, int &DEPflg,
			vector<vector<vector<vector<pair<int, string> > > > > &predefH, map<string, int> &predefE, vector<vector<vector<string> > > &substr, map<char, int> &n2i, const char* nucdef){

	InitRand();

	int s = 0;
	int b = 0;
	sector[++s].i = 1;
	sector[s].j = length;
	sector[s].Li = initL;
	sector[s].Rj = initR;
	sector[s].ml = 0;


	OUTLOOP:
	while (s>0) {
		//	    int fij, fi, ij, cij, traced, traced_Lk, i1, j1, k, p , q;
		int fij, fi, ij, cij, traced, i1, j1, k, p , q;

	    //ここでoptseqに値を反映させるべき？
	    int i  = sector[s].i;
	    int j  = sector[s].j;
	    int Li  = sector[s].Li;
	    int Rj  = sector[s].Rj;
	    int ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
	                              	  occur in the M- (1) or in the F-array (0) */
	    ij = getIndx(i,j,w, indx);
	    int Li_nuc = pos2nuc[i][Li];
	    int Rj_nuc = pos2nuc[j][Rj];

	    if(i + 1 == j && Dep1[ii2r[Li_nuc*10+Rj_nuc]][i] == 0){ continue;}
	    if(i + 2 == j && Dep2[ii2r[Li_nuc*10+Rj_nuc]][i] == 0){ continue;}

	    int type_LiRj = BP_pair[i2r[Li_nuc]][i2r[Rj_nuc]];

	    (*optseq)[i] = i2n[Li_nuc];
	    (*optseq)[j] = i2n[Rj_nuc];

	    if (ml==2) {
	      base_pair[++b].i = i;
	      base_pair[b].j   = j;
	      goto repeat1;
	    }

	    //	    if (j < i+TURN+1) continue; /* no more pairs in this interval */
	    if (j == i + 1) continue;

	    //	    fij = (ml == 1)? m[indx[j]+i][Li][Rj] : f[j][Li][Rj];
	    fij = (ml == 1)? m[getIndx(i,j,w,indx)][Li][Rj] : f2[ij][Li][Rj];
	    cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << fij << ")" << ":" << *optseq << ":" << s << endl;


	    // trace i,j from i,j-1 for multi-loop
	    if(ml == 1){
		    for(unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size(); Rj1++){
		    	int Rj1_nuc = pos2nuc[j-1][Rj1];
			    if(Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;}
			    int mi = m[getIndx(i,j-1,w,indx)][Li][Rj1] + P->MLbase;

			    if (fij == mi) {  /* 3' end is unpaired */
			    	sector[++s].i = i;
			    	sector[s].j   = j-1;
			    	sector[s].Li   = Li;
			    	sector[s].Rj   = Rj1;
			    	sector[s].ml  = ml;
			    	//continue;
			    	goto OUTLOOP;
			    }
		    }
	    }

	    if (ml == 0) { /* backtrack in f */
	    	//traced = 0;
	    	int label[4];
	    	for(int l = 0; l < 4; l++){
	    		label[l] = l;
	    	}
	    	shuffle(label, 4);

	    	for(int l = 0; l < 4; l++){
	    		cout << "go to label F" << label[l]+1 << endl;
	    		if(label[l] == 0){
	    			goto F1;
	    		}
	    		else if(label[l] == 1){
	    			goto F2;
	    		}
	    		else if(label[l] == 2){
	    			goto F3;
	    		}
	    		else{
	    			goto F4;
	    		}

	    		F1:
	    		// trace i,j from i,j-1
	    		for(unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size(); Rj1++){
	    			int Rj1_nuc = pos2nuc[j-1][Rj1];
	    			if(Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;}

	    			fi  = f2[getIndx(i,j-1,w,indx)][Li][Rj1];

	    			if (fij == fi) {  /* 3' end is unpaired */
	    				sector[++s].i = i;
	    				sector[s].j   = j-1;
	    				sector[s].Li   = Li;
	    				sector[s].Rj   = Rj1;
	    				sector[s].ml  = ml;
	    				//continue;
                		cout << "Traceback path found." << endl;
	    				goto OUTLOOP;
	    			}

	    		}
	    		continue;

	    		F2:
	    		// trace i,j from i+1,j
	    		for(unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size(); Li1++){
	    			int Li1_nuc = pos2nuc[i+1][Li1];
	    			if(Dep1[ii2r[Li_nuc*10+Li1_nuc]][i+1] == 0){ continue;}

	    			fi  = f2[getIndx(i+1,j,w,indx)][Li1][Rj];

	    			if (fij == fi) {  /* 5' end is unpaired */
	    				sector[++s].i = i+1;
	    				sector[s].j   = j;
	    				sector[s].Li   = Li1;
	    				sector[s].Rj   = Rj;
	    				sector[s].ml  = ml;
	    				//continue;
                		cout << "Traceback path found." << endl;
	    				goto OUTLOOP;
	    			}

	    		}
	    		continue;

	    		F3:
	    		// trace i,j from C(i,j)
	    		if(type_LiRj && j - i + 1 <= w){
	    			int en_c = TermAU(type_LiRj, P) +  c[getIndx(i,j,w,indx)][Li][Rj];
	    			int en_f = f2[ij][Li][Rj];
	    			cout << en_c << "," << en_f << endl;
	    			if(en_c ==  en_f){
//	    			k = i;
//	    			traced = j;
	    				base_pair[++b].i = i;
	    				base_pair[b].j   = j;
//	    			traced_Lk = Li;
	    	    	//	    			goto LABEL1;
                		cout << "Traceback path found." << endl;
	    				goto repeat1;
	    			}
	    		}
	    		continue;

	    		F4:
	    		//	    		for(k=j-TURN-1; k>= i + TURN + 2; k--){
	    		//	    		for(k=j-1; k>= i + 2; k--){
	    		for(k=j-1; k>= i + 2; k--){

	    			for(unsigned int Rk1 = 0; Rk1 < pos2nuc[k-1].size(); Rk1++){
	    				int Rk1_nuc = pos2nuc[k-1][Rk1];
	    				if(DEPflg && k == 3 && Dep1[ii2r[Li_nuc*10+Rk1_nuc]][i] == 0){ continue;} // dependency between 1(i) and 2(k-1)
	    				if(DEPflg && k == 4 && Dep2[ii2r[Li_nuc*10+Rk1_nuc]][i] == 0){ continue;} // dependency between 1(i) and 3(k-1)

//	    				if((k - 1) - i + 1 > w ||
//	    					j - k + 1 > w)
//	    						continue;

	    				for(unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++){
	    					int Lk_nuc = pos2nuc[k][Lk];
	    					//if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){continue;}
	    					if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k-1 and k

	    					int en_f1 = f2[getIndx(i,k-1,w,indx)][Li][Rk1];
	    					int en_f2 = f2[getIndx(k,j,w,indx)][Lk][Rj];
	    					if(fij ==  en_f1 + en_f2){
//	                    		traced = j;
//	                    		traced_Lk = Lk;
	                    		/* push back the remaining f portion */
	                    		sector[++s].i = i;
	                    		sector[s].j   = k-1;
	                    		sector[s].Li = Li;
	                    		sector[s].Rj = Rk1;
	                    		sector[s].ml  = ml;
	                    		sector[++s].i = k;
	                    		sector[s].j   = j;
	                    		sector[s].Li = Lk;
	                    		sector[s].Rj = Rj;
	                    		sector[s].ml  = ml;
	                    		cout << "Traceback path found in " << i << "," << k-1 << " and " << k << "," << j << endl;
	        	    			goto OUTLOOP;
	    					}
	    				}
	    			}
	    		}
	    		continue;

	    	}
//	    	LABEL1:
//	    	if (!traced){
	    		fprintf(stderr, "backtrack failed in f2\n");
	    		fprintf(stderr, "cannot trace f2[%d][%d][%d][%d] Lnuc=%c Rnuc=%c \n", i, j, Li, Rj, i2n[Li_nuc], i2n[Rj_nuc]);
	    		exit(0);
//	    	}


	    	/* trace back the base pair found */
	    	/*
	    	// [1]
	    	i=k;            // iを更新
	    	j=traced;       // この代入は多分必要ない。jに関しては不変だから
	    	Li = traced_Lk; // Liを更新
	    	//Rjは不変
	    	base_pair[++b].i = i;
	    	base_pair[b].j   = j;
	    	goto repeat1;
	    	*/
	    }
	    else { /* trace back in fML array */

		    for(unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size(); Li1++){
		    	int Li1_nuc = pos2nuc[i+1][Li1];
		    	if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){continue;}
    		    if(DEPflg && Dep1[ii2r[Li_nuc*10+Li1_nuc]][i] == 0){ continue;} // dependency between k-1 and k

    		    //	  	if (m[indx[j]+i+1][Li1][Rj]+P->MLbase == fij) { /* 5' end is unpaired */
    		    if (m[getIndx(i+1,j,w,indx)][Li1][Rj]+P->MLbase == fij) { /* 5' end is unpaired */
	    			sector[++s].i = i+1;
	    			sector[s].j   = j;
	    			sector[s].Li  = Li1;
	    			sector[s].Rj  = Rj;
	    			sector[s].ml  = ml;
	    			goto OUTLOOP;
		    	}

    		    //	ij  = indx[j]+i;
    		    ij  = getIndx(i,j,w,indx);

		    	if(fij == c[ij][Li][Rj] + TermAU(type_LiRj, P) + P->MLintern[type_LiRj]){
		    		base_pair[++b].i = i;
		    		base_pair[b].j   = j;
		    		goto repeat1;
		    	}

		    }

		    //		    for(k = i + 1 + TURN; k <= j - 2 - TURN; k++){
		    for(k = i + 2 + TURN; k <= j - 1 - TURN; k++){

	    		for(unsigned int Rk1 = 0; Rk1 < pos2nuc[k-1].size(); Rk1++){
	    	    	int Rk1_nuc = pos2nuc[k-1][Rk1];
	    	    	if(NCflg == 1 && i2r[Rk1_nuc] != NucConst[k-1]){continue;}
	    	    	// check dependency is not needed because i+2<k,k+2<J
	    	    	for(unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++){
		    			int Lk_nuc = pos2nuc[k][Lk];
		    	    	if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){continue;}
		    		    if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k-1 and k

		    		    //if(fij == (m[indx[k-1]+i][Li][Rk1]+m[indx[j]+k][Lk][Rj])){
		    		    if(fij == (m[getIndx(i,k-1,w,indx)][Li][Rk1]+m[getIndx(k,j,w,indx)][Lk][Rj])){
		    	    		sector[++s].i = i;
		    	    		sector[s].j   = k-1;
		    	    		sector[s].Li  = Li;
		    	    		sector[s].Rj  = Rk1;
		    	    		sector[s].ml  = ml;
		    	    		sector[++s].i = k;
		    	    		sector[s].j   = j;
		    	    		sector[s].Li  = Lk;
		    	    		sector[s].Rj  = Rj;
		    	    		sector[s].ml  = ml;
		    	    		goto OUTLOOP;
		    	    	}
		    		}
	    		}

		    }
    		if (k>j-1-TURN){
    			fprintf(stderr, "backtrack failed in fML\n");
    			exit(1);
    		}
	    }

	    repeat1: // いちいちスタックに積まずに、ここで部分的なトレースバックをしてしまう。
	    //continue;
	    /*----- begin of "repeat:" -----*/
	    if(j - i + 1 > w){
	    	cerr << "backtrack failed at " << i << "," << j << " : the length must at most << w << endl";
	    }
	    //	    ij = indx[j]+i; // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
	    ij = getIndx(i,j,w,indx); // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
	    Li_nuc = pos2nuc[i][Li]; // Liは更新されている。
	    Rj_nuc = pos2nuc[j][Rj]; //Rj_nucは更新されていない場合もある[1]。
	    type_LiRj = BP_pair[i2r[Li_nuc]][i2r[Rj_nuc]];
	    cij = c[ij][Li][Rj];
		(*optseq)[i] = i2n[Li_nuc]; //塩基対部分を記録
		(*optseq)[j] = i2n[Rj_nuc];

//		if(TB_CHK_flg == 1)
//			cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << cij << ")" << ":" << *optseq << endl;


	    // predefinedなヘアピンのトレースバック
//		if ((j-i+1 == 5 ||j-i+1 == 6 ||j-i+1 == 8) &&
 //   			cij == predefH[i][j-i+1][i2r[Li_nuc]][i2r[Rj_nuc]].first){
//			string s1 = predefH[i][j-i+1][i2r[Li_nuc]][i2r[Rj_nuc]].second;
//			for(unsigned int k = 0; k < s1.size(); k++){
//				(*optseq)[i+k] = s1[k]; //塩基を記録
//			}
//			cout << "Predefined Hairpin " << s1 << " at " << i << ":" << j << endl;
//			goto OUTLOOP;
 //   	}
		if (j-i+1 == 5 ||j-i+1 == 6 ||j-i+1 == 8){
			int l = j-i+1;
			for(unsigned int s = 0; s < substr[i][l].size(); s++){
				string hpn = substr[i][l][s];
				int hL_nuc  = n2i[hpn[0]];
				int hL2_nuc = n2i[hpn[1]];
				int hR2_nuc = n2i[hpn[l-2]];
				int hR_nuc  = n2i[hpn[l-1]];
				if(hL_nuc != i2r[Li_nuc]) continue;
				if(hR_nuc != i2r[Rj_nuc]) continue;

				if(NCflg == 1){
					string s1 = string(nucdef).substr(i, l);
					if(hpn != s1) continue;
				}


				if(DEPflg && Li_nuc > 4 && Dep1[ii2r[Li_nuc*10+hL2_nuc]][i] == 0){continue;} // Dependency is already checked.
				if(DEPflg && Rj_nuc > 4 && Dep1[ii2r[hR2_nuc*10+Rj_nuc]][j-1] == 0){continue;}// ただし、Li_nuc、Rj_nucがVWXYのときだけは、一つ内側との依存関係をチェックする必要がある。

				// predefinedなヘアピンとの比較
				if(predefE.count(hpn) > 0){
					if(c[ij][Li][Rj] == predefE[hpn]){
						cout << "Predefined Hairpin at " << i << "," << j << endl;
						for(unsigned int k = 0; k < hpn.size(); k++){
							(*optseq)[i+k] = hpn[k]; //塩基を記録
						}
						goto OUTLOOP;
					}

				}
				else{	// 普通のヘアピンとの比較
					//					int energy = HairpinE(j - i - 1, type_LiRj,
					//							i2r[hL2_nuc], i2r[hR2_nuc],
					//							"NNNNNNNNN");
					int energy = E_hairpin(j - i - 1, type_LiRj,
												i2r[hL2_nuc], i2r[hR2_nuc],
												"NNNNNNNNN", P);

					if(c[ij][Li][Rj] == energy){
						for(unsigned int k = 0; k < hpn.size(); k++){
							(*optseq)[i+k] = hpn[k]; //塩基を記録
						}
						goto OUTLOOP;
					}
				}
			}

		}
		else{
			// 普通のヘアピンのトレースバック
			for(unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size(); Li1++){
				int Li1_nuc = pos2nuc[i+1][Li1];
				if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){continue;}
				if(DEPflg && Dep1[ii2r[Li_nuc*10+Li1_nuc]][i] == 0){ continue;} // dependency between i and i+1

				for(unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size(); Rj1++){
					int Rj1_nuc = pos2nuc[j-1][Rj1];
					if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){continue;}
					if(DEPflg && Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;} // dependency between j-1 and j

					//if (cij == HairpinE(j-i-1, type_LiRj, i2r[Li1_nuc], i2r[Rj1_nuc], "NNNNNNNNN")){
					if (cij == E_hairpin(j-i-1, type_LiRj, i2r[Li1_nuc], i2r[Rj1_nuc], "NNNNNNNNN", P)){
						(*optseq)[i+1] = i2n[Li1_nuc]; //塩基対の内側のミスマッチ塩基を記録
						(*optseq)[j-1] = i2n[Rj1_nuc];
						goto OUTLOOP;
					}
					else{
						continue;
					}
				}
			}
		}
	    // Hairpinに該当がなければ、Internal loopのトレースバック。もっとも手強い。
	    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
		    for (unsigned int Lp = 0; Lp < pos2nuc[p].size() ; Lp++) {
		    	int Lp_nuc = pos2nuc[p][Lp];
		    	if(NCflg == 1 && i2r[Lp_nuc] != NucConst[p]){continue;}
		    	if(DEPflg && p == i + 1 && Dep1[ii2r[Li_nuc*10+Lp_nuc]][i] == 0){ continue;} // dependency between i and q
		    	if(DEPflg && p == i + 2 && Dep2[ii2r[Li_nuc*10+Lp_nuc]][i] == 0){ continue;} // dependency between i and q

		    	int minq = j-i+p-MAXLOOP-2;
		    	if (minq<p+1+TURN) minq = p+1+TURN;
		    	for (q = j-1; q >= minq; q--) {
				    for (unsigned int Rq = 0; Rq < pos2nuc[q].size() ; Rq++) {
				    	int Rq_nuc = pos2nuc[q][Rq];
				    	if(NCflg == 1 && i2r[Rq_nuc] != NucConst[q]){continue;}
				    	if(DEPflg && q == j - 1 && Dep1[ii2r[Rq_nuc*10+Rj_nuc]][q] == 0){ continue;} // dependency between q and j
				    	if(DEPflg && q == j - 2 && Dep2[ii2r[Rq_nuc*10+Rj_nuc]][q] == 0){ continue;} // dependency between q and j

				    	int type_LpRq = BP_pair[i2r[Lp_nuc]][i2r[Rq_nuc]];
				    	if (type_LpRq==0) continue;
				    	type_LpRq = rtype[type_LpRq];

					    for (unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size() ; Li1++) {
					    	int Li1_nuc = pos2nuc[i+1][Li1];
					    	if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){continue;}
					    	if(DEPflg && Dep1[ii2r[Li_nuc*10+Li1_nuc]][i] == 0){ continue;} // dependency between i and i+1
					    	if(i+1 == p && Li1_nuc != Lp_nuc){ continue; } // i,pの時は、i+1の塩基とpの塩基は一致していないといけない。(1)

				    		for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size() ; Rj1++) {
						    	int Rj1_nuc = pos2nuc[j-1][Rj1];
						    	if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){continue;}
						    	if(DEPflg && Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;} // dependency between j-1 and j
					    		if(q == j - 1 && Rj1_nuc != Rq_nuc){ continue;} // q,jの時は、qの塩基とj-1の塩基は一致していないといけない。(2)

						    	for (unsigned int Lp1 = 0; Lp1 < pos2nuc[p-1].size() ; Lp1++) {
							    	int Lp1_nuc = pos2nuc[p-1][Lp1];
							    	if(NCflg == 1 && i2r[Lp1_nuc] != NucConst[p-1]){continue;}

							    	if(DEPflg && Dep1[ii2r[Lp1_nuc*10+Lp_nuc]][p-1] == 0){ continue;} // dependency between p-1 and p
							    	if(DEPflg && i == p-2 && Dep1[ii2r[Li_nuc*10+Lp1_nuc]][i] == 0){ continue; }  // i,X,p: dependency between i and p-1
							    	if(DEPflg && i == p-3 && Dep1[ii2r[Li1_nuc*10+Lp1_nuc]][i+1] == 0){ continue; }  // i,X,X,p: dependency between i+1 and p-1

							    	if(i == p-1 && Li_nuc != Lp1_nuc){ continue; }   // i,pの時は、iの塩基とp-1の塩基は一致していないといけない。(1)の逆
						    		if(i == p-2 && Li1_nuc != Lp1_nuc){ continue; }  // i,X,pの時は、i+1の塩基とp-1の塩基(X)は一致していないといけない。

							    	for (unsigned int Rq1 = 0; Rq1 < pos2nuc[q+1].size() ; Rq1++) {
								    	int Rq1_nuc = pos2nuc[q+1][Rq1];
								    	if(NCflg == 1 && i2r[Rq1_nuc] != NucConst[q+1]){continue;}
								    	if(DEPflg && Dep1[ii2r[Rq_nuc*10+Rq1_nuc]][q] == 0){ continue;} // dependency between q and q+1

								    	if(DEPflg && j == q+2 && Dep1[ii2r[Rq1_nuc*10+Rj_nuc]][q+1] == 0){ continue; }   // q,X,j: dependency between j and q-1
								    	if(DEPflg && j == q+3 && Dep1[ii2r[Rq1_nuc*10+Rj1_nuc]][q+1] == 0){ continue; }  // q,X,X,j: dependency between j+1 and q-1


								    	if(q+1 == j && Rq1_nuc != Rj_nuc){ continue;} // q,jの時は、q+1の塩基とjの塩基は一致していないといけない。(2)の逆
								    	if(q+2 == j && Rq1_nuc != Rj1_nuc){ continue;} // q,X,jの時は、q+1の塩基とj-1の塩基は一致していないといけない。


								    	//int energy = LoopEnergy(p-i-1, j-q-1, type_LiRj, type_LpRq,
								    	// 			i2r[Li1_nuc], i2r[Rj1_nuc], i2r[Lp1_nuc], i2r[Rq1_nuc]);
								    	int energy = E_intloop(p-i-1, j-q-1, type_LiRj, type_LpRq,
								    			i2r[Li1_nuc], i2r[Rj1_nuc], i2r[Lp1_nuc], i2r[Rq1_nuc], P);

								    	//	int energy_new = energy+c[indx[q]+p][Lp][Rq];
								    	int energy_new = energy+c[getIndx(p,q,w,indx)][Lp][Rq];
								    	traced = (cij == energy_new);
								    	if (traced) {
								    		base_pair[++b].i = p;
								    		base_pair[b].j   = q;

								    		(*optseq)[p] = i2n[Lp_nuc];
								    		(*optseq)[q] = i2n[Rq_nuc];

								    		(*optseq)[i+1] = i2n[Li1_nuc];
								    		(*optseq)[p-1] = i2n[Lp1_nuc];
								    		(*optseq)[j-1] = i2n[Rj1_nuc];
								    		(*optseq)[q+1] = i2n[Rq1_nuc];

								    		i = p, j = q; // i,jの更新
								    		Li = Lp;
								    		Rj = Rq;

								    		goto repeat1;
								    	}

								    }
							    }
					    	}
					    }
				    }
		    	}
		    }

	    }
	    /* end of repeat: --------------------------------------------------*/

	    /* (i.j) must close a multi-loop */

	    int rtype_LiRj = rtype[type_LiRj];
	    i1 = i+1; j1 = j-1;

	    sector[s+1].ml  = sector[s+2].ml = 1;

	    int en = cij - TermAU(rtype_LiRj, P) - P->MLintern[rtype_LiRj] - P->MLclosing;
	    //	    for(k = i+2+TURN; k < j-2-TURN; k++){
	    int Li1_save, Rk1_save, Lk_save, Rj1_save;
	    Li1_save = Rk1_save = Lk_save = Rj1_save = -1;
	    for(k = i+3+TURN; k < j-1-TURN; k++){
		    for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k-1].size(); Rk1++) {
		    	int Rk1_nuc = pos2nuc[k-1][Rk1];
		    	if(NCflg == 1 && i2r[Rk1_nuc] != NucConst[k-1]){continue;}
		    	// i,k,jでの塩基の矛盾はチェックする必要がない。なぜなら、i+4<k, k+4<jだから。

		    	for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
			    	int Lk_nuc = pos2nuc[k][Lk];
			    	if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){continue;}
			    	if(DEPflg && Dep1[ii2r[Rk1_nuc*10+Lk_nuc]][k-1] == 0){ continue;} // dependency between k-1 and k

			    	for (unsigned int Li1 = 0; Li1 < pos2nuc[i+1].size(); Li1++) {
				    	int Li1_nuc = pos2nuc[i+1][Li1];
				    	if(NCflg == 1 && i2r[Li1_nuc] != NucConst[i+1]){continue;}
				    	if(DEPflg && Dep1[ii2r[Li_nuc*10+Li1_nuc]][i] == 0){ continue;} // dependency between i and i+1

				    	for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j-1].size(); Rj1++) {
					    	int Rj1_nuc = pos2nuc[j-1][Rj1];
					    	if(NCflg == 1 && i2r[Rj1_nuc] != NucConst[j-1]){continue;}
					    	if(DEPflg && Dep1[ii2r[Rj1_nuc*10+Rj_nuc]][j-1] == 0){ continue;} // dependency between j-1 and j

					    	//マルチループを閉じるところと、bifucationを同時に探している。
					    	//if(en == m[indx[k-1]+i+1][Li1][Rk1] + m[indx[j-1]+k][Lk][Rj1]){
					    	if(en == m[getIndx(i+1,k-1,w,indx)][Li1][Rk1] + m[getIndx(k,j-1,w,indx)][Lk][Rj1]){
					    		Li1_save = Li1;
					    		Rk1_save = Rk1;
					    		Lk_save = Lk;
					    		Rj1_save = Rj1;
					    		goto LABEL2;
					    	}
					    }
				    }
			    }
		    }
	    }
	    LABEL2:

	    if (k<=j-2-TURN) { /* found the decomposition successfully*/
	    	sector[++s].i = i1;
	    	sector[s].j   = k-1;
	    	sector[s].Li  = Li1_save;
	    	sector[s].Rj  = Rk1_save;

	    	sector[++s].i = k;
	    	sector[s].j   = j1;
	    	sector[s].Li  = Lk_save;
	    	sector[s].Rj  = Rj1_save;

	    } else {
	    	fprintf(stderr, "backtracking failed in repeat %d %d\n", i , j);
//	    	exit(1);
	    }

	}

	base_pair[0].i = b;    /* save the total number of base pairs */
//	cout << base_pair[0].i << endl;

}

/*
string init_string(int const &len){
	string s;
	s.resize(len+1);
	s[0] = ' ';
	for(int i = 1; i <= len; i++){
		s[i] = 'N';
	}
	return s;
}
*/

inline int TermAU(int const &type, paramT * const &P){
	if(type>2){
		return P->TerminalAU;
	}
	return 0;
}

inline int E_hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P){
	  int energy;
	  //fprintf(stderr, "ok\n");
	  energy = (size <= 30) ? P->hairpin[size] : P->hairpin[30]+(int)(P->lxc*log((size)/30.));
	  if (P->model_details.special_hp){
	    if (size == 4) { /* check for tetraloop bonus */
	      char tl[7]={0}, *ts;
	      strncpy(tl, string, 6);
	      if ((ts=strstr(P->Tetraloops, tl)))
	        return (P->Tetraloop_E[(ts - P->Tetraloops)/7]);
	    }
	    else if (size == 6) {
	      char tl[9]={0}, *ts;
	      strncpy(tl, string, 8);
	      if ((ts=strstr(P->Hexaloops, tl)))
	        return (energy = P->Hexaloop_E[(ts - P->Hexaloops)/9]);
	    }
	    else if (size == 3) {
	      char tl[6]={0,0,0,0,0,0}, *ts;
	      strncpy(tl, string, 5);
	      if ((ts=strstr(P->Triloops, tl))) {
	        return (P->Triloop_E[(ts - P->Triloops)/6]);
	      }
	      return (energy + (type>2 ? P->TerminalAU : 0));
	    }
	  }
	  energy += P->mismatchH[type][si1][sj1];

	  return energy;
	}


inline int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P){
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns, energy;
  energy = INF;
  int MAX_NINIO=300;
  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}

  if (nl == 0)
    return P->stack[type][type_2];  /* stack */

  if (ns==0) {                      /* bulge */
    energy = (nl<=MAXLOOP)?P->bulge[nl]:
      (P->bulge[30]+(int)(P->lxc*log(nl/30.)));
    if (nl==1) energy += P->stack[type][type_2];
    else {
      if (type>2) energy += P->TerminalAU;
      if (type_2>2) energy += P->TerminalAU;
    }
    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1)                    /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1];
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];
        return energy;
      }
      else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(P->internal_loop[nl+1]) : (P->internal_loop[30]+(int)(P->lxc*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);
        energy += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1];}
      else if (nl==3){              /* 2x3 loop */
        energy = P->internal_loop[5]+P->ninio[2];
        energy += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]) : (P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));
      energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}

int getMemoryUsage(const string &fname){
	//cout << fname << endl;
	 ifstream ifs(fname.c_str());
	 if(ifs){
	    string line;
	    while(getline(ifs, line)){ // $B9T$NFI$_9~$_(B
	    	string index = line.substr(0,6);
	    	if(index == "VmRSS:"){
	    		//cout << line << endl;
	    		string mem_str;
	    		for(unsigned int i = 6; i < line.size(); i++){
	    			if(line[i] == ' ' || line[i] == '\t') continue;
	    			if(isdigit(line[i])){
	    				mem_str.append(1, line[i]);
	    			}
	    			else if(line[i] == 'k' || line[i] == 'B'  ){
	    				continue;
	    			}
	    			else{
	    				cerr << "Unexpected letter found in " << fname << "(" << line[i] << ")" << endl;
	    				return -1;
	    			}
	    		}
	    		stringstream ss(mem_str);
	    		int val;
	    		ss >> val;
	    		return val;
	    	}
	    }

	 }
	 else{
	    cerr << "Error: cannot open file(" << fname << ")" << endl;
		return -1;
	 }

	 cerr << "VmRSS line was not found in " << fname << endl;
	 return -1;
}

void fill_optseq(string *optseq, int I, int J, vector <vector<int> > &pos2nuc, vector <vector<int> > Dep1){

	int i2r[20], ii2r[100];
	map<char, int> n2i = make_n2i();
	char i2n[20];
	make_i2n(i2n);
	make_i2r(i2r);
	make_ii2r(ii2r); // encoding a dinucleotide (each eight variation) to an integer from 11 - 88

	// main routine
//	cout << I << endl;
//	cout << pos2nuc[I][0]<< endl;
//	cout << i2n[pos2nuc[I][0]] << endl;
	(*optseq)[I] = i2n[pos2nuc[I][0]];
//	cout << "CHK:" << (*optseq)[I] << endl;
	for (int i = I+1; i <= J; i++) {
		for (unsigned int L = 0; L < pos2nuc[i].size(); L++) { // search possible nucleotides
			int L_nuc = pos2nuc[i][L];
			int L1_nuc;
			L1_nuc = n2i[(*optseq)[i-1]];
			if(Dep1[ii2r[L1_nuc*10+L_nuc]][i-1] == 0){continue;}
			(*optseq)[i] =i2n[L_nuc];
			break;
		}
	}

	//塩基V,W,X,Yの修正
	for (int i = I; i <= J; i++) {
		if((*optseq)[i] == 'V' || (*optseq)[i] == 'W'){
			(*optseq)[i] = 'U';
		}
		else if((*optseq)[i] == 'X' || (*optseq)[i] == 'Y'){
			(*optseq)[i] = 'G';
		}
	}

}

void fixed_init_matrix(const int &nuclen, const int &size, int *C, int *M, int *F, int *DMl, int *DMl1, int *DMl2){
	for(int i = 0; i <= nuclen; i++){
		F[i] = 0;
		DMl[i] = INF;
		DMl1[i] = INF;
		DMl2[i] = INF;
	}

	for(int i = 0; i < size; i++){
		C[i] = INF;
		M[i] = INF;
	}
}

void fixed_backtrack(string optseq, bond *base_pair, int *c, int *m, int *f,
		int *indx, paramT *P, int nuclen, int w, const int (&BP_pair)[5][5], map<string, int> predefE){
	int rtype[7] = { 0, 2, 1, 4, 3, 6, 5 };
	int s = 0;
	int b = 0;
	stack sector[500];
	sector[++s].i = 1;
	sector[s].j = nuclen;
	sector[s].ml = 0;

	map<char, int> n2i = make_n2i();
	int ioptseq[nuclen+1];
	ioptseq[0] = 0;
	for(int i = 1; i <= nuclen; i++){
		ioptseq[i] = n2i[optseq[i]];
	}

	OUTLOOP:
	while (s>0) {
	    int fij, fi, ij, cij, traced, i1, j1, k, p , q;

	    int i  = sector[s].i;
	    int j  = sector[s].j;
	    int ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
	                              	  occur in the M- (1) or in the F-array (0) */

	    int type = BP_pair[ioptseq[i]][ioptseq[j]];

	    if (ml==2) {
	      base_pair[++b].i = i;
	      base_pair[b].j   = j;
	      goto repeat1;
	    }

	    if (j == i) break;

	    fij = (ml == 1)? m[getIndx(i,j,w,indx)] : f[j];
	    cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << fij << ")" << endl;

	    fi  = (ml == 1)? m[getIndx(i,j-1,w,indx)] + P->MLbase: f[j-1];

	    if (fij == fi) {  /* 3' end is unpaired */
	    	sector[++s].i = i;
	    	sector[s].j   = j-1;
	    	sector[s].ml  = ml;
	    	//continue;
	    	goto OUTLOOP;
	    }

	    if (ml == 0) { /* backtrack in f */

	    	if(i != 1){
    			cerr << "Traceback failure: i must be 1 during bachtrack in f" << endl;
	    	}

	    	//f[j]とC[1][j]が一致している時の処理。Vieenaでは、次for文に統合されている。
	    	if(type && j <= w){
	    		int en_c = TermAU(type, P) +  c[getIndx(i,j,w,indx)];
                int en_f = f[j];
              	if(en_c ==  en_f){
              		k = i;
              		traced = j;
               		goto LABEL1;
              	}
	    	}

	    	for(k=j-TURN-1,traced=0; k>=MAX2(2,j-w+1); k--){

	    		int type_kj = BP_pair[ioptseq[k]][ioptseq[j]];
	    		if(type_kj){
	    			int en_c = TermAU(type_kj, P) +  c[getIndx(k,j,w,indx)];
	    			int en_f = f[k-1];
	    			if(fij ==  en_c + en_f){
	    				traced = j;
	    				/* push back the remaining f portion */
	    				sector[++s].i = i;
	    				sector[s].j   = k-1;
	    				sector[s].ml  = 0;

	    				goto LABEL1;
	    			}

	    		}
	    	}
	    	LABEL1:

	    	if (!traced){
	    		fprintf(stderr, "backtrack failed in f\n");
	    		fprintf(stderr, "cannot trace f[%d] \n", j);
	    		exit(0);
	    	}

	    	/* trace back the base pair found */
	    	// [1]
	    	i=k;            // iを更新
	    	j=traced;       // この代入は多分必要ない。jに関しては不変だから
	    	base_pair[++b].i = i;
	    	base_pair[b].j   = j;
	    	goto repeat1;
	    }
	    else { /* trace back in fML array */

	    	if (m[getIndx(i+1,j,w,indx)]+P->MLbase == fij) { /* 5' end is unpaired */
	    		sector[++s].i = i+1;
	    		sector[s].j   = j;
	    		sector[s].ml  = ml;
	    		goto OUTLOOP;
	    	}

	    	ij  = getIndx(i,j,w,indx);

	    	if(fij == c[ij] + TermAU(type, P) + P->MLintern[type]){
	    		base_pair[++b].i = i;
	    		base_pair[b].j   = j;
	    		goto repeat1;
	    	}

	    	for(k = i + 2 + TURN; k <= j - 1 - TURN; k++){
	    		if(fij == (m[getIndx(i,k-1,w,indx)]+m[getIndx(k,j,w,indx)])){
	    			sector[++s].i = i;
	    			sector[s].j   = k-1;
	    			sector[s].ml  = ml;
	    			sector[++s].i = k;
	    			sector[s].j   = j;
	    			sector[s].ml  = ml;
	    			goto OUTLOOP;
	    		}
	    	}

	    	if (k>j-1-TURN){
    			fprintf(stderr, "backtrack failed in fML\n");
    			exit(1);
    		}
	    }

	    repeat1: // いちいちスタックに積まずに、ここで部分的なトレースバックをしてしまう。
	    //continue;
	    /*----- begin of "repeat:" -----*/
	    if(j - i + 1 > w){
	    	cerr << "backtrack failed at " << i << "," << j << " : the length must at most << w << endl";
	    }
	    ij = getIndx(i,j,w,indx); // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
	    type = BP_pair[ioptseq[i]][ioptseq[j]];
	    cij = c[ij];

	    if (j-i+1 == 5 ||j-i+1 == 6 ||j-i+1 == 8){
			string hpn = "";
			for(int k = i; k <= j; k++){
				hpn += optseq[k];
			}

			// predefinedなヘアピンとの比較
			if(predefE.count(hpn) > 0){
				if(c[ij] == predefE[hpn]){
					cout << "Predefined Hairpin at " << i << "," << j << endl;
					goto OUTLOOP;
				}
			}
			else{	// 普通のヘアピンとの比較
				int energy = E_hairpin(j - i - 1, type,
						ioptseq[i+1], ioptseq[j-1],
						"NNNNNNNNN", P);

				if(c[ij] == energy){
					goto OUTLOOP;
				}
			}

		}
		else{
			// 普通のヘアピンのトレースバック
			if (cij == E_hairpin(j-i-1, type, ioptseq[i+1], ioptseq[j-1], "NNNNNNNNN", P)){
				goto OUTLOOP;
			}
		}
	    // Hairpinに該当がなければ、Internal loopのトレースバック。もっとも手強い。
	    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
	    	int minq = j-i+p-MAXLOOP-2;
	    	if (minq<p+1+TURN) minq = p+1+TURN;
	    	for (q = j-1; q >= minq; q--) {

	    		int type_pq = BP_pair[ioptseq[p]][ioptseq[q]];
	    		if (type_pq==0) continue;
	    		type_pq = rtype[type_pq];

	    		int energy = E_intloop(p-i-1, j-q-1, type, type_pq,
	    				ioptseq[i+1], ioptseq[j-1], ioptseq[p-1], ioptseq[q+1], P);

	    		int energy_new = energy+c[getIndx(p,q,w,indx)];
	    		traced = (cij == energy_new);
	    		if (traced) {
	    			base_pair[++b].i = p;
	    			base_pair[b].j   = q;

	    			i = p, j = q; // i,jの更新
	    			goto repeat1;
	    		}
	    	}
	    }
	    /* end of repeat: --------------------------------------------------*/

	    /* (i.j) must close a multi-loop */

	    int type_rev = rtype[type];
	    i1 = i+1; j1 = j-1;

	    sector[s+1].ml  = sector[s+2].ml = 1;

	    int en = cij - TermAU(type_rev, P) - P->MLintern[type_rev] - P->MLclosing;
	    for(k = i+3+TURN; k < j-1-TURN; k++){
	    	//マルチループを閉じるところと、bifucationを同時に探している。
	    	if(en == m[getIndx(i+1,k-1,w,indx)] + m[getIndx(k,j-1,w,indx)]){
	    		goto LABEL2;
	    	}
	    }
	    LABEL2:

	    if (k<=j-2-TURN) { /* found the decomposition successfully*/
	    	sector[++s].i = i1;
	    	sector[s].j   = k-1;

	    	sector[++s].i = k;
	    	sector[s].j   = j1;

	    } else {
	    	fprintf(stderr, "fixed_backtracking failed in repeat %d %d\n", i , j);
//	    	exit(1);
	    }

	}

	base_pair[0].i = b;    /* save the total number of base pairs */
//	cout << base_pair[0].i << endl;

}

void fixed_fold(string optseq, int *indx, const int &w, map<string, int> &predefE,
		const int (&BP_pair)[5][5], paramT *P, char *aaseq, codon codon_table){
	int nuclen = optseq.size() - 1;
	int aalen = (optseq.size() - 1)/3;
	int size = getMatrixSize(nuclen, w);
	int C[size];
	int M[size];
	int F[nuclen+1];
	int DMl[nuclen+1];
	int DMl1[nuclen+1];
	int DMl2[nuclen+1];
	bond base_pair[nuclen/2];

	map<char, int> n2i = make_n2i();
	int ioptseq[nuclen+1];
	ioptseq[0] = 0;
	for(int i = 1; i <= nuclen; i++){
		ioptseq[i] = n2i[optseq[i]];
//		cout << optseq[i] << endl;
	}
//	exit(0);

	fixed_init_matrix(nuclen, size, C, M, F, DMl, DMl1, DMl2);
	int rtype[7] = { 0, 2, 1, 4, 3, 6, 5 };

	// investigate mfe
	const char dummy_str[10] = "XXXXXXXXX";
	for (int l = 5; l <= nuclen; l++) {
		if(l > w) break;
		//cout << "process:" << l << endl;

		//	  for(int l = 5; l <= 5; l++){
		for (int i = 1; i <= nuclen - l + 1; i++) {
			int j = i + l - 1;
			int ij = getIndx(i,j,w,indx);
			C[ij] = INF;
			M[ij] = INF;
//			cout << "test:" << j << endl;
			int type = BP_pair[ioptseq[i]][ioptseq[j]];

			if (type) {
				// hairpin
				int energy = E_hairpin(j - i - 1, type,
						ioptseq[i+1], ioptseq[j-1],
						dummy_str, P);
				C[ij] = MIN2(energy, C[ij]);

				if(l == 5 || l ==6 || l == 8){
					string hpn = "";
					for(int k = i; k <= j; k++){
						hpn += optseq[k];
					}
					//cout << i << ":" << j << "=" << hpn << endl;
					if(predefE.count(hpn) > 0){
						C[ij] = predefE[hpn];
					}
				}

				// interior loop
				for (int p = i + 1;
					p <= MIN2(j-2-TURN, i+MAXLOOP+1); p++) { // loop for position q, p
					int minq = j - i + p - MAXLOOP - 2;
					if (minq < p + 1 + TURN)
						minq = p + 1 + TURN;
					for (int q = minq; q < j; q++) {

						int pq = getIndx(p,q,w, indx);

						int type_2 =
								BP_pair[ioptseq[p]][ioptseq[q]];

						if (type_2 == 0)
							continue;
						type_2 = rtype[type_2];

						int int_energy =
								E_intloop(p	- i	- 1,
										j - q   - 1,
										type,
										type_2,
										ioptseq[i+1],
										ioptseq[j-1],
										ioptseq[p-1],
										ioptseq[q+1],
										P);

						int energy =
									int_energy
										+ C[pq];
						C[ij] =
								MIN2(energy,
										C[ij]);

					} /* end q-loop */
				} /* end p-loop */
				//cout << i << "," << j << ":" << C[ij] << endl;

				// multi-loop
				energy = DMl2[i+1];
				int tt = rtype[type];

				energy += P->MLintern[tt];
				if(tt > 2)
					energy += P->TerminalAU;

				energy += P->MLclosing;
				C[ij] =
						MIN2(energy,
								C[ij]);
			}
			else C[ij] = INF;

			// fill M
			// create M[ij] from C[ij]
			if(type){
		        int energy_M = C[ij];
		        if(type > 2)
		          energy_M += P->TerminalAU;

		        energy_M += P->MLintern[type];
		        M[ij] = energy_M;
			}

			// create M[ij] from M[i+1][j]
			int energy_M = M[getIndx(i+1, j, w, indx)]+P->MLbase;
			M[ij] = MIN2(energy_M, M[ij]);

			// create M[ij] from M[i][j-1]
			energy_M = M[getIndx(i, j-1, w, indx)]+P->MLbase;
			M[ij] = MIN2(energy_M, M[ij]);

			/* modular decomposition -------------------------------*/
			for (int k = i + 2 + TURN; k <= j - TURN - 1; k++) { // Is this correct?
				int energy_M =  M[getIndx(i,k-1,w,indx)]+M[getIndx(k,j,w,indx)];
				DMl[i] = MIN2(energy_M, DMl[i]);
				M[ij] = MIN2(energy_M, M[ij]);
			}
		}
		// rotate DMl arrays
		int FF[nuclen+1];
		for(int j = 1; j <= nuclen; j++){
			FF[j] = DMl2[j]; DMl2[j] = DMl1[j]; DMl1[j] = DMl[j]; DMl[j] = FF[j];
		}
		for(int j = 1; j <= nuclen; j++){
			DMl[j] = INF;
		}
	}

	// Fill F matrix
	// Initialize F[1]
	F[1] = 0;
	for (int j = 2; j <= nuclen; j++) {
		F[j] = INF;
		int type = BP_pair[ioptseq[1]][ioptseq[j]];
		if (type) {
			int au_penalty = 0;
			if (type > 2)
				au_penalty = P->TerminalAU;
			if(j <= w)
				F[j] = MIN2(F[j], C[getIndx(1,j,w,indx)] + au_penalty); // recc 1
		}

		// create F[j] from F[j-1]
		F[j] = MIN2(F[j], F[j - 1]); // recc 2

		for (int k = MAX2(2, j-w+1); k <= j - TURN - 1; k++) { // Is this correct?
			int type_k =
					BP_pair[ioptseq[k]][ioptseq[j]];

			int au_penalty = 0;
			if (type_k > 2)
				au_penalty = P->TerminalAU;
			int kj = getIndx(k,j,w,indx);

			int energy = F[k - 1] + C[kj] + au_penalty; // recc 4
			F[j] = MIN2(F[j], energy);
		}
	}

	int MFE = F[nuclen];
	//	showFixedMatrix(C, indx, nuclen, w);
//	showFixedMatrix(M, indx, nuclen, w);
//	for(int i = 1; i <= nuclen; i++){
//		cout << "FMAT: " << i << " " << F[i] << "  " <<  w << endl;
//	}
	fixed_backtrack(optseq, base_pair, C, M, F,
			indx, P, nuclen, w, BP_pair, predefE);
//}

//	cout << "end" << endl;
//	exit(0);
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
	optseq.erase(0, 1);
	optstr.erase(0, 1);
	cout << optseq << endl;
	cout << optstr << endl;
	cout << "MFE:" << float(MFE)/100 << " kcal/mol" << endl;

}



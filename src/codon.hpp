#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cstdlib>

using namespace std;

/* tableとextendedTableはL、Rが異なる
 table	extendedTable
 L UUA <-> UVA
 UUG <-> UVG
 CUA <-> CVA
 CUC <-> CWC
 CUG <-> CVG
 CUU <-> CWU
 R AGA <-> AXA
 AGG <-> AXG
 CGU <-> CYU
 CGC <-> CYC
 CGA <-> CXA
 CGG <-> CXG

 V、X：次の塩基はAまたはG
 W、Y：次の塩基はCまたはU
 */

class codon {
public:
	codon();

	vector<string> getCodons(char c, string exceptedCodons) {
		if (table.count(c) == 0) {
			cerr << "ERR: table doesn't have " << c << "." << endl;
			exit(1);
		}

		char delim = ',';
		vector<string> exceptVector = split(exceptedCodons, delim);
		map<string,int> exceptMap;
		for(unsigned int i=0; i<exceptVector.size();++i){
			string tmp = exceptVector.at(i);
		if (exceptMap.find(tmp) == exceptMap.end()) {
			exceptMap[tmp] = 1;
			} else {
				int count = exceptMap[tmp] + 1;
				exceptMap[tmp] = count;
			}
		}

		vector<string> filterTable;
		for (unsigned int i = 0; i < table[c].size(); ++i) {
			string codon = table[c].at(i);
			if(exceptMap.find(codon) == exceptMap.end()){
				filterTable.push_back(codon);
			}
		}

		return filterTable;
	}

	vector<string> getExtendedCodons(char c, string exceptedCodons) {
		if (extendedTable.count(c) == 0) {
			cerr << "ERR: extended table doesn't have " << c << "." << endl;
			exit(1);
		}
		char delim = ',';
		vector<string> exceptVector = split(exceptedCodons, delim);
		map<string,int> exceptMap;
		for(unsigned int i=0; i<exceptVector.size();++i){
			string exceptCodon = exceptVector.at(i);
			string convertedExceptCodon = exceptCodon;
			map<string,string>::iterator itr;
			if((itr = expectedCodonOfCodon.find(exceptCodon)) != expectedCodonOfCodon.end()){
				convertedExceptCodon = itr->second;
			}

		if (exceptMap.find(convertedExceptCodon) == exceptMap.end()) {
			exceptMap[convertedExceptCodon] = 1;
			} else {
				int count = exceptMap[convertedExceptCodon] + 1;
				exceptMap[convertedExceptCodon] = count;
			}
		}

		// 除外コドンを除いたコドンテーブルを作成
		vector<string> filterTable;
		for (unsigned int i = 0; i < extendedTable[c].size(); ++i) {
			string codon = extendedTable[c].at(i);
			if(exceptMap.find(codon) == exceptMap.end()){
				filterTable.push_back(codon);
			}
		}
		return filterTable;
	}

	char c2a(int p1, int p2, int p3) {
		return table_rev[p1][p2][p3];
	}

	void showTable() {
		for (map<char, vector<string> >::iterator it = table.begin();
				it != table.end(); ++it) {
			char aa = it->first;
			vector<string> codons = getCodons(aa,"");
			for (unsigned int j = 0; j < codons.size(); j++) {
				cout << aa << " " << codons[j] << endl;
			}
		}
	}

private:
	map<char, vector<string> > table;
	map<char, vector<string> > extendedTable;
	map<string, string > expectedCodonOfCodon;
	char table_rev[5][5][5];

	vector<string> split(string &str, char delim) {
		istringstream iss(str);
		string tmp;
		vector<string> res;

		while (getline(iss, tmp, delim)) {
			res.push_back(tmp);
		}
		return res;
	}
};

codon::codon() {
	expectedCodonOfCodon.insert(make_pair("UUA", "UVA"));
	expectedCodonOfCodon.insert(make_pair("UUG", "UVG"));
	expectedCodonOfCodon.insert(make_pair("CUU", "CWU"));
	expectedCodonOfCodon.insert(make_pair("CUC", "CWC"));
	expectedCodonOfCodon.insert(make_pair("CUA", "CVA"));
	expectedCodonOfCodon.insert(make_pair("CUG", "CVG"));
	expectedCodonOfCodon.insert(make_pair("CGU", "CYU"));
	expectedCodonOfCodon.insert(make_pair("CGC", "CYC"));
	expectedCodonOfCodon.insert(make_pair("CGA", "CXA"));
	expectedCodonOfCodon.insert(make_pair("CGG", "CXG"));
	expectedCodonOfCodon.insert(make_pair("AGA", "AXA"));
	expectedCodonOfCodon.insert(make_pair("AGG", "AXG"));

	table['F'].push_back("UUU");
	table['F'].push_back("UUC");
	table['L'].push_back("UUA");
	table['L'].push_back("UUG");
	table['L'].push_back("CUU");
	table['L'].push_back("CUC");
	table['L'].push_back("CUA");
	table['L'].push_back("CUG");
	table['I'].push_back("AUU");
	table['I'].push_back("AUC");
	table['I'].push_back("AUA");
	table['M'].push_back("AUG");
	table['V'].push_back("GUU");
	table['V'].push_back("GUC");
	table['V'].push_back("GUA");
	table['V'].push_back("GUG");
	table['S'].push_back("UCU");
	table['S'].push_back("UCC");
	table['S'].push_back("UCA");
	table['S'].push_back("UCG");
	table['S'].push_back("AGU");
	table['S'].push_back("AGC");
	table['P'].push_back("CCU");
	table['P'].push_back("CCC");
	table['P'].push_back("CCA");
	table['P'].push_back("CCG");
	table['T'].push_back("ACU");
	table['T'].push_back("ACC");
	table['T'].push_back("ACA");
	table['T'].push_back("ACG");
	table['A'].push_back("GCU");
	table['A'].push_back("GCC");
	table['A'].push_back("GCA");
	table['A'].push_back("GCG");
	table['Y'].push_back("UAU");
	table['Y'].push_back("UAC");
	table['*'].push_back("UAA");
	table['*'].push_back("UAG");
	table['*'].push_back("UGA");
	table['H'].push_back("CAU");
	table['H'].push_back("CAC");
	table['Q'].push_back("CAA");
	table['Q'].push_back("CAG");
	table['N'].push_back("AAU");
	table['N'].push_back("AAC");
	table['K'].push_back("AAA");
	table['K'].push_back("AAG");
	table['D'].push_back("GAU");
	table['D'].push_back("GAC");
	table['E'].push_back("GAA");
	table['E'].push_back("GAG");
	table['C'].push_back("UGU");
	table['C'].push_back("UGC");
	table['W'].push_back("UGG");
	table['R'].push_back("CGU");
	table['R'].push_back("CGC");
	table['R'].push_back("CGA");
	table['R'].push_back("CGG");
	table['R'].push_back("AGA");
	table['R'].push_back("AGG");
	table['G'].push_back("GGU");
	table['G'].push_back("GGC");
	table['G'].push_back("GGA");
	table['G'].push_back("GGG");

	extendedTable['F'].push_back("UUU");
	extendedTable['F'].push_back("UUC");
	extendedTable['L'].push_back("UVA");
	extendedTable['L'].push_back("UVG");
	extendedTable['L'].push_back("CWU");
	extendedTable['L'].push_back("CWC");
	extendedTable['L'].push_back("CVA");
	extendedTable['L'].push_back("CVG");
	extendedTable['I'].push_back("AUU");
	extendedTable['I'].push_back("AUC");
	extendedTable['I'].push_back("AUA");
	extendedTable['M'].push_back("AUG");
	extendedTable['V'].push_back("GUU");
	extendedTable['V'].push_back("GUC");
	extendedTable['V'].push_back("GUA");
	extendedTable['V'].push_back("GUG");
	extendedTable['S'].push_back("UCU");
	extendedTable['S'].push_back("UCC");
	extendedTable['S'].push_back("UCA");
	extendedTable['S'].push_back("UCG");
	extendedTable['S'].push_back("AGU");
	extendedTable['S'].push_back("AGC");
	extendedTable['P'].push_back("CCU");
	extendedTable['P'].push_back("CCC");
	extendedTable['P'].push_back("CCA");
	extendedTable['P'].push_back("CCG");
	extendedTable['T'].push_back("ACU");
	extendedTable['T'].push_back("ACC");
	extendedTable['T'].push_back("ACA");
	extendedTable['T'].push_back("ACG");
	extendedTable['A'].push_back("GCU");
	extendedTable['A'].push_back("GCC");
	extendedTable['A'].push_back("GCA");
	extendedTable['A'].push_back("GCG");
	extendedTable['Y'].push_back("UAU");
	extendedTable['Y'].push_back("UAC");
	extendedTable['*'].push_back("UAA");
	extendedTable['*'].push_back("UAG");
	extendedTable['*'].push_back("UGA");
	extendedTable['H'].push_back("CAU");
	extendedTable['H'].push_back("CAC");
	extendedTable['Q'].push_back("CAA");
	extendedTable['Q'].push_back("CAG");
	extendedTable['N'].push_back("AAU");
	extendedTable['N'].push_back("AAC");
	extendedTable['K'].push_back("AAA");
	extendedTable['K'].push_back("AAG");
	extendedTable['D'].push_back("GAU");
	extendedTable['D'].push_back("GAC");
	extendedTable['E'].push_back("GAA");
	extendedTable['E'].push_back("GAG");
	extendedTable['C'].push_back("UGU");
	extendedTable['C'].push_back("UGC");
	extendedTable['W'].push_back("UGG");
	extendedTable['R'].push_back("CYU");
	extendedTable['R'].push_back("CYC");
	extendedTable['R'].push_back("CXA");
	extendedTable['R'].push_back("CXG");
	extendedTable['R'].push_back("AXA");
	extendedTable['R'].push_back("AXG");
	extendedTable['G'].push_back("GGU");
	extendedTable['G'].push_back("GGC");
	extendedTable['G'].push_back("GGA");
	extendedTable['G'].push_back("GGG");

	table_rev[4][4][4] = 'F';
	table_rev[4][4][2] = 'F';
	table_rev[4][4][1] = 'L';
	table_rev[4][4][3] = 'L';

	table_rev[2][4][4] = 'L';
	table_rev[2][4][2] = 'L';
	table_rev[2][4][1] = 'L';
	table_rev[2][4][3] = 'L';

	table_rev[1][4][4] = 'I';
	table_rev[1][4][2] = 'I';
	table_rev[1][4][1] = 'I';
	table_rev[1][4][3] = 'M';

	table_rev[3][4][4] = 'V';
	table_rev[3][4][2] = 'V';
	table_rev[3][4][1] = 'V';
	table_rev[3][4][3] = 'V';

	table_rev[4][2][4] = 'S';
	table_rev[4][2][2] = 'S';
	table_rev[4][2][1] = 'S';
	table_rev[4][2][3] = 'S';

	table_rev[2][2][4] = 'P';
	table_rev[2][2][2] = 'P';
	table_rev[2][2][1] = 'P';
	table_rev[2][2][3] = 'P';

	table_rev[1][2][4] = 'T';
	table_rev[1][2][2] = 'T';
	table_rev[1][2][1] = 'T';
	table_rev[1][2][3] = 'T';

	table_rev[3][2][4] = 'A';
	table_rev[3][2][2] = 'A';
	table_rev[3][2][1] = 'A';
	table_rev[3][2][3] = 'A';

	table_rev[4][1][4] = 'Y';
	table_rev[4][1][2] = 'Y';
	table_rev[4][1][1] = '*';
	table_rev[4][1][3] = '*';

	table_rev[2][1][4] = 'H';
	table_rev[2][1][2] = 'H';
	table_rev[2][1][1] = 'Q';
	table_rev[2][1][3] = 'Q';

	table_rev[1][1][4] = 'N';
	table_rev[1][1][2] = 'N';
	table_rev[1][1][1] = 'K';
	table_rev[1][1][3] = 'K';

	table_rev[3][1][4] = 'D';
	table_rev[3][1][2] = 'D';
	table_rev[3][1][1] = 'E';
	table_rev[3][1][3] = 'E';

	table_rev[4][3][4] = 'C';
	table_rev[4][3][2] = 'C';
	table_rev[4][3][1] = '*';
	table_rev[4][3][3] = 'W';

	table_rev[2][3][4] = 'R';
	table_rev[2][3][2] = 'R';
	table_rev[2][3][1] = 'R';
	table_rev[2][3][3] = 'R';

	table_rev[1][3][4] = 'S';
	table_rev[1][3][2] = 'S';
	table_rev[1][3][1] = 'R';
	table_rev[1][3][3] = 'R';

	table_rev[3][3][4] = 'G';
	table_rev[3][3][2] = 'G';
	table_rev[3][3][1] = 'G';
	table_rev[3][3][3] = 'G';

}


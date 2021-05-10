/*
 * AASeqConverter.hpp
 */

#ifndef AASEQCONVERTER_H_
#define AASEQCONVERTER_H_

#include "FileManager.hpp"
#include <string>
#include <vector>
#include <algorithm>
//#include "codon.hpp"
#include <time.h>
#include "Util.hpp"
#include <limits>

using namespace std;

class AASeqConverter {
public:
	AASeqConverter() {
		getBaseNumberMap();
		getPairNumberMap();
		getBaseEnergy();
	}

	~AASeqConverter() {
	}

	map<string, int> getBaseEnergy() {
		if (baseEnergy.empty()) {
			FileManager fm;
			fm.loadEnergyFile(baseEnergy);
		}
		return baseEnergy;
	}

	vector<vector<int> > countNeighborTwoBase(string aaseq,
			string exceptedCodons) {
		vector<vector<int> > result;

		// 結果格納マップを初期化
		int twoPairSize = aaseq.size() * 3 - 1;
		initBasePairMap(twoPairSize, result);

		// アミノ酸配列の文字列を配列に変換
		int position = 0;
		map<string, int> preBaseMap;
		for (unsigned int i = 0; i < aaseq.size(); ++i) {
			char aa = aaseq[i];

			// アミノ酸の塩基候補リストを取得
			vector<string> baseList = codonTable.getExtendedCodons(aa,
					exceptedCodons);

			// 塩基配列の候補配列を作成
			map<string, int> nextPreBaseMap;
			vector<string>::iterator itr;
			itr = baseList.begin();
			while (itr != baseList.end()) {
				string codon = *itr;

				if (preBaseMap.empty()) {

				} else {
					map<string, int>::iterator mapItr;
					for (mapItr = preBaseMap.begin();
							mapItr != preBaseMap.end(); mapItr++) {
						string preBase = mapItr->first;
						string basePair = preBase + codon.substr(0, 1);
						setBasePairMap(position, basePair, result);
					}
				}

				setBasePairMap(position + 1, codon.substr(0, 2), result);
				setBasePairMap(position + 2, codon.substr(1, 2), result);
				nextPreBaseMap.insert(make_pair(codon.substr(2, 1), 1));
				itr++;
			}
			position = position + 3;
			preBaseMap = nextPreBaseMap;

		}

		return result;
	}

	vector<vector<int> > countEveryOtherTwoBase(string aaseq,
			string exceptedCodons) {
		vector<vector<int> > result;

		// 結果格納マップを初期化
		int twoPairSize = aaseq.size() * 3 - 2;
		initBasePairMap(twoPairSize, result);

		// アミノ酸配列の文字列を配列に変換
		int position = 0;
		map<string, int> preOneBaseMap;
		map<string, int> preTwoBaseMap;
		for (unsigned int i = 0; i < aaseq.size(); ++i) {
			char aa = aaseq[i];

			// アミノ酸の塩基候補リストを取得
			vector<string> baseList = codonTable.getExtendedCodons(aa,
					exceptedCodons);

			// 塩基配列の候補配列を作成
			map<string, int> nextOneBaseMap;
			map<string, int> nextTwoBaseMap;
			vector<string>::iterator itr;
			itr = baseList.begin();
			while (itr != baseList.end()) {
				string codon = *itr;

				if (preOneBaseMap.empty()) {

				} else {
					map<string, int>::iterator mapItr;
					for (mapItr = preTwoBaseMap.begin();
							mapItr != preTwoBaseMap.end(); mapItr++) {
						string preBase = mapItr->first;
						string basePair = preBase + codon.substr(0, 1);
						setBasePairMap(position - 1, basePair, result);
					}
					for (mapItr = preOneBaseMap.begin();
							mapItr != preOneBaseMap.end(); mapItr++) {
						string preBase = mapItr->first;
						string basePair = preBase + codon.substr(1, 1);
						setBasePairMap(position, basePair, result);
					}
				}
				string twoBase = codon.substr(0, 1) + codon.substr(2, 1);
				setBasePairMap(position + 1, twoBase, result);
				nextTwoBaseMap.insert(make_pair(codon.substr(1, 1), 1));
				nextOneBaseMap.insert(make_pair(codon.substr(2, 1), 1));
				itr++;
			}
			position = position + 3;
			preTwoBaseMap = nextTwoBaseMap;
			preOneBaseMap = nextOneBaseMap;

		}

		return result;
	}

	vector<vector<vector<string> > >  getExtendedBases(
		string aminoAcid, string exceptedCodons) {
		// 0～(n-8)要素の部分塩基配列を取得
		unsigned int maxLength = 8;
		vector<vector<vector<string> > > bases = getSelectedLengthExtendedBases(aminoAcid, maxLength,
				exceptedCodons);

		// (n-8+1)以降の部分塩基配列を取得
		int preMaxLength = maxLength;
		for(int i=maxLength-1;i > 0;--i){
			maxLength = i;
			int startElement = aminoAcid.size() * 3 - preMaxLength + 1;
			addExtendedBases(aminoAcid, startElement, maxLength, bases,
					exceptedCodons);
			preMaxLength = maxLength;
		}

		return bases;
	}

	/*
	 * アミノ酸配列の部分塩基配列部位の自由エネルギーを取得する
	 *
	 * @param　aminoAcid　アミノ酸配列
	 * @return　各塩基配列における自由エネルギーのベクター
	 * 　要素1次元：開始位置（0～アミノ酸配列の塩基数、実際には要素番号1以降を使用）
	 * 　要素2次元：部分配列長（0～8、実際には要素番号は、5,6,8のみ使用）
	 * 　要素3次元：開始位置塩基番号（0～8、実際には要素番号1以降を使用）
	 * 　要素4次元：終了位置塩基番号（0～8、実際には要素番号1以降を使用）
	 * 　　（塩基番号：A=1, C=2, G=3, U=4, V=5, W=6, X=7, Y=8)
	 *
	 * （使用例）
	 *   string seq = "MLYF";
	 * 	 AASeqConverter conv;
	 *	 vector<vector<vector<vector<pair<int, string> > > > > result
	 *	   = conv.calcQueryBaseEnergy(seq);
	 *	 vector<int> baseLengths;
	 *	 baseLengths.push_back(5);
	 *	 baseLengths.push_back(6);
	 *	 baseLengths.push_back(8);
	 *	 for (unsigned int start = 1; start <= result.size(); start++) {
	 *	 	for (unsigned int size = 0; size < baseLengths.size(); size++) {
	 *	 		int baseSize = baseLengths[size];
	 *	 		for (int sBase = 1; sBase <= 8; sBase++) {
	 *	 			for (int eBase = 1; eBase <= 8; eBase++) {
	 *	 				int energy = result[start][baseSize][sBase][eBase].first;
	 *	 				string seq = result[start][baseSize][sBase][eBase].second;
	 *	 				cout << start << "\t" << baseSize << "\t" << sBase << "\t" << eBase
	 *	 						<< "\t" << energy << "\t" << seq << endl;
	 *	 			}
	 *	 		}
	 *	 	}
	 *	 }
	 */
	vector<vector<vector<vector<pair<int, string> > > > > calcQueryExtendedBaseEnergy(
			string aminoAcid, string exceptedCodons) {
		// 部分配列を取得
		vector<vector<vector<string> > > bases = getExtendedBases(aminoAcid,exceptedCodons);

		// 各部位の最小エネルギーを取得
		vector<vector<vector<vector<pair<int, string> > > > > result =
				calcEachExtendedBaseEnergy(aminoAcid, bases);

		return result;
	}

	vector<vector<vector<string> > > getOriginalBases(
			string aminoAcid, string exceptedCodons) {
		// 0～(n-8)要素の部分塩基配列を取得
		unsigned int maxLength = 8;
		vector<vector<vector<string> > > bases = getSelectedLengthOriginalBases(aminoAcid,
				maxLength,exceptedCodons);

		// (n-8+1)以降の部分塩基配列を取得
		int preMaxLength = maxLength;
		for(int i=maxLength-1;i > 0;--i){
			maxLength = i;
			int startElement = aminoAcid.size() * 3 - preMaxLength + 1;
			addExtendedBases(aminoAcid, startElement, maxLength, bases,
					exceptedCodons);
			preMaxLength = maxLength;
		}

//		// (n-8+1)～(n-6)要素の部分塩基配列を取得
//		int preMaxLength = 8;
//		int maxLength = 6;
//		int startElement = aminoAcid.size() * 3 - preMaxLength + 1;
//		addOriginalBases(aminoAcid, startElement, maxLength, bases,
//				exceptedCodons);
//
//		// (n-6+1)～(n-5)要素の部分塩基配列を取得
//		preMaxLength = maxLength;
//		maxLength = 5;
//		startElement = aminoAcid.size() * 3 - preMaxLength + 1;
//		addOriginalBases(aminoAcid, startElement, maxLength, bases,
//				exceptedCodons);

		return bases;
	}

	/*
	 * アミノ酸配列の部分塩基配列部位の自由エネルギーを取得する
	 *
	 * @param　aminoAcid　アミノ酸配列
	 * @return　各塩基配列における自由エネルギーのベクター
	 * 　要素1次元：開始位置（0～アミノ酸配列の塩基数、実際には要素番号1以降を使用）
	 * 　要素2次元：部分配列長（0～8、実際には要素番号は、5,6,8のみ使用）
	 * 　要素3次元：開始位置塩基番号（0～4、実際には要素番号1以降を使用）
	 * 　要素4次元：終了位置塩基番号（0～4、実際には要素番号1以降を使用）
	 * 　　（塩基番号：A=1, C=2, G=3, U=4)
	 *
	 * （使用例）
	 *   string seq = "MLYF";
	 * 	 AASeqConverter conv;
	 *	 vector<vector<vector<vector<pair<int, string> > > > > result
	 *	   = conv.calcQueryBaseEnergy(seq);
	 *	 vector<int> baseLengths;
	 *	 baseLengths.push_back(5);
	 *	 baseLengths.push_back(6);
	 *	 baseLengths.push_back(8);
	 *	 for (unsigned int start = 1; start <= result.size(); start++) {
	 *	 	for (unsigned int size = 0; size < baseLengths.size(); size++) {
	 *	 		int baseSize = baseLengths[size];
	 *	 		for (int sBase = 1; sBase <= 4; sBase++) {
	 *	 			for (int eBase = 1; eBase <= 4; eBase++) {
	 *	 				int energy = result[start][baseSize][sBase][eBase].first;
	 *	 				string seq = result[start][baseSize][sBase][eBase].second;
	 *	 				cout << start << "\t" << baseSize << "\t" << sBase << "\t" << eBase
	 *	 						<< "\t" << energy << "\t" << seq << endl;
	 *	 			}
	 *	 		}
	 *	 	}
	 *	 }
	 */
	vector<vector<vector<vector<pair<int, string> > > > > calcQueryOriginalBaseEnergy(
			string aminoAcid, string exceptedCodons) {
		// 部分塩基配列を取得
		vector<vector<vector<string> > > bases = getOriginalBases(aminoAcid,
				exceptedCodons);

		// 各部位の最小エネルギーを取得
		vector<vector<vector<vector<pair<int, string> > > > > result =
				calcEachOriginalBaseEnergy(aminoAcid, bases);

		return result;
	}

private:
	codon codonTable;
	map<string, int> baseNumberMap;
	map<string, int> pairNumberMap;
	map<string, int> baseEnergy;

	const static int BASE_ORIGINAL = 0;
	const static int BASE_EXTENDED = 1;

	int getBaseNumber(string base) {
		map<string, int>::iterator itr;
		itr = baseNumberMap.find(base);
		return itr->second;
	}

	int getPairNumber(string twoBases) {
		map<string, int>::iterator itr;
		itr = pairNumberMap.find(twoBases);
		return itr->second;
	}

	void getPairNumberMap() {
		pairNumberMap.insert(make_pair("AA", 1));
		pairNumberMap.insert(make_pair("AC", 2));
		pairNumberMap.insert(make_pair("AG", 3));
		pairNumberMap.insert(make_pair("AU", 4));
		pairNumberMap.insert(make_pair("AV", 5));
		pairNumberMap.insert(make_pair("AW", 6));
		pairNumberMap.insert(make_pair("AX", 7));
		pairNumberMap.insert(make_pair("AY", 8));
		pairNumberMap.insert(make_pair("CA", 9));
		pairNumberMap.insert(make_pair("CC", 10));
		pairNumberMap.insert(make_pair("CG", 11));
		pairNumberMap.insert(make_pair("CU", 12));
		pairNumberMap.insert(make_pair("CV", 13));
		pairNumberMap.insert(make_pair("CW", 14));
		pairNumberMap.insert(make_pair("CX", 15));
		pairNumberMap.insert(make_pair("CY", 16));
		pairNumberMap.insert(make_pair("GA", 17));
		pairNumberMap.insert(make_pair("GC", 18));
		pairNumberMap.insert(make_pair("GG", 19));
		pairNumberMap.insert(make_pair("GU", 20));
		pairNumberMap.insert(make_pair("GV", 21));
		pairNumberMap.insert(make_pair("GW", 22));
		pairNumberMap.insert(make_pair("GX", 23));
		pairNumberMap.insert(make_pair("GY", 24));
		pairNumberMap.insert(make_pair("UA", 25));
		pairNumberMap.insert(make_pair("UC", 26));
		pairNumberMap.insert(make_pair("UG", 27));
		pairNumberMap.insert(make_pair("UU", 28));
		pairNumberMap.insert(make_pair("UV", 29));
		pairNumberMap.insert(make_pair("UW", 30));
		pairNumberMap.insert(make_pair("UX", 31));
		pairNumberMap.insert(make_pair("UY", 32));
		pairNumberMap.insert(make_pair("VA", 33));
		pairNumberMap.insert(make_pair("VC", 34));
		pairNumberMap.insert(make_pair("VG", 35));
		pairNumberMap.insert(make_pair("VU", 36));
		pairNumberMap.insert(make_pair("VV", 37));
		pairNumberMap.insert(make_pair("VW", 38));
		pairNumberMap.insert(make_pair("VX", 39));
		pairNumberMap.insert(make_pair("VY", 40));
		pairNumberMap.insert(make_pair("WA", 41));
		pairNumberMap.insert(make_pair("WC", 42));
		pairNumberMap.insert(make_pair("WG", 43));
		pairNumberMap.insert(make_pair("WU", 44));
		pairNumberMap.insert(make_pair("WV", 45));
		pairNumberMap.insert(make_pair("WW", 46));
		pairNumberMap.insert(make_pair("WX", 47));
		pairNumberMap.insert(make_pair("WY", 48));
		pairNumberMap.insert(make_pair("XA", 49));
		pairNumberMap.insert(make_pair("XC", 50));
		pairNumberMap.insert(make_pair("XG", 51));
		pairNumberMap.insert(make_pair("XU", 52));
		pairNumberMap.insert(make_pair("XV", 53));
		pairNumberMap.insert(make_pair("XW", 54));
		pairNumberMap.insert(make_pair("XX", 55));
		pairNumberMap.insert(make_pair("XY", 56));
		pairNumberMap.insert(make_pair("YA", 57));
		pairNumberMap.insert(make_pair("YC", 58));
		pairNumberMap.insert(make_pair("YG", 59));
		pairNumberMap.insert(make_pair("YU", 60));
		pairNumberMap.insert(make_pair("YV", 61));
		pairNumberMap.insert(make_pair("YW", 62));
		pairNumberMap.insert(make_pair("YX", 63));
		pairNumberMap.insert(make_pair("YY", 64));
	}

	void getBaseNumberMap() {
		baseNumberMap.insert(make_pair("A", 1));
		baseNumberMap.insert(make_pair("C", 2));
		baseNumberMap.insert(make_pair("G", 3));
		baseNumberMap.insert(make_pair("U", 4));
		baseNumberMap.insert(make_pair("V", 5));
		baseNumberMap.insert(make_pair("W", 6));
		baseNumberMap.insert(make_pair("X", 7));
		baseNumberMap.insert(make_pair("Y", 8));
	}

	void initBasePairMap(int size, vector<vector<int> > &map) {
		// 各行に代入する列の配列を作成（実際に使用するのは、要素1以降）
		vector<int> rows;
		for (int i = 0; i <= size; i++) {
			rows.push_back(0);
		}

		// 各行に列を代入（実際に使用するのは、要素1以降）
		for (unsigned int i = 0; i <= pairNumberMap.size(); i++) {
			map.push_back(rows);
		}
	}

	void setBasePairMap(int position, string basePair,
			vector<vector<int> > &map) {
		int pair = getPairNumber(basePair);
		map[pair][position] = 1;
	}

	void initResultBaseVector(int inputBaseLength, int getBaseLength,
			vector<vector<vector<string> > > &result) {
		vector<vector<string> > result_length;
		vector<string> result_base;

		for (int i = 0; i <= getBaseLength; i++) {
			result_length.push_back(result_base);
		}

		for (int i = 0; i <= inputBaseLength; i++) {
			result.push_back(result_length);
		}
	}

	void initResultExtendedBaseEnergyVector(int inputBaseLength,
			int getBaseLength,
			vector<vector<vector<vector<pair<int, string> > > > > &result) {
		initResultBaseEnergyVectorBase(inputBaseLength, getBaseLength, result,
				BASE_EXTENDED);
	}

	void initResultOriginalBaseEnergyVector(int inputBaseLength,
			int getBaseLength,
			vector<vector<vector<vector<pair<int, string> > > > > &result) {
		initResultBaseEnergyVectorBase(inputBaseLength, getBaseLength, result,
				BASE_ORIGINAL);
	}

	void initResultBaseEnergyVectorBase(int inputBaseLength, int getBaseLength,
			vector<vector<vector<vector<pair<int, string> > > > > &result,
			int flag) {
		pair<int, string> baseEnergy;
		baseEnergy.first = numeric_limits<int>::max();
		baseEnergy.second = "";
		// 塩基数
		int baseNumber;
		if (flag == BASE_ORIGINAL) {
			baseNumber = 4;
		} else if (flag == BASE_EXTENDED) {
			baseNumber = 8;
		}
		vector<pair<int, string> > resEndBase;
		vector<vector<pair<int, string> > > resStartBase;
		vector<vector<vector<pair<int, string> > > > resLength;

		// resEndBaseに適切な数のbaseEnergyを追加（実際に使用するのは、要素1以降）
		resEndBase.assign(baseNumber + 1, baseEnergy);

		// resStartBaseに適切な数のresEndBaseを追加（実際に使用するのは、要素1以降）
		resStartBase.assign(baseNumber + 1, resEndBase);

		// resLengthに適切な数のresStartBaseを追加（実際に使用するのは、要素1以降）
		resLength.assign(getBaseLength + 1, resStartBase);

		// resultに適切な数のresLengthを追加（実際に使用するのは、要素1以降）
		result.assign(inputBaseLength + 1, resLength);
	}

	void getBasesBase(string aminoAcid, unsigned int maxLength,
			vector<vector<vector<string> > > &result, int flag,
			string exceptedCodons) {
		int baseLength = aminoAcid.size() * 3;

		// 結果格納ベクターを初期化
		initResultBaseVector(baseLength, maxLength, result);

		// maxLengthの長さを持つ塩基配列を作成
		for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
			vector<string> bases; // 作成した塩基配列を格納
			unsigned int startPosition = i + 1;
			for (unsigned int j = 0; j < aminoAcid.size(); j++) {
				if (startPosition > (j + 1) * 3) {
					// 開始位置より手前なら処理をスキップ
					continue;
				} else {
					// 作成中の塩基配列に連結する塩基配列をコドンから取得
					map<string, int> baseMap;
					int nowPosition = (j + 1) * 3 - 2;
					int addLength;
					if (bases.empty()) {
						addLength = maxLength;
					} else {
						addLength = maxLength - bases[0].length();
					}

					vector<string> codons;
					if (flag == BASE_ORIGINAL) {
						codons = codonTable.getCodons(aminoAcid.at(j),
								exceptedCodons);
					} else if (flag == BASE_EXTENDED) {
						codons = codonTable.getExtendedCodons(aminoAcid.at(j),
								exceptedCodons);
					} else {
						cout << "エラー：getBasesBase()のflagが無効な値です:" << flag
								<< endl;
						exit(1);
					}

					for (unsigned int k = 0; k < codons.size(); k++) {
						string codon = codons[k];
						int addNumber = 0;
						string addingBase;
						for (unsigned int l = 0; l < codon.size(); l++) {
							if (nowPosition + l >= startPosition) {
								// maxLengthを超える場合は処理をスキップ
								if (addLength - addNumber > 0) {
									addingBase += codon.at(l);
									addNumber++;
								}
							}
						}

						baseMap.insert(make_pair(addingBase, 1));
					}

					// コドンから取得した塩基配列を作成中の塩基配列に連結
					vector<string> newBases;
					map<string, int>::iterator itr;
					for (itr = baseMap.begin(); itr != baseMap.end(); itr++) {
						string base = itr->first;

						if (bases.empty()) {
							newBases.push_back(base);
						} else {
							for (unsigned int l2 = 0; l2 < bases.size(); l2++) {
								string newBase = bases[l2];
								newBase += base;
								newBases.push_back(newBase);
							}
						}
					}

					// 塩基配列を伸長した新しい結果に置き換える
					if (!newBases.empty()) {
						bases = newBases;
					}

					// 部分配列の長さmaxLengthに達したらループを抜ける
					if (bases[0].length() == maxLength) {
						break;
					}
				}
			}

			// 作成した塩基配列を結果格納ベクターにセット
			for (unsigned int m = 0; m < bases.size(); m++) {
				result[startPosition][maxLength].push_back(bases[m]);
			}
		}

		// maxLength以下の部分配列を作成
		for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
			int startPosition = i + 1;
			for (int j = maxLength; j > 1; j--) {
				map<string, int> baseMap;
				for (unsigned int k = 0; k < result[startPosition][j].size(); k++) {
					string base = result[startPosition][j][k];
					baseMap.insert(make_pair(base.substr(0, j - 1), 1));
				}

				map<string, int>::iterator itr;
				for (itr = baseMap.begin(); itr != baseMap.end(); itr++) {
					string base = itr->first;
					result[startPosition][j - 1].push_back(base);
				}
			}
		}
	}

	void addBasesBase(string aminoAcid, int startElement,
			unsigned int maxLength, vector<vector<vector<string> > > &result,
			int flag, string exceptedCodons) {
		int baseLength = aminoAcid.size() * 3;

		// maxLengthの長さを持つ塩基配列を作成
		vector<vector<vector<string> > > nowResult;
		initResultBaseVector(baseLength, maxLength, nowResult);
		for (unsigned int i = startElement; i <= baseLength - maxLength; i++) {
			vector<string> bases; // 作成した塩基配列を格納
			unsigned int startPosition = i + 1;
			for (unsigned int j = 0; j < aminoAcid.size(); j++) {
				if (startPosition > (j + 1) * 3) {
					// 開始位置より手前なら処理をスキップ
					continue;
				} else {
					// 作成中の塩基配列に連結する塩基配列をコドンから取得
					map<string, int> baseMap;
					int nowPosition = (j + 1) * 3 - 2;
					int addLength;
					if (bases.empty()) {
						addLength = maxLength;
					} else {
						addLength = maxLength - bases[0].length();
					}

					vector<string> codons;
					if (flag == BASE_ORIGINAL) {
						codons = codonTable.getCodons(aminoAcid.at(j),
								exceptedCodons);
					} else if (flag == BASE_EXTENDED) {
						codons = codonTable.getExtendedCodons(aminoAcid.at(j),
								exceptedCodons);
					} else {
						cout << "エラー：addBasesBase()のflagが無効な値です:" << flag
								<< endl;
						exit(1);
					}

					for (unsigned int k = 0; k < codons.size(); k++) {
						string codon = codons[k];
						int addNumber = 0;
						string addingBase;
						for (unsigned int l = 0; l < codon.size(); l++) {
							if (nowPosition + l >= startPosition) {
								// maxLengthを超える場合は処理をスキップ
								if (addLength - addNumber > 0) {
									addingBase += codon.at(l);
									addNumber++;
								}
							}
						}

						baseMap.insert(make_pair(addingBase, 1));
					}

					// コドンから取得した塩基配列を作成中の塩基配列に連結
					vector<string> newBases;
					map<string, int>::iterator itr;
					for (itr = baseMap.begin(); itr != baseMap.end(); itr++) {
						string base = itr->first;

						if (bases.empty()) {
							newBases.push_back(base);
						} else {
							for (unsigned int l2 = 0; l2 < bases.size(); l2++) {
								string newBase = bases[l2];
								newBase += base;
								newBases.push_back(newBase);
							}
						}
					}

					// 塩基配列を伸長した新しい結果に置き換える
					if (!newBases.empty()) {
						bases = newBases;
					}

					// 部分配列の長さmaxLengthに達したらループを抜ける
					if (bases[0].length() == maxLength) {
						break;
					}

				}
			}

			// 作成した塩基配列を結果格納ベクターにセット
			for (unsigned int m = 0; m < bases.size(); m++) {
				nowResult[startPosition][maxLength].push_back(bases[m]);
			}
		}

		// maxLength以下の部分配列を作成
		for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
			int startPosition = i + 1;
			for (int j = maxLength; j > 1; j--) {
				map<string, int> baseMap;
				for (unsigned int k = 0; k < nowResult[startPosition][j].size(); k++) {
					string base = nowResult[startPosition][j][k];
					baseMap.insert(make_pair(base.substr(0, j - 1), 1));
				}

				map<string, int>::iterator itr;
				for (itr = baseMap.begin(); itr != baseMap.end(); itr++) {
					string base = itr->first;
					nowResult[startPosition][j - 1].push_back(base);
				}
			}
		}

		// 取得した内容を結果ベクターに設定
		for (unsigned int i = 1; i < nowResult.size(); i++) {
			for (unsigned int j = 1; j < nowResult[i].size(); j++) {
				for (unsigned int k = 0; k < nowResult[i][j].size(); k++) {
					string base = nowResult[i][j][k];
					result[i][j].push_back(base);
				}
			}
		}
	}

	void addExtendedBases(string aminoAcid, int startElement,
			unsigned int maxLength, vector<vector<vector<string> > > &result,
			string exceptedCodons) {
		addBasesBase(aminoAcid, startElement, maxLength, result, BASE_EXTENDED,
				exceptedCodons);
	}

	void addOriginalBases(string aminoAcid, int startElement,
			unsigned int maxLength, vector<vector<vector<string> > > &result,
			string exceptedCodons) {
		addBasesBase(aminoAcid, startElement, maxLength, result, BASE_ORIGINAL,
				exceptedCodons);
	}

	vector<vector<vector<vector<pair<int, string> > > > > calcEachBaseEnergyBase(
			string aminoAcid, vector<vector<vector<string> > > &bases,
			int flag) {
		// 結果格納ベクターを初期化
		int baseLength = aminoAcid.size() * 3;
		int maxLength = 8;
		vector<vector<vector<vector<pair<int, string> > > > > result;
		initResultExtendedBaseEnergyVector(baseLength, maxLength, result);

		// 各ポジション（配列長=5, 6, 8）における部分配列の自由エネルギーを取得
		vector<int> baseLengths;
		baseLengths.push_back(5);
		baseLengths.push_back(6);
		baseLengths.push_back(8);
		for (unsigned int size = 0; size < baseLengths.size(); size++) {
			int baseSize = baseLengths[size];
			for (int i = 0; i < baseLength; i++) {
				int position = i + 1;
				for (unsigned int k = 0; k < bases[position][baseSize].size();
						k++) {
					string base = bases[position][baseSize][k];

					// 配列、エネルギーマップから自由エネルギーを検索
					string originalBase = base;
					if (flag == BASE_EXTENDED) {
						// 拡張した塩基を使用する場合は、V、WはUに変換、X、YはGに変換
						Util::baseReplace(originalBase, "V", "U");
						Util::baseReplace(originalBase, "W", "U");
						Util::baseReplace(originalBase, "X", "G");
						Util::baseReplace(originalBase, "Y", "G");
					}

					// エネルギー取得（未設定なら以降の処理スキップ）
					map<string, int>::iterator itr;
					itr = baseEnergy.find(originalBase);
					if (itr == baseEnergy.end()) {
						continue;
					}
					int energy = itr->second;

					// 結果ベクターに格納
					int startBase = getBaseNumber(base.substr(0, 1));
					int endBase = getBaseNumber(
							base.substr(base.length() - 1, 1));
					int resultEnergy =
							result[position][baseSize][startBase][endBase].first;
					string resultBase =
							result[position][baseSize][startBase][endBase].second;

					if (resultBase == "") {
						result[position][baseSize][startBase][endBase].first =
								energy;
						result[position][baseSize][startBase][endBase].second =
								base;
					} else {
						if (energy < resultEnergy) {
							result[position][baseSize][startBase][endBase].first =
									energy;
							result[position][baseSize][startBase][endBase].second =
									base;
						}
					}

				}
			}
		}
		return result;
	}

	vector<vector<vector<vector<pair<int, string> > > > > calcEachExtendedBaseEnergy(
			string aminoAcid, vector<vector<vector<string> > > &bases) {
		return calcEachBaseEnergyBase(aminoAcid, bases, BASE_EXTENDED);
	}

	vector<vector<vector<vector<pair<int, string> > > > > calcEachOriginalBaseEnergy(
			string aminoAcid, vector<vector<vector<string> > > &bases) {
		return calcEachBaseEnergyBase(aminoAcid, bases, BASE_ORIGINAL);
	}

	vector<vector<vector<string> > > getSelectedLengthExtendedBases(string aminoAcid,
			unsigned int maxLength, string exceptedCodons) {
		vector<vector<vector<string> > > result;
		getBasesBase(aminoAcid, maxLength, result, BASE_EXTENDED,
				exceptedCodons);
		return result;
	}

	vector<vector<vector<string> > > getSelectedLengthOriginalBases(string aminoAcid,
			unsigned int maxLength, string exceptedCodons) {
		vector<vector<vector<string> > > result;
		getBasesBase(aminoAcid, maxLength, result, BASE_ORIGINAL,
				exceptedCodons);
		return result;
	}

	/*
	 * デバッグ用
	 */
	void printTime(string msg, clock_t start, clock_t end) {
		double margin = (double) (end - start) / CLOCKS_PER_SEC;
		cout << msg << "\t" << margin << endl;
	}

};

#endif /* AASEQCONVERTER_H_ */

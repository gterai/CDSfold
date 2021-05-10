#include <iostream>
#include <stdlib.h>
//#include <string> // なくても良い。
#include <fstream>
#include <vector>
#include <string.h> // to use strcpy for gcc-4.3 or later
using namespace std;

class eachseq{
public:
  char *desc;
  char *seq;
  int seqlen;
};

class fasta{
public:
  fasta(const char *fname); //コンストラクタには戻り値
                            //が要らないという特別ルールがある。
  ~fasta();
  char *getDesc(){
    return data[p].desc;
  }
  char *getSeq(){
    return data[p].seq;
  }
  int getSeqLen(){
    return data[p].seqlen;
  }
  void initP(){
    p = 0;
  }
  void printP(){
    cout << p << endl;
  }
  
  int next(){
    p++;
    if(p > numSeq){
      p = 0;
      return 0;
    }
    return 1;
  }
  
private:
  vector<eachseq> data;    // いくつ配列があるか分からないので、Vector型とする。
  int p;      // ポインタ
  int numSeq; // 配列数
};


fasta::fasta(const char *fname){ //DPマトリクスの動的メモリ確保
  
  initP();
  numSeq = 0; 
  
  ifstream ifs(fname);
  if(ifs){
    string line;
    string tmp_desc;
    string tmp_seq;
    //int id = 0;
    getline(ifs, line);
    if(line[0] != '>'){
      cerr << "Invalid fasta format." << endl;
      exit(1);
    }
    else{
      line.erase(0,1); // ">"を削除
      tmp_desc = line;
    }
    while(getline(ifs, line)){ // 行の読み込み
      //cout << line << endl;
      if(line[0] == '>'){
	eachseq e;
	//cout << "ok1" << endl;
	e.seq  = new char[tmp_seq.length()+1];     // 配列の長さ分だけ領域を新たに確保する +1入れないとなぜかsegmentaionエラーに。 <- うーん。良く分からん。
	e.desc = new char[tmp_desc.length()+1];    // 配列の長さ分だけ領域を新たに確保する
	strcpy(e.seq, (char*)tmp_seq.c_str());   // 文字列をコピー（char型への変換が必要）
	strcpy(e.desc, (char*)tmp_desc.c_str()); // 文字列をコピー（char型への変換が必要）
	e.seqlen = tmp_seq.length();
	data.push_back(e); // 
	tmp_seq.erase(0); // 空にする
 	tmp_desc.erase(0);// 空にする

	// e.seq should be deleted?

	//id++;
	numSeq++; // 配列数
	
	line.erase(0,1); // ">"を削除
	tmp_desc = line;
      }
      else{
	tmp_seq += line;
      }
      
    }
    // 一番最後の配列
    eachseq e;
    e.seq  = new char[tmp_seq.length()+1];
    e.desc = new char[tmp_desc.length()+1];
    strcpy(e.seq, (char*)tmp_seq.c_str());
    strcpy(e.desc, (char*)tmp_desc.c_str());
	e.seqlen = tmp_seq.length();
    data.push_back(e);


    // e.seq should be deleted?
    //tmp_seq.erase(0); // this should be done?
    //tmp_desc.erase(0);// this should be done?
    
  }
  else{
	  if(fname == NULL){
		  cerr << "Error: no input file" << endl;
	  }
	  else{
		  cerr << "Error: cannot open file(" << fname << ")" << endl;
	  }
	  exit(1);

  }
  
}

fasta::~fasta(){
  //cout << "# elemetns :" << data.size() << endl;
  for(unsigned int i = 0; i < data.size(); i++){
    //cout << "deleting fasta: " << i << endl;
    delete [] data[i].desc;
    delete [] data[i].seq;
  }
}

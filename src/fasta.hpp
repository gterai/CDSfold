#include <iostream>
#include <stdlib.h>
//#include <string> // $B$J$/$F$bNI$$!#(B
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
  fasta(const char *fname); //$B%3%s%9%H%i%/%?$K$OLa$jCM(B
                            //$B$,MW$i$J$$$H$$$&FCJL%k!<%k$,$"$k!#(B
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
  vector<eachseq> data;    // $B$$$/$DG[Ns$,$"$k$+J,$+$i$J$$$N$G!"(BVector$B7?$H$9$k!#(B
  int p;      // $B%]%$%s%?(B
  int numSeq; // $BG[Ns?t(B
};


fasta::fasta(const char *fname){ //DP$B%^%H%j%/%9$NF0E*%a%b%j3NJ](B
  
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
      line.erase(0,1); // ">"$B$r:o=|(B
      tmp_desc = line;
    }
    while(getline(ifs, line)){ // $B9T$NFI$_9~$_(B
      //cout << line << endl;
      if(line[0] == '>'){
	eachseq e;
	//cout << "ok1" << endl;
	e.seq  = new char[tmp_seq.length()+1];     // $BG[Ns$ND9$5J,$@$1NN0h$r?7$?$K3NJ]$9$k(B +1$BF~$l$J$$$H$J$<$+(Bsegmentaion$B%(%i!<$K!#(B <- $B$&!<$s!#NI$/J,$+$i$s!#(B
	e.desc = new char[tmp_desc.length()+1];    // $BG[Ns$ND9$5J,$@$1NN0h$r?7$?$K3NJ]$9$k(B
	strcpy(e.seq, (char*)tmp_seq.c_str());   // $BJ8;zNs$r%3%T!<!J(Bchar$B7?$X$NJQ49$,I,MW!K(B
	strcpy(e.desc, (char*)tmp_desc.c_str()); // $BJ8;zNs$r%3%T!<!J(Bchar$B7?$X$NJQ49$,I,MW!K(B
	e.seqlen = tmp_seq.length();
	data.push_back(e); // 
	tmp_seq.erase(0); // $B6u$K$9$k(B
 	tmp_desc.erase(0);// $B6u$K$9$k(B

	// e.seq should be deleted?

	//id++;
	numSeq++; // $BG[Ns?t(B
	
	line.erase(0,1); // ">"$B$r:o=|(B
	tmp_desc = line;
      }
      else{
	tmp_seq += line;
      }
      
    }
    // $B0lHV:G8e$NG[Ns(B
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

/*
 * Util.h
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <string>
#include <map>
#include <iostream>

using namespace std;

class Util {
public:
	static void baseReplace(string &base, string from, string to){
		int pos = 0;
		while((pos = base.find(from,pos)) != string::npos){
			base.replace(pos,from.length(),to);
			pos += to.length();
		}
	}
};

#endif /* UTIL_H_ */

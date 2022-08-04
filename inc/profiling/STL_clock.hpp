#pragma once

#include <string>			// for C++ string
#include <vector>			// for vector container
#include <chrono>
using namespace std;

class stopWatch {
	chrono::high_resolution_clock::time_point t_start, t_stop;
	double total=0; 
public:
	void start() {
		t_start = chrono::high_resolution_clock::now();
	}
	void stop() {
		t_stop = chrono::high_resolution_clock::now();
		chrono::duration<double> d = t_stop - t_start;
		total += d.count();
	}

	double elapsedTime() {
		return total;
	}
	static double resolution() {
		auto tmp = chrono::high_resolution_clock::period();
		return (double)tmp.num / tmp.den;
	}
};



struct clockstruct {
	stopWatch beginNew;
	stopWatch beginMainloop;
	stopWatch convectionDifussion;
	stopWatch p;
	stopWatch updateT1toT3;
	stopWatch checkL2norm;
	stopWatch updateT3toT0;
	stopWatch BC;
	stopWatch set_MaA;
};
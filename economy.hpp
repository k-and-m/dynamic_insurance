/*
 * economy.hpp
 *
 *  Created on: May 28, 2015
 *      Author: arun
 */

#ifndef ECONOMY_HPP_
#define ECONOMY_HPP_

#include "BaseHeader.h"
#include "households.h"
#include <string>
#include <fstream>
#include <sstream>
#include <random>
using namespace std;

class Economy{
public:
	Economy(int numHouseholds, int seed);
	~Economy();
	void initialize(const EquilFns &pol, const StochProc& proc, const State& p_currentState);
	void simulateOnePeriod(int newPhiState, double r, double agg1, double agg2);
	void testOnePeriod(int newPhiState, double r, double agg1, double agg2, double randNum) const;
	double getAverage(unsigned int state);
	double getAverageTest(unsigned int state) const;
	double getAverageAssets();
	void printEconomy(std::string fileName);
protected:
	Household **myHHs;
private:
	int phiState;
	const unsigned int econSize;
	uniform_real_distribution<double> distr;
	mt19937 gener;

};

#endif /* ECONOMY_HPP_ */

#pragma once
#ifndef WORLD_ECONOMY_HPP_
#define WORLD_ECONOMY_HPP_

#include "BaseHeader.h"
#include "economy.hpp"
#include "brent.hpp"
#include <random>
using namespace std;
using namespace brent;

class WorldEconomy : public func_base
{
public:
	WorldEconomy(int p_numCountries, int currentState);
	~WorldEconomy();
	void initialize(int whichCountry, const EquilFns &pol, const StochProc& proc, const State& p_currentState);
	double distance(VecDoub targets);
	void printEconomies();
	void simulateToSS();
	void simulateNPeriods(int n);
	vector<VecDoub> getHistory();
	double operator() (double r) const;

private:
	Economy **e;
	int numCountries;
	int curSt;
	int currentPeriod;
	StochProc **myStoch;
	vector < VecDoub > history;

	uniform_real_distribution<double> distr;
	mt19937 gener;

	uniform_real_distribution<double> testdistr;
	mt19937 testgener;
	int testSeed;

	double distance(VecDoub targets, int targetPhi);
	void simulateOnePeriod(double r);
};

#endif
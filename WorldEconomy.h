#pragma once
#ifndef WORLD_ECONOMY_HPP_
#define WORLD_ECONOMY_HPP_

#include "BaseHeader.h"
#include "economy.hpp"
#include <random>
using namespace std;

class WorldEconomy
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

private:
	Economy **e;
	int numCountries;
	int curSt;
	int currentPeriod;
	StochProc **myStoch;
	vector < VecDoub > history;

	uniform_real_distribution<double> distr;
	mt19937 gener;

	double distance(VecDoub targets, int targetPhi);
	void simulateOnePeriod();
};

#endif
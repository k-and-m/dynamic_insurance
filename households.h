/*
 * households.h
 *
 *  Created on: May 28, 2015
 *      Author: arun
 */

#ifndef HOUSEHOLDS_H_
#define HOUSEHOLDS_H_

#include "BaseHeader.h"
#include "utilityFunctions.h"
#include "vfiMaxUtil.h"
#include "eulerUtility.h"
#include <random>
#include <string>
#include <sstream>
using namespace std;

class Household{
public:
	Household(int seed, const EquilFns &pol, const StochProc& proc);
	~Household();
	void setInitialState(const State& p_currentState);
	void setAggState(int newAggState, int newPhiState);
	void iterate(int newAggState, int newPhiState);
	void test(int newAggState, int newPhiState, double r);
	double getCurrentState(int whichState);
	double getPreviousState(int whichState);
	double getCurrentAsset(int whichAsset);
	string toString();

protected:
	int randSeed;
	const EquilFns *valAndPol;
	State oldState;
	VecDoub oldAssetDist;
	State currentState;
	VecDoub currentAssetDist;
	const StochProc *s_proc;
	uniform_real_distribution<double> distr;
	mt19937 gener;

private:
	vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>>& getReformPolicyFn();
	static vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>>& refPolInit(const EquilFns *param);
};

#endif /* HOUSEHOLDS_H_ */

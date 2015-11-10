#pragma once

#include "typedefs.h"

class StochProc
{
public:
	vector<vector<vector<vector<vector<VecDoub>>>>> shocks; // By defining this separately, we can parallelize the problem.
	vector<vector<vector<vector<vector<Mat5Doub>>>>> transition; // this transition is across all stochastic states
	vector<vector<vector<vector<vector<Mat5Doub>>>>> cdf;  // this cdf is conditional on aggregate states, i.e. phi, agg
	double phiTrans[PHI_STATES][PHI_STATES];
	double phiTransCDF[PHI_STATES][PHI_STATES];

#if ASSETSTATE
	VecDoub assets;
	VecDoub aggAssets;
#else
	double capital[CAPITAL_SIZE];
	double bonds[BOND_SIZE];
#endif
	StochProc(const StochProc& init);
	StochProc(VecDoub phis);
	~StochProc();
	int getCondNewPhi(const int curSt, const double randNum) const;
	VecInt getCondNewState(const VecInt& currentState, const int newAggState, const int newPhi, const double randNum) const;
	StochProc& operator=(const StochProc& fnSource);
};
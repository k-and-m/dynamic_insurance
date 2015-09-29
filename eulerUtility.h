#pragma once
#ifndef EULER_UTILITY_H_
#define EULER_UTILITY_H_

#include "BaseHeader.h"
#include "vfiMaxUtil.h"

class eulerUtility : public vfiMaxUtil
{
private:

protected:
	vector<vector<vector<vector<vector<VecDoub>>>>> valueFnDeriv;

public:
	eulerUtility(const State& current, const StochProc& stoch, const EquilFns& fns);
	~eulerUtility();
	double operator() (VecDoub state_prime) const;
};
#endif
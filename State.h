#pragma once
#include "typedefs.h"
#include "StochProc.h"

class State
{
public:
	VecDoub current_states;
	VecInt current_indices;
#if ASSETSTATE==0
	alglib::spline3dinterpolant s;
#endif
	VecDoub phi;
	VecDoub prices;
	vector<vector<VecDoub>> coefficients;

	State();
	State(const VecDoub& phi1, const VecDoub& prices);
	State(const State& orig);
	State& operator=(const State& fnSource);
	double getRecursiveVal(int whichVal) const;
	double getTau() const;

	void defaultInitialState(StochProc& stoch2);
};


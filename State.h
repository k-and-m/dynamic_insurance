#pragma once
#include "typedefs.h"
#include "StochProc.h"
#include "matrix.h"

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

	State();
	State(const VecDoub& phi1, const VecDoub& prices);
	State(const State& orig);
	State& operator=(const State& fnSource);
	double getTau() const;
	void defaultInitialState(StochProc& stoch2);
	double getNextR() const;
};
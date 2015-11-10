#pragma once
#include "typedefs.h"
#include "StochProc.h"
#include "matrix.h"

class State
{
private:
	Mat3Doub coefficients;

public:
	VecDoub current_states;
	VecInt current_indices;
#if ASSETSTATE==0
	alglib::spline3dinterpolant s;
#endif
	VecDoub phi;
	VecDoub prices;

	State();
	State(const VecDoub& phi1, const VecDoub& prices, const Mat3Doub& recursEst);
	State(const State& orig);
	State& operator=(const State& fnSource);
	double getRecursiveVal(int whichVal) const;
	double getTau() const;
	void defaultInitialState(StochProc& stoch2);
	double getNextR(double a1, double a2, int phi_state) const;

};


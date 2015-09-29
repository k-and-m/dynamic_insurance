#if 0
#pragma once
#include "typedefs.h"
class EndogenousGrid
{
private:
	const EquilFns *valAndPol;
	State currentState;
	StochProc *s_proc;

public:
	EndogenousGrid();
	~EndogenousGrid();
	void solve();
};

#endif
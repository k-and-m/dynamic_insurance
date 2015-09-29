#pragma once

#include "typedefs.h"
class EquilFns
{
public:
#if ASSETSTATE
	vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>> value_fn;
	vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>> consumption;
	vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>> policy_fn;
#else
	VecDoub *******value_fn; //Could be just a double and haven't changed yet
	VecDoub *******consumption;//Could be just a double and haven't changed yet
	VecDoub *******policy_fn;
#endif

	EquilFns();
	~EquilFns();
	EquilFns& operator=(const EquilFns& fnSource);
};


#include "utilityBase.h"


utilityBase::utilityBase(const State& current, const StochProc& stoch, const EquilFns& fns) :
curSt(&current), curStoch(stoch), curFns(fns)
{
}


utilityBase::~utilityBase()
{
}

double utilityBase::consUtil(double consumption){
	return utilityFunctions::integer_power(consumption, 1 - RRA) / (1 - RRA);
}

double utilityBase::marginalConsUtil(double consumption){
	return utilityFunctions::integer_power(consumption, -RRA);
}

void utilityBase::updateCurrent(const State& current){
	curSt = &current;
}


#pragma once
#ifndef UTILITY_BASE_H_
#define UTILITY_BASE_H_

#include "BaseHeader.h"
#include "utilityFunctions.h"

class utilityBase
{
protected:
	const State* curSt;
	const StochProc& curStoch;
	const EquilFns& curFns;

public:
	utilityBase(const State& current, const StochProc& stoch, const EquilFns& fns);
	~utilityBase();
	virtual double operator() (VecDoub state_prime) const = 0;
	virtual void updateCurrent(const State& current) final;
	virtual bool constraintBinds() const = 0;
	virtual	double getBoundBorrow() const = 0;
	virtual	double getBoundUtil() const = 0;
	static double consUtil(double consumption);
	static double marginalConsUtil(double consumption);
};

#endif
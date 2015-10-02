#pragma once
#include "BaseHeader.h"
#include "utilityBase.h"
class vfiMaxUtil :
	public utilityBase
{
protected:
	double get_wage(const State& current, const StochProc& stoch) const;
	double get_zshock(const State& current, const StochProc& stoch, COUNTRYID whichCountry) const;
	double prod_fn(const double k1, const double k2, /*const double c1mgmt,*/ const State& current,
		const StochProc& stoch) const;
	double getMaxBorrow(const double k1, const double k2/*, const double c1mgmt*/) const;

public:
	vfiMaxUtil(const State& current, const StochProc& stoch, const EquilFns& fns);
	~vfiMaxUtil();
	virtual double operator() (VecDoub state_prime) const;
	double getBoundUtil() const final;
	double getBoundBorrow() const final;
	bool constraintBinds() const final;
	double calculateCurrentAssets(const double k1, const double k2, const double bonds, double r) const;
	double calculateTestAssets(const double k1, const double k2, const double bonds, const double r) const;
	double getNextPeriodAggAssets(int whichCountry, const double currentAggAsst) const;
};


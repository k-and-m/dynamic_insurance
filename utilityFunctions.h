#pragma once
#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include "typedefs.h"

namespace utilityFunctions{
	double boundValue(const double origVal, const double lower, const double upper);
	double interpolate(const VecDoub &x, const VecDoub &f, double a_prime);
	double interpolate2d(const VecDoub &x1, const VecDoub &x2, const vector<VecDoub> &y, double x1_prime, double x2_prime);
	double interpolate3d(const VecDoub &x1, const VecDoub &x2, const VecDoub &x3, const vector<vector<VecDoub>> &y, double x1_prime, double x2_prime, double x3_prime);
	double integer_power(const double base, const int exponent);
	double distanceFn(VecDoub& x, VecDoub& y);
}

#endif
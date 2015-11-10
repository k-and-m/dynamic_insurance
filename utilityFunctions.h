#pragma once
#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include "typedefs.h"
#include <sys/stat.h>

namespace utilityFunctions{
	double boundValue(const double origVal, const double lower, const double upper);
	double interpolate(const VecDoub &x, const VecDoub::const_iterator f, double a_prime);
	double integer_power(const double base, const int exponent);
	double distanceFn(VecDoub& x, VecDoub& y);
	bool fileExists(const std::string& filename);
}

#endif
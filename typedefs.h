/*
 * typedefs.h
 *
 *  Created on: Mar 4, 2015
 *      Author: robmrk
 */

#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include "defines.h"
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cstddef>
#include <cctype>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include <algorithm>
#include "matrix.h"
#include "matrixIO.h"
#include <vector>
using std::vector;

//Numerical recipes type definitions for Amoeba/Golden

typedef std::vector<int> VecInt, VecInt_O, VecInt_IO;
typedef const std::vector<int> VecInt_I;
typedef std::vector<double> VecDoub, VecDoub_O, VecDoub_IO;
typedef const std::vector<double> VecDoub_I;

typedef Numeric_lib::Matrix<int, 2> MatInt, MatInt_O, MatInt_IO;
typedef const Numeric_lib::Matrix<int, 2> MatInt_I;
typedef Numeric_lib::Matrix<double, 2> MatDoub, MatDoub_O, MatDoub_IO;
typedef const Numeric_lib::Matrix<double, 2> MatDoub_I;

typedef Numeric_lib::Matrix<double, 3> Mat3Doub;
typedef Numeric_lib::Matrix<double, 4> Mat4Doub;
typedef Numeric_lib::Matrix<double, 5> Mat5Doub;

enum COUNTRYID {
	C1 = 0, C2 = 1
};

#endif /* TYPEDEFS_H_ */


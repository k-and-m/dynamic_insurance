/*
 * numericalrecipes.h
 *
 *  Created on: Feb 11, 2015
 *      Author: robmrk
 */

#ifndef NUMERICALRECIPES_H_
#define NUMERICALRECIPES_H_
#include "BaseHeader.h"

double golden(double ax, double bx, double cx, State& current /*InitialStruct& initial EquilsFns& fns */, double (*f)(VecDoub, State& /*,InitialStruct& ,EquilsFns& */), double tol, double *xmin);
void tauchen(double mean, double rho, double sigma, double m, double (&z)[CAP_SHOCK_SIZE], double (&zprob)[CAP_SHOCK_SIZE][CAP_SHOCK_SIZE]);
void tauchen2(double mean, double rho, double sigma, double m, double (&z)[WAGE_SHOCK_SIZE], double (&zprob)[WAGE_SHOCK_SIZE][WAGE_SHOCK_SIZE]);

#endif /* NUMERICALRECIPES_H_ */

/*
 * numericalrecipes.cc
 *
 *  Created on: Feb 11, 2015
 *      Author: robmrk
 */


#include <math.h>
#include "numericalrecipes.h"

#define R 0.61803399
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
double golden(double ax, double bx, double cx, State& current, //InitialStruct& initial EquilsFns& fns
		double (*f)(VecDoub, State&/*,InitialStruct& ,EquilsFns& */), double tol, double *xmin)
{
	double f1,f2,x0,x1,x2,x3;
	VecDoub tempParam(1);
	x0=ax;
	x3=cx;
	if (fabs(cx-bx) > fabs(bx-ax))
	{
		x1=bx;
		x2=bx+C*(cx-bx);
	}
	else
	{
		x2=bx;
		x1=bx-C*(bx-ax);
	}
	tempParam[0]=x1;
	f1=(*f)(tempParam, current /*,initial , fns */);
	tempParam[0]=x2;
	f2=(*f)(tempParam, current /*,initial , fns */);

	int counter=0;
	while ((counter<1000)&&(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))))
	{
		counter++;
		if (f2 < f1)
		{
			SHFT3(x0,x1,x2,R*x1+C*x3)
			tempParam[0]=x2;
			SHFT2(f1,f2,(*f)(tempParam, current /*,initial , fns */))
		}
		else
		{
			SHFT3(x3,x2,x1,R*x2+C*x0)
			tempParam[0]=x1;
			SHFT2(f2,f1,(*f)(tempParam, current /*,initial , fns */))
		}
	}
	if (counter==1000)
	{
		std::cout<<"Error: Golden exceeded max counter!!"<<std::endl;
		exit(-1);
	}
	if (f1 < f2)
	{
		*xmin=x1;
		return f1;
	}
	else
	{
		*xmin=x2;
		return f2;
	}
}



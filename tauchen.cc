/*
 * tauchen.cc
 *
 *  Created on: Mar 12, 2015
 *      Author: robmrk
 */

/*function [Z,Zprob] = tauchen(mu,rho,sigma,m)
%Function TAUCHEN
%
%Purpose:    Finds a Markov chain whose sample paths
%            approximate those of the AR(1) process
%                z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
%            where eps are normal with stddev sigma
%
%Format:     {Z, Zprob} = Tauchen(mu,rho,sigma,m)
%
%Input:
%            mu      scalar, unconditional mean of process
%            rho     scalar
%            sigma   scalar, std. dev. of epsilons
%            m       max +- std. devs.
%
%Output:     Z       SHOCK_SIZE*1 vector, nodes for Z
%            Zprob   SHOCK_SIZE*SHOCK_SIZE matrix, transition probabilities
%
%    Martin Floden
%    Fall 1996
%
%    This procedure is an implementation of George Tauchen's algorithm
%    described in Ec. Letters 20 (1986) 177-181.
%
*/

#include <iostream>
#include <cstdlib>
#include <cstddef>
#include <cctype>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <fstream>
#include <cassert>
#include "numericalrecipes.h"
#include <algorithm>
#include "typedefs.h"

double cdf_normal(double x);

void tauchen(double mean, double rho, double sigma, double m, double (&z)[CAP_SHOCK_SIZE], double (&zprob)[CAP_SHOCK_SIZE][CAP_SHOCK_SIZE])
{
	double a;
	double zstep=0;

	if (CAP_SHOCK_SIZE<0)
	{
		std::cout<<"Tauchen method, CAP_SHOCK_SIZE must be greater than zero. SHOCK_SIZE= "<<CAP_SHOCK_SIZE<<std::endl;
		exit(-1);
	}

	for(int i=0; i<CAP_SHOCK_SIZE; i++)
	{
		z[i]=0;
		for(int j=0; j<CAP_SHOCK_SIZE; j++)
		{
			zprob[i][j]=0;
		}
	}

	a=(1-rho)*mean;

	if(CAP_SHOCK_SIZE>1)
	{
		z[CAP_SHOCK_SIZE-1]= m*sqrt(pow(sigma,2)/(1-pow(rho,2)));
		z[0]= -z[CAP_SHOCK_SIZE-1];
		zstep=(z[CAP_SHOCK_SIZE-1]-z[0])/(CAP_SHOCK_SIZE-1);
	}
	else
	{
		z[0]=0;
	}

	for (int i=1; i<CAP_SHOCK_SIZE-1; i++)
	{
		z[i]=z[0]+zstep*i;
	}

	for (int i=0; i<CAP_SHOCK_SIZE; i++)
	{
		z[i]=z[i]+a/(1-rho);
	}

	if(CAP_SHOCK_SIZE==1)
	{
		zprob[0][0]=1;
		return;
	}
	for (int j=0; j<CAP_SHOCK_SIZE; j++)
	{
		for(int k=0; k<CAP_SHOCK_SIZE; k++)
		{
			if (k==0)
			{
				zprob[j][k]=cdf_normal((z[0]-a-rho*z[j]+zstep/2)/sigma);
			}
			else if (k==CAP_SHOCK_SIZE-1)
			{
				zprob[j][k]=1-cdf_normal((z[CAP_SHOCK_SIZE-1]-a-rho*z[j]-zstep/2)/sigma);
			}
			else
			{
				zprob[j][k]=cdf_normal((z[k]-a-rho*z[j]+zstep/2)/sigma)-cdf_normal((z[k]-a-rho*z[j]-zstep/2)/sigma);
			}
		}
	}
}

double cdf_normal(double x)
{
	return 0.5 * erfc(-x/sqrt(2));
}



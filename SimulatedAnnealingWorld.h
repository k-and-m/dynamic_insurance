#pragma once
#ifndef SIM_ANNEAL_WORLD_HPP
#define SIM_ANNEAL_WORLD_HPP 1

#include "typedefs.h"
#include "simAnneal.hpp"
#include "brent.hpp"

class SimulatedAnnealingWorld :
	public SimulatedAnnealing<VecDoub, double>
{
public:

	//      void printStatus(double solution, double temp);
	VecDoub* giveRandomNeighbour(const VecDoub& lastSolution) const;

	double calcDistanceToTarget(const VecDoub& solution) const;

	bool shouldStopHook(const VecDoub& solution);

	SimulatedAnnealingWorld(const VecDoub& startSolution, const double& target,
		double starttemp, double precision, brent::func_base& f, double p_maxstep, int max_mult);

	~SimulatedAnnealingWorld() {
	}
	;

	double getEnergy(const VecDoub& lastSolution) const {
		return (*func)(lastSolution[0]);
	}
	;


private:
	double unifRand() const {
		return rand() / double(RAND_MAX);
	};

	double unifRand(double a, double b) const {
		double x;
		x = (b - a) * unifRand() + a;
		return x;
	};

	brent::func_base *func;
	const double maxstep;
	long int numCalls=0;
	int maxMult;
};
#endif
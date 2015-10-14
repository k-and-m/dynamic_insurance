#include "SimulatedAnnealingWorld.h"

VecDoub* SimulatedAnnealingWorld::giveRandomNeighbour(
	const VecDoub& lastPoint) const {
	VecDoub *newPoint = new VecDoub(lastPoint.size());
	*newPoint = lastPoint;
	for (unsigned int i = 0; i < lastPoint.size(); i++) {
		double dist = (unifRand(-maxstep, maxstep));
		(*newPoint)[i] = lastPoint[i] + dist;
	}
	return newPoint;
}

double SimulatedAnnealingWorld::calcDistanceToTarget(const VecDoub& angle) const {
	double distance = getEnergy(angle) - *TARGET;
	if (distance > 0) {
		return distance;
	}
	else {
		return (-distance);
	}
}

bool SimulatedAnnealingWorld::shouldStopHook(const VecDoub& solution) {
	numCalls++;
	return numCalls > 10*maxMult;
}

SimulatedAnnealingWorld::SimulatedAnnealingWorld(const VecDoub& startSolution, const double& target,
	double starttemp, double precision, brent::func_base& f, double p_maxstep, int max_mult) :
	SimulatedAnnealing<VecDoub, double>(startSolution, target, starttemp, precision), maxstep(p_maxstep),
	maxMult(max_mult){
	static long int seed = 1;
	srand(seed++);
	func = &f;
}

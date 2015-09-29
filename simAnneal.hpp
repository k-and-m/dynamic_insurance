/*
 * simAnneal.hpp
 *
 *  Created on: Apr 27, 2015
 *      Author: arun
 */

#ifndef SIMANNEAL_HPP_
#define SIMANNEAL_HPP_

#include <iostream>	// needed for basic IO
#include <cmath>	// needed for chance calculation
#include <assert.h> // will use assert to check certain values
#include <cstdlib>
#include <cstddef>

#include "typedefs.h"

/***********************************************************************************************//**

 \brief Simulated Annealing Framework:
 A class to help users implement simulated annealing in an easy and flexible way.

 Simulated Annealing Class

 This class offers the user a flexibel framework to use for the implementation of applications
 of Simulated Annealing. To use it the user needs to create a child class that implements
 the 3 required functions listed below.

 Require overriding:
 - giveRandomNeighbour()
 - calcDistanceToTarget()

 Optionale overriding:
 - shouldStopHook()
 - printStatus()
 - calcProbability()
 - calcNewTemp()



 The type of the candidate solutions can be set using template parameters.
 The type of the target value we want to achieve can also be set using template parameters.
 The link between these two types is made by the calcDistanceToTarget function which is one of
 the functions that require overriding.

 ***************************************************************************************************/

template<class Solution, class Target>
class SimulatedAnnealing {

public:

	/**
	 Constructor
	 @param startSolution The solution in the domain where you want to start the search
	 @param TARGET The TARGET value you are looking for
	 @param starttemp The starttemperature, the higher this temperature, the longer the
	 algorythm will allow uphill moves
	 @param PRECISION The required PRECISION for a solution to be acceptable
	 */
	SimulatedAnnealing(const Solution& startSolution, const Target& target,
			double starttemp, double precision);

	/**
	 The public method that is called from a SimulatedAnnealing(or child class)-Object
	 */
	Solution* solve();

	/**
	 This function returns the energy of the given solution
	 @param lastSolution
	 */
	virtual double getEnergy(const Solution& lastSolution) const =0;

protected:
        SimulatedAnnealing(){};
	virtual ~SimulatedAnnealing();

	/***********************************************************************************************

	 Following functions HAVE to be implemented to adapt the algorythm to the specific problem

	 ***********************************************************************************************/

	/**
	 This function returns a new random neighbour, how the neihbour is calculated is entirely
	 up to the user (and can influence the result of the algorythm a lot)
	 @param lastSolution The previous solution which can optionally (and should) be used to
	 determine the new neighbour
	 @return A pointer to the newly chosen neighbour (you don't have to worry about deleting)
	 */
	virtual Solution* giveRandomNeighbour(
			const Solution& lastSolution) const =0;

	/**
	 This function determines how we can compare the target value and the current solution
	 their types aren't required to be the same, as long as this function can determine a
	 double which symbolises the distance.
	 @param solution The previous solution of which the distance to the target needs to be
	 determined
	 @return The value that symbolises the distance to the target
	 */
	virtual double calcDistanceToTarget(const Solution& solution) const =0;

	/***********************************************************************************************

	 Following functions CAN be adapted to your needs

	 ***********************************************************************************************/

	/**
	 This hook determines whether or not to go on with the search, it extends the already
	 present stop condition that stops execution when the distance of a solution to the
	 target is smaller than the precision required.

	 Standard implementation always returns false.
	 @param solution The solution that will be considered as final solution
	 @return A boolean value determining whether or not to stop
	 */
	virtual bool shouldStopHook(const Solution& solution);

	/**
	 This function can be used to print some info during execution of the algorythm.
	 @param solution The solution that was last considered
	 @param temp The current temperature
	 */
	virtual void printStatus(const Solution& solution, double temp);

	/**
	 This function can be used to change the way the probability for acceptance of uphill
	 moves is calculated. It's a function considering the improvement of the distance to the TARGET
	 of the old and the new solution and the temperature.

	 CAUTION: The probability value returned will be checked, probabilities adhere to the following
	 rule \f$ 0<=probability<=1 \f$

	 Standard implementation:
	 \f$
	 P(\Delta(f(s)), T) = exp(\frac{\Delta(f(s))}{T})
	 \f$

	 @param change The change in distance to TARGET between the last and the new solution
	 @param temp The current temperature
	 */
	virtual double calcProbability(double change, double temp) const;

	/**
	 This function calculates the new temperature based on the previous temperature.

	 CAUTION: The temperature value returned will be checked, it has to be a positive value.
	 @param lastTemp The previous temperature which now needs to be updated
	 @return The newly calculated temperature
	 */
	virtual double calcNewTemp(double lastTemp) const;

	/***********************************************************************************************

	 Following functions are standard and SHOULD NOT be overridden

	 ***********************************************************************************************/

	bool shouldStop(const Solution& solution);

	double calcDistanceChange(const Solution& oldSolution,
			const Solution& newSolution) const;

	bool accept(const Solution& oldSolution, const Solution& newSolution,
			double temp) const;

	const Target* TARGET;
	const double PRECISION;
	const double ALPHA = 0.9;

	Solution* solution;
	Solution* best;
	double temp;

};

template<class Solution, class Target>
bool SimulatedAnnealing<Solution, Target>::shouldStop(
		const Solution& solution) {

	if (shouldStopHook(solution)) {
//		std::cout << std::endl << "STOP REASON: ShouldStopHook" << std::endl;
		return true;
	} else if (calcDistanceToTarget(solution) < PRECISION) {
		std::cout << std::endl
				<< "STOP REASON: Distance to target is smaller than the required precision. Solution found."
				<< std::endl;
		return true;
	} else {
		return false;
	}

}

template<class Solution, class Target>
double SimulatedAnnealing<Solution, Target>::calcDistanceChange(
		const Solution& oldSolution, const Solution& newSolution) const {

	double oldDistance = calcDistanceToTarget(oldSolution);
	double newDistance = calcDistanceToTarget(newSolution);
	assert(oldDistance >= 0 && newDistance >= 0); // distances are always positive
	return (newDistance - oldDistance);

}

template<class Solution, class Target>
double SimulatedAnnealing<Solution, Target>::calcProbability(double change,
		double temp) const { // should be overwritten

	return exp(-1.0 * change / temp);

}

template<class Solution, class Target>
double SimulatedAnnealing<Solution, Target>::calcNewTemp(
		double lastTemp) const {
	return lastTemp * ALPHA;
}

template<class Solution, class Target>
bool SimulatedAnnealing<Solution, Target>::accept(const Solution& oldSolution,
		const Solution& newSolution, double temp) const {

	//std::cout << "Testing acceptance: Old Solution: " << oldSolution << " New Solution: " << newSolution << std::endl;
	double change = calcDistanceChange(oldSolution, newSolution);
	//std::cout << "\t ->Change: " << change << std::endl;
	if (change < 0) {
//		std::cout << "Result is better so: ";
		return true;
	} else {
		double probability = calcProbability(change, temp);
		//DEBUG:BEGIN
		if (!(probability >= 0 && probability <= 1)) {
			std::cout << "Probability wrong: " << probability << " temp: "
					<< temp << std::endl;
			std::cin.get();
		}
		//DEBUG:END
		assert(probability >= 0 && probability <= 1); // probability has to be checked
//		std::cout << "Probability we're gonna accept: " << probability << " result: ";
		return (rand() % ((int) (1.0 / probability)) == 0);
	}

}

template<class Solution, class Target>
Solution* SimulatedAnnealing<Solution, Target>::solve() {

	printStatus(*solution, temp);
	long int counter = 0;
	while (!shouldStop(*solution)) {
//		std::cout<<counter++<<std::endl<<std::flush;
		Solution *solPtr = solution;
		Solution* newSolution = giveRandomNeighbour(*solution);
		if (accept(*solution, *newSolution, temp)) {
			solution = newSolution;
//			std::cout << "ACCEPTED" << std::endl<<std::flush;
		} else {
//			std::cout << "REJECTED" << std::endl<<std::flush;
			delete newSolution;
		}
//		std::cout << "energy(solution)" << std::endl<<std::flush;
		double newEn = -getEnergy(*solution);
//		std::cout << "energy(best)" << std::endl<<std::flush;
		double bestEn = -getEnergy(*best);
		if (newEn > bestEn) {
			if (solPtr == best) {
//				std::cout << "delete solPtr" << std::endl<<std::flush;
				delete solPtr;
			} else {
//				std::cout << "delete solPtr2" << std::endl<<std::flush;
				delete solPtr;
//				std::cout << "delete best" << std::endl<<std::flush;
				delete best;
			}
			best = solution;
		} else {
			if ((solPtr != best) && (solPtr != solution)) {
//				std::cout << "delete solPtr3" << std::endl<<std::flush;
				delete solPtr;
			} else {
//				std::cout << "don't delete solPtr3" << std::endl<<std::flush;
			}
		}
//		std::cout << "pointer deletion done?" << std::endl<<std::flush;

		//std::cout << "\n*****************\n" << std::endl;
		printStatus(*solution, temp);
		temp = calcNewTemp(temp);
		assert(temp >= 0); // only positive temperatures are allowed!
	}
//	std::cout << "\n\n\n*************************************************************\n\n"
//		<< "We're done, Solution: \n\n" << *solution << "\n\n*************************************************************\n" << std::endl;

//	std::cout<<"in solve"<<(best==solution)<<std::endl<<std::flush;

	if (solution != best) {
		delete solution;
	}
	delete TARGET;

//	std::cout << "\n\n\nPress enter to continue..." << std::endl;

//	std::cin.get();
	return best;
}

template<class Solution, class Target>
bool SimulatedAnnealing<Solution, Target>::shouldStopHook(
		const Solution& solution) {
	return false;
}

template<class Solution, class Target>
void SimulatedAnnealing<Solution, Target>::printStatus(const Solution& solution,
		double temp) {
//	std::cout << "Current solution: " << solution << " at Temp: " << temp << std::endl;
}

template<class Solution, class Target>
SimulatedAnnealing<Solution, Target>::SimulatedAnnealing(
		const Solution& startSolution, const Target& target, double starttemp,
		double precision) :
		TARGET(new Target(target)), solution(new Solution(startSolution)), temp(
				starttemp), PRECISION(precision) {
	assert(starttemp >= 0); // only positive temperatures are allowed!
	best = solution;
}

template<class Solution, class Target>
SimulatedAnnealing<Solution, Target>::~SimulatedAnnealing() {
}

class SimulatedAnnealingWorld: public SimulatedAnnealing<VecDoub, double> {
public:

//      void printStatus(double solution, double temp);
	VecDoub* giveRandomNeighbour(const VecDoub& lastSolution) const;

	double calcDistanceToTarget(const VecDoub& solution) const;

	bool shouldStopHook(const VecDoub& solution);

	SimulatedAnnealingWorld(const VecDoub& startSolution, const double& target,
			double starttemp, double precision, double (*function)(double,double,double,double,double,double),
                        double p_maxstep, double p_phi1, double p_phi2, double p_prop);

	~SimulatedAnnealingWorld() {
	}
	;

	double getEnergy(const VecDoub& lastSolution) const {
		return (*func)(phi1,phi2,lastSolution[0],lastSolution[1],lastSolution[2],prop);
	}
	;

private:
	double unifRand() const {
		return rand() / double(RAND_MAX);
	}
	;

	double unifRand(double a, double b) const {
		double x;
		x = (b - a) * unifRand() + a;
		return x;
	}
	;

	double (*func)(double,double,double,double,double,double);
        const double maxstep;
        const double phi1;
        const double phi2;
        const double prop;
};

VecDoub* SimulatedAnnealingWorld::giveRandomNeighbour(
		const VecDoub& lastPoint) const {
	VecDoub *newPoint = new VecDoub;
	*newPoint = lastPoint;
	for (unsigned int i = 0; i < lastPoint.size(); i++) {
		(*newPoint)[i] = lastPoint[i] + unifRand(-maxstep, maxstep);
	}
	return newPoint;
}

double SimulatedAnnealingWorld::calcDistanceToTarget(
		const VecDoub& angle) const {
	static long int count = 0;
	count += 1;
	double distance = -getEnergy(angle) - *TARGET;
	if (distance > 0) {
		return distance;
	} else {
		return (-distance);
	}
}

bool SimulatedAnnealingWorld::shouldStopHook(const VecDoub& solution) {
	return temp < MINIMIZATION_TOL;
}

SimulatedAnnealingWorld::SimulatedAnnealingWorld(const VecDoub& startSolution,
					const double& target, double starttemp, double precision, 
                			double (*function)(double,double,double,double,double,double),
                                        double p_maxstep, double p_phi1, double p_phi2, double p_prop) :
		SimulatedAnnealing<VecDoub,double>(startSolution, target, starttemp, precision),maxstep(p_maxstep),
                phi1(p_phi1),phi2(p_phi2),prop(p_prop) {
		static long int seed = 1;
	        srand(seed++);
		func = function;
}

#endif /* SIMANNEAL_HPP_ */

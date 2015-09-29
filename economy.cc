/*
 * economy.cc
 *
 *  Created on: May 28, 2015
 *      Author: arun
 */
#include "economy.hpp"
//static int randSeed = 1;

Economy::Economy(int numHouseholds, int seed):myHHs(),econSize(numHouseholds){
	myHHs=new Household*[numHouseholds];

	distr = uniform_real_distribution<double>(0.0, 1.0);
	gener = mt19937(seed);
//	srand(randSeed++);
}

Economy::~Economy(){
	for (unsigned int i = 0; i < econSize; i++){
		delete myHHs[i];
	}
	delete[] myHHs;
}
void Economy::initialize(const EquilFns &pol, const StochProc& proc, const State& p_currentState){
	for (unsigned int i = 0; i < econSize; i++){
		myHHs[i] = new Household(floor(10*NUMHHS*distr(gener)), pol, proc);
//		myHHs[i] = new Household(rand(), pol, proc);
	}

	for (unsigned int i = 0; i < econSize; i++){
		myHHs[i]->setInitialState(p_currentState);
	}
}

double Economy::getAverage(unsigned int asset){
	double total=0;
	for (unsigned int i = 0; i < econSize; i++){
		total+=myHHs[i]->getCurrentAsset(asset);
	}
	return total/econSize;
}

double Economy::getAverageAssets(){
	double total = 0;
	for (unsigned int i = 0; i < econSize; i++){
		total += myHHs[i]->getCurrentState(ASTATE);
	}
	return total / econSize;
}

void Economy::simulateOnePeriod(int newPhiState){
	phiState = newPhiState;
#if OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(3)
#endif
	for (int i = 0; i < (int)econSize; i++){
		myHHs[i]->iterate(0,newPhiState);
	}
}

void Economy::printEconomy(std::string fileName){
	ofstream out_stream;
	out_stream.open(fileName.c_str());

	out_stream
			<< "hh,z1,z2,wage_shock,k1,k2,b"
			<< endl;

	for (unsigned int i = 0; i < econSize; i++){
		out_stream<<i<<","<<myHHs[i]->toString()<<endl;
	}
	out_stream.close();
}


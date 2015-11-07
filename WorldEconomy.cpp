#include "WorldEconomy.h"


WorldEconomy::WorldEconomy(int p_numCountries, int currentState) : numCountries(p_numCountries), curSt(currentState)
{
	e = new Economy*[numCountries];
	myStoch = new StochProc*[numCountries];

	distr = uniform_real_distribution<double>(0.0, 1.0);
	gener = mt19937(WORLD_ECON_SEED);

	testdistr = uniform_real_distribution<double>(0.0, 1.0);
	testgener = mt19937(WORLD_ECON_SEED);

	history.resize(TOTALPERIODS+1);
	for (int i = 0; i < TOTALPERIODS+1; i++){
		history[i].resize(p_numCountries + 3);
	}
	currentPeriod = 0;
}


WorldEconomy::~WorldEconomy()
{
	for (int i = 0; i < numCountries; i++){
		delete e[i];
	}
	delete[] e;

	delete myStoch;
}

void WorldEconomy::initialize(int whichCountry, const EquilFns &pol, const StochProc& proc, const State& p_currentState)
{
	stateVar = p_currentState;

	myStoch[whichCountry] = new StochProc(proc);
	curSt = p_currentState.current_indices[PHI_STATE];

	e[whichCountry] = new Economy(NUMHHS, floor(numCountries*10*distr(gener)));
	e[whichCountry]->initialize(pol,proc,p_currentState);
	currentPeriod = 0;

	history[currentPeriod][numCountries] = curSt;
	history[currentPeriod][whichCountry] = e[whichCountry]->getAverageAssets();

	history[currentPeriod][numCountries + 2] += e[whichCountry]->getAverage(BSTATE);
	if (whichCountry == 1) {
		history[currentPeriod][numCountries + 1] = p_currentState.getNextR();
	}
}

void WorldEconomy::simulateToSS(){
	simulateNPeriods(SSPERIODS);
}

void WorldEconomy::simulateNPeriods(int n)
{
	for (int index = 0; index < n; index++){
		testSeed = testdistr(testgener);
#if 1
		double r = history[currentPeriod][numCountries + 1];
#elif 0
		for (int j = 0; j < 1000; j++) {
			double r = j*3.0 / 1000 + 1;
			double netBonds = (*this)(r);
			std::cout << "NetBonds: " << netBonds << " , r=" << r << std::endl;
		}
		exit(-1);
#else
		const double targ = 0;
		VecDoub start = VecDoub(1);
		start[0] = r;
		std::cout << "Simulate again." << std::endl;
		SimulatedAnnealingWorld worldSolve = SimulatedAnnealingWorld(start, targ, 100, .01, *this, 1, currentPeriod+1);
		VecDoub *soln = worldSolve.solve();
		r = (*soln)[0];
#endif
		simulateOnePeriod(r);
	}
}
void WorldEconomy::simulateOnePeriod(double r)
{
	currentPeriod += 1;
	double randNum = distr(gener);
	curSt = myStoch[0]->getCondNewPhi(curSt, randNum);
	history[currentPeriod][numCountries] = curSt;

	double netBonds = 0;
	for (int i = 0; i < numCountries; i++){
		e[i]->simulateOnePeriod(curSt, r, e[0]->getAverageAssets(),e[1]->getAverageAssets());
		history[currentPeriod][i] = e[i]->getAverageAssets();
		netBonds += e[i]->getAverage(BSTATE);
	}
	history[currentPeriod][numCountries + 1] = stateVar.getNextR();
	history[currentPeriod][numCountries + 2] = netBonds;
	std::cout << "Period " << currentPeriod << ", phi: " << curSt << " , A1: " << history[currentPeriod][0] << ", A2: " << history[currentPeriod][1] << ", Net Bonds: " << netBonds << " , r: " << r << std::endl;
}

double WorldEconomy::operator() (double r) const
{
	int localcurSt = curSt;
	double retVal = 0;
	for (int i = 0; i < numCountries; i++) {
		e[i]->testOnePeriod(localcurSt, r, e[0]->getAverageAssets(), e[1]->getAverageAssets(), testSeed+i);
		retVal += e[i]->getAverageTest(BSTATE);
	}

	return retVal;
}

double WorldEconomy::netBonds()
{
	double dist = 0;
	for (int i = 0; i < numCountries; i++) {
		dist += e[i]->getAverage(BSTATE);
	}
	return dist;
}

void WorldEconomy::printEconomies(){
	for (int i = 0; i < numCountries; i++){
		ostringstream os;
		os << "economy_" << i << ".txt";
		e[i]->printEconomy(os.str());
	}
}

vector<VecDoub> WorldEconomy::getHistory(){
	return history;
}
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
		history[currentPeriod][numCountries + 1] = p_currentState.getNextR(history[0][0], history[0][1], curSt);
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
		e[i]->simulateOnePeriod(curSt, r, MAX(MIN_AGG_ASSETS,e[0]->getAverageAssets()), 
			MAX(MIN_AGG_ASSETS,e[1]->getAverageAssets()));
		history[currentPeriod][i] = MAX(MIN_AGG_ASSETS, e[i]->getAverageAssets());
		netBonds += e[i]->getAverage(BSTATE);
	}
	history[currentPeriod][numCountries + 1] = stateVar.getNextR(history[currentPeriod][0], history[currentPeriod][1],curSt);
	history[currentPeriod][numCountries + 2] = netBonds;
	std::cout << "Period " << currentPeriod << ", phi: " << curSt << " , A1: " << history[currentPeriod][0] << ", A2: " << history[currentPeriod][1] << ", Net Bonds: " << netBonds << " , r: " << r << std::endl;
}

double WorldEconomy::operator() (double r) const
{
	int localcurSt = curSt;
	double retVal = 0;
	for (int i = 0; i < numCountries; i++) {
		e[i]->testOnePeriod(localcurSt, r, MAX(MIN_AGG_ASSETS, e[0]->getAverageAssets()),
			MAX(MIN_AGG_ASSETS, e[1]->getAverageAssets()), testSeed+i);
		retVal += e[i]->getAverageTest(BSTATE);
	}

	return retVal;
}
/*
double WorldEconomy::distance(VecDoub targets)
{
	if (targets.size() != numCountries){
		std::cerr << "WorldEconomy.distance(targets): Different size between countries and target proportions" << std::endl;
		exit(-1);
	}

	double dist = 0;
	for (int i = 0; i < PHI_STATES; i++){
		dist += distance(targets, i);
	}

	return dist;
}

double WorldEconomy::distance(VecDoub targets, int targetPhi)
{
	curSt = targetPhi;
	for (int i = 0; i < numCountries; i++){
		testSeed = testdistr(testgener);
#if 1
		double r = zero(-1, 1, 0.00001, *this);
#else
		double r = 0;
		double netBonds = glomin(-1, 1, 0, 1000, 0.001, 0.0001, *this, r);
#endif
		e[i]->simulateOnePeriod(targetPhi,r);
	}

	//find distance from market clearing
	double totalK1 = 0;
	double totalK2 = 0;
	double totalB = 0;

	std::cout << "State = " << targetPhi << std::endl;
	std::cout << "Country                K1                 K2                  B                  TotalAssets" << std::endl;
	for (int i = 0; i < numCountries; i++){
		totalK1 += e[i]->getAverage(K1STATE);
		totalK2 += e[i]->getAverage(K2STATE);
		totalB += e[i]->getAverage(BSTATE);
		std::cout << i + 1 << "(" << myStoch[i]->shocks[targetPhi][0][0][0][0][EF_PHI] << ")               " << e[i]->getAverage(K1STATE)
			           << "          " << e[i]->getAverage(K2STATE) 
					   << "          " << e[i]->getAverage(BSTATE)
					   << "          " << e[i]->getAverageAssets()
					   << std::endl;
	}

	double c1Prop = targets[0];
	double c2Prop = targets[1];

	std::cout << "total                " << totalK1 << "          " << totalK2 << "          " << totalB << endl;
	std::cout << "target                  " << c1Prop << "              " << c2Prop << "                 0" << endl;

	double k1d = 10 * (utilityFunctions::integer_power(1 + (abs(totalK1 / numCountries - c1Prop)) / c1Prop, 2) - 1);
	double k2d = 10 * (utilityFunctions::integer_power(1 + (abs(totalK2 / numCountries - c2Prop)) / c2Prop, 2) - 1);
	double bd = utilityFunctions::integer_power(1 + abs(totalB) / 0.1, 2) - 1;

	std::cout << "distance              " << k1d << "          " << k2d << "         " << bd << std::endl << std::endl;

	return sqrt(k1d + k2d + bd);
}
*/
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
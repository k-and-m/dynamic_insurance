#include "WorldEconomy.h"


WorldEconomy::WorldEconomy(int p_numCountries, int currentState) : numCountries(p_numCountries), curSt(currentState)
{
	e = new Economy*[numCountries];
	myStoch = new StochProc*[numCountries];

	distr = uniform_real_distribution<double>(0.0, 1.0);
	gener = mt19937(WORLD_ECON_SEED);

	testdistr = uniform_real_distribution<double>(0.0, 1.0);
	testgener = mt19937(WORLD_ECON_SEED);

	history.resize(TOTALPERIODS);
	for (int i = 0; i < TOTALPERIODS; i++){
		history[i].resize(p_numCountries + 2);
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
	myStoch[whichCountry] = new StochProc(proc);
	curSt = p_currentState.current_indices[PHI_STATE];

	e[whichCountry] = new Economy(NUMHHS, floor(numCountries*10*distr(gener)));
	e[whichCountry]->initialize(pol,proc,p_currentState);
	currentPeriod = 0;

}

void WorldEconomy::simulateToSS(){
	simulateNPeriods(SSPERIODS);
}

void WorldEconomy::simulateNPeriods(int n)
{
	for (int index = 0; index < n; index++){
		testSeed = testdistr(testgener);
#if 1
		double r=zero(-1, 1, 0.00001, *this);
#else
		double r = 0;
		double netBonds = glomin(-1, 1, 0, 1000, 0.001, 0.0001, *this, r);
#endif
		simulateOnePeriod(r);
	}
}

void WorldEconomy::simulateOnePeriod(double r)
{
	double randNum = distr(gener);
	curSt = myStoch[0]->getCondNewPhi(curSt, randNum);
	history[currentPeriod][numCountries] = curSt;
	history[currentPeriod][numCountries+1] = r;
	double netBonds = 0;
	for (int i = 0; i < numCountries; i++){
		e[i]->simulateOnePeriod(curSt, r, e[0]->getAverageAssets(), e[1]->getAverageAssets());
		history[currentPeriod][i] = e[i]->getAverageAssets();
		netBonds += e[i]->getAverage(BSTATE);
	}
	std::cout << "Period " << currentPeriod << " , Net Bonds: " << netBonds << " , r: " << r << std::endl;
	currentPeriod += 1;
}

double WorldEconomy::operator() (double r) const
{
	int localcurSt = myStoch[0]->getCondNewPhi(curSt, testSeed);
	for (int i = 0; i < numCountries; i++) {
		e[i]->testOnePeriod(localcurSt, r, e[0]->getAverageAssets(), e[1]->getAverageAssets(), testSeed + i);
	}

	double retVal = 0;
	for (int i = 0; i < numCountries; i++) {
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
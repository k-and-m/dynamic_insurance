/*
 * households.cc
 *
 *  Created on: May 28, 2015
 *      Author: arun
 */
#include "households.h"
Household::Household(int seed, const EquilFns &pol, const StochProc& proc) :
randSeed(seed){
	valAndPol = &pol;
	s_proc = &proc;
	distr = uniform_real_distribution<double>(0.0, 1.0);
	gener = mt19937(seed);
	oldAssetDist = VecDoub(NUM_CHOICE_VARS);
	currentAssetDist = VecDoub(NUM_CHOICE_VARS);
	testAssetDist = VecDoub(NUM_CHOICE_VARS);

	refPolInit(valAndPol);
}

vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>>& Household::getReformPolicyFn(){
	return refPolInit(valAndPol);
}

vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>>& Household::refPolInit(const EquilFns *param,
	bool reset){
	static vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>> reform_policy_fn;
	static bool reform_pol_set = false;

	if ((reform_pol_set == false)&&(param!=NULL)){
		std::cout << "Re-ordering policy function arrays." << std::endl << std::flush;
		reform_pol_set = true;
		reform_policy_fn.resize(PHI_STATES);
		for (int f = 0; f < PHI_STATES; f++){
			reform_policy_fn[f].resize(AGG_SHOCK_SIZE);
			for (int g = 0; g < AGG_SHOCK_SIZE; g++){
				reform_policy_fn[f][g].resize(CAP_SHOCK_SIZE);
				for (int h = 0; h < CAP_SHOCK_SIZE; h++) {
					reform_policy_fn[f][g][h].resize(CAP_SHOCK_SIZE);
					for (int i = 0; i < CAP_SHOCK_SIZE; i++) {
						reform_policy_fn[f][g][h][i].resize(WAGE_SHOCK_SIZE);
						for (int ii = 0; ii < WAGE_SHOCK_SIZE; ii++) {
#if K2CHOICE
							reform_policy_fn[f][g][h][i][ii].resize(NUM_CHOICE_VARS);
							for (int l = 0; l < NUM_CHOICE_VARS; l++) {
#else
							reform_policy_fn[f][g][h][i][ii].resize(NUM_CHOICE_VARS-1);
							for (int l = 0; l < NUM_CHOICE_VARS-1; l++) {
#endif
								reform_policy_fn[f][g][h][i][ii][l].resize(AGG_ASSET_SIZE);
								for (int ll = 0; ll < AGG_ASSET_SIZE; ll++) {
									reform_policy_fn[f][g][h][i][ii][l][ll].resize(AGG_ASSET_SIZE);
									for (int m = 0; m < AGG_ASSET_SIZE; m++) {
										reform_policy_fn[f][g][h][i][ii][l][ll][m].resize(ASSET_SIZE);
										for (int n = 0; n < ASSET_SIZE; n++) {
#if K2CHOICE
											reform_policy_fn[f][g][h][i][ii][l][ll][m][n] = param->policy_fn[ll][m][f][g][n][h][i][ii][l];
#else
											if (l == K1STATE) {
												reform_policy_fn[f][g][h][i][ii][l][ll][m][n] = param->policy_fn[ll][m][f][g][n][h][i][ii][l];
											}
											else {
												reform_policy_fn[f][g][h][i][ii][l][ll][m][n] = param->policy_fn[ll][m][f][g][n][h][i][ii][l+1];
											}
#endif
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	if (reset == true) {
		reform_pol_set = false;
	}

	return reform_policy_fn;
}

Household::~Household(){
	refPolInit(NULL, true);
}

void Household::setInitialState(const State& p_currentState) {
	oldState = p_currentState;
	currentState = p_currentState;

	currentAssetDist[K1STATE] = p_currentState.current_states[ASTATE] / (4 * NOMINAL_PRICE * p_currentState.getTau());
	currentAssetDist[K2STATE] = p_currentState.current_states[ASTATE] / (4 * NOMINAL_PRICE);
	currentAssetDist[BSTATE] = p_currentState.current_states[ASTATE] / 4;

	oldAssetDist[K1STATE] = currentAssetDist[K1STATE];
	oldAssetDist[K2STATE] = currentAssetDist[K2STATE];
	oldAssetDist[BSTATE] = currentAssetDist[BSTATE];

}

void Household::setRandomInitialState(const State& p_currentState) {
	oldState = p_currentState;
	currentState = p_currentState;
	double randNum = distr(gener);
	currentState.current_indices[ASTATE] = floor(randNum*ASSET_SIZE);
	currentState.current_states[ASTATE] = s_proc->assets[currentState.current_indices[ASTATE]];
	currentAssetDist[K1STATE] = currentState.current_states[ASTATE] / (4*NOMINAL_PRICE*p_currentState.getTau());
	currentAssetDist[K2STATE] = currentState.current_states[ASTATE] / (4*NOMINAL_PRICE);
	currentAssetDist[BSTATE] = currentState.current_states[ASTATE] / 4;

	oldAssetDist[K1STATE] = currentAssetDist[K1STATE];
	oldAssetDist[K2STATE] = currentAssetDist[K2STATE];
	oldAssetDist[BSTATE] = currentAssetDist[BSTATE];
}

void Household::setAggState(int newAggState, int newPhiState) {
	cerr << "need to complete code for household.setAggState." << endl;
	exit(-1);
}

void Household::iterate(int newAggState, int newPhiState, double r, const double agg1, const double agg2) {
	oldState = currentState;
	oldAssetDist[K1STATE] = currentAssetDist[K1STATE];
	oldAssetDist[K2STATE] = currentAssetDist[K2STATE];
	oldAssetDist[BSTATE] = currentAssetDist[BSTATE];

	VecInt curSt = VecInt(NUM_SHOCK_VARS);
	curSt[EF_K1] = oldState.current_indices[CAP1_SHOCK_STATE];
	curSt[EF_K2] = oldState.current_indices[CAP2_SHOCK_STATE];
	curSt[EF_W] = oldState.current_indices[WAGE_SHOCK_STATE];
	curSt[EF_A] = oldState.current_indices[AGG_SHOCK_STATE];
	curSt[EF_PHI] = oldState.current_indices[PHI_STATE];

	double k1 = oldAssetDist[K1STATE];
	double k2 = oldAssetDist[K2STATE];
	double b = oldAssetDist[BSTATE];
	k1 = MAX(k1, MIN_CAPITAL);
	k2 = MAX(k2, MIN_CAPITAL);

	double randNum = distr(gener);
	VecInt newSt = s_proc->getCondNewState(curSt, newAggState, newPhiState, randNum);

	currentState.current_indices[CAP1_SHOCK_STATE] = newSt[EF_K1];
	currentState.current_indices[CAP2_SHOCK_STATE] = newSt[EF_K2];
	currentState.current_indices[WAGE_SHOCK_STATE] = newSt[EF_W];
	currentState.current_indices[AGG_SHOCK_STATE] = newSt[EF_A];
	currentState.current_indices[PHI_STATE] = newSt[EF_PHI];

	currentState.current_states[CAP1_SHOCK_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_K1];
	currentState.current_states[CAP2_SHOCK_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_K2];
	currentState.current_states[WAGE_SHOCK_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_W];
	currentState.current_states[AGG_SHOCK_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_A];
	currentState.current_states[PHI_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_PHI];

	vfiMaxUtil vfiMU = vfiMaxUtil(currentState,*s_proc,*valAndPol);

	currentState.current_states[ASTATE] = vfiMU.calculateCurrentAssets(k1, k2, b, r);
	if (currentState.current_states[ASTATE] != currentState.current_states[ASTATE]){
		cout << "ERROR! Household.cc, interate() : currentState.current_state[ASTATE]=NaN" << endl;
		exit(-1);
	}
	if (currentState.current_states[ASTATE] > MAX_ASSETS){
		currentState.current_states[ASTATE] = MAX_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (currentState.current_states[ASTATE] < MIN_ASSETS){
		currentState.current_states[ASTATE] = MIN_ASSETS;
	}

	oldState.current_states[AGG_ASSET_STATE] = agg1;
	oldState.current_states[AGG2_ASSET_STATE] = agg2;

	currentState.current_states[AGG_ASSET_STATE] = oldState.getRecursiveVal(AGG_ASSET_C1);
	if (currentState.current_states[AGG_ASSET_STATE] != currentState.current_states[AGG_ASSET_STATE]){
		cout << "ERROR! Household.cc, interate() : currentState.current_state[AGG_ASSET_STATE]=NaN" << endl;
		exit(-1);
	}
	if (currentState.current_states[AGG_ASSET_STATE] > MAX_AGG_ASSETS){
		currentState.current_states[AGG_ASSET_STATE] = MAX_AGG_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (currentState.current_states[AGG_ASSET_STATE] < MIN_AGG_ASSETS){
		currentState.current_states[AGG_ASSET_STATE] = MIN_AGG_ASSETS;
	}

	currentState.current_states[AGG_ASSET_STATE] = oldState.getRecursiveVal(AGG_ASSET_C2);
	if (currentState.current_states[AGG2_ASSET_STATE] != currentState.current_states[AGG2_ASSET_STATE]){
		cout << "ERROR! Household.cc, interate() : currentState.current_state[AGG2_ASSET_STATE]=NaN" << endl;
		exit(-1);
	}
	if (currentState.current_states[AGG2_ASSET_STATE] > MAX_AGG_ASSETS){
		currentState.current_states[AGG2_ASSET_STATE] = MAX_AGG_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (currentState.current_states[AGG2_ASSET_STATE] < MIN_AGG_ASSETS){
		currentState.current_states[AGG2_ASSET_STATE] = MIN_AGG_ASSETS;
	}

#if K2CHOICE
	double k1prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K1STATE]).begin(),
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double k2prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K2STATE]).begin(),
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double bprime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][BSTATE]).begin(),
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
#else
	double k1prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K1STATE]).begin(),
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double k2prime = k1prime * pow(oldState.getTau(), 1/(1 - ALPHA1));
	double bprime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][BSTATE-1]).begin(),
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
#endif

	currentAssetDist[K1STATE] = k1prime;
	currentAssetDist[K2STATE] = k2prime;
	currentAssetDist[BSTATE] = bprime;
}
void Household::test(int newAggState, int newPhiState, double r, const double agg1, const double agg2, double randNum){
	
	VecInt curSt = VecInt(NUM_SHOCK_VARS);
	curSt[EF_K1] = currentState.current_indices[CAP1_SHOCK_STATE];
	curSt[EF_K2] = currentState.current_indices[CAP2_SHOCK_STATE];
	curSt[EF_W] = currentState.current_indices[WAGE_SHOCK_STATE];
	curSt[EF_A] = currentState.current_indices[AGG_SHOCK_STATE];
	curSt[EF_PHI] = currentState.current_indices[PHI_STATE];

	double k1 = currentAssetDist[K1STATE];
	double k2 = currentAssetDist[K2STATE];
	double b = currentAssetDist[BSTATE];
	k1 = MAX(k1, MIN_CAPITAL);
	k2 = MAX(k2, MIN_CAPITAL);

	VecInt newSt = s_proc->getCondNewState(curSt, newAggState, newPhiState, randNum);
	State localCurrentState(currentState);

	localCurrentState.current_indices[CAP1_SHOCK_STATE] = newSt[EF_K1];
	localCurrentState.current_indices[CAP2_SHOCK_STATE] = newSt[EF_K2];
	localCurrentState.current_indices[WAGE_SHOCK_STATE] = newSt[EF_W];
	localCurrentState.current_indices[AGG_SHOCK_STATE] = newSt[EF_A];
	localCurrentState.current_indices[PHI_STATE] = newSt[EF_PHI];

	localCurrentState.current_states[CAP1_SHOCK_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_K1];
	localCurrentState.current_states[CAP2_SHOCK_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_K2];
	localCurrentState.current_states[WAGE_SHOCK_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_W];
	localCurrentState.current_states[AGG_SHOCK_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_A];
	localCurrentState.current_states[PHI_STATE] =
		s_proc->shocks[newSt[EF_PHI]][newSt[EF_A]][newSt[EF_K1]][newSt[EF_K2]][newSt[EF_W]][EF_PHI];

	vfiMaxUtil vfiMU = vfiMaxUtil(localCurrentState, *s_proc, *valAndPol);

	localCurrentState.current_states[ASTATE] = vfiMU.calculateCurrentAssets(k1, k2, b, r);
	if (localCurrentState.current_states[ASTATE] != localCurrentState.current_states[ASTATE]) {
		cout << "ERROR! Household.cc, interate() : currentState.current_state[ASTATE]=NaN" << endl;
		exit(-1);
	}
	if (localCurrentState.current_states[ASTATE] > MAX_ASSETS) {
		localCurrentState.current_states[ASTATE] = MAX_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (localCurrentState.current_states[ASTATE] < MIN_ASSETS) {
		localCurrentState.current_states[ASTATE] = MIN_ASSETS;
	}

	localCurrentState.current_states[AGG_ASSET_STATE] = agg1;
	localCurrentState.current_states[AGG2_ASSET_STATE] = agg2;

	double temp1 = localCurrentState.getRecursiveVal(AGG_ASSET_C1);
	double temp2 = localCurrentState.getRecursiveVal(AGG_ASSET_C2);
	localCurrentState.current_states[AGG_ASSET_STATE] = temp1;
	if (localCurrentState.current_states[AGG_ASSET_STATE] != localCurrentState.current_states[AGG_ASSET_STATE]) {
		cout << "ERROR! Household.cc, interate() : localCurrentState.current_state[AGG_ASSET_STATE]=NaN" << endl;
		exit(-1);
	}
	if (localCurrentState.current_states[AGG_ASSET_STATE] > MAX_AGG_ASSETS) {
		localCurrentState.current_states[AGG_ASSET_STATE] = MAX_AGG_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (localCurrentState.current_states[AGG_ASSET_STATE] < MIN_AGG_ASSETS) {
		localCurrentState.current_states[AGG_ASSET_STATE] = MIN_AGG_ASSETS;
	}

	localCurrentState.current_states[AGG2_ASSET_STATE] = temp2;
	if (localCurrentState.current_states[AGG2_ASSET_STATE] != localCurrentState.current_states[AGG2_ASSET_STATE]) {
		cout << "ERROR! Household.cc, interate() : currentState.current_state[AGG2_ASSET_STATE]=NaN" << endl;
		exit(-1);
	}
	if (localCurrentState.current_states[AGG2_ASSET_STATE] > MAX_AGG_ASSETS) {
		localCurrentState.current_states[AGG2_ASSET_STATE] = MAX_AGG_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (localCurrentState.current_states[AGG2_ASSET_STATE] < MIN_AGG_ASSETS) {
		localCurrentState.current_states[AGG2_ASSET_STATE] = MIN_AGG_ASSETS;
	}

#if K2CHOICE
	double k1prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K1STATE]).begin(),
		localCurrentState.current_states[AGG_ASSET_STATE], localCurrentState.current_states[AGG2_ASSET_STATE], localCurrentState.current_states[ASTATE]);
	double k2prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K2STATE]).begin(),
		localCurrentState.current_states[AGG_ASSET_STATE], localCurrentState.current_states[AGG2_ASSET_STATE], localCurrentState.current_states[ASTATE]);
	double bprime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][BSTATE]).begin(),
		localCurrentState.current_states[AGG_ASSET_STATE], localCurrentState.current_states[AGG2_ASSET_STATE], localCurrentState.current_states[ASTATE]);
#else
	double k1prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K1STATE]).begin(),
		localCurrentState.current_states[AGG_ASSET_STATE], localCurrentState.current_states[AGG2_ASSET_STATE], localCurrentState.current_states[ASTATE]);
	double k2prime = k1prime * pow(oldState.getTau(), 1/(1 - ALPHA1));
	double bprime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, (getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][BSTATE-1]).begin(),
		localCurrentState.current_states[AGG_ASSET_STATE], localCurrentState.current_states[AGG2_ASSET_STATE], localCurrentState.current_states[ASTATE]);
#endif
	testAssetDist[K1STATE] = k1prime;
	testAssetDist[K2STATE] = k2prime;
	testAssetDist[BSTATE] = bprime;
}

/* This returns the actual "state", not k1, k2, and b levels*/
double Household::getCurrentState(int whichState) {
	return currentState.current_states[whichState];
}

/* This returns the k1, k2, and b levels*/
double Household::getCurrentAsset(int whichAsset) {
	return currentAssetDist[whichAsset];
}

/* This returns the k1, k2, and b levels*/
double Household::getTestAsset(int whichAsset) {
	return testAssetDist[whichAsset];
}

/* This returns the actual "state", not k1, k2, and b levels*/
double Household::getPreviousState(int whichState) {
	return oldState.current_states[whichState];
}


std::string Household::toString(){
	std::ostringstream os;
	os << "z1=" << getCurrentState(CAP1_SHOCK_STATE) << ",";
	os << "z2=" << getCurrentState(CAP2_SHOCK_STATE) << ",";
	os << "w =" << getCurrentState(WAGE_SHOCK_STATE) << ",";
	os << "ph=" << getCurrentState(PHI_STATE) << ",";
	os << "k1=" << getCurrentAsset(K1STATE) << ",";
	os << "k2=" << getCurrentAsset(K2STATE) << ",";
	os << "b =" << getCurrentAsset(BSTATE);

	return os.str();
}

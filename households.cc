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

	refPolInit(valAndPol);
}

vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>>& Household::getReformPolicyFn(){
	return refPolInit(valAndPol);
}

vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>>& Household::refPolInit(const EquilFns *param){
	static vector<vector<vector<vector<vector<vector<vector<vector<VecDoub>>>>>>>> reform_policy_fn;
	static bool reform_pol_set = false;

	if (reform_pol_set == false){
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
							reform_policy_fn[f][g][h][i][ii].resize(NUM_CHOICE_VARS);
							for (int l = 0; l < NUM_CHOICE_VARS; l++) {
								reform_policy_fn[f][g][h][i][ii][l].resize(AGG_ASSET_SIZE);
								for (int ll = 0; ll < AGG_ASSET_SIZE; ll++) {
									reform_policy_fn[f][g][h][i][ii][l][ll].resize(AGG_ASSET_SIZE);
									for (int m = 0; m < AGG_ASSET_SIZE; m++) {
										reform_policy_fn[f][g][h][i][ii][l][ll][m].resize(ASSET_SIZE);
										for (int n = 0; n < ASSET_SIZE; n++) {
											reform_policy_fn[f][g][h][i][ii][l][ll][m][n] = param->policy_fn[ll][m][f][g][n][h][i][ii][l];
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

	return reform_policy_fn;
}

Household::~Household(){
	//delete s_proc;
}

void Household::setInitialState(const State& p_currentState) {
	oldState = p_currentState;
	currentState = p_currentState;
	currentAssetDist[K1STATE] = MIN_CAPITAL;
	currentAssetDist[K2STATE] = MIN_CAPITAL;
	//currentAssetDist[MGMT_C1_STATE] = 1;
	currentAssetDist[BSTATE] = p_currentState.current_states[ASTATE] - 2 * MIN_CAPITAL;
}

void Household::setAggState(int newAggState, int newPhiState) {
	cerr << "need to complete code for household.setAggState." << endl;
	exit(-1);
}

void Household::iterate(int newAggState, int newPhiState) {
	oldState = currentState;
	oldAssetDist[K1STATE] = currentAssetDist[K1STATE];
	oldAssetDist[K2STATE] = currentAssetDist[K2STATE];
	oldAssetDist[BSTATE] = currentAssetDist[BSTATE];
	//oldAssetDist[MGMT_C1_STATE] = currentAssetDist[MGMT_C1_STATE];

	/*
	vector<vector<VecDoub>> k1Grid,k2Grid,bGrid;
	k1Grid.resize(AGG_ASSET_SIZE);
	k2Grid.resize(AGG_ASSET_SIZE);
	bGrid.resize(AGG_ASSET_SIZE);
	for (int h = 0; h < AGG_ASSET_SIZE; h++){
		k1Grid[h].resize(AGG_ASSET_SIZE);
		k2Grid[h].resize(AGG_ASSET_SIZE);
		bGrid[h].resize(AGG_ASSET_SIZE);
		for (int i = 0; i < AGG_ASSET_SIZE; i++){
			k1Grid[h][i].resize(ASSET_SIZE);
			k2Grid[h][i].resize(ASSET_SIZE);
			bGrid[h][i].resize(ASSET_SIZE);
		}
	}
	//VecDoub mgmtc1Grid(ASSET_SIZE);
	*/

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
	//double mgmt = oldAssetDist[MGMT_C1_STATE];
	//mgmt = MAX(MIN(mgmt, MAX_MGMT), MIN_MGMT);

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

	currentState.current_states[ASTATE] = vfiMU.calculateCurrentAssets(k1, k2, b/*, mgmt*/);
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

	currentState.current_states[AGG_ASSET_STATE] = vfiMU.getNextPeriodAggAssets(1,currentState.current_states[AGG_ASSET_STATE]);
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

	currentState.current_states[AGG2_ASSET_STATE] = vfiMU.getNextPeriodAggAssets(2,currentState.current_states[AGG2_ASSET_STATE]);
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

#if 0
	for (int g = 0; g < AGG_ASSET_SIZE; g++) {
		for (int h = 0; h < AGG_ASSET_SIZE; h++) {
			for (int i = 0; i < ASSET_SIZE; i++) {
				k1Grid[g][h][i] =
					valAndPol->policy_fn[g][h][curSt[EF_PHI]][curSt[EF_A]][i][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K1STATE];
				k2Grid[g][h][i] =
					valAndPol->policy_fn[g][h][curSt[EF_PHI]][curSt[EF_A]][i][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K2STATE];
				bGrid[g][h][i] =
					valAndPol->policy_fn[g][h][curSt[EF_PHI]][curSt[EF_A]][i][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][BSTATE];
				/*
				mgmtc1Grid[i] =
				valAndPol->policy_fn[curSt[EF_PHI]][curSt[EF_A]][i][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][MGMT_C1_STATE];
				*/
			}
		}
	}

	double k1prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, k1Grid,
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double k2prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, k2Grid,
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double bprime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, bGrid,
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
#else
	double k1prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K1STATE],
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double k2prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K2STATE],
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double bprime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][BSTATE],
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
#endif


	/*
	double mgmtprime = utilityFunctions::interpolate(s_proc->assets, mgmtc1Grid,
		currentState.current_states[ASTATE]);
		*/

	currentAssetDist[K1STATE] = k1prime;
	currentAssetDist[K2STATE] = k2prime;
	currentAssetDist[BSTATE] = bprime;
	//currentAssetDist[MGMT_C1_STATE] = mgmtprime;
}
void Household::test(int newAggState, int newPhiState, double r) {
	
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
	//double mgmt = oldAssetDist[MGMT_C1_STATE];
	//mgmt = MAX(MIN(mgmt, MAX_MGMT), MIN_MGMT);

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

	vfiMaxUtil vfiMU = vfiMaxUtil(currentState, *s_proc, *valAndPol);

	currentState.current_states[ASTATE] = vfiMU.calculateCurrentAssets(k1, k2, b/*, mgmt*/);
	if (currentState.current_states[ASTATE] != currentState.current_states[ASTATE]) {
		cout << "ERROR! Household.cc, interate() : currentState.current_state[ASTATE]=NaN" << endl;
		exit(-1);
	}
	if (currentState.current_states[ASTATE] > MAX_ASSETS) {
		currentState.current_states[ASTATE] = MAX_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (currentState.current_states[ASTATE] < MIN_ASSETS) {
		currentState.current_states[ASTATE] = MIN_ASSETS;
	}

	currentState.current_states[AGG_ASSET_STATE] = vfiMU.getNextPeriodAggAssets(1, currentState.current_states[AGG_ASSET_STATE]);
	if (currentState.current_states[AGG_ASSET_STATE] != currentState.current_states[AGG_ASSET_STATE]) {
		cout << "ERROR! Household.cc, interate() : currentState.current_state[AGG_ASSET_STATE]=NaN" << endl;
		exit(-1);
	}
	if (currentState.current_states[AGG_ASSET_STATE] > MAX_AGG_ASSETS) {
		currentState.current_states[AGG_ASSET_STATE] = MAX_AGG_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (currentState.current_states[AGG_ASSET_STATE] < MIN_AGG_ASSETS) {
		currentState.current_states[AGG_ASSET_STATE] = MIN_AGG_ASSETS;
	}

	currentState.current_states[AGG2_ASSET_STATE] = vfiMU.getNextPeriodAggAssets(2, currentState.current_states[AGG2_ASSET_STATE]);
	if (currentState.current_states[AGG2_ASSET_STATE] != currentState.current_states[AGG2_ASSET_STATE]) {
		cout << "ERROR! Household.cc, interate() : currentState.current_state[AGG2_ASSET_STATE]=NaN" << endl;
		exit(-1);
	}
	if (currentState.current_states[AGG2_ASSET_STATE] > MAX_AGG_ASSETS) {
		currentState.current_states[AGG2_ASSET_STATE] = MAX_AGG_ASSETS;
		//cout << endl << "HH hits max assets" << endl;
	}
	if (currentState.current_states[AGG2_ASSET_STATE] < MIN_AGG_ASSETS) {
		currentState.current_states[AGG2_ASSET_STATE] = MIN_AGG_ASSETS;
	}

#if 0
	for (int g = 0; g < AGG_ASSET_SIZE; g++) {
		for (int h = 0; h < AGG_ASSET_SIZE; h++) {
			for (int i = 0; i < ASSET_SIZE; i++) {
				k1Grid[g][h][i] =
					valAndPol->policy_fn[g][h][curSt[EF_PHI]][curSt[EF_A]][i][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K1STATE];
				k2Grid[g][h][i] =
					valAndPol->policy_fn[g][h][curSt[EF_PHI]][curSt[EF_A]][i][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K2STATE];
				bGrid[g][h][i] =
					valAndPol->policy_fn[g][h][curSt[EF_PHI]][curSt[EF_A]][i][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][BSTATE];
				/*
				mgmtc1Grid[i] =
				valAndPol->policy_fn[curSt[EF_PHI]][curSt[EF_A]][i][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][MGMT_C1_STATE];
				*/
			}
		}
	}

	double k1prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, k1Grid,
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double k2prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, k2Grid,
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double bprime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, bGrid,
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
#else
	double k1prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K1STATE],
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double k2prime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][K2STATE],
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
	double bprime = utilityFunctions::interpolate3d(s_proc->aggAssets, s_proc->aggAssets, s_proc->assets, getReformPolicyFn()[curSt[EF_PHI]][curSt[EF_A]][curSt[EF_K1]][curSt[EF_K2]][curSt[EF_W]][BSTATE],
		currentState.current_states[AGG_ASSET_STATE], currentState.current_states[AGG2_ASSET_STATE], currentState.current_states[ASTATE]);
#endif


	/*
	double mgmtprime = utilityFunctions::interpolate(s_proc->assets, mgmtc1Grid,
	currentState.current_states[ASTATE]);
	*/

	currentAssetDist[K1STATE] = k1prime;
	currentAssetDist[K2STATE] = k2prime;
	currentAssetDist[BSTATE] = bprime;
	//currentAssetDist[MGMT_C1_STATE] = mgmtprime;
}

/* This returns the actual "state", not k1, k2, and b levels*/
double Household::getCurrentState(int whichState) {
	return currentState.current_states[whichState];
}

/* This returns the k1, k2, and b levels*/
double Household::getCurrentAsset(int whichAsset) {
	return currentAssetDist[whichAsset];
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

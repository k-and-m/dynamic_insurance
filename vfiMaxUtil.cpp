#include "vfiMaxUtil.h"


vfiMaxUtil::vfiMaxUtil(const State& current, const StochProc& stoch, const EquilFns& fns) : utilityBase(current, stoch, fns)
{
}


vfiMaxUtil::~vfiMaxUtil()
{
}

double vfiMaxUtil::calculateCurrentAssets(const double k1, const double k2, const double bonds/*, const double mgmtprime*/) const{
	int curAgg = curSt->current_indices[AGG_SHOCK_STATE];
	int curPhi = curSt->current_indices[PHI_STATE];

	double nextAssets = get_wage(*curSt, curStoch) + (curSt->getTau() * k1) + (k2) + (1 + curSt->getRecursiveVal(P_R)) * bonds
		+ prod_fn(k1, k2, /*mgmtprime,*/ *curSt, curStoch);

	if (nextAssets != nextAssets){
		std::cout << "ERROR! vfiMaxUtil, calculateCurrentAssets() : assets=NaN" << std::endl
			<< "Wage: " << get_wage(*curSt, curStoch) << std::endl
			<< "P1: " << curSt->getTau() << std::endl
			<< "K1: " << k1 << std::endl
			<< "P2: " << 1 << std::endl
			<< "K2: " << k2 << std::endl
			<< "r: " << curSt->getRecursiveVal(P_R) << std::endl
			<< "B: " << bonds << std::endl
			//<< "M1: " << mgmtprime << std::endl
			<< "prod: " << prod_fn(k1, k2, /*mgmtprime,*/ *curSt, curStoch) << std::endl;
		exit(-1);
	}

	return nextAssets;
}

double vfiMaxUtil::getNextPeriodAggAssets(int whichCountry, const double currentAggAsst) const{
	return currentAggAsst;
}

double vfiMaxUtil::operator() (const VecDoub state_prime) const
{
	double u = 0;
	double k1prime, k2prime;
	double bprime;
	//double c1mgmt;
	double penalty = 0;
	double consumption;

	k1prime = MIN_CAPITAL + utilityFunctions::integer_power(state_prime[K1STATE], 2);
	k2prime = MIN_CAPITAL + utilityFunctions::integer_power(state_prime[K2STATE], 2);
	bprime = MIN_BONDS + utilityFunctions::integer_power(state_prime[BSTATE], 2);
	
	if (bprime != bprime){
		std::cout << "ERROR! vfiMaxUtil.cc utility() - bprime " << bprime;
		exit(-1);
	}

	//c1mgmt = MIN_MGMT + utilityFunctions::integer_power(state_prime[MGMT_C1_STATE], 2);
	//c1mgmt = utilityFunctions::boundValue(c1mgmt, MIN_MGMT, MAX_MGMT);

	double x = getMaxBorrow(k1prime, k2prime/*, c1mgmt*/);
	bprime = utilityFunctions::boundValue(bprime, x, 2 * MAX_ASSETS);

	int curAgg = curSt->current_indices[AGG_SHOCK_STATE];
	int curPhi = curSt->current_indices[PHI_STATE];

	consumption = curSt->current_states[ASTATE] - (curSt->getTau() * k1prime) - (k2prime) - bprime;

	if (consumption <= 0) {
		penalty = utilityFunctions::integer_power(1 + abs(consumption) + 0.001, 10) - 1;
		consumption = .001;
	}

	u = consUtil(consumption) - penalty;

	if (u != u) {
		std::cout.precision(15);
		std::cout << std::scientific;
		std::cout << "ERROR! vfiMaxUtil.cc utility() - utility at start = " << u
			<< ". (A,a:z1,z2)=(" << curSt->current_states[AGG_ASSET_STATE] << curSt->current_states[ASTATE] << ","
			<< curSt->current_states[CAP1_SHOCK_STATE] << ","
			<< curSt->current_states[CAP2_SHOCK_STATE] << ")"
			<< " k1prime=" << k1prime << " k2prime=" << k2prime
			<< " bprime=" << bprime << " consumption=" << consumption
			<< std::endl << " k1_input=" << state_prime[K1STATE]
			<< " k2_input=" << state_prime[K2STATE] << " b_input="
			<< state_prime[BSTATE] << std::endl << std::flush;
		std::cout << "P1:" << curSt->getTau() << std::endl
			<< "P2:" << 1 << std::endl
			<< "r:" << curSt->getRecursiveVal(P_R) << std::endl;
		exit(-1);
	}

	State newState(*curSt);
	VecDoub valueFn(ASSET_SIZE);

	int f = curSt->current_indices[AGG_ASSET_STATE];
	int ff = curSt->current_indices[AGG2_ASSET_STATE];

	for (int g = 0; g < PHI_STATES; g++){
			for (int h = 0; h < AGG_SHOCK_SIZE; h++) {
				for (int j = 0; j < WAGE_SHOCK_SIZE; j++) {
#if IID_CAP
					for (int k = 0; k < ASSET_SIZE; k++) {
						valueFn[k] = curFns.value_fn[f][ff][g][h][k][0][0][j];
					}
#endif
					for (int i = 0; i < CAP_SHOCK_SIZE; i++) {
						for (int ii = 0; ii < CAP_SHOCK_SIZE; ii++) {
#if IID_CAP==0
							for (int k = 0; k < ASSET_SIZE; k++) {
								valueFn[k] = fns.value_fn[h][k][i][ii][j][0];
							}
#endif

							newState.current_states[CAP1_SHOCK_STATE] =
								curStoch.shocks[g][h][i][ii][j][EF_K1];
							newState.current_states[CAP2_SHOCK_STATE] =
								curStoch.shocks[g][h][i][ii][j][EF_K2];
							newState.current_states[WAGE_SHOCK_STATE] =
								curStoch.shocks[g][h][i][ii][j][EF_W];
							newState.current_states[AGG_SHOCK_STATE] =
								curStoch.shocks[g][h][i][ii][j][EF_A];
							newState.current_states[PHI_STATE] =
								curStoch.shocks[g][h][i][ii][j][EF_PHI];
							newState.current_indices[CAP1_SHOCK_STATE] = i;
							newState.current_indices[CAP2_SHOCK_STATE] = ii;
							newState.current_indices[WAGE_SHOCK_STATE] = j;
							newState.current_indices[AGG_SHOCK_STATE] = h;
							newState.current_indices[PHI_STATE] = g;

							double aprime = get_wage(newState, curStoch) + (curSt->getTau() * k1prime)
								+ (k2prime)
								+ (1 + curSt->getRecursiveVal(P_R)) * bprime
								+ prod_fn(k1prime, k2prime, /*c1mgmt,*/ newState, curStoch);

							if (aprime != aprime){
								std::cout.precision(15);
								std::cout << std::scientific;
								std::cout << "ERROR! vfiMaxUtil.cc utility() - aprime for next period is NaN"
									<< std::endl
									<< "Wage: " << get_wage(newState, curStoch)
									<< "k1: " << k1prime
									<< "k2: " << k2prime
									<< "B: " << bprime
									<< "prod: " << prod_fn(k1prime, k2prime, /*c1mgmt,*/ newState, curStoch)
									<< std::endl;
								exit(-1);
							}

							int phiSt = curSt->current_indices[PHI_STATE];
							int agSt = curSt->current_indices[AGG_SHOCK_STATE];
							int cap1St = curSt->current_indices[CAP1_SHOCK_STATE];
							int cap2St = curSt->current_indices[CAP2_SHOCK_STATE];
							int wgSt = curSt->current_indices[WAGE_SHOCK_STATE];

							double intp = utilityFunctions::interpolate(curStoch.assets, valueFn, aprime);

							u +=
								BETA
								* curStoch.transition[phiSt][agSt][cap1St][cap2St][wgSt][g][h][i][ii][j]
								* MIN(0, intp);
						}
					}
				}
			}
		}
	if (u != u) {
		std::cout.precision(15);
		std::cout << std::scientific;
		std::cout << "ERROR! vfiMaxUtil.cc utility() - utility at end is " << u
			<< ". (A,a:z1,z2)=(" << curSt->current_states[AGG_ASSET_STATE] << curSt->current_states[ASTATE] << ","
			<< curSt->current_states[CAP1_SHOCK_STATE] << ","
			<< curSt->current_states[CAP2_SHOCK_STATE] << ")"
			<< " k1prime=" << k1prime << " k2prime=" << k2prime
			<< " bprime=" << bprime << " consumption=" << consumption
			<< " p1 =" << curSt->getTau() << " p2 =" << 1 << " r = " << curSt->getRecursiveVal(P_R)
			<< std::endl << std::flush;
		exit(-1);
	}
	return -u;
}

double vfiMaxUtil::prod_fn(const double k1, const double k2, /*const double c1mgmt,*/ const State& current,
	const StochProc& stoch)  const{
	double c1 = 0;
	double c2 = 0;
	double c1mgmt = (k1+k2==0)?0.5:k1 / (k1 + k2);
	c1 = current.current_states[AGG_SHOCK_STATE] * THETA
		* get_zshock(current, stoch, C1) * pow(k1, ALPHA1) * pow(c1mgmt, 1 - ALPHA1);
	c2 = current.current_states[AGG_SHOCK_STATE] * THETA
		* get_zshock(current, stoch, C2) * pow(k2, ALPHA1) * pow(1 - c1mgmt, 1 - ALPHA1);
	return c1 + c2;
}

double vfiMaxUtil::get_wage(const State& current, const StochProc& stoch)  const{
	double w;
	double exp_w = 0;
#if 0
	w = exp(current.current_states[WAGE_SHOCK_STATE]);
#else
	w = current.current_states[WAGE_SHOCK_STATE];
#endif

#if IID_CAP
	if (current.current_indices[WAGE_SHOCK_STATE] == 0){
		exp_w = 0.391;
	}
	else{
		exp_w = 1.309;
	}
#else
	int agg_state, cap1_state, cap2_state, wage_state;
	double temp_var;
	agg_state = current.current_indices[AGG_STATE];
	cap1_state = current.current_indices[CAP1_SHOCK_STATE];
	cap2_state = current.current_indices[CAP2_SHOCK_STATE];
	wage_state = current.current_indices[WAGE_SHOCK_STATE];
	std::cout << "fix this: Should only calculate zshock once since for every input since it doesn't change" << std::endl;
	exit(-1);
	for (int m = 0; m < AGG_SHOCK_SIZE; m++) {
		for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
			for (int o = 0; o < CAP_SHOCK_SIZE; o++) {
				for (int i = 0; i < WAGE_SHOCK_SIZE; i++) {
					temp_var = 5 * exp(stoch.shocks[m][n][o][i][EF_W]);
					exp_w +=
						stoch.transition[agg_state][cap1_state][cap2_state][wage_state][m][n][o][i]
						* temp_var;
				}
			}
		}
	}
#endif
	exp_w = w - exp_w;
	w = w - current.current_states[PHI_STATE] * exp_w;
	return w;
}

double vfiMaxUtil::get_zshock(const State& current, const StochProc& stoch, COUNTRYID whichCountry)  const{
	double z;
	double exp_z = 0;

#if 0
	z = exp(current.current_states[CAP1_SHOCK_STATE + whichCountry]);
#else
	z = current.current_states[CAP1_SHOCK_STATE + whichCountry];
#endif

#if IID_CAP
#if 0
	exp_z = exp((whichCountry == C1) ? CAP1_SHOCK_MEAN : CAP2_SHOCK_MEAN);
#else
	exp_z = (whichCountry == C1) ? CAP1_SHOCK_MEAN : CAP2_SHOCK_MEAN;
#endif
#else
	double temp_var1;
	int agg_state, cap1_state, cap2_state, wage_state;
	agg_state = current.current_indices[AGG_STATE];
	cap1_state = current.current_indices[CAP1_SHOCK_STATE];
	cap2_state = current.current_indices[CAP2_SHOCK_STATE];
	wage_state = current.current_indices[WAGE_SHOCK_STATE];
	std::cout << "fix this: Should only calculate zshock once since for every input since it doesn't change" << std::endl;
	exit(-1);
	for (int m = 0; m < AGG_SHOCK_SIZE; m++) {
		for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
			for (int o = 0; o < CAP_SHOCK_SIZE; o++) {
				for (int i = 0; i < WAGE_SHOCK_SIZE; i++) {
					temp_var1 = exp(
						stoch.shocks[m][n][o][i][EF_K1 + whichCountry]);
					exp_z +=
						stoch.transition[agg_state][cap1_state][cap2_state][wage_state][m][n][o][i]
						* temp_var1;
				}
			}
		}
	}
#endif
	exp_z = z - exp_z;
	z = z - current.current_states[PHI_STATE] * exp_z;
	return z;
}


double vfiMaxUtil::getBoundBorrow() const{
	return getMaxBorrow(0, 0/*, 50*/);
}

double vfiMaxUtil::getBoundUtil() const{
	VecDoub temp2 = VecDoub(NUM_CHOICE_VARS);
	temp2[K1STATE] = 0;
	temp2[K2STATE] = 0;
	temp2[BSTATE] = sqrt(getBoundBorrow() - MIN_BONDS);
	//temp2[MGMT_C1_STATE] = 7;
	return (*this)(temp2);
}

double vfiMaxUtil::getMaxBorrow(const double k1, const double k2/*, const double c1mgmt*/) const{
	State newState(*curSt);

	newState.current_states[CAP1_SHOCK_STATE] = curStoch.shocks[0][0][0][0][0][EF_K1];
	newState.current_states[CAP2_SHOCK_STATE] = curStoch.shocks[0][0][0][0][0][EF_K2];
	newState.current_states[WAGE_SHOCK_STATE] = curStoch.shocks[0][0][0][0][0][EF_W];
	newState.current_states[AGG_SHOCK_STATE] = curStoch.shocks[0][0][0][0][0][EF_A];
	newState.current_states[PHI_STATE] = curSt->current_states[PHI_STATE];
//	newState.current_states[PHI_STATE] = curStoch.shocks[0][0][0][0][0][EF_PHI];
	newState.current_indices[CAP1_SHOCK_STATE] = 0;
	newState.current_indices[CAP2_SHOCK_STATE] = 0;
	newState.current_indices[WAGE_SHOCK_STATE] = 0;
	newState.current_indices[AGG_SHOCK_STATE] = 0;
	newState.current_indices[PHI_STATE] = curSt->current_indices[PHI_STATE];
//	newState.current_indices[PHI_STATE] = 0;

	int phiState = newState.current_indices[PHI_STATE];
#if 0
	if (phiState == 1){
		std::cout << newState.current_states[PHI_STATE] << std::endl;
		std::cout << newState.prices[phiState][P1] << " " << newState.prices[phiState][P2] << std::endl;
		exit(-1);
	}
#endif
	double x = get_wage(newState, curStoch) + (newState.getTau() * k1) + (k2) + prod_fn(k1, k2, /*c1mgmt,*/ newState, curStoch);
	return MAX(MIN_BONDS, MIN_ASSETS - x) / (1 + curSt->getRecursiveVal(P_R));
}

bool vfiMaxUtil::constraintBinds() const{
	//calculate marginal utility of consuming all assets
	double x = getMaxBorrow(0, 0/*, 50*/);
	double consUtil = marginalConsUtil(curSt->current_states[ASTATE] + (-x));

	//find derivative of expected value function at a=0 
	double valueFn[2];
	double deriv = 0;
	int phiSt = curSt->current_indices[PHI_STATE];
	int aggSt = curSt->current_indices[AGG_SHOCK_STATE];
	int cap1St = curSt->current_indices[CAP1_SHOCK_STATE];
	int cap2St = curSt->current_indices[CAP2_SHOCK_STATE];
	int wageSt = curSt->current_indices[WAGE_SHOCK_STATE];
	int f = curSt->current_indices[AGG_ASSET_STATE];
	int ff = curSt->current_indices[AGG2_ASSET_STATE];

	for (int g = 0; g < PHI_STATES; g++) {
		for (int h = 0; h < AGG_SHOCK_SIZE; h++) {
			for (int i = 0; i < CAP_SHOCK_SIZE; i++) {
				for (int ii = 0; ii < CAP_SHOCK_SIZE; ii++) {
					for (int j = 0; j < WAGE_SHOCK_SIZE; j++) {
						for (int k = 0; k < 2; k++) {
							valueFn[k] = curFns.value_fn[f][ff][g][h][k][i][ii][j];
						}
						deriv += curStoch.transition[phiSt][aggSt][cap1St][cap2St][wageSt][g][h][i][ii][j] * (valueFn[1] - valueFn[0]) / (curStoch.assets[1] - curStoch.assets[0]);
					}
				}
			}
		}
	}
	//if beta*derivative is less than marginal consumption utility, then constrained 
	deriv = BETA*deriv;
	return (consUtil > deriv);
}

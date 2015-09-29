#include "eulerUtility.h"

eulerUtility::eulerUtility(const State& current, const StochProc& stoch, const EquilFns& fns) : vfiMaxUtil(current, stoch, fns)
{
	//get derivative of value function and save
	valueFnDeriv.resize(PHI_STATES);
	for (int h = 0; h < AGG_SHOCK_SIZE; h++) {
		valueFnDeriv[h].resize(AGG_SHOCK_SIZE);
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			valueFnDeriv[h][i].resize(ASSET_SIZE);
			for (int ii = 0; ii < ASSET_SIZE; ii++) {
				valueFnDeriv[h][i][ii].resize(CAP_SHOCK_SIZE);
				for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
					valueFnDeriv[h][i][ii][l].resize(CAP_SHOCK_SIZE);
					for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
						valueFnDeriv[h][i][ii][l][ll].resize(WAGE_SHOCK_SIZE);
						new VecDoub[WAGE_SHOCK_SIZE];
						for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
							switch (ii){
							case 0:
								valueFnDeriv[h][i][ii][l][ll][m] = (fns.value_fn[h][i][ii + 1][l][ll][m] - fns.value_fn[h][i][ii][l][ll][m]) /
									(stoch.assets[ii + 1] - stoch.assets[0]);
								break;
							case ASSET_SIZE - 1:
								valueFnDeriv[h][i][ii][l][ll][m] = (fns.value_fn[h][i][ii][l][ll][m] - fns.value_fn[h][i][ii - 1][l][ll][m]) /
									(stoch.assets[ii] - stoch.assets[ii - 1]);
								break;
							default:
								double fulldistance = stoch.assets[ii + 1] - stoch.assets[ii - 1];
								double lhdistance = stoch.assets[ii] - stoch.assets[ii - 1];
								double lhderiv = (fns.value_fn[h][i][ii][l][ll][m] - fns.value_fn[h][i][ii - 1][l][ll][m]) / lhdistance;
								double rhderiv = (fns.value_fn[h][i][ii + 1][l][ll][m] - fns.value_fn[h][i][ii][l][ll][m]) / (fulldistance - lhdistance);
								valueFnDeriv[h][i][ii][l][ll][m] = (1 - lhdistance / fulldistance)*lhderiv + lhdistance*rhderiv;
								break;
							}
						}
					}
				}
			}
		}
	}
}


eulerUtility::~eulerUtility()
{
}

double eulerUtility::operator() (const VecDoub state_prime) const
{

	double k1prime, k2prime;
	double bprime, c1mgmt,tmpbprime;
	double penalty = 0;
	double consumption;

	k1prime = MIN_CAPITAL + utilityFunctions::integer_power(state_prime[K1STATE], 2);
	k2prime = MIN_CAPITAL + utilityFunctions::integer_power(state_prime[K2STATE], 2);
	bprime = MIN_BONDS + utilityFunctions::integer_power(state_prime[BSTATE], 2);
	c1mgmt = MIN_MGMT + utilityFunctions::integer_power(state_prime[MGMT_C1_STATE], 2);

	c1mgmt = utilityFunctions::boundValue(c1mgmt,MIN_MGMT,MAX_MGMT);

	double x = getMaxBorrow(k1prime, k2prime, c1mgmt);
	if (!std::isfinite(bprime)){
		std::cout << "ERROR! eulerUtility() - b for next period is infinite"
			<< bprime
			<< std::endl;
	}
	tmpbprime = bprime;
	bprime = utilityFunctions::boundValue(bprime, x, (double)(2 * MAX_ASSETS));
	if (!std::isfinite(bprime)){
		std::cout << "ERROR! eulerUtility() - b for next period is infinite after boundValue"<<std::endl
			<< "K1:" << k1prime << std::endl
			<< "K2:" << k2prime << std::endl
			<< "C1Mg:" << c1mgmt << std::endl
			<< "Min:" << x << std::endl 
			<< "Orig: " << tmpbprime << std::endl 
			<< "New: " 	<< bprime << std::endl;
		exit(-1);
	}

	consumption = curSt->current_states[ASTATE] - (curSt->prices[curSt->current_indices[PHI_STATE]][P1] * k1prime)
		- (curSt->prices[curSt->current_indices[PHI_STATE]][P2] * k2prime) - bprime;

	if (consumption <= 0) {
		consumption = .001;
	}

	//find marginal utility of consumption
	double margUtilC = marginalConsUtil(consumption);

	//calculate expected marginal utility of a'
	State newState(*curSt);
	double expMargU=0;
	double *valueFn = new double[ASSET_SIZE];
	for (int g = 0; g < PHI_STATES; g++) {
		for (int h = 0; h < AGG_SHOCK_SIZE; h++) {
			for (int j = 0; j < WAGE_SHOCK_SIZE; j++) {
				for (int i = 0; i < CAP_SHOCK_SIZE; i++) {
					for (int ii = 0; ii < CAP_SHOCK_SIZE; ii++) {
						for (int k = 0; k < ASSET_SIZE; k++) {
							valueFn[k] = valueFnDeriv[g][h][k][i][ii][j];
						}

						newState.current_states[CAP1_SHOCK_STATE] =
							curStoch.shocks[g][h][i][ii][j][EF_K1];
						newState.current_states[CAP2_SHOCK_STATE] =
							curStoch.shocks[g][h][i][ii][j][EF_K2];
						newState.current_states[WAGE_SHOCK_STATE] =
							curStoch.shocks[g][h][i][ii][j][EF_W];
						newState.current_states[AGG_STATE] =
							curStoch.shocks[g][h][i][ii][j][EF_A];
						newState.current_states[PHI_STATE] =
							curStoch.shocks[g][h][i][ii][j][EF_PHI];
						newState.current_indices[CAP1_SHOCK_STATE] = i;
						newState.current_indices[CAP2_SHOCK_STATE] = ii;
						newState.current_indices[WAGE_SHOCK_STATE] = j;
						newState.current_indices[AGG_STATE] = h;
						newState.current_indices[PHI_STATE] = g;

						double aprime = get_wage(newState, curStoch) + (newState.prices[g][P1] * k1prime)
							+ (newState.prices[g][P2] * k2prime)
							+ (1 + newState.prices[g][P_R]) * bprime
							+ prod_fn(k1prime, k2prime, c1mgmt, newState, curStoch);

						if (aprime != aprime){
							std::cout << "ERROR! eulerUtility() - aprime for next period is NaN"
								<< std::endl
								<< "Wage: " << get_wage(newState, curStoch) << std::endl
								<< "k1: " << k1prime << std::endl
								<< "k2: " << k2prime << std::endl
								<< "B: " << bprime << std::endl
								<< "prod: " << prod_fn(k1prime, k2prime, c1mgmt, newState, curStoch)
								<< std::endl;
							exit(-1);
						}
						int phiSt = curSt->current_indices[PHI_STATE];
						int agSt = curSt->current_indices[AGG_STATE];
						int cap1St = curSt->current_indices[CAP1_SHOCK_STATE];
						int cap2St = curSt->current_indices[CAP2_SHOCK_STATE];
						int wgSt = curSt->current_indices[WAGE_SHOCK_STATE];
						expMargU += BETA * curStoch.transition[phiSt][agSt][cap1St][cap2St][wgSt][g][h][i][ii][j] * MAX(0, utilityFunctions::interpolate(curStoch.assets, valueFn, aprime));
					}
				}
			}
		}
	}
	delete valueFn;
	if (expMargU != expMargU) {
		std::cout << "ERROR! eulerUtility() - utility is " << expMargU
			<< ". (a:z1,z2)=(" << curSt->current_states[ASTATE] << ","
			<< curSt->current_states[CAP1_SHOCK_STATE] << ","
			<< curSt->current_states[CAP2_SHOCK_STATE] << ")"
			<< " k1prime=" << k1prime << " k2prime=" << k2prime
			<< " bprime=" << bprime << " consumption=" << consumption
			<< std::endl << std::flush;
		exit(-1);
	}

	//return absolute value of difference
	return abs(expMargU-margUtilC);
}

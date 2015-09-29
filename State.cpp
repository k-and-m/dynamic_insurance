#include "State.h"

State::State() {
	phi = VecDoub(PHI_STATES);
	prices = VecDoub(1);
	prices[0] = 1;

	for (int i = 0; i < PHI_STATES; i++){
		phi[i] = 1;
	}
	coefficients.resize(NUM_RECURSIVE_FNS);
	for (int i = 0; i < NUM_RECURSIVE_FNS; i++){
		coefficients[i].resize(PHI_STATES);
		for (int j = 0; j < PHI_STATES; j++){
			coefficients[i][j].resize(3);
			coefficients[i][j][0] = 0;
			coefficients[i][j][1] = 1;
			coefficients[i][j][2] = 0;
		}
	}
}

State::State(const VecDoub& phi1, const VecDoub& p_prices) : phi(phi1),prices(p_prices){
	if (phi1.size() != PHI_STATES){
		std::cerr << "State(VecDoub,VecDoub): phi has different dimension than number of states" << std::endl;
		exit(-1);
	}

	if (prices.size() != 1){
		std::cerr << "State(VecDoub,VecDoub): only expect 1 price" << std::endl;
		exit(-1);
	}

	current_states = VecDoub(NUM_STATE_VARS); //capital and bonds
	current_indices = VecInt(NUM_STATE_VARS);
	int counter = 0;
	coefficients.resize(NUM_RECURSIVE_FNS);
	for (int i = 0; i < NUM_RECURSIVE_FNS; i++){
		coefficients[i].resize(PHI_STATES);
		for (int j = 0; j < PHI_STATES; j++){
			coefficients[i][j].resize(3);
			coefficients[i][j][0] = 0;
			coefficients[i][j][1] = 1;
			coefficients[i][j][2] = 0;
		}
	}
}

State::State(const State& orig) : phi(orig.phi), prices(orig.prices){
	current_states = VecDoub(orig.current_states); //capital and bonds
	current_indices = VecInt(orig.current_indices);
	coefficients.resize(NUM_RECURSIVE_FNS);
	for (int i = 0; i < NUM_RECURSIVE_FNS; i++){
		coefficients[i].resize(PHI_STATES);
		for (int j = 0; j < PHI_STATES; j++){
			coefficients[i][j].resize(3);
			coefficients[i][j][0] = orig.coefficients[i][j][0];
			coefficients[i][j][1] = orig.coefficients[i][j][1];
			coefficients[i][j][2] = orig.coefficients[i][j][2];
		}
	}
}

State& State::operator=(const State& fnSource) {
	current_states = fnSource.current_states; //capital and bonds
	current_indices = fnSource.current_indices;
	phi = fnSource.phi;
	coefficients.resize(fnSource.coefficients.size());
	for (int i = 0; i < coefficients.size(); i++){
		coefficients[i].resize(PHI_STATES);
		for (int j = 0; j < PHI_STATES; j++){
			coefficients[i][j].resize(3);
			coefficients[i][j][0] = fnSource.coefficients[i][j][0];
			coefficients[i][j][1] = fnSource.coefficients[i][j][1];
			coefficients[i][j][2] = fnSource.coefficients[i][j][2];
		}
	}
	return *this;
}

double State::getRecursiveVal(int whichVal) const{
	int phiState = current_indices[PHI_STATE];
	if (whichVal == P_R) {
		return 0.03;
	}
	return exp(coefficients[whichVal][phiState][0] + coefficients[whichVal][phiState][1] * log(current_states[AGG_ASSET_STATE]) + coefficients[whichVal][phiState][2] * log(current_states[AGG2_ASSET_STATE]));
}

double State::getTau() const{
	return prices[0];
}

void State::defaultInitialState(StochProc& stoch2){

	 current_states[ASTATE] = stoch2.assets[ASSET_SIZE / 4];
	 current_states[CAP1_SHOCK_STATE] =
		stoch2.shocks[0][0][0][0][0][EF_K1];
	 current_states[CAP2_SHOCK_STATE] =
		stoch2.shocks[0][0][0][0][0][EF_K2];
	 current_states[WAGE_SHOCK_STATE] =
		stoch2.shocks[0][0][0][0][0][EF_W];
	 current_states[AGG_SHOCK_STATE] =
		stoch2.shocks[0][0][0][0][0][EF_A];
	 current_states[PHI_STATE] =
		stoch2.shocks[0][0][0][0][0][EF_PHI];
	 current_states[AGG_ASSET_STATE] = stoch2.aggAssets[AGG_ASSET_SIZE / 2];
	 current_states[AGG2_ASSET_STATE] = stoch2.aggAssets[AGG_ASSET_SIZE / 2];

	 current_indices[ASTATE] = ASSET_SIZE / 4;
	 current_indices[CAP1_SHOCK_STATE] = 0;
	 current_indices[CAP2_SHOCK_STATE] = 0;
	 current_indices[WAGE_SHOCK_STATE] = 0;
	 current_indices[AGG_SHOCK_STATE] = 0;
	 current_indices[PHI_STATE] = 0;
	 current_indices[AGG_ASSET_STATE] = AGG_ASSET_SIZE / 2;
	 current_indices[AGG2_ASSET_STATE] = AGG_ASSET_SIZE / 2;
}
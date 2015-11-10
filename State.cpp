#include "State.h"

State::State(){
	phi = VecDoub(PHI_STATES);
	for (int i = 0; i < PHI_STATES; i++) {
		phi[i] = 1;
	}

	prices = VecDoub(2);
	prices[0] = 1;
	prices[1] = 1;

	current_states = VecDoub(NUM_STATE_VARS); //capital and bonds
	current_indices = VecInt(NUM_STATE_VARS);
}

State::State(const VecDoub& phi1, const VecDoub& p_prices) : phi(phi1),prices(p_prices) {
	if (phi1.size() != PHI_STATES){
		std::cerr << "State(VecDoub,VecDoub): phi has different dimension than number of states" << std::endl;
		exit(-1);
	}

	if (prices.size() != 2){
		std::cerr << "State(VecDoub,VecDoub): expect 2 prices, received " << prices.size() << std::endl;
		exit(-1);
	}

	current_states = VecDoub(NUM_STATE_VARS); //capital and bonds
	current_indices = VecInt(NUM_STATE_VARS);
}

State::State(const State& orig) : phi(orig.phi), prices(orig.prices){
	current_states = VecDoub(orig.current_states); //capital and bonds
	current_indices = VecInt(orig.current_indices);
}

State& State::operator=(const State& fnSource) {
	current_states = fnSource.current_states; //capital and bonds
	current_indices = fnSource.current_indices;
	phi = fnSource.phi;
	prices = fnSource.prices;
	return *this;
}

double State::getNextR() const{
	return prices[1];
}

double State::getTau() const{
	return prices[0];
}

void State::defaultInitialState(StochProc& stoch2){

	 current_states[ASTATE] = stoch2.assets[ASSET_SIZE / 2];
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
	 current_indices[ASTATE] = ASSET_SIZE / 2;
	 current_indices[CAP1_SHOCK_STATE] = 0;
	 current_indices[CAP2_SHOCK_STATE] = 0;
	 current_indices[WAGE_SHOCK_STATE] = 0;
	 current_indices[AGG_SHOCK_STATE] = 0;
	 current_indices[PHI_STATE] = 0;
}
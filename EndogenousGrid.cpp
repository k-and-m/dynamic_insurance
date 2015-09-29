#if 0
#include "EndogenousGrid.h"


EndogenousGrid::EndogenousGrid()
{
}


EndogenousGrid::~EndogenousGrid()
{
}

void EndogenousGrid::solve()
{
	extern double get_wage(const State& current, const StochProc& stoch);

	double G_y[WAGE_SHOCK_SIZE][ASSET_SIZE];
	double valueSlope[WAGE_SHOCK_SIZE][ASSET_SIZE];
	double xDim[ASSET_SIZE];
	double fx[ASSET_SIZE];
	double temp_v_np1[WAGE_SHOCK_SIZE][ASSET_SIZE];

	//preliminary setup
	currentState.current_indices[AGG_STATE] = 0;
	for (int i = 0; i < WAGE_SHOCK_SIZE; i++){
		currentState.current_indices[WAGE_SHOCK_STATE] = i;
		for (int j = 0; j < ASSET_SIZE; j++){
			G_y(i, j) = a_v(i)*(1 + currentState.r) + get_wage(currentState, *s_proc);
		}
	}

	for (int i = 0; i < WAGE_SHOCK_SIZE; i++){
		currentState.current_indices[WAGE_SHOCK_STATE] = i;
		valueSlope(i, 0) = (value_v(i, 2) - value_v(i, 1)) / (a_v(1) - a_v(0));
		valueSlope(i, ASSET_SIZE - 1) = (value_v(i, ASSET_SIZE - 1) - value_v(i, ASSET_SIZE - 2))
			/ (a_v(ASSET_SIZE - 1) - a_v(ASSET_SIZE - 2));
		for (int j = 0; j < ASSET_SIZE; j++){
			valueSlope(i, j) = (value_v(i, j) - value_v(i, j - 1)) / (a_v(j) - a_v(j - 1));
			maxC(i, j) = (1 + currentState.r)*a_v(j) + get_wage(currentState, *s_proc) + (-MIN_ASSETS);
			c_star(i, j) = temp_vk(i, j)**(-1.0 / RRA);
			if (c_star(i, j) > maxC(i, j)){
				c_star(i, j) = maxC(i, j);
			}
			y_end(i, j) = c_star(i, j) + a_v(j);
			temp_v_n(i, j) = (c_star(i, j)**(1 - RRA)) / (1.0 - RRA) + value_v(i, j);
		}
	}

	for (int i = 0; i < WAGE_SHOCK_SIZE; i++){
		xs = ;
		fxs = ;
		for (int j = 0; j < ASSET_SIZE; j++){
			temp_v_np1(i, j) = MIN(sub_myinterp1(xs, fxs, G_y(i, j)), 0);
		}
	}

	for (int i = 1; i < ASSET_SIZE; i++){
		for (int j = 0; j < WAGE_SHOCK_SIZE; j++){
				value_v_t(j, i) = beta_p*dot_product(transition(j, :), temp_v_np1(:, i))
				temp_k(j, i) = (y_end(j, i) - exp(stateSpace(j))) / (1.0_dp + r1)
		}
	}

}
#endif
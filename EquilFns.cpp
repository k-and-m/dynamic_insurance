#include "EquilFns.h"

EquilFns::EquilFns() {
	value_fn.resize(AGG_ASSET_SIZE);
	consumption.resize(AGG_ASSET_SIZE);
	policy_fn.resize(AGG_ASSET_SIZE);
	for (int f = 0; f < AGG_ASSET_SIZE; f++){
		value_fn[f].resize(AGG_ASSET_SIZE);
		consumption[f].resize(AGG_ASSET_SIZE);
		policy_fn[f].resize(AGG_ASSET_SIZE);
		for (int g = 0; g < AGG_ASSET_SIZE; g++){
			value_fn[f][g].resize(PHI_STATES);
			consumption[f][g].resize(PHI_STATES);
			policy_fn[f][g].resize(PHI_STATES);
			for (int h = 0; h < PHI_STATES; h++) {
				value_fn[f][g][h].resize(AGG_SHOCK_SIZE);
				consumption[f][g][h].resize(AGG_SHOCK_SIZE);
				policy_fn[f][g][h].resize(AGG_SHOCK_SIZE);
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					value_fn[f][g][h][i].resize(ASSET_SIZE);
					consumption[f][g][h][i].resize(ASSET_SIZE);
					policy_fn[f][g][h][i].resize(ASSET_SIZE);
					for (int ii = 0; ii < ASSET_SIZE; ii++) {
						value_fn[f][g][h][i][ii].resize(CAP_SHOCK_SIZE);
						consumption[f][g][h][i][ii].resize(CAP_SHOCK_SIZE);
						policy_fn[f][g][h][i][ii].resize(CAP_SHOCK_SIZE);
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							value_fn[f][g][h][i][ii][l].resize(CAP_SHOCK_SIZE);
							consumption[f][g][h][i][ii][l].resize(CAP_SHOCK_SIZE);
							policy_fn[f][g][h][i][ii][l].resize(CAP_SHOCK_SIZE);
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								value_fn[f][g][h][i][ii][l][ll].resize(WAGE_SHOCK_SIZE);
								consumption[f][g][h][i][ii][l][ll].resize(WAGE_SHOCK_SIZE);
								policy_fn[f][g][h][i][ii][l][ll].resize(WAGE_SHOCK_SIZE);
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									policy_fn[f][g][h][i][ii][l][ll][m].resize(NUM_CHOICE_VARS);
								}
							}
						}
					}
				}
			}
		}
	}
}

EquilFns::~EquilFns() {
}

EquilFns& EquilFns::operator=(const EquilFns& fnSource) {
	if (this != &fnSource) {
		for (int f = 0; f < AGG_ASSET_SIZE; f++){
			for (int g = 0; g < AGG_ASSET_SIZE; g++){
				for (int h = 0; h < PHI_STATES; h++){
					for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
						for (int ii = 0; ii < ASSET_SIZE; ii++) {
							for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
								for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
									for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
										value_fn[f][g][h][i][ii][l][ll][m] =
											fnSource.value_fn[f][g][h][i][ii][l][ll][m];
										consumption[f][g][h][i][ii][l][ll][m] =
											fnSource.consumption[f][g][h][i][ii][l][ll][m];
										for (int n = 0; n < NUM_CHOICE_VARS; n++) {
											policy_fn[f][g][h][i][ii][l][ll][m][n] =
												fnSource.policy_fn[f][g][h][i][ii][l][ll][m][n];
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
	return *this;
}

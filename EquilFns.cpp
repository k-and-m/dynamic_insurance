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
									value_fn[f][g][h][i][ii][l][ll][m] = 0;
									consumption[f][g][h][i][ii][l][ll][m] = 0;
									policy_fn[f][g][h][i][ii][l][ll][m].resize(NUM_CHOICE_VARS);
									for (int n = 0; n < NUM_CHOICE_VARS; n++) {
										policy_fn[f][g][h][i][ii][l][ll][m][n] = 0;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	reordered_value_fn.resize(PHI_STATES);
	for (int f = 0; f < PHI_STATES; f++) {
		reordered_value_fn[f].resize(AGG_SHOCK_SIZE);
		for (int g = 0; g < AGG_SHOCK_SIZE; g++) {
			reordered_value_fn[f][g].resize(CAP_SHOCK_SIZE);
			for (int h = 0; h < CAP_SHOCK_SIZE; h++) {
				reordered_value_fn[f][g][h].resize(CAP_SHOCK_SIZE);
				for (int i = 0; i < CAP_SHOCK_SIZE; i++) {
					reordered_value_fn[f][g][h][i].resize(WAGE_SHOCK_SIZE);
					for (int l = 0; l < WAGE_SHOCK_SIZE; l++) {
						reordered_value_fn[f][g][h][i][l].resize(AGG_ASSET_SIZE);
						for (int ll = 0; ll < AGG_ASSET_SIZE; ll++) {
							reordered_value_fn[f][g][h][i][l][ll].resize(AGG_ASSET_SIZE);
							for (int m = 0; m < AGG_ASSET_SIZE; m++) {
								reordered_value_fn[f][g][h][i][l][ll][m].resize(ASSET_SIZE);
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

		for (int f = 0; f < PHI_STATES; f++) {
			for (int g = 0; g < AGG_SHOCK_SIZE; g++) {
				for (int h = 0; h < CAP_SHOCK_SIZE; h++) {
					for (int i = 0; i < CAP_SHOCK_SIZE; i++) {
						for (int l = 0; l < WAGE_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < AGG_ASSET_SIZE; ll++) {
								for (int m = 0; m < AGG_ASSET_SIZE; m++) {
									for (int n = 0; n < ASSET_SIZE; n++) {
										reordered_value_fn[f][g][h][i][l][ll][m][n] =
											fnSource.reordered_value_fn[f][g][h][i][l][ll][m][n];
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

void EquilFns::setValueFn(VecInt_I param, const double val) {
	if (param.size() != 8) {
		std::cerr << "EquilFns.setValueFn(): Expect 8 indicies. Only got " << param.size() << std::endl;
		exit(-1);
	}

	value_fn[param[0]][param[1]][param[2]][param[3]][param[4]][param[5]][param[6]][param[7]] = val;
	reordered_value_fn[param[2]][param[3]][param[5]][param[6]][param[7]][param[0]][param[1]][param[4]] = val;

}

double EquilFns::getValueFn(VecInt_I param) const{
	if (param.size() != 8) {
		std::cerr << "EquilFns.getValueFn(): Expect 8 indicies. Only got " << param.size() << std::endl;
		exit(-1);
	}

	return value_fn[param[0]][param[1]][param[2]][param[3]][param[4]][param[5]][param[6]][param[7]];
}

int EquilFns::policyToArray(double *toAlloc) const {
	int counter = 0;
	for (int f = 0; f < AGG_ASSET_SIZE; f++) {
		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int h = 0; h < PHI_STATES; h++) {
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int ii = 0; ii < ASSET_SIZE; ii++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									for (int n = 0; n < NUM_CHOICE_VARS; n++) {
										toAlloc[counter] = policy_fn[f][g][h][i][ii][l][ll][m][n];
										counter++;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return counter;
}

int EquilFns::valueToArray(double *toAlloc) const {
	int counter = 0;
	for (int f = 0; f < AGG_ASSET_SIZE; f++) {
		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int h = 0; h < PHI_STATES; h++) {
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int ii = 0; ii < ASSET_SIZE; ii++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
										toAlloc[counter] = value_fn[f][g][h][i][ii][l][ll][m];
										counter++;
								}
							}
						}
					}
				}
			}
		}
	}
	return counter;
}

void EquilFns::setPolicyFromArray(double *values) {
	int counter = 0;
	for (int f = 0; f < AGG_ASSET_SIZE; f++) {
		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int h = 0; h < PHI_STATES; h++) {
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int ii = 0; ii < ASSET_SIZE; ii++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									for (int n = 0; n < NUM_CHOICE_VARS; n++) {
										policy_fn[f][g][h][i][ii][l][ll][m][n] = values[counter];
										counter++;
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

void EquilFns::setValueFromArray(double *values) {
	int counter = 0;
	for (int f = 0; f < AGG_ASSET_SIZE; f++) {
		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int h = 0; h < PHI_STATES; h++) {
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int ii = 0; ii < ASSET_SIZE; ii++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
										VecInt param(8);
										param[0] = f;
										param[1] = g;
										param[2] = h;
										param[3] = i;
										param[4] = ii;
										param[5] = l;
										param[6] = ll;
										param[7] = m;

										setValueFn(param, values[counter]);
										counter++;
								}
							}
						}
					}
				}
			}
		}
	}
}


vector<vector<VecDoub>>::const_iterator EquilFns::getReorderedValueFnVector(VecInt_I param) const{
	if (param.size() != 5) {
		std::cerr << "EquilFns.getReorderedValueFnVector(): Expect 7 indicies. Only got " << param.size() << std::endl;
		exit(-1);
	}

	return reordered_value_fn[param[0]][param[1]][param[2]][param[3]][param[4]].cbegin();
}

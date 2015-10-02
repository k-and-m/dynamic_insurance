#pragma once

#include "StochProc.h"

StochProc::StochProc(const StochProc& init){
	if (PHI_STATES == 1){
		phiTrans[0][0] = init.phiTrans[0][0];
		phiTransCDF[0][0] = init.phiTransCDF[0][0];
	}
	else if (PHI_STATES == 2){
		phiTrans[0][0] = init.phiTrans[0][0];
		phiTrans[0][1] = init.phiTrans[0][1];
		phiTrans[1][0] = init.phiTrans[1][0];
		phiTrans[1][1] = init.phiTrans[1][1];
		phiTransCDF[0][0] = init.phiTransCDF[0][0];
		phiTransCDF[0][1] = init.phiTransCDF[0][1];
		phiTransCDF[1][0] = init.phiTransCDF[1][0];
		phiTransCDF[1][1] = init.phiTransCDF[1][1];
	}
	else{
		std::cout << "expect only two phi states" << std::endl;
		exit(-1);
	}

	assets.resize(ASSET_SIZE);
	aggAssets.resize(AGG_ASSET_SIZE);
	transition.resize(PHI_STATES);
	cdf.resize(PHI_STATES);
	shocks.resize(PHI_STATES);
	for (int h = 0; h < PHI_STATES; h++){
		transition[h].resize(AGG_SHOCK_SIZE);
		cdf[h].resize(AGG_SHOCK_SIZE);
		shocks[h].resize(AGG_SHOCK_SIZE);
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			transition[h][i].resize(CAP_SHOCK_SIZE);
			cdf[h][i].resize(CAP_SHOCK_SIZE);
			shocks[h][i].resize(CAP_SHOCK_SIZE);
			for (int ii = 0; ii < CAP_SHOCK_SIZE; ii++) {
				transition[h][i][ii].resize(CAP_SHOCK_SIZE);
				cdf[h][i][ii].resize(CAP_SHOCK_SIZE);
				shocks[h][i][ii].resize(CAP_SHOCK_SIZE);
				for (int j = 0; j < CAP_SHOCK_SIZE; j++) {
					transition[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					cdf[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					shocks[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					for (int k = 0; k < WAGE_SHOCK_SIZE; k++) {
						shocks[h][i][ii][j][k].resize(NUM_SHOCK_VARS);
					}
				}
			}
		}
	}

	for (int h = 0; h < PHI_STATES; h++) {
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			for (int j = 0; j < CAP_SHOCK_SIZE; j++) {
				for (int jj = 0; jj < CAP_SHOCK_SIZE; jj++) {
					for (int l = 0; l < WAGE_SHOCK_SIZE; l++) {
						shocks[h][i][j][jj][l][EF_K1] = init.shocks[h][i][j][jj][l][EF_K1];
						shocks[h][i][j][jj][l][EF_K2] = init.shocks[h][i][j][jj][l][EF_K2];
						shocks[h][i][j][jj][l][EF_W] = init.shocks[h][i][j][jj][l][EF_W];
						shocks[h][i][j][jj][l][EF_A] = init.shocks[h][i][j][jj][l][EF_A];
						shocks[h][i][j][jj][l][EF_PHI] = init.shocks[h][i][j][jj][l][EF_PHI];

						for (int mm1 = 0; mm1 < PHI_STATES; mm1++) {
							for (int m = 0; m < AGG_SHOCK_SIZE; m++) {
								for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
									for (int nn = 0; nn < CAP_SHOCK_SIZE; nn++) {
										for (int o = 0; o < WAGE_SHOCK_SIZE; o++) {
											transition[h][i][j][jj][l][mm1][m][n][nn][o] = init.transition[h][i][j][jj][l][mm1][m][n][nn][o];
											cdf[h][i][j][jj][l][mm1][m][n][nn][o] = init.cdf[h][i][j][jj][l][mm1][m][n][nn][o];
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

	for (int i = 0; i < ASSET_SIZE; i++) {
		assets[i] = init.assets[i];
	}

	for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
		aggAssets[i] = init.aggAssets[i];
	}

}

StochProc::StochProc(VecDoub phis) {

	if (phis.size() != PHI_STATES){
		std::cerr << "StochProc(VecDoub): phis parameter is not of size PHI_STATES" << std::endl;
		exit(-1);
	}

	if (PHI_STATES == 1){
		phiTrans[0][0] = 1;
		phiTransCDF[0][0] = 1;
	}
	else if (PHI_STATES == 2){
#if 1
#if 1
		phiTrans[0][0] = 0.5;
		phiTrans[0][1] = 0.5;
		phiTrans[1][0] = 0.125;
		phiTrans[1][1] = 0.875;
		phiTransCDF[0][0] = 0.5;
		phiTransCDF[0][1] = 1;
		phiTransCDF[1][0] = 0.125;
		phiTransCDF[1][1] = 1;
#else
		phiTrans[0][0] = 0.99999;
		phiTrans[0][1] = 0.00001;
		phiTrans[1][0] = 0.00001;
		phiTrans[1][1] = 0.99999;
		phiTransCDF[0][0] = 0.99999;
		phiTransCDF[0][1] = 1;
		phiTransCDF[1][0] = 0.00001;
		phiTransCDF[1][1] = 1;
#endif
#else
		phiTrans[0][0] = 0.5;
		phiTrans[0][1] = 0.5;
		phiTrans[1][0] = 0.5;
		phiTrans[1][1] = 0.5;
		phiTransCDF[0][0] = 0.5;
		phiTransCDF[0][1] = 1;
		phiTransCDF[1][0] = 0.5;
		phiTransCDF[1][1] = 1;

#endif
	}
	else{
		std::cout << "expect only two phi states" << std::endl;
		exit(-1);
	}

	assets.resize(ASSET_SIZE);
	aggAssets.resize(AGG_ASSET_SIZE);
	transition.resize(PHI_STATES);
	cdf.resize(PHI_STATES);
	shocks.resize(PHI_STATES);
	for (int h = 0; h < PHI_STATES; h++){
		transition[h].resize(AGG_SHOCK_SIZE);
		cdf[h].resize(AGG_SHOCK_SIZE);
		shocks[h].resize(AGG_SHOCK_SIZE);
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			transition[h][i].resize(CAP_SHOCK_SIZE);
			cdf[h][i].resize(CAP_SHOCK_SIZE);
			shocks[h][i].resize(CAP_SHOCK_SIZE);
			for (int ii = 0; ii < CAP_SHOCK_SIZE; ii++) {
				transition[h][i][ii].resize(CAP_SHOCK_SIZE);
				cdf[h][i][ii].resize(CAP_SHOCK_SIZE);
				shocks[h][i][ii].resize(CAP_SHOCK_SIZE);
				for (int j = 0; j < CAP_SHOCK_SIZE; j++) {
					transition[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					cdf[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					shocks[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					for (int k = 0; k < WAGE_SHOCK_SIZE; k++) {
						shocks[h][i][ii][j][k].resize(NUM_SHOCK_VARS);
					}
				}
			}
		}
	}

	for (int i = 0; i < ASSET_SIZE; i++) {
		assets[i] = 1.0 / (ASSET_SIZE - 1) * i;
		assets[i] = pow(assets[i], CURVATURE);
		assets[i] = assets[i] * (MAX_ASSETS - MIN_ASSETS) + MIN_ASSETS;
	}

	for (int i = 0; i < AGG_ASSET_SIZE; i++) {
		if (AGG_ASSET_SIZE == 1){
			aggAssets[0] = (MAX_AGG_ASSETS + MIN_AGG_ASSETS) / 2;
		}else{
			aggAssets[i] = 1.0 / (MAX(1,AGG_ASSET_SIZE - 1)) * i;
			aggAssets[i] = pow(aggAssets[i], AGG_CURVATURE);
			aggAssets[i] = aggAssets[i] * (MAX_AGG_ASSETS - MIN_AGG_ASSETS) + MIN_AGG_ASSETS;
		}
	}

	for (int h = 0; h < PHI_STATES; h++) {
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			for (int j = 0; j < CAP_SHOCK_SIZE; j++) {
				for (int jj = 0; jj < CAP_SHOCK_SIZE; jj++) {
					for (int l = 0; l < WAGE_SHOCK_SIZE; l++) {
						double shocks1[CAP_SHOCK_SIZE];
						double shocks2[CAP_SHOCK_SIZE];
						double wageShocks[WAGE_SHOCK_SIZE];
						double transition1[CAP_SHOCK_SIZE][CAP_SHOCK_SIZE];
						double transition2[CAP_SHOCK_SIZE][CAP_SHOCK_SIZE];
						double wageTrans[WAGE_SHOCK_SIZE][WAGE_SHOCK_SIZE];
						if (CAP_SHOCK_SIZE == 1) {
							shocks1[j] = 0;
							shocks2[jj] = 0;
							transition1[j][jj] = 1;
							transition2[j][jj] = 1;
						}
						else {
#if 0
							extern void tauchen(double mean, double rho,
								double sigma, double m,
								double(&z)[CAP_SHOCK_SIZE],
								double(&zprob)[CAP_SHOCK_SIZE][CAP_SHOCK_SIZE]);
							tauchen(CAP1_SHOCK_MEAN, CAP1_SHOCK_PERSISTENCE,
								CAP1_SHOCK_STD, CAP1_SHOCK_SPREAD, shocks1,
								transition1);
							tauchen(CAP2_SHOCK_MEAN, CAP2_SHOCK_PERSISTENCE,
								CAP2_SHOCK_STD, CAP2_SHOCK_SPREAD, shocks2,
								transition2);
#else
							shocks1[0] = -0.225;
							shocks1[1] = 0.525;
							transition1[0][0] = 0.5;
							transition1[0][1] = 0.5;
							transition1[1][0] = 0.5;
							transition1[1][1] = 0.5;
							shocks2[0] = -0.225;
							shocks2[1] = 0.525;
							transition2[0][0] = 0.5;
							transition2[0][1] = 0.5;
							transition2[1][0] = 0.5;
							transition2[1][1] = 0.5;
#endif
						}
						if (WAGE_SHOCK_SIZE == 1) {
							wageShocks[l] = 0;
							wageTrans[l][l] = 1;
						}
						else {
#if 0
							extern void tauchen2(double mean, double rho,
								double sigma, double m,
								double(&z)[WAGE_SHOCK_SIZE],
								double(&zprob)[WAGE_SHOCK_SIZE][WAGE_SHOCK_SIZE]);
							tauchen2(WAGE_SHOCK_MEAN, WAGE_SHOCK_PERSISTENCE,
								WAGE_SHOCK_STD, WAGE_SHOCK_SPREAD,
								wageShocks, wageTrans);
#else
							if (WAGE_SHOCK_SIZE != 2){
								std::cout << "expect only two wage shocks" << std::endl;
								exit(-1);
							}
							wageShocks[0] = 0.34;
							wageShocks[1] = 1.36;
							wageTrans[0][0] = 0.95;
							wageTrans[0][1] = 0.05;
							wageTrans[1][0] = 0.05;
							wageTrans[1][1] = 0.95;

#endif
						}

						shocks[h][i][j][jj][l][EF_K1] = shocks1[j];
						shocks[h][i][j][jj][l][EF_K2] = shocks2[jj];
						shocks[h][i][j][jj][l][EF_W] = wageShocks[l];
						shocks[h][i][j][jj][l][EF_A] = 1;
						shocks[h][i][j][jj][l][EF_PHI] = phis[h];

						for (int mm1 = 0; mm1 < PHI_STATES; mm1++) {
							for (int m = 0; m < AGG_SHOCK_SIZE; m++) {
								for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
									for (int nn = 0; nn < CAP_SHOCK_SIZE; nn++) {
										for (int o = 0; o < WAGE_SHOCK_SIZE; o++) {
											transition[h][i][j][jj][l][mm1][m][n][nn][o] =
												transition1[j][n]
												* transition2[jj][nn]
												* wageTrans[l][o]
												* phiTrans[h][mm1];
										}
 									}
								}
							}
						}

						for (int mm1 = 0; mm1 < PHI_STATES; mm1++) {
							for (int m = 0; m < AGG_SHOCK_SIZE; m++) {
								double prevValue = 0;
								for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
									for (int nn = 0; nn < CAP_SHOCK_SIZE; nn++) {
										for (int o = 0; o < WAGE_SHOCK_SIZE; o++) {
											cdf[h][i][j][jj][l][mm1][m][n][nn][o] =
												prevValue +
												transition[h][i][j][jj][l][mm1][m][n][nn][o];
											prevValue = cdf[h][i][j][jj][l][mm1][m][n][nn][o];
										}
									}
								}
								for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
									for (int nn = 0; nn < CAP_SHOCK_SIZE; nn++) {
										for (int o = 0; o < WAGE_SHOCK_SIZE; o++) {
											cdf[h][i][j][jj][l][mm1][m][n][nn][o] = cdf[h][i][j][jj][l][mm1][m][n][nn][o] /
												cdf[h][i][j][jj][l][mm1][m][CAP_SHOCK_SIZE - 1][CAP_SHOCK_SIZE - 1][WAGE_SHOCK_SIZE - 1];
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
}

StochProc::~StochProc() {
}

int StochProc::getCondNewPhi(const int curSt, const double randNum) const{
	for (int i = 0; i < PHI_STATES; i++){
		if (randNum < phiTransCDF[curSt][i]){
			return i;
		}
	}
	std::cerr << "StochProc:getCondNewPhi - Should not reach here. rand="<<randNum << std::endl;
	std::cerr << "Transition CDF"<<std::endl;
	for (int i = 0; i < PHI_STATES; i++){
		std::cout << phiTransCDF[curSt][i] << " ";
	}
	exit(-1);
	return PHI_STATES - 1;
}

VecInt StochProc::getCondNewState(const VecInt& currentState, const int newAggState, const int newPhi, const double randNum) const{
	int h = currentState[EF_PHI];
	int i = currentState[EF_A];
	int j = currentState[EF_K1];
	int jj = currentState[EF_K2];
	int l = currentState[EF_W];

	VecInt newState = VecInt(NUM_SHOCK_VARS);

	for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
		for (int nn = 0; nn < CAP_SHOCK_SIZE; nn++) {
			for (int o = 0; o < WAGE_SHOCK_SIZE; o++) {
				if (randNum < cdf[h][i][j][jj][l][newPhi][newAggState][n][nn][o]){
					newState[EF_A] = newAggState;
					newState[EF_PHI] = newPhi;
					newState[EF_K1] = n;
					newState[EF_K2] = nn;
					newState[EF_W] = o;
					return newState;
				}
			}
		}
	}

	for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
		for (int nn = 0; nn < CAP_SHOCK_SIZE; nn++) {
			for (int o = 0; o < WAGE_SHOCK_SIZE; o++) {
				if (randNum < cdf[newPhi][i][j][jj][l][newPhi][newAggState][n][nn][o]){
					newState[EF_A] = newAggState;
					newState[EF_PHI] = newPhi;
					newState[EF_K1] = n;
					newState[EF_K2] = nn;
					newState[EF_W] = o;
					return newState;
				}
			}
		}
	}

	std::cerr << "Could not transition to new phi and new state." << std::endl;
	exit(-1);
	return newState;
}

StochProc& StochProc::operator=(const StochProc& init) {
	if (PHI_STATES == 1){
		phiTrans[0][0] = init.phiTrans[0][0];
		phiTransCDF[0][0] = init.phiTransCDF[0][0];
	}
	else if (PHI_STATES == 2){
		phiTrans[0][0] = init.phiTrans[0][0];
		phiTrans[0][1] = init.phiTrans[0][1];
		phiTrans[1][0] = init.phiTrans[1][0];
		phiTrans[1][1] = init.phiTrans[1][1];
		phiTransCDF[0][0] = init.phiTransCDF[0][0];
		phiTransCDF[0][1] = init.phiTransCDF[0][1];
		phiTransCDF[1][0] = init.phiTransCDF[1][0];
		phiTransCDF[1][1] = init.phiTransCDF[1][1];
	}
	else{
		std::cout << "expect only two phi states" << std::endl;
		exit(-1);
	}

	assets.resize(ASSET_SIZE);
	aggAssets.resize(AGG_ASSET_SIZE);
	transition.resize(PHI_STATES);
	cdf.resize(PHI_STATES);
	shocks.resize(PHI_STATES);
	for (int h = 0; h < PHI_STATES; h++){
		transition[h].resize(AGG_SHOCK_SIZE);
		cdf[h].resize(AGG_SHOCK_SIZE);
		shocks[h].resize(AGG_SHOCK_SIZE);
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			transition[h][i].resize(CAP_SHOCK_SIZE);
			cdf[h][i].resize(CAP_SHOCK_SIZE);
			shocks[h][i].resize(CAP_SHOCK_SIZE);
			for (int ii = 0; ii < CAP_SHOCK_SIZE; ii++) {
				transition[h][i][ii].resize(CAP_SHOCK_SIZE);
				cdf[h][i][ii].resize(CAP_SHOCK_SIZE);
				shocks[h][i][ii].resize(CAP_SHOCK_SIZE);
				for (int j = 0; j < CAP_SHOCK_SIZE; j++) {
					transition[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					cdf[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					shocks[h][i][ii][j].resize(WAGE_SHOCK_SIZE);
					for (int k = 0; k < WAGE_SHOCK_SIZE; k++) {
						shocks[h][i][ii][j][k].resize(NUM_SHOCK_VARS);
					}
				}
			}
		}
	}

	for (int h = 0; h < PHI_STATES; h++) {
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			for (int j = 0; j < CAP_SHOCK_SIZE; j++) {
				for (int jj = 0; jj < CAP_SHOCK_SIZE; jj++) {
					for (int l = 0; l < WAGE_SHOCK_SIZE; l++) {
						shocks[h][i][j][jj][l][EF_K1] = init.shocks[h][i][j][jj][l][EF_K1];
						shocks[h][i][j][jj][l][EF_K2] = init.shocks[h][i][j][jj][l][EF_K2];
						shocks[h][i][j][jj][l][EF_W] = init.shocks[h][i][j][jj][l][EF_W];
						shocks[h][i][j][jj][l][EF_A] = init.shocks[h][i][j][jj][l][EF_A];
						shocks[h][i][j][jj][l][EF_PHI] = init.shocks[h][i][j][jj][l][EF_PHI];

						for (int mm1 = 0; mm1 < PHI_STATES; mm1++) {
							for (int m = 0; m < AGG_SHOCK_SIZE; m++) {
								for (int n = 0; n < CAP_SHOCK_SIZE; n++) {
									for (int nn = 0; nn < CAP_SHOCK_SIZE; nn++) {
										for (int o = 0; o < WAGE_SHOCK_SIZE; o++) {
											transition[h][i][j][jj][l][mm1][m][n][nn][o] = init.transition[h][i][j][jj][l][mm1][m][n][nn][o];
											cdf[h][i][j][jj][l][mm1][m][n][nn][o] = init.cdf[h][i][j][jj][l][mm1][m][n][nn][o];
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

	for (int i = 0; i < ASSET_SIZE; i++) {
		assets[i] = init.assets[i];
	}
	for (int i = 0; i < AGG_ASSET_SIZE; i++) {
		aggAssets[i] = init.aggAssets[i];
	}
	return *this;
}

/*
 * akmodel.cc
 *
 *  Created on: Feb 11, 2015
 *      Author: robmrk
 */

#define POLICY_ITERATION 1
#define AMOEBA 1
#define SIMULATED_ANNEALING 0
#define POLICY_CONVERGENCE 0

#include "BaseHeader.h"
#include "utilityFunctions.h"
#include "matrixIO.h"
#include "matrix.h"
#include "amoeba.h"
#include "WorldEconomy.h"
#include "vfiMaxUtil.h"
#include <Eigen/Dense>

#if SIMULATED_ANNEALING
#include "simAnneal.hpp"
#endif

using namespace Eigen;

vector<MatrixXf> simulate(int numC, const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target);
double getNextParameters(const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target);
void printResults(const EquilFns& e, const VecDoub& phis, COUNTRYID whichCountry);
void initialize(EquilFns& fns, const VecDoub& phis);
void solve(const VecDoub& phis, const VecDoub& prices, const EquilFns& orig, EquilFns& final);
int policy_iteration(double difference, State& initial, const StochProc& stoch,
	EquilFns& lastFns, EquilFns& currentFns);
double maxValueDistance(const EquilFns &a, const EquilFns &b, VecInt& location);
double maxPolicyDistance(const EquilFns &a, const EquilFns &b,
	VecInt& location);
double solveProblem(const VecDoub& phis, const VecDoub& prices, double c1prop, int seqNo);

#if 1
int main(int argc, char *argv[]) {
	using namespace std;

	if (argc != 7){
		std::cout << "Incorrect number of arguments. Require: c1phiL c1phiH c2phiL c2phiH tau seqNo" << endl;
		//				exit(-1);
	}

#if 0
	VecDoub start = VecDoub(3);
	start[0] = atof(argv[3]);
	start[1] = atof(argv[4]);
	start[2] = atof(argv[5]);
	const double targ = 0;
	SimulatedAnnealingWorld worldSolve = SimulatedAnnealingWorld(start, targ, 100, 0.9, solveProblem, 0.01, atof(argv[1]), atof(argv[2]), 0.3);
	VecDoub *soln = worldSolve.solve();
	cout << (*soln)[0] << ":" << (*soln)[1] << ":" << (*soln)[2] << endl;
#endif

#if 0
	VecDoub phi = VecDoub(4);
	VecDoub prices = VecDoub(1);
	phi[0] = atof(argv[1]);
	phi[1] = atof(argv[2]);
	phi[2] = atof(argv[3]);
	phi[3] = atof(argv[4]);

	prices[0] = atof(argv[5]);
	double dis = solveProblem(phi, prices, 0.3, atoi(argv[argc-1]));
#else 
	VecDoub phi = VecDoub(4);
	VecDoub prices = VecDoub(1);
	phi[0] = 0.4;
	phi[1] = 0.8;
	phi[2] = 0;
	phi[3] = 0;

	prices[0] = 1.01;
	double dis = solveProblem(phi, prices, 0.3, 0);
#endif
	std::cout << "soln dist: " << dis << endl;

	return 0;
}
#endif

#if 0
extern "C"{
	double runonce_(double *c1p, double *c2p, double *c1Price, double *c2Price, double *intRate,
		double *c1prop)
	{
		return solveProblem(*c1p, *c2p, *c1Price, *c2Price, *intRate, *c1prop);
	}
}
#endif

double solveProblem(const VecDoub& phis, const VecDoub& prices, double c1prop, int seqNo)
{
	if (phis.size() != 2 * PHI_STATES){
		std::cerr << "solveProblem: incorrect number of phis " << phis.size() << endl;
		exit(-1);
	}
	if (prices.size() != 1){
		std::cerr << "solveProblem: incorrect number of prices " << prices.size() << endl;
		exit(-1);
	}


	EquilFns orig1, final1;
	EquilFns orig2, final2;

	VecDoub myphis = VecDoub(PHI_STATES);
	for (int i = 0; i < PHI_STATES; i++){
		myphis[i] = phis[i];
	}

	initialize(orig1, myphis);
	std::cout << "Solving Country 1" << std::endl;
	solve(myphis, prices, orig1, final1);
#if 1
	printResults(final1, myphis, C1);
#else
	std::cout << "Skipping printing" << std::endl;
#endif
	State current1(myphis, prices);
	StochProc stoch1(myphis);
	current1.defaultInitialState(stoch1);

	for (int i = PHI_STATES; i < 2 * PHI_STATES; i++){
		myphis[i - PHI_STATES] = phis[i];
	}

	initialize(orig2, myphis);
	std::cout << "Solving Country 2" << std::endl;
	solve(myphis, prices, orig2, final2);
#if 1
	printResults(final2, myphis, C2);
#else
	std::cout << "Skipping printing" << std::endl;
#endif

	State current2(myphis, prices);
	StochProc stoch2(myphis);
	current2.defaultInitialState(stoch2);

	std::cout << "Simulating world economies. " << NUMHHS << " households over " << NUMPERIODS << " periods." << std::endl;
	double xres = getNextParameters(final1, stoch1, current1, final2, stoch2, current2, c1prop);

	ostringstream os;
	ofstream out_stream;
	os << "intermediateResults.dat";
	out_stream.precision(15);
	out_stream << std::scientific;
	out_stream.open(os.str(), std::ofstream::out | std::ofstream::app);
	out_stream << xres << "          ";
	std::cout << "Tau: " << prices[0] << std::endl;
	out_stream << "Tau: " << prices[0] << std::endl;
	out_stream.close();

	ostringstream os2;
	ofstream out_stream2;
	os2 << "output" << seqNo << ".dat";
	out_stream2.precision(15);
	out_stream2 << std::scientific;
	out_stream2.open(os2.str(), std::ofstream::out | std::ofstream::trunc);
	out_stream2 << xres << endl;
	out_stream2.close();

	return xres;
}

double getNextParameters(const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target)
{
	int numC = 2;
	vector<MatrixXf> data = simulate(numC, policies1, stoch1, curSt1, policies2, stoch2, curSt2, c1target);

	MatrixXf badRHS = data[0].middleCols(0, numC+1);
	MatrixXf badLHS = data[0].middleCols(numC, 1);

	MatrixXf goodRHS = data[1].middleCols(0, numC+1);
	MatrixXf goodLHS = data[1].middleCols(numC, 1);

	ofstream out_stream;
	ostringstream os;
	os << "data.dat";
	out_stream.open(os.str());
	out_stream << "LHS" << endl;
	out_stream <<  badLHS << endl;
	out_stream << "RHS" << endl;
	out_stream << badRHS << endl;
	out_stream << "soln: " << endl << (badRHS.transpose() * badRHS).ldlt().solve(badRHS.transpose() * badLHS.col(0)) << endl;
	out_stream.close();
	exit(-1);
}

vector<MatrixXf> simulate(int numC, const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target)
{
	//setup
	WorldEconomy we = WorldEconomy(numC, 0);
	we.initialize(0, policies1, stoch1, curSt1);
	we.initialize(1, policies2, stoch2, curSt2);

	//run simulation
	std::cout << "Simulating to steady state " << std::endl << std::flush;
	we.simulateToSS();
	std::cout << "Simulating " << NUMPERIODS << " periods. " << std::endl << std::flush;
	we.simulateNPeriods(NUMPERIODS);
	std::cout << "Done simulating" << std::endl << std::flush;
	/*
	VecDoub targets = VecDoub(2);
	targets[0] = c1target;
	targets[1] = 1 - c1target;
	double mydist = we.distance(targets);
	*/
	vector<MatrixXf> mydata;
	std::cout << "resizing" << std::endl;
	mydata.resize(PHI_STATES);
	std::cout << "complete resize" << std::endl;
	vector<VecDoub> hist = we.getHistory();
	int numGoodStates = 0;
	int numBadStates = 0;
	for (int i = 0; i < NUMPERIODS-1; i++){
		std::cout << i << " ";
		if (hist[i+SSPERIODS][2] == 1){
			numGoodStates++;
		}
		else{
			numBadStates++;
		}
	}
	mydata[0] = MatrixXf(numBadStates,2*numC+1);
	mydata[1] = MatrixXf(numGoodStates, 2 * numC + 1);

	int currentGood = 0;
	int currentBad = 0;
	for (int i = 0; i < NUMPERIODS-1; i++){
		std::cout << std::endl << i << " ";
		if (hist[i + SSPERIODS][2] == 0){
			if (currentBad > numBadStates){
				std::cerr << "At bad state " << currentBad << " but only expect " << numBadStates;
				exit(-1);
			}
			mydata[0](currentBad, 0) = 1;
			for (int j = 0; j < numC; j++){
				mydata[0](currentBad, j+1) = hist[i+SSPERIODS][j];
				mydata[0](currentBad, j+numC+1) = hist[i+SSPERIODS+1][j];
			}
			currentBad++;
		}
		else{
			if (currentGood > numGoodStates){
				std::cerr << "At good state " << currentGood << " but only expect " << numGoodStates;
				exit(-1);
			}

			mydata[1](currentGood, 0) = 1;
			for (int j = 0; j < numC; j++){
				mydata[1](currentGood, j+1) = hist[i+SSPERIODS][j];
				mydata[1](currentGood, j + numC + 1) = hist[i + SSPERIODS + 1][j];
			}
			currentGood++;
		}
	}

	return mydata;
}

void printResults(const EquilFns& e, const VecDoub& phis, COUNTRYID whichCountry) {
	using namespace std;
	ofstream out_stream;

	ostringstream os;
	os << "policy_" << whichCountry << ".dat";

	StochProc stoch(phis);
	out_stream.open(os.str());

	out_stream
		<< "aggAsset,aggAsset2,phi,agg_shock,z1,z2,wage_shock,a,value_fn,consumption,k1_prime,k2_prime,b_prime,c1mgmt"
		<< endl;
	for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int f = 0; f < AGG_ASSET_SIZE; f++) {
				for (int j = 0; j < ASSET_SIZE; j++) {
#if IID_CAP==0
					for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
						for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
#else
					int l = 0;
					int ll = 0;
#endif
					for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
						for (int h = 0; h < PHI_STATES; h++) {
							out_stream
								<< stoch.aggAssets[g] << ","
								<< stoch.aggAssets[f] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_PHI] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_A] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_K1] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_K2] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_W] << ","
								<< stoch.assets[j] << ","
								<< e.value_fn[f][g][h][i][j][l][ll][m] << ","
								<< e.consumption[f][g][h][i][j][l][ll][m] << ","
								<< e.policy_fn[f][g][h][i][j][l][ll][m][K1STATE] << ","
								<< e.policy_fn[f][g][h][i][j][l][ll][m][K2STATE] << ","
								<< e.policy_fn[f][g][h][i][j][l][ll][m][BSTATE] //<< ","
								//						<< e.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE]
								<< endl;
						}
#if IID_CAP==0
					}
						}
#endif
					}
				}
			}
		}
	}
	out_stream.close();
}

void initialize(EquilFns& fns, const VecDoub& phis) {

	using namespace std;

	StochProc stoch(phis);
	for (int f = 0; f < AGG_ASSET_SIZE; f++){
		for (int g = 0; g < AGG_ASSET_SIZE; g++){
			for (int h = 0; h < PHI_STATES; h++){
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int j = 0; j < ASSET_SIZE; j++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									fns.value_fn[f][g][h][i][j][l][ll][m] = vfiMaxUtil::consUtil(stoch.assets[j] - MIN_ASSETS + 0.1);
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}

void solve(const VecDoub& phis, const VecDoub& prices, const EquilFns& orig, EquilFns& final) {
	using namespace std;
	EquilFns in_process_values, last_values;

	StochProc stoch(phis);

	int counter = 0;
	double diff = 100;
	double diff2 = 100;

	last_values = orig;

	while (counter < MAX_ITER && diff > VAL_TOL) {

		counter++;
#if OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(3)
#endif
		for (int j = 0; j < ASSET_SIZE; j++) {
			for (int g = 0; g < AGG_ASSET_SIZE; g++) {
				for (int gg = 0; gg < AGG_ASSET_SIZE; gg++) {
					for (int h = 0; h < PHI_STATES; h++){
						for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
							State current(phis, prices);
							vfiMaxUtil ub = vfiMaxUtil(current, stoch, last_values);
							//				eulerUtility ub = eulerUtility(current, stoch, last_values);

#if AMOEBA
							Amoeba am(MINIMIZATION_TOL);
#endif
							for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
								for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
									for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
										current.current_states[ASTATE] = stoch.assets[j];
										current.current_states[AGG_ASSET_STATE] = stoch.aggAssets[g];
										current.current_states[AGG2_ASSET_STATE] = stoch.aggAssets[gg];
										current.current_states[CAP1_SHOCK_STATE] =
											stoch.shocks[h][i][l][ll][m][EF_K1];
										current.current_states[CAP2_SHOCK_STATE] =
											stoch.shocks[h][i][l][ll][m][EF_K2];
										current.current_states[WAGE_SHOCK_STATE] =
											stoch.shocks[h][i][l][ll][m][EF_W];
										current.current_states[AGG_SHOCK_STATE] =
											stoch.shocks[h][i][l][ll][m][EF_A];
										current.current_states[PHI_STATE] = stoch.shocks[h][i][l][ll][m][EF_PHI];
										current.current_indices[ASTATE] = j;
										current.current_indices[AGG_ASSET_STATE] = g;
										current.current_indices[AGG2_ASSET_STATE] = gg;
										current.current_indices[CAP1_SHOCK_STATE] = l;
										current.current_indices[CAP2_SHOCK_STATE] = ll;
										current.current_indices[WAGE_SHOCK_STATE] = m;
										current.current_indices[AGG_SHOCK_STATE] = i;
										current.current_indices[PHI_STATE] = h;

										ub.updateCurrent(current);

										VecDoub temp1(NUM_CHOICE_VARS);
										temp1[K1STATE] =
											sqrt(
											MAX(MIN_CAPITAL, last_values.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE])
											- MIN_CAPITAL);
										temp1[K2STATE] =
											sqrt(
											MAX(MIN_CAPITAL, last_values.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE])
											- MIN_CAPITAL);
										temp1[BSTATE] =
											sqrt(
											MAX(MIN_BONDS, last_values.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE])
											- MIN_BONDS);
										/*
																		temp1[MGMT_C1_STATE] =
																		sqrt(
																		MAX(MIN_MGMT, last_values.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE])
																		- MIN_MGMT);
																		*/
										VecDoub_I initial_point(temp1);
										VecDoub temp2;
										double minVal = 0;
										//#if IID_CAP==0
										//								std::cerr<<"akmodel.cc:solve - Error. Expect iid capital shocks."<<std::endl;
										//								exit(-1);
										//#endif
										if (l == 0 && ll == 0){
											if (ub.constraintBinds()){
												temp2 = VecDoub(NUM_CHOICE_VARS);
												temp2[K1STATE] = 0;
												temp2[K2STATE] = 0;
												temp2[BSTATE] = sqrt(ub.getBoundBorrow() - MIN_BONDS);
												//temp2[MGMT_C1_STATE] = 7;
												minVal = ub.getBoundUtil();
											}
											else{
												//#if	SIMULATED_ANNEALING
												//										cerr<<"cannot set both amoeba and simulated annealing"<<endl;
												//										exit(-1);
												//#endif
#if AMOEBA
												const double delta = MAX(0.001, 0.1 * initial_point[K1STATE]);
												temp2 = am.minimize(initial_point, delta, ub);
												minVal = ((vfiMaxUtil)ub)(temp2);
#else
												SimulatedAnnealingWorld saw(initial_point, 100, MAX_ITER/counter, MINIMIZATION_TOL, 0.9,
													&current,&stoch,&last_values,&utility, 1.0/(counter/10+1));
												VecDoub *temp3=saw.solve();
												minVal=saw.getEnergy(*temp3);
												temp2=*temp3;
												delete temp3;
#endif
											}
											//#if IID_CAP==0
											//								std::cerr<<"akmodel.cc:solve - Error. Expect iid capital shocks."<<std::endl;
											//								exit(-1);
											//#endif
										}
										else {
											temp2 = VecDoub(NUM_CHOICE_VARS);
											temp2[K1STATE] =
												sqrt(
												in_process_values.policy_fn[g][gg][h][i][j][0][0][m][K1STATE]
												- MIN_CAPITAL);
											temp2[K2STATE] =
												sqrt(
												in_process_values.policy_fn[g][gg][h][i][j][0][0][m][K2STATE]
												- MIN_CAPITAL);
											temp2[BSTATE] =
												sqrt(
												in_process_values.policy_fn[g][gg][h][i][j][0][0][m][BSTATE]
												- MIN_BONDS);
											/*
											temp2[MGMT_C1_STATE] =
											sqrt(
											in_process_values.policy_fn[h][i][j][0][0][m][MGMT_C1_STATE]
											- MIN_MGMT);
											*/
											minVal =
												-in_process_values.value_fn[g][gg][h][i][j][0][0][m];
										}
										//#endif

										in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE] =
											MIN_CAPITAL + utilityFunctions::integer_power(temp2[K1STATE], 2);
										in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE] =
											MIN_CAPITAL + utilityFunctions::integer_power(temp2[K2STATE], 2);
										in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE] =
											MIN_BONDS + utilityFunctions::integer_power(temp2[BSTATE], 2);
										/*
										in_process_values.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE] =
										MIN_MGMT + utilityFunctions::integer_power(temp2[MGMT_C1_STATE], 2);
										in_process_values.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE] =
										utilityFunctions::boundValue(in_process_values.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE], MIN_MGMT, MAX_MGMT);
										if (temp2[K1STATE] == 0 && temp2[K2STATE] == 0){
										in_process_values.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE] = 50;
										}
										*/
										in_process_values.value_fn[g][gg][h][i][j][l][ll][m] = -minVal;
										in_process_values.consumption[g][gg][h][i][j][l][ll][m] =
											current.current_states[ASTATE]
											- (MIN_CAPITAL
											+ utilityFunctions::integer_power(temp2[K1STATE], 2))
											- (MIN_CAPITAL
											+ utilityFunctions::integer_power(temp2[K2STATE], 2))
											- (MIN_BONDS + utilityFunctions::integer_power(temp2[BSTATE], 2));
									}
								}
							}
						}
					}
				}
			}
		}

		{
			VecInt location1 = VecInt(NUM_STATE_VARS);
			VecInt location2 = VecInt(NUM_STATE_VARS);
#if	POLICY_CONVERGENCE
			diff = maxPolicyDistance(last_values, in_process_values, location1);
			diff2 = maxValueDistance(last_values, in_process_values, location2);
#else
			diff = maxValueDistance(last_values, in_process_values, location1);
			diff2 = maxPolicyDistance(last_values, in_process_values, location2);
#endif
			int max_iterations = 0;
#if POLICY_ITERATION
			State current(phis, prices);
			max_iterations = policy_iteration(diff, current, stoch, last_values, in_process_values);
#if	POLICY_CONVERGENCE
			diff = maxPolicyDistance(last_values, in_process_values, location1);
			diff2 = maxValueDistance(last_values, in_process_values, location2);
#else
			diff = maxValueDistance(last_values, in_process_values, location1);
			diff2 = maxPolicyDistance(last_values, in_process_values, location2);
#endif
#endif

			std::cout
				<< ".......................................................................................Iteration "
				<< counter << std::endl
#if POLICY_CONVERGENCE
				<< "policy diff=" << diff << "       Value diff=" << diff2 << std::endl;
#else
				<< "policy diff=" << diff2 << "       Value diff=" << diff << std::endl;
#endif
#if 0
			std::cout << std::scientific;
			std::cout << " Params: phi=" << phis[0] << ":" << phis[1]
				<< " p1=" << prices[0] << ":" << prices[3]
				<< " p2=" << prices[1] << ":" << prices[4]
				<< " r=" << prices[2] << ":" << prices[5]
				<< endl << diff << "         " << diff2;
			std::cout.unsetf(std::ios_base::floatfield);
			cout << " , " << location1[0] << ":"
				<< location1[1] << ":" << location1[2] << ":"
				<< location1[3] << "  "
				<< last_values.value_fn[0][location1[0]][location1[1]][location1[2]][location1[3]][0]
				<< "  "
				<< in_process_values.value_fn[0][location1[0]][location1[1]][location1[2]][location1[3]][0]

				<< endl << std::scientific << diff2;
			cout.unsetf(std::ios_base::floatfield);
			cout << " , " << location2[0] << ":"
				<< location2[1] << ":" << location2[2] << ":"
				<< location2[3] << "  "
				<< last_values.policy_fn[0][location2[0]][location2[1]][location2[2]][location2[3]][0]
				<< "  "
				<< last_values.policy_fn[0][location2[0]][location2[1]][location2[2]][location2[3]][1]
				<< "  "
				<< last_values.policy_fn[0][location2[0]][location2[1]][location2[2]][location2[3]][2]
				<< "  "
				<< last_values.policy_fn[0][location2[0]][location2[1]][location2[2]][location2[3]][3]
				<< "  "
				<< in_process_values.policy_fn[0][location2[0]][location2[1]][location2[2]][location2[3]][0]
				<< "  "
				<< in_process_values.policy_fn[0][location2[0]][location2[1]][location2[2]][location2[3]][1]
				<< "  "
				<< in_process_values.policy_fn[0][location2[0]][location2[1]][location2[2]][location2[3]][2]
				<< "  "
				<< in_process_values.policy_fn[0][location2[0]][location2[1]][location2[2]][location2[3]][3]
				std::cout << endl;
#endif
			std::cout << max_iterations << endl;
		}

		last_values = in_process_values;
	}

	final = last_values;

	return;
}

int policy_iteration(double difference, State& current, const StochProc& stoch,
	EquilFns& oldValue, EquilFns& latestValue) {

	int max_iterations = MIN(20, floor(log(1.0 / difference)));

	if (max_iterations < 1){
		return max_iterations;
	}

	EquilFns temp(latestValue);
	EquilFns altPol(latestValue);

	for (int g = 0; g < AGG_ASSET_SIZE; g++){
		for (int gg = 0; gg < AGG_ASSET_SIZE; gg++){
			for (int h = 0; h < PHI_STATES; h++){
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
						for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
							for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
								for (int j = 0; j < ASSET_SIZE; j++) {
									altPol.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE] = sqrt(
										temp.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE]
										- MIN_CAPITAL);
									altPol.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE] = sqrt(
										temp.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE]
										- MIN_CAPITAL);
									altPol.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE] = sqrt(
										temp.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE]
										- MIN_BONDS);
									/*
									altPol.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE] = sqrt(
									temp.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE]
									- MIN_MGMT);
									*/
								}
							}
						}
					}
				}
			}
		}
	}

	for (int count = 0; count < max_iterations; count++){
		vfiMaxUtil ub = vfiMaxUtil(current, stoch, latestValue);
		//eulerUtility vmu = eulerUtility(current,stoch, initial);

		for (int g = 0; g < AGG_ASSET_SIZE; g++){
			for (int gg = 0; gg < AGG_ASSET_SIZE; gg++){
				for (int h = 0; h < PHI_STATES; h++){
					for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									for (int j = 0; j < ASSET_SIZE; j++) {
										current.current_states[ASTATE] = stoch.assets[j];
										current.current_states[AGG_ASSET_STATE] = stoch.aggAssets[g];
										current.current_states[AGG2_ASSET_STATE] = stoch.aggAssets[gg];
										current.current_states[CAP1_SHOCK_STATE] =
											stoch.shocks[h][i][l][ll][m][EF_K1];
										current.current_states[CAP2_SHOCK_STATE] =
											stoch.shocks[h][i][l][ll][m][EF_K2];
										current.current_states[WAGE_SHOCK_STATE] =
											stoch.shocks[h][i][l][ll][m][EF_W];
										current.current_states[AGG_SHOCK_STATE] =
											stoch.shocks[h][i][l][ll][m][EF_A];
										current.current_states[PHI_STATE] = stoch.shocks[h][i][l][ll][m][EF_PHI];
										current.current_indices[ASTATE] = j;
										current.current_indices[AGG_ASSET_STATE] = g;
										current.current_indices[AGG2_ASSET_STATE] = gg;
										current.current_indices[CAP1_SHOCK_STATE] = l;
										current.current_indices[CAP2_SHOCK_STATE] = ll;
										current.current_indices[WAGE_SHOCK_STATE] = m;
										current.current_indices[AGG_SHOCK_STATE] = i;
										current.current_indices[PHI_STATE] = h;

										ub.updateCurrent(current);

										VecDoub policy = VecDoub(NUM_CHOICE_VARS);
										policy[K1STATE] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE];
										policy[K2STATE] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE];
										policy[BSTATE] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE];
										//policy[MGMT_C1_STATE] = altPol.policy_fn[h][i][j][l][ll][m][MGMT_C1_STATE];
										if (l == 0 && ll == 0){
											temp.value_fn[g][gg][h][i][j][l][ll][m] = -ub(policy);
										}
										else{
											temp.value_fn[g][gg][h][i][j][l][ll][m] = temp.value_fn[g][gg][h][i][j][0][0][m];
										}
									}
								}
							}
						}
					}
				}
			}
		}

		for (int g = 0; g < AGG_ASSET_SIZE; g++){
			for (int gg = 0; gg < AGG_ASSET_SIZE; gg++){
				for (int h = 0; h < PHI_STATES; h++){
					for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									for (int j = 0; j < ASSET_SIZE; j++) {
										if (count == (max_iterations - 2)){
											oldValue.value_fn[g][gg][h][i][j][l][ll][m] = temp.value_fn[g][gg][h][i][j][l][ll][m];
										}
										latestValue.value_fn[g][gg][h][i][j][l][ll][m] = temp.value_fn[g][gg][h][i][j][l][ll][m];
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return max_iterations;
}


double maxValueDistance(const EquilFns &a, const EquilFns &b,
	VecInt& location) {
	VecDoub val1(1);
	val1[0] = a.value_fn[0][0][0][0][0][0][0][0];
	VecDoub val2(1);
	val2[0] = b.value_fn[0][0][0][0][0][0][0][0];
	double temp_var = 0;
	double diff = utilityFunctions::distanceFn(val1, val2);
	location[ASTATE] = 0;
	location[AGG_ASSET_STATE] = 0;
	location[AGG2_ASSET_STATE] = 0;
	location[CAP1_SHOCK_STATE] = 0;
	location[CAP2_SHOCK_STATE] = 0;
	location[WAGE_SHOCK_STATE] = 0;
	location[AGG_SHOCK_STATE] = 0;
	location[PHI_STATE] = 0;
	for (int g = 0; g < AGG_ASSET_SIZE; g++) {
		for (int gg = 0; gg < AGG_ASSET_SIZE; gg++) {
			for (int h = 0; h < PHI_STATES; h++) {
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int j = 0; j < ASSET_SIZE; j++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									val1[0] = a.value_fn[g][gg][h][i][j][l][ll][m];
									val2[0] = b.value_fn[g][gg][h][i][j][l][ll][m];
									temp_var = utilityFunctions::distanceFn(val1, val2);
									if (temp_var > diff) {
										diff = temp_var;
										location[ASTATE] = j;
										location[CAP1_SHOCK_STATE] = l;
										location[CAP2_SHOCK_STATE] = ll;
										location[WAGE_SHOCK_STATE] = m;
										location[AGG_SHOCK_STATE] = i;
										location[PHI_STATE] = h;
										location[AGG_ASSET_STATE] = g;
										location[AGG2_ASSET_STATE] = gg;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return diff;
}

double maxPolicyDistance(const EquilFns &a, const EquilFns &b,
	VecInt& location) {
	VecDoub val1 = a.policy_fn[0][0][0][0][0][0][0][0];
	VecDoub val2 = b.policy_fn[0][0][0][0][0][0][0][0];
	//val1[MGMT_C1_STATE] = val1[MGMT_C1_STATE] / MAX_MGMT;
	//val2[MGMT_C1_STATE] = val2[MGMT_C1_STATE] / MAX_MGMT;
	double temp_var = 0;
	double diff = utilityFunctions::distanceFn(val1, val2);
	location[ASTATE] = 0;
	location[AGG_ASSET_STATE] = 0;
	location[AGG2_ASSET_STATE] = 0;
	location[CAP1_SHOCK_STATE] = 0;
	location[CAP2_SHOCK_STATE] = 0;
	location[WAGE_SHOCK_STATE] = 0;
	location[AGG_SHOCK_STATE] = 0;
	location[PHI_STATE] = 0;
	for (int g = 0; g < AGG_ASSET_SIZE; g++) {
		for (int gg = 0; gg < AGG_ASSET_SIZE; gg++) {
			for (int h = 0; h < PHI_STATES; h++) {
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int j = 0; j < ASSET_SIZE; j++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									val1 = a.policy_fn[g][gg][h][i][j][l][ll][m];
									val2 = b.policy_fn[g][gg][h][i][j][l][ll][m];
									//val1[MGMT_C1_STATE] = val1[MGMT_C1_STATE] / MAX_MGMT;
									//val2[MGMT_C1_STATE] = val2[MGMT_C1_STATE] / MAX_MGMT;
									temp_var = utilityFunctions::distanceFn(val1, val2);
									if (temp_var > diff) {
										diff = temp_var;
										location[ASTATE] = j;
										location[AGG_ASSET_STATE] = g;
										location[AGG2_ASSET_STATE] = gg;
										location[CAP1_SHOCK_STATE] = l;
										location[CAP2_SHOCK_STATE] = ll;
										location[WAGE_SHOCK_STATE] = m;
										location[AGG_SHOCK_STATE] = i;
										location[PHI_STATE] = h;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return diff;
}



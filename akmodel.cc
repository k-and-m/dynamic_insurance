/*
 * akmodel.cc
 *
 *  Created on: Feb 11, 2015
 *      Author: robmrk
 */

#define POLICY_ITERATION 1
#define AMOEBA 1
#define BFGS 0
#define SIMULATED_ANNEALING 0
#define POLICY_CONVERGENCE 0

#include "BaseHeader.h"
#include "utilityFunctions.h"
#include "matrixIO.h"
#include "matrix.h"
#if 0
#include "amoeba.h"
#else
#include "simplex.h"
#endif
#include "WorldEconomy.h"
#include "vfiMaxUtil.h"
#include <Eigen/Dense>
#if USE_MPI
#include "mpi.h"
#endif
#if BFGS
#include <dlib/optimization.h>
using namespace dlib;
#endif

#if 0
#if SIMULATED_ANNEALING
#include "simAnneal.hpp"
#endif
#endif

using namespace Eigen;

double simulate(int numC, const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target);
double getNextParameters(const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target);
void printResults(const EquilFns& e, const VecDoub& phis, const COUNTRYID whichCountry);
void readResults(EquilFns& e, COUNTRYID whichCountry);
void initialize(EquilFns& fns, const VecDoub& phis);
void solve(const VecDoub& phis, const VecDoub& prices, const EquilFns& orig, EquilFns& final, const COUNTRYID whichCountry);
int policy_iteration(double difference, State& initial, const StochProc& stoch,
	EquilFns& lastFns, EquilFns& currentFns);
double maxValueDistance(const EquilFns &a, const EquilFns &b, VecInt& location);
double maxPolicyDistance(const EquilFns &a, const EquilFns &b,
	VecInt& location);
double solveProblem(const VecDoub& phis, const VecDoub& prices, double c1prop, int seqNo);

#if 1
int main(int argc, char *argv[]) {
	using namespace std;

	if (argc != 7) {
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
	double dis = solveProblem(phi, prices, 0.3, atoi(argv[argc - 1]));
#else 
	VecDoub phi = VecDoub(2);
	VecDoub tau = VecDoub(2);
	phi[0] = 1;
	phi[1] = 0;

	tau[0] = 1.1;
	tau[1] = atof(argv[1]);
	double dis = solveProblem(phi, tau, 0.3, 0);
#endif
	std::cout << "soln dist: " << dis << endl;

	return 0;
}
#endif

#if 0
extern "C" {
	double runonce_(double *c1p, double *c2p, double *c1Price, double *c2Price, double *intRate,
		double *c1prop)
	{
		return solveProblem(*c1p, *c2p, *c1Price, *c2Price, *intRate, *c1prop);
	}
}
#endif

double solveProblem(const VecDoub& phis, const VecDoub& tau, double c1prop, int seqNo)
{
	if (phis.size() != 2 * PHI_STATES) {
		std::cerr << "solveProblem: incorrect number of phis " << phis.size() << endl;
		exit(-1);
	}
	if (tau.size() != 2) {
		std::cerr << "solveProblem: incorrect number of prices (taus)" << tau.size() << endl;
		exit(-1);
	}

	double avgDist = 10;
	for (int loopC = 0; (loopC < 1) && (avgDist>VAL_TOL); loopC++) {

		EquilFns orig1, final1;
		EquilFns orig2, final2;

		VecDoub myphis = VecDoub(PHI_STATES);
		for (int i = 0; i < PHI_STATES; i++) {
			myphis[i] = phis[i];
		}

		StochProc stoch1(myphis);
		initialize(orig1, myphis);
		std::cout << "Solving Country 1" << std::endl;
		solve(myphis, tau, orig1, final1, C1);

		State current1(myphis, tau);
		current1.defaultInitialState(stoch1);

		for (int i = PHI_STATES; i < 2 * PHI_STATES; i++) {
			myphis[i - PHI_STATES] = phis[i];
		}
		StochProc stoch2(myphis);
		State current2(myphis, tau);
		current2.defaultInitialState(stoch2);
		initialize(orig2, myphis);

		std::cout << "Solving Country 2" << std::endl;
		solve(myphis, tau, orig2, final2, C2);

		std::cout << "Simulating world economies. " << NUMHHS << " households over " << TOTALPERIODS << " periods." << std::endl;

		int numC = 2;
		double avgDist = getNextParameters(final1, stoch1, current1, final2, stoch2, current2, c1prop);
	}

	ostringstream os;
	ofstream out_stream;
	os << "intermediateResults.dat";
	out_stream.precision(15);
	out_stream << std::scientific;
	out_stream.open(os.str(), std::ofstream::out | std::ofstream::app);
	out_stream << avgDist << "          ";
	std::cout << "Tau: " << tau[0] << std::endl;
	out_stream << "Tau: " << tau[0] << std::endl;
	out_stream.close();

	ostringstream os2;
	ofstream out_stream2;
	os2 << "output" << seqNo << ".dat";
	out_stream2.precision(15);
	out_stream2 << std::scientific;
	out_stream2.open(os2.str(), std::ofstream::out | std::ofstream::trunc);
	out_stream2 << avgDist << endl;
	out_stream2.close();

	return avgDist;
}

double getNextParameters(const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target)
{
	int numC = 2;
	double results = simulate(numC, policies1, stoch1, curSt1, policies2, stoch2, curSt2, c1target);
	return results;
}

double simulate(int numC, const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
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

	return we.netBonds();
}

void printResults(const EquilFns& e, const VecDoub& phis, const COUNTRYID whichCountry) {
	using namespace std;
	ofstream out_stream;

	ostringstream os;
	os << "policy_" << whichCountry << ".dat";

	StochProc stoch(phis);
	out_stream.open(os.str());

	out_stream.precision(15);
	for (int i = 0; i < PHI_STATES; i++) {
		if (i == (PHI_STATES - 1)) {
			out_stream << phis[i] << std::endl;
		}
		else {
			out_stream << phis[i] << ",";
		}
	}

	for (int i = 0; i < CAP_SHOCK_SIZE; i++) {
		if (i == (CAP_SHOCK_SIZE - 1)) {
			out_stream << stoch.shocks[0][0][i][0][0][EF_K1] << std::endl;
		}
		else {
			out_stream << stoch.shocks[0][0][i][0][0][EF_K1] << ",";
		}
	}

	for (int i = 0; i < WAGE_SHOCK_SIZE; i++) {
		if (i == (WAGE_SHOCK_SIZE - 1)) {
			out_stream << stoch.shocks[0][0][0][0][i][EF_W] << std::endl;
		}
		else {
			out_stream << stoch.shocks[0][0][0][0][i][EF_W] << ",";
		}
	}

	out_stream
		<< "phi,agg_shock,z1,z2,wage_shock,a,value_fn,consumption,k1_prime,k2_prime,b_prime"
		<< endl;

	VecInt index(6);

	for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
		for (int j = 0; j < ASSET_SIZE; j++) {
			for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
				for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
					for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
						for (int h = 0; h < PHI_STATES; h++) {
							index[0] = h;
							index[1] = i;
							index[2] = j;
							index[3] = l;
							index[4] = ll;
							index[5] = m;
#if 1
							out_stream
								<< h << ","
								<< i << ","
								<< l << ","
								<< ll << ","
								<< m << ","
								<< j << ","
								<< e.getValueFn(index) << ","
								<< e.consumption[h][i][j][l][ll][m] << ","
								<< e.policy_fn[h][i][j][l][ll][m][K1STATE] << ","
								<< e.policy_fn[h][i][j][l][ll][m][K2STATE] << ","
								<< e.policy_fn[h][i][j][l][ll][m][BSTATE]
								<< endl;
#else
							out_stream
								<< stoch.aggAssets[g] << ","
								<< stoch.aggAssets[f] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_PHI] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_A] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_K1] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_K2] << ","
								<< stoch.shocks[h][i][l][ll][m][EF_W] << ","
								<< stoch.assets[j] << ","
								<< e.getValueFn(index) << ","
								<< e.consumption[f][g][h][i][j][l][ll][m] << ","
								<< e.policy_fn[f][g][h][i][j][l][ll][m][K1STATE] << ","
								<< e.policy_fn[f][g][h][i][j][l][ll][m][K2STATE] << ","
								<< e.policy_fn[f][g][h][i][j][l][ll][m][BSTATE]
								<< endl;
#endif
						}
					}
				}
			}
		}
	}
	out_stream.close();
}

void readResults(EquilFns& e, COUNTRYID whichCountry) {
	using namespace std;

	FILE *fp = NULL;
	if (whichCountry == C1) {
		fp = std::fopen("policy_0.dat", "r");
	}
	else {
		fp = std::fopen("policy_1.dat", "r");
	}

	int a1, a2, aggShock, phi, z1, z2, wage, a;
	double valuefn, cons, k1p, k2p, bp;
	double phi1, phi2;

	if (PHI_STATES != 2) {
		std::cerr << " ERROR! akmodel.cc-readResults(): only expect two phi states, not " << PHI_STATES << std::endl;
		std::fclose(fp);
		exit(-1);
	}

	VecDoub vphis(PHI_STATES);
	if (fscanf(fp, "%lf, %lf", &phi1, &phi2) != PHI_STATES) {
		std::cerr << "ERROR! akmodel.cc:readResults() - unable to read phis." << std::endl;
		std::fclose(fp);
		exit(-1);
	}
	bool doublePhi = (phi1 == phi2);
	vphis[0] = phi1;
	vphis[1] = phi2;
	StochProc stoch(vphis);

	if (CAP_SHOCK_SIZE != 2) {
		std::cerr << " ERROR! akmodel.cc-readResults(): only expect two capital states, not " << CAP_SHOCK_SIZE << std::endl;
		std::fclose(fp);
		exit(-1);
	}

	VecDoub vzs(CAP_SHOCK_SIZE);
	if (fscanf(fp, "%lf, %lf", &phi1, &phi2) != CAP_SHOCK_SIZE) {
		std::cerr << "ERROR! akmodel.cc:readResults() - unable to read z's." << std::endl;
		std::fclose(fp);
		exit(-1);
	}
	vzs[0] = phi1;
	vzs[1] = phi2;

	VecDoub vws(WAGE_SHOCK_SIZE);
	if (fscanf(fp, "%lf, %lf", &phi1, &phi2) != WAGE_SHOCK_SIZE) {
		std::cerr << "ERROR! akmodel.cc:readResults() - unable to read w's." << std::endl;
		std::fclose(fp);
		exit(-1);
	}
	vws[0] = phi1;
	vws[1] = phi2;

	if (1) {
		char str[256];
		fscanf(fp, "%s", str);
	}
	while (fscanf(fp, "%d, %d, %d, %d, %d, %d, %lf, %lf, %lf, %lf, %lf",
		&phi, &aggShock, &z1, &z2, &wage, &a, &valuefn, &cons, &k1p, &k2p, &bp) == 13) {

		VecInt index(6);

		index[0] = phi;
		index[1] = aggShock;
		index[2] = a;
		index[3] = z1;
		index[4] = z2;
		index[5] = wage;

		e.setValueFn(index, valuefn);
		e.consumption[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]] = cons;
		e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][K1STATE] = k1p;
		e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][K2STATE] = k2p;
		e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][BSTATE] = bp;

		if (doublePhi) {
			index[0] = 1;
			e.setValueFn(index, valuefn);
			e.consumption[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]] = cons;
			e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][K1STATE] = k1p;
			e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][K2STATE] = k2p;
			e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][BSTATE] = bp;
		}
	}
	fclose(fp);
}


void initialize(EquilFns& fns, const VecDoub& phis) {

	using namespace std;

	StochProc stoch(phis);
	VecInt vect(6);
	for (int h = 0; h < PHI_STATES; h++) {
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			for (int j = 0; j < ASSET_SIZE; j++) {
				for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
					for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
						for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
							vect[0] = h;
							vect[1] = i;
							vect[2] = j;
							vect[3] = l;
							vect[4] = ll;
							vect[5] = m;
							//fns.setValueFn(vect, vfiMaxUtil::consUtil(stoch.assets[j] - MIN_ASSETS + 0.1));
							fns.setValueFn(vect, 0);
						}
					}
				}
			}
		}
	}
	return;
}

void solve(const VecDoub& phis, const VecDoub& prices, const EquilFns& orig, EquilFns& final, const COUNTRYID whichCountry) {
	using namespace std;


	EquilFns last_values;
	StochProc stoch(phis);

	int counter = 0;
	double diff = 100;
	double diff2 = 100;

	last_values = orig;

	while (counter < MAX_ITER && diff > VAL_TOL) {

		counter++;
		EquilFns in_process_values;

#if OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(3)
#endif
		for (int j = 0; j < ASSET_SIZE; j++) {
			for (int h = 0; h < PHI_STATES; h++) {
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					State current(phis, prices);
					vfiMaxUtil ub = vfiMaxUtil(current, stoch, last_values);
#if AMOEBA
#if 0
					Amoeba am(MINIMIZATION_TOL);
#endif
#endif
					for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								VecInt vect(6);
								vect[0] = h;
								vect[1] = i;
								vect[2] = j;
								vect[3] = l;
								vect[4] = ll;
								vect[5] = m;

#if 0
                                                                std::cout<<h<<":"<<i<<":"<<j<<":"<<l<<":"<<ll<<":"<<m<<std::endl<<std::flush;
#endif

								current.current_states[ASTATE] = stoch.assets[j];
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
								current.current_indices[CAP1_SHOCK_STATE] = l;
								current.current_indices[CAP2_SHOCK_STATE] = ll;
								current.current_indices[WAGE_SHOCK_STATE] = m;
								current.current_indices[AGG_SHOCK_STATE] = i;
								current.current_indices[PHI_STATE] = h;

								ub.updateCurrent(current);

#if K2CHOICE
								VecDoub temp1(NUM_CHOICE_VARS);
#else
								VecDoub temp1(NUM_CHOICE_VARS - 1);
#endif
								temp1[K1STATE] =
									sqrt(
										MAX(MIN_CAPITAL+0.00001, last_values.policy_fn[h][i][j][l][ll][m][K1STATE])
										- MIN_CAPITAL);
#if K2CHOICE
								temp1[K2STATE] =
									sqrt(
										MAX(MIN_CAPITAL, last_values.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE])
										- MIN_CAPITAL);
								temp1[BSTATE] =
#else
								temp1[BSTATE - 1] =
#endif
									sqrt(
										MAX(MIN_BONDS+1, last_values.policy_fn[h][i][j][l][ll][m][BSTATE])
										- MIN_BONDS);
								VecDoub_I initial_point(temp1);
								VecDoub temp2;
								double minVal = 0;
									if (ub.constraintBinds()) {
#if K2CHOICE
										temp2 = VecDoub(NUM_CHOICE_VARS);
#else
										temp2 = VecDoub(NUM_CHOICE_VARS - 1);
#endif
										temp2[K1STATE] = 0;
#if K2CHOICE
										temp2[K2STATE] = 0;
										temp2[BSTATE] = sqrt(ub.getBoundBorrow() - MIN_BONDS);
#else
										temp2[BSTATE - 1] = sqrt(ub.getBoundBorrow() - MIN_BONDS);
#endif
										minVal = ub.getBoundUtil();
									}
									else {
#if AMOEBA
										const double delta = MAX(0.001, 0.1 * initial_point[K1STATE]);
										if (initial_point.size() != 2) {
											std::cerr << "akmodel.cc-solve(): initial point is not of size 2, but of size " << initial_point.size();
											exit(-1);
										}
#if 0
										temp2 = am.minimize(initial_point, delta, ub);
#else
temp2 = BT::Simplex(ub, initial_point,MINIMIZATION_TOL);
#endif
										if (temp2.size() != 2) {
											std::cerr << "ERROR! akmodel.cc - solve(): minimize returned std::vector of size " << temp2.size() << std::endl;
											exit(-1);
										}
										minVal = ((vfiMaxUtil)ub)(temp2);
#elif BFGS
										find_min_using_approximate_derivatives(bfgs_search_strategy(),
											objective_delta_stop_strategy(1e-7),
											ub, initial_point, -1);
#else
										SimulatedAnnealingWorld saw(initial_point, 100, MAX_ITER / counter, MINIMIZATION_TOL, 0.9,
											&current, &stoch, &last_values, &utility, 1.0 / (counter / 10 + 1));
										VecDoub *temp3 = saw.solve();
										minVal = saw.getEnergy(*temp3);
										temp2 = *temp3;
										delete temp3;
#endif
									}

								in_process_values.policy_fn[h][i][j][l][ll][m][K1STATE] =
									MIN_CAPITAL + utilityFunctions::integer_power(temp2[K1STATE], 2);
								in_process_values.policy_fn[h][i][j][l][ll][m][K2STATE] =
#if K2CHOICE
									MIN_CAPITAL + utilityFunctions::integer_power(temp2[K2STATE], 2);
#else
									in_process_values.policy_fn[h][i][j][l][ll][m][K1STATE]
									* pow(prices[0], 1 / (1 - ALPHA1));
#endif
								in_process_values.policy_fn[h][i][j][l][ll][m][BSTATE] =
#if K2CHOICE
									MIN_BONDS + utilityFunctions::integer_power(temp2[BSTATE], 2);
#else
									MIN_BONDS + utilityFunctions::integer_power(temp2[BSTATE - 1], 2);
#endif
								in_process_values.setValueFn(vect, MIN(-minVal, 0));
								in_process_values.consumption[h][i][j][l][ll][m] =
									current.current_states[ASTATE]
									- in_process_values.policy_fn[h][i][j][l][ll][m][K1STATE]
									- in_process_values.policy_fn[h][i][j][l][ll][m][K2STATE]
									- in_process_values.policy_fn[h][i][j][l][ll][m][BSTATE];
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
			for (int i = 0; i < NUM_STATE_VARS; i++) {
				if (i > 0) {
					std::cout << ",";
				}
				std::cout << location1[i];
			}
			std::cout << std::endl;
#else
				<< "policy diff=" << diff2 << "       Value diff=" << diff << std::endl;
			std::cout << "PolicyLoc: ";
			for (int i = 0; i < NUM_STATE_VARS; i++) {
				if (i > 0) {
					std::cout << ",";
				}
				std::cout << location2[i];
			}
			std::cout << std::endl;
			std::cout << "ValueLoc: ";
			for (int i = 0; i < NUM_STATE_VARS; i++) {
				if (i > 0) {
					std::cout << ",";
				}
				std::cout << location1[i];
			}
			std::cout << std::endl;
#endif
			std::cout << max_iterations << endl;
		}

		last_values = in_process_values;

		printResults(last_values, phis, whichCountry);
	}

	final = last_values;
	return;
}

int policy_iteration(double difference, State& current, const StochProc& stoch,
	EquilFns& oldValue, EquilFns& latestValue) {

	int max_iterations = MIN(20, floor(log(1.0 / difference)));

	if (max_iterations < 1) {
		return max_iterations;
	}

	EquilFns temp(latestValue);
	EquilFns altPol(latestValue);

	for (int h = 0; h < PHI_STATES; h++) {
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
				for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
					for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
						for (int j = 0; j < ASSET_SIZE; j++) {
							altPol.policy_fn[h][i][j][l][ll][m][K1STATE] = sqrt(
								temp.policy_fn[h][i][j][l][ll][m][K1STATE]
								- MIN_CAPITAL);
							altPol.policy_fn[h][i][j][l][ll][m][K2STATE] = sqrt(
								temp.policy_fn[h][i][j][l][ll][m][K2STATE]
								- MIN_CAPITAL);
							altPol.policy_fn[h][i][j][l][ll][m][BSTATE] = sqrt(
								temp.policy_fn[h][i][j][l][ll][m][BSTATE]
								- MIN_BONDS);
						}
					}
				}
			}
		}
	}

	for (int count = 0; count < max_iterations; count++) {
		vfiMaxUtil ub = vfiMaxUtil(current, stoch, latestValue);
		//eulerUtility vmu = eulerUtility(current,stoch, initial);

		for (int h = 0; h < PHI_STATES; h++) {
			for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
				for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
					for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
						for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
							for (int j = 0; j < ASSET_SIZE; j++) {
								current.current_states[ASTATE] = stoch.assets[j];
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
								current.current_indices[CAP1_SHOCK_STATE] = l;
								current.current_indices[CAP2_SHOCK_STATE] = ll;
								current.current_indices[WAGE_SHOCK_STATE] = m;
								current.current_indices[AGG_SHOCK_STATE] = i;
								current.current_indices[PHI_STATE] = h;

								ub.updateCurrent(current);

#if K2CHOICE
								VecDoub policy = VecDoub(NUM_CHOICE_VARS);
#else
								VecDoub policy = VecDoub(NUM_CHOICE_VARS - 1);
#endif
								policy[K1STATE] = altPol.policy_fn[h][i][j][l][ll][m][K1STATE];
#if K2CHOICE
								policy[K2STATE] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE];
								policy[BSTATE] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE];
#else
								policy[BSTATE - 1] = altPol.policy_fn[h][i][j][l][ll][m][BSTATE];
#endif
								VecInt vect(6);
								vect[0] = h;
								vect[1] = i;
								vect[2] = j;
								vect[3] = l;
								vect[4] = ll;
								vect[5] = m;

								temp.setValueFn(vect, -ub(policy));
							}
						}
					}
				}
			}
		}

		for (int h = 0; h < PHI_STATES; h++) {
			for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
				for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
					for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
						for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
							for (int j = 0; j < ASSET_SIZE; j++) {
								VecInt vect(6);
								vect[0] = h;
								vect[1] = i;
								vect[2] = j;
								vect[3] = l;
								vect[4] = ll;
								vect[5] = m;
								if (count == (max_iterations - 2)) {
									oldValue.setValueFn(vect, temp.getValueFn(vect));
								}
								latestValue.setValueFn(vect, temp.getValueFn(vect));
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

	VecInt vect(6);
	vect[0] = 0;
	vect[1] = 0;
	vect[2] = 0;
	vect[3] = 0;
	vect[4] = 0;
	vect[5] = 0;

	VecDoub val1(1);
	val1[0] = a.getValueFn(vect);
	VecDoub val2(1);
	val2[0] = b.getValueFn(vect);
	double temp_var = 0;
	double diff = utilityFunctions::distanceFn(val1, val2);
	location[ASTATE] = 0;
	location[CAP1_SHOCK_STATE] = 0;
	location[CAP2_SHOCK_STATE] = 0;
	location[WAGE_SHOCK_STATE] = 0;
	location[AGG_SHOCK_STATE] = 0;
	location[PHI_STATE] = 0;
	for (int h = 0; h < PHI_STATES; h++) {
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			for (int j = 0; j < ASSET_SIZE; j++) {
				for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
					for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
						for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
							vect[0] = h;
							vect[1] = i;
							vect[2] = j;
							vect[3] = l;
							vect[4] = ll;
							vect[5] = m;
							val1[0] = a.getValueFn(vect);
							val2[0] = b.getValueFn(vect);
							temp_var = utilityFunctions::distanceFn(val1, val2);
							if (temp_var > diff) {
								diff = temp_var;
								location[ASTATE] = j;
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
	return diff;
}

double maxPolicyDistance(const EquilFns &a, const EquilFns &b,
	VecInt& location) {
	VecDoub val1 = a.policy_fn[0][0][0][0][0][0];
	VecDoub val2 = b.policy_fn[0][0][0][0][0][0];
	double temp_var = 0;
	double diff = utilityFunctions::distanceFn(val1, val2);
	location[ASTATE] = 0;
	location[CAP1_SHOCK_STATE] = 0;
	location[CAP2_SHOCK_STATE] = 0;
	location[WAGE_SHOCK_STATE] = 0;
	location[AGG_SHOCK_STATE] = 0;
	location[PHI_STATE] = 0;
	for (int h = 0; h < PHI_STATES; h++) {
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			for (int j = 0; j < ASSET_SIZE; j++) {
				for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
					for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
						for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
							val1 = a.policy_fn[h][i][j][l][ll][m];
							val2 = b.policy_fn[h][i][j][l][ll][m];
							temp_var = utilityFunctions::distanceFn(val1, val2);
							if (temp_var > diff) {
								diff = temp_var;
								location[ASTATE] = j;
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
	return diff;
}

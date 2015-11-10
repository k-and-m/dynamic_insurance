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
#define USE_MPI 1

#include "BaseHeader.h"
#include "utilityFunctions.h"
#include "matrixIO.h"
#include "matrix.h"
#include "amoeba.h"
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

std::vector<MatrixXd> simulate(int numC, const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target);
Mat3Doub getNextParameters(const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target,
	VecDoub& r_squared);
void printResults(const EquilFns& e, const Mat3Doub& betas, const VecDoub& phis, const COUNTRYID whichCountry);
void readSimpleResults(EquilFns& e, COUNTRYID whichCountry);
void readResults(EquilFns& e, Mat3Doub& betas, COUNTRYID whichCountry);
void initialize(EquilFns& fns, const VecDoub& phis);
void solve(const VecDoub& phis, const VecDoub& prices, const EquilFns& orig, EquilFns& final, const Mat3Doub& recEst, const COUNTRYID whichCountry);
int policy_iteration(double difference, State& initial, const StochProc& stoch,
	EquilFns& lastFns, EquilFns& currentFns);
double maxValueDistance(const EquilFns &a, const EquilFns &b, VecInt& location);
double maxPolicyDistance(const EquilFns &a, const EquilFns &b,
	VecInt& location);
double solveProblem(const VecDoub& phis, const VecDoub& prices, double c1prop, int seqNo);

#if 1

int main(int argc, char *argv[]) {
	using namespace std;
#if USE_MPI
	std::cout << "init" << std::endl << std::flush;
	MPI_Init(&argc, &argv);
	std::cout << "init done" << std::endl << std::flush;
#endif

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
	VecDoub phi = VecDoub(4);
	VecDoub tau = VecDoub(1);
	phi[0] = 0.4;
	phi[1] = 0.8;
	phi[2] = 0;
	phi[3] = 0;

	tau[0] = 1.1;
	double dis = solveProblem(phi, tau, 0.3, 0);
#endif
	std::cout << "soln dist: " << dis << endl;
#if USE_MPI
	MPI_Finalize();
#endif

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
	if (tau.size() != 1) {
		std::cerr << "solveProblem: incorrect number of prices (taus)" << tau.size() << endl;
		exit(-1);
	}

	Mat3Doub recursEst(NUM_RECURSIVE_FNS, PHI_STATES, 3);
	for (int i = 0; i < NUM_RECURSIVE_FNS; i++) {
		for (int j = 0; j < PHI_STATES; j++) {
			switch (i) {
			case P_R:
				recursEst[i][j][0] = 0.01502;
				recursEst[i][j][1] = 0;
				recursEst[i][j][2] = 0;
				break;
			case AGG_ASSET_C1:
				recursEst[i][j][0] = 0;
				recursEst[i][j][1] = 1;
				recursEst[i][j][2] = 0;
				break;
			case AGG_ASSET_C2:
				recursEst[i][j][0] = 0;
				recursEst[i][j][1] = 0;
				recursEst[i][j][2] = 1;
				break;
			default:
				std::cerr << "ERROR! Unknown recursive function: " << i << std::endl;
				exit(-1);
				break;
			}
		}
	}

	double avgDist = 0;
	for (int loopC = 0; (loopC < 200) && (avgDist < 0.99); loopC++) {

		EquilFns orig1, final1;
		EquilFns orig2, final2;

		VecDoub myphis = VecDoub(PHI_STATES);
		for (int i = 0; i < PHI_STATES; i++) {
			myphis[i] = phis[i];
		}

		StochProc stoch1(myphis);
		initialize(orig1, myphis);
		if (utilityFunctions::fileExists("policy_simple_0.dat")) {
			readSimpleResults(orig1, C1);
		}
		if (utilityFunctions::fileExists("policy_0.dat")) {
			readResults(orig1, recursEst, C1);
		}

		std::cout << "Solving Country 1" << std::endl;
		solve(myphis, tau, orig1, final1, recursEst, C1);

		State current1(myphis, tau, recursEst);
		current1.defaultInitialState(stoch1);

		for (int i = PHI_STATES; i < 2 * PHI_STATES; i++) {
			myphis[i - PHI_STATES] = phis[i];
		}
		StochProc stoch2(myphis);
		State current2(myphis, tau, recursEst);
		current2.defaultInitialState(stoch2);
		initialize(orig2, myphis);
		if (utilityFunctions::fileExists("policy_simple_1.dat")) {
			readSimpleResults(orig1, C2);
		}
		if (utilityFunctions::fileExists("policy_1.dat")) {
			readResults(orig2, recursEst, C2);
		}

		std::cout << "Solving Country 2" << std::endl;
		solve(myphis, tau, orig2, final2, recursEst, C2);

		std::cout << "Simulating world economies. " << NUMHHS << " households over " << TOTALPERIODS << " periods." << std::endl;

		VecDoub r_squareds(NUM_RECURSIVE_FNS);
		Mat3Doub xres = getNextParameters(final1, stoch1, current1, final2, stoch2, current2, c1prop, r_squareds);

		avgDist = 0;
		for (int i = 0; i < NUM_RECURSIVE_FNS; i++) {
			avgDist += r_squareds[i];
		}
		avgDist = avgDist / NUM_RECURSIVE_FNS;
		double non_r_R2 = (avgDist * NUM_RECURSIVE_FNS - r_squareds[P_R]) / (NUM_RECURSIVE_FNS - 1);
		bool updateR = non_r_R2 > 0.95;
		if (updateR) {
			std::cout << "Updating R estimate: " << non_r_R2 << std::endl;
		}
		else {
			std::cout << "Skipping R update" << std::endl;
		}
		for (int dim1 = 0; dim1 < xres.dim1(); dim1++) {
			for (int dim2 = 0; dim2 < xres.dim2(); dim2++) {
				for (int dim3 = 0; dim3 < xres.dim3(); dim3++) {
					if (dim1 == P_R) {
						if (updateR) {
							recursEst[dim1][dim2][dim3] = 0.9*recursEst[dim1][dim2][dim3] + 0.1*xres[dim1][dim2][dim3];
						}
					}
					else {
						std::cout << "Updating beta estimates" << std::endl;
						recursEst[dim1][dim2][dim3] = 0.3*recursEst[dim1][dim2][dim3] + 0.7*xres[dim1][dim2][dim3];
					}
				}
			}
		}
		for (int dim1 = 0; dim1 < recursEst.dim1(); dim1++) {
			for (int dim2 = 0; dim2 < recursEst.dim2(); dim2++) {
				for (int dim3 = 0; dim3 < recursEst.dim3(); dim3++) {
					if (dim3 > 0) {
						std::cout << ",";
					}
					std::cout << recursEst[dim1][dim2][dim3];
				}
				std::cout << std::endl;
			}
		}
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

Mat3Doub getNextParameters(const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
	const EquilFns& policies2, const StochProc& stoch2, const State& curSt2, const double c1target,
	VecDoub& r_squared)
{

	if (r_squared.size() != NUM_RECURSIVE_FNS) {
		std::cerr << "akmodel.cc.getNextParameters(): Expect a std::vector of size " << NUM_RECURSIVE_FNS << " for r-squareds. Received std::vector of size " << r_squared.size() << std::endl;
		exit(-1);
	}

	int numC = 2;
	std::vector<MatrixXd> data = simulate(numC, policies1, stoch1, curSt1, policies2, stoch2, curSt2, c1target);

	std::cout << "phi,1,A1,A2,R',A1',A2',NB" << std::endl;
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data[i].rows(); j++) {
			std::cout << i << "," << data[i](j, 0)
				<< "," << data[i](j, 1)
				<< "," << data[i](j, 2)
				<< "," << data[i](j, 3)
				<< "," << data[i](j, 4)
				<< "," << data[i](j, 5)
				<< "," << data[i](j, 6)
				<< std::endl;
		}
	}
	MatrixXd badRHS = data[0].middleCols(0, 4);
	std::vector<VectorXd> badLHS;
	badLHS.resize(NUM_RECURSIVE_FNS);
	badLHS[AGG_ASSET_C1] = data[0].middleCols(4, 1);
	badLHS[AGG_ASSET_C2] = data[0].middleCols(5, 1);
	badLHS[P_R] = data[0].middleCols(6, 1);
	for (int i = 1; i < NUM_RECURSIVE_FNS; i++) {
		for (int j = 0; j < badLHS[i].rows(); j++) {
			badLHS[i][j] = log(badLHS[i][j]);
		}
	}
	for (int i = 0; i < badRHS.cols(); i++) {
		for (int j = 0; j < badRHS.rows(); j++) {
			if (i == 0) {
				if (badRHS(j, i) != 1) {
					std::cerr << "ERROR!: akmodel.cc-getNextParameters(): incorrect constant for badRHS regression. How?";
					badRHS(j, i) = 1;
				}
			}
			else {
				badRHS(j, i) = log(badRHS(j, i));
			}
		}
	}

	MatrixXd goodRHS = data[1].middleCols(0, 4);
	std::vector<VectorXd> goodLHS;
	goodLHS.resize(NUM_RECURSIVE_FNS);
	goodLHS[AGG_ASSET_C1] = data[1].middleCols(4, 1);
	goodLHS[AGG_ASSET_C2] = data[1].middleCols(5, 1);
	goodLHS[P_R] = data[1].middleCols(6, 1);
	for (int i = 1; i < NUM_RECURSIVE_FNS; i++) {
		for (int j = 0; j < goodLHS[i].rows(); j++) {
			goodLHS[i][j] = log(goodLHS[i][j]);
		}
	}
	for (int i = 0; i < goodRHS.cols(); i++) {
		for (int j = 0; j < goodLHS[0].rows(); j++) {
			if (i == 0) {
				if (goodRHS(j, i) != 1) {
					std::cerr << "ERROR!: akmodel.cc-getNextParameters(): incorrect constant for goodRHS regression. How?";
					goodRHS(j, i) = 1;
				}
			}
			else {
				goodRHS(j, i) = log(goodRHS(j, i));
			}
		}
	}

	Mat3Doub results(NUM_RECURSIVE_FNS, PHI_STATES, 3);
	std::vector<VectorXd> bondBetas(2);
	bondBetas[0].resize(4);
	bondBetas[1].resize(4);

	/* Note: This order reflects the order stored in State.h*/
	for (int i = 0; i < NUM_RECURSIVE_FNS; i++) {
		for (int j = 0; j < PHI_STATES; j++) {
			MatrixXd tempRHS;
			VectorXd tempLHS;
			if (j == 0) {
				tempRHS = (i == P_R) ? badRHS : badRHS.middleCols(0, 3);
				tempLHS = badLHS[i];
			}
			else {
				tempRHS = (i == P_R) ? goodRHS : goodRHS.middleCols(0, 3);
				tempLHS = goodLHS[i];
			}
			//			VectorXd temp = (tempRHS.transpose() * tempRHS).ldlt().solve(tempRHS.transpose() * tempLHS);
			//			VectorXd temp = tempRHS.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(tempLHS);
			VectorXd temp = (tempRHS.transpose()*tempRHS).inverse()*tempRHS.transpose()*tempLHS;
			std::cout << "Betas (" << i << "," << j << ")" << temp << std::endl;
			if (i == P_R) {
				bondBetas[j] = temp;
			}
			else {
				for (int k = 0; k < 3; k++) {
					results[i][j][k] = temp[k];
				}
			}
		}
	}

	//NOW CALCULATE New R to get Net Bonds = 0
	VectorXd newGoodRs(goodLHS[0].rows()), newBadRs(badLHS[0].rows());
	if (abs(bondBetas[1][3]) < 0.0001) {
		std::cerr << "ERROR! akmodel.cc-getNextParameter(): goodNBs have almost zero alpha on R: " << bondBetas[1][3] << std::endl;
		exit(-1);
	}
	for (int i = 0; i < newGoodRs.size(); i++) {
		newGoodRs[i] = goodRHS(i, 3) - goodLHS[0][i] / bondBetas[1][3];
	}
	if (abs(bondBetas[0][3]) < 0.0001) {
		std::cerr << "ERROR! akmodel.cc-getNextParameter(): badNBs have almost zero alpha on R: " << bondBetas[0][3] << std::endl;
		exit(-1);
	}
	for (int i = 0; i < newBadRs.size(); i++) {
		newBadRs[i] = badRHS(i, 3) - badLHS[0][i] / bondBetas[0][3];
	}
	goodLHS[P_R] = newGoodRs;
	badLHS[P_R] = newBadRs;

	MatrixXd tempRHS(goodRHS.rows(), goodRHS.cols() - 1);
	for (int i = 0; i < tempRHS.cols(); i++) {
		for (int j = 0; j < tempRHS.rows(); j++) {
			tempRHS(j, i) = goodRHS(j, i);
		}
	}
	goodRHS = tempRHS;

	MatrixXd tempRHS2(badRHS.rows(), badRHS.cols() - 1);
	for (int i = 0; i < tempRHS2.cols(); i++) {
		for (int j = 0; j < tempRHS2.rows(); j++) {
			tempRHS2(j, i) = badRHS(j, i);
		}
	}
	badRHS = tempRHS2;

	//	VectorXd tempBetas =(badRHS.transpose() * badRHS).ldlt().solve(badRHS.transpose() * badLHS[P_R]);
	//	VectorXd tempBetas = badRHS.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(badLHS[P_R]);
	//	VectorXd tempBetas = badRHS.fullPivHouseholderQr().solve(badLHS[P_R]);
	VectorXd tempBetas = (badRHS.transpose()*badRHS).inverse()*badRHS.transpose()*badLHS[P_R];
	for (int i = 0; i < 3; i++) {
		results[P_R][0][i] = tempBetas[i];
	}

	//	tempBetas = (goodRHS.transpose() * goodRHS).ldlt().solve(goodRHS.transpose() * goodLHS[P_R]);
	//	tempBetas = goodRHS.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(goodLHS[P_R]);
	//	tempBetas = goodRHS.fullPivHouseholderQr().solve(goodLHS[P_R]);
	tempBetas = (goodRHS.transpose()*goodRHS).inverse()*goodRHS.transpose()*goodLHS[P_R];
	for (int i = 0; i < 3; i++) {
		results[P_R][1][i] = tempBetas[i];
	}

	/* Calculate the r-squared*/
	VecDoub TSS(NUM_RECURSIVE_FNS);
	VecDoub RSS(NUM_RECURSIVE_FNS);

	for (int h = 0; h < NUM_RECURSIVE_FNS; h++) {

		double avg_y = 0;
		int totalObs = 0;
		for (int hh = 0; hh < PHI_STATES; hh++) {
			VectorXd tempLHS = (hh == 0) ? badLHS[h] : goodLHS[h];
			MatrixXd tempRHS = (hh == 0) ? badRHS : goodRHS;
			int NUMROWS = tempLHS.rows();
			totalObs += NUMROWS;
			for (int i = 0; i < NUMROWS; i++) {
				avg_y += tempLHS[i];
			}
		}
		avg_y = avg_y / totalObs;

		TSS[h] = 0;
		RSS[h] = 0;
		for (int hh = 0; hh < PHI_STATES; hh++) {
			VectorXd tempLHS = (hh == 0) ? badLHS[h] : goodLHS[h];
			MatrixXd tempRHS = (hh == 0) ? badRHS : goodRHS;
			int NUMROWS = tempLHS.rows();
			for (int i = 0; i < NUMROWS; i++) {
				if (tempRHS(i, 0) != 1) {
					std::cerr << "ERROR! akmodel.cc-getNextParameter(): tempRHS[" << i << ",0]!=1. Why not?" << std::endl;
					exit(-1);
				}
				TSS[h] += pow(tempLHS[i] - avg_y, 2);
				double pred = results[h][hh][0] * tempRHS(i, 0) + results[h][hh][1] * tempRHS(i, 1) + results[h][hh][2] * tempRHS(i, 2);
				RSS[h] += pow(tempLHS[i] - pred, 2);
			}
		}
		r_squared[h] = 1 - RSS[h] / TSS[h];

		std::cout << "akmodel.getNextParams(): R-squared(" << h << ")=" << r_squared[h] << std::endl;
	}

	//Update to return the beta estimates
	return results;
}

std::vector<MatrixXd> simulate(int numC, const EquilFns& policies1, const StochProc& stoch1, const State& curSt1,
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
	std::vector<MatrixXd> mydata;
	mydata.resize(PHI_STATES);
	std::vector<VecDoub> hist = we.getHistory();
	int numGoodStates = 0;
	int numBadStates = 0;
	for (int i = 0; i < NUMPERIODS - 1; i++) {
		if (hist[i + SSPERIODS][2] == 1) {
			numGoodStates++;
		}
		else {
			numBadStates++;
		}
	}
	mydata[0] = MatrixXd(numBadStates, 2 * numC + 3);
	mydata[1] = MatrixXd(numGoodStates, 2 * numC + 3);

	int currentGood = 0;
	int currentBad = 0;
	for (int i = 0; i < NUMPERIODS - 1; i++) {
		if (hist[i + SSPERIODS][2] == 0) {
			if (currentBad > numBadStates) {
				std::cerr << "At bad state " << currentBad << " but only expect " << numBadStates;
				exit(-1);
			}
			mydata[0](currentBad, 0) = 1;
			mydata[0](currentBad, 3) = hist[i + SSPERIODS][3];
			for (int j = 0; j < numC; j++) {
				mydata[0](currentBad, j + 1) = hist[i + SSPERIODS][j];
				mydata[0](currentBad, j + 4) = hist[i + SSPERIODS + 1][j];
			}
			mydata[0](currentBad, 6) = hist[i + SSPERIODS][4];
			currentBad++;
		}
		else {
			if (currentGood > numGoodStates) {
				std::cerr << "At good state " << currentGood << " but only expect " << numGoodStates;
				exit(-1);
			}

			mydata[1](currentGood, 0) = 1;
			mydata[1](currentGood, 3) = hist[i + SSPERIODS][3];
			for (int j = 0; j < numC; j++) {
				mydata[1](currentGood, j + 1) = hist[i + SSPERIODS][j];
				mydata[1](currentGood, j + 4) = hist[i + SSPERIODS + 1][j];
			}
			mydata[1](currentGood, 6) = hist[i + SSPERIODS][4];
			currentGood++;
		}
	}

	return mydata;
}

void printResults(const EquilFns& e, const Mat3Doub& betas, const VecDoub& phis, const COUNTRYID whichCountry) {
	using namespace std;
	ofstream out_stream;

	ostringstream os;
	os << "policy_" << whichCountry << ".dat";

	StochProc stoch(phis);
	out_stream.open(os.str());

	for (int dim1 = 0; dim1 < betas.dim1(); dim1++) {
		for (int dim2 = 0; dim2 < betas.dim2(); dim2++) {
			for (int dim3 = 0; dim3 < betas.dim3(); dim3++) {
				if (dim3 > 0) {
					out_stream << ",";
				}
				out_stream << betas[dim1][dim2][dim3];
			}
			out_stream << std::endl;
		}
	}

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
		<< "aggAsset,aggAsset2,phi,agg_shock,z1,z2,wage_shock,a,value_fn,consumption,k1_prime,k2_prime,b_prime"
		<< endl;

	VecInt index(8);

	for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
		for (int f = 0; f < AGG_ASSET_SIZE; f++) {
			for (int g = 0; g < AGG_ASSET_SIZE; g++) {
				for (int j = 0; j < ASSET_SIZE; j++) {
					for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
						for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
							for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
								for (int h = 0; h < PHI_STATES; h++) {
									index[0] = f;
									index[1] = g;
									index[2] = h;
									index[3] = i;
									index[4] = j;
									index[5] = l;
									index[6] = ll;
									index[7] = m;
#if 1
									out_stream
										<< f << ","
										<< g << ","
										<< h << ","
										<< i << ","
										<< l << ","
										<< ll << ","
										<< m << ","
										<< j << ","
										<< e.getValueFn(index) << ","
										<< e.consumption[f][g][h][i][j][l][ll][m] << ","
										<< e.policy_fn[f][g][h][i][j][l][ll][m][K1STATE] << ","
										<< e.policy_fn[f][g][h][i][j][l][ll][m][K2STATE] << ","
										<< e.policy_fn[f][g][h][i][j][l][ll][m][BSTATE]
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
		}
	}
	out_stream.close();
}

void readSimpleResults(EquilFns& e, COUNTRYID whichCountry) {
	using namespace std;

	FILE *fp = NULL;
	if (whichCountry == C1) {
		fp = std::fopen("policy_simple_0.dat", "r");
	}
	else {
		fp = std::fopen("policy_simple_1.dat", "r");
	}

	{
		double temp = 0;
		if (fscanf(fp, "%lf", &temp) == 1) {
		}
		else {
			std::cerr << " ERROR! akmodel.cc-readSimpleResults(): could not read expected phi " << std::endl;
			std::fclose(fp);
			exit(-1);
		}
	}
	{
		double temp1 = 0, temp2 = 0;
		if (fscanf(fp, "%lf, %lf", &temp1, &temp2) == 2) {
		}
		else {
			std::cerr << " ERROR! akmodel.cc-readSimpleResults(): could not read capital shocks " << std::endl;
			std::fclose(fp);
			exit(-1);
		}
	}
	{
		double temp1 = 0, temp2 = 0;
		if (fscanf(fp, "%lf, %lf", &temp1, &temp2) == 2) {
		}
		else {
			std::cerr << " ERROR! akmodel.cc-readSimpleResults(): could not read wage shocks " << std::endl;
			std::fclose(fp);
			exit(-1);
		}
	}

	int phi, aggShock, z1, z2, wage, a;
	double valuefn, cons, k1p, k2p, bp;

	if (1) {
		char str[256];
		fscanf(fp, "%s", str);
	}
	while (fscanf(fp, "%d, %d, %d, %d, %d, %d, %lf, %lf, %lf, %lf, %lf",
		&phi, &aggShock, &z1, &z2, &wage, &a, &valuefn, &cons, &k1p, &k2p, &bp) == 11) {

		VecInt index(8);
		for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
			for (int j = 0; j < AGG_SHOCK_SIZE; j++) {
				for (int k = 0; k < PHI_STATES; k++) {
					index[0] = i;
					index[1] = j;
					index[2] = k;
					index[3] = aggShock;
					index[4] = a;
					index[5] = z1;
					index[6] = z2;
					index[7] = wage;

					e.setValueFn(index, valuefn);
					e.consumption[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]] = cons;
					e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][K1STATE] = k1p;
					e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][K2STATE] = k2p;
					e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][BSTATE] = bp;
				}
			}
		}
	}
	fclose(fp);
}

void readResults(EquilFns& e, Mat3Doub& betas, COUNTRYID whichCountry) {
	using namespace std;

	FILE *fp = NULL;
	if (whichCountry == C1) {
		fp = std::fopen("policy_0.dat", "r");
	}
	else {
		fp = std::fopen("policy_1.dat", "r");
	}

	for (int dim1 = 0; dim1 < betas.dim1(); dim1++) {
		for (int dim2 = 0; dim2 < betas.dim2(); dim2++) {
			double val1, val2, val3;
			if (fscanf(fp, "%lf,%lf,%lf", &val1, &val2, &val3) == 3) {
				betas[dim1][dim2][0] = val1;
				betas[dim1][dim2][1] = val2;
				betas[dim1][dim2][2] = val3;
			}
			else {
				std::cerr << " ERROR! akmodel.cc-readResults(): could not read beta values for dim1=" << dim1 << " dim2=" << dim2 << std::endl;
				std::fclose(fp);
				exit(-1);
			}
		}
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
	while (fscanf(fp, "%d, %d, %d, %d, %d, %d, %d, %d, %lf, %lf, %lf, %lf, %lf",
		&a1, &a2, &phi, &aggShock, &z1, &z2, &wage, &a, &valuefn, &cons, &k1p, &k2p, &bp) == 13) {

		VecInt index(8);

		index[0] = a1;
		index[1] = a2;
		index[2] = phi;
		index[3] = aggShock;
		index[4] = a;
		index[5] = z1;
		index[6] = z2;
		index[7] = wage;

		e.setValueFn(index, valuefn);
		e.consumption[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]] = cons;
		e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][K1STATE] = k1p;
		e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][K2STATE] = k2p;
		e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][BSTATE] = bp;

		if (doublePhi) {
			index[2] = 1;
			e.setValueFn(index, valuefn);
			e.consumption[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]] = cons;
			e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][K1STATE] = k1p;
			e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][K2STATE] = k2p;
			e.policy_fn[index[0]][index[1]][index[2]][index[3]][index[4]][index[5]][index[6]][index[7]][BSTATE] = bp;
		}
	}
	fclose(fp);
}


void initialize(EquilFns& fns, const VecDoub& phis) {

	using namespace std;

	StochProc stoch(phis);
	VecInt vect(8);
	for (int f = 0; f < AGG_ASSET_SIZE; f++) {
		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int h = 0; h < PHI_STATES; h++) {
				for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
					for (int j = 0; j < ASSET_SIZE; j++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									vect[0] = f;
									vect[1] = g;
									vect[2] = h;
									vect[3] = i;
									vect[4] = j;
									vect[5] = l;
									vect[6] = ll;
									vect[7] = m;
									fns.setValueFn(vect, vfiMaxUtil::consUtil(stoch.assets[j] - MIN_ASSETS + 0.1));
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

void solve(const VecDoub& phis, const VecDoub& prices, const EquilFns& orig, EquilFns& final, const Mat3Doub& recEst, const COUNTRYID whichCountry) {
	using namespace std;

#if USE_MPI
	int numprocs, rank, namelen;
	static int firstCall = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	std::cout << "numprocs:" << numprocs << std::endl << std::flush;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::cout << "rank:" << rank << std::endl << std::flush;
#endif

	EquilFns last_values;
	StochProc stoch(phis);

	int counter = 0;
	double diff = 100;
	double diff2 = 100;

	last_values = orig;

	while (counter < MAX_ITER && diff > VAL_TOL) {

		counter++;
		EquilFns in_process_values;

#if USE_MPI
		int mystart = (AGG_ASSET_SIZE / numprocs) * rank;
		int myend;
		if (AGG_ASSET_SIZE % numprocs > rank) {
			mystart += rank;
			myend = mystart + (AGG_ASSET_SIZE / numprocs) + 1;
		}
		else {
			mystart += AGG_ASSET_SIZE % numprocs;
			myend = mystart + (AGG_ASSET_SIZE / numprocs);
		}


		int SIZE = AGG_ASSET_SIZE*AGG_ASSET_SIZE*PHI_STATES*AGG_SHOCK_SIZE*ASSET_SIZE*CAP_SHOCK_SIZE*CAP_SHOCK_SIZE*WAGE_SHOCK_SIZE*NUM_CHOICE_VARS;
		int SIZE2 = AGG_ASSET_SIZE*AGG_ASSET_SIZE*PHI_STATES*AGG_SHOCK_SIZE*ASSET_SIZE*CAP_SHOCK_SIZE*CAP_SHOCK_SIZE*WAGE_SHOCK_SIZE;
		double *policyArray = new double[SIZE];
		double *valueArray = new double[SIZE2];
		if (rank == 0) {
			last_values.policyToArray(policyArray);
			last_values.valueToArray(valueArray);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(policyArray, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(valueArray, SIZE2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (rank != 0) {
			std::cout << rank << " set" << std::endl << std::flush;
			last_values.setPolicyFromArray(policyArray);
			last_values.setValueFromArray(valueArray);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int g = mystart; g < myend; g++) {
			std::cout << rank << ":" << g << std::endl << std::flush;
#if OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(8)
#endif
			for (int j = 0; j < ASSET_SIZE; j++) {
#else
#if OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(8)
#endif
		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int j = 0; j < ASSET_SIZE; j++) {
#endif
				for (int gg = 0; gg < AGG_ASSET_SIZE; gg++) {
					for (int h = 0; h < PHI_STATES; h++) {
						for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
							State current(phis, prices, recEst);
							vfiMaxUtil ub = vfiMaxUtil(current, stoch, last_values);
#if AMOEBA
							Amoeba am(MINIMIZATION_TOL);
#endif
							for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
								for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
									for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
										VecInt vect(8);
										vect[0] = g;
										vect[1] = gg;
										vect[2] = h;
										vect[3] = i;
										vect[4] = j;
										vect[5] = l;
										vect[6] = ll;
										vect[7] = m;

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

#if K2CHOICE
										VecDoub temp1(NUM_CHOICE_VARS);
#else
										VecDoub temp1(NUM_CHOICE_VARS - 1);
#endif
										temp1[K1STATE] =
											sqrt(
												MAX(MIN_CAPITAL, last_values.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE])
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
												MAX(MIN_BONDS, last_values.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE])
												- MIN_BONDS);
										VecDoub_I initial_point(temp1);
										VecDoub temp2;
										double minVal = 0;
										if (l == 0 && ll == 0) {
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
												temp2 = am.minimize(initial_point, delta, ub);
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
										}
										else {
#if K2CHOICE
											temp2 = VecDoub(NUM_CHOICE_VARS);
#else
											temp2 = VecDoub(NUM_CHOICE_VARS - 1);
#endif
											temp2[K1STATE] =
												sqrt(
													in_process_values.policy_fn[g][gg][h][i][j][0][0][m][K1STATE]
													- MIN_CAPITAL);
#if K2CHOICE
											temp2[K2STATE] =
												sqrt(
													in_process_values.policy_fn[g][gg][h][i][j][0][0][m][K2STATE]
													- MIN_CAPITAL);
											temp2[BSTATE] =
#else
											temp2[BSTATE - 1] =
#endif
												sqrt(
													in_process_values.policy_fn[g][gg][h][i][j][0][0][m][BSTATE]
													- MIN_BONDS);

											minVal =
												-in_process_values.getValueFn(vect);
										}

										in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE] =
											MIN_CAPITAL + utilityFunctions::integer_power(temp2[K1STATE], 2);
										in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE] =
#if K2CHOICE
											MIN_CAPITAL + utilityFunctions::integer_power(temp2[K2STATE], 2);
#else
											in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE]
											* pow(prices[0], 1 / (1 - ALPHA1));
#endif
										in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE] =
#if K2CHOICE
											MIN_BONDS + utilityFunctions::integer_power(temp2[BSTATE], 2);
#else
											MIN_BONDS + utilityFunctions::integer_power(temp2[BSTATE - 1], 2);
#endif
										in_process_values.setValueFn(vect, MIN(-minVal, 0));
										in_process_values.consumption[g][gg][h][i][j][l][ll][m] =
											current.current_states[ASTATE]
											- in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE]
											- in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE]
											- in_process_values.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE];
									}
								}
							}
						}
					}
				}
			}
		}
#if USE_MPI
		in_process_values.policyToArray(policyArray);
		in_process_values.valueToArray(valueArray);
		double *combinedPolicy = new double[SIZE];
		double *combinedValue = new double[SIZE2];
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(policyArray, combinedPolicy, SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		delete[] policyArray;
		policyArray = NULL;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(valueArray, combinedValue, SIZE2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		delete[] valueArray;
		valueArray = NULL;
		if (rank == 0) {
			in_process_values.setPolicyFromArray(combinedPolicy);
			in_process_values.setValueFromArray(combinedValue);
#else
		{
#endif
			for (int f = 0; f < AGG_ASSET_SIZE; f++) {
				for (int g = 0; g < AGG_ASSET_SIZE; g++) {
					for (int h = 0; h < PHI_STATES; h++) {
						for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
							for (int ii = 0; ii < ASSET_SIZE; ii++) {
								for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
									for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
										for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
											in_process_values.consumption[f][g][h][i][ii][l][ll][m] =
												stoch.assets[ii]
												- in_process_values.policy_fn[f][g][h][i][ii][l][ll][m][K1STATE]
												- in_process_values.policy_fn[f][g][h][i][ii][l][ll][m][K2STATE]
												- in_process_values.policy_fn[f][g][h][i][ii][l][ll][m][BSTATE];
										}
									}
								}
							}
						}
					}
				}
			}

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
			State current(phis, prices, recEst);
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
			std::cout << max_iterations << endl;
		}

		last_values = in_process_values;
		printResults(last_values, recEst, phis, whichCountry);

		delete[] combinedPolicy;
		delete[] combinedValue;
		combinedPolicy = NULL;
		combinedValue = NULL;

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

	for (int g = 0; g < AGG_ASSET_SIZE; g++) {
		for (int gg = 0; gg < AGG_ASSET_SIZE; gg++) {
			for (int h = 0; h < PHI_STATES; h++) {
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
								}
							}
						}
					}
				}
			}
		}
	}

	for (int count = 0; count < max_iterations; count++) {
		vfiMaxUtil ub = vfiMaxUtil(current, stoch, latestValue);
		//eulerUtility vmu = eulerUtility(current,stoch, initial);

		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int gg = 0; gg < AGG_ASSET_SIZE; gg++) {
				for (int h = 0; h < PHI_STATES; h++) {
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

#if K2CHOICE
										VecDoub policy = VecDoub(NUM_CHOICE_VARS);
#else
										VecDoub policy = VecDoub(NUM_CHOICE_VARS - 1);
#endif
										policy[K1STATE] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][K1STATE];
#if K2CHOICE
										policy[K2STATE] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][K2STATE];
										policy[BSTATE] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE];
#else
										policy[BSTATE - 1] = altPol.policy_fn[g][gg][h][i][j][l][ll][m][BSTATE];
#endif
										VecInt vect(8);
										vect[0] = g;
										vect[1] = gg;
										vect[2] = h;
										vect[3] = i;
										vect[4] = j;
										vect[5] = l;
										vect[6] = ll;
										vect[7] = m;

										if (l == 0 && ll == 0) {
											temp.setValueFn(vect, -ub(policy));
										}
										else {
											VecInt vect2 = vect;
											vect2[5] = 0;
											vect2[6] = 0;
											temp.setValueFn(vect, temp.getValueFn(vect2));
										}
									}
								}
							}
						}
					}
				}
			}
		}

		for (int g = 0; g < AGG_ASSET_SIZE; g++) {
			for (int gg = 0; gg < AGG_ASSET_SIZE; gg++) {
				for (int h = 0; h < PHI_STATES; h++) {
					for (int i = 0; i < AGG_SHOCK_SIZE; i++) {
						for (int l = 0; l < CAP_SHOCK_SIZE; l++) {
							for (int ll = 0; ll < CAP_SHOCK_SIZE; ll++) {
								for (int m = 0; m < WAGE_SHOCK_SIZE; m++) {
									for (int j = 0; j < ASSET_SIZE; j++) {
										VecInt vect(8);
										vect[0] = g;
										vect[1] = gg;
										vect[2] = h;
										vect[3] = i;
										vect[4] = j;
										vect[5] = l;
										vect[6] = ll;
										vect[7] = m;
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
		}
	}
	return max_iterations;
}


double maxValueDistance(const EquilFns &a, const EquilFns &b,
	VecInt& location) {

	VecInt vect(8);
	vect[0] = 0;
	vect[1] = 0;
	vect[2] = 0;
	vect[3] = 0;
	vect[4] = 0;
	vect[5] = 0;
	vect[6] = 0;
	vect[7] = 0;

	VecDoub val1(1);
	val1[0] = a.getValueFn(vect);
	VecDoub val2(1);
	val2[0] = b.getValueFn(vect);
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
									vect[0] = g;
									vect[1] = gg;
									vect[2] = h;
									vect[3] = i;
									vect[4] = j;
									vect[5] = l;
									vect[6] = ll;
									vect[7] = m;
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

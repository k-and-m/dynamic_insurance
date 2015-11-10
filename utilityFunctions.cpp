#include "utilityFunctions.h"

namespace utilityFunctions{

	double integer_power(const double base, const int exponent){

		unsigned int exp = abs(exponent);
		if (base == 0 && exponent < 0){
			std::cerr << "base = 0, exponent is negative. This is not a valid combination." << std::endl;
			exit(-1);
		}
		double newBase = (exponent < 0) ? 1 / base : base;
		double retval = 1;

		switch (exp){

		case 0:
			return 1;
		case 1:
			return newBase;
		case 2:
			return newBase*newBase;
		default:
			break;
		}

		while (exp > 0){
			newBase = newBase * newBase;
			if (exp % 2 == 0){
				exp = exp / 2;
			}
			else{
				exp = (exp - 1) / 2;
				retval = retval*newBase;
			}
		}
		return retval;
	}

	double distanceFn(VecDoub& x, VecDoub& y) {
		double retVal = 0;
		for (unsigned int i = 0; i < x.size(); i++) {
			retVal += integer_power(x[i] - y[i], 2);
		}
		retVal = sqrt(retVal);
		return retVal;
	}

	double interpolate(const VecDoub &x, const VecDoub::const_iterator y, double a_prime) {
		using namespace std;

		int right = 0;
		int left = 0;
		double slope;
		double distance;
		double interp_val;
		int size = x.size();

		if (size == 1) {
			std::cerr << "UtilityFunctions.cpp-interpolate():Should never pass a size 1 vector to interpolate" << std::endl << std::flush;
			exit(-1);
		}

		if (a_prime <= x[0]) {
			right = 1;
			left = 0;
		}
		else if (a_prime >= x[size-1]) {
			right = size - 1;
			left = right - 1;
		}
		else {
			int i;
			for (i = 0; i < size; i++) {
				if (x[i] >= a_prime) {
					right = i;
					left = i - 1;
					break;
				}
			}
			if (i == size) {
				cout << "interpolate: Error! Couldn't find aprime= " << a_prime << endl;
				for (i = 0; i < size; i++) {
					cout << x[i] << ",";
				}
				exit(-1);
			}
		}
		slope = (y[right] - y[left]) / (x[right] - x[left]);
		distance = a_prime - x[left];
		interp_val = y[left] + distance * slope;

		if (interp_val != interp_val) {
			std::cout << "ERROR! utilityFunctions-interpolate() - return value is NaN. "
				<< y[right] << "-" << y[left] << "/(" << x[right] << "-"
				<< x[left] << "), distance=" << distance << " slope=" << slope
				<< " a_prime=" << a_prime << std::endl << std::flush;
			std::cout << "left: " << left << "  right: " << right << "size: " << size << std::endl;
			for (int i = 0; i < size; i++) {
				std::cout << "x[" << i << "]: " << x[i] << std::endl;
			}
			exit(-1);
		}

		return interp_val;
	}

	double interpolate2d(const VecDoub &x1, const VecDoub &x2, const vector<VecDoub>::const_iterator y, double x1_prime, double x2_prime) {
		using namespace std;

		if (AGG_ASSET_SIZE == 1){
			return interpolate(x2, y[0].begin(), x2_prime);
		}

		VecDoub int1(AGG_ASSET_SIZE);

		for (int i = 0; i < AGG_ASSET_SIZE; i++){
			int1[i] = interpolate(x2, y[i].begin(), x2_prime);
		}

		double interp_val = interpolate(x1, int1.begin(), x1_prime);
		return interp_val;
	}

	double interpolate3d(const VecDoub &x1, const VecDoub &x2, const VecDoub &x3, const vector<vector<VecDoub>>::const_iterator y, double x1_prime, double x2_prime, double x3_prime){
		using namespace std;

		if (AGG_ASSET_SIZE == 1){
			return interpolate2d(x2, x3, y[0].begin(), x2_prime, x3_prime);
		}

		VecDoub int1(AGG_ASSET_SIZE);

		for (int i = 0; i < AGG_ASSET_SIZE; i++){
			int1[i] = interpolate2d(x2, x3, y[i].begin(), x2_prime, x3_prime);
		}

		double interp_val = interpolate(x1, int1.begin(), x1_prime);
		return interp_val;
	}

	std::string tostr(double t) {
		std::ostringstream os;
		os << t;
		return os.str();
	}

	double boundValue(const double origVal, const double lower, const double upper){
		if (MAX(MIN(origVal, upper), lower) == origVal){
			return origVal;
		}

		double newVal = origVal;

		if (origVal < lower){
			newVal = 2 * lower - origVal;
		}

		double distFromLower = newVal - lower;
		double maxDist = upper - lower;

		int multiple = floor(distFromLower / maxDist);
		if (multiple % 2 == 1){
			newVal = upper - (distFromLower - multiple*maxDist);
		}
		else{
			newVal = lower + (distFromLower - multiple*maxDist);
		}
		return newVal;
	}

	// Function: fileExists
	/**
	Check if a file exists
	@param[in] filename - the name of the file to check

	@return    true if the file exists, else false

	*/
	bool fileExists(const std::string& filename)
	{
		struct stat buf;
		if (stat(filename.c_str(), &buf) != -1)
		{
			return true;
		}
		return false;
	}
}
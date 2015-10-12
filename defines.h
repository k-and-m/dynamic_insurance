#ifndef DEFINES_H_
#define DEFINES_H_
#define OPENMP 1
#define ASSETSTATE 1
#define IID_CAP 1
#define WORLD_ECON_SEED 1234

#define MIN(a,b) ((a<b)?a:b)
#define MAX(a,b) ((a>b)?a:b)

#define RRA 2

#define AGG_SHOCK_SIZE 1

#define AGG_ASSET_SIZE 5
#define MAX_AGG_ASSETS 10
#define MIN_AGG_ASSETS 1 

#define PHI_STATES 2

#if ASSETSTATE
#define ASSET_SIZE 150
#else
#define BOND_SIZE 30
#define CAPITAL_SIZE 30
#endif

#define CAP_SHOCK_SIZE	2
#define WAGE_SHOCK_SIZE 2

#define CAP1_SHOCK_MEAN  0.15
#define CAP1_SHOCK_STD .1
#define CAP1_SHOCK_PERSISTENCE 0
#define CAP1_SHOCK_SPREAD 3 //Number of standard deviations per side for Tauchen

#define CAP2_SHOCK_MEAN  0.15
#define CAP2_SHOCK_STD   .1
#define CAP2_SHOCK_PERSISTENCE 0
#define CAP2_SHOCK_SPREAD 3 //Number of standard deviations per side for Tauchen


#if ASSETSTATE
#define MAX_ASSETS 40
#define MIN_ASSETS 0 

//#define MIN_MGMT 0
//#define MAX_MGMT 100
#endif

#define MIN_CAPITAL 0
#define MIN_BONDS -MAX_ASSETS 

#define WAGE_SHOCK_MEAN 0.85
#define WAGE_SHOCK_STD .2 //Assuming log wages
#define WAGE_SHOCK_PERSISTENCE 0.9
#define WAGE_SHOCK_SPREAD 3
#define CURVATURE 3
#define AGG_CURVATURE 1
#define MINIMIZATION_TOL pow(10,-8)
#define VAL_TOL 3*pow(10,-4)
#define THETA 1
#define ALPHA1 0.25
#define BETA 0.925
#define MAX_ITER 100

//Definitions for state.current_states, which is a VecDoub of five dimensions.
//Also, first two dimensions work for EquilFns.policy_fn.
#define ASTATE 0

#define K1STATE 0 //Capital is first state index
#define K2STATE K1STATE+1 //capital in country 2 is second state index
#define BSTATE K2STATE+1 //Bonds are second state index
//#define MGMT_C1_STATE BSTATE+1

#define NUM_STATE_CHOICE_VARS BSTATE+1
#define NUM_PERIOD_CHOICE_VARS 0
#define NUM_CHOICE_VARS NUM_STATE_CHOICE_VARS+NUM_PERIOD_CHOICE_VARS

#define CAP1_SHOCK_STATE ASTATE+1 //Shocks to capital are the third state index
#define CAP2_SHOCK_STATE CAP1_SHOCK_STATE+1 //Shocks to capital are the third state index
#define WAGE_SHOCK_STATE CAP2_SHOCK_STATE+1// Shocks to wages are the fourth state index
#define AGG_SHOCK_STATE WAGE_SHOCK_STATE+1 //Aggregate shocks are the fifth state index
#define PHI_STATE AGG_SHOCK_STATE + 1
#define AGG_ASSET_STATE PHI_STATE+1
#define AGG2_ASSET_STATE PHI_STATE+1
#define NUM_STATE_VARS AGG2_ASSET_STATE+1

//Definitions for within "StochProc.shocks"
#define EF_K1 0
#define EF_K2 1
#define EF_W 2
#define EF_A 3
#define EF_PHI 4
#define NUM_SHOCK_VARS 5


//definitions for order of prices
#define P_R 0
#define AGG_ASSET_C1 1
#define AGG_ASSET_C2 2
#define NUM_RECURSIVE_FNS AGG_ASSET_C2+1
#endif

//definitions for Krusell-Smith economy simulation
#define SSPERIODS 100
#define NUMPERIODS 100
#define TOTALPERIODS SSPERIODS+NUMPERIODS
#define NUMHHS 10000

//random price
#define NOMINAL_PRICE 3.0
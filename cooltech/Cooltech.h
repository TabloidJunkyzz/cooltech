// "Cooltech.h": Header for "Cooltech.cpp" and "Cooltech_main.cpp"

//Header Guard as per usual
#ifndef COOLTECH_H
#define COOLTECH_H

#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
//#define NDEBUG
#include <cassert>
#define M_PI   3.14159265358979323846


//constants for length of arrays and Epsilon !
const int NPIPE = 6;
const int NPUMP = 9;
const int NCOST = 8;
const double EPSILON = 1e-6;
const double TOTALFLOW = 0.432; //[m³/s]



//Creates a struct for all of our Variables
struct data{
	//data pipe
	double length[NPIPE];
	double flow[NPIPE];
	double deltaP[NPIPE];
	double k;

	//data cost
	double insideDiameter[NCOST];
	double wallThickness[NCOST];
	double outsideDiameter[NCOST];
	int    dIN[NCOST];
	double operationTime;

	//data water
	double density;
	double viscosity;

	//data pump
	double powerMech[NPUMP];
	double powerEL[NPUMP];
	double efficiency[NPUMP];

	//data result
	double pressureLossFriction[NPUMP][NPIPE];
	double diameter[NPUMP][NPIPE];
	double reynolds[NPUMP][NPIPE];
	double insideDiameterChart[NPUMP][NPIPE];
	double outsideDiameterChart[NPUMP][NPIPE];
	double pressureLossValve[NPUMP][NPIPE];
	double pumpCost[NPUMP];
	double pipeCost[NPUMP][NPIPE];
	double powerCost[NPUMP];
	double totalCost[NPUMP];


};

//Extern Variables for Header of read and output
extern std::ifstream file;
extern double gelesen;
extern std::string skip;

extern std::ofstream out;
extern double temp;

//String for Versions
extern std::string  version;

//Extern Variables for Header of functions in pipe measurements and main
extern double lambda;
extern double lambdaStart;
extern double lambdaTemp;

//Inwiefern brauch man das noch ??
extern double newton;
extern double newtonHelp;

//temporarily stores Pipe costs for each Pump
extern double pipeCostStorage;

//declarations for polynomfit
extern int degree;
extern double* coefficents;

//Struct Data declarations!
extern data pipe;
extern data water;
extern data pump;
extern data cost;
extern data result;

//Function Prototypes!

//READ
void read();
void printRead();
void printResults();
void output();

//PIPE MEASUREMENTS
double pressureLossFrictionFunction (int i,int j);
double nikuradse(double k);
double diameterFunction (int i,int j);
double reynoldsFunction (int i,int j);
double pressureLossValveFunction(int i,int j);


//NEWTON
double colebrookFunction (int i, int j, double lambda);
double colebrookDerivativeFunction (int i, int j, std::function<double(int, int, double)>func);
double newtonFunction (int i, int j, double lambda, std::function<double(int,int,double)>func);

//Chart Diameter
void pickDiameterFromChart ();

//Information for user  --> Still needs to be written
void userInformation();

//Costfunctions

double pipeCostFunction(int i,int j);
double pumpCostFunction(int i);
double powerCostFunction(int i);
double totalCostFunction(int i);

void polynomialFit();

#endif //Header Guard END

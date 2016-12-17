// "Cooltech.h": Header for "Cooltech.cpp" and "Cooltech_main.cpp"

//Header Guard as per usual
#ifndef COOLTECH_H
#define COOLTECH_H

#include <fstream>
#include <iostream>
#include <string>
#include <functional>


//constants for length of arrays and Epsilon !
const int NPIPE = 6;
const int NPUMP = 9;
const int NCOST = 8;
const double EPSILON = 1e-6;



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
	double pipeCost[NPIPE];
	double powerCost[NCOST];
	double totalCost[NPUMP][NPIPE];


};

//Extern Variables for Header of read and output
extern std::ifstream file;
extern double gelesen;
extern std::string skip;

extern std::ofstream out;
extern double temp;

//Extern Variables for Header of functions in pipe measurements and main
extern double lambda;
extern double lambdaStart;
extern double lambdaTemp;

extern double newton;
extern double newtonHelp;

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

//Information for user
void userInformation();

//Costfunctions

double pipeCostFunction(int i,int j);
double pumpCostFunction(int i);
double powerCostFunction(int i);
double totalCostFunction(int i,int j);


#endif //Header Guard END

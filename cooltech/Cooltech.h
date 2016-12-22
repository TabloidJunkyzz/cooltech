
//++++++++++++++++++++++ "Cooltech.h": Header for "Cooltech.cpp" and "Cooltech_main.cpp" ++++++++++++++++++++++

//Header Guard as per usual
#ifndef COOLTECH_H
#define COOLTECH_H

//########################################### Includes ###########################################

#include <fstream>    //Ifstream and Ofstream to read and output Data in .txt
#include <functional> //Function pointers for general Newton function
#include <iostream>   //cin, cout, cerr etc.
#include <string>     //To Create strings
#include <math.h>     //pow(), fabs(), sqrt() --> Trivial Mathematical functions
#include <iomanip>    //For setw() To make the console and .txt look nice
//#define NDEBUG      //Ignores assertions !But not if commented!
#include <cassert>    //Assertions
#define M_PI   3.14159265358979323846

//########################################### Constants ###########################################

//constants for length of arrays & Epsilon & Total flow V(point) !
const int NPIPE = 6;
const int NPUMP = 9;
const int NCOST = 8;
const double EPSILON = 1e-6;
const double TOTALFLOW = 0.432; //[in m^3/s]

//########################################### Variables and Structs ###########################################

//Creates a struct for all of our Variables
struct data{
	//data pipe
	double length[NPIPE];
	double flow  [NPIPE];
	double deltaP[NPIPE];
	double k;

	//data cost
	double insideDiameter [NCOST];
	double wallThickness  [NCOST];
	double outsideDiameter[NCOST];
	int    dIN            [NCOST];
	double operationTime;

	//data water
	double density;
	double viscosity;

	//data pump
	double powerMech [NPUMP];
	double powerEL   [NPUMP];
	double efficiency[NPUMP];

	//data result
	double pressureLossFriction[NPUMP][NPIPE];
	double diameter            [NPUMP][NPIPE];
	double reynolds            [NPUMP][NPIPE];
	double insideDiameterChart [NPUMP][NPIPE];
	double outsideDiameterChart[NPUMP][NPIPE];
	double pressureLossValve   [NPUMP][NPIPE];
	double pumpCost            [NPUMP];
	double pipeCost            [NPUMP][NPIPE];
	double powerCost           [NPUMP];
	double totalCost           [NPUMP];


};

//Extern Variables for Header of read and output
extern std::ifstream file; //read --> Creates file operator to read things from a .txt file
extern double gelesen;     //read --> Values temporary variable that gets overwritten constantly
extern std::string skip;   //read --> Gets written but is never used ... only to skip last "words" in .txt files

extern std::ofstream out;  //output --> Creates file operator to write stuff in a .txt file

//String for Versions
extern std::string  version; //Important for switching between Dev and User programm

//Extern Variables for Header of functions in pipe measurements and main
extern double lambda;       //lambda which will be calculated through newton
extern double lambdaStart;  //Start value for lambda with nikuradse
extern double lambdaTemp;   //Helping variable to store lambda (while-loop of main)

//Newton Variables
extern double newton;       //"result" of newton --> lambda gets set to newton in newton function
extern double newtonHelp;   //Helping variable to store lambda (in newton-function)

//temporarily stores Pipe costs for each Pump
extern double pipeCostStorage;

//declarations for polynomfit
extern int degree;
extern double* coefficents;

//Minimum Calculation Variables
extern double xMin;                 //x-value of our global minimum
extern double minCost;              //both xmin and minCost compose the coordinates of our global minimum

//Mistake bool Variable
extern bool mistake;    //Variable to validate mistake with pump 1

//Struct Data declarations!
extern data pipe;
extern data water;
extern data pump;
extern data cost;
extern data result;

//########################################### Function Prototypes ###########################################

//READ
void read();
void printRead();
void printResults();
void output();

//Mistake
void checkForMistake(int i,int j);

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
double totalCostFunction(int i);

//Polynomfit
void polynomialFit();

//Calculation of Minimum
void minCostFunction();

#endif //Header Guard END

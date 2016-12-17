// "Cooltech.cpp" Implementation for our header "Cooltech.h"

#include "Cooltech.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <math.h>
//#define NDEBUG
#include <cassert>
#define M_PI   3.14159265358979323846
using namespace std;


//Initialisation of our variables again to make them local

//File Operators and variables for read and output
ifstream file;
string skip;
double gelesen;

ofstream out;
double temp;


//Variables for functions of pipe measurements and main
double lambda;
double lambdaStart;
double lambdaTemp;

double newton;
double newtonHelp;

//Struct variables
data pipe, water, pump, cost, result;

//Declaration of our "read"-function

void read()//Reads in Data from charts
{

	//Setup for Reading Pipe Data(Creating file operator)
file.open("Rohrdaten.txt");

	//Reads Flow of each pipe [First line of "Rohrdaten.txt"]; Skips end of line
	for(int i=0;i<NPIPE;i++)
	{
		file>>gelesen;
		pipe.flow[i] = 0.001*gelesen;
				assert(pipe.flow[i] > 0);
	}
	file >> skip;

	//Reads Pressureloss of each pipe [Second Line of "Rohrdaten.txt"]; Skips end of line
	for(int i=0;i<NPIPE;i++)
	{
		file>>gelesen;
		pipe.deltaP[i] = 1e5*gelesen;
				assert(pipe.deltaP[i] > 0);

	}
	file >> skip;

	//Reads Pipelength of each pipe [Third line of "Rohrdaten.txt"]; Skips end of line
	for(int i=0;i<NPIPE;i++)
	{
		file>>gelesen;
		pipe.length[i] = gelesen;
				assert(pipe.length[i] > 0);
	}
	file >> skip;

	#if 0
	Reads k [Fourth Line of "Rohrdaten.txt"]; Skips end of line
	Then closes file, cause every pipedata is in the system!
	#endif
	file>>gelesen;
	pipe.k = gelesen/1000;
		assert(pipe.k > 0);
file.close();


	//Setup for Reading water data
file.open("Wasserdaten.txt");

	//Reads Water density [First line of "Wasserdaten.txt"]; Skips end of line
	file >> gelesen;
	water.density = gelesen;
		assert(water.density > 0);
	file >> skip;

	//Reads Water viscosity [Second line of "Wasserdaten.txt"]; Skips end of line; Closes file
	file >> gelesen;
	water.viscosity = gelesen;
			assert(water.viscosity > 0);
file.close();


	//Setup for reading pump data
file.open("Pumpendaten.txt");

	//Reads Pumppowers [First line of "Pumpendaten.txt"]; Skips end of line
	for(int i=0;i<NPUMP;i++)
	{
		file>>gelesen;
		pump.powerMech[i] = 1000*gelesen;
				assert(pump.powerMech[i] > 0);
	}
	file >> skip;

	//Reading in Pumpefficiencies [Second line of "Pumpendaten.txt"]; skips end of line; Closes file
	for(int i=0;i<NPUMP;i++)
	{
		file>>gelesen;
		pump.efficiency[i] = gelesen;
				assert(pump.efficiency[i] > 0);
	}
	file >> skip;
file.close();


	// Setup for Reading Costparameters
file.open("Kosten.txt");

	//Reads DIN�s for naming available pipes [First line of "Kosten.txt"]; Skips end of line
	for(int i=0; i<NCOST; i++)
	{
		file >> gelesen;
		cost.dIN[i] = gelesen;
			assert(cost.dIN[i] > 0);
	}

	file >> skip;

	//Reads Inside Diameter of available pipes [Second line of "Kosten.txt"]; Skips end of line
	for(int i=0; i<NCOST; i++)
	{
		file >> gelesen;
		cost.insideDiameter[i] = gelesen/1000;
				assert(cost.insideDiameter[i] > 0);
	}

	file >> skip;
	cout <<endl;

	//Reads Wall Thickness of available pipes [Third line of "Kosten.txt"]; Skips end of line
	for(int i=0; i<NCOST; i++)
	{
		file >> gelesen;
		cost.wallThickness[i] = gelesen/1000;
				assert(cost.wallThickness[i] > 0);
	}

	file >> skip;

	//Reads CALCULATED outside Diameter by adding the inner one with the wall Thickness
	for(int i=0;i<NCOST;i++)
	{
		cost.outsideDiameter[i] = cost.insideDiameter[i] + (2*cost.wallThickness[i]);
				assert(cost.outsideDiameter[i] > 0);
	}

	// Reads CALCULATED powerEL with the Quotient of powerMech & efficiency
	for(int i=0;i<NPUMP;i++)
	{
		pump.powerEL[i] = (pump.powerMech[i]/1000) / pump.efficiency[i];
				assert(pump.powerEL[i] > 0);
	}

	//Reads production Time [Fourth Line of "Kosten.txt"]; (Skips end of line); Closes file
	file >> gelesen;
	cost.operationTime = gelesen;
			assert(cost.operationTime > 0);
	file >> skip;
file.close();

} //END OF VOID read




void printRead()//Gives read data in Console
{

//User Information
cout << "Gruppe 1.5 Numerik Praktikum [Lukas Samuel Jeck]" << "          --- Einlesen der Daten aus Tabelle 1 - 4 ---"
	 << endl << endl << "Beachte: Teilweise in untersch. Einheiten eingelesen! [Angegebene Einheiten]" << endl;

cout << endl << "                  --- Rohrdaten ---" << endl<<endl;

	//Gives flow of each pipe in Console
	for(int i=0;i<NPIPE;i++)
	{
		cout << "Durchstrom Rohr " << i+1 << ": " << pipe.flow[i] << " [m�/s]"<<endl;
	}
cout<<endl;

	//Gives Deltap of each pipe in Console
	for(int i=0;i<NPIPE;i++)
	{
		cout << "Delta p Rohr " << i+1 << ": " << pipe.deltaP[i] << " [Pa]" << endl;
	}
cout<<endl;

	//Gives length of each pipe in Console
	for(int i=0;i<NPIPE;i++)
	{
		cout << "Laenge Rohr " << i+1 << ": " << pipe.length[i] << " [m]" << endl;
	}
cout << endl;

//Gives k of each pipe in Console
cout << "Rohrrauigkeit:" << pipe.k/1000 << " [m]" << endl <<endl;


cout<< "                  --- Wasserdaten ---"<<endl<<endl;

//Gives density and viscosity ind Console
cout << "Dichte Wasser: " << water.density << " [kg/m^3]" << endl << endl;
cout << "Viskositaet Wasser: " << water.viscosity << " [Pa s]" << endl << endl;


cout<<"                  --- Pumpendaten ---"<<endl<<endl;

	//Gives powerMech of each Pump in Console
	for(int i=0;i<NPUMP;i++)
	{
		cout << "Mech. Pumpenleistung von Pumpe " << i+1 << ": " << pump.powerMech[i] << " [W]" << endl;

	}
cout<<endl;

	//Gives Pumpefficiency of each Pump in Console
	for(int i=0;i<NPUMP;i++)
	{
	cout << "Pumpenwirkungsgrad von Pumpe " << i+1 << ": " << 100*pump.efficiency[i]
					<< " [%] (Mit Faktor 1/100 eingelesen)" << endl;
	}
cout <<endl;


cout << "                  --- Kostenparameter ---" << endl << endl;

	//Gives insideDiameter of each available pipe in Console
	for(int i=0; i<NCOST; i++)
	{
		cout << "Innendurchmesser des lieferbaren Rohrs mit DIN " << cost.dIN[i]
			 << ": " << cost.insideDiameter[i]<< " [m]" <<endl;
	}
	cout<<endl;

	//Gives wallThickness of each available pipe in Console
	for(int i=0; i<NCOST; i++)
	{
		cout << "Wandstaerke des lieferbaren Rohrs mit DIN " << cost.dIN[i] << ": " << cost.wallThickness[i]
						 << " [m]" <<endl;
	}
	cout<<endl;

	//Gives calculated outsideDiameter of each available pipe in Console
	for(int i=0; i<NCOST; i++)
	{
		cout << "Aussendurchmesser des lieferbaren Rohrs mit DIN " << cost.dIN[i] << ": "
							 << cost.outsideDiameter[i] << " [m]" << endl;
	}
	cout<<endl;

	//Gives calculated PowerPEL of each Pump in Console
	for(int i=0; i<NPUMP; i++)
	{
		cout << "Pumpenleistung Pel von Pumpe " << i+1 << ": " << pump.powerEL[i] << " [kW]" << endl;
	}
	cout<<endl;

//Gives operationTime in Console
cout << "Betriebszeit: " << cost.operationTime << " [h]" << endl << endl;
}//END OF VOID printread

#if 0
+Hier muss Lukas die Arrays noch in die KONSOLE ausgeben lassen
        +Evtl ausklammern
        +Auf Formatierung aufpassen [Auch Einheiten nicht vergessen!]
#endif // 0


void printResults()
{

//Console Prints: It's mentioned in the cout's what is printed below
cout<<endl<<endl << "		~~~ Konsolenausgabe der Arrays die in Funktionen berechnet werden ~~~"<<endl<<endl;
cout << "--- PressureLossFriction ---"<<endl<<endl;
	for(int i=0; i<NPUMP;i++)
	{
		for(int j=0; j<NPIPE;j++)
		{
			cout << result.pressureLossFriction[i][j] << "[Pa]" << "\t"; //i:Zeilen j:Spalten
		}
		cout<<endl;
	}
cout << endl;

cout << "Lambda von Nikuradse: " << lambdaStart <<endl;
cout<<endl;

cout << "--- Calculated Diameter ---"<<endl<<endl;
	for(int i=0; i<NPUMP;i++)
		{
			for(int j=0; j<NPIPE;j++)
			{
				cout << result.diameter[i][j] << "[m]" << "\t"; //i:Zeilen j:Spalten
			}
			cout<<endl;
		}
cout<<endl;

cout << "--- Reynolds ---"<<endl<<endl;

	for(int i=0; i<NPUMP;i++)
		{
			for(int j=0; j<NPIPE;j++)
			{
				cout << result.reynolds[i][j] << "[Einheitslos]" << "\t"; //i:Zeilen j:Spalten
			}
			cout<<endl;
		}
cout<<endl;

cout << "--- Ideal (inside)Diameters from Charts ---"<<endl<<endl;

	for(int i=0; i<NPUMP;i++)
		{
			for(int j=0; j<NPIPE;j++)
			{
				cout << result.insideDiameterChart[i][j] << "[m]" << "\t"; //i:Zeilen j:Spalten
			}
			cout<<endl;
		}
cout<<endl;

cout << "--- Ideal (outside)Diameters from Charts ---"<<endl<<endl;

	for(int i=0; i<NPUMP;i++)
		{
			for(int j=0; j<NPIPE;j++)
			{
				cout << result.outsideDiameterChart[i][j] << "[m]" << "\t"; //i:Zeilen j:Spalten
			}
			cout<<endl;
		}
cout<<endl;

cout << "--- Pressure Loss Valve ---"<<endl<<endl;

	for(int i=0; i<NPUMP;i++)
		{
			for(int j=0; j<NPIPE;j++)
			{
				cout << result.pressureLossValve[i][j] << "[Pa]" << "\t"; //i:Zeilen j:Spalten
			}
			cout<<endl;
		}
cout<<endl;

}//END OF VOID printResults




void userInformation()
{

	//COUT WITH INFORMATION WHERE DATA IS SAVED AND STUFF

}
//Functions with which calculations are done




double pressureLossFrictionFunction (int i,int j)//Calculates pressureLossFriction 1DIM [Initialisation of an array in main]
{

	return (pump.powerMech[i]/pipe.flow[j])-pipe.deltaP[j];

}//END OF double pressureLossFrictionFunction




double nikuradse(double k)//Calculates a value for lambdaStart [Initialisation of double lambdaStart in main]
{

    double D_start = 0.5; // assume the start diameter is 0.5[m]
	return 1 / pow((2 * log10(3.71 * D_start / pipe.k)), 2.0);

}//END OF double nikuradse




double diameterFunction (int i,int j)//Calculates Diameter 1DIM [Initialisation of an array in main]
{

	return pow(((8.0*lambda*pipe.length[j]*water.density*pipe.flow[j]*pipe.flow[j])/(M_PI*M_PI*result.pressureLossFriction[i][j])),0.2);

}//END OF double** diameterFunction




double reynoldsFunction (int i,int j)//Calculates reynolds number 1DIM [Initialisation of an array in main]
{

	return (4*pipe.flow[j]*water.density)/(water.viscosity*diameterFunction(i,j)*M_PI);

}//END OF double reynoldsFunction




//NEWTON PROCESS Declarations
double colebrookFunction (int i, int j, double lambda)
{

    return 1.74 - 2.0*log10((2.0*pipe.k/diameterFunction(i,j))+18.7/(reynoldsFunction(i,j)*sqrt(lambda))) - 1/sqrt(lambda);

}//END OF colebrookfunction




double colebrookDerivativeFunction (int i, int j, std::function<double(int, int, double)>func)
{
       return (func(i, j, lambda+EPSILON) - func(i ,j ,lambda-EPSILON))/(2*EPSILON);

}




double newtonFunction (int i, int j, double lambda, std::function<double(int,int,double)>func)
{
    do{
        newtonHelp  = lambda;
        newton      = lambda - (func(i, j, lambda)/colebrookDerivativeFunction(i, j, colebrookFunction));
        lambda      = newton;

    } while (fabs(newton - newtonHelp) > EPSILON);

    return lambda;
}




void pickDiameterFromChart () //insideDiameter from chart out of read data, idealDiameter from diameterFunction (do this in main)
{

    for(int i=0; i<NPUMP; i++){
        for(int j=0; j<NPIPE; j++){

            int m=0;

            while(result.diameter[i][j] > cost.insideDiameter[m])//go through all the data of insideDiameter until you find a diameter bigger than the calculated one
            {
                m++;
                	assert(m < NCOST);
            }

            result.insideDiameterChart[i][j]        = cost.insideDiameter[m]; //save the chart inside Diameters

            result.outsideDiameterChart[i][j] = cost.insideDiameter[m]+2*cost.wallThickness[m]; //save the chart outside Diameters
        }
    }
}//END OF double pickDiameterFromChart




double pressureLossValveFunction(int i,int j)
{

	return (pump.powerMech[i]/pipe.flow[j]) - pipe.deltaP[j] - ((8*lambda*pipe.length[j]*water.density*pipe.flow[j]*pipe.flow[j]) / (M_PI*M_PI*pow(result.insideDiameterChart[i][j],5.0)));

}//END OF double pressureLossValveFunction

#if 0
+Hier muss Lukas sich die Formatierung noch anschauen [Einheiten nicht vergessen!]
#endif // 0



void output(){

//Setting up output(Creating file operator)
out.open("Cooltech_Daten.txt");

//Information (Headline of "Cooltech_Daten.txt")
	out << "Gruppe 1.5 Numerik Praktikum [Lukas Samuel Jeck]" << "            --- Ausgeben der Dateien in eine Textdatei ---";
	out << endl<<endl<<endl;

//Defining what data is shown
out << "~~~ Daten der Rohrauslegung ~~~"<<endl<<endl<<endl;


//Putting out ideal Diameters of pipes for each Pipe
out << "--- Durchmesser[Pumpe][Rohr] ---"<<endl<<endl;
for(int i=0;i<NPUMP;i++)
{

	for(int j=0;j<NPIPE;j++)
	{
      out << "Durchmesser[" << i+1 << "]" << "[" << j+1 << "] " << result.insideDiameterChart[i][j] << "\t";
	}
	out<<endl;
}
out<<endl;


//Putting out Pressureloss of Valve
out << "--- Druckverlust f�r Ausgleichsventil[Pumpe][Rohr] ---"<<endl<<endl;
for(int i=0;i<NPUMP;i++)
{

	for(int j=0;j<NPIPE;j++)
	{
      out << "Ventildruck[" << i+1 << "]" << "[" << j+1 << "] " << result.pressureLossValve[i][j] << "\t";
	}
	out<<endl;
}
out<<endl<<endl;


//Defining what Data is shown
out << "~~~ Daten zur Kostenberechnung ~~~"<<endl<<endl<<endl;


//Putting out the Total Cost of each Pump
out << "--- Kosten Pumpen[Pumpe] ---"<<endl<<endl;
for(int i=0;i<NPUMP;i++)
{
	out << "Kosten Pumpe[" << i+1 << "] " << result.pumpCost[i]<< "\t";
}
out<<endl;


out << "--- Kosten Rohre[Rohre] ---"<<endl<<endl;
for(int j=0;j<NPIPE;j++)
{
	out << "Kosten Rohr[" << j+1 << "] " << result.pipeCost[j]<< "\t";
}
out<<endl;


out << "--- Kosten Strom[Pumpe] ---"<<endl<<endl;
for(int i=0;i<NPUMP;i++)
{
	out << "Kosten Strom[" << i+1 << "] " << result.pumpCost[i]<< "\t";
}
out<<endl;


out << "--- Gesamtkosten[Pumpe][Rohr] ---"<<endl<<endl;
for(int i=0;i<NPUMP;i++)
{

	for(int j=0;j<NPIPE;j++)
	{
      out << "Gesamtkosten[" << i+1 << "]" << "[" << j+1 << "] " << result.totalCost[i][j] << "\t";
	}
	out<<endl;
}
out<<endl<<endl;

out.close();
} //END of Void output

/*------------------------------------------------------- END OF PUMPCALCULATION -----------------------------------------------------------*/
#if 0
Hier muss Lukas noch Arrays f�r 1DIM Arrays in main belegen
#endif // 0

double pumpCostFunction(int i)					//powerEl in kW
{

    return pump.powerEL[i] * 406 + 4011 * (5 - pump.efficiency[i]);

}




double pipeCostFunction(int i, int j)					//diameter in m, lenght in m
{

	return (result.outsideDiameterChart[i][j] * result.outsideDiameterChart[i][j] * 16458 - result.outsideDiameterChart[i][j] * 2109 + 151) * pipe.length[j];

}




double powerCostFunction(int i)						//powerEl in kW, time in h
{

	 return cost.operationTime * pump.powerEL[i] * 0,05;

}




double totalCostFunction(int i,int j)
{

	return pumpCostFunction(i) + pipeCostFunction(i,j) + powerCostFunction(i);

}




/*--------------------------------------------------------- POLYNOMIAL FIT -----------------------------------------------------------------*/


double* polynomialFit(int degree, int NPUMP, double* x, double* y)	//add cin>> degree for user
														//degree is the degree of the polynom, x[] and y[] the valuepairs
{
    int i, j, k;										//counting Variables
    double a[degree + 1];								//Array that stores the coefficents of the polynom
    double X[2 * degree + 1];							//Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    for (i = 0 ; i < 2 * degree + 1 ; i++)
    {
        X[i]=0;
        for (j = 0 ; j < NPUMP ; j++)
        {
			X[i] = X[i] + pow(x[j],i);
		}
    }
    double B[degree + 1][degree + 2];					//B is the Normal matrix that will store the equations
    for (i = 0 ; i <= degree ; i++)
    {
        for (j = 0 ; j <= degree ; j++)
        {
            B[i][j] = X[i + j];
        }
    }													//Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
    double Y[degree + 1];                    			//Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    for (i = 0 ; i < degree + 1 ; i++)
    {
        Y[i] = 0;
        for (j = 0 ; j < NPUMP ; j++)
        {
			Y[i] = Y[i] + pow(x[j],i) * y[j];
		}
    }
    for (i = 0 ; i <= degree ; i++)
    {
        B[i][degree+1] = Y[i];
    }               									//load the values of Y as the last column of B


    for (i = 0 ; i < degree + 1 ; i++)                 	//From now Gaussian Elimination starts to solve the set of linear equations (Pivotisation)
     {
        for (k = i + 1 ; k < degree + 1 ; k++)
        {
            if (B[i][i] < B[k][i])
            {
                for (j = 0 ; j <= degree + 1 ; j++)
                {
                    double temp = B[i][j];
                    B[i][j] = B[k][j];
                    B[k][j] = temp;
                }
			}
		}
	}
    for (i = 0 ; i < degree ; i++)
    {         											//loop to perform the gauss elimination
        for (k = i + 1 ; k < degree + 1 ; k++)
            {
                double t = B[k][i] / B[i][i];
                for (j = 0 ; j <= degree + 1 ; j++)
                {
                    B[k][j] = B[k][j] - t * B[i][j];    //make the elements below the pivot elements
				}
            }
    }
    for (i = degree ; i >= 0 ; i--)               		//back-substitution
    {
        a[i] = B[i][degree + 1];                		//make the variable to be calculated equal to the rhs of the last equation
        for (j = 0 ; j < degree + 1 ; j++)
        {
            if (j != i)
            {            								//then subtract all the lhs values except the coefficient of the variable whose value is being calculated
                a[i] = a[i] - B[i][j] * a[j];
			}
		}
        a[i] = a[i] / B[i][i];           				//now finally divide the rhs by the coefficient of the variable to be calculated
    }
    return a;
}
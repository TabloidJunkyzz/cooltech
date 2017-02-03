//++++++++++++++++++++++ "Cooltech.cpp" Implementation for our header "Cooltech.h" ++++++++++++++++++++++

//Includes Header which includes all necessary libraries/headers
#include "Cooltech.h"
using namespace std;


//########################################### Variables and Structs ###########################################

//Variables Declarations here to make them local ! ##Variables are commented on in our Header##

//Variables for read, output functions and all console commands
ifstream file;
string   skip;
double   reading;
double   totalFlow = 0; // Set to 0 for Summation

ofstream out;

string version = "";

//Variables for functions of pipe measurements
double lambda;
double lambdaStart;
double lambdaTemp;

double newton;
double newtonHelp;

//Array for checkForMistake --> Saves outcome of every line checked into an array
bool boolMistake[NPUMP];

//Variables for determineIntervallForMinCost
int nFunctionalPumps             = 0;  //Counts Functional Pipes
int nSettingFunctionalPumpsArray = 0;  //Proceedure Variable to set Values in New Array

//Polynomialfit
double  pipeCostStorage;
double* coefficents;
int     degree;

//Minimum calculation Variables
double xMin;
double minCost;
double iBegin;
double iEnd;
double precision;

//Struct variables
data pipe, water, pump, cost, result;

//########################################### Reading Routine ###########################################

#if 0
By Lukas Samuel Jeck
#endif // 0
//Saves data from charts in arrays
void read()
{
//Setup for Reading Pipe Data
file.open("Rohrdaten.txt");

        //Reads Flow of each pipe [First line of "Rohrdaten.txt"]; Skips end of line
        for(int iPipes = 0;iPipes < NPIPE; iPipes++)
        {
            file >> reading;
            pipe.flow[iPipes] = 0.001*reading;
                    assert(pipe.flow[iPipes] >= 0);

            totalFlow += pipe.flow[iPipes]; //Adds up all flows for "V(point)"
        }
	file >> skip;

        //Reads Pressureloss of each pipe [Second Line of "Rohrdaten.txt"]; Skips end of line
        for(int iPipes = 0;iPipes < NPIPE; iPipes++)
        {
            file >> reading;
            pipe.deltaP[iPipes] = 1e5*reading;
                    assert(pipe.deltaP[iPipes] >= 0);
        }
	file >> skip;

        //Reads Pipelength of each pipe [Third line of "Rohrdaten.txt"]; Skips end of line
        for(int iPipes = 0;iPipes < NPIPE; iPipes++)
        {
            file >> reading;
            pipe.length[iPipes] = reading;
                    assert(pipe.length[iPipes] >= 0);
        }
	file >> skip;

	/*
	Reads k [Fourth Line of "Rohrdaten.txt"]; Skips end of line
	Then closes file, cause every pipedata is in the system!
	*/
	file >> reading;
	pipe.k = reading/1000;
		assert(pipe.k > 0);
file.close();






	//Setup for Reading water data
file.open("Wasserdaten.txt");

	//Reads Water density [First line of "Wasserdaten.txt"]; Skips end of line
	file >> reading;
	water.density = reading;
		assert(water.density > 0);
	file >> skip;

	//Reads Water viscosity [Second line of "Wasserdaten.txt"]; Skips end of line; Closes file
	file >> reading;
	water.viscosity = reading;
			assert(water.viscosity > 0);
file.close();






	//Setup for reading pump data
file.open("Pumpendaten.txt");

        //Reads Pumppowers [First line of "Pumpendaten.txt"]; Skips end of line
        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            file >> reading;
            pump.powerMech[iPumps] = 1000*reading;
                    assert(pump.powerMech[iPumps] > 0);
        }
	file >> skip;

        //Reading in Pumpefficiencies [Second line of "Pumpendaten.txt"]; skips end of line; Closes file
        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            file >> reading;
            pump.efficiency[iPumps] = reading;
                    assert(pump.efficiency[iPumps] > 0);
        }
	file >> skip;
file.close();







	// Setup for Reading Costparameters
file.open("Kosten.txt");

        //Reads DINs for naming available pipes [First line of "Kosten.txt"]; Skips end of line
        for(int iCost = 0; iCost < NCOST; iCost++)
        {
            file >> reading;
            cost.dIN[iCost] = reading;
                assert(cost.dIN[iCost] > 0);
        }

	file >> skip;

        //Reads Inside Diameter of available pipes [Second line of "Kosten.txt"]; Skips end of line
        for(int iCost = 0; iCost < NCOST; iCost++)
        {
            file >> reading;
            cost.insideDiameter[iCost] = reading/1000;
                    assert(cost.insideDiameter[iCost] > 0);
        }

	file >> skip;

        //Reads Wall Thickness of available pipes [Third line of "Kosten.txt"]; Skips end of line
        for(int iCost = 0; iCost < NCOST; iCost++)
        {
            file >> reading;
            cost.wallThickness[iCost] = reading/1000;
                    assert(cost.wallThickness[iCost] > 0);
        }

	file >> skip;

        //Reads outside Diameter by adding the inner one with the wall Thickness
        for(int iCost = 0;iCost < NCOST; iCost++)
        {
            cost.outsideDiameter[iCost] = cost.insideDiameter[iCost] + (2*cost.wallThickness[iCost]);
                    assert(cost.outsideDiameter[iCost] > 0);
        }

        // Reads CALCULATED powerEl with the Quotient of powerMech & efficiency
        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            pump.powerEl[iPumps] = (pump.powerMech[iPumps] / 1000) / pump.efficiency[iPumps];
                    assert(pump.powerEl[iPumps] > 0);
        }

	//Reads production Time [Fourth Line of "Kosten.txt"]; (Skips end of line); Closes file
	file >> reading;
	cost.operationTime = reading;
			assert(cost.operationTime > 0);
	file >> skip;
file.close();

} //END OF VOID read

//########################################### Console ###########################################

//Prints out saved arrays of "read"-function in console
void printRead()
{
    //General Informations
    cout << "Group 1.5 Numerical Project" << "          --- Console Outputs in Dev-Version ---"<<endl<<endl;
    cout << "		~~~ Reading Routine: Caution! Certain datas are in a different unity! [Unity is shown] ~~~"<<endl<<endl;

    cout <<endl<< "                  --- Pipe Data ---"<<endl<<endl;

        //Gives flow of each pipe in Console
        for(int iPipes = 0; iPipes < NPIPE; iPipes++)
        {
            cout << "Pipe Flow " << iPipes+1 << ": " << setw(7) << left << pipe.flow[iPipes] << " [m^3/s]"<<endl;
        }
    cout <<endl;

        //Gives Deltap of each pipe in Console
        for(int iPipes = 0; iPipes < NPIPE; iPipes++)
        {
            cout << "Delta P Pipe " << iPipes+1 << ": " << setw(7) << left << pipe.deltaP[iPipes] << " [Pa]"<<endl;
        }
    cout <<endl;

        //Gives length of each pipe in Console
        for(int iPipes = 0;iPipes < NPIPE; iPipes++)
        {
            cout << "Pipe Length " << iPipes+1 << ": " << setw(5) << left << pipe.length[iPipes] << " [m]"<<endl;
        }
    cout <<endl;

    //Gives k of each pipe in Console
    cout << "Pipe Roughness:" << pipe.k / 1000 << " [m]"<<endl<<endl;


    cout << "                  --- Water Data ---"<<endl<<endl;

    //Gives density and viscosity ind Console
    cout << "Density Water: " << water.density << " [kg/m^3]"<<endl<<endl;
    cout << "Viscosity Water: " << water.viscosity << " [Pa s]"<<endl<<endl;


    cout << "                  --- Pump Data ---"<<endl<<endl;

        //Gives powerMech of each Pump in Console
        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            cout << "Mech. Power of Pump " << iPumps+1 << ": " << setw(8) << left << pump.powerMech[iPumps] << " [W]"<<endl;

        }
    cout <<endl;

        //Gives Pumpefficiency of each Pump in Console
        for(int iPumps = 0;iPumps < NPUMP;iPumps++)
        {
            cout << "Efficiency for Pump " << iPumps+1 << ": " << 100*pump.efficiency[iPumps]
                 << " [%] (Read with  Factor 1/100)"<<endl;
        }
    cout <<endl;


    cout << "                  --- Cost Parameter ---"<<endl<<endl;

        //Gives insideDiameter of each available pipe in Console
        for(int iCost = 0; iCost < NCOST; iCost++)
        {
            cout << "Inside Diameter of available Pipes by DIN " << setw(3) << right << cost.dIN[iCost]
                 << ": " << setw(6) << left << cost.insideDiameter[iCost]<< " [m]"<<endl;
        }
    cout <<endl;

        //Gives wallThickness of each available pipe in Console
        for(int iCost = 0; iCost < NCOST; iCost++)
        {
            cout << "Wall Thickness of available Pipe by DIN " << setw(3) << right << cost.dIN[iCost] << ": " << setw(8) << left
                 << cost.wallThickness[iCost] << " [m]"<<endl;
        }
    cout <<endl;

        //Gives calculated outsideDiameter of each available pipe in Console
        for(int iCost = 0; iCost < NCOST; iCost++)
        {
            cout << "Outside Diameter of available Pipe by DIN " << setw(3) << right << cost.dIN[iCost] << ": "
                 << setw(8) << left << cost.outsideDiameter[iCost] << " [m]"<<endl;
        }
    cout <<endl;

        //Gives calculated PowerPEL of each Pump in Console
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            cout << "powerEl of Pump " << iPumps+1 << ": " << setw(8) << left << pump.powerEl[iPumps] << " [kW]"<<endl;
        }
    cout <<endl;

    //Gives operationTime in Console
    cout << "Operation Time: " << cost.operationTime << " [h]"<<endl<<endl;
}//END OF VOID printread




//Prints calculated data in Console
void printResults(double* powerElSort)
{
    //It's mentioned in the cout's what is printed below
    cout <<endl << "		~~~ Consoleprint of Arrays that are calculated within the Functions ~~~"<<endl<<endl;
    cout << "		~~~ Rows are for Pumps and Columns are for Pipes ~~~"<<endl<<endl;

    cout << "--- PressureLossFriction ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            for(int iPipes = 0; iPipes < NPIPE; iPipes++)
            {
                cout << setw(12) << left << result.pressureLossFriction[iPumps][iPipes] << "[Pa]" << "    "; //iPumps:Rows iPipes:Columns
            }
        cout <<endl<<endl;
        }

    cout << "--- Reynolds Number ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            for(int iPipes = 0; iPipes < NPIPE; iPipes++)
            {
                cout << setw(15) << left << result.reynolds[iPumps][iPipes] << "    "; //iPumps:Rows iPipes:Columns
            }
        cout <<endl<<endl;
        }
    cout <<endl;

    cout << "############# Lambda Nikuradse: " << lambdaStart << " #############"<<endl<<endl;

    cout << "--- Calculated Diameter ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            for(int iPipes = 0; iPipes < NPIPE; iPipes++)
            {
                cout << setw(12) << left << result.diameter[iPumps][iPipes] << "[m]" << "    "; //iPumps:Rows iPipes:Colums
            }
        cout <<endl<<endl;
        }
    cout <<endl;

    cout << "--- Ideal (inside)Diameters from Chart ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            for(int iPipes = 0; iPipes < NPIPE; iPipes++)
            {
                cout << setw(12) << left << result.insideDiameterChart[iPumps][iPipes] << "[m]" << "    "; //iPumps:Rows iPipes:Columns
            }
        cout <<endl<<endl;
        }
    cout <<endl;

    cout << "--- Ideal (outside)Diameters (from Chart) ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            for(int iPipes = 0; iPipes < NPIPE; iPipes++)
            {
                cout << setw(12) << left << result.outsideDiameterChart[iPumps][iPipes] << "[m]" << "    "; //iPumps:Rows iPipes:Columns
            }
        cout <<endl<<endl;
        }
    cout <<endl;

    cout << "--- Pressure Loss Valve ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            for(int iPipes = 0; iPipes < NPIPE; iPipes++)
            {
                cout << setw(12) << left << result.pressureLossValve[iPumps][iPipes] << "[Pa]" << "    "; //iPumps:Rows iPipes:Columns
            }
        cout <<endl<<endl;
        }
    cout <<endl;

    cout << "--- Cost Pipes[Pipe] ---"<<endl<<endl;
        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            for(int iPipes = 0;iPipes < NPIPE; iPipes++)
            {
                cout << "Cost Pipe[" << iPipes+1 << "] " << result.pipeCost[iPumps][iPipes] << "    ";
            }
        cout <<endl<<endl;
        }
    cout <<endl;

    cout << "--- Cost Pumps[Pump] ---"<<endl<<endl;
        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            cout << "Cost Pump[" << iPumps+1 << "] " << result.pumpCost[iPumps]<<endl;
        }
	cout <<endl<<endl;


    cout << "--- Cost Power[Pump] ---"<<endl<<endl;
        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            cout << "Cost Power[" << iPumps+1 << "] " << result.powerCost[iPumps]<<endl;
        }
	cout <<endl<<endl;


    cout << "--- Total Costs[Pumpe] ---"<<endl<<endl;
        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            cout << "Total cost[" << iPumps+1 << "] " << setw(6) << left << result.totalCost[iPumps]<<endl;
        }
	cout <<endl<<endl;

	cout << "--- Mistake Routine (1 Equals true / 0 Equals false) ---"<<endl<<endl;

        for(int iPumps = 0;iPumps < NPUMP; iPumps++)
        {
            cout << "Mistake[" << iPumps+1 << "] " << setw (6) << left << boolMistake[iPumps]<<endl;
        }
    cout <<endl<<endl;

	cout << "--- Coefficients of our Polynom ---"<<endl<<endl;
        for(int iPolynomial = 0; iPolynomial <= degree; iPolynomial++)
            {
                cout << "Coefficient " << iPolynomial+1 << ": " << coefficents[iPolynomial] << "   ";
            }
    cout <<endl<<endl;

    cout << "--- Minimum Calculation ---"<<endl<<endl;

    cout << "Amount of working Pumps: " << nFunctionalPumps<<endl<<endl;

    cout << "powerElSort-Array: "<<endl;

        for(int iWorkingPumps = 0; iWorkingPumps < nFunctionalPumps; iWorkingPumps++)
        {
            cout << "[" << iWorkingPumps+1 << "] " << powerElSort[iWorkingPumps]<<endl;
        }
    cout<<endl;

    cout << "Least Cost: " << minCost << " " << "In regard to Pumppower: " << xMin << "[kW]"<<endl<<endl;

}//END OF VOID printResults




//Prints a text for the User in the Console
void userInformation()
{

#if 0
Searching for ways to clear the Console after the User input Maybe ??
#endif // 0

  //system("clear");            //--> Cool easy way but "should be avoided with user interactive console management" ??
  cout << "\033[2J\033[1;1H"; //--> Could not be compatible with every Terminal (Works with Linux Terminal)

	cout << "Thank you for working with Cooltech_Industries"<<endl<<endl;
	cout << "Your Data is being Considered and Calculations are in Queue ..."<<endl<<endl;
	cout << "Pls open Cooltech_Data.txt to look at the Calculation Results"<<endl<<endl;

}//END OF VOID userInformation

//########################################### Output ###########################################

//Takes calculated data and writes it in "Cooltech_Daten.txt"
void output()
{
//Setting up output
out.open("Cooltech_Daten.txt");

	//Information (Headline of "Cooltech_Daten.txt")
	out <<endl;
	out << "Group 1.5 Numerical Project" << "            --- Writing results into a .txt File ---";
	out <<endl<<endl<<endl;


	//CALCULATIONS
	out << "~~~ Data of Pipecalculations ~~~"<<endl<<endl<<endl;


	//Output ideal Diameters of pipes for each Pipe
	out << "--- Ideal InsideDiameter[Pump][Pipe] ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            if(boolMistake[iPumps] == true)
            {
                out << "Pump " << iPumps+1 << " is not working with the given data";
            }

            else
            {
                for(int iPipes = 0; iPipes < NPIPE; iPipes++)
                {
                    out << "InsideDiameter[" << iPumps+1 << "]" << "[" << iPipes+1 << "] " << setw(6) << left << result.insideDiameterChart[iPumps][iPipes] << " [m]" << "    ";
                }
            }

        out <<endl;
        }
    out <<endl<<endl;

    //Output Pressureloss of Valve
    out << "--- Pressureloss Valve[Pump][Pipe] ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            if(boolMistake[iPumps] == true)
            {
                out << "Pump " << iPumps+1 << " is not working with the given data";
            }

            else
            {
                for(int iPipes = 0; iPipes < NPIPE; iPipes++)
                {
                    out << "Valve Pressureloss[" << iPumps+1 << "]" << "[" << iPipes+1 << "] " << setw(7) << left << result.pressureLossValve[iPumps][iPipes] << " [Pa]" << "    ";
                }
            }

        out <<endl;
        }
    out <<endl<<endl;

	//COSTS
	out << "~~~ Data of Costcalculations ~~~"<<endl<<endl<<endl;

    //Output Cost of the Pipes
    out << "--- Cost Pipes[Pump][Pipe] ---"<<endl<<endl;
        for(int iPumps = 0; iPumps < NPUMP; iPumps++)
        {
            if(boolMistake[iPumps] == true)
            {
                out << "Pump " << iPumps+1 << " is not working with the given data";
            }

            else
            {
                for(int iPipes = 0; iPipes < NPIPE; iPipes++)
                {
                    out << "Cost Pipe[" << iPumps+1 << "]" << "[" << iPipes+1 << "] " << setw(7) << left << result.pipeCost[iPumps][iPipes] << " [€]" << "    ";
                }
            }

        out <<endl;
        }
    out <<endl<<endl;

    //Output Cost of Pumps
	out << "--- Cost Pumps[Pump] ---"<<endl<<endl;
        for(int iPumps = 0;iPumps < NPUMP;iPumps++)
        {
            if(boolMistake[iPumps] == true)
            {
                out << "Pump " << iPumps+1 << " is not working with the given data"<<endl;
            }

            else
            {
                out << "Cost Pump[" << iPumps+1 << "] " << setw(8) << left << result.pumpCost[iPumps] << "[€]"<<endl;
            }
        }
	out <<endl<<endl;

	//Output Cost of Power for Pumps
	out << "--- Cost Power[Pump] ---"<<endl<<endl;
        for(int iPumps = 0;iPumps < NPUMP;iPumps++)
        {
            if(boolMistake[iPumps] == true)
            {
               out << "Pump " << iPumps+1 << " is not working with the given data"<<endl;
            }

            else
            {
            out << "Cost Power[" << iPumps+1 << "] " << setw(8) << left << result.powerCost[iPumps]<< "[€]"<<endl;
            }
        }
	out <<endl<<endl;

	//Output Total Costs
	out << "--- Total Costs[Pumpe] ---"<<endl<<endl;
        for(int iPumps = 0;iPumps < NPUMP;iPumps++)
        {
            if(boolMistake[iPumps] == true)
            {
                out << "Pump " << iPumps+1 << " is not working with the given data"<<endl;
            }

            else
            {
            out << "Total Cost[" << iPumps+1 << "] " << setw(7) << left << result.totalCost[iPumps]<< "[€]"<<endl;
            }
        }
	out <<endl<<endl;

	//Output Coefficients of our Polynom
	out << "--- Coefficients of our Polynom ---"<<endl<<endl;
        for(int iPolynomial = 0; iPolynomial <= degree; iPolynomial++)
        {
            out << "Coefficient " << iPolynomial+1 << ": " << coefficents[iPolynomial] << "   ";
        }
    out <<endl<<endl;

    //Output Polynomial
    out << "--- The Polynom that we fitted would then look like this: ---"<<endl<<endl;

        out << "Polynomial of degree " << degree << ":   ";
        for(int iPolynomial = 0; iPolynomial <= degree; iPolynomial++)
        {
            //When Degree = 0 only !
            if(degree < 1 && degree >= 0)
            {
                out << "(" << coefficents[0] << ")";
            }

            else
            {
                //When Degree is higher than 0 --> Prints "firstCoefficent + "
                if(iPolynomial < 1 && iPolynomial >= 0 && degree > 0)
                {
                    out << "(" << coefficents[iPolynomial] << ")" << " + ";
                }

                //When Degree = 1 only! --> Prints out "secondCoefficent*x"
                if(iPolynomial < 2 && iPolynomial >= 1   && degree <= 1 && degree > 0)
                {
                    out << "(" << coefficents[iPolynomial] << ")" << "*x";
                }

                //When Degree is higher than 1 --> Prints "secondCoefficent*x + "
                if(iPolynomial < 2 && iPolynomial >= 1 && degree > 1)
                {
                    out << "(" << coefficents[iPolynomial] << ")" << "*x + ";
                }

                //When Degree is higher than 1 --> Prints "thirdCoefficent*x^2 + ... nCoefficent*x^n-1"
                if(iPolynomial > 1  && iPolynomial < degree && degree > 1)
                {
                    out << "(" << coefficents[iPolynomial] << ")" << "*x^" << iPolynomial << " + ";
                }

                //When Degree is higher than 1 --> Prints "n+1Coefficent*x^n"
                if(iPolynomial > degree - 1 && iPolynomial <= degree && degree > 1)
                {
                    out << "(" << coefficents[iPolynomial] << ")" << "*x^" << iPolynomial;
                }
            }
        }
    out <<endl<<endl;

    out << "--- Minimum Calculation ---"<<endl<<endl;

    out << "Least Cost: " << minCost << "[€]" << " " << "In regard to Pumppower: " << xMin << "[kW]"<<endl<<endl;

out.close();

} //END of Void output

//########################################### Mistake Routine ###########################################

//Checks PressurelossFriction Array for negative Values and saves those outcomes by line in boolMistake
void checkForMistake()
{
    //First sets the Array to false in every entry to not interrupt my boolean commands and make sure that the array has no strange values in it
    for(int iPumps = 0; iPumps < NPUMP; iPumps++)
    {
        boolMistake[iPumps] = false;
    }

    //Checks Every line for a negative Value and then saves "true" or "false" for that in our !boolMistake!
    for(int iPumps = 0; iPumps < NPUMP; iPumps++)
    {
        for(int iPipes = 0; iPipes < NPIPE; iPipes++)
        {
            if(result.pressureLossFriction[iPumps][iPipes] < 0)
            {
                boolMistake[iPumps] = true;
            }

            /*Important here is the second boolean request --> If negative Value in a line --> outcome stays true and does not get set to false once
            there is a value after the negative one thats positive!*/
            if(result.pressureLossFriction[iPumps][iPipes] > 0 && boolMistake[iPumps] != true)
            {
                boolMistake[iPumps] = false;
            }
        }
    }

} //END OF VOID checkForMistake
#if 0
By Lukas Samuel Jeck
#endif // 0

//########################################### Calculations for Pipedimensioning ###########################################

#if 0
By Hermann Gegel & Boda Yang [Lukas Samuel Jeck as a Helping Hand]
#endif // 0
//Calculates pressureLossFriction 1DIM [Initialisation of an array in main]
double pressureLossFrictionFunction(int iPumps, int iPipes)
{
	return (pump.powerMech[iPumps] / totalFlow) - pipe.deltaP[iPipes];

}//END OF double pressureLossFrictionFunction




//Calculates a value for lambdaStart [Initialisation of double lambdaStart in main]
double nikuradse(double k)
{
    double D_start = 0.5; // assume the start diameter is 0.5[m]
	return 1 / pow((2*log10(3.71*D_start / pipe.k)), 2.0);

}//END OF double nikuradse




//Calculates Diameter 1DIM [Initialisation of an array in main]
double diameterFunction(int iPumps, int iPipes)
{
	return pow(((8.0*lambda*pipe.length[iPipes]*water.density*pipe.flow[iPipes]*pipe.flow[iPipes]) / (M_PI*M_PI*result.pressureLossFriction[iPumps][iPipes])), 0.2);

}//END OF double diameterFunction




//Calculates reynolds number 1DIM [Initialisation of an array in main]
double reynoldsFunction(int iPumps, int iPipes)
{
	return (4*pipe.flow[iPipes]*water.density) / (water.viscosity*diameterFunction(iPumps,iPipes)*M_PI);

}//END OF double reynoldsFunction

//########################################### Newton ###########################################

//Colebrookfunction
double colebrookFunction(int iPumps, int iPipes, double lambda)
{
    return 1.74 - 2.0*log10((2.0*pipe.k / diameterFunction(iPumps,iPipes)) + 18.7 / (reynoldsFunction(iPumps,iPipes)*sqrt(lambda))) - 1 / sqrt(lambda);

}//END OF colebrookfunction



//Derivating Colebrook numerically
double colebrookDerivativeFunction(int iPumps, int iPipes, std::function<double(int, int, double)>func)
{
       return (func(iPumps, iPipes, lambda+EPSILON) - func(iPumps, iPipes, lambda-EPSILON)) / (2*EPSILON);

}//END OF colebrookDerivativeFunction



//Solves Colebrook for lambda
double newtonFunction(int iPumps, int iPipes, double lambda, std::function<double(int,int,double)>func)
{
    do
        {
            newtonHelp  = lambda;
            newton      = lambda - (func(iPumps, iPipes, lambda) / colebrookDerivativeFunction(iPumps, iPipes, colebrookFunction));
            lambda      = newton;

        } while (fabs(newton - newtonHelp) > EPSILON);

    return lambda;

}//END OF double newtonFunction

//########################################### Pick from Chart ###########################################

//Picks insideDiameter from chart, Calculates outsideDiameter with insideDiameter from chart and corresponding Wall Thickness
void pickDiameterFromChart()
{
    for(int iPumps = 0; iPumps < NPUMP; iPumps++)
    {
        for(int iPipes = 0; iPipes < NPIPE; iPipes++)
        {

            int m = 0;

            while(result.diameter[iPumps][iPipes] > cost.insideDiameter[m])//go through all the data of insideDiameter until you find a diameter bigger than the calculated one
            {
                m++;
                	assert(m < NCOST);
            }

            result.insideDiameterChart[iPumps][iPipes]  = cost.insideDiameter[m]; //save the picked Diameter in Array

            result.outsideDiameterChart[iPumps][iPipes] = cost.insideDiameter[m] + 2*cost.wallThickness[m]; //Calculate outsideDiameter
        }
    }
}//END OF VOID pickDiameterFromChart
#if 0
By Hermann Gegel & Boda Yang [Lukas Samuel Jeck as a Helping Hand]
#endif // 0

//########################################### Valve ###########################################

#if 0
By Lina Lepp
#endif // 0
//Calculates the, to be regulated, pressure Loss by the Valve
double pressureLossValveFunction(int iPumps, int iPipes)
{
	return (8*lambda*pipe.length[iPipes]*water.density*pipe.flow[iPipes]*pipe.flow[iPipes]) / (M_PI*M_PI*pow(result.diameter[iPumps][iPipes], 5.0)) - ((8*lambda*pipe.length[iPipes]*water.density*pipe.flow[iPipes]*pipe.flow[iPipes]) / (M_PI*M_PI*pow(result.insideDiameterChart[iPumps][iPipes],5.0)));

}//END OF double pressureLossValveFunction
#if 0
By Lina Lepp
#endif // 0

//########################################### Cost Calculations ###########################################

#if 0
By Lukas Samuel Jeck
#endif // 0
//Calculates Pumpcost
double pumpCostFunction(int iPumps)
{
    return pump.powerEl[iPumps]*406 + 4011*(5 - pump.efficiency[iPumps]);

}//END OF pumpCostFunction




//Calculates pipeCost
double pipeCostFunction(int iPumps, int iPipes)
{
	return (result.outsideDiameterChart[iPumps][iPipes]*result.outsideDiameterChart[iPumps][iPipes]*16458 - result.outsideDiameterChart[iPumps][iPipes]*2109 + 151)*pipe.length[iPipes];

}//END OF pipeCostFunction




////Calculates powerCost
double powerCostFunction(int iPumps)
{
	 return cost.operationTime*pump.powerEl[iPumps]*0.05;

}//END OF powerCostFunction




////Calculates totalCost
double totalCostFunction(int iPumps)
{
    pipeCostStorage = 0;
        for (int iPipes = 0; iPipes < NPIPE; iPipes++)
        {
            pipeCostStorage += pipeCostFunction(iPumps,iPipes);
        }
	return pumpCostFunction(iPumps) + pipeCostStorage + powerCostFunction(iPumps);

}//END OF totalCostFunction
#if 0
By Lukas Samuel Jeck
#endif // 0

//########################################### Polynomial Fit ###########################################

#if 0
By Jakob Mangold .... Did not get his Part of the Programm done in time and did not respond, so Lukas Samuel Jeck had to fix it
#endif // 0
//PonynomialFit for our Total Cost Function
void polynomialFit()
{
    //pump.powerEl[iPumps] is considered the x value
    //result.totalCost[iPumps] is considered the y value
    //Those are used as setting points for our PolynomialFit

    //counting Variables
    int iPumps, kDegree, lDegree, mDegree;

    //coefficients array gets set to 0 for starting purpose
    coefficents = new double[degree+1];
    for (kDegree = 0; kDegree < degree; kDegree++)
    {
        coefficents[kDegree] = 0;
    }

    //powerElArray will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    double powerElArray[2 * degree + 1];
    for (kDegree = 0 ; kDegree < 2 * degree + 1 ; kDegree++)
    {
        powerElArray[kDegree] = 0;
        for (iPumps = 0 ; iPumps < NPUMP ; iPumps++)
        {
            if (boolMistake[iPumps] == false)
            {
                powerElArray[kDegree] = powerElArray[kDegree] + pow(pump.powerEl[iPumps],kDegree);
            }
		}
    }

    //Normal matrix that will store the equations
    double normalMatrix[degree + 1][degree + 2];
    for (kDegree = 0 ; kDegree <= degree ; kDegree++)
    {
        for (lDegree = 0 ; lDegree <= degree ; lDegree++)
        {
            normalMatrix[kDegree][lDegree] = powerElArray[kDegree + lDegree]; //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
        }
    }

    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    double resultArray[degree + 1];

    for (kDegree = 0 ; kDegree < degree + 1 ; kDegree++)
    {
        resultArray[kDegree] = 0;

            for (iPumps = 0 ; iPumps < NPUMP ; iPumps++)
            {
                if (boolMistake[iPumps] == false)
                {
                    resultArray[kDegree] = resultArray[kDegree] + pow(pump.powerEl[iPumps],kDegree) * result.totalCost[iPumps];
                }
            }
    }

    for (kDegree = 0 ; kDegree <= degree ; kDegree++)
    {
        normalMatrix[kDegree][degree+1] = resultArray[kDegree]; //load the values of resultArray as the last column of normalMatrix
    }

    //From now Gaussian Elimination starts to solve the set of linear equations (Pivotisation) --> We need "degree + 1" Equations for Gauss
    for (lDegree = 0 ; lDegree < degree + 1 ; lDegree++)
     {
        for (kDegree = lDegree + 1 ; kDegree < degree + 1 ; kDegree++)
        {
            if (normalMatrix[lDegree][lDegree] < normalMatrix[kDegree][lDegree])
            {
                for (mDegree = 0 ; mDegree <= degree + 1 ; mDegree++)
                {
                    double temp = normalMatrix[lDegree][mDegree];
                    normalMatrix[lDegree][mDegree] = normalMatrix[kDegree][mDegree];
                    normalMatrix[kDegree][mDegree] = temp;
                }
			}
		}
	}

    //loop to perform the gauss elimination
    for (lDegree = 0 ; lDegree < degree ; lDegree++)
    {
        for (kDegree = lDegree + 1 ; kDegree < degree + 1 ; kDegree++)
            {
                double matrixFactor = normalMatrix[kDegree][lDegree] / normalMatrix[lDegree][lDegree];
                    for (mDegree = 0 ; mDegree <= degree + 1 ; mDegree++)
                    {
                        normalMatrix[kDegree][mDegree] = normalMatrix[kDegree][mDegree] - matrixFactor * normalMatrix[lDegree][mDegree];   //make the elements below the pivot elements to zero or eliminate the variables
                    }
            }
    }

    //back-substitution
    for (kDegree = degree ; kDegree >= 0 ; kDegree--)
    {
        //make the variable to be calculated equal to the rhs of the last equation
        coefficents[kDegree] = normalMatrix[kDegree][degree + 1];
        for (lDegree = 0 ; lDegree < degree + 1 ; lDegree++)
        {
            if (lDegree != kDegree)
            {
                //then subtract all the lhs values except the coefficient of the variable whose value is being calculated
                coefficents[kDegree] = coefficents[kDegree] - normalMatrix[kDegree][lDegree] * coefficents[lDegree];
			}
		}
		//now finally divide the rhs by the coefficient of the variable to be calculated
        coefficents[kDegree] = coefficents[kDegree] / normalMatrix[kDegree][kDegree];
    }

}//END OF VOID polynomialFit
#if 0
By Jakob Mangold .... Did not get his Part of the Programm done in time and did not respond, so Lukas Samuel Jeck had to fix it
#endif // 0

//########################################### Minimum of Polynom ###########################################

#if 0
By Lukas Samuel Jeck [Idea came from Frederik Heberle to implement a dynamic display of the Intervall for MinCostFunction]
#endif // 0
//Counts FunctionalPumps
void countFunctionalPumps()
{
    //Counts Functionalpumps
	  for(int iPumps = 0; iPumps < NPUMP; iPumps++)
	  {
            if(boolMistake[iPumps] == false)
            {
                nFunctionalPumps++;
            }
	  }

}//END OF VOID countFunctionalPumps

//Routine to determine how many pumps are working and then saving them in another Array to use them sorted as the Intervall for MinCost
void determineIntervallMinCost(double* powerElSort)
{
	  //Saving functionalPumpEl's in the new Array
	  for(int iPumps = 0; iPumps < NPUMP; iPumps++)
	  {
	      if(boolMistake[iPumps] == false)
          {
            powerElSort[nSettingFunctionalPumpsArray] = pump.powerEl[iPumps];
            nSettingFunctionalPumpsArray++;
          }
	  }

}//END O VOID determineIntervallForMinCost

//Sorts an array[length] according to size of its values
void bubblesort(double* powerElSort, int length)
{
    //counting Variables
    int iSort, jSort;


    for (iSort = 0; iSort < length - 1; ++iSort)
    {

        for (jSort = 0; jSort < length - iSort - 1; ++jSort)
        {
            if (powerElSort[jSort] > powerElSort[jSort + 1])
            {
                double tmp = powerElSort[jSort];
                powerElSort[jSort] = powerElSort[jSort + 1];
                powerElSort[jSort + 1] = tmp;
            }
        }
    }

}//END OF VOID bubblesort
#if 0
By Lukas Samuel Jeck [Idea came from Frederik Heberle to implement a dynamic display of the Intervall for MinCostFunction and Improvement(determineIntervallMinCost-Function) after thinking about it]
#endif // 0




#if 0
By Frederik Heberle
#endif // 0
//Calculates values f(x) for any x
long double polyCalculation(double* coefficents, double minX, int degree)
{
	long double result = 0;
	for(int iPolynomial = 0; iPolynomial <= degree ; iPolynomial++)
    {
        result += coefficents[iPolynomial]*pow(minX,iPolynomial);
	}

	return result;

}//End OF polyCalculation

//Calculates a minimum inside of iBegin and iEnd
void minCostFunction()
{
	//minCost resembles y value of polynomial, xMin resembles x value

	double step;
	step = fabs(iBegin - iEnd) / precision;		//divide our range of values into (*precision)many parts
	xMin = iBegin;                              //set first x value for minimum

    //compute first value for minimum
	minCost = polyCalculation(coefficents, xMin, degree);

    //calculates the first step further
	double range = iBegin + step;

	//Calculating minCost and xMin
	while(range <= iEnd)
    {
        //check if value of one step further is lower than before
        if(polyCalculation(coefficents, range, degree) < minCost)
        {
            //saves the x and y value of that position if thats the case
            minCost = polyCalculation(coefficents, range, degree);
            xMin    = range;
        }

        //Afterwards jumping one step further
        range = range + step;
    }

}//END OF VOID minCostFunction
#if 0
By Frederik Heberle
#endif //0

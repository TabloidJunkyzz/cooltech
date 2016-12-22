
//++++++++++++++++++++++ "Cooltech.cpp" Implementation for our header "Cooltech.h" ++++++++++++++++++++++

//Includes Header which includes all necessary libraries/headers
#include "Cooltech.h"
using namespace std;


//Initialisation of our variables again to make them local [Commented on in "Cooltech.h"]

//File Operators and variables for read and output
ifstream file;
string skip;
double gelesen;

ofstream out;

//String for Versions
string version = "";

//Variables for functions of pipe measurements
double lambda;
double lambdaStart;
double lambdaTemp;

double newton;
double newtonHelp;

//Bool Variable mistake
bool mistake;

//Polynomfit
double pipeCostStorage;
double* coefficents;
int degree;

//Minimum calculation Variables
double xMin;
double minCost;

//Struct variables
data pipe, water, pump, cost, result;

//########################################### Reading Routine ###########################################

#if 0
By Lukas Samuel Jeck
#endif // 0
//Saves data from charts in arrays
void read()
{
	//Setup for Reading Pipe Data(Creating file operator)
file.open("Rohrdaten.txt");

        //Reads Flow of each pipe [First line of "Rohrdaten.txt"]; Skips end of line
        for(int i = 0;i < NPIPE;i++)
        {
            file>>gelesen;
            pipe.flow[i] = 0.001*gelesen;
                    assert(pipe.flow[i] > 0);
        }
	file >> skip;

        //Reads Pressureloss of each pipe [Second Line of "Rohrdaten.txt"]; Skips end of line
        for(int i = 0;i < NPIPE;i++)
        {
            file>>gelesen;
            pipe.deltaP[i] = 1e5*gelesen;
                    assert(pipe.deltaP[i] > 0);
        }
	file >> skip;

        //Reads Pipelength of each pipe [Third line of "Rohrdaten.txt"]; Skips end of line
        for(int i = 0;i < NPIPE;i++)
        {
            file>>gelesen;
            pipe.length[i] = gelesen;
                    assert(pipe.length[i] > 0);
        }
	file >> skip;

	/*
	Reads k [Fourth Line of "Rohrdaten.txt"]; Skips end of line
	Then closes file, cause every pipedata is in the system!
	*/
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
        for(int i = 0;i < NPUMP;i++)
        {
            file>>gelesen;
            pump.powerMech[i] = 1000*gelesen;
                    assert(pump.powerMech[i] > 0);
        }
	file >> skip;

        //Reading in Pumpefficiencies [Second line of "Pumpendaten.txt"]; skips end of line; Closes file
        for(int i = 0;i < NPUMP;i++)
        {
            file>>gelesen;
            pump.efficiency[i] = gelesen;
                    assert(pump.efficiency[i] > 0);
        }
	file >> skip;
file.close();


	// Setup for Reading Costparameters
file.open("Kosten.txt");

        //Reads DINs for naming available pipes [First line of "Kosten.txt"]; Skips end of line
        for(int i = 0; i < NCOST; i++)
        {
            file >> gelesen;
            cost.dIN[i] = gelesen;
                assert(cost.dIN[i] > 0);
        }

	file >> skip;

        //Reads Inside Diameter of available pipes [Second line of "Kosten.txt"]; Skips end of line
        for(int i = 0; i < NCOST; i++)
        {
            file >> gelesen;
            cost.insideDiameter[i] = gelesen/1000;
                    assert(cost.insideDiameter[i] > 0);
        }

	file >> skip;

        //Reads Wall Thickness of available pipes [Third line of "Kosten.txt"]; Skips end of line
        for(int i = 0; i < NCOST; i++)
        {
            file >> gelesen;
            cost.wallThickness[i] = gelesen/1000;
                    assert(cost.wallThickness[i] > 0);
        }

	file >> skip;

        //Reads CALCULATED outside Diameter by adding the inner one with the wall Thickness
        for(int i = 0;i < NCOST;i++)
        {
            cost.outsideDiameter[i] = cost.insideDiameter[i] + (2*cost.wallThickness[i]);
                    assert(cost.outsideDiameter[i] > 0);
        }

        // Reads CALCULATED powerEL with the Quotient of powerMech & efficiency
        for(int i = 0;i < NPUMP;i++)
        {
            pump.powerEL[i] = (pump.powerMech[i] / 1000) / pump.efficiency[i];
                    assert(pump.powerEL[i] > 0);
        }

	//Reads production Time [Fourth Line of "Kosten.txt"]; (Skips end of line); Closes file
	file >> gelesen;
	cost.operationTime = gelesen;
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

    cout <<endl<< "                  --- Pipe Data ---" <<endl<<endl;

        //Gives flow of each pipe in Console
        for(int i = 0;i < NPIPE;i++)
        {
            cout << "Pipe Flow " << i+1 << ": " << setw(7) << left << pipe.flow[i] << " [m^3/s]"<<endl;
        }
    cout<<endl;

        //Gives Deltap of each pipe in Console
        for(int i = 0;i < NPIPE;i++)
        {
            cout << "Delta P Pipe " << i+1 << ": " << setw(7) << left << pipe.deltaP[i] << " [Pa]"<<endl;
        }
    cout<<endl;

        //Gives length of each pipe in Console
        for(int i = 0;i < NPIPE;i++)
        {
            cout << "Pipe Length " << i+1 << ": " << setw(5) << left << pipe.length[i] << " [m]"<<endl;
        }
    cout<<endl;

    //Gives k of each pipe in Console
    cout << "Pipe Roughness:" << pipe.k / 1000 << " [m]" <<endl<<endl;


    cout<< "                  --- Water Data ---"<<endl<<endl;

    //Gives density and viscosity ind Console
    cout << "Density Water: " << water.density << " [kg/m^3]" <<endl<<endl;
    cout << "Viscosity Water: " << water.viscosity << " [Pa s]"<<endl<<endl;


    cout<<"                  --- Pump Data ---"<<endl<<endl;

        //Gives powerMech of each Pump in Console
        for(int i = 0;i < NPUMP;i++)
        {
            cout << "Mech. Power of Pump " << i+1 << ": " << setw(8) << left << pump.powerMech[i] << " [W]"<<endl;

        }
    cout<<endl;

        //Gives Pumpefficiency of each Pump in Console
        for(int i = 0;i < NPUMP;i++)
        {
            cout << "Efficiency for Pump " << i+1 << ": " << 100*pump.efficiency[i]
                 << " [%] (Read with  Factor 1/100)"<<endl;
        }
    cout<<endl;


    cout << "                  --- Cost Parameter ---" <<endl<<endl;

        //Gives insideDiameter of each available pipe in Console
        for(int i = 0; i < NCOST; i++)
        {
            cout << "Inside Diameter of available Pipes by DIN " << setw(3) << right << cost.dIN[i]
                 << ": " << setw(6) << left << cost.insideDiameter[i]<< " [m]"<<endl;
        }
    cout<<endl;

        //Gives wallThickness of each available pipe in Console
        for(int i = 0; i < NCOST; i++)
        {
            cout << "Wall Thickness of available Pipe by DIN " << setw(3) << right << cost.dIN[i] << ": " << setw(8) << left
                 << cost.wallThickness[i] << " [m]"<<endl;
        }
    cout<<endl;

        //Gives calculated outsideDiameter of each available pipe in Console
        for(int i = 0; i < NCOST; i++)
        {
            cout << "Outside Diameter of available Pipe by DIN " << setw(3) << right << cost.dIN[i] << ": "
                 << setw(8) << left << cost.outsideDiameter[i] << " [m]"<<endl;
        }
    cout<<endl;

        //Gives calculated PowerPEL of each Pump in Console
        for(int i = 0; i < NPUMP; i++)
        {
            cout << "PowerEL of Pump " << i+1 << ": " << setw(8) << left << pump.powerEL[i] << " [kW]"<<endl;
        }
    cout<<endl;

    //Gives operationTime in Console
    cout << "Operation Time: " << cost.operationTime << " [h]"<<endl<<endl;
}//END OF VOID printread




//Prints calculated data in Console
void printResults()
{
    //It's mentioned in the cout's what is printed below
    cout <<endl<< "		~~~ Consoleprint of Arrays that are calculated within the Functions ~~~"<<endl<<endl;
    cout << "		~~~ Rows are for Pumps and Columns are for Pipes ~~~"<<endl<<endl;

    cout << "--- PressureLossFriction ---"<<endl<<endl;
        for(int i = 0; i < NPUMP;i++)
        {
            for(int j = 0; j < NPIPE;j++)
            {
                cout << setw(12) << left << result.pressureLossFriction[i][j] << "[Pa]" << "    "; //i:Zeilen j:Spalten
            }
        cout<<endl<<endl;
        }

    cout << "--- Reynolds Number ---"<<endl<<endl;
        for(int i = 0; i < NPUMP;i++)
        {
            for(int j = 0; j < NPIPE;j++)
            {
                cout << setw(15) << left << result.reynolds[i][j] << "    "; //i:Zeilen j:Spalten
            }
        cout<<endl<<endl;
        }
    cout<<endl;

    cout << "############# Lambda Nikuradse: " << lambdaStart << " #############" <<endl<<endl;

    cout << "--- Calculated Diameter ---"<<endl<<endl;
        for(int i = 0; i < NPUMP;i++)
        {
            for(int j = 0; j < NPIPE;j++)
            {
                cout << setw(12) << left << result.diameter[i][j] << "[m]" << "    "; //i:Zeilen j:Spalten
            }
        cout<<endl<<endl;
        }
    cout<<endl;

    cout << "--- Ideal (inside)Diameters from Chart ---"<<endl<<endl;
        for(int i = 0; i < NPUMP;i++)
        {
            for(int j = 0; j < NPIPE;j++)
            {
                cout << setw(12) << left << result.insideDiameterChart[i][j] << "[m]" << "    "; //i:Zeilen j:Spalten
            }
        cout<<endl<<endl;
        }
    cout<<endl;

    cout << "--- Ideal (outside)Diameters (from Chart) ---"<<endl<<endl;
        for(int i = 0; i < NPUMP;i++)
        {
            for(int j = 0; j < NPIPE;j++)
            {
                cout << setw(12) << left << result.outsideDiameterChart[i][j] << "[m]" << "    "; //i:Zeilen j:Spalten
            }
        cout<<endl<<endl;
        }
    cout<<endl;

    cout << "--- Pressure Loss Valve ---"<<endl<<endl;
        for(int i = 0; i < NPUMP;i++)
        {
            for(int j = 0; j < NPIPE;j++)
            {
                cout << setw(12) << left << result.pressureLossValve[i][j] << "[Pa]" << "    "; //i:Zeilen j:Spalten
            }
        cout<<endl<<endl;
        }
    cout<<endl;

    cout << "--- Cost Pipes[Pipe] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            for(int j = 0;j < NPIPE;j++)
            {
                cout << "Cost Pipe[" << j+1 << "] " << result.pipeCost[i][j] << "    ";
            }
        cout<<endl<<endl;
        }
    cout<<endl;

    cout << "--- Cost Pumps[Pump] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            cout << "Cost Pump[" << i+1 << "] " << result.pumpCost[i]<<endl;
        }
	cout<<endl<<endl;


    cout << "--- Cost Power[Pump] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            cout << "Cost Power[" << i+1 << "] " << result.powerCost[i]<<endl;
        }
	cout<<endl<<endl;


    cout << "--- Total Costs[Pumpe] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            cout << "Total cost[" << i+1 << "]" << setw(6) << left << result.totalCost[i]<<endl;
        }
	cout<<endl<<endl;

    cout << "--- Minimum Calculation ---"<<endl<<endl;

    cout << "Least Cost: " << minCost << " " << "In regard to Pumppower: " << xMin << "[kW]"<<endl<<endl;

}//END OF VOID printResults




//Prints a text for the User in the Console
void userInformation()
{
    cout <<endl<<endl;
	cout << "Thank you for working with Cooltech_Industries"<<endl<<endl;
	cout << "Your Data is being Considered and Calculations are in Queue ..."<<endl<<endl;
	cout << "Pls open Cooltech_Data.txt to look at the Calculation Results"<<endl<<endl;

	if(mistake == true)
        {
            cout << "The Goddamn first Pump isn´t Working :( !"<<endl<<endl;
        }

}//END OF VOID userInformation




//Takes calculated data and writes it in Cooltech_Daten.txt
void output()
{
    //Setting up output(Creating file operator)
out.open("Cooltech_Daten.txt");

        if(mistake == true)
            {
            out << "The Goddamn first Pump isn't Working :( ! --> Consider the results Shown for the first Pump as False !"<<endl<<endl;
            }


	//Information (Headline of "Cooltech_Daten.txt")
	out << endl;
	out << "Group 1.5 Numerical Project" << "            --- Writing results into a .txt File ---";
	out << endl<<endl<<endl;


	//Defining what data is shown
	out << "~~~ Data of Pipecalculations ~~~"<<endl<<endl<<endl;


	//Putting out ideal Diameters of pipes for each Pipe
	out << "--- Ideal InsideDiameter[Pump][Pipe] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            for(int j = 0;j < NPIPE;j++)
            {
                    out << "InsideDiameter[" << i+1 << "]" << "[" << j+1 << "] " << setw(6) << left << result.insideDiameterChart[i][j] << " [m]" << "    ";
            }
        out<<endl;
        }
	out<<endl;

	//Putting out Pressureloss of Valve
	out << "--- Pressureloss Valve[Pump][Pipe] ---"<<endl<<endl;

        for(int i=0;i<NPUMP;i++)
        {
            for(int j=0;j<NPIPE;j++)
            {
	      		out << "Valve Pressureloss[" << i+1 << "]" << "[" << j+1 << "] " << setw(7) << left << result.pressureLossValve[i][j] << " [Pa]" << "    ";
            }
		out<<endl;
        }
	out<<endl<<endl;


	//Data that is getting put out is shown above
	out << "~~~ Data of Costcalculations ~~~"<<endl<<endl<<endl;


	out << "--- Cost Pipes[Pipe] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            for(int j = 0;j < NPIPE;j++)
            {
                out << "Cost Pipe[" << j+1 << "] " << setw(7) << left << result.pipeCost[i][j] << " [€]" << "    ";
            }
        out<<endl;
        }
    out<<endl<<endl;

	out << "--- Cost Pumps[Pump] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            out << "Cost Pump[" << i+1 << "] " << setw(8) << left << result.pumpCost[i] << "[€]"<<endl;
        }
	out<<endl<<endl;

	out << "--- Cost Power[Pump] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            out << "Cost Power[" << i+1 << "] " << setw(8) << left << result.powerCost[i]<< "[€]"<<endl;
        }
	out<<endl<<endl;

	out << "--- Total Costs[Pumpe] ---"<<endl<<endl;
        for(int i = 0;i < NPUMP;i++)
        {
            out << "Total Cost[" << i+1 << "] " << setw(7) << left << result.totalCost[i]<< "[€]"<<endl;
        }
	out<<endl<<endl;

out.close();
} //END of Void output




//Checks PressurelossFriction Array for negative Values (WE KNOW THAT THEY ARE FROM PUMP 1)
void checkForMistake(int i ,int j)
{
    if(result.pressureLossFriction[i][j]<0)
        {
            mistake = true;
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
double pressureLossFrictionFunction (int i, int j)
{
	return (pump.powerMech[i]/TOTALFLOW)-pipe.deltaP[j];

}//END OF double pressureLossFrictionFunction




//Calculates a value for lambdaStart [Initialisation of double lambdaStart in main]
double nikuradse(double k)
{
    double D_start = 0.5; // assume the start diameter is 0.5[m]
	return 1 / pow((2 * log10(3.71 * D_start / pipe.k)), 2.0);

}//END OF double nikuradse




//Calculates Diameter 1DIM [Initialisation of an array in main]
double diameterFunction (int i, int j)
{
	return pow(((8.0*lambda*pipe.length[j]*water.density*pipe.flow[j]*pipe.flow[j])/(M_PI*M_PI*result.pressureLossFriction[i][j])),0.2);

}//END OF double diameterFunction




//Calculates reynolds number 1DIM [Initialisation of an array in main]
double reynoldsFunction (int i, int j)
{
	return (4*pipe.flow[j]*water.density)/(water.viscosity*diameterFunction(i,j)*M_PI);

}//END OF double reynoldsFunction




//NEWTON PROCESS Declarations

//Colebrookfunction
double colebrookFunction (int i, int j, double lambda)
{
    return 1.74 - 2.0*log10((2.0*pipe.k/diameterFunction(i,j))+18.7/(reynoldsFunction(i,j)*sqrt(lambda))) - 1/sqrt(lambda);

}//END OF colebrookfunction



//Derivating Colebrook numerically
double colebrookDerivativeFunction (int i, int j, std::function<double(int, int, double)>func)
{
       return (func(i, j, lambda+EPSILON) - func(i ,j ,lambda-EPSILON))/(2*EPSILON);

}//END OF colebrookDerivativeFunction



//Solves Colebrook for lambda
double newtonFunction (int i, int j, double lambda, std::function<double(int,int,double)>func)
{
    do
        {
            newtonHelp  = lambda;
            newton      = lambda - (func(i, j, lambda)/colebrookDerivativeFunction(i, j, colebrookFunction));
            lambda      = newton;

        } while (fabs(newton - newtonHelp) > EPSILON);

    return lambda;

}//END OF double newtonFunction



//Picks insideDiameter from chart, Calculates outsideDiameter with insideDiameter from chart and corresponding Wall Thicknes
void pickDiameterFromChart ()
{
    for(int i = 0; i < NPUMP; i++)
    {
        for(int j = 0; j < NPIPE; j++)
        {

            int m = 0;

            while(result.diameter[i][j] > cost.insideDiameter[m])//go through all the data of insideDiameter until you find a diameter bigger than the calculated one
            {
                m++;
                	assert(m < NCOST);
            }

            result.insideDiameterChart[i][j]        = cost.insideDiameter[m]; //save the picked Diameter in Array

            result.outsideDiameterChart[i][j]       = cost.insideDiameter[m] + 2*cost.wallThickness[m]; //Calculate outsideDiameter
        }
    }
}//END OF VOID pickDiameterFromChart
#if 0
By Hermann Gegel & Boda Yang [Lukas Samuel Jeck as a Helping Hand]
#endif // 0




#if 0
By Lina Lepp
#endif // 0
//Calculates the, to be regulated, pressure Loss by the Valve
double pressureLossValveFunction(int i, int j)
{
	return (8*lambda*pipe.length[j]*water.density*pipe.flow[j]*pipe.flow[j]) / (M_PI*M_PI*pow(result.diameter[i][j],5.0)) - ((8*lambda*pipe.length[j]*water.density*pipe.flow[j]*pipe.flow[j]) / (M_PI*M_PI*pow(result.insideDiameterChart[i][j],5.0)));

}//END OF double pressureLossValveFunction
#if 0
By Lina Lepp
#endif // 0

//########################################### Cost Calculations ###########################################

#if 0
By Jakob Mangold
#endif // 0
//Calculates Pumpcost
double pumpCostFunction(int i)
{
    return pump.powerEL[i]*406 + 4011*(5 - pump.efficiency[i]);

}//END OF pumpCostFunction




//Calculates pipeCost
double pipeCostFunction(int i, int j)
{
	return (result.outsideDiameterChart[i][j]*result.outsideDiameterChart[i][j]*16458 - result.outsideDiameterChart[i][j]*2109 + 151)*pipe.length[j];

}//END OF pipeCostFunction




////Calculates powerCost
double powerCostFunction(int i)
{
	 return cost.operationTime*pump.powerEL[i]*0.05;

}//END OF powerCostFunction




////Calculates totalCost
double totalCostFunction(int i)
{
    pipeCostStorage = 0;
        for (int j = 0; j < NPIPE; j++)
        {
            pipeCostStorage += pipeCostFunction(i,j);
        }
	return pumpCostFunction(i) + pipeCostStorage + powerCostFunction(i);

}//END OF totalCostFunction

//########################################### Polynomial Fit ###########################################

//PonynomialFit for our Total Cost Function
void polynomialFit()	                    //add cin>> degree for user														//degree is the degree of the polynom, x[] and y[] the valuepairs
{                                                       // ATTENTION: ADD ASSERT FOR DEGREE > NPUMP ---> ABORT
    int i, k, l, m;										//counting Variables
    coefficents = new double[degree+1];	                //Array that stores the coefficents of the polynom
    for (k = 0; k < degree; k++)
    {
        coefficents[k] = 0;
    }
    double powerElArray[2 * degree + 1];							//Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    for (k = 0 ; k < 2 * degree + 1 ; k++)
    {
        powerElArray[k]=0;
        for (i = 0 ; i < NPUMP ; i++)
        {
			powerElArray[i] = powerElArray[i] + pow(pump.powerEL[i],k);
		}
    }
    double normalMatrix[degree + 1][degree + 2];					//B is the Normal matrix that will store the equations
    for (k = 0 ; k <= degree ; k++)
    {
        for (l = 0 ; l <= degree ; l++)
        {
            normalMatrix[k][l] = powerElArray[k + l];
        }
    }													//Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
    double resultArray[degree + 1];                    			//Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    for (k = 0 ; k < degree + 1 ; k++)
    {
        resultArray[i] = 0;
        for (i = 0 ; i < NPUMP ; i++)
        {
			resultArray[k] = resultArray[k] + pow(pump.powerEL[i],k) * result.totalCost[i];
		}
    }
    for (k = 0 ; k <= degree ; k++)
    {
        normalMatrix[k][degree+1] = resultArray[k];
    }               									//load the values of Y as the last column of B


    for (l = 0 ; l < degree + 1 ; l++)                 	//From now Gaussian Elimination starts to solve the set of linear equations (Pivotisation)
     {
        for (k = i + 1 ; k < degree + 1 ; k++)
        {
            if (normalMatrix[l][l] < normalMatrix[k][l])
            {
                for (m = 0 ; m <= degree + 1 ; m++)
                {
                    double temp = normalMatrix[l][m];
                    normalMatrix[l][m] = normalMatrix[k][m];
                    normalMatrix[k][m] = temp;
                }
			}
		}
	}
    for (l = 0 ; l < degree ; l++)
    {         											//loop to perform the gauss elimination
        for (k = i + 1 ; k < degree + 1 ; k++)
            {
                double t = normalMatrix[k][l] / normalMatrix[l][l];
                for (m = 0 ; m <= degree + 1 ; m++)
                {
                    normalMatrix[k][m] = normalMatrix[k][m] - t * normalMatrix[l][m];    //make the elements below the pivot elements
				}
            }
    }
    for (k = degree ; k >= 0 ; k--)               		//back-substitution
    {
        coefficents[k] = normalMatrix[k][degree + 1];                		//make the variable to be calculated equal to the rhs of the last equation
        for (l = 0 ; l < degree + 1 ; l++)
        {
            if (l != k)
            {            								//then subtract all the lhs values except the coefficient of the variable whose value is being calculated
                coefficents[k] = coefficents[k] - normalMatrix[k][l] * coefficents[l];
			}
		}
        coefficents[k] = coefficents[k] / normalMatrix[k][k];           				//now finally divide the rhs by the coefficient of the variable to be calculated
    }

}//END OF VOID polynomialFit
#if 0
By Jakob Mangold
#endif // 0

//########################################### Minimum of Polynom ###########################################

#if 0
By Frederik Heberle
#endif // 0
//Calculates the Minimum of the Polynom
void  minCostFunction()
{
	double tempMinCost = 0.0;
    minCost = 1e10;                                    //setting initial minCost value high to ensure that it will be replaced later
	double firstGuess;										//first value to start numerical evaluation later taken from powerEL

	double xValues[NPUMP];                                  //creating a new array to sort
	for(int i = 0; i < NPUMP; i++) {                        //copying all necessary elements of powerEL to new array
        xValues[i] = pump.powerEL[i];
	}

	for(int j = xValues[0]; j < xValues[NPUMP-1]; j++) {        //starting brute-force search for minimum at lowest value on power-scale of our pumps !!NOTIZ WENN xVALUES double dann muss j double!!
		firstGuess = j;                            			//firstGuess as variable to represent x in our polynomial

		for(int i = 0; i < degree; i++) {			//using Horner's scheme to compute result of polynomial
			tempMinCost += pow(firstGuess, i) * coefficents[i];      //getting tempMinCost as temp-variable for our computed y-value of the polynomial
		}

		if(minCost >= tempMinCost && tempMinCost > 0) {                        //comparing our newly computed y-value to the last given y-value
			minCost = tempMinCost;                          //if newer y-value is lower (cheaper) than one before, overwrite old value with new one
			xMin = firstGuess;                              //save corresponding x-value to new y-value for output
		}
	}

}//END OF VOID minCostFucntion
#if 0
By Frederik Heberle
#endif // 0

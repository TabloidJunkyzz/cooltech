
//++++++++++++++++++++++ "Cooltech_main.cpp" [Hauptprogramm #unfinished] ++++++++++++++++++++++

//Includes Header which includes all necessary libraries/headers
#include "Cooltech.h"
using namespace std;


//Here Happens the good Stuff
int main()
{

//########################################### General Information ###########################################

    //Question for User or Dev-Version
    cout << "############# Launch in Developer - Version? #############"<<endl<< "Yes[y], No[n]" <<endl;
    //Expects to get a letter in Console ---> Throws it into "version"
    getline(cin,version);
    cout<<endl<<endl;

//########################################### Reading Routine ###########################################

    //Puts Data in System
    read();

//############################################### Pipe Dimensioning #######################################

	//Here happens most of the stuff: 1.) Calculating our needed datas  2.) Saving those in arrays with [NPUMP][NPIPE]
    for (int i = 0;i < NPUMP;i++)
    {
        for(int j = 0;j < NPIPE;j++)
        {
            result.pressureLossFriction[i][j] = pressureLossFrictionFunction(i, j);
            lambdaStart                       = nikuradse(pipe.k);
            lambda                            = lambdaStart;

                //Recalculating diameter and reynolds with newly aquired lambda as long as these newly aquired and temporarily saved lambdas are small enough in difference
            do
            {
                lambdaTemp                      = lambda;
                result.diameter[i][j] 	        = diameterFunction(i, j);
                result.reynolds[i][j]           = reynoldsFunction(i, j);
                lambda                          = newtonFunction(i, j, lambda, colebrookFunction) ;

            } while (fabs(lambdaTemp - lambda) > 0.0001);
        }
    }

    //Picks Diametere from Chart --> Saved in "result.insideDiameterChart[i][j]" & "result.outsideDiameterChart[i][j]"
    pickDiameterFromChart();

    //Loop for saving pipeCost & pressureLossValve
    for (int i = 0;i < NPUMP;i++)
	{
		for(int j = 0;j < NPIPE;j++)
		{
            result.pressureLossValve[i][j] = pressureLossValveFunction(i, j);
		}
	}

//########################################### Cost Functions ###########################################

    //loop for saving pipeCost
    for (int i = 0;i < NPUMP;i++)
        {
            for(int j = 0;j < NPIPE;j++)
            {
                result.pipeCost[i][j] = pipeCostFunction(i, j);
            }
        }

	//Loop for saving pumpCost and powerCost, then Calculating totalCost
    for(int i = 0;i < NPUMP;i++)
        {
            result.pumpCost [i]  = pumpCostFunction(i);
            result.powerCost[i]  = powerCostFunction(i);
            result.totalCost[i]  = totalCostFunction(i);
        }

//########################################### Polynomial Fit & Minimum Cost ###########################################

    //Asks for a a degree to polynomialfit !! Degree cant be higher than NPUMP !!
	cout << "Please enter degree of the polynom (has to be smaller than the amount of pumps): " << endl;
	cin >> degree;
        assert(degree < NPUMP);
    cout<<endl<<endl;

    //Does Polnomialfit and Calculates Minimum
    polynomialFit();
    minCostFunction();

//########################################### Results and Console ###########################################

     //Prints out Usertext in Console [User-Version]
	if(version == "n")
    {
        userInformation();
    }

    //Prints read Values in console [Dev-Version]
    if(version == "y")
        {
            printRead();
        }

    //Prints Results of various calculations in Console [Dev-Version]
    if (version == "y")
    {
        printResults();
    }

    //Puts out results in "Cooltech_Daten.txt"
    output();

    //deletes Pointer cause coefficients is a dynamic Field!
    delete [] coefficents;

} //INT MAIN END

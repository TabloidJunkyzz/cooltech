
//++++++++++++++++++++++ "Cooltech_main.cpp" ++++++++++++++++++++++

//Includes Header which includes all necessary libraries/headers
#include "Cooltech.h"
using namespace std;


//All the Tasks happen here! Every function we have, can be called and all the Goodness can be done!!
int main()
{

//########################################### General Information ###########################################

    //Question for User or Dev-Version
    cout << "############# Launch in Developer - Version? #############"<<endl<< "Yes[y], No[n]" <<endl;

    //Expects to get a letter in Console ---> Throws it into "version"
    //If something else than "n" or "y" is entered it just asks again
    while(version != "n" && version != "y")
    {
        getline(cin,version);
    }
    cout<<endl<<endl;

//########################################### Reading Routine ###########################################

    //Puts Data in System
    read();

//########################################### Pipe Dimensioning ###########################################

	//Here happens most of the stuff: 1.) Calculating our needed datas  2.) Saving those in arrays with [NPUMP][NPIPE]
    for (int iPumps = 0;iPumps < NPUMP;iPumps++)
    {
        for(int iPipes = 0;iPipes < NPIPE;iPipes++)
        {
            result.pressureLossFriction[iPumps][iPipes] = pressureLossFrictionFunction(iPumps, iPipes);
            lambdaStart                                 = nikuradse(pipe.k);
            lambda                                      = lambdaStart;

                //Recalculating diameter and reynolds with newly aquired lambda as long as these newly aquired and temporarily saved lambdas are small enough in difference
            do
            {
                lambdaTemp                      = lambda;
                result.diameter[iPumps][iPipes] = diameterFunction(iPumps, iPipes);
                result.reynolds[iPumps][iPipes] = reynoldsFunction(iPumps, iPipes);
                lambda                          = newtonFunction(iPumps, iPipes, lambda, colebrookFunction) ;

            } while (fabs(lambdaTemp - lambda) > 0.0001);
        }
    }

    //Picks Diametere from Chart --> Saved in "result.insideDiameterChart[iPumps][iPipes]" & "result.outsideDiameterChart[iPumps][iPipes]"
    pickDiameterFromChart();

    //Loop for saving pipeCost & pressureLossValve
    for (int iPumps = 0;iPumps < NPUMP;iPumps++)
	{
		for(int iPipes = 0;iPipes < NPIPE;iPipes++)
		{
            result.pressureLossValve[iPumps][iPipes] = pressureLossValveFunction(iPumps, iPipes);
		}
	}

//########################################### Cost Functions ###########################################

    //loop for saving pipeCost
    for (int iPumps = 0;iPumps < NPUMP;iPumps++)
    {
        for(int iPipes = 0;iPipes < NPIPE;iPipes++)
        {
            result.pipeCost[iPumps][iPipes] = pipeCostFunction(iPumps, iPipes);
        }
    }

	//Loop for saving pumpCost and powerCost, then Calculating totalCost
    for(int iPumps = 0;iPumps < NPUMP;iPumps++)
        {
            result.pumpCost [iPumps]  = pumpCostFunction(iPumps);
            result.powerCost[iPumps]  = powerCostFunction(iPumps);
            result.totalCost[iPumps]  = totalCostFunction(iPumps);
        }

    //Checks pressureLossFriction Array for Mistakes
    checkForMistake();

//########################################### Polynomial Fit ###########################################

    //Asks for a a degree to polynomialfit !! Degree cant be higher than NPUMP !!
	cout << "Please enter the degree of the polynom (has to be smaller than the amount of pumps): " << endl;

	//Stops the User from getting "wrong" values in the system
	do
    {
       	cin  >> degree;

    } while(degree < 0 || degree > NPUMP);
    cout <<endl<<endl;

    //Polnomialfits
    polynomialFit();

//########################################### Minimum Cost ###########################################

    countFunctionalPumps();

    double powerElSort[nFunctionalPumps]; //Creating New Array which gets set after and sorted later

    determineIntervallMinCost(powerElSort);

    bubblesort(powerElSort, nFunctionalPumps);

    //User Information and input for minCost
    cout << "Please enter the range in which the calculation of a minimum should take place (cannot be lower/higher than lowest and highest pump values)." << endl;
    cout << "Lowest pump value = " << powerElSort[0] << ", highest pump value = "<< powerElSort[nFunctionalPumps-1]<<endl<<endl;

    //Stops the User from getting "wrong" values in the system
    do
    {
        cout << "Please enter the lower end of your desired range: ";
        cin  >> iBegin;

    } while(iBegin < powerElSort[0] || iBegin > powerElSort[nFunctionalPumps-1]);
    cout <<endl;

    //Stops the User from getting "wrong" values in the system
    do
    {
        cout << "Please enter the high end of your desired range: ";
        cin  >> iEnd;

    } while(iEnd > powerElSort[nFunctionalPumps-1] || iEnd < iBegin);
    cout <<endl;

    cout << "Please enter a value for the precision of our minimum calculation"<<endl;
    cout << "(higher values for precision result in more accurate calculations it is also advisable to choose a value greater than 1 or equal to it):";
    cout <<endl;
    do
    {
        cin  >> precision;

    } while(precision < 1.0);
    cout <<endl<<endl;

    //Calculates Minimum
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
        printResults(powerElSort);
    }

    //Puts out results in "Cooltech_Daten.txt"
    output();

    //deletes Pointer
    delete [] coefficents;

} //INT MAIN END

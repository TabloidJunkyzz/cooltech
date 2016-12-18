//"Cooltech_main.cpp" [Hauptprogramm #unfinished] : --- Doing Stuff ---

#include "Cooltech.h"
using namespace std;


#if 0
+Kostenfunktion aufrufen fehlt noch
#endif

int main(){

cout << "############# Wollen sie die Entwicklerversion starten?? #############"<<endl<< "Ja[j], Nein[n]" <<endl;
getline(cin,version);

//Puts Data in System
read();
//Prints Data in console

if(version == "j")
    {
        printRead();
    }

//for loop for saving pumpCost
for(int i=0;i<NPUMP;i++)
    {
      result.pumpCost[i] = pumpCostFunction(i);
    }

	//Here happens most of the stuff: 1.) Calculating our needed datas  2.) Saving those in arrays with [NPUMP][NPIPE]
	for (int i=0;i<NPUMP;i++)
	{
		for(int j=0;j<NPIPE;j++)
		{
			result.pressureLossFriction[i][j] = pressureLossFrictionFunction(i,j);
			lambdaStart = nikuradse(pipe.k);
			lambda =lambdaStart;

				do
				{
					lambdaTemp = lambda;
					result.diameter[i][j] 	        = diameterFunction(i,j);
					result.reynolds[i][j]           = reynoldsFunction(i,j);
					lambda = newtonFunction(i,j, lambda, colebrookFunction) ;//Functioncalling of newton

				} while (fabs(lambdaTemp - lambda) > 0.0001);

			result.pressureLossValve[i][j] = pressureLossValveFunction(i,j);
			result.pipeCost[j]             = pipeCostFunction(i,j);
			result.totalCost[i][j]         = totalCostFunction(i,j);
			checkForMistake(i,j);

		}
	}

	if(version == "n")
    {
        userInformation();
    }

	pickDiameterFromChart(); //Saved in "result.insideDiameterChart[i][j]" & "result.outsideDiameterChart[i][j]"

if (version == "j")
{
	printResults();
}
	output();


	cout << "Please enter degree of the polynom (has to be smaller than the amount of pumps): " << endl;
	cin >> degree;
	assert(degree < NPUMP);


    polynomialFit();


} //INT MAIN END

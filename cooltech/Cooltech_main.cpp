//"Cooltech_main.cpp" [Hauptprogramm #unfinished] : --- Doing Stuff ---

#include "Cooltech.h"
using namespace std;


#if 0
+Kostenfunktion aufrufen fehlt noch
#endif

int main(){

cout << "############# Wollen sie die Entwicklerversion starten?? #############"<<endl<< "Ja[j], Nein[n]" <<endl;
getline(cin,version);

if(version == "n")
{
    userInformation();
}
	//Puts Data in System
	read();
	//Prints Data in console

if(version == "j")
{
	printRead();
}
	//Here we have a for-loop to save data in arrays with [NPUMP]
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
			result.pipeCost[i][j] 		   = pipeCostFunction(i,j);
			result.totalCost[i][j]         = totalCostFunction(i,j);

		}
	}

	pickDiameterFromChart(); //Saved in "result.insideDiameterChart[i][j]" & "result.outsideDiameterChart[i][j]"

if (version == "j")
{
	printResults();
}
	output();

} //INT MAIN END

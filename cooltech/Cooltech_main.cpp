//"Cooltech_main.cpp" [Hauptprogramm #unfinished] : --- Doing Stuff ---

#include "Cooltech.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <functional>
#include <math.h>
using namespace std;


int main(){

	//Puts Data in System
	read();
	//Prints Data in console
	printRead();


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

            //cout <<lambda <<" zu kombi " <<i <<j <<endl;
			result.pressureLossValve[i][j] = pressureLossValveFunction(i,j);
			result.totalCost[i][j]         = totalCostFunction(i,j);
		}
	}

	pickDiameterFromChart();//Saved in result.insideDiameterChart[i][j]

	printResults();

	output();

} //INT MAIN END

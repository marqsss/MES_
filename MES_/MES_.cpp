#include <iostream>
#include <armadillo>
#include <random>
#include "MESGrid.h"
#include "MESCalc.h"

int main()
{
	int choice = 1;
	arma::mat globH;
	arma::mat globC;
	arma::mat temp;
	mes::Grid grid;
	mes::Calc calc;
	while (choice != 0)
	{
		grid.reset();
		std::cout << "0: exit\n1: excel test\n2: pdf test (test case)" << std::endl;
		std::cin >> choice;
		switch (choice)
		{
		case 0:
			break;
		case 1:

			grid.loadFromFile("testGrid.txt");
			grid.print();

			calc.getLocalHMatrix(grid, 0).print("H:");

			calc.getLocalCMatrix(grid, 0).print("Macierz C:");

			calc.getHBCMatrix(grid, 24).print("Macierz HBC:");

			calc.getPVector(grid, 24).print("P vector:");

			break;
		case 2:
			grid.loadFromFile("TestCaseGrid.txt");
			calc.printExtremeTemp(grid); // step 0

			for (unsigned int i = 0; i < grid.getTime() / grid.getDeltaTau(); i++)
			{
				calc.applyGauss(grid);
				calc.printExtremeTemp(grid);
			}

			break;
		}
	}



	//system("pause");
	return 0;
}
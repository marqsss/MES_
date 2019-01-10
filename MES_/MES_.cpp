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

			calc.getLocalCMatrix(grid, 0, true).print("Macierz C:");

			calc.getHBCMatrix(grid, 24).print("Macierz HBC:");

			calc.getPVector(grid, 24, true).print("P vector:");

			break;
		case 2:
			grid.loadFromFile("TestCaseGrid.txt");
			grid.print();
			globH = calc.getGlobalMatrix(grid, true);
			globH.print("Global H:");

			globC = calc.getGlobalMatrix(grid, false);
			globC.print("Global C:");
			std::cout << "Max temp: " << calc.getMaxTemp(grid) << ", min temp: " << calc.getMinTemp(grid) << std::endl;
			break;
		}
	}



	//system("pause");
	return 0;
}
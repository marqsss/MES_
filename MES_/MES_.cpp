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

			temp = calc.getLocalHMatrix(grid, 0);
			for (unsigned int i = 0; i < temp.n_cols; i++)
			{
				for (unsigned int j = 0; j < temp.n_rows; j++)
				{
					std::cout << "\t" << temp(i, j);
				}
				std::cout << std::endl;
			}
			temp.print("H:");
			calc.getLocalCMatrix(grid, 0).print("Macierz C:");

			globH = calc.getGlobalMatrix(grid);

			globC = calc.getGlobalMatrix(grid, false);
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
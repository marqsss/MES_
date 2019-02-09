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
		std::cout << "0: exit\n1: pdf test case 1\n2: pdf test case 2" << std::endl;
		std::cin >> choice;
		switch (choice)
		{
		case 0:
			break;
		case 1:
			grid.loadFromFile("TestCaseGrid.txt");
			std::cout << "#0: ";
			calc.printExtremeTemp(grid); // step 0

			for (unsigned int i = 0; i < grid.getTime() / grid.getDeltaTau(); i++)
			{
				calc.applyGauss(grid);
				std::cout << "#"<<i+1<<": ";
				calc.printExtremeTemp(grid);
			}
			break;
		case 2:
			grid.loadFromFile("TestCase_2.txt");
			std::cout << "#0: ";
			calc.printExtremeTemp(grid); // step 0

			for (unsigned int i = 0; i < grid.getTime() / grid.getDeltaTau(); i++)
			{
				calc.applyGauss(grid);
				std::cout << "#"<<i+1<<": ";
				calc.printExtremeTemp(grid);
			}
			break;
		}
	}



	//system("pause");
	return 0;
}
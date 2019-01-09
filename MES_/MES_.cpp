#include <iostream>
#include <armadillo>
#include <random>
#include "MESGrid.h"
#include "MESCalc.h"

int main()
{
	mes::Grid grid("testGrid.txt");
	grid.print();

	mes::Calc calc;
	
	arma::mat temp = calc.getLocalHMatrix(grid, 0);
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
	
	arma::mat globH = calc.getGlobalMatrix(grid);

	arma::mat globC = calc.getGlobalMatrix(grid, false);



	system("pause");
	return 0;
}
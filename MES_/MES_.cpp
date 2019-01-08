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
	calc.getLocalHMatrix(grid, 0).print("Macierz H:");
	calc.getLocalCMatrix(grid, 0).print("Macierz C:");
	calc.getGlobalHMatrix(grid, true).print("Global H:");


	system("pause");
	return 0;
}
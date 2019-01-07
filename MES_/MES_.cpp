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
	//calc.getLocalHMatrix(grid.getElement(0)->getK(), arma::dvec("0, 0.025, 0.025, 0"), arma::dvec("0, 0, 0.025, 0.025")).print("Macierz H:");
	calc.getLocalHMatrix(grid, 0).print("Macierz H:");
	//calc.getLocalCMatrix(arma::dvec("0, 0.025, 0.025, 0"), arma::dvec("0, 0, 0.025, 0.025"), 700, 7800).print("Macierz C:"); 
	calc.getLocalCMatrix(grid, 0).print("Macierz C:");


	system("pause");
	return 0;
}
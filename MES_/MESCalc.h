#ifndef MESCALC_H
#define MESCALC_H

#include "MESGrid.h"

namespace mes
{
	enum MatrixType {
		H,
		C,
		HBC,
		H_central // bez warunkow brzegowych
	};

	class Calc
	{
	public:
		Calc();
		//Calc(Grid& GRID);
		static const double ksi;
		static const double eta;
		static const arma::dvec etaCol;
		static const arma::dvec ksiCol;
		arma::dvec transformKsiEta(const arma::dvec &pos);
		arma::mat getDNDKsi();
		arma::mat getDNDEta();
		arma::mat getJacobian(const arma::dvec &x, const arma::dvec &y);
		arma::mat getLocalHMatrix(double k, const arma::dvec &x, const arma::dvec &y, bool debug = false);
		arma::mat getLocalHMatrix(Grid& GRID, unsigned int index, bool debug = false);
		arma::mat getLocalHMatrix(Element& e, bool debug = false);
		arma::mat getDetJ(const arma::mat &Jac);
		arma::mat getJacInv(const arma::mat &Jac, const arma::mat &detJ);
		arma::mat getDndy(const arma::mat &JacInv, const arma::mat &dKsi, const arma::mat &dEta);
		arma::mat getDndx(const arma::mat &JacInv, const arma::mat &dKsi, const arma::mat &dEta);
		arma::mat getLocalCMatrix(const arma::dvec &x, const arma::dvec &y, double c, double ro, bool debug = false);
		arma::mat getLocalCMatrix(Grid& GRID, unsigned int index, bool debug = false);
		// HorC: true for H_matrix, false for C_matrix; default true
		arma::mat getGlobalMatrix(Grid& grid, MatrixType type, bool debug = false);
		arma::mat getHBCMatrix(Grid& grid, unsigned int index, bool debug = false);
		arma::dvec getPVector(Grid& grid, unsigned int index, bool debug = false);
		arma::dvec getGlobalPVector(Grid& grid, bool debug = false);
		arma::mat getHCdTMatrix(Grid& grid, bool debug = false);
		arma::dvec getTemperaturesVector(Grid& grid, bool debug = false);
		arma::dvec dummy(Grid& grid, unsigned int index);
		void dummy2(Grid& grid, unsigned int index);
		arma::dvec gauss(Grid& grid, bool debug = false);
		void applyGauss(Grid& grid, unsigned int iterations = 1, bool debug = false);
		void applyGaussImproved(Grid& grid, unsigned int iterations = 1, bool debug = false); // does not work.
		arma::mat getHCdTGP(Grid& grid, bool debug = false);
		void printExtremeTemp(Grid& grid, unsigned int iterations = 0);


		double N1(double a, double b) { return 0.25*((1 - a)*(1 - b)); }
		double N2(double a, double b) { return 0.25*((1 + a)*(1 - b)); }
		double N3(double a, double b) { return 0.25*((1 + a)*(1 + b)); }
		double N4(double a, double b) { return 0.25*((1 - a)*(1 + b)); }

	private:
		//mes::Grid& grid;
		arma::mat nMatrix;
	};
}

#endif
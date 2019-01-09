#ifndef MESCALC_H
#define MESCALC_H

#include "MESGrid.h"

namespace mes
{
	class Calc
	{
	public:
		Calc();
		//Calc(Grid& GRID);
		static double ksi;
		static double eta;
		static arma::dvec etaCol;
		static arma::dvec ksiCol;
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
		arma::mat getGlobalMatrix(Grid& grid, bool HorC = true, bool debug = false);
		double getMinTemp(Grid& grid);
		double getMaxTemp(Grid& grid);

	private:
		//mes::Grid& grid;
		arma::mat nMatrix;
	};
}

#endif
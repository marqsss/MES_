#include <armadillo>
#include <string>
#include "MESCalc.h"

double mes::Calc::ksi = 1.0 / sqrt(3);
double mes::Calc::eta = 1.0 / sqrt(3);
arma::dvec mes::Calc::ksiCol = {-1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3), -1.0 / sqrt(3)};
arma::dvec mes::Calc::etaCol = {-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3)};

//mes::Calc::Calc(Grid& GRID): grid(GRID)
mes::Calc::Calc()
{
	// -- +- ++ -+
	nMatrix.insert_cols(0, (1 - ksiCol) % (1 - etaCol));
	nMatrix.insert_cols(1, (1 + ksiCol) % (1 - etaCol));
	nMatrix.insert_cols(2, (1 + ksiCol) % (1 + etaCol));
	nMatrix.insert_cols(3, (1 - ksiCol) % (1 + etaCol));
	nMatrix *= 0.25;
}

arma::dvec mes::Calc::transformKsiEta(const arma::dvec &pos)
{
	arma::dvec res = {arma::sum(pos%nMatrix.row(0)), arma::sum(pos%nMatrix.row(1)),
		arma::sum(pos%nMatrix.row(2)), arma::sum(pos%nMatrix.row(3))};
	return res;
}

arma::mat mes::Calc::getDNDKsi()
{
	arma::mat res;
	res.insert_rows(0, (-0.25*(1 - ksiCol.t())));
	res.insert_rows(1, (-0.25*(1 + ksiCol.t())));
	res.insert_rows(2, (0.25*(1 + ksiCol.t())));
	res.insert_rows(3, (0.25*(1 - ksiCol.t())));
	return res;
}

arma::mat mes::Calc::getDNDEta()
{
	arma::mat res;
	res.insert_rows(0, (-0.25*(1 - etaCol.t())));
	res.insert_rows(1, (0.25*(1 - etaCol.t())));
	res.insert_rows(2, (0.25*(1 + etaCol.t())));
	res.insert_rows(3, (-0.25*(1 + etaCol.t())));
	return res;
}

arma::mat mes::Calc::getJacobian(const arma::dvec &x, const arma::dvec &y)
{
	arma::mat res(4, 4), deta = getDNDEta(), dksi = getDNDKsi();
	res.row(0) = arma::drowvec{arma::sum(deta.col(0) % x), arma::sum(deta.col(1) % x),
		arma::sum(deta.col(2) % x), arma::sum(deta.col(3) % x)};
	res.row(1) = arma::drowvec{arma::sum(deta.col(0) % y), arma::sum(deta.col(1) % y),
		arma::sum(deta.col(2) % y), arma::sum(deta.col(3) % y)};
	res.row(2) = arma::drowvec{arma::sum(dksi.col(0) % x), arma::sum(dksi.col(1) % x),
		arma::sum(dksi.col(2) % x), arma::sum(dksi.col(3) % x)};
	res.row(3) = arma::drowvec{arma::sum(dksi.col(0) % y), arma::sum(dksi.col(1) % y),
		arma::sum(dksi.col(2) % y), arma::sum(dksi.col(3) % y)};
	return res;
}

arma::mat mes::Calc::getDetJ(const arma::mat &Jac)
{
	arma::dvec detJ(4);
	for (int i = 0; i < 4; i++)
	{
		detJ.row(i) = (Jac.col(i)(0)*Jac.col(i)(3) - Jac.col(i)(1)*Jac.col(i)(2));
	}
	return detJ;
}

arma::mat mes::Calc::getJacInv(const arma::mat &Jac, const arma::mat &detJ)
{
	arma::mat JacInv;
	for (int i = 0; i < 4; i++)
		JacInv.insert_cols(i, arma::dvec{Jac.col(i)(3) / detJ(i),-Jac.col(i)(1) / detJ(i),
			-Jac.col(i)(2) / detJ(i),Jac.col(i)(0) / detJ(i)});
	return JacInv;
}

arma::mat mes::Calc::getDndy (const arma::mat &JacInv, const arma::mat &dKsi, const arma::mat &dEta)
{
	arma::mat dndy;
	for (int i = 0; i < 4; i++)
		dndy.insert_rows(i, JacInv.row(0)(i)*dKsi.col(i).t() + JacInv.row(1)(i)*dEta.col(i).t());
	return dndy;
}

arma::mat mes::Calc::getDndx(const arma::mat &JacInv, const arma::mat &dKsi, const arma::mat &dEta)
{
	arma::mat dndx;
	for (int i = 0; i < 4; i++)
		dndx.insert_rows(i, JacInv.row(2)(i)*dKsi.col(i).t() + JacInv.row(3)(i)*dEta.col(i).t());
	return dndx;
}

arma::mat mes::Calc::getLocalHMatrix(double k, const arma::dvec &x, const arma::dvec &y, bool debug)
{
	arma::mat Jac = getJacobian(x, y);
	if (debug)
		Jac.print("Jac:");
	arma::dvec detJ(4);
	detJ = getDetJ(Jac);
	if (debug)
		detJ.print("detJ:");
	arma::mat JacInv;
	JacInv = getJacInv(Jac, detJ);
	if(debug)
		JacInv.print("JacInv:");
	arma::mat dKsi = getDNDKsi();
	arma::mat dEta = getDNDEta();
	if (debug)
	{
		dKsi.print("dKsi:");
		dEta.print("dEta:");
	}
	arma::mat dndx, dndy;
	dndx = getDndx(JacInv, dKsi, dEta);
	dndy = getDndy(JacInv, dKsi, dEta);
	if (debug)
	{
		dndx.print("dndx:");
		dndy.print("dndy:");
	}
	arma::mat dxdxT[4], dydyT[4];
	arma::mat matK[4];
	for (int i = 0; i < 4; i++)
	{
		dxdxT[i] = (dndx.row(i).t()*dndx.row(i))*detJ(i);
		dydyT[i] = (dndy.row(i).t()*dndy.row(i))*detJ(i);
		matK[i] = k * (dxdxT[i] + dydyT[i]);
	}
	if (debug)
	{
		dxdxT[0].print("dxdxT:");
		dydyT[0].print("dydyT:");
		matK[0].print("matK:");
	}
	arma::mat matH;
	for (int i = 0; i < 4; i++)
		matH.insert_cols(i, matK[0].col(i) + matK[1].col(i) + matK[2].col(i) + matK[3].col(i));
	return matH;
}

arma::mat mes::Calc::getLocalHMatrix(Grid & GRID, unsigned int index, bool debug)
{
	Element *e = GRID.getElement(index);
	arma::dvec x = { e->getNodes().at(0).x , e->getNodes().at(1).x , e->getNodes().at(2).x , e->getNodes().at(3).x };
	arma::dvec y = { e->getNodes().at(0).y , e->getNodes().at(1).y , e->getNodes().at(2).y , e->getNodes().at(3).y };
	return getLocalHMatrix(e->getK(), x, y, debug);
}

arma::mat mes::Calc::getLocalCMatrix(const arma::dvec &x, const arma::dvec &y, double c, double ro, bool debug)
{
	arma::mat Jac, detJ, JacInv;
	Jac = getJacobian(x, y);
	if (debug)
		Jac.print("Jac:");
	detJ = getDetJ(Jac);
	if (debug)
		detJ.print("deJ:");
	JacInv = getJacInv(Jac, detJ);
	if (debug)
		JacInv.print("JacInv:");
	arma::mat Ni;
	for (int i = 0; i < 4; i++)
		Ni.insert_rows(i, arma::drowvec{0.25*(1.0 - ksiCol[i])*(1.0 - etaCol[i]), 0.25*(1.0 + ksiCol[i])*(1.0 - etaCol[i]),
			0.25*(1.0 + ksiCol[i])*(1.0 + etaCol[i]), 0.25*(1.0 - ksiCol[i])*(1.0 + etaCol[i])});
	if (debug)
		Ni.print("Ni:");
	arma::mat temp[4];
	for (int i = 0; i < 4; i++)
	{
		temp[i] = Ni.col(i) * Ni.col(i).t()*detJ[i]*c*ro;
		if (debug)
			temp[i].print("temp[" + std::to_string(i) + "]:");
	}
	arma::mat C = temp[0] + temp[1] + temp[2] + temp[3];
	return C;
}

arma::mat mes::Calc::getLocalCMatrix(Grid & GRID, unsigned int index, bool debug)
{
	Element *e = GRID.getElement(index);
	arma::dvec y = { e->getNodes().at(0).x , e->getNodes().at(1).x , e->getNodes().at(2).x , e->getNodes().at(3).x };
	arma::dvec x = { e->getNodes().at(0).y , e->getNodes().at(1).y , e->getNodes().at(2).y , e->getNodes().at(3).y };
	return getLocalCMatrix(x, y, GRID.getC(), GRID.getRo(), debug);
}

arma::mat mes::Calc::getGlobalHMatrix()
{
	return arma::mat();
}

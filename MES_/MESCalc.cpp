#include <armadillo>
#include <string>
#include "MESCalc.h"

const double mes::Calc::ksi = 1.0 / sqrt(3);
const double mes::Calc::eta = 1.0 / sqrt(3);
const arma::dvec mes::Calc::ksiCol = { -1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3), -1.0 / sqrt(3) };
const arma::dvec mes::Calc::etaCol = { -1.0 / sqrt(3), -1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3) };

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
	arma::dvec res = { arma::sum(pos%nMatrix.row(0)), arma::sum(pos%nMatrix.row(1)),
		arma::sum(pos%nMatrix.row(2)), arma::sum(pos%nMatrix.row(3)) };
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
	res.row(0) = arma::drowvec{ arma::sum(deta.col(0) % x), arma::sum(deta.col(1) % x),
		arma::sum(deta.col(2) % x), arma::sum(deta.col(3) % x) };
	res.row(1) = arma::drowvec{ arma::sum(deta.col(0) % y), arma::sum(deta.col(1) % y),
		arma::sum(deta.col(2) % y), arma::sum(deta.col(3) % y) };
	res.row(2) = arma::drowvec{ arma::sum(dksi.col(0) % x), arma::sum(dksi.col(1) % x),
		arma::sum(dksi.col(2) % x), arma::sum(dksi.col(3) % x) };
	res.row(3) = arma::drowvec{ arma::sum(dksi.col(0) % y), arma::sum(dksi.col(1) % y),
		arma::sum(dksi.col(2) % y), arma::sum(dksi.col(3) % y) };
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
		JacInv.insert_cols(i, arma::dvec{ Jac.col(i)(3) / detJ(i),-Jac.col(i)(1) / detJ(i),
			-Jac.col(i)(2) / detJ(i),Jac.col(i)(0) / detJ(i) });
	return JacInv;
}

arma::mat mes::Calc::getDndy(const arma::mat &JacInv, const arma::mat &dKsi, const arma::mat &dEta)
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
	if (debug)
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
	if (debug)
		std::cout << *e << std::endl;
	arma::dvec x = { e->getNodes().at(0)->x , e->getNodes().at(1)->x , e->getNodes().at(2)->x , e->getNodes().at(3)->x };
	arma::dvec y = { e->getNodes().at(0)->y , e->getNodes().at(1)->y , e->getNodes().at(2)->y , e->getNodes().at(3)->y };
	return getLocalHMatrix(e->getConductivity(), x, y, debug);
}

arma::mat mes::Calc::getLocalHMatrix(Element& e, bool debug)
{
	arma::dvec x = { e.getNodes().at(0)->x , e.getNodes().at(1)->x , e.getNodes().at(2)->x , e.getNodes().at(3)->x };
	arma::dvec y = { e.getNodes().at(0)->y , e.getNodes().at(1)->y , e.getNodes().at(2)->y , e.getNodes().at(3)->y };
	return getLocalHMatrix(e.getConductivity(), x, y, debug);

}

arma::mat mes::Calc::getLocalCMatrix(const arma::dvec &x, const arma::dvec &y, double c, double ro, bool debug)
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
	if (debug)
		JacInv.print("JacInv:");
	arma::mat Ni;
	for (int i = 0; i < 4; i++)
		Ni.insert_rows(i, arma::drowvec{ 0.25*(1.0 - ksiCol[i])*(1.0 - etaCol[i]), 0.25*(1.0 + ksiCol[i])*(1.0 - etaCol[i]),
			0.25*(1.0 + ksiCol[i])*(1.0 + etaCol[i]), 0.25*(1.0 - ksiCol[i])*(1.0 + etaCol[i]) });
	if (debug)
		Ni.print("Ni:");
	arma::mat temp[4];
	for (int i = 0; i < 4; i++)
	{
		temp[i] = Ni.col(i) * Ni.col(i).t()*detJ[i] * c*ro;
		if (debug)
			temp[i].print("temp[" + std::to_string(i) + "]:");
	}
	arma::mat C = temp[0] + temp[1] + temp[2] + temp[3];
	return C;
}

arma::mat mes::Calc::getLocalCMatrix(Grid & GRID, unsigned int index, bool debug)
{
	Element *e = GRID.getElement(index);
	arma::dvec x = { e->getNodes().at(0)->x , e->getNodes().at(1)->x , e->getNodes().at(2)->x , e->getNodes().at(3)->x };
	arma::dvec y = { e->getNodes().at(0)->y , e->getNodes().at(1)->y , e->getNodes().at(2)->y , e->getNodes().at(3)->y };
	return getLocalCMatrix(x, y, GRID.getSpecificHeat(), GRID.getDensity(), debug);
}

arma::mat mes::Calc::getGlobalMatrix(Grid& grid, MatrixType type, bool debug)
{
	if (debug) // first and last element
		std::cout << *grid.getElement(0) << std::endl << *grid.getElement(static_cast<unsigned int>(grid.getSize()) - 1) << std::endl;
	arma::mat res(grid.getCols() * grid.getRows(), grid.getCols() * grid.getRows(), arma::fill::zeros);
	if (debug)
		std::cout << sizeof(double) << "*" << res.size() << " = " << sizeof(double)*res.size() << "b = "
		<< sizeof(double)*res.size() / 1024.0 << "kb = " << sizeof(double)*res.size() / 1024.0 / 1024.0 << "mb" << std::endl;

	for (unsigned int k = 0; k < grid.getSize(); k++)
	{
		Element *e = grid.getElement(k);
		arma::mat temp;

		switch (type)
		{
		case mes::H:
			temp = getLocalHMatrix(*e, debug);
			break;
		case mes::C:
			temp = getLocalCMatrix(grid, k, debug);
			break;
		case mes::HBC:
			temp = getHBCMatrix(grid, k, debug);
			break;
		case mes::H_central:
			temp = getLocalHMatrix(*e, debug) + getHBCMatrix(grid, k, debug);
			break;
		default:
			std::cout << "Error determining type of global matrix" << std::endl;
			break;
		}

		if (debug)
		{
			std::string s = "Local matrix ";
			s.append(type ? (type - 1) ? (type - 2) ? "H_central" : "HBC" : "C" : "H").append(" #").append(std::to_string(k)).append(":");
			temp.print(s);
		}

		if (debug && !(k % grid.getSize() / 10))
			std::cout << "inserting " << *e << " into (";
		for (unsigned int i = 0; i < temp.n_cols; i++)
			for (unsigned int j = 0; j < temp.n_rows; j++)
			{
				res(e->getNodes().at(j)->index, e->getNodes().at(i)->index) += temp(i, j);
				if (debug && !(k % 100))
					std::cout << e->getNodes().at(j)->index << ":" << e->getNodes().at(i)->index << ", ";
			}
		if (debug && !(k % 100))
			std::cout << ")" << std::endl;
	}
	return res;
}

arma::mat mes::Calc::getHBCMatrix(Grid& grid, unsigned int index, bool debug)
{
	if (debug)
		std::cout << *grid.getElement(index) << std::endl;
	std::vector<Node*> n = grid.getElement(index)->getNodes();
	arma::mat HBC(4, 4, arma::fill::zeros);
	std::vector<arma::mat> PcMat(4);
	for (int i = 0; i < 4; i++)
		PcMat.at(i) = { 4, 4, arma::fill::zeros };

	arma::mat Ni(2, 4, arma::fill::zeros);

	for (unsigned int i = 0; i < PcMat.size(); i++) // for edge i
	{
		arma::dvec ksi(2, arma::fill::zeros), eta(2, arma::fill::zeros), ones(2, arma::fill::ones);
		if (i == 0) // "if segment" could be remade as a vector of 2x2 matrixes  
		{//-- +- ++ -+
			ksi(0) -= ksi(1) = mes::Calc::ksi;
			eta.fill(-1);
		}
		if (i == 1)
		{
			ksi.fill(1);
			eta(0) -= eta(1) = mes::Calc::eta;
		}
		if (i == 2)
		{
			ksi(0) -= ksi(1) = -mes::Calc::ksi;
			eta.fill(1);
		}
		if (i == 3)
		{
			ksi.fill(-1);
			eta(0) -= eta(1) = -mes::Calc::eta;
		}
		// filling Ni
		Ni.col(0) = (ones - ksi) % (ones - eta);
		Ni.col(1) = (ones + ksi) % (ones - eta);
		Ni.col(2) = (ones + ksi) % (ones + eta);
		Ni.col(3) = (ones - ksi) % (ones + eta);
		Ni = Ni / 4.0;
		if (debug)
		{
			Ni.print("Ni:");
			(Ni.row(0).t() * Ni.row(0)*grid.getAlpha()).print("pc#1:");
			(Ni.row(1).t() * Ni.row(1)*grid.getAlpha()).print("pc#2:");

		}
		// filling PcMat
		PcMat.at(i) += (Ni.row(0).t() * Ni.row(0)*grid.getAlpha());
		PcMat.at(i) += (Ni.row(1).t() * Ni.row(1)*grid.getAlpha());
		// detJ from local length
		double detJ = (sqrt(pow(abs(n.at(i)->x - n.at((i + 1) % 4)->x), 2) +
			pow(abs(n.at(i)->y - n.at((i + 1) % 4)->y), 2))) / 2.0;
		PcMat.at(i) *= detJ;
		if (debug)
		{
			std::cout << detJ << std::endl;
			PcMat.at(i).print("PcMat:");
		}
		// filling HBC
		HBC += PcMat.at(i) * (static_cast<int>(grid.checkEdge(index) / pow(2, i)) % 2);
	}
	return HBC;
}

arma::dvec mes::Calc::getPVector(Grid& grid, unsigned int index, bool debug)
{
	Element *e = grid.getElement(index);
	if (debug)
		std::cout << *e << std::endl;
	std::vector<Node*> n = e->getNodes();
	arma::dvec P(4, arma::fill::zeros);

	arma::dvec Ni(4, arma::fill::zeros);
	arma::mat ksieta({
		{ -ksi, -1, ksi, -1 },
		{ 1, -ksi, 1, ksi },
		{ ksi, 1, -ksi, 1 },
		{ -1, ksi, -1, -ksi }
		});

	double length[4] = { 0, 0, 0, 0 };
	double detJ[4] = { 0, 0, 0, 0 };
	std::vector<Node*> nodess = grid.getElement(index)->getNodes();
	length[0] = sqrt(pow(nodess.at(1)->x - nodess.at(0)->x, 2) + pow(nodess.at(1)->y - nodess.at(0)->y, 2));
	length[1] = sqrt(pow(nodess.at(1)->x - nodess.at(2)->x, 2) + pow(nodess.at(1)->y - nodess.at(2)->y, 2));
	length[2] = sqrt(pow(nodess.at(2)->x - nodess.at(3)->x, 2) + pow(nodess.at(2)->y - nodess.at(3)->y, 2));
	length[3] = sqrt(pow(nodess.at(0)->x - nodess.at(3)->x, 2) + pow(nodess.at(0)->y - nodess.at(3)->y, 2));
	for (int i = 0; i < 4; i++)
		detJ[i] = length[i] / 2;

	for (int i = 0; i < 4; i++)
	{
		// filling Ni
		Ni(0) = ((1 - ksieta.row(i)(0))*(1 - ksieta.row(i)(1)) + (1 - ksieta.row(i)(2))*(1 - ksieta.row(i)(3))) / 4.0;
		Ni(1) = ((1 + ksieta.row(i)(0))*(1 - ksieta.row(i)(1)) + (1 + ksieta.row(i)(2))*(1 - ksieta.row(i)(3))) / 4.0;
		Ni(2) = ((1 + ksieta.row(i)(0))*(1 + ksieta.row(i)(1)) + (1 + ksieta.row(i)(2))*(1 + ksieta.row(i)(3))) / 4.0;
		Ni(3) = ((1 - ksieta.row(i)(0))*(1 + ksieta.row(i)(1)) + (1 - ksieta.row(i)(2))*(1 + ksieta.row(i)(3))) / 4.0;
		if (debug)
			Ni.print("Ni:");

		for (unsigned int j = 0; j < n.size(); j++)
			P(i) -= (static_cast<int>(grid.checkEdge(index) / pow(2, i)) % 2)*Ni(j)*grid.getAmbientTemperature()*grid.getAlpha()*detJ[j];
		if (debug)
		{
			std::string s = "Partial P # ";
			s.append(std::to_string(i));
			P.print(s);
		}
	}
	if (debug)
		P.print(std::string("Local P#").append(std::to_string(index)).append(":"));
	return P;
}

arma::dvec mes::Calc::getGlobalPVector(Grid& grid, bool debug)
{
	arma::dvec res(grid.getNodes().size(), arma::fill::zeros);
	std::vector<Element>& e = grid.getElements();
	for (unsigned int i = 0; i < grid.getSize(); i++)
	{
		arma::dvec temp = getPVector(grid, i); // localP
		std::vector<Node*> n = grid.getElement(i)->getNodes();
		for (unsigned int k = 0; k < 4; k++)
			res(n.at(k)->index) -= temp(k);
	}
	return res + ((getGlobalMatrix(grid, mes::C) / grid.getDeltaTau())*getTemperaturesVector(grid, debug));
}

arma::mat mes::Calc::getHCdTMatrix(Grid& grid, bool debug)
{
	return getGlobalMatrix(grid, mes::H_central, debug) + getGlobalMatrix(grid, mes::C, debug) / grid.getDeltaTau();
}

arma::dvec mes::Calc::getTemperaturesVector(Grid& grid, bool debug)
{
	arma::dvec res(grid.getNodes().size(), arma::fill::zeros);
	for (unsigned int i = 0; i < grid.getNodes().size(); i++)
		res(i) = grid.getNode(i)->t;
	return res;
}

arma::dvec mes::Calc::gauss(Grid& grid, bool debug)
{
	unsigned int n = grid.getNodes().size();
	arma::dvec res(n, arma::fill::zeros);
	arma::mat HCdT = getHCdTMatrix(grid);
	arma::dvec GP = getGlobalPVector(grid);

	arma::mat temp = HCdT;
	temp.reshape(n, n + 1);
	for (unsigned int i = 0; i < n; i++)
		temp(i, n) = GP(i);
	// [temp] = [HCdT][GP]

	if (debug)
	{
		temp.print("[HCdT][GP]:");
	}
	double m = 0, s = 0;

	for (unsigned int i = 0; i < n - 1; i++)
		for (unsigned int j = i + 1; j < n; j++)
			if (temp(i, i))
			{
				m = -temp(j, i) / temp(i, i);
				for (unsigned int k = 0; k < n + 1; k++)
					temp(j, k) += m * temp(i, k);
			}

	for (int i = n - 1; i >= 0; i--)
	{
		s = temp(i, n);
		for (int j = n - 1; j >= 0; j--)
			s -= temp(i, j) * res(j);
		if (temp(i, i))
			res(i) = s / temp(i, i);
	}
	return res;
}

void mes::Calc::applyGauss(Grid& grid, unsigned int iterations, bool debug)
{
	for (; iterations > 0; iterations--)
	{
		unsigned int n = grid.getNodes().size();
		arma::dvec res(n, arma::fill::zeros);
		arma::mat HCdT = getHCdTMatrix(grid);
		arma::dvec GP = getGlobalPVector(grid);

		arma::mat temp = HCdT;
		temp.reshape(n, n + 1);
		for (unsigned int i = 0; i < n; i++)
			temp(i, n) = GP(i);
		// [temp] = [HCdT][GP]

		if (debug)
			temp.print("[HCdT][GP]:");
		double m = 0, s = 0;

		for (unsigned int i = 0; i < n - 1; i++)
			for (unsigned int j = i + 1; j < n; j++)
				if (temp(i, i))
				{
					m = -temp(j, i) / temp(i, i);
					for (unsigned int k = 0; k < n + 1; k++)
						temp(j, k) += m * temp(i, k);
				}

		for (int i = n - 1; i >= 0; i--)
		{
			s = temp(i, n);
			for (int j = n - 1; j >= 0; j--)
				s -= temp(i, j) * res(j);
			if (temp(i, i))
				res(i) = s / temp(i, i);
		}

		for (unsigned int i = 0; i < n; i++)
			grid.getNode(i)->t = res(i);
	}
}

void mes::Calc::printExtremeTemp(Grid& grid, unsigned int iteration)
{
	if (!iteration)
		std::cout << "T_min = " << grid.getMinTemp() << ", T_max = " << grid.getMaxTemp() << std::endl;
	else
	{
		applyGauss(grid, iteration);
		std::cout << "T_min = " << grid.getMinTemp() << ", T_max = " << grid.getMaxTemp() << std::endl;
	}
}
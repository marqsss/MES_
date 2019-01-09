#ifndef MESGRID_H
#define MESGRID_H

#include <iostream>
#include <string>
#include <vector>
#include "MESNode.h"
#include "MESElement.h"
#include <armadillo>

namespace mes
{
	class Grid
	{
	public:
		Grid(std::string filename);
		void print(bool verbose = false);
		double getSpecificHeat() { return specific_heat; }
		void setSpecificHeat(double C) { specific_heat = C; }
		double getDensity() { return density; }
		unsigned int getSize() { return elements.size(); }
		unsigned int getCols() { return n_cols; }
		unsigned int getRows() { return n_rows; }
		Element* getElement(unsigned int index) { return &(elements.at(index)); }

	private:
		std::vector<Node> nodes;
		double gridHeight; // unit/px
		double gridWidth; // unit/px
		unsigned int n_cols;
		unsigned int n_rows;
		std::vector<Element> elements;
		double specific_heat;
		double density;
	};
}

#endif
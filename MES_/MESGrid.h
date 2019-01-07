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
		double getC() { return c; }
		void setC(double C) { c = C; }
		double getRo() { return ro; }
		Element* getElement(unsigned int index) { return &(elements.at(index)); }

	private:
		std::vector<Node> nodes;
		double height; // unit/px
		double width; // unit/px
		unsigned int cols;
		unsigned int rows;
		std::vector<Element> elements;
		double c;
		double ro;
	};
}

#endif
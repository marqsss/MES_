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
	enum Edge {
		Bottom = 1,
		Right=2,
		Top=4,
		Left=8
	};

	class Grid
	{
	public:
		Grid() {}
		// parameters in the text file must be separated by new lines and must be in this order:
		// specific heat (c), density(ro), width of grid in units(pixels), height of grid in units(pixels),
		// number of Node columns, number of Node rows, initial temperature of all Nodes, conductivity,
		// time of simulation(s), time step(s), ambient temperature, alpha
		Grid(std::string filename) { loadFromFile(filename); }
		void reset();
		void loadFromFile(std::string filename);
		void print(bool verbose = false);
		double getSpecificHeat() { return specific_heat; }
		void setSpecificHeat(double C) { specific_heat = C; }
		double getDensity() { return density; }
		size_t getSize() { return elements.size(); }
		unsigned int getCols() { return n_cols; }
		unsigned int getRows() { return n_rows; }
		Element* getElement(unsigned int index) { return &(elements.at(index)); }
		std::vector<Element>& getElements() { return elements; }
		std::vector<Node>& getNodes() { return nodes; }
		Node* getNode(unsigned int index) { return &nodes.at(index); }
		unsigned int checkEdge(unsigned int index);
		double getAmbientTemperature() { return ambientTemperature; }
		double getAlpha() { return alpha; }
		double getDeltaTau() { return timeStep; }
		double getTime() { return totalTime; }
		double getConductivity() { return getElement(0)->getConductivity(); }
		std::vector<double> getExtremeTemp();
		double getMinTemp();
		double getMaxTemp();

	private:
		std::vector<Node> nodes;
		double gridHeight; // unit/px
		double gridWidth; // unit/px
		unsigned int n_cols;
		unsigned int n_rows;
		std::vector<Element> elements;
		double specific_heat;
		double density;
		double totalTime;
		double timeStep;
		double ambientTemperature;
		double alpha;
	};
}

#endif
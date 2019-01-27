#include <iostream>
#include <string>
#include <fstream>
#include "MESGrid.h"


void mes::Grid::loadFromFile(std::string filename)
{
	std::fstream file;
	file.open(filename, std::ios::in);
	if (file.is_open())
	{
		std::string line = "";
		double temperature = 0, conduct = 0;
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#'); // read lines omitting empty ones and #comments
		specific_heat = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		density = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		gridWidth = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		gridHeight = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		n_cols = std::stoul(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		n_rows = std::stoul(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		temperature = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		conduct = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		totalTime = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		timeStep = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		ambientTemperature = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		alpha = std::stod(line);

		// create nodes and elements from the data
		for (unsigned int col = 0; col < n_cols; col++)
			for (unsigned int row = 0; row < n_rows; row++)
				nodes.emplace_back(col*(gridWidth / (n_cols - 1.0)), row*(gridHeight / (n_rows - 1.0)), col*n_rows + row, temperature);
		for (unsigned int col = 0; col < n_cols - 1; col++)
			for (unsigned int row = 0; row < n_rows - 1; row++)
				elements.emplace_back(&nodes.at((col*n_rows) + row), &nodes.at(((col + 1)*n_rows) + row),
					&nodes.at(((col + 1) *n_rows) + 1 + row), &nodes.at((col*n_rows) + row + 1),
					col*(n_rows - 1) + row, conduct);
	}
	else
		printf("Could not open source file\n");
	file.close();
}

void mes::Grid::reset()
{
	nodes.clear();
	elements.clear();
}

void mes::Grid::print(bool verbose)
{
	if (!verbose)
		std::cout << "Grid of (" << gridWidth << "px:" << gridHeight << "px) size with c=" << specific_heat << ";\n"
		<< nodes.size() << " (" << n_cols << "*" << n_rows << " = " << n_cols * n_rows << ") nodes;\n"
		<< elements.size() << " (" << (n_cols - 1) << "*" << (n_rows - 1) << " = " << (n_cols - 1) * (n_rows - 1) << ") elements." << std::endl;
	else
	{
		std::cout << "Grid of (" << gridWidth << ":" << gridHeight << ") size. C=" << specific_heat << ", ro=" << density << ". Nodes:" << std::endl;
		for (unsigned int i = 0; i < n_cols; i++)
		{
			for (unsigned int j = 0; j < n_rows; j++)
				std::cout << nodes.at(i*n_rows + j);
			std::cout << std::endl;
		}
		std::cout << "Elements:" << std::endl;
		for (unsigned int i = 0; i < n_cols - 1; i++)
		{
			for (unsigned int j = 0; j < n_rows - 1; j++)
				std::cout << elements.at(i*(n_rows - 1) + j);
			std::cout << std::endl;
		}
		std::cout << "Element size is " << gridWidth / (n_cols - 1) << " x " << gridHeight / (n_rows - 1) << " units." << std::endl;
	}
}

unsigned int mes::Grid::checkEdge(unsigned int i)
{
	unsigned int res = 0;
	if (!(i % (getRows() - 1)))
		res += Edge::Bottom;
	if (i < getRows() - 1)
		res += Edge::Left;
	if (i > (getCols() - 2)*(getRows() - 1) - 1)
		res += Edge::Right;
	if (!((i + 1) % (getRows() - 1)))
		res += Edge::Top;
	return res;
}

std::vector<double> mes::Grid::getExtremeTemp()
{
	std::vector<double> res;
	res.resize(2);

	double min = nodes.at(0).t, max = nodes.at(0).t;
	for (unsigned int j = 1; j < nodes.size(); j++)
	{
		if (min > nodes.at(j).t)
			min = nodes.at(j).t;
		if (max < nodes.at(j).t)
			max = nodes.at(j).t;
	}
	res.at(0) = min;
	res.at(1) = max;
	return res;
}

double mes::Grid::getMinTemp()
{
	double res = nodes.at(0).t;
	for (unsigned int j = 1; j < nodes.size(); j++)
		if (res > nodes.at(j).t)
			res = nodes.at(j).t;
	return res;
}

double mes::Grid::getMaxTemp()
{
	double res = nodes.at(0).t;
	for (unsigned int j = 1; j < nodes.size(); j++)
		if (res < nodes.at(j).t)
			res = nodes.at(j).t;
	return res;
}
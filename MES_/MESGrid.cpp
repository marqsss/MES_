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
		double temperature = 0, conduct=0;
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
		n_cols = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		n_rows = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		temperature = std::stod(line);
		do
			std::getline(file, line);
		while (line.empty() || line.at(0) == '#');
		conduct = std::stod(line);

		// create nodes and elements from the data
		for (unsigned int col = 0; col < n_cols; col++)
			for (unsigned int row = 0; row < n_rows; row++)
				nodes.emplace_back(col*(gridWidth / (n_cols - 1.0)), row*(gridHeight / (n_rows - 1.0)), col*n_rows + row, temperature);
		for (unsigned int col = 0; col < n_cols - 1; col++)
			for (unsigned int row = 0; row < n_rows - 1; row++)
				elements.emplace_back(&nodes.at((col*n_rows) + row), &nodes.at(((col + 1)*n_rows) + row + 1),
					&nodes.at((col *n_rows) + 1 + row), &nodes.at(((col + 1)*n_rows) + row), col*(n_rows - 1) + row, conduct);
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

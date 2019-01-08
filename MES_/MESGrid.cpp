#include <iostream>
#include <string>
#include <fstream>
#include "MESGrid.h"

mes::Grid::Grid(std::string filename)
{
	std::fstream file;
	file.open(filename, std::ios::in);
	if (file.is_open())
	{
		std::string line = "";
		std::getline(file, line);
		if (!line.empty())
			c = std::stod(line);
		std::getline(file, line);
		if (!line.empty())
			ro = std::stod(line);
		std::getline(file, line);
		if (!line.empty())
			width = std::stoul(line);
		std::getline(file, line);
		if (!line.empty())
			height = std::stoul(line);
		std::getline(file, line);
		if (!line.empty())
			cols = std::stoul(line);
		std::getline(file, line);
		if (!line.empty())
			rows = std::stoul(line);
		for (unsigned int col = 0; col < cols; col++)
			for (unsigned int row = 0; row < rows; row++)
				nodes.emplace_back(col*(width / (cols - 1.0)), row*(height / (rows - 1.0)), col*rows+row);
		for (unsigned int col = 0; col < cols - 1; col++)
			for (unsigned int row = 0; row < rows - 1; row++)
				elements.emplace_back(nodes.at((col*rows) + row), nodes.at(((col + 1)*rows) + row + 1),
					nodes.at((col *rows) + 1 + row), nodes.at(((col + 1)*rows) + row), col*(rows-1)+row);
	}
	else
		printf("Could not open source file\n");
	file.close();
}

void mes::Grid::print(bool verbose)
{
	if (!verbose)
		std::cout << "Grid of (" << width << "px:" << height << "px) size with c="<<c<<";\n"
		<< nodes.size() << " (" << cols << "*" << rows << " = " << cols*rows << ") nodes;\n"
		<< elements.size() << " (" << (cols - 1)<<"*"<<(rows - 1) <<" = "<< (cols - 1) * (rows - 1) << ") elements." << std::endl;
	else
	{
		std::cout << "Grid of (" << width << ":" << height << ") size. K="<<c<<". Nodes:" << std::endl;
		for (unsigned int i = 0; i < cols; i++)
		{
			for (unsigned int j = 0; j < rows; j++)
				std::cout << nodes.at(i*rows + j);
			std::cout << std::endl;
		}
		std::cout << "Elements:" << std::endl;
		for (unsigned int i = 0; i < cols-1; i++)
		{
			for (unsigned int j = 0; j < rows-1; j++)
				std::cout << elements.at(i*(rows-1) + j);
			std::cout << std::endl;
		}
		std::cout << "Element size is " << width / (cols - 1) << " x " << height / (rows - 1) << " units." << std::endl;
	}
}

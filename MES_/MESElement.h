#ifndef MESELEMENT_H
#define MESELEMENT_H

#include <vector>
#include <fstream>

namespace mes
{
	class Element
	{
	public:
		Element(std::vector<Node*> Nodes, unsigned int idx, double K = 30);
		Element(Node *n1, Node *n2, Node *n3, Node *n4, unsigned int idx, double K = 30);
		double getConductivity() { return conductivity; }
		std::vector<Node*> getNodes() { return nodes; }
		unsigned int getIndex() { return index; }
		friend std::ostream& operator<<(std::ostream& s, Element e)
		{
			s << "e#" << e.index << " \t(" << *e.nodes.at(0) << ", " << *e.nodes.at(1) << ", "
				<< *e.nodes.at(2) << ", " << *e.nodes.at(3) << " @ " << e.conductivity << ")";
			return s;
		}
	private:
		std::vector<Node*> nodes;
		double conductivity;
		unsigned int index;
	};
}

#endif
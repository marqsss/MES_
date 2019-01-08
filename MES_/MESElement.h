#ifndef MESELEMENT_H
#define MESELEMENT_H

#include <vector>
#include <fstream>

namespace mes
{
	class Element
	{
	public:
		Element();
		Element(std::vector<Node> Nodes, unsigned int idx, double K = 30);
		Element(Node n1, Node n2, Node n3, Node n4, unsigned int idx, double K = 30);
		double getK() { return k; }
		std::vector<Node> getNodes() { return nodes; }
		friend std::ostream& operator<<(std::ostream& s, Element e)
		{
			s << "# " << e.index << "(" << e.nodes.at(0).index << ", " << e.nodes.at(1).index << ", "
				<< e.nodes.at(2).index << ", " << e.nodes.at(3).index << " @ " << e.k << ")";
			return s;
		}
	private:
		std::vector<Node> nodes;
		double k;
		unsigned int index;
	};
}

#endif
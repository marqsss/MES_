#include <vector>
#include "MESNode.h"
#include "MESElement.h"

mes::Element::Element(std::vector<Node> Nodes, double K) : k(K)
{
	nodes = Nodes;
}

mes::Element::Element(Node n1, Node n2, Node n3, Node n4, double K) : k(K)
{
	nodes.emplace_back(n1);
	nodes.emplace_back(n3);
	nodes.emplace_back(n2);
	nodes.emplace_back(n4);
}
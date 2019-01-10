#include <iostream>
#include <vector>
#include "MESNode.h"
#include "MESElement.h"

mes::Element::Element(std::vector<Node*> Nodes, unsigned int idx, double K) : conductivity(K)
{
	nodes = Nodes;
	index = idx;
}

mes::Element::Element(Node *n1, Node *n2, Node *n3, Node *n4, unsigned int idx, double K) : conductivity(K)
{
	index = idx;
	nodes.emplace_back(n1);
	nodes.emplace_back(n3);
	nodes.emplace_back(n2);
	nodes.emplace_back(n4);
}

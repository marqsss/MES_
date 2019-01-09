#ifndef MESNODE_H
#define MESNODE_H

#include <fstream>

namespace mes
{
	class Node
	{
	public:
		Node(double X, double Y, unsigned int idx, double T = 0) : x(X), y(Y), index(idx), t(T){};
		friend std::ostream& operator<<(std::ostream& s, const Node n)
		{
			s << "#" << n.index << "(" << n.x << ":" << n.y << "@" << n.t << ")";
			return s;
		}
		unsigned int index;
		double x;
		double y;
		double t;
	private:		
	};
}

#endif
#ifndef MESNODE_H
#define MESNODE_H

#include <fstream>

namespace mes
{
	class Node
	{
	public:
		Node();
		Node(double X, double Y, unsigned int NR, double T = 0) : x(X), y(Y), nr(NR), t(T){};
		unsigned int nr;
		friend std::ostream& operator<<(std::ostream& s, const Node n)
		{
			s << "#" << n.nr << "(" << n.x << ":" << n.y << "@" << n.t << ")";
			return s;
		}
		double x;
		double y;
		double t;
	private:		
	};
}

#endif
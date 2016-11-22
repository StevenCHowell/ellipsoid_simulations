// Hailiang Zhang

#include <coor.h>

// stdio
std::ostream & operator<<(std::ostream & os, const Coor & c)
{
	os<<c._x<<" "<<c._y<<" "<<c._z;
    return os;
}

// destructor
Coor::
~Coor()
{
}

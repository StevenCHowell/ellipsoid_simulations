// Hailiang Zhang

#ifndef COOR_H
#define COOR_H

// externels
#include <cmath>
#include <iostream>
#include <vector>

class Coor
{
	// data
	private:
		float _x;
		float _y;
		float _z;

	// interface
	public:
		inline float distance(const Coor &) const; // distance
		inline Coor operator+(const Coor &) const; // add operator
		inline Coor operator-(const Coor &) const; // minus operator
		inline float operator*(const Coor &) const; // multiply operator
		inline Coor operator*(const float) const; // multiply by a value
		inline Coor & operator+=(const Coor &); // add/assign operator
		inline Coor & operator-=(const Coor &); // minus/assign operator
		inline Coor & operator*=(const float); // *= operator
		inline Coor & operator/=(const float); // /= operator
		inline float x() const; // get the x coordinate
		inline float y() const; // get the y coordinate
		inline float z() const; // get the z coordinate

	// stdio
	public:
		friend std::ostream & operator<<(std::ostream &, const Coor &);

	// constructor/destructor
	public:
		inline Coor();
		inline Coor(const float x, const float y, const float z);
		inline Coor(const std::vector<float> & coor);
		virtual ~Coor();
};


#define COOR_ICC
#include "coor.icc"
#undef COOR_ICC 

#endif

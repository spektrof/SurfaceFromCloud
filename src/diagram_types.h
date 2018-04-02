#pragma once

#include "cgal_types.h"

struct Pole
{
	Point* center;
	float radius;
	Pole(Point* c = NULL, const float& r = 0.0f) : center(c), radius(r) {}

	bool isNull() const { return center == NULL; }
};
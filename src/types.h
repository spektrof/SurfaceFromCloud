#pragma once

#include <iostream>
#include <boost/sort/spreadsort/spreadsort.hpp>		//TODO just for rightshift - so probably wont need this

// Coords & Textures later-------------- -
struct Coord
{
	float x, y, z;
	Coord(const float& _x = 0.0f, const float& _y = 0.0f, const float& _z = 0.0f) : x(_x), y(_y), z(_z) {}

	friend class XYZOrder;
	friend class YZXOrder;
	friend class ZXYOrder;

	class XYZOrder {
	public:
		bool operator()(const Coord& c, const Coord& _c) {
			return c.x < _c.x || (c.x == _c.x &&  c.y < _c.y) || (c.x == _c.x && c.y == _c.y && c.z < _c.z);
		}
	};
	class YZXOrder {
	public:
		bool operator()(const Coord& c, const Coord& _c) {
			return c.y < _c.y || (c.y == _c.y &&  c.z < _c.z) || (c.y == _c.y && c.z == _c.z && c.x < _c.x);
		}
	};
	class ZXYOrder {
	public:
		bool operator()(const Coord& c, const Coord& _c) {
			return c.z < _c.z || (c.z == _c.z &&  c.x < _c.x) || (c.z == _c.z && c.x == _c.x && c.y < _c.y);
		}
	};

	friend std::ostream& operator << (std::ostream& out, const Coord& c)
	{
		out << c.x << " " << c.y << " " << c.z << "\n";
		return out;
	}
};

inline bool operator < (const Coord& c, const Coord& _c)
{
	return c.x < _c.x || (c.x == _c.x &&  c.y < _c.y) || (c.x == _c.x && c.y == _c.y && c.z < _c.z);
}

inline bool operator <= (const Coord& c, const Coord& _c)
{
	return c.x <= _c.x || (c.x == _c.x &&  c.y <= _c.y) || (c.x == _c.x && c.y == _c.y && c.z <= _c.z);
}

inline bool operator >= (const Coord& c, const Coord& _c)
{
	return c.x >= _c.x || (c.x == _c.x &&  c.y >= _c.y) || (c.x == _c.x && c.y == _c.y && c.z >= _c.z);
}

inline bool operator > (const Coord& c, const Coord& _c)
{
	return c.x > _c.x || (c.x == _c.x &&  c.y > _c.y) || (c.x == _c.x && c.y == _c.y && c.z > _c.z);
}

//---------------------------------------
//Boxes for nodes - just for testing kdTree
struct Section
{
	float min, max;
	Section(const float& mi = 0.0f, const float& ma = 20.0f) : min(mi), max(ma) {}

	bool isBetween(const float& val) { return val >= min && val <= max; }
	bool isLower(const float& val) { return val < min; }
	bool isHigher(const float& val) { return val > max; }
};

struct Box
{
	Section x, y, z;

	Box(const Section& s) : x(s), y(s), z(s) {}
	Box(const float& min = 0.0f, const float& max = 20.0f) : x(Section(min, max)), y(Section(min, max)), z(Section(min, max)) {}
	Box(const Section& _x, const Section& _y, const Section& _z) : x(_x), y(_y), z(_z) {}
	Box(const float& xmi, const float& xma, const float& ymi, const float& yma, const float& zmi, const float& zma) : x(Section(xmi, xma)), y(Section(ymi, yma)), z(Section(zmi, zma)) {}

	bool isInside(const Coord& c) { return x.isBetween(c.x) && y.isBetween(c.y) && z.isBetween(c.z); }

	float getXMin() { return x.min; }
	float getXMax() { return x.max; }
	float getYMin() { return y.min; }
	float getYMax() { return y.max; }
	float getZMin() { return z.min; }
	float getZMax() { return z.max; }

	float getXLength() { return x.max - x.min; }
	float getYLength() { return y.max - y.min; }
	float getZLength() { return z.max - z.min; }
};

//-----------------------------------------------

//rossz de alapból rossz speciális compare funcra
struct rightshiftCoord {
	int operator()(const Coord &c, const unsigned offset) const {
		return boost::sort::spreadsort::float_mem_cast<float, int>(c.x) >> offset;
	}
};

#pragma once

#include <iostream>

enum ProcessReturns
{
	CLOUD_LOAD_SUCCESS,
	PROCESS_DONE,
	NO_INPUT_DATA,
	PROCESS_STOPPED,
	EMPTY_POINTS_SET,
	NOT_ENABLED_COMPONENT,
	UNDEFINED_ERROR
};

enum DrawPossibilites
{
	POINTSS,
	BOX,
	SURFACE,
	UNDEFINED
};

static float sumOfsquare(const float& x, const float& y, const float& z)
{
	return x * x + y * y + z * z;
}

/*
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

inline Coord& operator - (const Coord& lhs, const Coord& rhs)
{
	return Coord(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

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
*/

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

	//bool isInside(const Point& c) { return x.isBetween(c.x()) && y.isBetween(c.y()) && z.isBetween(c.z()); }

	float getXMin() const { return x.min; }
	void setXMin(const float& v) {		x.min = v; }
	float getXMax() const { return x.max; }
	void setXMax(const float& v) {		x.max = v; ; }
	float getYMin() const { return y.min; }
	void setYMin(const float& v) {		y.min = v; ; }
	float getYMax() const { return y.max; }
	void setYMax(const float& v) {		y.max = v; ; }
	float getZMin() const { return z.min; }
	void setZMin(const float& v) {		z.min = v; ; }
	float getZMax() const { return z.max; }
	void setZMax(const float& v) {		z.max = v; ; }

	float getXLength() const { return x.max - x.min; }
	float getYLength() const { return y.max - y.min; }
	float getZLength() const { return z.max - z.min; }

	float getVolume() const { return getXLength() * getYLength() * getZLength(); }
};

//-----------------------------------------------
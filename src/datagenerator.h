#pragma once
#include "types.h"
#include "defines.h"
#include <vector>

class DataGenerator
{
public:
	DataGenerator() : box(Box(-1.0f,1.0f)), n_point(0) {}

	~DataGenerator() {}

	void setNumberOfPoint(const unsigned int& pn) { n_point = pn; }
	unsigned int getNumberOfPoint() const { return n_point; }

	void setBox(const Section& s) { box = Box(s); }
	void setBox(const float& min, const float& max) { box = Box(min, max); }
	void setBox(const Section& _x, const Section& _y, const Section& _z) { box = Box(_x, _y, _z); }
	void setBox(const float& xmi, const float& xma, const float& ymi, const float& yma, const float& zmi, const float& zma) { box = Box(xmi, xma, ymi, yma, zmi, zma); }

	//to get starter box in cloudcontainer but we need to calculate it when we read from file... TODO
	Box* getBoxObject() { return &box;  }

	//---------------------------

	std::vector<Coord> getCloudRandom();
	std::vector<Coord>  getCloudSpecific();
	//std::set<Coord> getCloudFromFile();
	//std::set<Coord> getCloudFromClient();

private:
	Box box;
	unsigned int n_point;

	//TODO: udp/tcp server
};
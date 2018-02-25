#include "datagenerator.h"
#include <cstdlib>
#include <ctime>

#include <set>		//for random point generator

std::vector<Coord> DataGenerator::getCloudRandom()
{
	srand(static_cast <unsigned> (time(0)));

	std::set<Coord> res;

	float xmin = box.getXMin();
	float xmax = box.getXMax();
	float ymin = box.getYMin();
	float ymax = box.getYMax();
	float zmin = box.getZMin();
	float zmax = box.getZMax();

	while (res.size() != n_point)
	{
		Coord tmp;
		tmp.x = xmin + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (xmax - xmin)));
		tmp.y = ymin + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (ymax - ymin)));
		tmp.z = zmin + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (zmax - zmin)));

		res.insert(tmp);
	}

	return std::vector<Coord>(res.begin(), res.end());
}

std::vector<Coord> DataGenerator::getCloudSpecific()
{
	std::vector<Coord> tmp;

	tmp.push_back(Coord(2, 3, 3));
	tmp.push_back(Coord(5, 4, 2));
	tmp.push_back(Coord(9, 6, 7));
	tmp.push_back(Coord(4, 7, 9));
	tmp.push_back(Coord(8, 1, 5));
	tmp.push_back(Coord(7, 2, 6));
	tmp.push_back(Coord(9, 4, 1));
	tmp.push_back(Coord(8, 4, 2));
	tmp.push_back(Coord(9, 7, 8));
	tmp.push_back(Coord(6, 3, 1));
	tmp.push_back(Coord(3, 4, 5));
	tmp.push_back(Coord(1, 6, 8));
	tmp.push_back(Coord(9, 5, 3));
	tmp.push_back(Coord(2, 1, 3));
	tmp.push_back(Coord(8, 7, 6));

	return tmp;
}
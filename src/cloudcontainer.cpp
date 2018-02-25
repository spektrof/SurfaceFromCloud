#include "cloudcontainer.h"

CloudContainer::CloudContainer(const int& n_poi, const Section& s) :
	pC(RANDOM), dP(POINTS),
	tmp(0)
{
	box = dataGen.getBoxObject();
	dataGen.setNumberOfPoint(n_poi);
	dataGen.setBox(s);
}

std::vector<Coord> CloudContainer::getDrawingData() const
{
	switch (dP)
	{
	case POINTS:
		return data;
	case BOX:
		return drawNode();
	case SURFACE:

	default:
		return std::vector<Coord>();
	}
}

CloudContainer::DATAFLOWRTN CloudContainer::getDataFromClient()
{
	if (tmp > 0) return UNDEFINED_ERROR;

	if (pC == RANDOM)
	{
		releaseNodes();
		data = dataGen.getCloudRandom();
		return DATA_LOAD_SUCCESS;
	}
	else
	{
		releaseNodes();
		data = dataGen.getCloudSpecific();
		return DATA_LOAD_SUCCESS;
	}

	return UNDEFINED_ERROR;
}
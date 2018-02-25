#pragma once

#include "types.h"
#include "datagenerator.h"

#include <vector>

class CloudContainer
{
public:
	CloudContainer(const int& n_poi = 20000, const Section& s = Section(-20.0f, 20.0f));

	virtual ~CloudContainer() {}

	enum DRAWPOSSIBILITIES
	{
		POINTS,
		BOX,
		SURFACE,
		UNDEFINED
	};

	//Virtuals
	virtual void buildAllNode() = 0;
	virtual void releaseNodes() = 0;
	virtual bool process() = 0;

	virtual std::vector<Coord> drawNode() const = 0;
	virtual void goLeftNode() = 0;
	virtual void goRightNode() = 0;
	virtual void goParentNode() = 0;
	
	virtual void changeThreadCapacity(const unsigned int&) = 0;
	virtual unsigned int getThreadCapacity() const = 0;
	virtual unsigned int getPointNumbFromBox() const = 0;
	virtual void getTimes(float& join, float& sort, float& all) const = 0;

	//Setters&Getters
	void setNumberOfPoints(const int& n_p) { dataGen.setNumberOfPoint(n_p); }
	void setDataRange(const float& min, const float& max) { dataGen.setBox(min, max); }
	void setDataRange(const Section& s) { dataGen.setBox(s); }

	//TODO: wait until it finish the actual session
	void changeDrawType(const unsigned int& dp) { dP = DRAWPOSSIBILITIES(dp); }

	std::vector<Coord> getDrawingData() const;
	DRAWPOSSIBILITIES getDrawType() const { return dP; }

	void resetCloud() { tmp = 0; }
	void changeComputation() { tmp = 1 - tmp; }

protected:

	enum POINTCLOUD
	{
		RANDOM,
		SPECIFIC,
		FILE,
		UDP
	};

	enum DATAFLOWRTN
	{
		DATA_LOAD_SUCCESS,
		FAILED_ON_OPEN_FILE,
		FAILED_ON_TIMEOUT,
		UNDEFINED_ERROR
	};

	std::vector<Coord> data;
	Box* box;

	DATAFLOWRTN getDataFromClient();

private:
	DataGenerator dataGen;

	POINTCLOUD pC;
	DRAWPOSSIBILITIES dP;

	int tmp;	//TODO change logic name or delete
};
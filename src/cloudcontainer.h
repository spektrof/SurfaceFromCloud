#pragma once

#include "types.h"
#include "c_io.h"

#include <vector>
#include <queue>
#include <limits>

template<typename Data>
class CloudContainer
{
	enum CloudSource;

protected:
	struct node;
	typedef std::pair<float, node*> DistanceTuple;

	struct LargestOnTop 
	{
		bool operator()(const DistanceTuple &a, const DistanceTuple &b) const { return a.first < b.first; }
	};
	typedef std::priority_queue<DistanceTuple, std::vector<DistanceTuple>, LargestOnTop > kNearest_queue;

	struct nearest_node
	{
		node* neighb;
		float distance;
		nearest_node(node* n = NULL, const float& d = FLT_MAX /*std::numeric_limits<float>::max()*/ ) : neighb(n), distance(d) {}
	};

	struct node
	{
		node* left, *right;
		node* parent;

		int depth, length;
		Box box;			//TODO: just pointer and make an other node for boxes -> compute_bounding_boxes

		Data* point;	//pointer to node val
		Data* start;	//pointer to starter point in the vector

		nearest_node nn;
		kNearest_queue knn;
		float spacing;		//TODO: some functions for compute this, i.e. centroid , compute_average_spacing()

		node(node* p, Data* s, const int& d, const int& lg, Box b, Data* po = NULL, node* l = NULL, node* r = NULL, const nearest_node& _nn = nearest_node())
			: parent(p), start(s), depth(d), length(lg), box(b), point(po), left(l), right(r), nn(_nn) {}
	};

	//------------------------------------
	std::vector<Data> data;
	Box* box;

	int meshaccuracy;

public:
	//maybe we dont need the pcloud if we separate the consturctors
	CloudContainer(const unsigned int& pCloud, const int& ma, bool repeat, const int& n_poi, const Section& s = Section(-20.0f, 20.0f))	:
		meshaccuracy(ma),
		cloud_source(CloudSource(pCloud)),
		run_count(0), started(false), repeatable(repeat)
	{
		box = c_io.getBoxObject();
		c_io.setNumberOfPoint(n_poi);
		c_io.setBox(s);
	}

	CloudContainer(const unsigned int& pCloud, const int& ma, bool repeat, char* filename) :
		meshaccuracy(ma),
		cloud_source(CloudSource(pCloud)),
		run_count(0), started(false), repeatable(repeat)
	{
		c_io.setFileName(filename);
	}

	virtual ~CloudContainer() {}

	//Virtuals
	virtual void build() = 0;
	virtual void release() = 0;
	virtual ProcessReturns load_process() = 0;
	virtual ProcessReturns build_process(const std::vector<Data>&) = 0;

	virtual void goLeftNode() = 0;
	virtual void goRightNode() = 0;
	virtual void goParentNode() = 0;
	
	virtual void changeThreadCapacity(const unsigned int&) = 0;
	virtual unsigned int getThreadCapacity() const = 0;
	virtual unsigned int getPointNumbFromBox() const = 0;
	virtual void getTimes(float& join, float& sort, float& all) const = 0;
	virtual void getActivePoints(std::vector<Data>& active_points, node* r = NULL) = 0;

	virtual node* getRootPtr() const = 0;
	virtual GLPaintFormat drawNode() = 0;

	//Setters&Getters
	void setNumberOfPoints(const int& n_p) { c_io.setNumberOfPoint(n_p); }
	void setDataRange(const float& min, const float& max) { c_io.setBox(min, max); }
	void setDataRange(const Section& s) { c_io.setBox(s); }
	void setMeshAccuracy(const int& ma) { if (data.size() > ma ) meshaccuracy = ma; }

	//TODO Refactor this
	void resetCloud() { run_count = 0; }
	void swapStartFlag() 
	{ 
		started = 1 -  started; 
		if (cloud_source == FILE && !c_io.checkFileType()) c_io.resetFileIndex();
	}

	bool checkFileType() { return c_io.checkFileType(); }

	std::vector<Data> getData() const { return data;  }

	/*void nearest() { 
		bt = BOXDRAWTYPES(1 - bt); 
	}
	void knearest() 
	{ 
		if ( bt % 2 != 0 )	bt = BOXDRAWTYPES(2);
		else bt = BOXDRAWTYPES(2 - bt);
	}*/

protected:

	enum CloudSource
	{
		RANDOM,
		FILE,
		SPECIFIC,
		UDP
	};

	/*enum BOXDRAWTYPES
	{
		POINT,
		NEAREST,
		KNEAREST,
		ALL
	};

	BOXDRAWTYPES bt;*/

	ProcessReturns getDataFromClient()
	{
		if (!started) return PROCESS_STOPPED;

		switch (cloud_source)
		{
			case RANDOM:
			{
				if (run_count > 0 && !repeatable) return NO_INPUT_DATA;
				run_count = 1;

				release();
				data = c_io.getCloudRandom();
				return CLOUD_LOAD_SUCCESS;
			}
			case FILE:
			{
				if (!c_io.nextFileExist())
				{
					if (!repeatable) return NO_INPUT_DATA;
					c_io.resetFileIndex();
				}

				release();
				data = c_io.getCloudFromFile();
				box = c_io.getBoxObject();
				return CLOUD_LOAD_SUCCESS;
			}
			case SPECIFIC:
			{
				release();
				data = c_io.getCloudFromFile();
				return CLOUD_LOAD_SUCCESS;
			}
			default:
				return UNDEFINED_ERROR;
		}
			
	}

private:
	CloudSource cloud_source;
	C_IO<Data> c_io;

	int run_count;	
	bool started;
	bool repeatable;
};

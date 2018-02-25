#pragma once

#include "types.h"
#include "cloudcontainer.h"
#include "taskscheduler.h"

#include <set>	// only neighbours right now
#include <cmath> //for log2
#include <boost/chrono.hpp>

//TODO: do not test with boost sort, its wrong
#include <boost/sort/sort.hpp>
#include <boost/sort/spreadsort/float_sort.hpp>

class kdTree : public CloudContainer
{
	struct Node
	{
		Node* left, *right;
		Node* parent;

		int depth, length;
		Box box;

		Coord* point;	//pointer to node val
		Coord* start;	//pointer to starter point in the vector

		Node(Node* p, Coord* s, const int& d, const int& lg, Box b, Coord* po = NULL, Node* l = NULL, Node* r = NULL)
			: parent(p), start(s), depth(d), length(lg), box(b), point(po), left(l), right(r) {}
	};

public:
	kdTree(const int& thr_c = 4, const int& n_poi = 20000, const Section& s = Section(-20.0f, 20.0f));

	~kdTree() {
		releaseNodes();
		delete t_schedular;
	}

	//Overrides
	void buildAllNode() override;
	void releaseNodes() override;
	bool process() override;

	void changeThreadCapacity(const unsigned int& tc) override;

	void getTimes(float& join, float& sort, float& all) const override;
	unsigned int getThreadCapacity() const override { return t_schedular->getCloudThreadCapacity(); }
	unsigned int getPointNumbFromBox() const override;

	std::vector<Coord> drawNode() const override;
	void goLeftNode() override    { if (_draw && _draw->left)	_draw = _draw->left; }
	void goRightNode() override   { if (_draw && _draw->right)	_draw = _draw->right; }
	void goParentNode() override  { if (_draw && _draw->parent) _draw = _draw->parent; }

	//--------------
	void inOrderPrint(Node* r, const int& level);

protected:
	std::set<Node*> getNeighbour(Node* n) { std::set<Node*> tmp; return tmp; }

	void doBuildWithAutomech();
	void doBuildRecursive(Node* r, bool stop = false, const int& level = 0);

	void sortPartVector(const int& axis, Coord* c, const int& length);
	void releaseNode(Node*);

private:
	Node * root;
	Node* _draw;

	TaskSchedular<kdTree, Node>* t_schedular;	//TODO: in cloudcontainer???
												//REFACTOR mybe
	boost::chrono::duration<float> time_sortMS;
	boost::chrono::milliseconds calc_time;

	unsigned int stopBuildDepth;

};

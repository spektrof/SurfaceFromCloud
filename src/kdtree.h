#pragma once

#include "types.h"
#include "cloudcontainer.h"
#include "taskscheduler.h"

#include <set>
#include <cmath> //for log2
#include <boost/chrono.hpp>
#include <QDebug>

template <typename Data>
class kdTree : public CloudContainer<Data>
{
	//task struct, task* megy a taskmanagerbe - hátrány nekünk kell megszüntetni, elõny: eltünik a start és a length

	typedef boost::chrono::high_resolution_clock::time_point time_p;
	typedef boost::chrono::milliseconds milliseconds;
	typedef boost::chrono::duration<float> duration_f;
	typedef boost::chrono::high_resolution_clock hrclock;
	typedef boost::chrono::high_resolution_clock::time_point time_p;

public:
	typedef typename node* kdnode_ptr;

	kdTree(const unsigned int& pCloudType, const int& ma, bool repeat, const int& thr_c = 4, const int& n_poi = 20000, const Section& s = Section(-0.1f, 0.1f))
		: CloudContainer<Data>(pCloudType, ma, repeat, n_poi, s),
		root(NULL), _draw(NULL),
		stopBuildDepth(STOPLEVEL),
		kn(1),
		sort_time(duration_f(0.0f)), calc_time(milliseconds(0))
	{
		t_schedular = new TaskSchedular<kdTree, node>();
		changeThreadCapacity(thr_c);
	}

	kdTree(const unsigned int& pCloudType, const int& ma, bool repeat, const int& thr_c, char* filename) 
		: CloudContainer<Data>(pCloudType, ma, repeat, filename),
		root(NULL), _draw(NULL),
		stopBuildDepth(STOPLEVEL),
		kn(1),
		sort_time(duration_f(0.0f)), calc_time(milliseconds(0))
	{
		t_schedular = new TaskSchedular<kdTree, node>();
		t_schedular->setCloudThreadCapacity(thr_c);
	}

	~kdTree() {
		release();
		delete t_schedular;
	}

	//Overrides
	void build() override
	{
		calc_time = milliseconds(0);
		sort_time = duration_f(0.0f);

		root = new node(NULL, &data[0], 0, data.size(), *box, NULL, NULL, NULL);	//TODO CHECK DATA
		_draw = root;

		//TODO: test the performance
		if (data.size() > CLOUDRECURSIVEPOINTCAP && t_schedular->getCloudThreadCapacity() >= 2)
			buildWithAutomech();
		else
			buildRecursive(root);		//with few points the simple recursion is faster than syncing with main thread

		if (PRINT) inOrderPrint(root, 0);
	}

	void release() override
	{
		_draw = NULL;
		if (root) releaseNode(root);

		data.clear();
		//TODO: IF WE ADD DIFFERENT POINT -> new variable with one more clear here
//		active_data.clear();
	}

	ProcessReturns load_process() override
	{
		ProcessReturns res;
		if (res = getDataFromClient()) return res;

		return res;
	}

	//REFACTOR
	//EZ MÁS DE TODO nem tudom nem megcsinálni, mert nem akarok adatot másolni
	//nem tudom kiszedni több szálon csak lockfree vel, de annak van korlátja
	//viszont meg kell csinálni, ha nem akarok minden pontot kirajzolni
	// - -> 1x végig kell mennem rekurzívan a kapott fán, ha 1 volt az accuracy ha nem

	ProcessReturns build_process(const std::vector<Data>& input) override
	{
		data = input;

		time_p t1 = hrclock::now();
		time_p t2;

		build();
		if (!root) return UNDEFINED_ERROR;		//TODO

		//----------------------------------------
		t = 0;
		//addNeighbours();
		addKNeighbours();

		//----------------------------------------
		t2 = hrclock::now();
		calc_time = boost::chrono::duration_cast<milliseconds>(t2 - t1);

		return PROCESS_DONE;
	}

	//------------------------------
	void changeThreadCapacity(const unsigned int& tc) override;

	void getTimes(float& join, float& sort, float& all) const override;
	unsigned int getThreadCapacity() const override { return t_schedular->getCloudThreadCapacity(); }
	unsigned int getPointNumbFromBox() const override;

	void goLeftNode() override { if (_draw && _draw->left)	_draw = _draw->left; }
	void goRightNode() override { if (_draw && _draw->right)	_draw = _draw->right; }
	void goParentNode() override { if (_draw && _draw->parent) _draw = _draw->parent; }

	node* getRootPtr() const override { return root; }
	GLPaintFormat drawNode(/*const unsigned int&*/) override;

	//--------------------------------------------------------------
	void nearestAll(node* query, bool stop = false, const int& stopDepth = 0)
	{
		//test = 0;
		nearestRecursive(query);
		//t++;

		if (stop && stopDepth == query->depth + 1)
		{
			if (query->left)   t_schedular->addTask(query->left);
			if (query->right)  t_schedular->addTask(query->right);
			return;
		}
		else
		{
			if (query->left) nearestAll(query->left);
			if (query->right) nearestAll(query->right);
		}
	}

	void knearestAll(node* query, bool stop = false, const int& stopDepth = 0)
	{
		//test = 0;
		kNearestRecursive(query);
		//t++;

		if (stop && stopDepth == query->depth + 1)
		{
			if (query->left)   t_schedular->addTask(query->left);
			if (query->right)  t_schedular->addTask(query->right);
			return;
		}
		else
		{
			if (query->left) knearestAll(query->left);
			if (query->right) knearestAll(query->right);
		}
	}

	//------------------------------
	void inOrderPrint(node* r, const int& level);

protected:
	//Tree buildings
	void buildWithAutomech();
	void buildRecursive(node* r, bool stop = false, const int& level = 0);
	void sortPartVector(const int& axis, Data* c, const int& length);

	//----------------
	// Neighbours calculators
	void addNeighbours()
	{
		if (data.size() > CLOUDRECURSIVEPOINTCAP && t_schedular->getCloudThreadCapacity() >= 2)
		{
			if (root->depth != stopBuildDepth) nearestAll(root, true, stopBuildDepth);
			else t_schedular->addTask(root);

			while (!t_schedular->isEmptyTask())
				t_schedular->addSubscribeShit(this, &kdTree::nearestAll);

			t_schedular->joinAll();		//wait our threads
		}
		else
			nearestAll(root);
	}

	void addKNeighbours()
	{
		if (data.size() > CLOUDRECURSIVEPOINTCAP && t_schedular->getCloudThreadCapacity() >= 2)
		{
			if (root->depth != stopBuildDepth) knearestAll(root, true, stopBuildDepth);
			else t_schedular->addTask(root);

			while (!t_schedular->isEmptyTask())
				t_schedular->addSubscribeShit(this, &kdTree::knearestAll);

			t_schedular->joinAll();		//wait our threads
		}
		else
			knearestAll(root);
	}

	void nearestRecursive(node* query) {
		if (!root) {
			return;
		}
		nearest(query->point, root, query->nn);
	}

	void kNearestRecursive(node* query) {
		if (!root) {
			return;
		}
		knearest(query->point, root, query->knn);
	}

	void nearest(Data *query, node* currentNode, nearest_node& best);

	void knearest(Data *query, node* currentNode, kNearest_queue& best);

	//----------------
	//Others
	void releaseNode(node*);

	void getActivePoints(std::vector<Data>& active_points, node* r = NULL);

private:
	//node* _nearest;
	node * root;
	node* _draw;

	unsigned int stopBuildDepth;
	TaskSchedular<kdTree, node>* t_schedular;	//TODO: in cloudcontainer???
												//REFACTOR mybe
	int kn;										//TODO: better name?? k nearest neighbours

	duration_f sort_time;
	milliseconds calc_time;
//----------
	int test;
	boost::atomic_int t;
};

//--------------------------
//---		-----		----
//million point under 0,4 sec with 4 thread
template <typename Data>
void kdTree<Data>::buildWithAutomech()
{
	if (root->depth != stopBuildDepth) buildRecursive(root, true, stopBuildDepth);
	else t_schedular->addTask(root);

	while (!t_schedular->isEmptyTask())
		t_schedular->addSubscribeWithoutAutomechanism(this, &kdTree::buildRecursive);

	t_schedular->joinAll();		//wait our threads
}

template <typename Data>
void kdTree<Data>::buildRecursive(node* r, bool stop, const int& stopDepth)
{
	if (r->length == 1 || r->length <= meshaccuracy || r->box.getVolume() < MINVOLUME)
	{
		r->point = r->start + (r->length >> 1);
		return;
	}

	int axis = r->depth % DIMENSION;
	int med = r->length >> 1;
	//---------------------
	time_p t1 = hrclock::now();

	sortPartVector(axis, r->start, r->length);

	time_p t2 = hrclock::now();
	sort_time += (boost::chrono::duration_cast<duration_f>(t2 - t1));
	//---------------------

	r->point = r->start + med;

	//---------------------
	Box Lbox = r->box;
	Box Rbox = r->box;

	if (!axis)
	{
		Lbox.setXMax(r->point->x());
		Rbox.setXMin(r->point->x());
	}
	else if (axis == 1)
	{
		Lbox.setYMax(r->point->y());
		Rbox.setYMin(r->point->y());
	}
	else
	{
		Lbox.setZMax(r->point->z());
		Rbox.setZMin(r->point->z());
	}

	r->left = r->point - r->start == 0 ? NULL : new node(r, r->start, r->depth + 1, r->point - r->start, Lbox);
	r->right = r->length - r->left->length - 1 == 0 ? NULL : new node(r, r->point + 1, r->depth + 1, r->length - r->left->length - 1, Rbox);

	r->length = 1;

	if (stop && stopDepth == r->depth + 1)
	{
		if (r->left)  t_schedular->addTask(r->left);
		if (r->right) t_schedular->addTask(r->right);
		return;
	}
	else
	{
		if (r->left)  buildRecursive(r->left, stop, stopDepth);
		if (r->right) buildRecursive(r->right, stop, stopDepth);
	}
}

template <typename Data>
void kdTree<Data>::releaseNode(node* r)
{
	node* _l = r->left;
	node* _r = r->right;

	delete r;

	if (_l) releaseNode(_l);
	if (_r) releaseNode(_r);
}

template <typename Data>
void kdTree<Data>::sortPartVector(const int& axis, Data* start, const int& length)
{
	if (axis  == 2)
	std::nth_element(start, start + length /2, start + length, CGALTool::ZXYOrder());
	else if (axis == 1)
	std::nth_element(start, start + length / 2, start + length, CGALTool::YZXOrder());
	else
	std::nth_element(start, start + length / 2, start + length, CGALTool::XYZOrder());
}

//---		-----		----
//--------------------------
//---		-----		----
// Neighbours

template <typename Data>
void kdTree<Data>::nearest(Data *query, node* currentNode, nearest_node& best) {
	if (!currentNode) return;
	//test++;

	int axis = currentNode->depth % DIMENSION;

	float d = sumOfsquare(currentNode->point->x - query->x, currentNode->point->y - query->y, currentNode->point->z - query->z);

	if (d == 0.0f)	//megtaláltuk a pontunkat, mindkét gyerek irányba menni kell mert a pont rajta van a vágósíkon
	{
		nearest(query, currentNode->left, best);
		nearest(query, currentNode->right, best);
		return;
	}

	//adott vágósíktól vett távolság
	float dx = axis == 0 ? currentNode->point->x - query->x : axis == 1 ? currentNode->point->y - query->y : currentNode->point->z - query->z;

	if (d < best.distance) {
		best.neighb = currentNode;
		best.distance = d;
	}

	node* _near = dx <= 0 ? currentNode->right : currentNode->left;
	node* _far = dx <= 0 ? currentNode->left : currentNode->right;
	nearest(query, _near, best);
	if (dx * dx >= best.distance) return;	//pitagoras
	nearest(query, _far, best);
}

template <typename Data>
void kdTree<Data>::knearest(Data *query, node* currentNode, kNearest_queue& best) {
	if (!currentNode) return;

	int axis = currentNode->depth % DIMENSION;

	float d = sumOfsquare(currentNode->point->x() - query->x(), currentNode->point->y() - query->y(), currentNode->point->z() - query->z());

	if (d == 0.0f)	//megtaláltuk a pontunkat, mindkét gyerek irányba menni kell mert a pont rajta van a vágósíkon
	{
		knearest(query, currentNode->left, best);
		knearest(query, currentNode->right, best);
		return;
	}

	//adott vágósíktól vett távolság
	float dx = axis == 0 ? currentNode->point->x() - query->x() : axis == 1 ? currentNode->point->y() - query->y() : currentNode->point->z() - query->z();

	if (best.size() < kn || d <= best.top().first) {
		best.push(DistanceTuple(d, currentNode));
		/*kNearest_queue tmp = best;
		while (!tmp.empty())
		{
		qDebug() << tmp.top() << " ";
		tmp.pop();
		}
		qDebug() << "\n";*/
		if (best.size() > kn) {
			best.pop();
			/*tmp = best;
			while (!tmp.empty())
			{
			qDebug() << tmp.top() << " ";
			tmp.pop();
			}
			qDebug() << "\n";*/
		}
	}

	node* _near = dx <= 0 ? currentNode->right : currentNode->left;
	node* _far = dx <= 0 ? currentNode->left : currentNode->right;
	knearest(query, _near, best);
	if (dx * dx >= best.top().first) return;	//pitagoras
	knearest(query, _far, best);

}

//---		-----		----
//--------------------------

template <typename Data>
GLPaintFormat kdTree<Data>::drawNode(/*const unsigned int& bt*/)
{
	GLPaintFormat gpf;
	//node* draw = bt == 0 ? _draw : bt == 1 ? _draw->nn.neighb : NULL;
	if (!_draw) return gpf;

	Box act_box = _draw->box;

	gpf.box_points.push_back(Data(act_box.getXMin(), act_box.getYMin(), act_box.getZMin()));
	gpf.box_points.push_back(Data(act_box.getXMax(), act_box.getYMin(), act_box.getZMin()));
	gpf.box_points.push_back(Data(act_box.getXMin(), act_box.getYMax(), act_box.getZMin()));
	gpf.box_points.push_back(Data(act_box.getXMax(), act_box.getYMax(), act_box.getZMin()));

	gpf.box_points.push_back(Data(act_box.getXMin(), act_box.getYMin(), act_box.getZMax()));
	gpf.box_points.push_back(Data(act_box.getXMax(), act_box.getYMin(), act_box.getZMax()));
	gpf.box_points.push_back(Data(act_box.getXMin(), act_box.getYMax(), act_box.getZMax()));
	gpf.box_points.push_back(Data(act_box.getXMax(), act_box.getYMax(), act_box.getZMax()));

	
	gpf.ix.push_back(0);	gpf.ix.push_back(1);	gpf.ix.push_back(2);
	gpf.ix.push_back(2);	gpf.ix.push_back(1);	gpf.ix.push_back(3);
							
	gpf.ix.push_back(4);	gpf.ix.push_back(0);	gpf.ix.push_back(6);
	gpf.ix.push_back(6);	gpf.ix.push_back(0);	gpf.ix.push_back(2);

	gpf.ix.push_back(4);	gpf.ix.push_back(5);	gpf.ix.push_back(0);
	gpf.ix.push_back(0);	gpf.ix.push_back(5);	gpf.ix.push_back(1);
	
	gpf.ix.push_back(1);	gpf.ix.push_back(5);	gpf.ix.push_back(3);
	gpf.ix.push_back(3);	gpf.ix.push_back(5);	gpf.ix.push_back(7);

	gpf.ix.push_back(5);	gpf.ix.push_back(4);	gpf.ix.push_back(7);
	gpf.ix.push_back(7);	gpf.ix.push_back(4);	gpf.ix.push_back(6);

	gpf.ix.push_back(2);	gpf.ix.push_back(3);	gpf.ix.push_back(6);
	gpf.ix.push_back(6);	gpf.ix.push_back(3);	gpf.ix.push_back(7);
	
	//gpf.ix.push_back(1);	gpf.ix.push_back(5);	gpf.ix.push_back(7);
	
	return gpf;
}

//---		-----		----
//--------------------------

template <typename Data>
void kdTree<Data>::getActivePoints(std::vector<Data>& active_points, node* r = NULL)
{
	if (!r) r = root;
	if (!r) return;

	active_points.emplace_back(*(r->point));

	if (r->left) getActivePoints(active_points, r->left);
	if (r->right) getActivePoints(active_points, r->right);
}

template <typename Data>
void kdTree<Data>::changeThreadCapacity(const unsigned int& tc)
{
	t_schedular->setCloudThreadCapacity(tc);
	stopBuildDepth = (unsigned int)log2(tc);
}

template <typename Data>
void kdTree<Data>::getTimes(float& join, float& sort, float& all) const
{
	//	std::cout << "Sorting time: " << time_sortMS << "\nJoin time: " << t_schedular->getJoinTime() << "\n";
	join = t_schedular->getJoinTime();
	sort = sort_time.count();
	all = calc_time.count();
}

//TODO: CHANGE IF NECESSARY - depends on remove from kdtree
template <typename Data>
unsigned int kdTree<Data>::getPointNumbFromBox() const
{
	return t;
	return test;
	return 0;
}

template <typename Data>
void kdTree<Data>::inOrderPrint(node* r, const int& level)
{
	if (r == NULL) return;
	if (r->point == NULL)
	{
		std::cout << "Empty Node\n";
		return;
	}

	inOrderPrint(r->left, level + 1);
	std::cout << level << ": " << *(r->point);
	inOrderPrint(r->right, level + 1);
}

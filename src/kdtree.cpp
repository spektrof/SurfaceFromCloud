#include "kdtree.h"

#include <algorithm>

kdTree::kdTree(const int& thr_c, const int& n_poi , const Section& s) : CloudContainer(n_poi, s),
	calc_time( boost::chrono::milliseconds(0) ),
	time_sortMS( boost::chrono::duration<float>(0.0f) ),
	root( NULL ),
	_draw( NULL ),
	stopBuildDepth( STOPLEVEL )
{
	t_schedular = new TaskSchedular<kdTree, Node>();
	t_schedular->setCloudThreadCapacity(thr_c);
}

//--------------------------
//---		-----		----

void kdTree::buildAllNode()
{
	calc_time = boost::chrono::milliseconds(0);
	time_sortMS = boost::chrono::duration<float>(0.0f);

	boost::chrono::high_resolution_clock::time_point t1 = boost::chrono::high_resolution_clock::now();;
	boost::chrono::high_resolution_clock::time_point t2;

	root = new Node(NULL, &data[0], 0, data.size(), *box , NULL, NULL, NULL);
	_draw = root;

   //TODO: test the performance
	if (data.size() > CLOUDRECURSIVEPOINTCAP && t_schedular->getCloudThreadCapacity() >= 2)
		doBuildWithAutomech();
	else
		doBuildRecursive(root);		//with few points the simple recursion is faster than syncing with main thread

	t2 = boost::chrono::high_resolution_clock::now();
	if (PRINT) inOrderPrint(root, 0);

	calc_time = boost::chrono::duration_cast<boost::chrono::milliseconds>(t2 - t1);
}

//million point under 0,4 sec with 4 thread
void kdTree::doBuildWithAutomech()
{
	if (root->depth != stopBuildDepth) doBuildRecursive(root, true, stopBuildDepth);
	else t_schedular->addTask(root);

	while (!t_schedular->isEmptyTask())
		t_schedular->addSubscribeWithoutAutomechanism(this, &kdTree::doBuildRecursive);

	t_schedular->joinAll();		//wait our threads
}

void kdTree::doBuildRecursive(kdTree::Node* r, bool stop, const int& stopDepth)
{
	if (r->length == 1)
	{
		r->point = r->start;
		return;
	}

	int axis = r->depth % DIMENSION;
	int med = r->length >> 1;

	//---------------------
	boost::chrono::high_resolution_clock::time_point t1 = boost::chrono::high_resolution_clock::now();

	sortPartVector(axis, r->start, r->length);

	boost::chrono::high_resolution_clock::time_point t2 = boost::chrono::high_resolution_clock::now();
	time_sortMS += (boost::chrono::duration_cast<boost::chrono::duration<float>>(t2 - t1));

	r->point = r->start + med;

	//---------------------
	float xMin = r->box.getXMin();
	float xMax = r->box.getXMax();
	float yMin = r->box.getYMin();
	float yMax = r->box.getYMax();
	float zMin = r->box.getZMin();
	float zMax = r->box.getZMax();

	//do set X coord for Box maybe, TODO - refactor
	Box Lbox = !axis ? Box(xMin, r->point->x, yMin, yMax, zMin, zMax) :
		   axis == 1 ? Box(xMin, xMax, yMin, r->point->y, zMin, zMax) :
		               Box(xMin, xMax, yMin, yMax, zMin, r->point->z);

	Box Rbox = !axis ? Box(r->point->x, xMax, yMin, yMax, zMin, zMax) :
		   axis == 1 ? Box(xMin, xMax, r->point->y, yMax, zMin, zMax) :
		               Box(xMin, xMax, yMin, yMax, r->point->z, zMax);


	r->left = r->point - r->start == 0 ? NULL : new Node(r, r->start, r->depth + 1, r->point - r->start, Lbox);
	r->right = r->length - r->left->length - 1 == 0 ? NULL : new Node(r, r->point + 1, r->depth + 1, r->length - r->left->length - 1, Rbox);

	if (stop && stopDepth == r->depth + 1)
	{
		if (r->left)  t_schedular->addTask(r->left);
		if (r->right) t_schedular->addTask(r->right);
		return;
	}
	else
	{
		if (r->left)  doBuildRecursive(r->left, stop, stopDepth);
		if (r->right) doBuildRecursive(r->right, stop, stopDepth);
	}
}

bool kdTree::process()
{
	if (getDataFromClient() == DATA_LOAD_SUCCESS)
	{
		buildAllNode();
		//getNeighbour
		//getSimplifiedCloud
		return true;
	}

	return false;
}

void kdTree::releaseNodes()
{
	if (root) releaseNode(root);

	data.clear();
}

void kdTree::releaseNode(Node* r)
{
	Node* _l = r->left;
	Node* _r = r->right;

	delete r;

	if (_l) releaseNode(_l);
	if (_r) releaseNode(_r);
}

void kdTree::sortPartVector(const int& axis, Coord* start, const int& length)
{
	//switch extract as only 3 options are available
	if (axis /*& 0x3)*/ == 2)
		BOOST ? boost::sort::spreadsort::float_sort(start, start + length, rightshiftCoord(), Coord::ZXYOrder()) : std::sort(start, start + length, Coord::ZXYOrder());
	else if (axis /*& 0x3)*/ == 1)
		BOOST ? boost::sort::spreadsort::float_sort(start, start + length, rightshiftCoord(), Coord::YZXOrder()) : std::sort(start, start + length, Coord::YZXOrder());
	else
		BOOST ? boost::sort::spreadsort::float_sort(start, start + length, rightshiftCoord(), Coord::XYZOrder()) : std::sort(start, start + length, Coord::XYZOrder());
}

//---		-----		----
//--------------------------

std::vector<Coord> kdTree::drawNode() const
{
	std::vector<Coord> tmp;
	if (!_draw) return tmp;

	Box act_box = _draw->box;

	//TODO INIT LIST??
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMin(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMin(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMax(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMax(), act_box.getZMin()));
	
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMax(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMax(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMax(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMax(), act_box.getZMax()));
	
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMin(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMin(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMax(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMax(), act_box.getZMax()));
	
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMin(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMin(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMax(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMax(), act_box.getZMin()));
	
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMin(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMin(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMin(), act_box.getZMin()));
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMin(), act_box.getZMin()));
		
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMin(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMin(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMax(), act_box.getYMax(), act_box.getZMax()));
	tmp.emplace_back(Coord(act_box.getXMin(), act_box.getYMax(), act_box.getZMax()));
	
	return tmp;
}

void kdTree::changeThreadCapacity(const unsigned int& tc)
{
	t_schedular->setCloudThreadCapacity(tc);
	stopBuildDepth = (unsigned int)log2(tc);
}

void kdTree::getTimes(float& join, float& sort, float& all) const
{
	//	std::cout << "Sorting time: " << time_sortMS << "\nJoin time: " << t_schedular->getJoinTime() << "\n";
	join = t_schedular->getJoinTime();
	sort = time_sortMS.count();
	all = calc_time.count();
}

//TODO: CHANGE IF NECESSARY - depends on remove from kdtree
unsigned int kdTree::getPointNumbFromBox() const
{
	return 0;
}

void kdTree::inOrderPrint(kdTree::Node* r, const int& level)
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

#pragma once
#include "defines.h"
#include "types.h"
#include "cgal_types.h"
#include "filters.h"
#include "normalestimators.h"
#include "kdtree.h"

//#include <boost/tuple/tuple.hpp>
#ifdef ENABLE_CGAL_SURFACE
	#include <CGAL/Poisson_reconstruction_function.h>
	#include <CGAL/Surface_mesh_default_triangulation_3.h>
	#include <CGAL/Implicit_surface_3.h>
	#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#endif

#ifdef ENABLE_DEBUG
	#include <QDebug>
#endif

class SurfaceGenerator
{
	enum Algorithm
	{
		POISSON,
		MINE,
		PARTIALDIFF
	};

	//TODO: after kdtree refactor put these to types.h
	typedef boost::chrono::high_resolution_clock::time_point time_p;
	typedef boost::chrono::milliseconds milliseconds;
	typedef boost::chrono::duration<float> duration_f;
	typedef boost::chrono::high_resolution_clock hrclock;
	typedef boost::chrono::high_resolution_clock::time_point time_p;

	typedef typename kdTree<Point>::kdnode_ptr kdnode_ptr;

#ifdef ENABLE_CGAL_SURFACE
	typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
	typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;
	typedef CGAL::Surface_mesh_default_triangulation_3 STr;
	typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
#endif

public:
	SurfaceGenerator() :
		alg(POISSON), draw_object(POINTSS)
		//,root(NULL)
		, normal_estimator(NULL)
	{
		cloud = new kdTree<Point>(1, 2, false, 3, "D:/ELTE/Diplomamunka/sfc/sfc/clouds/sanyi/00008_127.0.0.1.ply");	//TODO: file path from config

		filters.clear();

		last_simplif = last_outlier = last_plane_outliner = 0;

		all_calculation_time       = milliseconds(0);
		surface_calculation_time   = milliseconds(0);
		cloud_load_time			   = milliseconds(0);

#ifdef ENABLE_CGAL_SURFACE
		normal_estimator = new PcaNormal<PwN_vector>();	//TODO: adding this as well
#endif
	}

	~SurfaceGenerator()
	{
		//std::for_each(filters.begin(), filters.end(), [](Filter<Point>* fil) { delete fil; });	//for simplification / smooth algorithms

		//active_points.clear();
		//active_nodes.clear();
	//	delete cloud;
	}
	
	/*************************************
	*	Process and Paint helper functions	*
	**************************************/

	ProcessReturns process()
	{
		time_p start = hrclock::now();
		time_p tmp;

		// STEP 0. : cloud import
		switch (cloud->load_process())
		{
			case CLOUD_LOAD_SUCCESS:
			{
				points.clear();
				points = cloud->getData();
				break;
			}
			case NO_INPUT_DATA:
				return NO_INPUT_DATA;
			default:
				return UNDEFINED_ERROR;
		}

		tmp = hrclock::now();
		cloud_load_time = boost::chrono::duration_cast<milliseconds>(tmp - start);

		//STEP 1. : remove outliers: a) plane based, b) average space based

		auto filter_part = last_plane_outliner;
		auto it = 0;

		//a)
		for (; it < filter_part; ++it)
		{
			filters[it]->filter_process(points);
		}

		qDebug() << "a : " << points.size() << "\n";
#ifdef ENABLE_CGAL_FILTER
		//b)
		filter_part += last_outlier;
		compute_average_spacing();

		for (; it < filter_part; ++it)
		{
			filters[it]->SetThirdProperty(average_spacing);
			filters[it]->filter_process(points);
			compute_average_spacing();
		}
#endif
		qDebug() << "b : " << points.size() << "\n";

		//STEP 2. : USE MINE KDTREE
		//REFACTOR -> put into filter
		/*switch (cloud->build_process(points))
		{
			case PROCESS_DONE:
			{
				points.clear();
				cloud->getActivePoints(points);
				break;
			}
			case UNDEFINED_ERROR:
			default:
			{

			#ifdef ENABLE_DEBUG
				qDebug() << "UNDEFINED ERROR at MINE KDTREE\n";
			#endif
			}
		}*/

		//STEP 3. : do pointset simplification and smooth
#ifdef ENABLE_CGAL_FILTER
		bool changed = false;

		//simplifiers
		filter_part += last_simplif;
		for (; it < filter_part; ++it)
		{
			filters[it]->filter_process(points);
			changed = true;
		}

		if (changed) compute_average_spacing();

		//smoothers
		for (; it < filters.size(); ++it)
		{
			filters[it]->filter_process(points);
		}

#endif
		qDebug() << "filtering done : " << points.size() << "\n";

		tmp = hrclock::now();
		all_calculation_time = boost::chrono::duration_cast<milliseconds>(tmp - start);

		return PROCESS_DONE;
	}

	ProcessReturns surface_process()
	{
		points_normals.clear();
		if (points.empty()) return EMPTY_POINTS_SET;

		time_p start = hrclock::now();
		time_p tmp;

		switch (alg)
		{
		case POISSON:
		{
		#ifdef ENABLE_CGAL_SURFACE
			for (auto it : points)
				points_normals.push_back(Point_with_Normal(std::make_pair(it, Vector())));

			normal_estimator->normal_c_process(points_normals);		//nullcheck!!

			PoissonByCGAL();

			tmp = hrclock::now();
			surface_calculation_time = boost::chrono::duration_cast<milliseconds>(tmp - start);
			qDebug() << surface_calculation_time.count();
			return PROCESS_DONE;
		#else
			return NOT_ENABLED_COMPONENT;
		#endif
		}
		case MINE:
			return NOT_ENABLED_COMPONENT;
		case PARTIALDIFF:
			return NOT_ENABLED_COMPONENT;
		default:
			return UNDEFINED_ERROR;
		}
	}

	GLPaintFormat getPaintData()
	{
		switch (draw_object)
		{
			case POINTSS:
			{
				return GLPaintFormat(points);
			}
			case BOX:
			{
				return cloud->drawNode();
			}
			case SURFACE:
			{
				return getPaintFromPolyhedron();
			}
			default:
				return GLPaintFormat();
		}
	}


	/***********************
	*	New Objectums      *
	***********************/

	bool newFileInput(char* filename, const int& ma, bool repeat)
	{
		for (auto& it : filters) delete it;
		position_link.clear();
		filters.clear();
		last_outlier = last_plane_outliner = last_simplif = 0;

		points_normals.clear();
		points.clear();
		mesh.clear();

		all_calculation_time = milliseconds(0);
		surface_calculation_time = milliseconds(0);
		cloud_load_time = milliseconds(0);

		const unsigned int core = cloud->getThreadCapacity();
		if (cloud) delete cloud;

		cloud = new kdTree<Point>(1, ma, repeat, core, filename);

		return cloud->checkFileType();
	}

	bool newRandomInput(const unsigned int& n, const int& ma, bool repeat)
	{
		for (auto& it : filters) delete it;
		position_link.clear();
		filters.clear();
		last_outlier = last_plane_outliner = last_simplif = 0;

		points_normals.clear();
		points.clear();
		mesh.clear();

		all_calculation_time = milliseconds(0);
		surface_calculation_time = milliseconds(0);
		cloud_load_time = milliseconds(0);

		const unsigned int core = cloud->getThreadCapacity();
		if (cloud) delete cloud;

		cloud = new kdTree<Point>(0, ma, repeat, core, n);

		return true;
	}
	//-------------------------

	bool addNewOutlierFilter(const unsigned int& type, const float& kn_px, const float& limit_py = 0.0f, const float& pz = 0.0f, const float& nx = 0.0f, const float& ny = 0.0f, const float& nz = 0.0f)
	{
		switch (type)
		{
			case 1:
			{
			#ifdef ENABLE_CGAL_FILTER
				int index = last_plane_outliner + last_outlier;

				filters.insert(filters.begin() + index, new OutlierRemoval<Point>(kn_px, limit_py));
				last_outlier++;

				qDebug() << " Before \n";
				for (auto it : position_link) qDebug() << it << "\n";

				for (auto& it : position_link) { if (it >= index) it++; }
				position_link.push_back(index);
				qDebug() << " After \n";
				for (auto it : position_link) qDebug() << it << "\n";
				
				return true;
			#else
				return false;
			#endif
			}
			case 2:
			{
				int index = last_plane_outliner;

				filters.insert(filters.begin() + index,
					new OutlierComponentRemoval<Point>(Vector(kn_px, limit_py, pz), Vector(nx, ny, nz)));
				last_plane_outliner++;

				qDebug() << " Before \n";
				for (auto it : position_link) qDebug() << it << "\n";

				for (auto& it : position_link) { if (it >= index) it++; }
				position_link.push_back(index);
				qDebug() << " After \n";
				for (auto it : position_link) qDebug() << it << "\n";

				return true;
			}
		}
		return false;
	}

	bool addNewSimplifierFilter(const unsigned int& type, const float& first, const float& second = 0.0f)
	{
		qDebug() << "NEW SIMPLIF " << type << " " << first <<" " << second << "\n";
		switch (type)
		{
		case 1:
		{
		#ifdef ENABLE_CGAL_FILTER
			int index = last_plane_outliner + last_outlier + last_simplif;

			filters.insert(filters.begin() + index, new GridSimplification<Point>(first));
			last_simplif++;

			qDebug() << " Before \n";
			for (auto it : position_link) qDebug() << it << "\n";

			for (auto& it : position_link) { if (it >= index) it++; }
			position_link.push_back(index);
			qDebug() << " After \n";
			for (auto it : position_link) qDebug() << it << "\n";

			return true;
		#else
			return false;
		#endif
		}
		case 2:
		{
			//TODO MINE
			//filters.push_back(new OutlierComponentRemoval<Point>());
			return false;
		}
		case 3:
		{
#ifdef ENABLE_CGAL_FILTER
			int index = last_plane_outliner + last_outlier + last_simplif;

			filters.insert(filters.begin() + index, new HiearchySimplification<Point>(first, second));
			last_simplif++;

			qDebug() << " Before \n";
			for (auto it : position_link) qDebug() << it << "\n";

			for (auto& it : position_link) { if (it >= index) it++; }
			position_link.push_back(index);
			qDebug() << " After \n";
			for (auto it : position_link) qDebug() << it << "\n";
			return true;
#else
			return false;
#endif
		}
	}
/*
#ifdef ENABLE_CGAL_SURFACE
		filters.push_back(new WLOPSimplification<Point>());
#endif*/
		return false;
	}

	bool addNewSmootherFilter(const unsigned int& type, const float& first, const float& second = 0.0f)
	{
		switch (type)
		{
			case 1:
			{
			#ifdef ENABLE_CGAL_FILTER
				#ifdef CGAL_EIGEN3_ENABLED
	
					filters.push_back(new JetSmoothing<Point>(first));

					position_link.push_back(filters.size()-1);
					for (auto it : position_link) qDebug() << it << "\n";

					return true;
				#else
					return false;
				#endif
			#else	
				return false;
			#endif
			}
		}
		return false;
	}

	bool addNewSurfaceReconstructor()
	{
	#ifdef ENABLE_CGAL_SURFACE
		normal_estimator = new PcaNormal<PwN_vector>();

		//TODO
		return true;
	#else
		return false;
	#endif
	}

//-----------------------------------------------------------------------------

	void SetFilterProperty(const unsigned int& id, const double& first, const double& second = 0.0f, const double& third = 0.0f)
	{
		filters[position_link[id]]->SetFirstProperty(first);
		filters[position_link[id]]->SetSecondProperty(second);
		filters[position_link[id]]->SetThirdProperty(third);
	}

	void SetFilterProperty(const unsigned int& id, const double& x, const double& y, const double& z, const double& nx, const double& ny, const double& nz)
	{
		filters[position_link[id]]->SetFirstProperty(Vector(x,y,z));
		filters[position_link[id]]->SetSecondProperty(Vector(nx,ny,nz));
	}

	//TODO: replace to shared_ptr
	void deleteFilter(const unsigned int& id)
	{
		int rm = position_link[id];
		int index = last_plane_outliner;

		delete filters[rm];
		filters.erase(filters.begin() + rm, filters.begin() + rm + 1);

		position_link.erase(position_link.begin() + id, position_link.begin() + id + 1);

		for (auto& it : position_link) { if (it > rm) it--; }
		for (auto it : position_link) qDebug() << it << "\n";

		if (index > rm)
		{
			last_plane_outliner--;
			qDebug() << "Plane outlier remov. deleted\n";
			return;
		}

		index += last_outlier;

		if (index > rm)
		{	
			last_outlier--;
			qDebug() << "Outlier remov. deleted\n";
			return;
		}

		index += last_simplif;

		if (index > rm)
		{	
			last_simplif--;
			qDebug() << "Simplifier deleted\n";
			return;
		}

		qDebug() << "Smoother deleted\n";
		return;
		
	}

//-----------------------------------------------------------------------------
	//Setters & Getters
	//TODO: wait until it finish the actual session
	void setDrawObject(const unsigned int& dp) { draw_object = DrawPossibilites(dp); }
	DrawPossibilites getDrawObject() const { return draw_object; }

	//-------------------------
	// Routers for cloud
	void swapStartFlag() { cloud->swapStartFlag(); }
	void setMeshAccuracy(const int& val) { cloud->setMeshAccuracy(val); }
	void resetCloud() { cloud->resetCloud(); }
	void goRightNode() { cloud->goRightNode(); }
	void goLeftNode() { cloud->goLeftNode(); }
	void goParentNode() { cloud->goParentNode(); }
	int getPointNumbFromBox() const { return cloud->getPointNumbFromBox(); }
	void getTimes(float& c_load_time, float& s_cal_time, float& all_time, float& join_time, float& sort_time) {
		all_time = all_calculation_time.count();
		s_cal_time = surface_calculation_time.count();
		c_load_time = cloud_load_time.count();
		qDebug() << s_cal_time;

		float dummy;
		cloud->getTimes(join_time, sort_time, dummy); 
	}
	void changeThreadCapacity(const unsigned int& cloud_t) { cloud->changeThreadCapacity(cloud_t); }
	unsigned int getThreadCapacity() const { return cloud->getThreadCapacity(); }

	/************************************************************************
	*																		*
	*************************************************************************/
protected:
/*
2. minden cgallal
3. TBB vel linkelni a CGALT

	Voronoi -> convex hull (CGAL) -> polyhedron (CGAL)
*/

/*
	2. Beolvasás + kn szomszéd alapján outlier eltávolítás + CGAL SImplifi. + CGAL Poisson-> összehasonlítás -> MI LEGYEN???
	2. b) TBB vel CGAL
	3. Voronoi
	4. Döntés + Vornoi based surface rec.
	5. Rajzolásuk
	----------------------------
	6. Saját simplifierek, saját architektúrában
*/

	//-----------------------------------------
	/***********************
	*	Filter componenet   *
	***********************/
#ifdef ENABLE_CGAL_FILTER
	void compute_average_spacing();
#endif

	//-----------------------------------------
	/***********************
	*	Surface componenet   *
	***********************/
#ifdef ENABLE_CGAL_SURFACE

	Polyhedron mesh;

	void PoissonByCGAL();

#endif

	//-----------------------------------------

	std::vector<Point> getResult() {}
	std::vector<Point> getVoronoy() {}

	GLPaintFormat getPaintFromPolyhedron();
	
private:
	/***********************
	*	Private componenet   *
	***********************/

	Algorithm alg;
	NormalEstimator<PwN_vector>* normal_estimator;
	DrawPossibilites draw_object;

	CloudContainer<Point>* cloud;			//TODO: delete the points after everything is done, rethink how can we use our kdtree effectively, related to the points variable
	std::vector< Filter<Point>* > filters;	//for simplification / smooth algorithms
	std::vector<int> position_link;			//position link between the UI and this

	unsigned int last_plane_outliner;
	unsigned int last_outlier;
	unsigned int last_simplif;

	// [Outlier, ... , Simplif, Simplif ,... ,Simplif, Smooth, ..., Smooth]
	
	double average_spacing;		//TODO for filtering... mybe dont need here

	std::vector<Point> points;
	PwN_vector points_normals;

	//kdnode_ptr root;
	//std::vector<kdnode_ptr> active_nodes;

	//-----------------------------
	milliseconds all_calculation_time;
	milliseconds surface_calculation_time;
	milliseconds cloud_load_time;
};

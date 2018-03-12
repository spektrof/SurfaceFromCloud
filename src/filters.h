#pragma once
#include "cgal_types.h"		//tag and vector, mybe refactor

#include <QDebug>

#ifdef ENABLE_CGAL_FILTER

#include <CGAL/remove_outliers.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/hierarchy_simplify_point_set.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>

#include <CGAL/jet_smooth_point_set.h>
//#include <CGAL/bilateral_smooth_point_set.h>

#endif

template< typename T>
class Filter
{
public:
	Filter() {}
	virtual ~Filter() = 0;

	virtual void filter_process(std::vector<T>& points) = 0;
	virtual void SetFirstProperty(const double&) {}
	virtual void SetSecondProperty(const double&) {}
	virtual void SetFirstProperty(const Vector&) {}
	virtual void SetSecondProperty(const Vector&) {}
	virtual void SetThirdProperty(const double&) {}

};
	
template< typename T>
Filter<T>::~Filter() {
	qDebug() << "DELETED FILTER\n";
}	//provided body for pure virtual destr

template< typename T>
class Simplifier : public Filter<T>
{
public:
	Simplifier() {}
	virtual ~Simplifier() = 0;
	virtual void filter_process(std::vector<T>& points) = 0;
	virtual void SetFirstProperty(const double&) {}
	virtual void SetSecondProperty(const double&) {}
	virtual void SetThirdProperty(const double&) {}
	virtual void SetFirstProperty(const Vector&) {}
	virtual void SetSecondProperty(const Vector&) {}
	
};
template< typename T>
Simplifier<T>::~Simplifier() {
	qDebug() << "DELETED Simplifier\n";
}	//provided body for pure virtual destr

template< typename T>
class Smoother : public Filter<T>
{
public:
	Smoother() {}
	virtual ~Smoother() = 0;
	virtual void filter_process(std::vector<T>& points) = 0;
	virtual void SetFirstProperty(const double&) {}
	virtual void SetSecondProperty(const double&) {}
};
template< typename T>
Smoother<T>::~Smoother() { }	//provided body for pure virtual destr

//----------------------------------
//----------	------	-----------

/******************************
*	Outlier removal componenet   *
******************************/

template< typename T>
class OutlierComponentRemoval : public Simplifier<T>
{
public:
	OutlierComponentRemoval(const Vector& _p = Vector(0, 0, 0), const Vector& _n = Vector(1, 0, 0)) : p(_p), n(_n) {}
	~OutlierComponentRemoval() {
		qDebug() << "DELETED OUTLIER\n";
	}

	void filter_process(std::vector<T>& points) override
	{
		float md = p.x()*n.x() + p.y()*n.y() + p.z()*n.z();		//calculate d in plane equation
		std::vector<size_t> del_pos;

		for (size_t it = 0; it < points.size(); ++it)
		{
			float s_d = points[it].x()*n.x() + points[it].y()*n.y() + points[it].z()*n.z() - md;
			if (s_d < 0) del_pos.push_back(it);
		}
		//qDebug() << md << " , " << p.x() << "-" << p.y() << "-" << p.z() << " " << n.x() << "-" << n.y() << "-" << n.z() << " " << del_pos.size() << "\n";
		if (del_pos.size() == 0) return;

		for (auto it = del_pos.end() - 1; it >= del_pos.begin(); --it)
		{
			points.erase(points.begin() + *it, points.begin() + *it + 1);
		}

	}

	void SetFirstProperty(const Vector& first) override {
		qDebug() << first.x() << " " << first.y() << " " <<  first.z() << "\n";
		p = first; 
	}
	void SetSecondProperty(const Vector& second)override { n = second; }

private:
	Vector p;
	Vector n;
};

#ifdef ENABLE_CGAL_FILTER
template< typename T>		//kn neighbours, outlier_limit, average spacing
class OutlierRemoval : public Simplifier<T>
{
public:
	OutlierRemoval(const int& kn = 6, const double& o_lim = 2.0f, const double& avg_sp = 1.0f) : nb_neighbors(kn), outlier_limit(o_lim), average_spacing(avg_sp) {}
	~OutlierRemoval() {}

	void filter_process(std::vector<T>& points) override
	{
		qDebug() << "Outlier removal in progress with " << nb_neighbors << " " << outlier_limit << "\n";
		//FIRST: I dont know the ratio
		std::vector<T>::iterator first_to_remove
			= CGAL::remove_outliers(points.begin(), points.end(),
				CGAL::Identity_property_map<T>(),
				nb_neighbors,
				100.,                  // No limit on the number of outliers to remove
				outlier_limit * average_spacing); // Point with distance above 2*average_spacing are considered outliers

		/*if (_DEBUG)
		{
			qDebug() << (100. * std::distance(first_to_remove, points.end()) / (double)(points.size()))
				<< "% of the points are considered outliers when using a distance threshold of "
				<< 2. * average_spacing;
		}*/

		points.erase(first_to_remove, points.end());
		
		// SECOND OPTION //
		// I know the ratio of outliers present in the point set

		/*const double removed_percentage = 5.0; // percentage of points to remove

		points.erase(CGAL::remove_outliers(points.begin(), points.end(),
			CGAL::Identity_property_map<Point>(),
			nb_neighbors,
			removed_percentage, // Minimum percentage to remove
			0.), // No distance threshold (can be omitted)
			points.end());
			*/
		std::vector<T>(points).swap(points);
	}

	void SetFirstProperty(const double& first) override { nb_neighbors = first; }
	void SetSecondProperty(const double& second) override { outlier_limit = second; }
	void SetThirdProperty(const double& third) override { average_spacing = third; }

private:
	int nb_neighbors;
	double outlier_limit;	//TODO: better name?
	double average_spacing;
};

//----------------------------------
//----------	------	-----------

/******************************
*	Simplifier componenet     *
******************************/

template <typename T>		//cell size
class GridSimplification : public Simplifier<T>
{
public:
	GridSimplification(const double& cs = 0.001f) :cell_size(cs) {}
	~GridSimplification() {}

	void filter_process(std::vector<T>& points)	override	//if we use mine kd tree then this is not relevant
	{
		qDebug() << "GridSimplification process with " << cell_size << "\n";
		points.erase(CGAL::grid_simplify_point_set(points.begin(), points.end(), cell_size),
			points.end());
		
		std::vector<Point>(points).swap(points);
	}

	void SetFirstProperty(const double& first) { cell_size = first; }

private:
	double cell_size;
};

template <typename T>	//max cluster size, max surface variation
class HiearchySimplification : public Simplifier<T>		//quite good
{
public:
	HiearchySimplification(const double& mcs = 100.0f, const double& msv = 0.01f) : max_cluster_size(mcs), max_surface_variation(msv){}
	~HiearchySimplification() {}

	void filter_process(std::vector<T>& points) override
	{
		qDebug() << "HiearchySimplification process with " << max_cluster_size <<" " <<  max_surface_variation  << "\n";

		points.erase(CGAL::hierarchy_simplify_point_set(points.begin(), points.end(),
			max_cluster_size, 
			max_surface_variation), 
			points.end());

		std::vector<Point>(points).swap(points);
	}

	void SetFirstProperty(const double& mcs) { max_cluster_size = mcs; }
	void SetSecondProperty(const double& msv) { max_surface_variation = msv; }

private:
	double max_cluster_size;
	double max_surface_variation;
};

template <typename T>	//2 double
class WLOPSimplification : public Simplifier<T>			//wow
{
public:
	WLOPSimplification(const double& rp = 2.0f, const double& nr = 0.5f) : retain_percentage(rp), neighbor_radius(nr){}
	~WLOPSimplification() {}

	void filter_process(std::vector<T>& points) override
	{
		qDebug() << "WLOPSimplification process with " << rp << " " << nr << "\n";

		std::vector<Point> output;

		CGAL::wlop_simplify_and_regularize_point_set
			<Concurrency_tag>
			(points.begin(),
				points.end(),
				std::back_inserter(output),
				retain_percentage,
				neighbor_radius
				);

		std::vector<Point>(output).swap(points);
	}

	void SetFirstProperty(const double& rp) { retain_percentage = rp; }
	void SetSecondProperty(const double& nr) { neighbor_radius = nr; }

private:
	double retain_percentage;
	double neighbor_radius;
};

//----------------------------------
//----------	------	-----------

/******************************
*	Smoother componenet      *
******************************/

#ifdef CGAL_EIGEN3_ENABLED
template <typename T>		//neighbours
class JetSmoothing : public Smoother<T>		//Need Eigen, and eigen flag
{
public:
	JetSmoothing(const double& kn = 8) : nb_neighbors(kn){}
	~JetSmoothing() {}

	void filter_process(std::vector<T>& points)
	{
		qDebug() << "JetSmoothing process with " << nb_neighbors << "\n";

		CGAL::jet_smooth_point_set<Concurrency_tag>(points.begin(), points.end(), nb_neighbors);
	}

	void SetFirstProperty(const double& first) { nb_neighbors = first; }

private:
	double nb_neighbors; // default is 24 for real-life point sets
};
#endif

/*template <typename T>		//TODO
class BilateralSmoothing : public Smoother<T>
{
public:
	BilaterSmoothing() {}
	~BilateralSmoothing() {}

	void filter_process(std::vector<T>& points, const double& average_spacing)	//need the normals!!!
	{
		int k = 120;                 // size of neighborhood. The bigger the smoother the result will be.
									 // This value should bigger than 1.
		double sharpness_angle = 25; // control sharpness of the result.
								 // The bigger the smoother the result will be
		int iter_number = 3;         // number of times the projection is applied

		for (int i = 0; i < iter_number; ++i)
		{
			// double error 
			CGAL::bilateral_smooth_point_set <Concurrency_tag>(
				points.begin(),
				points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>(),
				k,
				sharpness_angle);
	}
};*/

#endif

/*
void filter_process(std::vector<T>& points) override
{
std::sort(points.begin(), points.end(), CGALTool::XYZOrder());

float tol = 0.2f;	//0.2f sanyi example sort XYZ , nthelement más

auto it = points.begin() + 1;

float max = 0;

while (it != points.end())
{
float xDif = (it->x() - (it - 1)->x());
if (xDif > max) max = xDif;
if (xDif * xDif > tol)
			{
			break;
			}

++it;
		}

		if (it != points.end())
		{
			if ((it - points.begin()) > (points.end() - it))
				points.erase(it, points.end());
			else
				points.erase(points.begin(), it - 1);

			std::vector<T>(points).swap(points);
		}

	}
*/
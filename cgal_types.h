#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "defines.h"

#ifdef ENABLE_CGAL_SURFACE
#include <CGAL/Polyhedron_3.h>
#include <CGAL/property_map.h>
#endif
//#include <CGAL/Point_with_normal_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;

typedef Kernel::Point_3 Point;
#ifdef ENABLE_CGAL_SURFACE
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef typename Polyhedron::Vertex_const_iterator VCI;
typedef typename Polyhedron::Facet_const_iterator FCI;
typedef typename Polyhedron::Halfedge_around_facet_const_circulator HFCC;
typedef CGAL::Inverse_index<VCI> Index;
#endif

typedef Kernel::Vector_3 Vector;
typedef Kernel::Sphere_3 Sphere;
typedef std::pair<Point, Vector>         Point_with_Normal;
typedef std::vector<Point_with_Normal> PwN_vector;

#ifdef ENABLE_CGAL_SURFACE
typedef CGAL::First_of_pair_property_map<Point_with_Normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_Normal> Normal_map;
#endif

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

struct Color
{
	int r, b, g, a;
};

struct GLPaintFormat
{
	std::vector<Point> p;
	std::vector<unsigned int> ix;
	std::vector<Color> col;
	std::vector<Point> box_points;

	GLPaintFormat() {}
	GLPaintFormat(const std::vector<Point>& _p) :
		p(_p) {}
	GLPaintFormat(const std::vector<Point>& _p, const std::vector<Color>& c) :
		p(_p), col(c) {}
	GLPaintFormat(const std::vector<Point>& _p, const std::vector<unsigned int>& _i, const std::vector<Color>& c) :
		p(_p), ix(_i), col(c) {}
};

namespace CGALTool
{
	struct XYZOrder {
		bool operator()(const Point& c, const Point& _c) {
			return c.x() < _c.x() || (c.x() == _c.x() && c.y() < _c.y()) || (c.x() == _c.x() && c.y() == _c.y() && c.z() < _c.z());
		}
	};
	struct YZXOrder {
		bool operator()(const Point& c, const Point& _c) {
			return c.y() < _c.y() || (c.y() == _c.y() && c.z() < _c.z()) || (c.y() == _c.y() && c.z() == _c.z() && c.x() < _c.x());
		}
	};
	struct ZXYOrder {
		bool operator()(const Point& c, const Point& _c) {
			return c.z() < _c.z() || (c.z() == _c.z() && c.x() < _c.x()) || (c.z() == _c.z() && c.x() == _c.x() && c.y() < _c.y());
		}
	};
};

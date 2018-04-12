#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "defines.h"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Segment_3.h>
#include <QVector3D>
#include <QDebug>

#ifdef ENABLE_CGAL_SURFACE
#include <CGAL/property_map.h>
#endif
//#include <CGAL/Point_with_normal_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;

typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Triangle_3<Kernel> Triangle;
typedef CGAL::Segment_3<Kernel> Segment;

//#ifdef ENABLE_CGAL_SURFACE
typedef typename Polyhedron::Vertex_const_iterator VCI;
typedef typename Polyhedron::Facet_const_iterator FCI;
typedef typename Polyhedron::Halfedge_around_facet_const_circulator HFCC;
typedef CGAL::Inverse_index<VCI> Index;
//#endif

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

template<typename T>
struct triple
{
	triple(const unsigned int& f, const unsigned int& s, const unsigned int& t) : first(f), second(s), third(t) {}

	T first, second, third;

	template<typename T>
	bool operator < (const triple<T>& t) const
	{
		return (first < t.first) || (first == t.first && second < t.second) || (first == t.first && second == t.second && third < t.third);
	}

	template<typename T>
	bool operator <= (const triple<T>& t) const
	{
		return (first <= t.first) || (first == t.first && second <= t.second) || (first == t.first && second == t.second && third <= t.third);
	}
};

typedef std::pair<uint16_t, triple<uint16_t>> RGB;
typedef std::pair<int16_t, std::pair<float, float>> UV;

typedef std::vector<RGB> rbg_map;
typedef std::vector<UV> uv_map;

typedef std::pair<Point, rbg_map> point_rbg_tuple;
typedef std::vector<point_rbg_tuple> point_rbg_map;
typedef std::pair<Point, uv_map> point_uv_tuple;
typedef std::vector<point_uv_tuple> point_uv_map;

typedef CGAL::First_of_pair_property_map<point_uv_tuple>  filter_point_map;

struct GLPaintFormat
{
	std::vector<Point> points;
	std::vector<QVector3D> points_col;
	std::vector<QVector3D> col;
	std::vector<unsigned int> ix;
	std::vector<std::pair<Point, float>> centers_with_radius;
	std::vector<unsigned int> center_part_lengths;
	std::vector<unsigned int> point_part_lengths;
	std::vector<std::pair<float, float>> uv_coords;

	GLPaintFormat() {}
	GLPaintFormat(const std::vector<Point>& _p) :
		points(_p) {
		col.push_back(QVector3D(1, 0, 0));
		points_col.push_back(QVector3D(1, 0, 0));
	}
	GLPaintFormat(const std::vector<Point>& _p, const std::vector<QVector3D>& c) :
		points(_p), col(c) {}
	GLPaintFormat(const std::vector<Point>& _p, const std::vector<unsigned int>& _i, const std::vector<QVector3D>& c) :
		points(_p), ix(_i), col(c) {}
};

namespace CGALTool
{
	struct XYZOrder {
		bool operator()(const point_uv_tuple& c, const point_uv_tuple& _c) {
			return c.first.x() < _c.first.x() || (c.first.x() == _c.first.x() && c.first.y() < _c.first.y()) || (c.first.x() == _c.first.x() && c.first.y() == _c.first.y() && c.first.z() < _c.first.z());
		}
	};
	struct YZXOrder {
		bool operator()(const point_uv_tuple& c, const point_uv_tuple& _c) {
			return c.first.y() < _c.first.y() || (c.first.y() == _c.first.y() && c.first.z() < _c.first.z()) || (c.first.y() == _c.first.y() && c.first.z() == _c.first.z() && c.first.x() < _c.first.x());
		}
	};
	struct ZXYOrder {
		bool operator()(const point_uv_tuple& c, const point_uv_tuple& _c) {
			return c.first.z() < _c.first.z() || (c.first.z() == _c.first.z() && c.first.x() < _c.first.x()) || (c.first.z() == _c.first.z() && c.first.x() == _c.first.x() && c.first.y() < _c.first.y());
		}
	};
};

inline std::ostream& operator << (std::ostream& out, const Point& p)
{
	out << p.x() << ", " << p.y() << ", " << p.z();
	return out;
}

inline std::ostream& operator << (std::ostream& out, const point_uv_tuple& puv_m)
{
	out << puv_m.first << " : ";
	for (auto & tup : puv_m.second)
	{
		out << tup.first << " [ " << tup.second.first << ", " << tup.second.second << " ] ,";
	}
	return out;
}

inline QDebug operator << (QDebug& out, const Point& p)
{
	out << p.x() << ", " << p.y() << ", " << p.z();
	return out;
}

inline std::pair<float, float>& operator * (std::pair<float, float>& lhs, float& rhs)
{
	lhs = std::pair<float, float>(lhs.first * rhs, lhs.second * rhs);
	return lhs;
}

inline std::pair<float, float>& operator + (std::pair<float, float>& lhs, std::pair<float, float>& rhs)
{
	lhs = std::pair<float, float>(lhs.first + rhs.first, lhs.second + rhs.second);
	return lhs;
}

inline std::pair<float, float>& operator / (std::pair<float, float>& lhs, std::pair<float, float>& rhs)
{
	lhs = std::pair<float, float>(lhs.first / rhs.first, lhs.second / rhs.second);
	return lhs;
}

inline Point& operator +(const Point& lhs, const Point& rhs)
{
	Point res = Point(lhs.x() + rhs.x(), lhs.y() + rhs.y(), lhs.z() + rhs.z());
	return res;
}

inline Point& operator / (Point& lhs, const float& rhs)
{
	lhs = Point(lhs.x() / rhs, lhs.y() / rhs, lhs.z() / rhs);
	return lhs;
}

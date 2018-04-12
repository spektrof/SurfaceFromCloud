#pragma once

#include "types.h"
#include "cgal_types.h"
#include "diagram_types.h"

#include <QDebug>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

typedef CGAL::Triangulation_cell_base_with_info_3<cell_inf, Kernel, CGAL::Triangulation_cell_base_3<Kernel>> delaunay_triangulation_cell_base_with_info;

typedef CGAL::Triangulation_data_structure_3<
	CGAL::Triangulation_vertex_base_3<Kernel>,
	delaunay_triangulation_cell_base_with_info,
	Concurrency_tag>                          delaunay_triangulation_data_structure;
typedef CGAL::Delaunay_triangulation_3<Kernel, delaunay_triangulation_data_structure> delaunay_triangulation;

typedef delaunay_triangulation::Finite_vertices_iterator delaunay_finite_vertices_iterator;
typedef delaunay_triangulation::Finite_cells_iterator delaunay_finite_cells_iterator;
typedef delaunay_triangulation::Cell_handle delaunay_cell_handle;
typedef delaunay_triangulation::Vertex_handle delaunay_vertex_handle;

//1 Cellahoz 1 surf pont és 2 pole tartozik (másik pont, sugárral)

class VoronoiCell
{
public:
	VoronoiCell(Point& sp) : surf_point(sp) {
		cell_edges.clear();
	}
	~VoronoiCell() {
		cell_edges.clear();
	}

	GLPaintFormat getPaintData()
	{
		GLPaintFormat res;

		int index = 0;
		std::map<Point*, unsigned int> tmp_points;
		for (auto& it : cell_edges)
		{ 
			std::pair<std::map<Point*, unsigned int>::iterator, bool> pos;
			pos = tmp_points.insert(std::pair<Point*, unsigned int>(it->first.second, index));
			if (pos.second == true)
			{
				res.points.push_back(*pos.first->first);
				index++;
			}

			res.ix.push_back(pos.first->second);

			pos = tmp_points.insert(std::pair<Point*, unsigned int>(it->second.second, index));
			if (pos.second == true)
			{
				res.points.push_back(*pos.first->first);
				index++;
			}

			res.ix.push_back(pos.first->second);
		}

		/*for (auto it : cell_edges)
		{
			res.ix.push_back(it.first.first);
			res.ix.push_back(it.second.first);
		}*/
		res.centers_with_radius.push_back(std::pair<Point, float>(surf_point, 0.05f));

		if (!poles[0].isNull())	res.centers_with_radius.push_back(std::pair<Point, float>(*poles[0].center, poles[0].radius));
		if (!poles[1].isNull()) res.centers_with_radius.push_back(std::pair<Point, float>(*poles[1].center, poles[1].radius));

		res.center_part_lengths.push_back(1);
		res.center_part_lengths.push_back(res.centers_with_radius.size() - 1);
		res.col.push_back(QVector3D(148.0f / 255.0f, 0, 211.0f / 255.0f));
		res.col.push_back(QVector3D(0, 0, 1));
		res.points_col.push_back(QVector3D(1, 0, 0));

		return res;
	}

	GLPaintFormat getPaintDataWithCellD(const unsigned int& ind)
	{
		return GLPaintFormat();
	}

	void getPolesPaintData(std::vector<std::pair<Point, float>>& _poles) const
	{
		if (!poles[0].isNull())	_poles.push_back(std::pair<Point, float>(*poles[0].center, poles[0].radius));
		if (!poles[1].isNull()) _poles.push_back(std::pair<Point, float>(*poles[1].center, poles[1].radius));
	}

	std::vector<Point> getPoints()
	{
		std::set<Point> tmp;
		for (auto& it : cell_edges)
		{
			tmp.insert(*it->first.second);
			tmp.insert(*it->second.second);
		}
		return std::vector<Point>(tmp.begin(),tmp.end());
	}

	std::pair<Pole, Point> getFirstPole()
	{
		return std::pair<Pole, Point>(poles[0], surf_point);
	}

	std::pair<Pole, Point> getSecondPole()
	{
		return std::pair<Pole, Point>(poles[1], surf_point);
	}

	Point getSurfPoint() const { return surf_point; }

	void calcPoles()
	{
		std::set<Point*> tmp;
		for (auto& it : cell_edges)
		{
			tmp.insert(it->first.second);
			tmp.insert(it->second.second);
		}
		auto vertices = std::vector<Point*>(tmp.begin(), tmp.end());
		//-----------------
		std::vector<Point*>::iterator farthest = vertices.begin();
		std::vector<Point*>::iterator it = vertices.begin()++;
		float distance = CGAL::squared_distance(surf_point, **farthest);

		for (; it != vertices.end(); ++it)
		{
			float tmp_distance = CGAL::squared_distance(surf_point, **it);
			if (tmp_distance > distance)
			{
				distance = tmp_distance;
				farthest = it;
			}
		}

		//qDebug() << "\tFirst pole is (" << farthest->second << "): " << (farthest->first)->x() << "," << (farthest->first)->y() << "," << (farthest->first)->z() << " radius_square: " << distance << "\n";

		poles[0].center = *farthest;
		poles[0].radius = sqrt(distance);

		CGAL::Vector_3<Kernel> sp1 = **farthest - surf_point;
		farthest = vertices.end();
		float scalar_dot = 1.0f;

		//choose the most negative product
		it = vertices.begin();
		for (; it != vertices.end(); ++it)
		{
			CGAL::Vector_3<Kernel> sp2 = **it - surf_point;
			float tmp_scalar_dot = CGAL::scalar_product(sp1, sp2);
			if (scalar_dot > tmp_scalar_dot)
			{
				farthest = it;
				scalar_dot = tmp_scalar_dot;
			}
		}

		/*distance = 0.0f;
		it = vertices_map.begin();
		for (; it != vertices_map.end(); ++it)
		{
			CGAL::Vector_3<Kernel> sp2 = *it->first - surf_point;
			float tmp_scalar_dot = CGAL::scalar_product(sp1, sp2);
			float tmp_distance = CGAL::squared_distance(surf_point, *it->first);
			if  ((tmp_scalar_dot < 0 && distance < tmp_distance) || (scalar_dot > tmp_scalar_dot && distance == tmp_distance))
			{
				farthest = it;
				scalar_dot = tmp_scalar_dot;
				distance = tmp_distance;
			}
		}*/

		if (scalar_dot >= 0)
		{
			if (scalar_dot > 0) qDebug() << "ERROR: COULDNT FIND second pole" << scalar_dot << "\n";
			else qDebug() << "ERROR: COULDNT FIND second pole because scalar dot is 0\n";

			it = vertices.begin();
			for (; it != vertices.end(); ++it)
			{
				CGAL::Vector_3<Kernel> sp2 = **it - surf_point;
				qDebug() << "\t" << CGAL::scalar_product(sp1, sp2);
			}

			poles[1].center = NULL;
		}
		else
		{
			distance = CGAL::squared_distance(surf_point, **farthest);
			//qDebug() << "\tSecond pole is (" << farthest->second << "): " << (farthest->first)->x() << "," << (farthest->first)->y() << "," << (farthest->first)->z() << " radius_square: " << distance << "\n";
			poles[1].center = *farthest;
			poles[1].radius = sqrt(distance);
		}
	}

	void addVoronoiSegment(edge_t* edge)
	{
		cell_edges.insert(edge);
	}

private:
	Point surf_point;

	Pole poles[2];

	std::set<edge_t*> cell_edges;
};

class VoronoiDiagram
{
public:
	VoronoiDiagram() : bounded(Box(-1.0f,1.0f))
	{
		cells.clear();
		voronoi_edges.clear();
		voronoi_vertices.clear();
		indicies_of_delauney_segments.clear();
	}

	~VoronoiDiagram() {
		cells.clear();
		voronoi_edges.clear();
		voronoi_vertices.clear();
		indicies_of_delauney_segments.clear();
	}

	typedef std::vector<VoronoiCell>::iterator Voronoi_cell_iterator;
	typedef std::vector<VoronoiCell>::const_iterator const_Voronoi_cell_iterator;

	Voronoi_cell_iterator cells_begin()
	{
		return cells.begin();
	}

	Voronoi_cell_iterator cells_end()
	{
		return cells.end();
	}

	VoronoiCell* getCell(const int& ind)
	{
		return &cells[ind];
	}

	void addCell(const VoronoiCell& vc)
	{
		cells.push_back(vc);
	}

	void clear() {
		cells.clear();
		voronoi_edges.clear();
		voronoi_vertices.clear();
		indicies_of_delauney_segments.clear();
	}

	size_t size() const { return cells.size(); }

	void setBox(const Box& b) { bounded = b; }

	std::vector<std::pair<Pole, Point> > getPoles()
	{
		std::vector<std::pair<Pole, Point>> tmp;
		for (auto& it : cells)
		{
			tmp.push_back(it.getFirstPole());
			auto tmp_pole_check = it.getSecondPole();
			if (!tmp_pole_check.first.isNull())
				tmp.push_back(tmp_pole_check);
			else
				qDebug() << "Couldnt add the second pole!";
		}
		return tmp;
	}

	GLPaintFormat getCellPaintData(const int& ind);
	GLPaintFormat getCellWDualPaintData(const int& ind, const int& ind2);
	GLPaintFormat getPolesSurfPaintData();
	GLPaintFormat getDelauneySegmentPaintData(GLPaintFormat& res);
	GLPaintFormat getVoronoiDiagramBySegmentPaintData();

	/*Precondition: set the bounding box, default: -1.0f,1.0f range on every axis*/
	void calc_diagram(std::vector<Point>& points);

	/*Precondition: calc_diagram*/
	void calc_cell_poles()
	{
		if (cells.size() == 0)
		{
			qDebug() << "Call the calc_diagram first! There is no voronoi cell.";
			return;
		}

		for (auto &it : cells)
			it.calcPoles();

		qDebug() << "Poles are calculated\n";
	}

	void calc_diagram_with_poles(std::vector<Point>& points)
	{
		calc_diagram(points);

		if (cells.size() == 0) return;
	
		calc_cell_poles();
	}

	void calc_test()
	{
		qDebug() << "point generatin started";
		std::vector<Point> points_w;

		int counter = 0;
		for (int i = 0; i < 440000; ++i)
		{
			Point tmp = Point((-1.0f / 140000.0f) * (float)i, 0.5f + (-1.0f / 140000.0f) * (float)i, (1.0f / 140000.0f)* (float)i);
			//	qDebug() << tmp.x() << ", " << tmp.y() << " , " << tmp.z();
			if (counter % 10000 == 0) qDebug() << "10k";
			points_w.push_back(tmp);
			counter++;
		}
		qDebug() << "Dela started";
		test = delaunay_triangulation(points_w.begin(), points_w.end());
		qDebug() << "TEST DELA have " << test.number_of_cells() << " " << test.number_of_finite_cells() << " cell " << test.number_of_edges() << " edge " << test.number_of_facets() << " face " << test.number_of_vertices();
		//0   0  cell  1319997  edge  879998  face  440000
		//170mb
	}

protected:
	void addBoxPoints(std::set<Point>& box_points, std::vector<Point>& points);
	void calcDelauneySegments(const delaunay_triangulation& T, const std::vector<Point>&);
	
	edge_t* insert_edge(std::map<Point*, edge_t*>& edge_list_by_cell, edge_t& edge_candidate);
	void set_cell_info(const delaunay_triangulation& T, delaunay_cell_handle& cell, cell_inf* info, unsigned int& vor_index);

private:
	delaunay_triangulation test;

	typedef typename _vertex::v_ptr v_ptr;
	typedef typename edge_p::e_ptr e_ptr;

	std::vector<VoronoiCell> cells;

	std::vector<v_ptr> voronoi_vertices;
	std::vector<e_ptr> voronoi_edges;
	//std::set<edge_t> voronoi_edges;

	std::vector<unsigned int> indicies_of_delauney_segments;

	Box bounded;
};


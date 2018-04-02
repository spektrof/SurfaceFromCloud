#pragma once

#include "types.h"
#include "cgal_types.h"
#include "diagram_types.h"

#include <QDebug>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
//#include <CGAL/convex_hull_3.h>

//typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, Kernel>    Vb;

typedef CGAL::Triangulation_data_structure_3<
	CGAL::Triangulation_vertex_base_3<Kernel>,
	CGAL::Triangulation_cell_base_3<Kernel>,
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
	VoronoiCell(const Point& sp, const unsigned int& spi) : surf_point(sp), surf_point_index(spi) {
		vertices_map.clear();
	//	voronoi_vertices.clear();
		indicies_of_cell_segments.clear();
	}
	~VoronoiCell() {
	}

	void addVoronoiVertex(Point* p, const unsigned int& ind)
	{
		vertices_map.insert(std::pair<Point*, unsigned int>(p, ind));
	}

	void setSegmentIndicies(const std::set<std::pair<Point*, Point*>>& segments)
	{
		std::map<Point*, unsigned int> tmp_vert_index_pairs;
		int index = 0;

		for (auto& it : vertices_map)
		{
		//	voronoi_vertices.push_back(it.first);
			tmp_vert_index_pairs.insert(std::pair <Point*, unsigned int>( it.first, index));
			index++;
		}

		for (auto it : segments)
		{
			Point* left = it.first;
			Point* right = it.second;

			auto left_pos = tmp_vert_index_pairs[left];
			auto right_pos = tmp_vert_index_pairs[right];
			
			indicies_of_cell_segments.push_back(left_pos);
			indicies_of_cell_segments.push_back(right_pos);
		}
	}

	GLPaintFormat getPaintData()
	{
		GLPaintFormat res;

		for (auto& it : vertices_map)
			res.points.push_back(*it.first);

		res.ix = indicies_of_cell_segments;
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
		if (inc_cell_w_dual.size() == 0) return GLPaintFormat();

		GLPaintFormat res;

		for (auto& it : vertices_map)
			res.points.push_back(*it.first);

		res.ix = indicies_of_cell_segments;
		res.centers_with_radius.push_back(std::pair<Point, float>(surf_point, 0.05f));

		res.center_part_lengths.push_back(1);
		res.center_part_lengths.push_back(1);

		res.col.push_back(QVector3D(148.0f / 255.0f, 0, 211.0f / 255.0f));
		res.col.push_back(QVector3D(0, 1, 1));

		res.points_col.push_back(QVector3D(1, 0, 0));
		res.points_col.push_back(QVector3D(0, 0, 0));

		auto pair = inc_cell_w_dual[ind % inc_cell_w_dual.size()];
		res.centers_with_radius.push_back(std::pair<Point, float>(pair.second, 0.05f));


		for (auto it : pair.first)
		{
			res.points.push_back(it);
		}

		int i1 = res.points.size() - 4;
		int i2 = res.points.size() - 3;
		int i3 = res.points.size() - 2;
		int i4 = res.points.size() - 1;

		res.point_part_lengths.push_back(res.ix.size());

		res.ix.push_back(i1); res.ix.push_back(i2);
		res.ix.push_back(i2); res.ix.push_back(i3);
		res.ix.push_back(i3); res.ix.push_back(i1);
		res.ix.push_back(i1); res.ix.push_back(i4);
		res.ix.push_back(i2); res.ix.push_back(i4);
		res.ix.push_back(i3); res.ix.push_back(i4);

		res.point_part_lengths.push_back(res.ix.size() - res.point_part_lengths[0]);

		return res;
	}

	void getPolesPaintData(std::vector<std::pair<Point, float>>& _poles) const
	{
		if (!poles[0].isNull())	_poles.push_back(std::pair<Point, float>(*poles[0].center, poles[0].radius));
		if (!poles[1].isNull()) _poles.push_back(std::pair<Point, float>(*poles[1].center, poles[1].radius));
	}

	std::vector<Point> getPoints()
	{
		std::vector<Point> tmp;
		for (auto& it : vertices_map)
			tmp.push_back(*it.first);
		return tmp;
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

	size_t vert_size() { return vertices_map.size(); }

	void calcPoles()
	{
		std::map<Point*, unsigned int>::iterator farthest = vertices_map.begin();
		std::map<Point*, unsigned int>::iterator it = vertices_map.begin()++;
		float distance = CGAL::squared_distance(surf_point, *farthest->first);

		for (; it != vertices_map.end(); ++it)
		{
			float tmp_distance = CGAL::squared_distance(surf_point, *it->first);
			if (tmp_distance > distance)
			{
				distance = tmp_distance;
				farthest = it;
			}
		}

		//qDebug() << "\tFirst pole is (" << farthest->second << "): " << (farthest->first)->x() << "," << (farthest->first)->y() << "," << (farthest->first)->z() << " radius_square: " << distance << "\n";

		poles[0].center = farthest->first;
		poles[0].radius = sqrt(distance);

		CGAL::Vector_3<Kernel> sp1 = *farthest->first - surf_point;
		farthest = vertices_map.end();
		float scalar_dot = 1.0f;

		//choose the most negative product
		it = vertices_map.begin();
		for (; it != vertices_map.end(); ++it)
		{
			CGAL::Vector_3<Kernel> sp2 = *it->first - surf_point;
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

			it = vertices_map.begin();
			for (; it != vertices_map.end(); ++it)
			{
				CGAL::Vector_3<Kernel> sp2 = *it->first - surf_point;
				qDebug() << "\t" << CGAL::scalar_product(sp1, sp2);
			}

			poles[1].center = NULL;
		}
		else
		{
			distance = CGAL::squared_distance(surf_point, *farthest->first);
			//qDebug() << "\tSecond pole is (" << farthest->second << "): " << (farthest->first)->x() << "," << (farthest->first)->y() << "," << (farthest->first)->z() << " radius_square: " << distance << "\n";
			poles[1].center = farthest->first;
			poles[1].radius = sqrt(distance);
		}
	}

	void addIncCell(std::vector<Point>& po, Point& d)
	{
		inc_cell_w_dual.push_back(std::pair<std::vector<Point>, Point>(po,d));
	}

private:
	Point surf_point;
	unsigned surf_point_index;		//mybe this is unneceserry

	std::map<Point*, unsigned int> vertices_map;	//with global indicies, only if we want to draw all
//	std::vector<Point*> voronoi_vertices;		//TODO delete if there is no error after TEST

	Pole poles[2];

	std::vector<unsigned int> indicies_of_cell_segments;
	std::vector<std::pair<std::vector<Point>, Point>> inc_cell_w_dual;	//TODO remove later, only for check
};

class VoronoiDiagram
{
public:
	VoronoiDiagram() : bounded(Box(-1.0f,1.0f))
	{
		cells.clear();
		point_map.clear();
		indicies_of_delauney_segments.clear();
		indicies_of_vornoi_diagram_segments.clear();
		points_with_box_points.clear();
	}

	~VoronoiDiagram() {
		cells.clear();
		point_map.clear();
		indicies_of_delauney_segments.clear();
		indicies_of_vornoi_diagram_segments.clear();
		points_with_box_points.clear();
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
		point_map.clear();
		indicies_of_delauney_segments.clear();
		indicies_of_vornoi_diagram_segments.clear();
		points_with_box_points.clear();
	}

	size_t size() const { return cells.size(); }
	size_t size_of_points() const { return point_map.size(); }

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

	std::vector<Point> getCellPoints(const int& ind)
	{
		return cells[ind].getPoints();
	}

	GLPaintFormat getCellPaintData(const int& ind);
	GLPaintFormat getCellWDualPaintData(const int& ind, const int& ind2);
	GLPaintFormat getPolesSurfPaintData();
	GLPaintFormat getDelauneySegmentPaintData();
	GLPaintFormat getVoronoiDiagramBySegmentPaintData();

	/*Precondition: set the bounding box, default: -1.0f,1.0f range on every axis*/
	void calc_diagram(std::vector<Point>& points)
	{
		qDebug() << "Voronoi diagram calculation";
		/*************************************
		*	Clear	*
		**************************************/
		cells.clear();
		point_map.clear();
		indicies_of_delauney_segments.clear();
		indicies_of_vornoi_diagram_segments.clear();
		points_with_box_points.clear();

		qDebug() << "\tOUR BOX IS: " << bounded.getXMin() << " " << bounded.getXMax() << " " << bounded.getYMin() << " " << bounded.getYMax() << " " << bounded.getZMin() << " " << bounded.getZMax();

		/*************************************
		*	Triangulation	*
		**************************************/

		delaunay_triangulation::Lock_data_structure locking_ds(CGAL::Bbox_3(bounded.getXMin(), bounded.getYMin(), bounded.getZMin(), bounded.getXMax(), bounded.getYMax(), bounded.getZMax()), 50);

		points_with_box_points = points;
		std::set<Point> box_points;

		addBoxPoints(box_points, points_with_box_points);

		delaunay_triangulation T(points_with_box_points.begin(), points_with_box_points.end(), &locking_ds);
		assert(T.is_valid());
		qDebug() << "\t" << T.number_of_cells() << " cells " << T.number_of_finite_cells() << " finite cells";
	
		calcDelauneySegments(T);
		/*************************************
		*	Voronoi Diagram calculation	*
		**************************************/

		unsigned int surf_index = 0;
		unsigned int vor_index = 0;
		std::set<std::pair<Point*, Point*>> voronoi_diagram_segments;

		for (delaunay_triangulation::Finite_vertices_iterator vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++)
		{
			if (box_points.count(vit->point()) == 1)	//its a bounding box point, dont need their cells
				continue;

			//Getting Incident cells
			std::vector<delaunay_cell_handle> inc_cells;
			std::vector<delaunay_cell_handle> o_inc_cells;
			T.incident_cells(vit, std::back_inserter(inc_cells));
			T.finite_incident_cells(vit, std::back_inserter(o_inc_cells));

			if (inc_cells.size() != o_inc_cells.size())
				qDebug() << "INCIDENT CELL SIZES DIFF " << o_inc_cells.size()- inc_cells.size();
			
			if (inc_cells.empty())
			{
				qDebug() << "EMPPTY CELLS FOR SRUF VERT\n";
				continue;
			}

			//Create the current voronoi cell
			VoronoiCell vc = VoronoiCell(vit->point(), surf_index);
			surf_index++;
		
			//calculate voronoi vertices
			std::map<delaunay_cell_handle, Point*> cell_dual_pairs;			//calculate the duals for every incident cells

			for (auto cit = inc_cells.begin(); cit != inc_cells.end(); ++cit)
			{
				auto related_dual =/* getBoundedDual(*/T.dual(*cit)/*)*/;

				std::pair<std::map<Point, unsigned int>::iterator, bool> res;
				res = point_map.insert(std::pair<Point, unsigned int >(related_dual, vor_index));	//put dual into the point map
				if (res.second == true)	//not in the map
					vor_index++;	

									//add the vertex to the voronoi cell
				vc.addVoronoiVertex(const_cast<Point*>(&(res.first->first)), res.first->second);
				cell_dual_pairs.insert(std::pair<delaunay_cell_handle, Point*>(*cit, const_cast<Point*>(&(res.first->first))));

				std::vector<Point> po;
				po.push_back((*cit)->vertex(0)->point());
				po.push_back((*cit)->vertex(1)->point());
				po.push_back((*cit)->vertex(2)->point());
				po.push_back((*cit)->vertex(3)->point());
				vc.addIncCell(po, related_dual);

			}
			//-----------------------------------------------
			//Calculate Segments
			std::set<std::pair<Point*, Point*>> act_cell_segments;
			
			for (auto cit = inc_cells.begin(); cit != inc_cells.end(); ++cit)
			{
				for (int i = 0; i < 4; ++i)
				{
					delaunay_cell_handle tmp_cell = (*cit)->neighbor(i);
					auto is_related = cell_dual_pairs.find(tmp_cell);

					if (is_related != cell_dual_pairs.end())	//both cell_handle are related to the actual point
					{
						//segments of cell
						auto already_in = act_cell_segments.find(std::pair<Point*, Point*>(cell_dual_pairs[tmp_cell], cell_dual_pairs[*cit]));
						if (already_in == act_cell_segments.end())
							act_cell_segments.insert(std::pair<Point*, Point*>(cell_dual_pairs[*cit], cell_dual_pairs[tmp_cell]));

						//ALL SEGMENTS WHAT ARE RELATED TO THE CELLS
						already_in = voronoi_diagram_segments.find(std::pair<Point*, Point*>(cell_dual_pairs[tmp_cell], cell_dual_pairs[*cit]));
						if (already_in == voronoi_diagram_segments.end())
							voronoi_diagram_segments.insert(std::pair<Point*, Point*>(cell_dual_pairs[*cit], cell_dual_pairs[tmp_cell]));
					}
					//So we wont put those which go to infinite, TODO mybe
				}
			}

			vc.setSegmentIndicies(act_cell_segments);	//cell indicies for drawing
			addCell(vc);
		}

		/*************************************
		*	Calculate The Segment Indicies for Drawing	*
		**************************************/

		setDiagramIndicies(voronoi_diagram_segments);

		/*************************************
		*	END - calc_diagram	
		
			Ready for calcule poles*
		**************************************/
	}

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

protected:

	Point getBoundedDual(Point& dual);
	void addBoxPoints(std::set<Point>& box_points, std::vector<Point>& points);
	void calcDelauneySegments(const delaunay_triangulation& T);

	void setDiagramIndicies(std::set<std::pair<Point*, Point*>> voronoi_diagram_segments)
	{
		std::map<Point*, unsigned int> tmp_vertex_order; //TODO: rework with the private one
		int index = 0;
		for (auto& it : point_map)
		{
			tmp_vertex_order.insert(std::pair<Point*, unsigned int>(const_cast<Point*>(&it.first),index));
			index++;
		}

		indicies_of_vornoi_diagram_segments.clear();

		for (auto it : voronoi_diagram_segments)
		{
			auto left = tmp_vertex_order[it.first];
			auto right = tmp_vertex_order[it.second];

			indicies_of_vornoi_diagram_segments.push_back(left);
			indicies_of_vornoi_diagram_segments.push_back(right);
		}
	}

private:
	std::vector<VoronoiCell> cells;
	std::map<Point, unsigned int> point_map;	//voronoi vertex points
	std::vector<Point> points_with_box_points;	//surface points with its bound box

	std::vector<unsigned int> indicies_of_delauney_segments;
	std::vector<unsigned int> indicies_of_vornoi_diagram_segments;

	Box bounded;
};


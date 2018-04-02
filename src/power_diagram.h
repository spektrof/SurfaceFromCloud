#pragma once

#include "types.h"
#include "cgal_types.h"
#include "diagram_types.h"
#include "priority_queue.h"

#include <QDebug>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
//#include <CGAL/convex_hull_3.h>
#include <CGAL/Linear_algebraHd.h>

typedef Kernel::Weighted_point_3                            weighted_point;

typedef CGAL::Regular_triangulation_vertex_base_3<Kernel>        regular_triangulation_vertex_base;
typedef CGAL::Triangulation_vertex_base_with_info_3<std::vector<Point>, Kernel, regular_triangulation_vertex_base> regular_triangulation_vertex_base_with_info;
typedef CGAL::Regular_triangulation_cell_base_3<Kernel>          regular_triangulation_cell_base;
typedef CGAL::Triangulation_data_structure_3<
	regular_triangulation_vertex_base_with_info,
	regular_triangulation_cell_base,
	Concurrency_tag
> regular_triangulation_data_structure;

typedef CGAL::Regular_triangulation_3<Kernel, regular_triangulation_data_structure> regular_triangulation;

typedef regular_triangulation::Cell_handle regular_cell_handle;

typedef CGAL::Linear_algebraHd<float> l_algebra;
typedef CGAL::Linear_algebraHd<float>::Matrix Matrix;

class PolarBall
{
	enum flag
	{
		INNER,
		OUTER,
		UNKNOWN
	};

public:
	PolarBall() {}
	PolarBall(Pole& _pole, const std::vector<Point>& sur_p) : pole(_pole), surf_points(sur_p), label(UNKNOWN), in (0.0f), out(0.0f) {}
	~PolarBall() {}

	void setPole(Point* _p, const float& r) { pole.center = _p; pole.radius = r; }

	Pole getPole() const { return pole; }
	Point getPoint() const { return *pole.center; }
	float getRadius() const { return pole.radius; }
	std::vector<Point> getSurfPoints() const { return surf_points; }
	bool isEmptyPole() const { return pole.isNull(); }
	bool isInnerPole() const { return label == INNER; }
	bool isOuterPole() const { return label == OUTER; }
	bool isUnkown() const { return label == UNKNOWN; }

	float getInValue() const { return in; }
	float getOutValue() const { return out; }
	void  setInValue(const float& _in) { in = _in; }
	void  setOutValue(const float& _out) { out = _out; }

	void setFlag(const unsigned& ind) { if (ind > 2) return;  label = flag(ind); }
	unsigned int getFlag() const { return label;  }

	unsigned int getOppositeFlag() const { return label == 2 ? 2 : 1 - label; }
	float getValueByFlag(const unsigned int& flag) const 
	{
		if (flag == 2) return -1.0f;

		return flag == 0 ? getInValue() : getOutValue();
	}

	void setValueByFlag(const unsigned int& flag, const float& new_val)
	{
		if (flag == 2) return;

		if (flag == 0)
			setInValue(new_val);
		else
			setOutValue(new_val);
	}

	float getAlphaWeight(PolarBall* rhs);					//intersect angle of shallowly intersect ball
	float getBetaWeight(PolarBall* rhs, Point& surfpoint);	//angle related to a common surfpoint

private:
	Pole pole;
	std::vector<Point> surf_points;
	flag label;

	float in, out;	//should be between 0.0f and 1.0f
};

class PowerCell
{
	
public:
	PowerCell(Point* pol_p, const float& r, const std::vector<Point>& sur_p) : polar_ball(PolarBall(Pole(pol_p, r), sur_p)) {
		vertices_map.clear();
		neighbour_candidates.clear();
		proper_neighbours.clear();
		indicies_of_cell_segments.clear();
	}
	~PowerCell() {
		vertices_map.clear();
		neighbour_candidates.clear();
		proper_neighbours.clear();
		indicies_of_cell_segments.clear();
	}

	void addPowerVertex(Point* p, const unsigned int& ind)
	{
		vertices_map.insert(std::pair<unsigned int, Point*>(ind, p));
	}

	void setSegmentIndicies(const std::set<std::pair<Point*, Point*>>& segments)
	{
		std::map<Point*, unsigned int> tmp_vert_index_pairs;
		int index = 0;

		for (auto& it : vertices_map)
		{
		//	power_vertices.push_back(it.second);
			tmp_vert_index_pairs.insert(std::pair <Point*, unsigned int>(it.second, index));
			index++;
		}

		indicies_of_cell_segments.clear();

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

	std::vector<unsigned int> getSegmentIndicies() const {
		return indicies_of_cell_segments;
	}

	GLPaintFormat getPaintData()
	{
		GLPaintFormat res;

		for (auto& it : vertices_map)
			res.points.push_back(*it.second);

		//res.ix = local_draw_indicies;
		res.ix = indicies_of_cell_segments;
		res.centers_with_radius.push_back(std::pair<Point, float>(polar_ball.getPoint(), polar_ball.getRadius()));
		//res.surf_centers.push_back(polar_ball.getSurfPoint());

		res.col.push_back(QVector3D(0, 0, 1));
		res.center_part_lengths.push_back(res.centers_with_radius.size());

	//	qDebug() << "Voronoi vertices: " << power_vertices.size() << "\n Indicies: " << local_draw_indicies.size() << "\n";
	//	qDebug() << "\t 1 polar ball" << res.centers_with_radius.size() << "\n";
		return res;
	}

	int getNeighbourIndex(const int& ind) const
	{
		if (proper_neighbours.empty() ||  ind < 0) return -1;
		return proper_neighbours[ind % proper_neighbours.size()];
	}
	
	void appendPaintDataAsNeighbour(GLPaintFormat& res, const unsigned int& related)
	{
	//	if (related < 0 || neighbour_candidates.size() <= related) return;

		const size_t pre_size = res.points.size();
		for (auto& it : vertices_map)
			res.points.push_back(*it.second);

		for (auto it : indicies_of_cell_segments)
			res.ix.push_back(it + pre_size);

		qDebug() << "Neighbour vertices: " << vertices_map.size() /*<< "\n Indicies: " << local_draw_indicies.size() */<< "\n";
		for (auto it2 : neighbour_candidates[related])
		{
			qDebug() << it2 << ": " << qSetRealNumberPrecision(10) << vertices_map[it2]->x() <<", " << vertices_map[it2]->y() << ", " << vertices_map[it2]->z();
			res.centers_with_radius.push_back(std::pair<Point, float>(*vertices_map[it2], 0.005f));
		}
		res.col.push_back(QVector3D(0, 0, 0));
		res.center_part_lengths.push_back(neighbour_candidates[related].size());
	}

	Pole getPolePoint() const { return polar_ball.getPole(); }

	std::vector<Point> getPoints()
	{
		std::vector<Point> tmp;
		for (auto& it : vertices_map)
			tmp.push_back(*it.second);
		return tmp;
	}

	Point getPointFromMap(const unsigned int& ind)
	{
		return *vertices_map[ind];
	}

	Point* getPointByIndex(const unsigned int& ind)
	{
		if (ind >= vertices_map.size() || ind < 0)
			return NULL;
		return vertices_map[ind];
	}

	bool hasInnerPole() const { return polar_ball.isInnerPole(); }
	PolarBall* getPolarBallPtr() { return &polar_ball; }
	std::vector<unsigned int> getProperNeighbours() const { return proper_neighbours; }
	std::vector<Point> getRelatedSurfPoints() const { return polar_ball.getSurfPoints(); }

	size_t vert_size() { return vertices_map.size(); }

	void printNeighbourCandidate()
	{
		qDebug() << "NEIGHBOUR CANDIDATE\n";
		for (auto it : neighbour_candidates)
		{
			qDebug() << "-----------------------" <<
					 it.first << ":\n";
			for (auto it2 : it.second)
				qDebug() << "\t" << it2;

		}
		qDebug() << "\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\n";

	}

	void addNeighbourCandidate(const unsigned int& related_power_cell_index, const unsigned int& actual_power_vertex_index)
	{
		std::pair<std::map<unsigned int, std::set<unsigned int>>::iterator, bool> res;
		res = neighbour_candidates.insert(std::pair< unsigned int, std::set<unsigned int> >(related_power_cell_index, std::set<unsigned int>()));
		res.first->second.insert(actual_power_vertex_index);
	}

	std::set<unsigned int> getCommonNeighbourPoints(const unsigned int& neigh)
	{
		if (neigh < 0) return std::set<unsigned int>();

		return neighbour_candidates[neigh];
	}

	void printProperNeightbours()
	{
		std::stringstream out;
		out << "\t" << " My neighbours (" << proper_neighbours.size() << "): ";
		for (auto it : proper_neighbours)
		{
			out << " " << it;
		}
		qDebug() << QString(out.str().c_str());
	}

	void calcProperNeighbours();

private:
	PolarBall polar_ball;

	std::map<unsigned int, Point*> vertices_map;	//power verticies with their indicies
	std::map<unsigned int, std::set<unsigned int>> neighbour_candidates;
			// global Power Cell index, related common global vertex index
	std::vector<unsigned int> proper_neighbours;

	//std::vector<Point*> power_vertices;		//local order, bcuse convex hull's verices order is different than in map!! - refactor if we can							
	//std::vector<unsigned int> local_draw_indicies; 

	std::vector<unsigned int> indicies_of_cell_segments; 
};

class PowerDiagram
{
			//		súly , polar
	typedef std::pair<float, PolarBall*> priority_ball_tuple;

	struct HighestOrderFirst
	{
		bool operator()(const priority_ball_tuple &a, const priority_ball_tuple &b) const { return a.first < b.first; }
	};

	typedef power_crust_utils::priority_queue_for_pole<priority_ball_tuple, std::vector<priority_ball_tuple>, PolarBall*, HighestOrderFirst > priority_ball_queue;
	
public:
	PowerDiagram() : bounded(Box(-1.0f, 1.0f))
	{
		cells.clear();
		point_map.clear();
		pole_map.clear();
		regular_triangles.clear();
		indicies_of_power_diagram_segments.clear();
		inner_outer_related_cell_pairs.clear();
		power_crust_indicies.clear();

	}
	~PowerDiagram() 
	{
		cells.clear();
		point_map.clear();
		pole_map.clear();
		regular_triangles.clear();
		indicies_of_power_diagram_segments.clear();
		inner_outer_related_cell_pairs.clear();
		power_crust_indicies.clear();
	}

	typedef std::vector<PowerCell>::iterator Power_cell_iterator;
	typedef std::vector<PowerCell>::const_iterator const_Power_cell_iterator;

	Power_cell_iterator cells_begin() {	return cells.begin(); }

	Power_cell_iterator cells_end() {	return cells.end(); }

	PowerCell* getCell(const int& ind) { return &cells[ind]; }

	void addCell(const PowerCell& vc) { cells.push_back(vc); }

	void clear() {
		cells.clear();
		point_map.clear();
		pole_map.clear();
		regular_triangles.clear();
		indicies_of_power_diagram_segments.clear();
		inner_outer_related_cell_pairs.clear();
		power_crust_indicies.clear();
	}

	size_t size() const { return cells.size(); }
	size_t size_of_points() const { return point_map.size(); }
	void setBox(const Box& b) { bounded = b; }

	std::vector<Point> getCellPoints(const int& ind) { return cells[ind].getPoints(); }

	std::vector<Point> getAllPoints()
	{
		std::vector<Point> tmp;
		for (auto it : point_map)
			tmp.push_back(it.first);

		return tmp;
	}

	GLPaintFormat getCellPaintData(const int& ind, const int& nid);
	GLPaintFormat getAllPaintData();
	GLPaintFormat getPowerDiagramBySegmentPaintData();
	GLPaintFormat getInnerPolesPaintData();
	GLPaintFormat getOuterPolesPaintData();
	GLPaintFormat getUnknownPolesPaintData();
	GLPaintFormat getInnerOuterPairsPaintData();
	GLPaintFormat getPowerCrustPaintData();
	GLPaintFormat getPowerShapePaintData();

	//------------------------------------------
	//surf pontok a polejaikkal
	void calc_diagram(const std::vector<std::pair<Pole, Point> >& weighted_points);
	void label_poles();
	void calc_medial_axis()
	{
		if (inner_outer_related_cell_pairs.size() == 0) 
			calc_inner_outer_pairs();

		qDebug() << "Medial axis calculation";
		//medial_axis_point_map.clear();
		//unsigned int medial_axis_p_index = 0;

	
		//qDebug() << "We have " << medial_axis_point_map.size() << " medial axis points";
	}

protected:
	Point getBoundedDual(Point& dual);
	void addBoxPoints(const Box& box, std::set<Point>& box_points, std::vector<std::pair<weighted_point, std::vector<Point>>>& points);
	Box UpdateBox(const std::vector<std::pair<Pole, Point> >& weighted_points);

	void setDiagramIndicies(std::set<std::pair<Point*, Point*>> power_diagram_segments)
	{
		std::map<Point*, unsigned int> tmp_vertex_order; //TODO: rework with the private one
		int index = 0;
		for (auto& it : point_map)
		{
			tmp_vertex_order.insert(std::pair<Point*, unsigned int>(const_cast<Point*>(&it.first), index));
			index++;
		}

		for (auto it : power_diagram_segments)
		{
			auto left = tmp_vertex_order[it.first];
			auto right = tmp_vertex_order[it.second];

			indicies_of_power_diagram_segments.push_back(left);
			indicies_of_power_diagram_segments.push_back(right);
		}
	}

	void calc_inner_outer_pairs();
	void calc_power_crust();

private:
	typedef std::pair<unsigned int, std::set<unsigned int>> Vertex_PowerCell_Ids;
						// global vertex index ; related Power Cell index
	std::vector<PowerCell> cells;
	std::map<Point, Vertex_PowerCell_Ids> point_map;	//Power verticies with index and related Power Cell indicies
	std::map<Point, unsigned int> pole_map;				//Pole verticies with index
	
	std::set<triple<unsigned int>> regular_triangles;

	std::vector<unsigned int> indicies_of_power_diagram_segments;
	std::set<std::pair<unsigned int, unsigned int >> inner_outer_related_cell_pairs;
	std::vector<unsigned int> power_crust_indicies;
	//-----------------------
	
	Box bounded;

//	std::map<Point, unsigned int> medial_axis_point_map;
//	std::vector<unsigned int> medial_axis_segment_indicies;

};


/*
TODO:
Occur "Problem":
As we have close points each other their duals also will be close to each other which means sometimes it seems 2 cells
have 1 common points but in real its more than 3 just have 10^5 error between them.

Resolve:
Mybe use some accuracy depends on the density
*/

#pragma once

#include "types.h"
#include "cgal_types.h"
#include "diagram_types.h"
#include "priority_queue.h"

#include <QDebug>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <CGAL/Linear_algebraHd.h>

typedef Kernel::Weighted_point_3                            weighted_point;

typedef CGAL::Regular_triangulation_vertex_base_3<Kernel>        regular_triangulation_vertex_base;
typedef CGAL::Triangulation_vertex_base_with_info_3<std::vector<Point>, Kernel, regular_triangulation_vertex_base> regular_triangulation_vertex_base_with_info;
typedef CGAL::Regular_triangulation_cell_base_3<Kernel>          regular_triangulation_cell_base;
typedef CGAL::Triangulation_cell_base_with_info_3<cell_inf, Kernel, regular_triangulation_cell_base>  regular_triangulation_cell_base_with_info;

typedef CGAL::Triangulation_data_structure_3<
	regular_triangulation_vertex_base_with_info,
	regular_triangulation_cell_base_with_info,
	Concurrency_tag
> regular_triangulation_data_structure;

typedef CGAL::Regular_triangulation_3<Kernel, regular_triangulation_data_structure> regular_triangulation;

typedef regular_triangulation::Cell_handle regular_cell_handle;

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
	PowerCell(Point* pol_p, const float& r, const std::vector<Point>& sur_p, const unsigned int& c_i) : polar_ball(PolarBall(Pole(pol_p, r), sur_p)), my_index(c_i){
		cell_faces.clear();
	}
	~PowerCell() {
		cell_faces.clear();
	}

	GLPaintFormat getPaintData()
	{
		GLPaintFormat res;

		std::map<Point*, unsigned int> points;
		std::map<std::pair<unsigned int, unsigned int>, unsigned int> edge_ix;
		int index = 0;
		int edge_index = 0;
		for (auto& it : cell_faces)
		{
			//qDebug() << it;
			//qDebug() << "Face have : " << it->edges.size();
			for (auto& edge : it->edges)
			{
				//qDebug() << "\t" <<  edge;

				auto first_edge = edge->edge.first;
				auto second_edge = edge->edge.second;

				auto first_poi = points.insert(std::pair<Point*, unsigned int>(first_edge.second, index));
				if (first_poi.second == true)
				{
					res.points.push_back(*first_edge.second);
					index++;
				}

				auto second_poi = points.insert(std::pair<Point*, unsigned int>(second_edge.second, index));
				if (second_poi.second == true)
				{
					res.points.push_back(*second_edge.second);
					index++;
				}

				auto edge_n = edge_ix.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(std::pair<unsigned int, unsigned int>(first_edge.first, second_edge.first), edge_index));
				if (edge_n.second == true)
				{
					res.ix.push_back(first_poi.first->second);
					res.ix.push_back(second_poi.first->second);
					edge_index++;
				}

			}
		}
		qDebug() << "We got everything: " << res.points.size() << " - " << res.ix.size();

		res.centers_with_radius.push_back(std::pair<Point, float>(polar_ball.getPoint(), polar_ball.getRadius()));
		//res.surf_centers.push_back(polar_ball.getSurfPoint());

		res.col.push_back(QVector3D(0, 0, 1));
		res.center_part_lengths.push_back(res.centers_with_radius.size());

		return res;
	}

	Pole getPolePoint() const { return polar_ball.getPole(); }

	bool hasInnerPole() const { return polar_ball.isInnerPole(); }
	PolarBall* getPolarBallPtr() { return &polar_ball; }
	std::vector<Point> getRelatedSurfPoints() const { return polar_ball.getSurfPoints(); }

	void appendPaintDataAsNeighbour(GLPaintFormat& res, const unsigned int& related)
	{
		//	if (related < 0 || neighbour_candidates.size() <= related) return;

		const size_t pre_size = res.points.size();

		std::map<Point*, unsigned int> points;
		std::map<std::pair<unsigned int, unsigned int>, unsigned int> edge_ix;
		int index = 0;
		int edge_index = 0;
		face_t* related_face;

		for (auto& it : cell_faces)
		{
			qDebug() << it;
			qDebug() << "Face have : " << it->edges.size();
			for (auto& edge : it->edges)
			{
				qDebug() << "\t" << edge;

				auto first_edge = edge->edge.first;
				auto second_edge = edge->edge.second;

				auto first_poi = points.insert(std::pair<Point*, unsigned int>(first_edge.second, index));
				if (first_poi.second == true)
				{
					res.points.push_back(*first_edge.second);
					index++;
				}

				auto second_poi = points.insert(std::pair<Point*, unsigned int>(second_edge.second, index));
				if (second_poi.second == true)
				{
					res.points.push_back(*second_edge.second);
					index++;
				}

				auto edge_n = edge_ix.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(std::pair<unsigned int, unsigned int>(first_edge.first, second_edge.first), edge_index));
				if (edge_n.second == true)
				{
					res.ix.push_back(first_poi.first->second + pre_size);
					res.ix.push_back(second_poi.first->second + pre_size);
					edge_index++;
				}

			}

			if (it->get_neighbour_cell(my_index) == related)
				related_face = it;
		}

		qDebug() << "Neighbour vertices: " << "\n";

		auto face_points = related_face->get_face_points();
		for (auto& point : face_points)
		{
			qDebug() << "\t " << qSetRealNumberPrecision(10) << point->x() << ", " << point->y() << ", " << point->z();
			res.centers_with_radius.push_back(std::pair<Point, float>(*point, 0.005f));
		}

		res.col.push_back(QVector3D(0, 0, 0));
		res.center_part_lengths.push_back(face_points.size());
	}

	void add_face(face_t * face)
	{
		cell_faces.push_back(face);
	}

	std::vector<unsigned int> get_neighbours()
	{
		if (cell_faces.empty())	return std::vector<unsigned int>();

		std::vector<unsigned int> my_neighbours;

		for (auto& face : cell_faces)
		{
			int neigh_candidate = face->get_neighbour_cell(my_index);
			if (neigh_candidate == -1 || neigh_candidate == -3) qDebug() << "ERR: BAD getting neighbours";
			if (neigh_candidate < 0) continue;

			my_neighbours.push_back(neigh_candidate);
		}

		return my_neighbours;
	}

	int get_neighbour(const int& face_id)
	{
		if (face_id < 0) return -1;	//important because of drawings

		int correct_face_id = face_id % cell_faces.size();

		return cell_faces[correct_face_id]->get_neighbour_cell(my_index);
	}

	face_t* get_common_face(const unsigned int& neighbour)
	{
		for (auto& face : cell_faces)
		{
			int neigh_candidate = face->get_neighbour_cell(my_index);
			if (neigh_candidate == -1 || neigh_candidate == -3) qDebug() << "ERR: BAD neighbours";
			if (neigh_candidate < 0) continue;

			if (neigh_candidate == neighbour) return face;
		}

		return nullptr;
	}

	std::vector<face_t*> faces() const
	{
		return cell_faces;
	}

private:
	PolarBall polar_ball;

	std::vector<face_t*> cell_faces;
	unsigned int my_index;		//TODO: this is redundant
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
		power_faces.clear();
		power_edges.clear();
		power_vertices.clear();
		pole_map.clear();
		regular_triangles.clear();
	}
	~PowerDiagram() 
	{
		cells.clear();
		power_edges.clear();
		power_faces.clear();
		power_vertices.clear();
		pole_map.clear();
		regular_triangles.clear();
	}

	typedef std::vector<PowerCell>::iterator Power_cell_iterator;
	typedef std::vector<PowerCell>::const_iterator const_Power_cell_iterator;

	Power_cell_iterator cells_begin() {	return cells.begin(); }

	Power_cell_iterator cells_end() {	return cells.end(); }

	PowerCell* getCell(const int& ind) { return &cells[ind]; }

	void addCell(const PowerCell& vc) { cells.push_back(vc); }

	void clear() {
		cells.clear();
		pole_map.clear();
		regular_triangles.clear();
	}

	size_t size() const { return cells.size(); }
	void setBox(const Box& b) { bounded = b; }

	GLPaintFormat getCellPaintData(const int& ind, const int& nid);
	GLPaintFormat getPowerDiagramBySegmentPaintData();
	GLPaintFormat getInnerPolesPaintData();
	GLPaintFormat getOuterPolesPaintData();
	GLPaintFormat getUnknownPolesPaintData();
	GLPaintFormat getInnerOuterPairsPaintData();
	GLPaintFormat getPowerCrustPaintData();
	GLPaintFormat getPowerShapePaintData();
	GLPaintFormat get_triangled_mesh_without_texture();

	//------------------------------------------
	//surf pontok a polejaikkal
	void calc_diagram(const std::vector<std::pair<Pole, Point> >& weighted_points);
	void label_poles();
	void calc_power_crust();

protected:
	typedef std::pair<Point*, unsigned int> point_identifier;
	typedef std::pair<edge_tt*, std::set<unsigned int>> edge_identifier;
	typedef std::pair<int, int> neighbour;

	void addBoxPoints(const Box& box, std::set<Point>& box_points, std::vector<std::pair<weighted_point, std::vector<Point>>>& points);
	Box UpdateBox(const std::vector<std::pair<Pole, Point> >& weighted_points);

	std::pair<edge_tt, edge_identifier> add_edge_as_power_edge(const edge_tt&, std::map<edge_tt, edge_identifier>&);
	void update_neighbour_map(std::map<neighbour, std::set<edge_tt*>>&, edge_identifier&, const int&);

	void set_cell_info(const regular_triangulation& r, regular_cell_handle& cell, unsigned int& pow_index);
	void correct_unoriented_power_crust_faces(std::vector<face_t*>& unoriented_face_ptrs);

private:
	typedef typename _vertex::v_ptr v_ptr;
	typedef typename face_tt::f_ptr f_ptr;

	std::vector<PowerCell> cells;

	std::vector<v_ptr> power_vertices;
	std::set<edge_tt> power_edges;
	std::vector<f_ptr> power_faces;

	std::map<Point, unsigned int> pole_map;				//Pole verticies with index
	
	std::set<triple<unsigned int>> regular_triangles;
	//-----------------------
	
	Box bounded;
};


/*
TODO:
Occur "Problem":
As we have close points each other their duals also will be close to each other which means sometimes it seems 2 cells
have 1 common points but in real its more than 3 just have 10^5 error between them.

Resolve:
Mybe use some accuracy depends on the density
*/

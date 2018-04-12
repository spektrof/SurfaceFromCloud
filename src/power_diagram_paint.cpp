#include "power_diagram.h"

GLPaintFormat PowerDiagram::getPowerCrustPaintData()
{
	GLPaintFormat res;

	std::map<Point*, unsigned int> point_map;
	int p_index = 0;
	std::map<neighbour, unsigned int> edge_map;
	int e_index = 0;

	for (auto& face : power_faces)
	{
		if (!face->face.is_power_crust_face()) continue;

		std::set<Point*> face_points = face->face.get_face_points();

		for (auto& point : face_points)
		{
			auto pos = point_map.insert(std::pair<Point*, unsigned int>(point, p_index));
			if (pos.second == true)
			{
				res.points.push_back(*point);
				p_index++;
			}
		}

		auto face_edges = face->face.edges;

		auto orientation = face_edges[0]->related_faces[&face->face];
		Point* p1 = orientation == 0 ? face_edges[0]->edge.first.second : face_edges[0]->edge.second.second;
		auto first_pos = point_map[p1];

		std::vector<edge_tt*> _face_edges(face_edges.begin() + 1, face_edges.end() - 1);

		for (auto& edge : _face_edges)
		{
			res.ix.push_back(first_pos);

			orientation = edge->related_faces[&face->face];
			p1 = orientation == 0 ? edge->edge.first.second : edge->edge.second.second;
			auto pos = point_map[p1];
			res.ix.push_back(pos);
			res.ix.push_back(pos);

			p1 = orientation == 0 ? edge->edge.second.second : edge->edge.first.second;
			pos = point_map[p1];
			res.ix.push_back(pos);
			res.ix.push_back(pos);

			res.ix.push_back(first_pos);
		}
	}

	res.points_col.push_back(QVector3D(0, 1, 0));
	res.point_part_lengths.push_back(res.ix.size());
	qDebug() << "Power Crust indicies: " << res.ix.size();

	return res;
}

GLPaintFormat PowerDiagram::get_triangled_mesh_without_texture()
{
	GLPaintFormat res;

	std::map<Point*, unsigned int> point_map;	// to get indicies of triangles
	int point_index = 0;

	/*for (auto& face : power_faces)
	{
		if (!face->face.is_power_crust_face()) continue;

		// store the new points of face

		std::set<Point*> face_points = face->face.get_face_points();

		for (auto& point : face_points)
		{
			auto pos = point_map.insert(std::pair<Point*, unsigned int>(point, point_index));
			if (pos.second == true)
			{
				res.points.push_back(*point);
				point_index++;
			}
		}

		// partition the face into triangles
		auto face_edges = face->face.edges;
		auto edge_orientation = face_edges[0]->related_faces[&face->face];

		Point* edge_point = edge_orientation == 0 ? face_edges[0]->edge.first.second : face_edges[0]->edge.second.second;
		auto first_triangle_index = point_map[edge_point];

		face_edges = std::vector<edge_tt*>(face_edges.begin() + 1, face_edges.end() - 1);	// we dont need the first and last edge

		if (face_edges.size() == 0) qDebug() << "ERR: Face contains only 2 edge";

		for (auto& edge : face_edges)
		{
			res.ix.push_back(first_triangle_index);

			edge_orientation = edge->related_faces[&face->face];
			edge_point = edge_orientation == 0 ? edge->edge.first.second : edge->edge.second.second;
			auto pos = point_map[edge_point];
			res.ix.push_back(pos);

			edge_point = edge_orientation == 0 ? edge->edge.second.second : edge->edge.first.second;
			pos = point_map[edge_point];
			res.ix.push_back(pos);
		}
	}
	*/

	face_t* start_face;

	for (auto& face : power_faces)
	{
		if (face->face.is_power_crust_face())
		{
			start_face = &face->face;
			break;
		}
	}

	std::queue<face_t*> face_list;
	std::set<face_t*> done_faces;

	done_faces.insert(start_face);

	face_list.push(start_face);

	while (!face_list.empty())
	{
		face_t* act_face = face_list.front();
		face_list.pop();

		std::set<Point*> face_points = act_face->get_face_points();

		for (auto& point : face_points)
		{
			auto pos = point_map.insert(std::pair<Point*, unsigned int>(point, point_index));
			if (pos.second == true)
			{
				res.points.push_back(*point);
				point_index++;
			}
		}

		auto face_edges = act_face->edges;

		for (auto& face_edge : face_edges)
		{
			face_t* next_boundary_face = nullptr;
			int orientation_on_actual_face = -1;
			int next_face_orientation = -1;
			for (auto& related_face : face_edge->related_faces)
			{
				if (!related_face.first->is_power_crust_face()) continue;
				if (related_face.first == act_face)
				{
					orientation_on_actual_face = related_face.second;
					continue;
				}

				if (next_boundary_face == nullptr)
				{
					next_boundary_face = related_face.first;
					next_face_orientation = related_face.second;
				}
				else
					qDebug() << "CORRECT ORIENTATION ON BOUND IS BAAD";

			}

			if (next_boundary_face != nullptr)
			{
				if (done_faces.find(next_boundary_face) != done_faces.end())
				{
					if (next_face_orientation == orientation_on_actual_face)
					{
						qDebug() << "ERR: Bad orientation!";
					}
					continue;
				}

				if (next_face_orientation < 0 || orientation_on_actual_face < 0)
				{
					qDebug() << "ERR: orientation didnt set properly";
					continue;
				}

				face_list.push(next_boundary_face);
				done_faces.insert(next_boundary_face);
			}
		}

		// partition the face into triangles
		auto edge_orientation = face_edges[0]->related_faces[act_face];

		Point* edge_point = edge_orientation == 0 ? face_edges[0]->edge.first.second : face_edges[0]->edge.second.second;
		auto first_triangle_index = point_map[edge_point];

		face_edges = std::vector<edge_tt*>(face_edges.begin() + 1, face_edges.end() - 1);	// we dont need the first and last edge

		if (face_edges.size() == 0) qDebug() << "ERR: Face contains only 2 edge";

		for (auto& edge : face_edges)
		{
			res.ix.push_back(first_triangle_index);

			edge_orientation = edge->related_faces[act_face];
			edge_point = edge_orientation == 0 ? edge->edge.first.second : edge->edge.second.second;
			auto pos = point_map[edge_point];
			res.ix.push_back(pos);

			edge_point = edge_orientation == 0 ? edge->edge.second.second : edge->edge.first.second;
			pos = point_map[edge_point];
			res.ix.push_back(pos);
		}
	}

	qDebug() << "Power Crust indicies: " << res.ix.size();

	return res;
}

GLPaintFormat PowerDiagram::getPowerShapePaintData()
{
	GLPaintFormat res;

	std::map<unsigned int, unsigned int> global_local_ind_of_poles;
	int index = 0;
	for (auto& it : pole_map)
	{
		res.points.push_back(it.first);
		global_local_ind_of_poles.insert(std::pair<unsigned int, unsigned int>(it.second, index));
		index++;
	}

	for (auto& _triple : regular_triangles)
	{
		int ct = cells[_triple.first].hasInnerPole() + cells[_triple.second].hasInnerPole() + cells[_triple.third].hasInnerPole();
		if (ct < 3) continue;
		//if (!cells[_triple.first].hasInnerPole() || !cells[_triple.second].hasInnerPole() || !cells[_triple.third].hasInnerPole() ) continue;

		int tmp = global_local_ind_of_poles.size();

		res.ix.push_back(global_local_ind_of_poles[_triple.first]);
		res.ix.push_back(global_local_ind_of_poles[_triple.second]);
		res.ix.push_back(global_local_ind_of_poles[_triple.second]);
		res.ix.push_back(global_local_ind_of_poles[_triple.third]);
		res.ix.push_back(global_local_ind_of_poles[_triple.third]);
		res.ix.push_back(global_local_ind_of_poles[_triple.first]);

		if (tmp != global_local_ind_of_poles.size()) qDebug() << "BAAAD";
	}

	res.point_part_lengths.push_back(res.ix.size());
	res.points_col.push_back(QVector3D(0, 1, 0));

	return res;
}

GLPaintFormat PowerDiagram::getCellPaintData(const int& ind, const int& nid)
{
	if (size() == 0) return GLPaintFormat();

	GLPaintFormat res = cells[ind].getPaintData();
	res.point_part_lengths.push_back(res.ix.size());

	int n = cells[ind].get_neighbour(nid);
	if (n >= 0)
	{
		cells[n].appendPaintDataAsNeighbour(res, ind);
		res.point_part_lengths.push_back(res.ix.size() - res.point_part_lengths[0]);
		res.points_col.push_back(QVector3D(0, 1, 0));
		res.points_col.push_back(QVector3D(1, 0, 0));
	}
	else
		res.points_col.push_back(QVector3D(1, 0, 0));

	return res;
}

GLPaintFormat PowerDiagram::getInnerPolesPaintData()
{
	GLPaintFormat res;

	for (auto& it : cells)
	{
		PolarBall* tmp = it.getPolarBallPtr();

		if (tmp->isInnerPole())
			res.centers_with_radius.push_back(std::pair<Point, float>(tmp->getPoint(), tmp->getRadius()));
	}

	res.col.push_back(QVector3D(0, 1, 0));
	res.center_part_lengths.push_back(res.centers_with_radius.size());

	return res;
}

GLPaintFormat PowerDiagram::getOuterPolesPaintData()
{
	GLPaintFormat res;

	for (auto& it : cells)
	{
		PolarBall* tmp = it.getPolarBallPtr();
		if (tmp->isOuterPole())
			res.centers_with_radius.push_back(std::pair<Point, float>(tmp->getPoint(), tmp->getRadius()));
	}
	res.col.push_back(QVector3D(1, 0, 0));
	res.center_part_lengths.push_back(res.centers_with_radius.size());

	return res;
}

GLPaintFormat PowerDiagram::getUnknownPolesPaintData()
{
	GLPaintFormat res;

	for (auto& it : cells)
	{
		PolarBall* tmp = it.getPolarBallPtr();
		if (tmp->isUnkown())
			res.centers_with_radius.push_back(std::pair<Point, float>(tmp->getPoint(), tmp->getRadius()));
	}

	res.col.push_back(QVector3D(0, 0, 0));
	res.center_part_lengths.push_back(res.centers_with_radius.size());

	return res;
}

GLPaintFormat PowerDiagram::getPowerDiagramBySegmentPaintData()
{
	GLPaintFormat res;
	Point max = Point(0, 0, 0);
	for (auto it : power_vertices)
	{
		res.points.push_back(it->point);
		Vector a = it->point - Point(0, 0, 0);
		if (a.squared_length() > (max - Point(0, 0, 0)).squared_length()) max = it->point;
	}
	qDebug() << max.x() << ", " << max.y() << ", " << max.z();

	for (auto it : power_edges)
	{
		res.ix.push_back(it.edge.first.first);
		res.ix.push_back(it.edge.second.first);
	}

	//	res.ix = indicies_of_power_diagram_segments;
	qDebug() << "points " << res.points.size();
	qDebug() << "indicies_of_power_diagram_segments " << res.ix.size();
	res.col.push_back(QVector3D(1, 0, 0));
	return res;
}

GLPaintFormat PowerDiagram::getInnerOuterPairsPaintData()
{
	GLPaintFormat res;

	//TODO: remove edge mltiplications
	std::map<Point, unsigned int> io_points;
	unsigned int point_index = 0;

	for (auto& face : power_faces)
	{
		if (!face->face.is_power_crust_face()) continue;

		unsigned int cell_one = face->face.related_cells[0];
		unsigned int cell_two = face->face.related_cells[1];

		unsigned int inner_cell = cells[cell_one].hasInnerPole() ? cell_one : cell_two;

		auto cell_faces = cells[inner_cell].faces();
		for (auto& cell_face : cell_faces)
		{
			auto face_points = cell_face->get_face_points();

			for (auto& point : face_points)
			{
				auto is_in = io_points.insert(std::pair<Point, unsigned int>(*point, point_index));
				if (is_in.second == true)
				{
					res.points.push_back(*point);
					point_index++;
				}
			}

			auto face_edges = cell_face->edges;
			for (auto& edge : face_edges)
			{
				Point* p1 = edge->edge.first.second;
				Point* p2 = edge->edge.second.second;

				auto point_it = io_points.find(*p1);
				res.ix.push_back(point_it->second);
				point_it = io_points.find(*p2);
				res.ix.push_back(point_it->second);
			}
		}
	}

	res.point_part_lengths.push_back(res.ix.size());
	qDebug() << "Inner cell indexes: " << res.point_part_lengths[0];


	res.points_col.push_back(QVector3D(0, 1, 0));

	return res;
}

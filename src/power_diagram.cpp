#include "power_diagram.h"

/* ****************************
	Diagram calculation main functions
**************************** */

void PowerDiagram::calc_diagram(const std::vector<std::pair<Pole, Point> >& weighted_points)
{
	qDebug() << "Power Diagram calculation";
	/*************************************
	*	Clear	*
	**************************************/

	cells.clear();
	power_vertices.clear();
	power_edges.clear();
	power_faces.clear();
	pole_map.clear();
	regular_triangles.clear();

	/*************************************
	*	Box update & prepare for Weighted Points	*
	**************************************/

	qDebug() << "\tOUR BOX IS: " << bounded.getXMin() << " " << bounded.getXMax() << " " << bounded.getYMin() << " " << bounded.getYMax() << " " << bounded.getZMin() << " " << bounded.getZMax();

	qDebug() << "\tInput: " << weighted_points.size() << " poles";

	std::map<Point, std::pair<float, std::vector<Point>>> pole_point_longest_radius_pairs;

	auto start = hrclock::now();
	auto process = hrclock::now();

	//max radius
	/*for (auto& it : weighted_points)
	{
		std::pair<std::map<Point, std::pair<float, std::vector<Point>>>::iterator, bool> pair_it;
		pair_it = pole_point_longest_radius_pairs.insert(std::pair<Point, std::pair<float, std::vector<Point>> >(*it.first.center, std::pair<float, std::vector<Point>>(it.first.radius, std::vector<Point>())));
		if (pair_it.second == false && pair_it.first->second.first < it.first.radius)
			pair_it.first->second.first = it.first.radius;

		pair_it.first->second.second.push_back(it.second);
	}*/

	//min radius
	/*for (auto& it : weighted_points)
	{
		std::pair<std::map<Point, std::pair<float, std::vector<Point>>>::iterator, bool> pair_it;
		pair_it = pole_point_longest_radius_pairs.insert(std::pair<Point, std::pair<float, std::vector<Point>> >(*it.first.center, std::pair<float, std::vector<Point>>(it.first.radius, std::vector<Point>())));
		if (pair_it.second == false && pair_it.first->second.first > it.first.radius)
			pair_it.first->second.first = it.first.radius;

		pair_it.first->second.second.push_back(it.second);
	}
	*/

	//atlag
	std::map<Point, std::vector<float>> tmp_rads;
	for (auto& it : weighted_points)
	{
		std::pair<std::map<Point, std::vector<float>>::iterator, bool> rads;
		rads = tmp_rads.insert(std::pair<Point, std::vector<float> >(*it.first.center, std::vector<float>()));
		rads.first->second.push_back(it.first.radius);

		std::pair<std::map<Point, std::pair<float, std::vector<Point>>>::iterator, bool> pair_it;
		pair_it = pole_point_longest_radius_pairs.insert(std::pair<Point, std::pair<float, std::vector<Point>> >(*it.first.center, std::pair<float, std::vector<Point>>(it.first.radius, std::vector<Point>())));
		if (pair_it.second == false)
		{
			float sum = 0.0f;
			for (auto s : rads.first->second) sum += s;

			pair_it.first->second.first = sum / (float)rads.first->second.size();
		}

		pair_it.first->second.second.push_back(it.second);
	}
	tmp_rads.clear();

	qDebug() << "\t\texact: " << pole_point_longest_radius_pairs.size() << " different poles!";
	qDebug() << "\t\t\tResolve redundant poles under : " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - process).count();

	/*************************************
	*	Checks	*
	**************************************/

	//Checking surf points for poles (added more than once)
	int ct = 0;
	for (auto it : pole_point_longest_radius_pairs)
	{
		std::set<Point> tmp;
		for (auto po_it : it.second.second)
			tmp.insert(po_it);

		ct += it.second.second.size();
		if (tmp.size() != it.second.second.size())
			qDebug() << "\tERR: Surf point added more than once!";
	}

	if (ct != weighted_points.size())
		qDebug() << "\tERR: We lost some surf point pole pair";

	/*************************************
	*	Make weighted points and do Triangulation	*
	**************************************/
	qDebug() << "\t\tMaking weighted points";

	std::vector<std::pair<weighted_point, std::vector<Point>>> w_p;
	for (auto it : pole_point_longest_radius_pairs)
		w_p.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(it.first, it.second.first* it.second.first), it.second.second));	//d*d - r -> az alap sugár négyzete kell

	pole_point_longest_radius_pairs.clear();

	qDebug() << "\t\tBox stuff";

	std::set<Point> box_points;
	Box pole_bound = UpdateBox(weighted_points);
	addBoxPoints(pole_bound, box_points, w_p);

	qDebug() << "\tPreparation ended in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - start).count();
	qDebug() << "\tTriangulation";
	process = hrclock::now();

	regular_triangulation::Lock_data_structure locks(CGAL::Bbox_3(pole_bound.getXMin(), pole_bound.getYMin(), pole_bound.getZMin(), pole_bound.getXMax(), pole_bound.getYMax(), pole_bound.getZMax()), 50);
	regular_triangulation R(w_p.begin(), w_p.end(), &locks);
	w_p.clear();

	qDebug() << "\t\t" << R.number_of_cells() << " cells " << R.number_of_finite_cells() << " finite cells";
	qDebug() << "\t\t" << R.number_of_vertices() << " vertices ";

	assert(R.is_valid());

	/*************************************
	*	Calculate Power Cellas	*
	**************************************/
	qDebug() << "\tTriangulation ended in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - process).count();

	process = hrclock::now();

	unsigned int pow_index = 0;
	unsigned int pol_index = 0;

	std::map<edge_tt, edge_identifier> edge_cell_ids_map;
	std::map<neighbour, std::set<edge_tt*>> neighbour_map;

	std::vector<regular_triangulation::Finite_vertices_iterator> box_v_its;

	for (regular_triangulation::Finite_vertices_iterator vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); vit++)
	{
		//skip box points
		if (box_points.count(vit->point().point()) == 1)
		{
			box_v_its.push_back(vit);
			continue;
		}

		//-----------------------------------------------
		//get incident cells
		std::vector<regular_cell_handle> inc_cells;
		R.incident_cells(vit, std::back_inserter(inc_cells));
	
		//check for empty
		if (inc_cells.empty())
			continue;

		//-----------------------------------------------
		//Precalculations
		auto inc_c_size = inc_cells.size();

		std::set<regular_cell_handle> inc_cell_map;	//for finding cell in logM

		for (auto& it : inc_cells)
			inc_cell_map.insert(it);

		//-----------------------------------------------
		//store pole point for reference

		std::pair<std::map<Point, unsigned int>::iterator, bool> pol_res;
		pol_res = pole_map.insert(std::pair<Point, unsigned int >(vit->point().point(), pol_index));
		if (pol_res.second == true)
			pol_index++;	//not in the map
		else
			qDebug() << "ERR: BAAAD, double POLE";

		//-----------------------------------------------
		//creating powercell								pole ref				get radius			related surf points
		PowerCell pc = PowerCell(const_cast<Point*>(&(pol_res.first->first)), sqrt(vit->point().weight()), vit->info(), cells.size());		//vissza konverzió
		unsigned int act_pc_index = cells.size();

		std::set<edge_tt*> related_edges_of_actual_inc_cells;		//for neighbours calculation
		
		//-----------------------------------------------

		for (size_t cit = 0; cit < inc_c_size; ++cit)
		{
			regular_cell_handle act_cell = inc_cells[cit];
			
			set_cell_info(R, act_cell, pow_index);
			unsigned int act_index = act_cell->info().index;
			Point* act_vert_ptr = act_cell->info().dual;

			//-----------------------------------------------
			//calculate edges

			for (int i = 0; i < 4; ++i)
			{
				regular_cell_handle tmp_cell = act_cell->neighbor(i);
				auto other_cell_candidate = inc_cell_map.find(tmp_cell);

				if (other_cell_candidate != inc_cell_map.end())	//both cell_handle are related to the actual point
				{
					//calculate the dual if its not ready yet
					set_cell_info(R, tmp_cell, pow_index);
					unsigned int tmp_index = tmp_cell->info().index;
					Point* tmp_vert_ptr = tmp_cell->info().dual;

					//calculate the edge
					edge_tt edge;

					if (act_index < tmp_index)
						edge = edge_tt(edge_t(edge_point(act_index, act_vert_ptr), edge_point(tmp_index, tmp_vert_ptr)));
					else
						edge = edge_tt(edge_t(edge_point(tmp_index, tmp_vert_ptr), edge_point(act_index, act_vert_ptr)));

					//put the edge into the map
					std::pair<edge_tt, edge_identifier> edge_element_from_map = add_edge_as_power_edge(edge, edge_cell_ids_map);
					
					//save the edges for the actual cell
					related_edges_of_actual_inc_cells.insert(edge_element_from_map.second.first);

					//create relations between the cells based on the edges
					update_neighbour_map(neighbour_map, edge_element_from_map.second, act_pc_index);
				}
			}
			
		}

		//add the new cell index to the edges list
		for (auto& it : related_edges_of_actual_inc_cells)
		{
			edge_cell_ids_map[*it].second.insert(act_pc_index);
		}

		//---------------------------

		addCell(pc);
	}

	//to get those face which belongs to only 1 cell!
	int box_ind = 0;
	for (auto& it : box_v_its)
	{
		box_ind--;
		//get incident cells
		std::vector<regular_cell_handle> inc_cells;
		R.incident_cells(it, std::back_inserter(inc_cells));
	
		//-----------------------------------------------
		//Precalculations
		auto inc_c_size = inc_cells.size();

		std::set<regular_cell_handle> inc_cell_map;	//for finding cell in logM

		for (auto& it : inc_cells)
			inc_cell_map.insert(it);

		unsigned int act_pc_index = box_ind;

		for (size_t cit = 0; cit < inc_c_size; ++cit)
		{
			regular_cell_handle act_cell = inc_cells[cit];
			int act_index = act_cell->info().index;

			if (act_index == -1) continue;

			Point* act_dual = act_cell->info().dual;

			//-----------------------------------------------
			//calculate edges
			for (int i = 0; i < 4; ++i)
			{
				regular_cell_handle tmp_cell = act_cell->neighbor(i);
				auto other_cell_candidate = inc_cell_map.find(tmp_cell);

				if (other_cell_candidate != inc_cell_map.end())	//both cell_handle are related to the actual point
				{
					int tmp_index = tmp_cell->info().index;

					if (tmp_index == -1) continue;

					Point* tmp_dual = tmp_cell->info().dual;

					//calculate the edge
					edge_tt edge;

					if (act_index < tmp_index)
						edge = edge_tt(edge_t(edge_point(act_index, act_dual), edge_point(tmp_index, tmp_dual)));
					else
						edge = edge_tt(edge_t(edge_point(tmp_index, tmp_dual), edge_point(act_index, act_dual)));
					
					auto edg_p = edge_cell_ids_map.find(edge);
					if (edg_p == edge_cell_ids_map.end()) continue;

					//create relations between the cells based on the edges
					update_neighbour_map(neighbour_map, edg_p->second, act_pc_index);
				}
			}

		}
	}

	int tt = 0;
	for (auto& face_candidate : neighbour_map)
	{
		if (face_candidate.second.size() < 3) continue;

		std::vector<unsigned int> related_cells;
		if (face_candidate.first.first >= 0) related_cells.push_back(face_candidate.first.first);
		if (face_candidate.first.second >= 0) related_cells.push_back(face_candidate.first.second);

		face_tt face;
		power_faces.push_back(std::make_shared<face_tt>(face));
		power_faces.back()->face.add_edges(face_candidate.second, 0);
		power_faces.back()->face.set_related_cell(related_cells);

		//qDebug() << "Neigh test";
		for (auto& ind : related_cells)
		{
			cells[ind].add_face(const_cast<face_t*>(&(power_faces.back()->face)));
		//	qDebug() << "\t" << ind << " my neigh is : " << power_faces.back()->face.get_neighbour_cell(ind);
		}
	}
	
	/*qDebug() << i << " vs " << power_faces.size();
	for (auto& it : power_faces)
	{
		qDebug() << &it->face << " face have " << it->face.edges.size() << " edge";
	}*/

	qDebug() << pole_map.size() << " poles";
	qDebug() << "Power Diagram base is created";
	qDebug() << "\t" << cells.size() << " cella";
	qDebug() << "\t" << power_vertices.size() << " poi";
	qDebug() << "\t" << power_edges.size() << " edge";
	qDebug() << "\t" << power_faces.size() << " faces";
	/*
	Power Diagram base is created
	 136359  cella
	 1025970  poi
	 2051934  edge
	 1162324  faces
	 2GB
	*/
	qDebug() << " *********************************";

	/*
	First implementation was over 10 seconds with voronoi diag poles of 24k points
	
	*/

	qDebug() << "Power diagram is calculated under ended in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - process).count();
	qDebug() << "calc_diagram ended in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - start).count();
	//	9 sec on 21k poi
	/*for (auto it : neighbour_map)
	{
		qDebug() << it.first.first << " - " << it.first.second << " : " ;
		for (auto& it2 : it.second)
			qDebug() << it2;
	}*/

	//qDebug() << "Faces: " << power_faces.size();
	/*for (auto & it : power_faces)
	{
		qDebug() <<"\t"<< it.face.edges.size();
	}*/
	/*qDebug() << "Edges: " << power_edges.size();

	for (auto & it : power_edges)
	{
		qDebug() << it.edge.first.first << " - " << it.edge.second.first;
		qDebug() << "\t" << it.related_faces.size();
	}*/
	//calculate those neighbours which have common face
	/*for (auto& it : cells)
	{
		it.calcProperNeighbours();
	}
	
	*/
	//checking önmagának szomszédja e
	int index = 0;
	for (auto it : cells)
	{
		auto ne = it.get_neighbours();
		std::set<unsigned int> tmp = std::set<unsigned int>(ne.begin(), ne.end());
		if (tmp.find(index) != tmp.end())
		{
			qDebug() << " Own neighbour";
		}

		index++;
	}

	qDebug() << "Power Diagram cells's neighbours are calculated";

	/*int ct2 = 0;
	for (auto& it : weighted_points)
	{
		if (pole_map.find(*it.first.center) == pole_map.end())
		{
			qDebug() << "LOSTED POLE: " << it.first.center->x() << ", " << it.first.center->y() << ", " << it.first.center->z();
			for (auto& it2 : w_p)
			{
				if (it2.first.point() == *it.first.center)
				{
					qDebug() << it2.first.weight() << " weight " << it2.second.size() << " surf_p";
					if (it2.second.size() != 0) qDebug() << it2.second[0].x() <<", " << it2.second[0].y() << ", " << it2.second[0].z() << " surf_p";
				}
			}
		}
	}
	*/

	for (regular_triangulation::Finite_cells_iterator cit = R.finite_cells_begin(); cit != R.finite_cells_end(); cit++)
	{
		Point p1 = cit->vertex(0)->point();
		Point p2 = cit->vertex(1)->point();
		Point p3 = cit->vertex(2)->point();
		Point p4 = cit->vertex(3)->point();

		if (pole_map.find(p1) == pole_map.end() || pole_map.find(p2) == pole_map.end() || pole_map.find(p3) == pole_map.end() || pole_map.find(p4) == pole_map.end()) continue;	//box point

		int presize = pole_map.size();
		
		unsigned int i1 = pole_map[p1];
		unsigned int i2 = pole_map[p2];
		unsigned int i3 = pole_map[p3];
		unsigned int i4 = pole_map[p4];

		if (presize != pole_map.size()) qDebug() << "Bad index calc for poles";

		regular_triangles.insert(triple<unsigned int>(i1, i2, i3));
		regular_triangles.insert(triple<unsigned int>(i1, i2, i4));
		regular_triangles.insert(triple<unsigned int>(i2, i3, i4));
		regular_triangles.insert(triple<unsigned int>(i3, i1, i4));
	}
	
}

void PowerDiagram::set_cell_info(const regular_triangulation& R, regular_cell_handle& cell, unsigned int& pow_index)
{
	auto info = &cell->info();

	if (info->index == -1)
	{
		info->index = pow_index;
		Point dual = R.dual(cell);
		typename _vertex::v_ptr node = std::make_shared<_vertex>(dual);
		power_vertices.push_back(node);
		info->dual = &power_vertices.back()->point;
		pow_index++;
	}
}

std::pair<edge_tt, PowerDiagram::edge_identifier> PowerDiagram::add_edge_as_power_edge(const edge_tt& e, std::map<edge_tt, edge_identifier>& edge_map)
{
	auto edge_it_ref = edge_map.insert(std::pair<edge_tt, edge_identifier>(e, edge_identifier()));
	if (edge_it_ref.second == true)
	{
		auto orig_edge_ref = power_edges.insert(edge_it_ref.first->first);
		edge_it_ref.first->second.first = const_cast<edge_tt*>(&(*orig_edge_ref.first));
	}
	return *edge_it_ref.first;
}

void PowerDiagram::update_neighbour_map(std::map<neighbour, std::set<edge_tt*>>& neighbour_map, edge_identifier& e_i, const int& powercell_index)
{
	for (auto& related_power_cell : e_i.second)
	{
		neighbour act_n;
		if (powercell_index < related_power_cell) act_n = neighbour(powercell_index, related_power_cell);
		else act_n = neighbour(related_power_cell, powercell_index);

		auto res = neighbour_map.insert(std::pair<neighbour, std::set<edge_tt*>>(act_n, std::set<edge_tt*>()));
		res.first->second.insert(e_i.first);
	}
}

/* ****************************
	Labeling functions
**************************** */

float PolarBall::getAlphaWeight(PolarBall* rhs)
{
	float R = this->pole.radius;
	float r = rhs->getRadius();
	float d = (this->getPoint() - rhs->getPoint()).squared_length();

	/**************************
	Law of cosinus for calculate intersection deepness
	*************************************************/

	float lhs_cosB = (r * r + R * R - d) / (2 * r*R);
	float B_angle = acosf(lhs_cosB);

	float A_angle = 3.14159265359f - B_angle;
	float result_attempt = -1 * cosf(A_angle);

	if (result_attempt < 0)		//only 0.0f - 1.0f range is valid, here the intersection was more deeply -> they have more the same label
	{
		//		qDebug() << "\t\tERR: ALPHA val:";
		/*		qDebug() << "\t\t\tR: " << R << ", r: " <<r <<", d: "<< d;
		qDebug() << "\t\t\tB_angle: " << B_angle;
		qDebug() << "\t\t\tA_angle: " << A_angle;
		qDebug() << "\t\t\tlhs_cosB " << lhs_cosB;
		qDebug() << "\t\t\tgetAlphaWeight return: " << -1 * cosf(A_angle);*/
		return 0.0f;
	}

	return result_attempt;
}

float PolarBall::getBetaWeight(PolarBall* rhs, Point& surfpoint)
{
	/**************************
	Check for common surf point
	**********************/

	std::vector<Point> lhs_surface_points = getSurfPoints();
	std::vector<Point> rhs_surface_points = rhs->getSurfPoints();

	std::set<Point> lhs_surface_points_s(lhs_surface_points.begin(), lhs_surface_points.end());
	std::set<Point> rhs_surface_points_s(rhs_surface_points.begin(), rhs_surface_points.end());

	if (lhs_surface_points_s.find(surfpoint) == lhs_surface_points_s.end() ||
		rhs_surface_points_s.find(surfpoint) == rhs_surface_points_s.end())
	{
		qDebug() << "\t\tERR: There is no surf point for beta calculation!";
		return 0.0f;
	}

	/**************************
	Cos of Angle between the 2 vector : from common points towards the center of pole
	*************************/
	Vector vec_one = this->getPoint() - surfpoint;
	Vector vec_two = rhs->getPoint() - surfpoint;

	float lhs = CGAL::scalar_product(vec_one, vec_two);
	float cos_val = lhs / sqrt(vec_one.squared_length() * vec_two.squared_length());

	float result_attempt = -1 * cos_val;

	if (result_attempt < 0)		// this should never happen
	{
		qDebug() << "\t\tERR: BETA val:";
		/*	qDebug() << "\t\t\tvec1 " << vec_one.x() << ", " << vec_one.y() << ", " << vec_one.z();
		qDebug() << "\t\t\tvec2 " << vec_two.x() << ", " << vec_two.y() << ", " << vec_two.z();
		qDebug() << "\t\t\tscalar dot " << lhs;
		qDebug() << "\t\t\tBetaWeight return: " << -1 * cos_val;*/
		return 0.0f;
	}

	return result_attempt;
}

void PowerDiagram::label_poles()
{
	qDebug() << "\nLabeling method";
	qDebug() << "\tMINIMAL BOUND BOX IS: " << bounded.getXMin() << " " << bounded.getXMax() << " " << bounded.getYMin() << " " << bounded.getYMax() << " " << bounded.getZMin() << " " << bounded.getZMax();
	qDebug() << "\n\tPrepare steps";
	auto start = hrclock::now();

	/*************************************
	*	Clear	*
	**************************************/
	priority_ball_queue pri_ball_q;
	std::map<PolarBall*, std::set<PolarBall*> > neigh_graph;
	neigh_graph.clear();

	bool error_during_init = false;
	/*************************************
	*	0. step:
		get Pole pair by surface point	*
	**************************************/
	std::map<Point, std::vector<PolarBall*>> neigh_by_surf;
	neigh_by_surf.clear();

	int ct = 0;
	for (auto& it : cells)
	{
		PolarBall* actual_polarball = it.getPolarBallPtr();
		std::vector<Point> related_surf_points = actual_polarball->getSurfPoints();
		for (auto surf_point : related_surf_points)
		{
			std::pair<std::map<Point, std::vector<PolarBall*>>::iterator, bool> res;
			res = neigh_by_surf.insert(std::pair<Point, std::vector<PolarBall*>>(surf_point, std::vector<PolarBall*>()));
			res.first->second.push_back(actual_polarball);
		}
		ct += related_surf_points.size();
	}
	qDebug() << "All surf_points: " << ct;
	/*************************************
	*	Checks:
			kül. polarball oknak a surf pointjai különbözõek!!!
			Nem 2 polarball tartozik 1 surf pointhoz *
	**************************************/
 
	for (auto& it : cells)
	{
		PolarBall* actual_polarball = it.getPolarBallPtr();
		std::vector<Point> related_surf_points = actual_polarball->getSurfPoints();
		std::set<Point> tmp_rel_surf(related_surf_points.begin(), related_surf_points.end());

		if (related_surf_points.size() != tmp_rel_surf.size())
		{
			error_during_init = true;
			qDebug() << "\t\tERR: One/More surf point added more than once to pole";
		}
	}

	for (auto it : neigh_by_surf)
	{
		//qDebug() << it.first.x() << ", " << it.first.y() << ", " << it.first.z() << " - " << it.second.size() << " " << it.second;
		if (it.second.size() != 2)
		{
		//	error_during_init = true;
			qDebug() << "\t\tERR: We lost a pole! - there is no pair " << it.second.size();
		}
	}

	qDebug() << "\n\tInit step";
	/*************************************
	*	First step:
			Creathe the neighbour graph	*
	**************************************/

	for (auto& cell : cells)
	{
		PolarBall* cell_ptr = cell.getPolarBallPtr();
		std::set<PolarBall*> actual_neighs_ptr;

		//adding neighbours by common face
		auto actual_neighs = cell.get_neighbours();
		for (auto& neigh : actual_neighs)
			actual_neighs_ptr.insert(cells[neigh].getPolarBallPtr());

		//adding neighbours by common surf point
		std::vector<Point> related_surf_points = cell.getRelatedSurfPoints();
		for (auto surf_point : related_surf_points)
		{
			auto pole_pairs = neigh_by_surf[surf_point];	//fix size ( 2 )

			if (pole_pairs[0] != cell_ptr && pole_pairs.size() == 2 && pole_pairs[1] != cell_ptr)
			{
				error_during_init = true;
				qDebug() << "\t\tERR: BAD NEIGHBOUR CALCULATION";
			}

			//surf point base edge
			if (pole_pairs[0] != cell_ptr)
				actual_neighs_ptr.insert(pole_pairs[0]);
			else if (pole_pairs.size() == 2)
				actual_neighs_ptr.insert(pole_pairs[1]);
		}

		neigh_graph.insert(std::pair<PolarBall*, std::set<PolarBall*>>(cell_ptr, actual_neighs_ptr));
	}

	qDebug() << "\t\tNeighbour graph created";
	/*************************************
	*	Second step:
				init the priority queue	*
	**************************************/

	for (auto &it : neigh_graph)
	{
		it.first->setInValue(0);
		it.first->setOutValue(0);
		pri_ball_q.push(priority_ball_tuple(0.0f, it.first));
	}
	qDebug() << "\t\tPriority queue initalized with " << pri_ball_q.size() << " element";

	/*************************************
	*	Third step:
	Setting outer poles	based on the bound box*
	**************************************/
	qDebug() << "\n\tSetting outlier poles";

	priority_ball_queue tmp_q = pri_ball_q;

	int out_p = 0;
	while (!tmp_q.empty())
	{
		priority_ball_tuple tmp = tmp_q.top();
		tmp_q.pop();

		Point pole_point = tmp.second->getPoint();
		if (!bounded.isInside(pole_point.x(), pole_point.y(), pole_point.z()))
		{
			auto element_in_priq = pri_ball_q.find(tmp.second);
			if (pri_ball_q.count(tmp.second) > 1)
			{
				error_during_init = true;
				qDebug() << "\t\tERR: Invalid Priority queue";
			}

			element_in_priq->second->setOutValue(1.0f);
			pri_ball_q.update(element_in_priq);
			out_p++;

			//	qDebug() << "\tFind " << tmp.second << " count " << coun;
			//	pri_ball_q.push(tmp);
			//	element_in_priq = pri_ball_q.find(tmp.second);
			//	coun = pri_ball_q.count(tmp.second);
			//	qDebug() << "\tUPDATE " << element_in_priq->second << " count " << coun;;
		}

	}
	qDebug() << "\t\tOuter pole count (init): " << out_p;

	qDebug() << "\n\\/\\/\\/\\/\\/\\/\\/\\/\\/";

	if (error_during_init)
	{
		qDebug() << "\tERROR: Resolve the errors before go for labeling";
		return;
	}

	/*************************************
	*	Forth step:
			Labeling	*
	**************************************/
	qDebug() << "\n\tStart labeling";

	while (!pri_ball_q.empty())
	{
		/*************************************
		*	Pop high priority element	: */

		auto p = pri_ball_q.top();
		//qDebug() << "\t\tPopped: " << p.second << " with prio : " << p.first;

		if (pri_ball_q.find_higher_prior(p.first))
			qDebug() << "\tPriority queue FAIL ";

		pri_ball_q.pop();

		if (pri_ball_q.count(p.second) > 0)
			qDebug() << "\tMore than one occur ";


		float tmp_p;
		float in_val = p.second->getInValue();
		float out_val = p.second->getOutValue();

		/*************************************
		*	Set its flag	: */

		if (in_val > out_val)
		{
			p.second->setFlag(0);
			tmp_p = in_val;
		//	qDebug() << "\t\t\tIts an Inner Pole";
		}
		else
		{
			p.second->setFlag(1);
			tmp_p = out_val;
		//	qDebug() << "\t\t\tIts an Outer Pole";
		}

		/*************************************
		*	Calculate new weights with Beta weight	: */

		std::vector<Point> tmp_surf_points = p.second->getSurfPoints();
		//qDebug() << "\t\tNumber of related Surf p: " << tmp_surf_points.size();

		for (std::vector<Point>::iterator related_surf_p = tmp_surf_points.begin();
										  related_surf_p != tmp_surf_points.end();		++related_surf_p)
		{
			std::vector<PolarBall*> tmp_related_polarballs = neigh_by_surf[*related_surf_p];

			if (tmp_related_polarballs.size() != 2)
			{
				qDebug() << "\tERR: Wrong related polarball ptrs";
				if (tmp_related_polarballs[0] != p.second) qDebug() << "\tERR: STRANGE - pole count thing to surface point";
				continue;
			}
			PolarBall* opposite_pole = tmp_related_polarballs[0] != p.second ? tmp_related_polarballs[0] : tmp_related_polarballs[1];
			// We got the opposite pole based on the current surface_point

			//	qDebug() << "\t\tFINDING: " << opposite_pole;
			auto opp_position = pri_ball_q.find(opposite_pole);
			if (opp_position != pri_ball_q.end()) { /*qDebug() << "\t\t\tFOUND: " << opp_position->second; */}
			else
			{
				/*qDebug() << "\t\t\tNOT FOUND: ";*/
				continue;
			}

			float w_pq = p.second->getBetaWeight(opp_position->second, *related_surf_p);
			unsigned int opp_flag = p.second->getOppositeFlag();
			float new_val = max(tmp_p * w_pq, opp_position->second->getValueByFlag(opp_flag));

			/*	qDebug() << "\t\tPopped items's flag:" << p.second->getFlag();
				qDebug() << "\t\tPopped items's Val for flag:" << tmp_p;
				qDebug() << "\t\t\tValues of opposite pole -  Out:" << opp_position->second->getValueByFlag(1) << " in: " << opp_position->second->getValueByFlag(0);
				qDebug() << "\t\t\tPre Val:" << opp_position->second->getValueByFlag(opp_flag) << " NEW VAL: " << new_val;
*/
			opp_position->second->setValueByFlag(opp_flag, new_val);

			pri_ball_q.update(opp_position);
			opp_position = pri_ball_q.find(opposite_pole);
			/*
				qDebug() << "\t\tAfter:" << opp_position->second;
				qDebug() << "\t\t\tOut:" << opp_position->second->getValueByFlag(1) << " in: " << opp_position->second->getValueByFlag(0);
				qDebug() << "\t\t\tNew Priority: " << opp_position->first;
*/
		}
		//Beta seems fine!

		//continue;

		std::set<PolarBall*> related_neighbours = neigh_graph[p.second];
		//qDebug() << "Neighbours ";
		for (std::set<PolarBall*>::iterator it = related_neighbours.begin(); it != related_neighbours.end(); ++it)
		{
			//we got the neighbours, including opposite

			//	qDebug() << "\t\t\tFINDING: " << *it;
			auto neighour_position = pri_ball_q.find(*it);
			if (neighour_position != pri_ball_q.end()) { /*qDebug() << "\t\t\t\tFOUND: " << neighour_position->second; */}
			else
			{
			//	qDebug() << "\t\t\t\tNOT FOUND: ";
				continue;
			}

			float w_pq = p.second->getAlphaWeight(neighour_position->second);
			unsigned int p_flag = p.second->getFlag();
			float new_val = max(tmp_p * w_pq, neighour_position->second->getValueByFlag(p_flag));

		/*	qDebug() << "\t\t\tPopped items's flag:" << p_flag;
			qDebug() << "\t\t\tPopped items's Val for flag:" << tmp_p;
			qDebug() << "\t\t\tValues of opposite pole -  Out:" << neighour_position->second->getValueByFlag(1) << " in: " << neighour_position->second->getValueByFlag(0);
			qDebug() << "\t\t\tPre Val:" << neighour_position->second->getValueByFlag(p_flag) << " NEW VAL: " << new_val;
*/
			neighour_position->second->setValueByFlag(p_flag, new_val);

			pri_ball_q.update(neighour_position);
			neighour_position = pri_ball_q.find(*it);
		/*	qDebug() << "\t\t\tAfter:" << neighour_position->second;
			qDebug() << "\t\t\tOut:" << neighour_position->second->getValueByFlag(1) << " in: " << neighour_position->second->getValueByFlag(0);
			qDebug() << "\t\t\tNew Priority: " << neighour_position->first;
*/
		}
	
	}

	qDebug() << "Update tims is " << pri_ball_q.upd_time.count();
	/*************************************
	*	Last step:
			Labeling	done	*
	**************************************/
	neigh_by_surf.clear();
	neigh_graph.clear();

	qDebug() << "Labeling finished";
	qDebug() << "Labeling ended in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - start).count();
	//currently 88,5 sec 21k poi
}

/* ****************************
	Power Crust calculation functions
**************************** */

void PowerDiagram::calc_power_crust()
{
	qDebug() << "Power Crust calculation";

	std::vector<face_t*> unoriented_power_crust_faces;

	for (auto& face : power_faces)
	{
		auto related_cells = face->face.related_cells;
		if (related_cells.size() > 2) qDebug() << "ERR: BAD FACE";
		if (related_cells.size() != 2)	continue;

		unsigned int cell_one = related_cells[0];
		unsigned int cell_two = related_cells[1];

		int border_face_check = (cells[cell_one].hasInnerPole() + cells[cell_two].hasInnerPole());

		if (border_face_check == 1)		//this will be a boundary face, correct its orientation
		{
			face->face.set_power_crust_flag(true);
			face->face.correct_edge_orientations();
			unoriented_power_crust_faces.push_back(&face->face);
		}
	}

	qDebug() << "We have " << unoriented_power_crust_faces.size() << " inner outer face!";

	correct_unoriented_power_crust_faces(unoriented_power_crust_faces);

	qDebug() << "Power Crust calculation finished";
}

void PowerDiagram::correct_unoriented_power_crust_faces(std::vector<face_t*>& unoriented_face_ptrs)
{
	qDebug() << "\tCorrect face orientations";

	face_t* starter_correct_face = nullptr;
	float distance_from_correct_face = FLT_MAX;
	Point outer_point = Point(3.0f*bounded.getXMax(), 3.0f*bounded.getYMax(), 3.0f*bounded.getZMax());
	Vector surface_normal;
	Vector outer_to_face;

	for (auto& power_face : unoriented_face_ptrs)
	{
		// we already oriented the edges one way
		auto edge1 = power_face->edges[0];
		auto edge2 = power_face->edges[1];

		// fill the first 3 point depends on the orienation
		Point* face_first_point;
		Point* face_second_point;
		Point* face_third_point;

		if (edge1->related_faces[power_face] == 0)
		{
			face_first_point = edge1->edge.first.second;
			face_second_point = edge1->edge.second.second;
		}
		else
		{
			face_first_point = edge1->edge.second.second;
			face_second_point = edge1->edge.first.second;
		}

		if (edge2->related_faces[power_face] == 0)
			face_third_point = edge2->edge.second.second;
		else
			face_third_point = edge2->edge.first.second;

		// calculate the actual normal
		Vector act_normal = CGAL::cross_product(*face_second_point - *face_first_point, *face_third_point - *face_first_point);
		Vector act_outer_to_face_v = *face_first_point - outer_point;
		float length_of_act_outer_to_face_v = act_outer_to_face_v.squared_length();

		if (length_of_act_outer_to_face_v < distance_from_correct_face)
		{
			distance_from_correct_face = length_of_act_outer_to_face_v;
			outer_to_face = act_outer_to_face_v;
			surface_normal = act_normal;
			starter_correct_face = power_face;
		}
	}

	// start the correction from the closest face to bound_point
	if (starter_correct_face != nullptr)
	{
		if (CGAL::scalar_product(outer_to_face, surface_normal) < 0)	// if thats obtuse angle then swap actual face orientation
			starter_correct_face->swap_edge_orientations();

		starter_correct_face->correct_boundary_orientation();			// we can go around the whole boundary because of the power_crust_flag
	}
	else
		qDebug() << "ERR: during correct unoriented normals";
}

/* ****************************
	Others
**************************** */

void PowerDiagram::addBoxPoints(const Box& box, std::set<Point>& box_points, std::vector<std::pair<weighted_point, std::vector<Point>>>& points)
{
	Point p1 = Point(box.getXMin(), box.getYMin(), box.getZMin());
	Point p2 = Point(box.getXMin(), box.getYMin(), box.getZMax());
	Point p3 = Point(box.getXMin(), box.getYMax(), box.getZMin());
	Point p4 = Point(box.getXMin(), box.getYMax(), box.getZMax());
	Point p5 = Point(box.getXMax(), box.getYMin(), box.getZMin());
	Point p6 = Point(box.getXMax(), box.getYMin(), box.getZMax());
	Point p7 = Point(box.getXMax(), box.getYMax(), box.getZMin());
	Point p8 = Point(box.getXMax(), box.getYMax(), box.getZMax());

	//TODO: 6 és 8 miatt ír warningot!

	box_points.insert(p1); points.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(p1, 0.000001f), std::vector<Point>()));
	box_points.insert(p2); points.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(p2, 0.000001f), std::vector<Point>()));
	box_points.insert(p3); points.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(p3, 0.000001f), std::vector<Point>()));
	box_points.insert(p4); points.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(p4, 0.000001f), std::vector<Point>()));
	box_points.insert(p5); points.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(p5, 0.000001f), std::vector<Point>()));
	box_points.insert(p6); points.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(p6, 0.000001f), std::vector<Point>()));
	box_points.insert(p7); points.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(p7, 0.000001f), std::vector<Point>()));
	box_points.insert(p8); points.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(p8, 0.000001f), std::vector<Point>()));
}

Box PowerDiagram::UpdateBox(const std::vector<std::pair<Pole, Point> >& weighted_points)
{
	Box res = Box(-1.0f, 1.0f);

	float minX = FLT_MAX, maxX = FLT_MIN;
	float minY = FLT_MAX, maxY = FLT_MIN;
	float minZ = FLT_MAX, maxZ = FLT_MIN;
	float x, y, z;

	for (auto it : weighted_points)
	{
		x = it.first.center->x();
		y = it.first.center->y();
		z = it.first.center->z();
	
		minX = minX > x ? x : minX;
		maxX = maxX < x ? x : maxX;
		minY = minY > -y ? -y : minY;
		maxY = maxY < -y ? -y : maxY;
		minZ = minZ > z ? z : minZ;
		maxZ = maxZ < z ? z : maxZ;	
	}

	if (!res.isInside(minX, minY, minZ) || !res.isInside(maxX, maxY, maxZ))
	{
		res = Box(	  minX < res.getXMin() ? minX : res.getXMin()
					, maxX > res.getXMax() ? maxX : res.getXMax()
					, minY < res.getYMin() ? minY : res.getYMin()
					, maxY > res.getYMax() ? maxY : res.getYMax()
					, minZ < res.getZMin() ? minZ : res.getZMin()
					, maxZ > res.getZMax() ? maxZ : res.getZMax());

		res = Box(	  res.getXMin() * 1.2f,
				      res.getXMax() * 1.2f,
				      res.getYMin() * 1.2f,
				      res.getYMax() * 1.2f,
				      res.getZMin() * 1.2f,
				      res.getZMax() * 1.2f);

	}

	return res;
}

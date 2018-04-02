#include "power_diagram.h"

float PolarBall::getAlphaWeight(PolarBall* rhs)
{
	float R = this->pole.radius;
	float r = rhs->getRadius();
	float d = (this->getPoint() - rhs->getPoint()).squared_length();

	//TODO: how measure this?
//	if (r + R - sqrt(d) > 0.000001f) return 0.0f;		//instersect shallowly

	/**************************
	Law of cosinus for calculate intersection deepness
	*************************************************/

	float lhs_cosB = (r * r + R * R - d) / (2 * r*R);
	float B_angle = acosf(lhs_cosB);
	
	float A_angle = 3.14159265359f - B_angle;
		
	if (-1 * cosf(A_angle) < 0)		//only 0.0f - 1.0f range is valid, here the intersection was more deeply
	{
	//		qDebug() << "\t\tERR: ALPHA val:";
	/*		qDebug() << "\t\t\tR: " << R << ", r: " <<r <<", d: "<< d;
			qDebug() << "\t\t\tB_angle: " << B_angle;
			qDebug() << "\t\t\tA_angle: " << A_angle;
			qDebug() << "\t\t\tlhs_cosB " << lhs_cosB;
			qDebug() << "\t\t\tgetAlphaWeight return: " << -1 * cosf(A_angle);*/
			return 0.0f;
	}

//	qDebug() << " ALPHA val calculated!";

	return -1 * cosf(A_angle);
}

float PolarBall::getBetaWeight(PolarBall* rhs, Point& surfpoint)
{
	Vector vec_one = this->getPoint() - surfpoint;
	Vector vec_two = rhs->getPoint() - surfpoint;

	/**************************
	Check for common surf point
	**********************/

	std::vector<Point> act = getSurfPoints();
	std::vector<Point> rhs_s = rhs->getSurfPoints();

	std::set<Point> act_s(act.begin(), act.end());
	std::set<Point> rhs_ss(rhs_s.begin(), rhs_s.end());

	if (act_s.find(surfpoint) == act_s.end() || rhs_ss.find(surfpoint) == rhs_ss.end())
	{
		qDebug() << "\t\tERR: There is no surf point for beta calculation!";
		return 0.0f;
	}

	/**************************
		Cos of Angle between the 2 vector
		**********************/

	float lhs = CGAL::scalar_product(vec_one, vec_two);

	float cos_val = lhs / sqrt(vec_one.squared_length() * vec_two.squared_length());

	if (-1 * cos_val < 0)		//this was never happen
	{
		qDebug() << "\t\tERR: BETA val:";
	/*	qDebug() << "\t\t\tvec1 " << vec_one.x() << ", " << vec_one.y() << ", " << vec_one.z();
		qDebug() << "\t\t\tvec2 " << vec_two.x() << ", " << vec_two.y() << ", " << vec_two.z();
		qDebug() << "\t\t\tscalar dot " << lhs;
		qDebug() << "\t\t\tBetaWeight return: " << -1 * cos_val;*/
		return 0.0f;
	}

	return -1 * cos_val;
}

//TODO: numerikus hiba miatt nem lesz teljesen sík, a pararelepipedon térfogat viszont túl pici így fölöslegesnek tartom a checket a coplanarityra
void PowerCell::calcProperNeighbours()
{
	/*for (auto it : neighbour_candidates)
	if (it.second.size() >= 3)	//van leg. 1 közös lapja, s nem egyenesre esnek a pontok
	proper_neighbours.push_back(it.first);
	*/
	proper_neighbours.clear();

	for (auto& it : neighbour_candidates)
	{
		if (it.second.size() < 3) continue;

		l_algebra::Vector center(3, 0.0f);

		for (auto& p_index : it.second)
		{
			center[0] += vertices_map[p_index]->x();
			center[1] += vertices_map[p_index]->y();
			center[2] += vertices_map[p_index]->z();
		}
		center /= it.second.size();

		std::vector<l_algebra::Vector> rows;
		for (auto& p_index : it.second)
		{
			l_algebra::Vector row(3);
			row[0] = vertices_map[p_index]->x() - center[0];
			row[1] = vertices_map[p_index]->y() - center[1];
			row[2] = vertices_map[p_index]->z() - center[2];
			rows.push_back(row);
		}
		Matrix tmp = Matrix(rows.begin(), rows.end());
		//	qDebug() << "Matrix dim: " << tmp.row_dimension() << " " << tmp.column_dimension();

		int rank = l_algebra::rank(tmp);
		//qDebug() << "Matrix rank is " << rank << " det " << l_algebra::determinant(l_algebra::transpose(tmp) * tmp );;
		proper_neighbours.push_back(it.first);

		for (int i = 0; i < rows.size() - 2; ++i)
		{
			Matrix test = Matrix(rows.begin() + i, rows.begin() + i + 3);
		//	qDebug() << "Determinants: " << l_algebra::determinant(test);
			if (l_algebra::determinant(test) > 0.0001f) qDebug() << "ERROR: too big determinant!";
		}

		if (rank == 2)
		{ 
		//	qDebug() << "Good Neigh";
		}
		else
		{
		//	qDebug() << "\tERR: Not 2D common face! :" << rank;
			/*if (it.second.size() == 3)
			{
				qDebug() << "Common Points";
				for (auto it : it.second)
					qDebug() << vertices_map[it]->x() << vertices_map[it]->y() << vertices_map[it]->z();
			}	*/				
		}						
	}
}

void PowerDiagram::calc_diagram(const std::vector<std::pair<Pole, Point> >& weighted_points)
{
	qDebug() << "Power Diagram calculation";
	/*************************************
	*	Clear	*
	**************************************/

	cells.clear();
	point_map.clear();
	pole_map.clear();
	indicies_of_power_diagram_segments.clear(); 
	inner_outer_related_cell_pairs.clear();
	regular_triangles.clear();
	power_crust_indicies.clear();

	/*************************************
	*	Box update & prepare for Weighted Points	*
	**************************************/

	qDebug() << "\tOUR BOX IS: " << bounded.getXMin() << " " << bounded.getXMax() << " " << bounded.getYMin() << " " << bounded.getYMax() << " " << bounded.getZMin() << " " << bounded.getZMax();

	qDebug() << "We got " << weighted_points.size() << " poles";

	std::map<Point, std::pair<float, std::vector<Point>>> pole_point_longest_radius_pairs;

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

	qDebug() << "\t exact: " << pole_point_longest_radius_pairs.size() << " different poles!";

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

	/*for (auto it : pole_point_longest_radius_pars)
	{
	qDebug() << "for a point " << it.second.first <<","  << " radius " << it.second.second.size() << " size surface points";
	}*/

	/*************************************
	*	Make weighted points and do Triangulation	*
	**************************************/
	qDebug() << "\tMaking weighted points";
   //locking?
	std::vector<std::pair<weighted_point, std::vector<Point>>> w_p;
	for (auto it : pole_point_longest_radius_pairs)
		w_p.push_back(std::pair<weighted_point, std::vector<Point>>(weighted_point(it.first, it.second.first* it.second.first), it.second.second));	//d*d - r -> az alap sugár négyzete kell

	pole_point_longest_radius_pairs.clear();

	qDebug() << "\tBox stuff";

	std::set<Point> box_points;
	Box pole_bound = UpdateBox(weighted_points);
	addBoxPoints(pole_bound, box_points, w_p);

	qDebug() << "\tTriangulation";

	regular_triangulation::Lock_data_structure locks(CGAL::Bbox_3(pole_bound.getXMin(), pole_bound.getYMin(), pole_bound.getZMin(), pole_bound.getXMax(), pole_bound.getYMax(), pole_bound.getZMax()), 50);
	regular_triangulation R(w_p.begin(), w_p.end(), &locks);
	//w_p.clear();

	qDebug() << R.number_of_cells() << " cells " << R.number_of_finite_cells() << " finite cells";
	qDebug() << R.number_of_vertices() << " vertices ";

	assert(R.is_valid());

	/*************************************
	*	Calculate Power Cellas	*
	**************************************/

	unsigned int pow_index = 0;
	unsigned int pol_index = 0;

	int count = 0;

	std::set<std::pair<Point*, Point*>> power_diagram_segments;

	for (regular_triangulation::Finite_vertices_iterator vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); vit++)
	{
		//skip box points
		if (box_points.count(vit->point().point()) == 1)
		{
			count++;
			continue;
		}

		//get incident cells
		std::vector<regular_cell_handle> inc_cells;
		R.incident_cells(vit, std::back_inserter(inc_cells));
	
		//check for empty
		if (inc_cells.empty())
		{
			qDebug() << "\tEMPPTY CELLS FOR SRUF VERT";
			continue;
		}

		//store pole point for reference
		std::pair<std::map<Point, unsigned int>::iterator, bool> pol_res;
		pol_res = pole_map.insert(std::pair<Point, unsigned int >(vit->point().point(), pol_index));
		if (pol_res.second == true)
		{
			pol_index++;	//not in the map
		}
		else
			qDebug() << "ERR: BAAAD, double POLE";

		//creating powercell								pole ref				get radius			related surf points
		PowerCell pc = PowerCell(const_cast<Point*>(&(pol_res.first->first)), sqrt(vit->point().weight()), vit->info());		//vissza konverzió

		std::set<Point> related_points_of_actual_inc_cells;		//for neighbours calculation
		std::map<regular_cell_handle, Point*> cell_dual_pairs;	//for segment calculation

		for (auto& cit = inc_cells.begin(); cit != inc_cells.end(); ++cit)
		{
			auto tmp = /*getBoundedDual(*/R.dual(*cit)/*)*/;

			std::pair<std::map<Point, Vertex_PowerCell_Ids>::iterator, bool> res;
			res = point_map.insert(std::pair<Point, Vertex_PowerCell_Ids >(tmp, std::pair<unsigned int, std::set<unsigned int>>(pow_index, std::set<unsigned int>())));
			if (res.second == true)
				pow_index++;	//not in the map

			const unsigned int actual_power_vertex_index = res.first->second.first;

			pc.addPowerVertex(const_cast<Point*>(&(res.first->first)), actual_power_vertex_index);		//hozzáadom a vertexet a powercellahoz
			cell_dual_pairs.insert(std::pair<regular_cell_handle, Point*>(*cit, const_cast<Point*>(&(res.first->first))));

			for (auto& related_power_cell : res.first->second.second)
			{
				//	qDebug() << related_power_cell << " < " << cells.size();
				cells[related_power_cell].addNeighbourCandidate(cells.size(), actual_power_vertex_index);		//új powercella szomszéd esély hozzáadás az aktuális vertex indexel
			}

			for (auto related_power_cell : res.first->second.second)
			{
				pc.addNeighbourCandidate(related_power_cell, actual_power_vertex_index);			//az újba berakjuk az eddigi powercella indexekhez az aktuális vertexet
			}

			related_points_of_actual_inc_cells.insert(res.first->first);
			//	res.first->second.second.insert(cells.size());		//adott vertexhez berakom az aktuális power cella indexét!
			//FAIL: beszúrtuk s újra elõjön a vertex de még nincs benn a cella!, ezért kell kintre rakni
		}

		for (auto& it : related_points_of_actual_inc_cells)		//az érintett vertexekhez, hozzáadjuk az aktuális cella indexet
			point_map[it].second.insert(cells.size());

		//---------------------------

		std::set<std::pair<Point*, Point*>> segments;

		for (auto& cit = inc_cells.begin(); cit != inc_cells.end(); ++cit)
		{
			for (int i = 0; i < 4; ++i)
			{
				regular_cell_handle tmp_cell = (*cit)->neighbor(i);
				auto is_related = cell_dual_pairs.find(tmp_cell);

				if (is_related != cell_dual_pairs.end())	//neighbour is incident cell as well
				{
					//segments of cell
					auto already_in = segments.find(std::pair<Point*, Point*>(cell_dual_pairs[tmp_cell], cell_dual_pairs[*cit]));
					if (already_in == segments.end())
						segments.insert(std::pair<Point*, Point*>(cell_dual_pairs[*cit], cell_dual_pairs[tmp_cell]));

					already_in = power_diagram_segments.find(std::pair<Point*, Point*>(cell_dual_pairs[tmp_cell], cell_dual_pairs[*cit]));
					if (already_in == power_diagram_segments.end())
						power_diagram_segments.insert(std::pair<Point*, Point*>(cell_dual_pairs[*cit], cell_dual_pairs[tmp_cell]));
				}
			}
		}

		pc.setSegmentIndicies(segments);

		addCell(pc);
	}

	setDiagramIndicies(power_diagram_segments);

	qDebug() << pole_map.size() << " poles";
	qDebug() << count << " left points";
	qDebug() << "Power Diagram base is created";
	qDebug() << "\t" << cells.size() << " cella";

	//calculate those neighbours which have common face
	for (auto& it : cells)
	{
		it.calcProperNeighbours();
	}
	int index = 0;
	
	//checking önmagának szomszédja e
	for (auto it : cells)
	{
		auto ne = it.getProperNeighbours();
		std::set<unsigned int> tmp = std::set<unsigned int>(ne.begin(), ne.end());
		if (tmp.find(index) != tmp.end())
		{
			qDebug() << " Own neighbour";
		}

		index++;
	}

	/*for (auto& it : cells)
	{
		std::stringstream ss;
		std::vector<unsigned int> tmp = it.getProperNeighbours();
		qDebug() << &it << " size " << tmp.size();
		for (auto it2 : tmp) qDebug() << "\t " << &cells[it2];
	}*/

	qDebug() << "Power Diagram cells's neighbours are calculated";

	int ct2 = 0;
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

void PowerDiagram::label_poles()
{
	qDebug() << "\nLabeling method";
	qDebug() << "\tMINIMAL BOUND BOX IS: " << bounded.getXMin() << " " << bounded.getXMax() << " " << bounded.getYMin() << " " << bounded.getYMax() << " " << bounded.getZMin() << " " << bounded.getZMax();
	qDebug() << "\n\tPrepare steps";
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
		auto actual_neighs = cell.getProperNeighbours();
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
	//	qDebug() << "\t\tPopped: " << p.second << " with prio : " << p.first;

		if (pri_ball_q.find_higher_prior(p.first))
			qDebug() << "\tPriority queue FAIL ";

		pri_ball_q.pop();

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

	/*************************************
	*	Last step:
			Labeling	done	*
	**************************************/
	neigh_by_surf.clear();
	neigh_graph.clear();

	qDebug() << "Labeling finished";
}

GLPaintFormat PowerDiagram::getCellPaintData(const int& ind, const int& nid)
{
	//cells[ind].printNeighbourCandidate();
	//cells[ind].printProperNeightbours();

	if (size() == 0) return GLPaintFormat();

	//cells[ind].calcConvexHull();
	/*return*/
	GLPaintFormat res = cells[ind].getPaintData();
	res.point_part_lengths.push_back(res.ix.size());

	int n = cells[ind].getNeighbourIndex(nid);
	if (n >= 0)
	{
	//	cells[n].calcConvexHull();
		cells[n].appendPaintDataAsNeighbour(res, ind);
		res.point_part_lengths.push_back(res.ix.size() - res.point_part_lengths[0]);
		res.points_col.push_back(QVector3D(0, 1, 0));
		res.points_col.push_back(QVector3D(1, 0, 0));
	}
	else
		res.points_col.push_back(QVector3D(1, 0, 0));

	return res;
}

GLPaintFormat PowerDiagram::getAllPaintData()
{
	GLPaintFormat res;
	for (auto it : point_map)
		res.points.push_back(it.first);

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
	for (auto it : point_map)
		res.points.push_back(it.first);

	res.ix = indicies_of_power_diagram_segments;
	qDebug() << "points " << res.points.size();
	qDebug() << "indicies_of_power_diagram_segments " << indicies_of_power_diagram_segments.size();
	res.col.push_back(QVector3D(1, 0, 0));
	return res;
}

Point PowerDiagram::getBoundedDual(Point& dual)
{
	if (bounded.isInside(dual.x(), dual.y(), dual.z()))
		return dual;

	QVector3D direction = QVector3D(dual.x(), dual.y(), dual.z());
	direction.normalize();

	float disMinX = bounded.getXMin() / direction.x();
	float disMaxX = bounded.getXMax() / direction.x();
	float disMinY = bounded.getYMin() / direction.y();
	float disMaxY = bounded.getYMax() / direction.y();
	float disMinZ = bounded.getZMin() / direction.z();
	float disMaxZ = bounded.getZMax() / direction.z();

	disMinX = disMinX  < 0 ? FLT_MAX : disMinX;
	disMaxX = disMaxX  < 0 ? FLT_MAX : disMaxX;
	disMinY = disMinY  < 0 ? FLT_MAX : disMinY;
	disMaxY = disMaxY  < 0 ? FLT_MAX : disMaxY;
	disMinZ = disMinZ  < 0 ? FLT_MAX : disMinZ;
	disMaxZ = disMaxZ  < 0 ? FLT_MAX : disMaxZ;

	float scaleX = disMinX >= disMaxX ? disMaxX : disMinX;
	float scaleY = disMinY >= disMaxY ? disMaxY : disMinY;
	float scaleZ = disMinZ >= disMaxZ ? disMaxZ : disMinZ;

	float scale = scaleX >= scaleY ? scaleY : scaleX;
	scale = scale >= scaleZ ? scaleZ : scale;
	Point tmp = Point(direction.x() * scale, direction.y() * scale, direction.z() *scale);

	//qDebug() << " NEW DUAL: " << tmp.x() << ", " << tmp.y() << ", " << tmp.z();
	return tmp;
}

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

void PowerDiagram::calc_inner_outer_pairs()
{
	int index = 0;
	for (auto it : cells)
	{
		index++;

		auto act_neighbours = it.getProperNeighbours();
		auto act_flag = it.getPolarBallPtr()->getFlag();

	//	if (act_flag == 0) qDebug() << index - 1 << " is an inner cell";

		if (act_flag == 2) continue;

		unsigned int opp_flag = 1 - act_flag;

		for (auto neigh_it : act_neighbours)
		{
			if (cells[neigh_it].getPolarBallPtr()->getFlag() != opp_flag) continue;

			std::pair<unsigned int, unsigned int> inner_outer_pair;
			if (act_flag == 0) inner_outer_pair = std::pair<unsigned int, unsigned int>(index-1, neigh_it);
			else inner_outer_pair = std::pair<unsigned int, unsigned int>(neigh_it, index - 1 );

			inner_outer_related_cell_pairs.insert(inner_outer_pair);
		}
	}

	qDebug() << "We have " << inner_outer_related_cell_pairs.size() << " related inner - outer pairs";
/*	for (auto it : inner_outer_related_cell_pairs)
		qDebug() << it;
*/

	calc_power_crust();
}

GLPaintFormat PowerDiagram::getInnerOuterPairsPaintData()
{
	GLPaintFormat res;
	std::map<Point, unsigned int> io_points;
	int io_index = 0;
	// maybe store these
	for (auto pair_it : inner_outer_related_cell_pairs)
	{
		std::vector<Point> rel_points = cells[pair_it.first].getPoints();
		for (auto poi_it : rel_points)
		{
			std::pair<std::map<Point, unsigned int>::iterator, bool> res_i;
			res_i = io_points.insert(std::pair<Point, unsigned int>(poi_it, io_index));
			if (res_i.second == true)
			{
				res.points.push_back(poi_it);
				io_index++;
			}
		}

		rel_points = cells[pair_it.second].getPoints();
		for (auto poi_it : rel_points)
		{
			std::pair<std::map<Point, unsigned int>::iterator, bool> res_i;
			res_i = io_points.insert(std::pair<Point, unsigned int>(poi_it, io_index));
			if (res_i.second == true)
			{
				res.points.push_back(poi_it);
				io_index++;
			}
		}
	}

	for (auto pair_it : inner_outer_related_cell_pairs)
	{
		std::vector<Point> rel_points = cells[pair_it.first].getPoints();
		std::vector<unsigned int> rel_points_ind = cells[pair_it.first].getSegmentIndicies();

		for (auto ind : rel_points_ind)
			res.ix.push_back(io_points[rel_points[ind]]);
	}
	
	res.point_part_lengths.push_back(res.ix.size());
	qDebug() << "Inner indexes: " << res.point_part_lengths[0];
	/*for (auto pair_it : inner_outer_related_cell_pairs)
	{
		std::vector<Point> rel_points = cells[pair_it.second].getPoints();
		std::vector<unsigned int> rel_points_ind = cells[pair_it.second].getSegmentIndicies();

		for (auto ind : rel_points_ind)
			res.ix.push_back(io_points[rel_points[ind]]);
	}
	*/
	res.point_part_lengths.push_back(res.ix.size() - res.point_part_lengths[0]);
	qDebug() << "Outer indexes: " << res.point_part_lengths[1];


	res.points_col.push_back(QVector3D(0, 1, 0));
	res.points_col.push_back(QVector3D(1, 0, 0));

	return res;
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

void PowerDiagram::calc_power_crust()
{
	qDebug() << "Power Crust calculation";

	int ct = 0;

	power_crust_indicies.clear();
	for (auto& inner_outer_pair : inner_outer_related_cell_pairs)
	{
		std::set<unsigned int> common_points_s = cells[inner_outer_pair.first].getCommonNeighbourPoints(inner_outer_pair.second);
		std::vector<unsigned int> common_points(common_points_s.begin(), common_points_s.end());

		std::vector<Point> points_in_new_coordsystem;
		Vector common_plane_normal;
		Point common_plane_barycenter = Point(0,0,0);
		
		for (auto& common_point_ind : common_points)
		{
			auto p = cells[inner_outer_pair.first].getPointFromMap(common_point_ind);
			common_plane_barycenter += Vector(p.x(), p.y(), p.z());
		//	qDebug() << p.x() << ", " << p.y() << ", " << p.z();
		}
		size_t com_p_siz = common_points.size();
		common_plane_barycenter = Point(common_plane_barycenter.x() / com_p_siz, common_plane_barycenter.y() / com_p_siz, common_plane_barycenter.z() / com_p_siz);

		common_plane_normal = CGAL::cross_product(cells[inner_outer_pair.first].getPointFromMap(common_points[0]) - common_plane_barycenter,
			cells[inner_outer_pair.first].getPointFromMap(common_points[1]) - common_plane_barycenter);

		Vector common_plane_y_axis = CGAL::cross_product(Vector(1, 0, 0), common_plane_normal);
		Vector common_plane_x_axis = CGAL::cross_product(common_plane_normal, common_plane_y_axis);

		common_plane_normal = common_plane_normal / CGAL::sqrt(common_plane_normal.squared_length());
		common_plane_y_axis = common_plane_y_axis / CGAL::sqrt(common_plane_y_axis.squared_length());
		common_plane_x_axis = common_plane_x_axis / CGAL::sqrt(common_plane_x_axis.squared_length());

	//	qDebug() << "\tNormal: " << common_plane_normal.x() << ", " << common_plane_normal.y() << ", " << common_plane_normal.z();
	//	qDebug() << "\tNormal: " << common_plane_y_axis.x() << ", " << common_plane_y_axis.y() << ", " << common_plane_y_axis.z();
	//	qDebug() << "\tNormal: " << common_plane_x_axis.x() << ", " << common_plane_x_axis.y() << ", " << common_plane_x_axis.z();
		
		std::map<unsigned int, Point> new_p;

		for (auto common_point_ind : common_points)
		{
			auto p = cells[inner_outer_pair.first].getPointFromMap(common_point_ind);
			Vector tmp = p - common_plane_barycenter;
			Point tmp_p = Point(CGAL::scalar_product(tmp, common_plane_x_axis), CGAL::scalar_product(tmp, common_plane_y_axis), 0);
		//	qDebug() << "New Point: " << tmp_p.x() << ", " << tmp_p.y() << ", " << tmp_p.z();
		//	qDebug() << "Point: " << p.x() << ", " << p.y() << ", " << p.z();
			new_p.insert(std::pair<unsigned int, Point>(common_point_ind, tmp_p));
		}

		Vector ref = new_p[new_p.begin()->first] - Point(0,0,0);
		ref = ref / CGAL::sqrt(ref.squared_length());
		auto it = new_p.begin();

		std::set<std::pair<float, unsigned int>> angle_id_s;
		for (; it != new_p.end(); ++it)
		{
			Vector tmp = it->second - Point(0, 0, 0);
			float dot = CGAL::scalar_product(ref, tmp);
			float tmp_y = (tmp.x() / ref.x())*ref.y();

			float angle;
			if ( tmp_y > tmp.y())
				angle = 2 * 3.14159265359f - acosf(dot / CGAL::sqrt((ref.squared_length() * tmp.squared_length())));
			else
				angle = acosf(dot / CGAL::sqrt((ref.squared_length() * tmp.squared_length())));

		//	qDebug() << acosf(dot / CGAL::sqrt((ref.squared_length() * tmp.squared_length())));

			angle_id_s.insert(std::pair<float, unsigned int>(angle, it->first));
		}

		
		std::vector<std::pair<float, unsigned int>> angle_id(angle_id_s.begin(), angle_id_s.end());
		for (int i = 1; i < angle_id.size() - 1; ++i)
		{
			//for triangs
		//	power_crust_indicies.push_back(0);
		//	power_crust_indicies.push_back(angle_id[i].second);
		//	power_crust_indicies.push_back(angle_id[i+1].second);

			power_crust_indicies.push_back(angle_id[0].second);
			power_crust_indicies.push_back(angle_id[i].second);
			power_crust_indicies.push_back(angle_id[i].second);
			power_crust_indicies.push_back(angle_id[i + 1].second);
			power_crust_indicies.push_back(angle_id[i + 1].second);
			power_crust_indicies.push_back(angle_id[0].second);

		}
		
	//	if( ct == 100)  return;
	//	ct++;
	}
	

	/*for (auto& inner_outer_pair : inner_outer_related_cell_pairs)
	{
		auto inner_cell = cells[inner_outer_pair.first];
		std::set<unsigned int> common_points_ix_s = inner_cell.getCommonNeighbourPoints(inner_outer_pair.second);
		std::vector<unsigned int> common_points_ix(common_points_ix_s.begin(), common_points_ix_s.end());

		for (size_t it = 1; it < common_points_ix.size()-1; ++it)
		{
			power_crust_indicies.push_back(common_points_ix[0]);
			power_crust_indicies.push_back(common_points_ix[it]);
			power_crust_indicies.push_back(common_points_ix[it]);
			power_crust_indicies.push_back(common_points_ix[it + 1]);
			power_crust_indicies.push_back(common_points_ix[it + 1]);
			power_crust_indicies.push_back(common_points_ix[0]);
		}

	}*/

	qDebug() << "Power Crust calculation finished" ;
}

GLPaintFormat PowerDiagram::getPowerCrustPaintData()
{
	GLPaintFormat res;

			//global index,  local order
	std::map<unsigned int, unsigned int> index_pair_of_vertices;
	std::map<unsigned int, Point> tmps;
	int index = 0;
	for (auto& it : point_map)
	{
		index_pair_of_vertices.insert(std::pair<unsigned int, unsigned int>(it.second.first, index));
		tmps.insert(std::pair<unsigned int, Point>(index, it.first));
		index++;
		res.points.push_back(it.first);
	}
	//qDebug() << "index vert";

	for (auto& it : power_crust_indicies)
	{ 
	//	qDebug() << index_pair_of_vertices[it];
	//	qDebug() << tmps[index_pair_of_vertices[it]].x() << ", " << tmps[index_pair_of_vertices[it]].y() << ", " << tmps[index_pair_of_vertices[it]].z();
		int a = index_pair_of_vertices.size();
		res.ix.push_back(index_pair_of_vertices[it]);
		if (a != index_pair_of_vertices.size()) qDebug() << "FAIL";
	}

	res.points_col.push_back(QVector3D(0, 1, 0));
	res.point_part_lengths.push_back(res.ix.size());
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


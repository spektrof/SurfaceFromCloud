#include "voronoi_diagram.h"

//TODO: enable multithread
void VoronoiDiagram::calc_diagram(std::vector<Point>& points)
{
	qDebug() << "Voronoi diagram calculation";
	/*************************************
	*	Clear	*
	**************************************/
	cells.clear();
	voronoi_edges.clear();
	voronoi_vertices.clear();
	indicies_of_delauney_segments.clear();

	qDebug() << "\t OUR BOX IS: " << bounded.getXMin() << " " << bounded.getXMax() << " " << bounded.getYMin() << " " << bounded.getYMax() << " " << bounded.getZMin() << " " << bounded.getZMax();
	auto start = hrclock::now();

	/*************************************
	*	Triangulation	*
	**************************************/

	delaunay_triangulation::Lock_data_structure locking_ds(CGAL::Bbox_3(bounded.getXMin(), bounded.getYMin(), bounded.getZMin(), bounded.getXMax(), bounded.getYMax(), bounded.getZMax()), 50);

	std::vector<Point> points_with_box_points = points;
	std::set<Point> box_points;

	addBoxPoints(box_points, points_with_box_points);

	auto start_2 = hrclock::now();

	delaunay_triangulation T(points_with_box_points.begin(), points_with_box_points.end(), &locking_ds);
	assert(T.is_valid());
	qDebug() << "\t" << T.number_of_cells() << " cells " << T.number_of_finite_cells() << " finite cells";
	points_with_box_points.clear();

	qDebug() << "\t\ttriangulation ended in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - start_2).count();

	start_2 = hrclock::now();
	//calcDelauneySegments(T, points_with_box_points);
	qDebug() << "\t\tdel segments got in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - start_2).count();

	/*************************************
	*	Voronoi Diagram calculation	*
	**************************************/
	/*TEST:
	GRID: 0,01
	OUTL: 4
	TIME: 2457 ms on 23707 point	-	first implementation

	Time include the triangulation calc time as well!

	with this implementation its 841 ms				-> this needs more space
							and  786 ms with set	-> this will be slower with more points
			now this is 720 ms (all)
	68k pont		2000 ms (all), triangulation ~ 300ms
					505 MB, ~ 310MB a voronoi diagram
	*/

	unsigned int vor_index = 0;

	start_2 = hrclock::now();
	
	for (delaunay_triangulation::Finite_vertices_iterator vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++)
	{
		if (box_points.count(vit->point()) == 1)	//its a bounding box point, dont need their cells
			continue;

		//-----------------------------------------------
		//Getting Incident cells
		std::vector<delaunay_cell_handle> inc_cells;
		T.incident_cells(vit, std::back_inserter(inc_cells));
	
		if (inc_cells.empty())
			continue;

		//-----------------------------------------------
		//Precalculations
		auto inc_c_size = inc_cells.size();

		std::set<delaunay_cell_handle> inc_cell_map;	//for finding cell in logM

		for (auto& it : inc_cells)		
			inc_cell_map.insert(it);
		
		//-----------------------------------------------
		//Create the current voronoi cell
		VoronoiCell vc = VoronoiCell(vit->point());

		for (auto &cell : inc_cells)
		{
			delaunay_cell_handle act_cell = cell;

			auto act_info = &act_cell->info();
			set_cell_info(T, act_cell, act_info, vor_index);

			unsigned int act_index = act_info->index;
			Point* act_vert_ptr = act_info->dual;

			//-----------------------------------------------
			//Calculate Segments
			for (int i = 0; i < 4; ++i)
			{
				delaunay_cell_handle tmp_cell = act_cell->neighbor(i);

				//if (tmp_info->an_neigh.size() > 4) qDebug() << "rlly bad";
				/*if (tmp_info->an_neigh.find(act_index) != tmp_cell->info().an_neigh.end())
				{
					vc.addVoronoiSegment(tmp_info->an_neigh[act_index]);
					continue;
				}*/

				auto nn = inc_cell_map.find(tmp_cell);			// < log 400 , with 70k point

				if (nn != inc_cell_map.end())	//both cell_handle are related to the actual point
				{
					auto tmp_info = &tmp_cell->info();

					set_cell_info(T, tmp_cell, tmp_info, vor_index);

					unsigned int tmp_index = tmp_info->index;
					Point* tmp_vert_ptr = tmp_info->dual;

					if (act_index < tmp_index)
					{
						auto edge_it_ref = act_info->edge.insert(std::pair<Point*, edge_t*>(tmp_vert_ptr, nullptr));
						if (edge_it_ref.second == true)
						{
							typename edge_p::e_ptr node = std::make_shared<edge_p>(edge_t(edge_point(act_index, act_vert_ptr), edge_point(tmp_index, tmp_vert_ptr)));
							voronoi_edges.push_back(node);
							edge_it_ref.first->second = &voronoi_edges.back()->edge;
						}
							
						vc.addVoronoiSegment(edge_it_ref.first->second);

						/*auto edge_it_ref = insert_edge(act_info->edge, edge_t(edge_point(act_index, act_vert_ptr), edge_point(tmp_index, tmp_vert_ptr)));

						vc.addVoronoiSegment(edge_it_ref);*/

						//tmp_info->an_neigh.insert(std::pair<unsigned int, edge_t*>(act_index, edge_it_ref.first->second));
						//act_info->an_neigh.insert(std::pair<unsigned int, edge_t*>(tmp_index, edge_it_ref.first->second));

						/*auto edge_it_ref = voronoi_edges.insert(edge_t(edge_point(act_cell->info().index, act_vert_ptr),
							edge_point(tmp_cell->info().index, tmp_vert_ptr)));

						vc.addVoronoiSegment(const_cast<edge_t*>(&(*edge_it_ref.first)));*/
						
					}
					else
					{
						auto edge_it_ref = tmp_info->edge.insert(std::pair<Point*, edge_t*>(act_vert_ptr, nullptr));
						if (edge_it_ref.second == true)
						{
							typename edge_p::e_ptr node = std::make_shared<edge_p>(edge_t(edge_point(tmp_index, tmp_vert_ptr), edge_point(act_index, act_vert_ptr)));
							voronoi_edges.push_back(node);
							edge_it_ref.first->second = &voronoi_edges.back()->edge;
						}
						vc.addVoronoiSegment(edge_it_ref.first->second);
						
						//too much function is slower here
						/*auto edge_it_ref = insert_edge(tmp_info->edge, edge_t(edge_point(tmp_index, tmp_vert_ptr), edge_point(act_index, act_vert_ptr)));

						vc.addVoronoiSegment(edge_it_ref);*/

						//tmp_info->an_neigh.insert(std::pair<unsigned int, edge_t*>(act_index, edge_it_ref.first->second));
						//act_info->an_neigh.insert(std::pair<unsigned int, edge_t*>(tmp_index, edge_it_ref.first->second));

						/*auto edge_it_ref = voronoi_edges.insert(edge_t(edge_point(tmp_cell->info().index, tmp_vert_ptr),
							edge_point(act_cell->info().index, act_vert_ptr)));

						vc.addVoronoiSegment(const_cast<edge_t*>(&(*edge_it_ref.first)));*/

						/*auto edge_it_ref = insert_edge(edge_map[tmp_index], edge_t(edge_point(tmp_index, tmp_vert_ptr),
																	  edge_point(act_index, act_vert_ptr)));

						vc.addVoronoiSegment(edge_it_ref);*/
					}
				}
			}
		}

		//-----------------------------------------------
		//Done, add to cells
		addCell(vc);
	}

	//T.clear();

	qDebug() << "\t\titerate over vertex ended in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - start_2).count();
	qDebug() << "\tVORONOI ENDED in: " << boost::chrono::duration_cast<milliseconds>(hrclock::now() - start).count();
	qDebug() << "\tVoronoi vertex number: " << voronoi_vertices.size();

	/*************************************
	*	END - calc_diagram

	Ready for calcule poles*
	**************************************/
}

void VoronoiDiagram::set_cell_info(const delaunay_triangulation& T, delaunay_cell_handle& cell, cell_inf* info, unsigned int& vor_index)
{
	if (info->index == -1)
	{
		info->index = vor_index;
		Point dual = T.dual(cell);
		typename _vertex::v_ptr node = std::make_shared<_vertex>(dual);
		voronoi_vertices.push_back(node);
		info->dual = &voronoi_vertices.back()->point;
		vor_index++;
	}
}

edge_t* VoronoiDiagram::insert_edge(std::map<Point*, edge_t*>& edge_list, edge_t& edge_candidate)
{
	auto is_in = edge_list.insert(std::pair<Point*, edge_t*>(edge_candidate.second.second , nullptr));
	
	if (is_in.second == true)
	{
		typename edge_p::e_ptr node = std::make_shared<edge_p>(edge_candidate);
		voronoi_edges.push_back(node);
		is_in.first->second = &voronoi_edges.back()->edge;
	}

	return is_in.first->second;
}

GLPaintFormat VoronoiDiagram::getCellPaintData(const int& ind)
{
	if (size() == 0) return GLPaintFormat();
	//cells[ind].calcConvexHull();
	return cells[ind].getPaintData();
}

GLPaintFormat VoronoiDiagram::getCellWDualPaintData(const int& ind, const int& ind2)
{
	return GLPaintFormat();
}

GLPaintFormat VoronoiDiagram::getPolesSurfPaintData()
{
	if (cells.size() == 0) return GLPaintFormat();

	GLPaintFormat res;

	for (auto it : cells)
		res.centers_with_radius.push_back(std::pair<Point, float>(it.getSurfPoint(), 0.05f));

	res.col.push_back(QVector3D(148.0f / 255.0f, 0, 211.0f / 255.0f));
	res.center_part_lengths.push_back(res.centers_with_radius.size());

	for (auto it : cells)
		it.getPolesPaintData(res.centers_with_radius);

	res.col.push_back(QVector3D(0, 0, 1));
	res.center_part_lengths.push_back(res.centers_with_radius.size() - res.center_part_lengths[0]);

	return res;
}

GLPaintFormat VoronoiDiagram::getDelauneySegmentPaintData(GLPaintFormat& res)
{
	if (indicies_of_delauney_segments.size() == 0) res;

	res.points.push_back(Point(bounded.getXMin(), bounded.getYMin(), bounded.getZMin()));
	res.points.push_back(Point(bounded.getXMin(), bounded.getYMin(), bounded.getZMax()));
	res.points.push_back(Point(bounded.getXMin(), bounded.getYMax(), bounded.getZMin()));
	res.points.push_back(Point(bounded.getXMin(), bounded.getYMax(), bounded.getZMax()));
	res.points.push_back(Point(bounded.getXMax(), bounded.getYMin(), bounded.getZMin()));
	res.points.push_back(Point(bounded.getXMax(), bounded.getYMin(), bounded.getZMax()));
	res.points.push_back(Point(bounded.getXMax(), bounded.getYMax(), bounded.getZMin()));
	res.points.push_back(Point(bounded.getXMax(), bounded.getYMax(), bounded.getZMax()));

	res.ix = indicies_of_delauney_segments;
	res.col.push_back(QVector3D(1, 0, 0));
	return res;
}

GLPaintFormat VoronoiDiagram::getVoronoiDiagramBySegmentPaintData()
{
	if (voronoi_vertices.empty()) return GLPaintFormat();

	GLPaintFormat res;
	
	for (auto it : voronoi_vertices)
		res.points.push_back(it->point);

	for (auto it : voronoi_edges)
	{
		res.ix.push_back(it->edge.first.first);
		res.ix.push_back(it->edge.second.first);
	}

	res.col.push_back(QVector3D(1, 0, 0));

	qDebug() << "Voronoi Diagram indicies: " << res.ix.size();

	return res;
}

void VoronoiDiagram::calcDelauneySegments(const delaunay_triangulation& T, const std::vector<Point>& points_with_box_points)
{
	std::map<Point, unsigned int> tmp_points;
	int index = 0;
	for (auto it : points_with_box_points)
	{
		std::pair<std::map<Point, unsigned int>::iterator, bool> res;
		res = tmp_points.insert(std::pair<Point, unsigned int>(it, index));
		if (res.second == true) index++;		//just for test, we know that this if condition wont be true
	}

	for (delaunay_triangulation::All_cells_iterator fit = T.all_cells_begin(); fit != T.all_cells_end(); ++fit)
	{
		Point vertex0 = fit->vertex(0)->point();
		Point vertex1 = fit->vertex(1)->point();
		Point vertex2 = fit->vertex(2)->point();
		Point vertex3 = fit->vertex(3)->point();

		unsigned int a = tmp_points[vertex0];
		unsigned int b = tmp_points[vertex1];
		unsigned int c = tmp_points[vertex2];
		unsigned int d = tmp_points[vertex3];

		indicies_of_delauney_segments.push_back(a);
		indicies_of_delauney_segments.push_back(b);

		indicies_of_delauney_segments.push_back(b);
		indicies_of_delauney_segments.push_back(c);

		indicies_of_delauney_segments.push_back(c);
		indicies_of_delauney_segments.push_back(a);

		indicies_of_delauney_segments.push_back(a);
		indicies_of_delauney_segments.push_back(d);

		indicies_of_delauney_segments.push_back(b);
		indicies_of_delauney_segments.push_back(d);

		indicies_of_delauney_segments.push_back(c);
		indicies_of_delauney_segments.push_back(d);
	}
}

void VoronoiDiagram::addBoxPoints(std::set<Point>& box_points, std::vector<Point>& points)
{
	Point p1 = Point(bounded.getXMin(), bounded.getYMin(), bounded.getZMin());
	Point p2 = Point(bounded.getXMin(), bounded.getYMin(), bounded.getZMax());
	Point p3 = Point(bounded.getXMin(), bounded.getYMax(), bounded.getZMin());
	Point p4 = Point(bounded.getXMin(), bounded.getYMax(), bounded.getZMax());
	Point p5 = Point(bounded.getXMax(), bounded.getYMin(), bounded.getZMin());
	Point p6 = Point(bounded.getXMax(), bounded.getYMin(), bounded.getZMax());
	Point p7 = Point(bounded.getXMax(), bounded.getYMax(), bounded.getZMin());
	Point p8 = Point(bounded.getXMax(), bounded.getYMax(), bounded.getZMax());

	box_points.insert(p1); points.push_back(p1);
	box_points.insert(p2); points.push_back(p2);
	box_points.insert(p3); points.push_back(p3);
	box_points.insert(p4); points.push_back(p4);
	box_points.insert(p5); points.push_back(p5);
	box_points.insert(p6); points.push_back(p6);
	box_points.insert(p7); points.push_back(p7);
	box_points.insert(p8); points.push_back(p8);
}

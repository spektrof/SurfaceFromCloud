//#include "kdtree.h"
#include "surfacegenerator.h"

GLPaintFormat SurfaceGenerator::getPaintFromPolyhedron(point_uv_map& points)
{
	GLPaintFormat res;

	std::map<unsigned int, Point> tmps;
	int ix = 0;

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		res.points.push_back(it->point());
		tmps.insert(std::pair<unsigned int, Point>(ix, it->point()));
		ix++;
	}

	Index index = Index(mesh.vertices_begin(), mesh.vertices_end());

	for (FCI fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi)
	{
		HFCC hc = fi->facet_begin();
		HFCC hc_end = hc;

		do {
			unsigned int draw_index = index[VCI(hc->vertex())];
			res.ix.push_back(draw_index);
			++hc;
			res.points_col.push_back(QVector3D(0, 0, 0));

		} while (hc != hc_end);

	}

	res.col.push_back(QVector3D(1, 0, 0));

	return res;

}

#ifdef ENABLE_CGAL_FILTER
#include <CGAL/compute_average_spacing.h>

void SurfaceGenerator::compute_average_spacing()
{
	const unsigned int nb_neighbors = 6; // 1 ring
	average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(
		points.begin(), points.end(),
		filter_point_map(),
		//CGAL::Nth_of_tuple_property_map<1, IndexedPointWithColorTuple>(),
		nb_neighbors);

#ifdef ENABLE_DEBUG
	qDebug() << "\tAverage spacing: " << average_spacing;
#endif
}
#endif

void SurfaceGenerator::add_tex_or_color_info(GLPaintFormat& res)
{
	kdTree tree(1, 2);
	ProcessReturns build_p = tree.build_process(points, cloud.getBox());
	tree.set_kn(9);

	if (build_p != PROCESS_DONE) return;
	qDebug() << "Success kdtree";

	//----------------------------------
	// calculate triangle's uv coordianetes
	//textúra pixelhez van nem indexhez!
	cam_coords camera_uv_coords_in_order;
	camera_uv_coords_in_order.resize(res.points.size());

	std::vector<uv_map> uv_cache;
	uv_cache.resize(res.points.size());
	qDebug() << "WE started with " << res.points.size() << " point";

	std::vector<index_point> one_triangle_points;
	one_triangle_points.clear();
	std::vector<std::vector<std::pair<float, float>>> camera_separeted_uv_coords;
	std::vector<std::vector<unsigned int>> camera_separeted_index_lists;
	camera_separeted_index_lists.resize(4);		//TODO: replace 4 with variable!
	camera_separeted_uv_coords.resize(4);

	int tt = 0;

	for (auto& index : res.ix)
	{
		one_triangle_points.push_back(std::pair<unsigned int, Point>(index, res.points[index]));

		if (one_triangle_points.size() == 3)
		{
			calc_triangle_details(tree, one_triangle_points, camera_uv_coords_in_order, uv_cache, camera_separeted_index_lists, res.points, tt);
		}

	}


	// fill uv coord data
	res.point_part_lengths.clear();

	for (auto& uv_p : camera_uv_coords_in_order)
	{
		int ca = uv_p.cur.first;
		auto uv = uv_p.cur.second;

		if (uv_p.cur.first == -1) qDebug() << "BAAAAAAAAAAAAAAAAAAAAAAAAAD";
		res.uv_coords.push_back(uv_p.cur.second);

		for (auto& it : uv_p.other_inds)
		{
			auto a = camera_uv_coords_in_order[it];
			if (a.cur.first == ca)	qDebug() << "SAME CAMERA - err";
			if (a.cur.second == uv)	qDebug() << "SAME uc - err??";
			if (a.other_inds.size() != 0)	qDebug() << "NOT NULL others";
		}
	//	qDebug() << "Texture from cam " <<  camera_related_uv.size();
	}

	// fill indicies data by cameras
	res.ix.clear();

	for (auto& camera_related_index : camera_separeted_index_lists)
	{
		for (auto& ind : camera_related_index)
			res.ix.push_back(ind);

		res.point_part_lengths.push_back(camera_related_index.size());
		qDebug() << "Texture from cam " << camera_related_index.size();

	}

	qDebug() << tt << " triangle point was correct from " << res.ix.size() / 3;
	qDebug() << "texture added" << res.uv_coords.size() << " to " << res.points.size() << " point";

}

/*****************
calculation new uv coordinates
*****************/

void SurfaceGenerator::calc_triangle_details(kdTree& tree, std::vector<index_point>& one_triangle_points, cam_coords& camera_uv_coords_in_order, std::vector<uv_map>& uv_cache,
											 std::vector<std::vector<unsigned int>>& camera_separeted_index_lists, std::vector<Point>& res_points, int& tt)
{
	std::pair<texture_tuple, bool> texture_uv_result = calculate_uv_coords_for_triangle(tree, one_triangle_points, camera_uv_coords_in_order, uv_cache, tt);

	if (texture_uv_result.second == true)
	{
	//	qDebug() << "OK";
		fill_uv_list(camera_uv_coords_in_order, texture_uv_result.first, res_points);
		fill_index_list(camera_separeted_index_lists, texture_uv_result.first.first, one_triangle_points);
		one_triangle_points.clear();
	}
	else
	{
	//	qDebug() << "Partition triangle";
		Point p1 = one_triangle_points[0].second;
		Point p2 = one_triangle_points[1].second;
		Point p3 = one_triangle_points[2].second;
		unsigned int i1 = one_triangle_points[0].first;
		unsigned int i2 = one_triangle_points[1].first;
		unsigned int i3 = one_triangle_points[2].first;

		Point p12 = (p1 + p2) / 2.0f;
		Point p23 = (p2 + p3) / 2.0f;
		Point p13 = (p1 + p3) / 2.0f;

	//	qDebug() << "New poi1: " << p12;
	//	qDebug() << "New poi2: " << p23;
	//	qDebug() << "New poi3: " << p13;

		unsigned int i12 = res_points.size();
		res_points.push_back(p12);
		unsigned int i23 = res_points.size();
		res_points.push_back(p23);
		unsigned int i13 = res_points.size();
		res_points.push_back(p13);
	//	qDebug() << "\tpart trip1 " << res_points.size()-3 << " " << p12;
	//	qDebug() << "\tpart trip2 " << res_points.size()-2 << " " << p23;
	//	qDebug() << "\tpart trip3 " << res_points.size()-1 << " " << p13;

		camera_uv_coords_in_order.push_back(UV(-1, std::pair<float, float>(0.0f, 0.0f)));
		camera_uv_coords_in_order.push_back(UV(-1, std::pair<float, float>(0.0f, 0.0f)));
		camera_uv_coords_in_order.push_back(UV(-1, std::pair<float, float>(0.0f, 0.0f)));

		uv_cache.push_back(uv_map());
		uv_cache.push_back(uv_map());
		uv_cache.push_back(uv_map());

		one_triangle_points.clear();
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i1, p1));
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i12, p12));
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i13, p13));
		calc_triangle_details(tree, one_triangle_points, camera_uv_coords_in_order, uv_cache, camera_separeted_index_lists, res_points, tt);
		one_triangle_points.clear();
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i12, p12));
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i2, p2));
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i23, p23));
		calc_triangle_details(tree, one_triangle_points, camera_uv_coords_in_order, uv_cache, camera_separeted_index_lists, res_points, tt); 
		one_triangle_points.clear();
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i12, p12));
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i23, p23));
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i13, p13));
		calc_triangle_details(tree, one_triangle_points, camera_uv_coords_in_order, uv_cache, camera_separeted_index_lists, res_points, tt); 
		one_triangle_points.clear();
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i23, p23));
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i3, p3));
		one_triangle_points.push_back(std::pair<unsigned int, Point>(i13, p13));
		calc_triangle_details(tree, one_triangle_points, camera_uv_coords_in_order, uv_cache, camera_separeted_index_lists, res_points, tt);
	//	qDebug() << "Partition end";

	}
	
}


std::pair<SurfaceGenerator::texture_tuple, bool> SurfaceGenerator::calculate_uv_coords_for_triangle(kdTree& tree, std::vector<index_point>& triangle, cam_coords& camera_uv_coords_in_order, std::vector<uv_map>& uv_cache, int& tt)
{
	std::vector<std::pair<unsigned int, std::pair<float, float>>> uvs;

	//remove these
	kdTree::nearest_node closest;
	tree.nearest_p(&triangle[0].second, tree.getRootPtr(), closest);
	uv_map possibilites = closest.neighb->point->second;
	int16_t cam = possibilites[0].first;
	////////////////////
							//point index  related uvs
	std::vector<std::pair<index_point, uv_map>> paired_point_uvs;

	/*auto pq = closest_p;
	qDebug() << "Nearest: " << closest_p2.neighb->point->first << " to " << *triangle[0].second;
	while(!pq.empty())
	{
		auto el = pq.top();
		pq.pop();
		qDebug() << el.first << " " <<  el.second->point->first << " - " << el.second->point->second;
	}
	qDebug() << "-----------";
	*/
	for (auto& poi : triangle)
	{
		//kdTree::nearest_node closest;
		//tree.nearest_p(poi.second, tree.getRootPtr(), closest);

		auto cache = &uv_cache[poi.first];
		if (cache->size() != 0)
		{
		//	qDebug() << "Read cached value";
			paired_point_uvs.push_back(std::pair<index_point, uv_map>(index_point(poi.first, poi.second), *cache));
			continue;
		}

		kdTree::kNearest_queue closest_points;
		tree.knearest_p(&poi.second, tree.getRootPtr(), closest_points);
	
		uv_map modified_uvs = calc_new_uv_coords_from_knearest(closest_points, poi.second);
		std::pair<index_point, uv_map> index_modified_uvs = std::pair<index_point, uv_map>(index_point(poi.first, poi.second), modified_uvs);

		paired_point_uvs.push_back(index_modified_uvs);
		*cache = modified_uvs;

		//---------------------------
		//uv_map possibilites = closest.neighb->point->second;
		//uvs.push_back(std::pair<int, std::pair<float, float>>(poi.first, possibilites[0].second));
	}

	std::pair<texture_tuple, bool> result = calc_camera_based_uv(paired_point_uvs, camera_uv_coords_in_order, tt);

	//cache map size correction
	for (auto& it : result.first.second)
	{
		unsigned int ind = it.first.first;
		while (ind >= uv_cache.size())
		{
			uv_cache.push_back(uv_map());
		}
	}

	return result;

	//texture_tuple uv_coords(cam, uvs);
	//return uv_coords;
}

uv_map SurfaceGenerator::calc_new_uv_coords_from_knearest(kdTree::kNearest_queue& closest_points, Point& triangle_point)
{
	uv_map uv_separated_by_source;

				// camera			distance			uv with its point
	std::vector< std::set<std::pair<float,			point_uv>> > separated_uvs;
	separated_uvs.resize(4);

	kdTree::kNearest_queue closest_points_t = closest_points;

	//getting all the closest points
	int test_uv_c = 0;

	while (!closest_points_t.empty())
	{
		auto pri_element = closest_points_t.top();
		closest_points_t.pop();

		float distance = pri_element.first;
		Point related_surf_point = pri_element.second->point->first;
		auto related_uvs_from_cameras = pri_element.second->point->second;

		//getting all the related uvs
		for (auto& uv : related_uvs_from_cameras)
		{
			unsigned int camera = uv.first;
			separated_uvs[camera].insert(std::pair<float, point_uv>(distance, point_uv(related_surf_point, uv.second)));
			test_uv_c++;
		}
	}
	if (test_uv_c > 12) qDebug() << "**********************************we are cool";

	/*for (int i = 0; i < 4; ++i)
	{
		qDebug() << i << " camera: ";
		for (auto& it : separated_uvs[i])
		{
			qDebug() << "\t " << it;
		}
	}
	qDebug() << "----------------------\n";
	*/
	
	//getting the 3 closest point to the triangle point
	for (int i = 0; i < 4; ++i)
	{
		auto uv_set = separated_uvs[i];

		if (uv_set.size() < 3) continue;

		std::vector<point_uv> closest_three;
		int j = 0;
		for (auto& it = uv_set.begin(); j < 3; ++it, ++j)
			closest_three.push_back(it->second);

		//calculate the new possible uv coordinate for our triangle point
		std::pair<float, float> new_uv = calc_new_uv_from_closest_three(closest_three, triangle_point);
		//qDebug() << "\t\tNew UV is " << new_uv;
		uv_separated_by_source.push_back( UV(i, new_uv) );
	}

	return uv_separated_by_source;
}

std::pair<float, float> SurfaceGenerator::calc_new_uv_from_closest_three(std::vector<point_uv>& closest_three, Point triangle_point)
{
	//qDebug() << "\t Triangle point" << *triangle_point;
	/*for (auto& it : closest_three)
	qDebug() << "\t\t" << it;
	qDebug() << "-**************";
	*/
	Vector da = closest_three[0].first - triangle_point;
	Vector db = closest_three[1].first - triangle_point;
	Vector dc = closest_three[2].first - triangle_point;

	float d_da = sqrt(da.squared_length());
	float d_db = sqrt(db.squared_length());
	float d_dc = sqrt(dc.squared_length());

	if (d_da == 0) return closest_three[0].second;
	if (d_db == 0) return closest_three[1].second;
	if (d_dc == 0) return closest_three[2].second;

	QVector3D closest_distances_from_point (1.0f / d_da, 1.0f / d_db, 1.0f / d_dc);
	closest_distances_from_point.normalize();

	std::pair<float, float> c_a = closest_three[0].second;
	std::pair<float, float> c_b = closest_three[1].second;
	std::pair<float, float> c_c = closest_three[2].second;
//	qDebug() << "\t\t" << c_a << " - " << c_b << " - " << c_c;

	float nd_a = closest_distances_from_point.x();
	float nd_b = closest_distances_from_point.y();
	float nd_c = closest_distances_from_point.z();
//	qDebug() << "\t\t" << nd_a << " - " << nd_b << " - " << nd_c;

	float sum = nd_a + nd_b + nd_c;
	std::pair<float, float> new_uv = (c_a * nd_a + c_b * nd_b + c_c * nd_c) / std::pair<float,float>(sum, sum);
	return new_uv;
}

/*****************
	calculation which camera data what we will use
*****************/

std::pair<SurfaceGenerator::texture_tuple, bool> SurfaceGenerator::calc_camera_based_uv(const pointuvs_v& paired_point_uvs, cam_coords& camera_uv_coords_in_order, int& tt)
{
	/*for (auto& it : paired_point_uvs)
	{
		qDebug() << it.first;

		for (auto& uv : it.second)
		{
			qDebug() << "\t" << qSetRealNumberPrecision(15) << uv.first << " - " << uv.second;
		}
	}
	qDebug() << "------------------------";
	*/
	std::vector< std::vector< std::pair<index_point, std::pair<float, float> > > > uvs_by_texture;
	uvs_by_texture.resize(4);
	std::vector<index_point> point_indicies;

	for (auto& point_uv : paired_point_uvs)
	{
		unsigned int point_index = point_uv.first.first;
		point_indicies.push_back(index_point(point_index, point_uv.first.second));

		for (auto& uv : point_uv.second)
		{
			uvs_by_texture[uv.first].push_back(std::pair<index_point, std::pair<float, float> >(index_point(point_index, point_uv.first.second), uv.second));
		}
	}

	/*for (int i = 0; i < 4; ++i)
		qDebug() << "\t" << i << " cam : " << uvs_by_texture[i].size() << " - " <<  uvs_by_texture[i];*/
	//////////////////////////////////
	//We choose the same coord if we can!
	//////////////////////////////////

	std::set< count_index_link, std::greater<count_index_link>> prefered_tex_file;		//based on previous calculations
	prefered_tex_file = get_tex_file_chances(point_indicies, camera_uv_coords_in_order);

	for (auto& tex : prefered_tex_file)
	{
		if (tex.first > 3) qDebug() << "TEX FILE CALC BAAAAD" << tex.first << " " << point_indicies.size();	//... this occurs
		auto related_uv_map = uvs_by_texture[tex.second];
		if (related_uv_map.size() == 3)
		{
			tt++;
			index_correction(related_uv_map, camera_uv_coords_in_order, tex.second);		//calculate that we need a second uv for the same point or not
			return std::pair<texture_tuple, bool>(texture_tuple(tex.second, related_uv_map), true);
		}
	}

	/*int i = 1;
	auto& tmp = uvs_by_texture[0];
	while (tmp.size() < 3)
	{
		for (auto& it : uvs_by_texture[i])
		{
			tmp.push_back(it);
			if (tmp.size() == 3) break;
		}
		++i;
	}

	return std::pair<texture_tuple, bool>(texture_tuple(0, tmp), true);*/
	return std::pair<texture_tuple, bool>(texture_tuple(), false);
}

std::set<SurfaceGenerator::count_index_link, std::greater<SurfaceGenerator::count_index_link>> SurfaceGenerator::get_tex_file_chances(std::vector<index_point>& point_indicies, cam_coords& camera_uv_coords)
{
	//több pontnak van 0-sa, 2 0 van alapból és 2 jön az added_ bõl... ? mért
	std::vector<int> number_of_tex;
	number_of_tex.resize(4);
	for (auto& it : number_of_tex) it = 0;

	for (auto& ind : point_indicies)
	{
		int cam = camera_uv_coords[ind.first].cur.first;
		int actc = cam;
		if (cam != -1)
		{
			number_of_tex[cam]++;
		//	qDebug() << "\t\tOk " << cam << " inc";
		}

		std::vector<unsigned int> other_inds = camera_uv_coords[ind.first].other_inds;
		if (other_inds.size() > 2) { qDebug() << "BAAAD, TOO MUCH INDS"; }

		for (auto& other_ind : other_inds)
		{
			int cam = camera_uv_coords[other_ind].cur.first;

			if (camera_uv_coords[other_ind].other_inds.size() != 0) qDebug() << "Got it more then once!!! BAAD";
			if (cam == actc) qDebug() << "ERR: MULTIPLE STORE";
			if (cam != -1)
			{
				number_of_tex[cam]++;
			//	qDebug() << "\t\tOk ADDED tex" << cam << " inc";
			}
		}
	}

	std::set<count_index_link, std::greater<count_index_link>> res;

	for (int i = 0; i < 4; ++i)
	{
		res.insert(count_index_link(number_of_tex[i], i));
	}

	/*qDebug() << "hajlandosag: ";
	for (auto& it : res)
		qDebug() << "\t" << it;
*/
	return res;
}

void SurfaceGenerator::index_correction(std::vector<std::pair<index_point, std::pair<float, float>>>& related_uv_map, cam_coords& camera_uv_coords_in_order, const int& camera)
{
	unsigned int act_size = camera_uv_coords_in_order.size();
	//qDebug() << "Finding indexes for camera " << camera << " act siz " << camera_uv_coords_in_order.size();
	for (auto& ind_uv : related_uv_map)
	{
		unsigned int index = ind_uv.first.first;						
		int act_camera = camera_uv_coords_in_order[index].cur.first;		
																		
		if (act_camera == -1) continue;
		if (act_camera == camera)	continue;	//index is fine
	
		std::vector<unsigned int> other_inds = camera_uv_coords_in_order[index].other_inds;

		if (other_inds.size() > 2) qDebug() << "BAAAD";

		bool found = false;
		for (auto& other_ind : other_inds)
		{
			int other_cam = camera_uv_coords_in_order[other_ind].cur.first;
			if (other_cam == act_camera)
			{
				qDebug() << "BAAD";
				continue;
			}

			if (camera_uv_coords_in_order[other_ind].other_inds.size() != 0) qDebug() << "SZAAAR";

			if (other_cam == camera)
			{
				ind_uv.first.first = other_cam;			//index correction
				found = true;
				break;
			}
		}

		if (found) continue;
		//we couldnt find any good camera, we should push it

		camera_uv_coords_in_order[index].other_inds.push_back(act_size);
		ind_uv.first.first = act_size;					//index correction
		act_size++;
	}

//	qDebug() << "*************************";

}

/*****************
	other main texturere functions
*****************/

void SurfaceGenerator::fill_uv_list(cam_coords& camera_uv_coords, const texture_tuple& texture_uv_result, std::vector<Point>& r_points)
{
	int16_t cam_index = texture_uv_result.first;

	for (auto& cam_uv_pair : texture_uv_result.second)
	{
		unsigned int poi_index = cam_uv_pair.first.first;
	//	qDebug() << " res - poi index: " << poi_index << " uv " << cam_uv_pair.first.second << " " << cam_index << " cam " <<  cam_uv_pair.second;
		if (poi_index >= camera_uv_coords.size())	//new texture coord for same point
		{
			camera_uv_coords.push_back(uv_list(UV(cam_index, cam_uv_pair.second)));
		//	qDebug() << "finally added poi " << r_points.size() << " " << cam_uv_pair.first.second;
			r_points.push_back(cam_uv_pair.first.second);
		//	qDebug() << "added additional poi";
			continue;
		}

		if (camera_uv_coords[poi_index].cur.first == -1)
		{
			camera_uv_coords[poi_index] = uv_list(UV(cam_index, cam_uv_pair.second));
		}
	}
}

void SurfaceGenerator::fill_index_list(index_by_cameras_v& camera_separeted_index_lists, const int16_t& camera_index, const std::vector<index_point>& triangle)
{
	for (auto& point : triangle)
	{
		camera_separeted_index_lists[camera_index].push_back(point.first);
	}
}

void SurfaceGenerator::update_point_set(std::vector<Point>& points, cam_coords& camera_uv_coords_in_order, const index_pointptr_v&)
{

}
/*void calcConvexHull()
{
if (local_draw_indicies.size() != 0) return;

local_draw_indicies.clear();	//unnecessary
global_draw_indicies.clear();

Polyhedron tmp_mesh;
std::vector<Point> tmp_vertex;

for (auto it : vertices_map)
tmp_vertex.push_back(*it.first);

CGAL::convex_hull_3(tmp_vertex.begin(), tmp_vertex.end(), tmp_mesh);

//std::vector<unsigned int> p_index;

for (auto it = tmp_mesh.vertices_begin(); it != tmp_mesh.vertices_end(); it++)
{
auto locate = vertices_map.begin();
while (*(locate->first) != it->point())
{
locate++;
}

voronoi_vertices.push_back(locate->first);
//p_index.push_back(locate->second);
}
//qDebug() << "MESH have: " << tmp_mesh.size_of_vertices() << " : " << tmp_mesh.size_of_facets() << "\n";

Index index = Index(tmp_mesh.vertices_begin(), tmp_mesh.vertices_end());

for (FCI fi = tmp_mesh.facets_begin(); fi != tmp_mesh.facets_end(); ++fi)
{
HFCC hc = fi->facet_begin();
HFCC hc_end = hc;

do {
local_draw_indicies.push_back(index[VCI(hc->vertex())]);
//global_draw_indicies.push_back(p_index[index[vertices_map[hc->vertex()]]]);
++hc;
} while (hc != hc_end);

}

tmp_mesh.clear();
tmp_vertex.clear();
}
*/


/*
void calcConvexHull()
{
if (local_draw_indicies.size() != 0) return;

local_draw_indicies.clear();	//unnecessary

Polyhedron tmp_mesh;
std::vector<Point> tmp_vertex;

for (auto it : vertices_map)
tmp_vertex.push_back(*it.second);

CGAL::convex_hull_3(tmp_vertex.begin(), tmp_vertex.end(), tmp_mesh);

//std::vector<unsigned int> p_index;

for (auto it = tmp_mesh.vertices_begin(); it != tmp_mesh.vertices_end(); it++)
{
auto locate = vertices_map.begin();
while (*(locate->second) != it->point())
{
locate++;
}

power_vertices.push_back(locate->second);
//p_index.push_back(locate->second);
}
qDebug() << "POWER CELL have: " << tmp_mesh.size_of_vertices() << " : " << tmp_mesh.size_of_facets() << "\n";

Index index = Index(tmp_mesh.vertices_begin(), tmp_mesh.vertices_end());

for (FCI fi = tmp_mesh.facets_begin(); fi != tmp_mesh.facets_end(); ++fi)
{
HFCC hc = fi->facet_begin();
HFCC hc_end = hc;

do {
local_draw_indicies.push_back(index[VCI(hc->vertex())]);
//global_draw_indicies.push_back(p_index[index[vertices_map[hc->vertex()]]]);
++hc;
} while (hc != hc_end);

}

tmp_mesh.clear();
tmp_vertex.clear();
}
*/


/*for (auto& it : cells)
{
std::pair<std::map<PolarBall*, std::vector<Point>>::iterator, bool> res;
res = surfs_by_poles.insert(std::pair<PolarBall*, std::vector<Point>>(it.getPolarBallPtr(), std::vector<Point>()));
res.first->second.push_back(it.getRelatedSurfPoint());
}*/

/*qDebug() << "\\/\\/\\/\\/\\/\\/\\/\\/\\/\nPrinting priority queue";

tmp_q = pri_ball_q;

while (!tmp_q.empty())
{
auto tmp = tmp_q.top();
tmp_q.pop();
qDebug() << tmp.first << " - " << tmp.second;
}








Point VoronoiDiagram::getBoundedDual(Point& dual)
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


void VoronoiDiagram::add_dual_as_voronoi_vertex(const Point& dual, std::vector<Point*>& dual_pointers, std::vector<unsigned int>& dual_places, const unsigned int& dual_act_ind, std::map<Point, point_identifier>& points_map, unsigned int& vor_index)
{
std::pair<std::map<Point, point_identifier>::iterator, bool> already_in;		//this inc cell related to more than 1 surface point
already_in = points_map.insert(std::pair<Point, point_identifier>(dual, point_identifier()));

if (already_in.second == true)
{
typename _vertex::v_ptr node = std::make_shared<_vertex>(dual);

//	last_voronoi_vertices->next = node;
//	last_voronoi_vertices = last_voronoi_vertices->next;

voronoi_vertices.push_back(node);

//	already_in.first->second = point_identifier(&(last_voronoi_vertices->point), vor_index);
already_in.first->second = point_identifier(&(voronoi_vertices.back()->point), vor_index);
vor_index++;
	}

	dual_pointers[dual_act_ind] = already_in.first->second.first;
	dual_places[dual_act_ind] = already_in.first->second.second;
}

void VoronoiDiagram::add_cell_voronoi_vertex(VoronoiCell& vc, Point* new_vert, std::map<Point*, unsigned int>& cell_verts, unsigned int& cell_voronoi_v_index, unsigned int& act_vor_ind_of_cell)
{
std::pair<std::map<Point*, unsigned int>::iterator, bool> cell_vert_pos;
cell_vert_pos = cell_verts.insert(std::pair<Point*, unsigned int>(new_vert, cell_voronoi_v_index));
if (cell_vert_pos.second == true)
{
//	vc.addVoronoiVertex(new_vert);		//should be different! -> TEST
cell_voronoi_v_index++;
}

act_vor_ind_of_cell = cell_vert_pos.first->second;
}

*/

/*unsigned int triangle_index = 0;
int cam = 0;

for (auto& index : res.ix)
{
kdTree::nearest_node closest_p;

tree.nearest_p(&res.points[index], tree.getRootPtr(), closest_p);
//	qDebug() << "Closest Point to " << res.points[index] << " is " << closest_p.neighb->point->first << " uv coords are : " << closest_p.neighb->point->second[0];

camera_uv_coords_in_order[index].first = closest_p.neighb->point->second[0].first;
camera_uv_coords_in_order[index].second = closest_p.neighb->point->second[0].second;
if (one_triangle_points.size() == 3)
{

texture_tuple texture_uv_result = calculate_uv_coords_for_triangle(tree, one_triangle_points);
qDebug() << "\tcamera is : " << texture_uv_result.first;
int i = 0;
for (auto& it : texture_uv_result.second)
{
qDebug() <<  "\t\t" << it.first << " ind " << " poi " << *one_triangle_points[i].second << " s " <<it.second;
i++;
}

one_triangle_points.clear();

switch (camera_uv_coords_in_order[index].first)
{
case 0:
{
cam = 0;
break;
}
case 1:
{
cam = 1;
break;
}
case 2:
{
cam = 2;
break;
}
case 3:
{
cam = 3;
break;
}

}

}
qDebug() << "Camera: " << cam << "index" << index << " poi " << res.points[index] << "uv " << closest_p.neighb->point->second[0].second;


camera_separeted_index_lists[cam].push_back(index);

one_triangle_points.push_back(std::pair<unsigned int, Point*>(index, &res.points[index]));
triangle_index++;

}*/

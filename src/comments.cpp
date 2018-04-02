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
*/
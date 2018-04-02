#include "voronoi_diagram.h"

GLPaintFormat VoronoiDiagram::getCellPaintData(const int& ind)
{
	if (size() == 0) return GLPaintFormat();
	//cells[ind].calcConvexHull();
	return cells[ind].getPaintData();
}

GLPaintFormat VoronoiDiagram::getCellWDualPaintData(const int& ind, const int& ind2)
{
	if (size() == 0) return GLPaintFormat();
	//cells[ind].calcConvexHull();
	return cells[ind].getPaintDataWithCellD(ind2);
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

GLPaintFormat VoronoiDiagram::getDelauneySegmentPaintData()
{
	if (indicies_of_delauney_segments.size() == 0) return GLPaintFormat();

	GLPaintFormat res;
	for (auto it : points_with_box_points)
		res.points.push_back(it);

	res.ix = indicies_of_delauney_segments;
	res.col.push_back(QVector3D(1, 0, 0));
	return res;
}

GLPaintFormat VoronoiDiagram::getVoronoiDiagramBySegmentPaintData()
{
	if (point_map.size() == 0) return GLPaintFormat();

	GLPaintFormat res;
	for (auto it : point_map)
		res.points.push_back(it.first);

	res.ix = indicies_of_vornoi_diagram_segments;
	//	qDebug() << "points " << res.points.size();
	//	qDebug() << "indicies_of_vornoi_diagram_segments " << indicies_of_vornoi_diagram_segments.size();
	res.col.push_back(QVector3D(1, 0, 0));
	return res;
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

void VoronoiDiagram::calcDelauneySegments(const delaunay_triangulation& T)
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

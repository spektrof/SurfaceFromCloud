#include "surfacegenerator.h"


GLPaintFormat SurfaceGenerator::getPaintFromPolyhedron()
{
	GLPaintFormat res;

//#ifdef ENABLE_CGAL_SURFACE

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		res.points.push_back(it->point());
	}

	Index index = Index(mesh.vertices_begin(), mesh.vertices_end());

	for (FCI fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi)
	{
		HFCC hc = fi->facet_begin();
		HFCC hc_end = hc;

		do {
			res.ix.push_back(index[VCI(hc->vertex())]);
			++hc;
		} while (hc != hc_end);

	}

//#endif

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
		//CGAL::Nth_of_tuple_property_map<1, IndexedPointWithColorTuple>(),
		nb_neighbors);

	#ifdef ENABLE_DEBUG
		qDebug() << "\tAverage spacing: " << average_spacing;
	#endif
}
#endif

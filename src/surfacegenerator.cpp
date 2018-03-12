#include "surfacegenerator.h"

GLPaintFormat SurfaceGenerator::getPaintFromPolyhedron()
{
	GLPaintFormat res;

#ifdef ENABLE_CGAL_SURFACE

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		res.box_points.push_back(it->point());
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

#endif

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
		qDebug() << "Average spacing: " << average_spacing;
	#endif
}
#endif

#ifdef ENABLE_CGAL_SURFACE

void SurfaceGenerator::PoissonByCGAL()
{
	FT sm_angle = 20.0; // Min triangle angle in degrees.
	FT sm_radius = 30; // Max triangle size w.r.t. point set average spacing.
	FT sm_distance = 0.375; // Surface Approximation error w.r.t. point set average spacing.

	Poisson_reconstruction_function function(points_normals.begin(), points_normals.end(), Point_map(), Normal_map());
	// Computes the Poisson indicator function f()
	// at each vertex of the triangulation.

	#ifdef ENABLE_DEBUG
		qDebug() << "Poisson_reconstruction_function rdy\n";
	#endif

	if (!function.compute_implicit_function())	return;

	// Gets one point inside the implicit surface
	// and computes implicit function bounding sphere radius.
	Point inner_point = function.get_inner_point();

	Sphere bsphere = function.bounding_sphere();
	FT radius = std::sqrt(bsphere.squared_radius());
	// Defines the implicit surface: requires defining a
	// conservative bounding sphere centered at inner point.
	FT sm_sphere_radius = 5.0 * radius;
	FT sm_dichotomy_error = sm_distance * average_spacing / 1000.0; // Dichotomy error must be << sm_distance

	Surface_3 surface(function,
		Sphere(inner_point, sm_sphere_radius*sm_sphere_radius),
		sm_dichotomy_error / sm_sphere_radius);

	// Defines surface mesh generation criteria
	CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
		sm_radius*average_spacing,  // Max triangle size
		sm_distance*average_spacing); // Approximation error
									  // Generates surface mesh with manifold option
	#ifdef ENABLE_DEBUG
		qDebug() << "Rdy to make surface\n";
	#endif

	STr tr; // 3D Delaunay triangulation for surface mesh generation
	C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
	CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
		surface,                              // implicit surface
		criteria,                             // meshing criteria
		CGAL::Manifold_with_boundary_tag());  // require manifold mesh

	if (tr.number_of_vertices() == 0)	return;

	// saves reconstructed surface mesh
	#ifdef ENABLE_DEBUG
		qDebug() << "Rdy to make polyhedron\n";
	#endif

	/*if (!c2t3.is_valid())
	qDebug() << "c2t3 not valid!!\n";
	*/

	mesh.clear();
	CGAL::output_surface_facets_to_polyhedron(c2t3, mesh);

	#ifdef ENABLE_DEBUG
		qDebug() << "Mesh created\n";
	#endif
}

#endif
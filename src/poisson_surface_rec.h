#pragma once

#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>

#include "cgal_types.h"
#include "defines.h"

typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;

static void calc_poisson_reconstruction_by_cgal(const PwN_vector& points_normals, FT& sm_angle, FT& sm_radius, FT& sm_distance, const float& average_spacing, Polyhedron& output_mesh)
{
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

	output_mesh.clear();
	CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);

#ifdef ENABLE_DEBUG
	qDebug() << "Mesh created\n";
#endif
}

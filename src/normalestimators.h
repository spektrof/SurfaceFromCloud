#pragma once
#include "cgal_types.h"	//only for concurency tag mybe refactor

#ifdef ENABLE_CGAL_SURFACE

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

#endif

template <typename T>
class NormalEstimator 
{
public:
	NormalEstimator() {}
	virtual ~NormalEstimator() = 0;

	virtual void normal_c_process(T& points) = 0;

protected:
	void delete_unoriented_normals(T& points, const unsigned int& k)
	{
#ifdef ENABLE_CGAL_SURFACE
		auto unoriented_points_begin =
			CGAL::mst_orient_normals(points.begin(), points.end(),
				Point_map(),
				Normal_map(),
				k);
		// Optional: delete points with an unoriented normal
		// if you plan to call a reconstruction algorithm that expects oriented normals.
		points.erase(unoriented_points_begin, points.end());
#endif
	}
};

template <typename T>
NormalEstimator<T>::~NormalEstimator() {}

#ifdef ENABLE_CGAL_SURFACE

template <typename T>
class PcaNormal : public NormalEstimator<T>
{
public:
	PcaNormal() {}
	~PcaNormal() {}

	void normal_c_process(T& points)
	{
		const unsigned int k = 18;
		CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(), Point_map(), Normal_map(), k);

		delete_unoriented_normals(points, k);
	}
};

#endif

//jetestimator use lapack!!
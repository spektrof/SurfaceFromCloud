#pragma once
#include "defines.h"
#include "types.h"
#include "cgal_types.h"

#include "filters.h"
#include "normalestimators.h"

#include "voronoi_diagram.h"
#include "power_diagram.h"
#include "kdtree.h"
#include "c_io.h"

#include <map>

//#include <boost/tuple/tuple.hpp>
#ifdef ENABLE_CGAL_SURFACE
	#include "poisson_surface_rec.h"

#endif

#ifdef ENABLE_DEBUG
	#include <QDebug>
#endif

class SurfaceGenerator
{

	enum Algorithm
	{
		POISSON,
		POWERCRUST,
		PARTIALDIFF
	};

	struct poisson_component
	{
		/*Parts*/
		NormalEstimator<PwN_vector>* normal_estimator;

		/*Properties*/
		FT sm_angle;			// Min triangle angle in degrees.
		FT sm_radius;			// Max triangle size w.r.t. point set average spacing.
		FT sm_distance;		// Surface Approximation error w.r.t. point set average spacing.

		/*Methods*/
		poisson_component(FT s_a = 20.0, FT s_r = 30, FT s_d = 0.375f) : sm_angle(s_a), sm_radius(s_r), sm_distance(s_d), normal_estimator(NULL)
		{
#ifdef ENABLE_CGAL_SURFACE
			normal_estimator = new PcaNormal<PwN_vector>();	//TODO: adding this as well
#endif
		}

		void SetProperties(FT s_a, FT s_r, FT s_d)
		{
			qDebug() << "New Poisson property: " << s_a << ", " << s_r << ", " << s_d;
			sm_angle = s_a;
			sm_radius = s_r;
			sm_distance = s_d;
		}

		void clear()	//dummy
		{

		}

		void calc_surface(const PwN_vector& points_normals, const float& average_spacing, Polyhedron& mesh)
		{
			calc_poisson_reconstruction_by_cgal(points_normals, sm_angle, sm_radius, sm_distance, average_spacing, mesh);
		}
	};

	struct power_crust_component
	{
		/*Parts*/
		VoronoiDiagram vd;
		PowerDiagram pd;

		int vc_index;
		int pc_index;
		int neighb_index;
		int c_index;

		/*Propeties*/

		/*Methods*/
		power_crust_component(const int& v_i = 0, const int& p_i = 0, const int& n_i = 0, const int& c_i = 0) : vc_index(v_i), pc_index(p_i), neighb_index(n_i), c_index(c_i) { }

		void SetProperties(const float& a, const float& b, const float& c)
		{
			qDebug() << "New PowerCrust property: " << a << ", " << b << ", " << c;
		}

		void clear()
		{
			pd.clear();
			vd.clear();
		}

		void calc_surface(std::vector<Point>& points, Box minimal_bound_box)
		{
			vd.setBox(minimal_bound_box * 2);	//TODO 5 ki configba
		//	vd.calc_diagram(points);
			vd.calc_diagram_with_poles(points);

			pd.setBox(minimal_bound_box);	//TODO 5 ki configba
			pd.calc_diagram(vd.getPoles());
			pd.label_poles();
			pd.calc_power_crust();
		}
	};

public:
	SurfaceGenerator() :
		alg(POWERCRUST), draw_object(POINTSS), draw_pc_object(RESULT)
		//,root(NULL)
	{
		//TODO: file path from config
		cloud.setFileName("D:/ELTE/Diplomamunka/sfc/sfc/clouds/sanyi/00008_127.0.0.1.ply");

		filters.clear();

		last_simplif = last_outlier = last_plane_outliner = 0;

		all_calculation_time = milliseconds(0);
		surface_calculation_time = milliseconds(0);
		cloud_load_time = milliseconds(0);
	}

	~SurfaceGenerator()
	{
		//std::for_each(filters.begin(), filters.end(), [](Filter<Point>* fil) { delete fil; });	//for simplification / smooth algorithms

		//active_points.clear();
		//active_nodes.clear();
	//	delete cloud;
	}

	/*************************************
	*	Process and Paint helper functions	*
	**************************************/

	ProcessReturns process()
	{
		time_p start = hrclock::now();

		// STEP 0. : cloud import

		switch (cloud.load_process(points))
		{
		case CLOUD_LOAD_SUCCESS:
		{
			//points_o.clear();
			contain_uv = false;
			if (points.size() == 0) return NO_INPUT_DATA;

			contain_uv = points[0].second[0].first != -1;

			/*int ct = 0;

			for (auto& data : points)
			{
				if (data.second.size() != 1)
				{
					qDebug() << "\t siz " << data.second.size();
					ct++;
				}
				points_o.push_back(data.first);

				if (data.second[0].first != -1) contain_uv = true;
			}
			qDebug() << "There are " << ct << " points with more than 1 color";*/

			qDebug() << "does it contain uv coords ?" << contain_uv;
			break;
		}
		case NO_INPUT_DATA:
			return NO_INPUT_DATA;
		default:
			return UNDEFINED_ERROR;
		}

		cloud_load_time = boost::chrono::duration_cast<milliseconds>(hrclock::now() - start);

		//STEP 1. : remove outliers: a) plane based, b) average space based

		auto filter_part = last_plane_outliner;
		auto it = 0;

		//a)
		for (; it < filter_part; ++it)
		{
			filters[it]->filter_process(points);
		}

		qDebug() << "After Plane Outlier removals : " << points.size() << " point remain";
#ifdef ENABLE_CGAL_FILTER
		//b)
		filter_part += last_outlier;
		compute_average_spacing();

		for (; it < filter_part; ++it)
		{
			filters[it]->SetThirdProperty(average_spacing);
			filters[it]->filter_process(points);
			compute_average_spacing();
		}
#endif
		qDebug() << "After Outlier removals : " << points.size() << " point remain";

		//STEP 2. : USE MINE KDTREE
		//REFACTOR -> put into filter
		/*switch (cloud->build_process(points))
		{
			case PROCESS_DONE:
			{
				points.clear();
				cloud->getActivePoints(points);
				break;
			}
			case UNDEFINED_ERROR:
			default:
			{

			#ifdef ENABLE_DEBUG
				qDebug() << "UNDEFINED ERROR at MINE KDTREE\n";
			#endif
			}
		}*/

		//STEP 3. : do pointset simplification and smooth
#ifdef ENABLE_CGAL_FILTER
		bool changed = false;

		//simplifiers
		filter_part += last_simplif;
		for (; it < filter_part; ++it)
		{
			filters[it]->filter_process(points);
			changed = true;
		}

		qDebug() << "After Simplifiers : " << points.size() << " point remain";

		if (changed) compute_average_spacing();

		//smoothers
		for (; it < filters.size(); ++it)
		{
			filters[it]->filter_process(points);
		}

#endif
		qDebug() << "\n------------------------------\nfiltering done : " << points.size() << "\n";

		all_calculation_time = boost::chrono::duration_cast<milliseconds>(hrclock::now() - start);
		/*kk.clear();
		for (auto& it : points)
		{
			Point a = it.first;
			auto& tt = std::find_if(points_o.begin(), points_o.end(), [a](const Point& p) { return a == p; });
			if (tt != points_o.end()) kk.push_back(it);
		}*/

		return PROCESS_DONE;
	}

	ProcessReturns surface_process()
	{
		//points_normals.clear();
		if (points.empty()) return EMPTY_POINTS_SET;

		time_p start = hrclock::now();

		switch (alg)
		{
		case POISSON:
		{
			PwN_vector points_normals;		//pontok normálissal

			qDebug() << points.size();
#ifdef ENABLE_CGAL_SURFACE
			for (auto& it : points)
				points_normals.push_back(Point_with_Normal(std::make_pair(it.first, Vector())));

			if (poisson.normal_estimator == NULL) return NOT_ENABLED_COMPONENT;

			poisson.normal_estimator->normal_c_process(points_normals);
			qDebug() << points_normals.size();

			poisson.calc_surface(points_normals, average_spacing, mesh);

			surface_calculation_time = boost::chrono::duration_cast<milliseconds>(hrclock::now() - start);

			qDebug() << surface_calculation_time.count();
			return PROCESS_DONE;
#else
			return NOT_ENABLED_COMPONENT;
#endif
		}
		case POWERCRUST:
		{
			std::vector<Point> points_s;
			for (auto& it : points)
			{
				points_s.push_back(it.first);
			}

			power_crust.calc_surface(points_s, cloud.getBox());

			return PROCESS_DONE;
		}
		case PARTIALDIFF:
			return NOT_ENABLED_COMPONENT;
		default:
			return UNDEFINED_ERROR;
		}
	}

	GLPaintFormat getPaintData()
	{
		switch (draw_object)
		{
		case POINTSS:
		{
			std::vector<Point> points_s;
			for (auto& it : points)
			{
				points_s.push_back(it.first);
			}
			return GLPaintFormat(points_s);
		}
		case BOX:
			//	return cloud->drawNode();
			return GLPaintFormat();
		case SURFACE:
		{
			GLPaintFormat res = getPaintFromPolyhedron(points);
			add_tex_or_color_info(res);
			return res;
		}
		case ALGORITHM:
		{
			switch (draw_pc_object)
			{
			case RESULT:
			{
				GLPaintFormat res = power_crust.pd.get_triangled_mesh_without_texture();
				if (contain_uv) add_tex_or_color_info(res);
				else
				{
					for (size_t i = 0; i < res.points.size(); ++i)
						res.points_col.push_back(QVector3D(0, 0, 0));

					res.point_part_lengths.push_back(res.ix.size());
				}

				return res;
			}
			case DELAUNEY:
			{
				GLPaintFormat res;
				for (auto& it : points)
					res.points.push_back(it.first);

				power_crust.vd.getDelauneySegmentPaintData(res);
				//	res.points = points;
				return res;
			}
			case VORONOI_WITH_SURF_POINT:
			case VORONOI_BY_CELL_FULL:
				return power_crust.vd.getCellPaintData(power_crust.vc_index);
			case POWER_DIAGRAM_BY_CELL:
				return power_crust.pd.getCellPaintData(power_crust.pc_index, -1);
			case VORONOI_DIAGRAM:
				return power_crust.vd.getVoronoiDiagramBySegmentPaintData();
			case POWER_DIAGRAM:
				return power_crust.pd.getPowerDiagramBySegmentPaintData();
			case POLES:
				return power_crust.vd.getPolesSurfPaintData();
			case POWER_DIAGRAM_CELL_WITH_NEIGHBOURS:
				return power_crust.pd.getCellPaintData(power_crust.pc_index, power_crust.neighb_index);
			case INNER_POLES:
				return power_crust.pd.getInnerPolesPaintData();
			case OUTER_POLES:
				return power_crust.pd.getOuterPolesPaintData();
			case UNKNOWN_POLES:
				return power_crust.pd.getUnknownPolesPaintData();
			case INNER_OUTER_POLES:
				return power_crust.pd.getInnerOuterPairsPaintData();
			case VORONOI_WITH_CELLDUAL:
				return power_crust.vd.getCellWDualPaintData(power_crust.vc_index, power_crust.c_index);
			case MEDIAL_AXIS:
				return power_crust.pd.getPowerCrustPaintData();
			case POWER_SHAPE:
				return power_crust.pd.getPowerShapePaintData();
			default:
				return GLPaintFormat();
			}
		}
		default:
			return GLPaintFormat();
		}
	}

	/***********************
	*	New Objectums      *
	***********************/

	bool newFileInput(char* filename, const int& ma, bool repeat)
	{
		for (auto& it : filters) delete it;
		position_link.clear();
		filters.clear();
		last_outlier = last_plane_outliner = last_simplif = 0;

		points.clear();
#ifdef	ENABLE_CGAL_SURFACE
		mesh.clear();
#endif
		all_calculation_time = milliseconds(0);
		surface_calculation_time = milliseconds(0);
		cloud_load_time = milliseconds(0);

		cloud.setCloudSource(1);
		cloud.setFileName(filename);

		return cloud.checkFileType();
	}

	bool newRandomInput(const unsigned int& n, const int& ma, bool repeat)
	{
		for (auto& it : filters) delete it;
		position_link.clear();
		filters.clear();
		last_outlier = last_plane_outliner = last_simplif = 0;

		points.clear();
#ifdef ENABLE_CGAL_SURFACE
		mesh.clear();
#endif
		all_calculation_time = milliseconds(0);
		surface_calculation_time = milliseconds(0);
		cloud_load_time = milliseconds(0);

		cloud.setCloudSource(0);
		cloud.setNumberOfPoint(n);

		return true;
	}
	//-------------------------

	bool addNewOutlierFilter(const unsigned int& type, const float& kn_px, const float& limit_py = 0.0f, const float& pz = 0.0f, const float& nx = 0.0f, const float& ny = 0.0f, const float& nz = 0.0f)
	{
		qDebug() << "NEW OUTLIER " << type;

		switch (type)
		{
		case 1:
		{
#ifdef ENABLE_CGAL_FILTER
			int index = last_plane_outliner + last_outlier;

			filters.insert(filters.begin() + index, new OutlierRemoval<point_uv_tuple>(kn_px, limit_py));
			last_outlier++;

			//qDebug() << " Before \n";
			//for (auto it : position_link) qDebug() << it << "\n";

			for (auto& it : position_link) { if (it >= index) it++; }
			position_link.push_back(index);
			//qDebug() << " After \n";
			//for (auto it : position_link) qDebug() << it << "\n";

			return true;
#else
			return false;
#endif
		}
		case 2:
		{
			int index = last_plane_outliner;

			filters.insert(filters.begin() + index,
				new OutlierComponentRemoval<point_uv_tuple>(Vector(kn_px, limit_py, pz), Vector(nx, ny, nz)));
			last_plane_outliner++;

			//	qDebug() << " Before \n";
			//	for (auto it : position_link) qDebug() << it << "\n";

			for (auto& it : position_link) { if (it >= index) it++; }
			position_link.push_back(index);
			//	qDebug() << " After \n";
			//	for (auto it : position_link) qDebug() << it << "\n";

			return true;
		}
		}
		return false;
	}

	bool addNewSimplifierFilter(const unsigned int& type, const float& first, const float& second = 0.0f)
	{
		qDebug() << "NEW SIMPLIF " << type << " " << first <<" " << second;
		switch (type)
		{
		case 1:
		{
		#ifdef ENABLE_CGAL_FILTER
			int index = last_plane_outliner + last_outlier + last_simplif;

			filters.insert(filters.begin() + index, new GridSimplification<point_uv_tuple>(first));
			last_simplif++;

		//	qDebug() << " Before \n";
		//	for (auto it : position_link) qDebug() << it << "\n";

			for (auto& it : position_link) { if (it >= index) it++; }
			position_link.push_back(index);
		//	qDebug() << " After \n";
		//	for (auto it : position_link) qDebug() << it << "\n";

			return true;
		#else
			return false;
		#endif
		}
		case 2:
		{
			//TODO MINE
			//filters.push_back(new OutlierComponentRemoval<Point>());
			return false;
		}
		case 3:
		{
#ifdef ENABLE_CGAL_FILTER
			int index = last_plane_outliner + last_outlier + last_simplif;

			filters.insert(filters.begin() + index, new HiearchySimplification<point_uv_tuple>(first, second));
			last_simplif++;

			//qDebug() << " Before \n";
		//	for (auto it : position_link) qDebug() << it << "\n";

			for (auto& it : position_link) { if (it >= index) it++; }
			position_link.push_back(index);
		//	qDebug() << " After \n";
		//	for (auto it : position_link) qDebug() << it << "\n";
			return true;
#else
			return false;
#endif
		}
		case 4:
		{
			int index = last_plane_outliner + last_outlier + last_simplif;

			filters.insert(filters.begin() + index, new WLOPSimplification<point_uv_tuple>(first, second));
			last_simplif++;

			qDebug() << " Before \n";
				for (auto it : position_link) qDebug() << it << "\n";

			for (auto& it : position_link) { if (it >= index) it++; }
			position_link.push_back(index);
				qDebug() << " After \n";
				for (auto it : position_link) qDebug() << it << "\n";
			return true;
		}
	}
/*
#ifdef ENABLE_CGAL_SURFACE
		filters.push_back(new WLOPSimplification<Point>());
#endif*/
		return false;
	}

	bool addNewSmootherFilter(const unsigned int& type, const float& first, const float& second = 0.0f)
	{
		qDebug() << "NEW SMOOTHER " << type << " " << first << " " << second;

		switch (type)
		{
			case 1:
			{
			#ifdef ENABLE_CGAL_FILTER
				#ifdef CGAL_EIGEN3_ENABLED
	
					filters.push_back(new JetSmoothing<point_uv_tuple>(first));

					position_link.push_back(filters.size()-1);
				//	for (auto it : position_link) qDebug() << it << "\n";

					return true;
				#else
					return false;
				#endif
			#else	
				return false;
			#endif
			}
		}
		return false;
	}

	bool addNewSurfaceReconstructor(const unsigned int& type, const float& first, const float& second = 0.0f, const float& third = 0.0f)
	{
		power_crust.clear();
		poisson.clear();

		switch (type)
		{
		case 0:
		{
			qDebug() << "Surface rec algorithm is Poission";
		#ifdef ENABLE_CGAL_SURFACE

			alg = Algorithm(type);
			poisson = poisson_component(first, second, third);
			return true;
		#else
			return false;
		#endif
		}
		case 1:
		{
			qDebug() << "Surface rec algorithm is Power Crust";
			alg = Algorithm(type);
			power_crust = power_crust_component();
			return true;
		}
		default:
			return false;
		}

		return false;
	}

//-----------------------------------------------------------------------------

	void SetFilterProperty(const unsigned int& id, const double& first, const double& second = 0.0f, const double& third = 0.0f)
	{
		filters[position_link[id]]->SetFirstProperty(first);
		filters[position_link[id]]->SetSecondProperty(second);
		filters[position_link[id]]->SetThirdProperty(third);
	}

	void SetFilterProperty(const unsigned int& id, const double& x, const double& y, const double& z, const double& nx, const double& ny, const double& nz)
	{
		filters[position_link[id]]->SetFirstProperty(Vector(x,y,z));
		filters[position_link[id]]->SetSecondProperty(Vector(nx,ny,nz));
	}

	void SetSurfaceProperty(const unsigned int& id, const float& first, const float& second = 0.0f, const float& third = 0.0f)
	{
		switch (id)
		{
		case 0:
		{
			poisson.SetProperties(first, second, third);
			break;
		}
		case 1:
		{
			power_crust.SetProperties(first, second, third);
			break;
		}
		}
	}

	//TODO: replace to shared_ptr
	void deleteFilter(const unsigned int& id)
	{
		int rm = position_link[id];
		int index = last_plane_outliner;

		delete filters[rm];
		filters.erase(filters.begin() + rm, filters.begin() + rm + 1);

		position_link.erase(position_link.begin() + id, position_link.begin() + id + 1);

		for (auto& it : position_link) { if (it > rm) it--; }
	//	for (auto it : position_link) qDebug() << it;

		if (index > rm)
		{
			last_plane_outliner--;
			qDebug() << "Plane outlier remov. deleted\n";
			return;
		}

		index += last_outlier;

		if (index > rm)
		{	
			last_outlier--;
			qDebug() << "Outlier remov. deleted\n";
			return;
		}

		index += last_simplif;

		if (index > rm)
		{	
			last_simplif--;
			qDebug() << "Simplifier deleted\n";
			return;
		}

		qDebug() << "Smoother deleted\n";
		return;
		
	}

	void export_triangle_mesh_with_tex_to_obj()
	{

	}

	void export_triangle_mesh_to_obj()
	{
		GLPaintFormat res = power_crust.pd.get_triangled_mesh_without_texture();
		cloud.export_to_obj(res.points, res.ix);
	}
//-----------------------------------------------------------------------------
	//Setters & Getters
	void setDrawObject(const unsigned int& dp) { draw_object = DrawPossibilites(dp); /*qDebug() << draw_object;*/ }
	DrawPossibilites getDrawObject() const { return draw_object; }
	void setPowerCrustDrawObject(const unsigned int& dp) { draw_pc_object = PowerCrustDrawPossibilites(dp); /*qDebug() << draw_pc_object;*/
	}
	PowerCrustDrawPossibilites getPowerCrustDrawObject() const { return draw_pc_object; }

	//-------------------------
	// Routers for cloud
	void swapStartFlag() { cloud.swapStartFlag(); }
	void setMeshAccuracy(const int& val) { /*cloud->setMeshAccuracy(val);*/ }
	void resetCloud() { /*cloud->resetCloud();*/ }
	void goRightNode() { /*cloud->goRightNode(); */}
	void goLeftNode() {  /*cloud->goLeftNode(); */}
	void goNextCell() {
		if (power_crust.vc_index < power_crust.vd.size() - 1) power_crust.vc_index++;
		if (power_crust.pc_index < power_crust.pd.size() - 1) power_crust.pc_index++;
	}
	void goPreviousCell() {
		if (power_crust.vc_index > 0) power_crust.vc_index--;
		if (power_crust.pc_index > 0) power_crust.pc_index--;
	}

	void nextNeighb() { power_crust.neighb_index++; }
	void nextCInd() {  power_crust.c_index++; }

	void goParentNode() { /*cloud->goParentNode();*/ }
	int getPointNumbFromBox() const { return 0;/*cloud.getPointNumbFromBox();*/ }
	void getTimes(float& c_load_time, float& s_cal_time, float& all_time, float& join_time, float& sort_time) {
		all_time = all_calculation_time.count();
		s_cal_time = surface_calculation_time.count();
		c_load_time = cloud_load_time.count();
		//qDebug() << s_cal_time;

		float dummy;
	//	cloud->getTimes(join_time, sort_time, dummy); 
	}
	void changeThreadCapacity(const unsigned int& cloud_t) { /*cloud->changeThreadCapacity(cloud_t);*/ }
	unsigned int getThreadCapacity() const { return 2; /*return cloud->getThreadCapacity();*/ }

	std::vector<std::string> get_file_name_list() { return cloud.get_file_name_list(); }

	/************************************************************************
	*																		*
	*************************************************************************/
protected:
/*
	6. Saját simplifierek, saját architektúrában
*/
	//-----------------------------------------
	/***********************
	*	Filter componenet   *
	***********************/
#ifdef ENABLE_CGAL_FILTER
	void compute_average_spacing();
#endif
	//-----------------------------------------
	/***********************
	*	Surface componenet   *
	***********************/
	Polyhedron mesh;
	//-----------------------------------------

	GLPaintFormat getPaintFromPolyhedron(point_uv_map&);
	void add_tex_or_color_info(GLPaintFormat& res);

	typedef std::pair<unsigned int, Point*> index_pointptr;
	typedef std::pair<unsigned int, Point> index_point;
	//					camera						point index		uv
	typedef std::pair<int16_t, std::vector<std::pair<index_point, std::pair<float, float>>>> texture_tuple;
	typedef std::pair<Point, std::pair<float, float>> point_uv;
	typedef std::vector<std::pair<unsigned int, Point*>> index_pointptr_v;
	typedef std::pair<int, int> count_index_link;
	typedef std::vector<std::vector<unsigned int>> index_by_cameras_v;
	typedef std::vector<std::pair<index_point, uv_map>> pointuvs_v;
						//point			cam			index for multiple textured points
	typedef std::map<Point, std::map<unsigned int, unsigned int>> p_cam_ind;

	struct uv_list
	{
		UV cur;
		std::vector<unsigned int> other_inds;
		uv_list(UV& c = UV(-1, std::pair<float, float>(0.0f, 0.0f))) : cur(c) { other_inds.clear(); }
	};
	typedef std::vector<uv_list> cam_coords;


	void calc_triangle_details(kdTree&, std::vector<index_point>&, cam_coords&, std::vector<uv_map>&,	std::vector<std::vector<unsigned int>>&,  std::vector<Point>&, int&);

	std::pair<texture_tuple, bool> calculate_uv_coords_for_triangle( kdTree& tree, std::vector<index_point>&, cam_coords&, std::vector<uv_map>& uv_cache,  int& tt);
	
	uv_map calc_new_uv_coords_from_knearest(kdTree::kNearest_queue&, Point&);
	std::pair<float, float> calc_new_uv_from_closest_three(std::vector<point_uv>& closest_three,Point triangle_point);

	std::pair<texture_tuple, bool> calc_camera_based_uv(const pointuvs_v& paired_point_uvs, cam_coords& camera_uv_coords_in_order, int& tt);
	std::set<count_index_link, std::greater<count_index_link>> get_tex_file_chances(std::vector<index_point>&, cam_coords&);
	void index_correction(std::vector<std::pair<index_point, std::pair<float, float>>>& related_uv_map, cam_coords& camera_uv_coords_in_order, const int& camera);

	void fill_uv_list(cam_coords& camera_separeted_uv_coords, const texture_tuple&texture_uv_result, std::vector<Point>&);
	void fill_index_list(index_by_cameras_v& camera_separeted_index_lists, const int16_t& camera_index, const std::vector<index_point>& );
	void update_point_set(std::vector<Point>& points, cam_coords& camera_uv_coords_in_order, const index_pointptr_v&);

private:
	/***********************
	*	Private componenet   *
	***********************/

	Algorithm alg;

	power_crust_component power_crust;
	poisson_component poisson;

	DrawPossibilites draw_object;
	PowerCrustDrawPossibilites draw_pc_object;

	C_IO<Point, UV> cloud;

	std::vector< Filter<point_uv_tuple>* > filters;	//for simplification / smooth algorithms
	std::vector<int> position_link;			//position link between the UI and this

	unsigned int last_plane_outliner;
	unsigned int last_outlier;
	unsigned int last_simplif;

	// [Outlier, ... , Simplif, Simplif ,... ,Simplif, Smooth, ..., Smooth]
	
	double average_spacing;		//TODO for filtering... mybe dont need here

	point_uv_map points;		//pontok
	//std::vector<Point> points_o;		//pontok_2
	//point_uv_map kk;

	bool contain_uv;

	//-----------------------------
	milliseconds all_calculation_time;
	milliseconds surface_calculation_time;
	milliseconds cloud_load_time;
};

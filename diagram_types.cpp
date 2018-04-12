#include "diagram_types.h"
#include <queue>

void face_t::correct_edge_orientations()
{
	if (edges.empty()) return;

	std::map<face_t*, unsigned int>::iterator location;

	location = edges[0]->related_faces.find(this);

	std::vector<edge_tt*> edges_o;
	edges_o.push_back(edges[0]);
	std::vector<edge_tt*> edges_cpy(edges.begin() + 1, edges.end());

	//todo: mivan ha rendezve tárolnánk õket vhogy...
	for (size_t i = 0; i < edges.size() -1 ; ++i)
	{
		Point* actual_last_point = location->second == 0 ? edges_o[i]->edge.second.second : edges_o[i]->edge.first.second;

		unsigned int new_orientation;
		auto next = std::find_if(edges_cpy.begin(), edges_cpy.end(), [&actual_last_point, &new_orientation](const edge_tt* e)->bool
		{
			if (e->edge.first.second == actual_last_point)
			{
				new_orientation = 0;
				actual_last_point = e->edge.second.second;
				return true;
			}
			else if (e->edge.second.second == actual_last_point)
			{
				new_orientation = 1;
				actual_last_point = e->edge.first.second;
				return true;
			}
			return false;
		});

		if (next == edges_cpy.end())
		{
			qDebug() << "FAIL not founded pair";
			return;
		}

		location = (*next)->related_faces.find(this);
		location->second = new_orientation;

		edges_o.push_back(*next);
		edges_cpy.erase(next);
	}

	std::vector<edge_tt*>(edges_o).swap(edges);

	//-------------------- TEST

	/*for (auto& it : edges)
	{
		if (it->related_faces[this] == 0)
			qDebug() << it->edge.first.first << " -> " << it->edge.second.first;
		else
			qDebug() << it->edge.second.first << " -> " << it->edge.first.first;
	}*/
}

void face_t::correct_boundary_orientation()
{
	qDebug() << "Correct Boundary orientation";
	std::queue<face_t*> face_list;
	std::set<face_t*> progress_faces;
	std::set<face_t*> done_faces;

	progress_faces.insert(this);
	done_faces.insert(this);

	face_list.push(this);

	while (!face_list.empty())
	{
		face_t* act_face = face_list.front();
		face_list.pop();

		auto act_face_edges = act_face->edges;

		for (auto& face_edge : act_face_edges)
		{
			face_t* next_boundary_face = nullptr;
			int orientation_on_actual_face = -1;
			int next_face_orientation = -1;
			for (auto& related_face : face_edge->related_faces)
			{
				if (!related_face.first->is_power_crust_face()) continue;
				if (related_face.first == act_face)
				{
					orientation_on_actual_face = related_face.second;
					continue;
				}

				if (next_boundary_face == nullptr)
				{
					next_boundary_face = related_face.first;
					next_face_orientation = related_face.second;
				}
				else
					qDebug() << "CORRECT ORIENTATION ON BOUND IS BAAD";

			}

			if (next_boundary_face != nullptr)
			{
				if (done_faces.find(next_boundary_face) != done_faces.end())
				{
					if (next_face_orientation == orientation_on_actual_face)
					{
						qDebug() << "ERR: Bad orientation!";
					}
					continue;
				}

				if (next_face_orientation < 0 || orientation_on_actual_face < 0)
				{
					qDebug() << "ERR: orientation didnt set properly";
					continue;
				}

				if ( next_face_orientation == orientation_on_actual_face )
					next_boundary_face->swap_edge_orientations();

				face_list.push(next_boundary_face);
				done_faces.insert(next_boundary_face);
			}
		}

	}
	qDebug() << "Correct Boundary orientation END";

}
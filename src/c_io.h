#pragma once
//#include "cgal_types.h"
#include "types.h"
#include "defines.h"
#include <vector>
#include <fstream>
#include <string.h>

#include <cstdlib>
#include <stdio.h>

//#include <set>
#include <map>
//#include <boost\thread\mutex.hpp>
//#include <boost\thread\thread.hpp>

#include <boost\filesystem\fstream.hpp>

#include <QDebug>

					//point			texture		uv coord -> not exist third!
template <typename P, typename T, bool UV = true>
class C_IO
{
	typedef std::pair<P, std::vector<T>> tuple;

	enum cloud_sources
	{
		RANDOM,
		FILE
	};

	enum file_extension
	{
		PLY,
		OBJ,
		UNDEFINED
	};

public:
	typedef typename std::vector<tuple> point_t_map;

	C_IO() : box(Box(-0.1f, 0.1f)), n_point(0), file_index(1), cloud_s(FILE), ext(UNDEFINED)
	{
		loadCalibrations();
		started = false;
		run_count = 0;
	}

	~C_IO() {}

	void setNumberOfPoint(const unsigned int& pn) { n_point = pn; }
	unsigned int getNumberOfPoint() const { return n_point; }

	void setBox(const Section& s) { box = Box(s); }
	void setBox(const float& min, const float& max) { box = Box(min, max); }
	void setBox(const Section& _x, const Section& _y, const Section& _z) { box = Box(_x, _y, _z); }
	void setBox(const float& xmi, const float& xma, const float& ymi, const float& yma, const float& zmi, const float& zma) { box = Box(xmi, xma, ymi, yma, zmi, zma); }

	Box* getBoxObject() { return &box; }
	Box getBox() const { return box; }

	void swapStartFlag() { started = 1 - started; run_count = 0; }

	//---------------------------
	void set_ply(char* fn)
	{
		ext = PLY;

		file_name = strrchr(file_path, '/');
		file_name++;

		adress_c = strrchr(file_name, '_');
		adress_c++;

		memset(counter_prefix, 0, sizeof(counter_prefix));
		strncpy(counter_prefix, file_name, adress_c - file_name - 1);

		file_index = file_index_start = atoi(counter_prefix);

		qDebug() << file_name << " " << adress_c << " " << file_index;
	}

	void set_obj(char* fn)
	{
		ext = OBJ;

		file_name = strrchr(file_path, '/');
		file_name++;

		qDebug() << file_name << " " << file_path;
	}

	void setFileName(char* fn)
	{
		strcpy(file_path, fn);
		postfix = strrchr(fn,'.');
		qDebug() << postfix;

		if (strcmp(postfix,".ply") == 0)
		{
			qDebug() << "ply load";
			set_ply(fn);
		}
		else if (strcmp(postfix, ".obj") == 0)
		{
			qDebug() << "obj load" <<fn << " " << file_path;
			qDebug() << postfix;
			set_obj(fn);
			qDebug() << postfix;

		}
		else
		{
			ext = UNDEFINED;
			qDebug() << "undefined";
		}
	
	}

	char* assemble_ply_filename(const unsigned int& i)
	{
		char c_file_index[50];
		itoa(file_index, c_file_index, 10);
		char* file_index_c = c_file_index;

		strncpy(counter_prefix + strlen(counter_prefix) - strlen(c_file_index), c_file_index, strlen(file_index_c));
		strncpy(file_path + (file_name - file_path), counter_prefix, strlen(counter_prefix));
		strncpy(file_path + (adress_c - file_path), ips[i], strlen(ips[i]));

		char tmp[100];
		memset(tmp, 0, sizeof(tmp));
		strncpy(tmp, file_path, adress_c - file_path + strlen(ips[i]));
		strcat(tmp, ".ply");
		strcpy(file_path, tmp);
		return file_path;
	}

	char* assemble_obj_filename()
	{
		return file_path;
	}

	void resetFileIndex() { file_index = file_index_start; }

	bool checkFileType() 
	{
		bool existOne = false;
		switch (ext)
		{
			case PLY:
			{
				for (int i = 0; i < kinectSources; ++i)
				{
					std::ifstream file(getFile(i), std::ios::binary);
					if (!file.fail()) existOne = true;
				}
				break;
			}
			case OBJ:
			{
				std::ifstream file(getFile(), std::ios::binary);
				if (!file.fail()) existOne = true;
				break;
			}
			default:
				return false;
		}
		
		return existOne;
	}

	bool nextFileExist()
	{
		return checkFileType();
	}

	void setCloudSource(const unsigned int& type)
	{
		started = false;
		run_count = 0;
		cloud_s = cloud_sources(type);
	}

	ProcessReturns load_process(point_t_map& points)
	{
		if (!started) return PROCESS_STOPPED;

		switch (cloud_s)
		{
			case RANDOM:
			{
				if (run_count > 0 /*&& !repeatable*/) return NO_INPUT_DATA;
				run_count = 1;

				points.clear();
				points = getCloudRandom();
				return CLOUD_LOAD_SUCCESS;
			}
			case FILE:
			{
				if (!nextFileExist() || run_count > 0)
				{
					//qDebug() << "Next FIle doesnt exist!";
					resetFileIndex();
					/*if (!repeatable)*/ return NO_INPUT_DATA;
				}
				run_count = 1;

				points.clear();
				points = getCloudFromFile();

				return CLOUD_LOAD_SUCCESS;
			}
			default:
				return UNDEFINED_ERROR;
		}
	}

	std::vector<std::string> get_file_name_list();

	void getCloudFromPly(std::set<tuple>& res)
	{
		//boost::filesystem::ifstream files[kinectSources];
		time_p t1 = hrclock::now();

		for (int i = 0; i<kinectSources; ++i)
		{
			boost::filesystem::ifstream file(getFile(i));

			if (file.fail()) continue;

			if (!parseHeader(file, i)) { file.close(); continue; }

			boost::filesystem::ifstream tex_file(getTexFile(i));

			//thread_ptrs.push_back(thread_pool.create_thread( boost::bind(&C_IO::parseContext, this, &files[i], &tmp, i)));
			//if (!parseContext(file,tmp)) { file.close(); continue; }
			parseContext(&file, &tex_file, &res, i);
			file.close();
			tex_file.close();
		}

		/*thread_pool.join_all();
		for (auto it : thread_ptrs)
		thread_pool.remove_thread(it);
		thread_ptrs.clear();
		*/

		qDebug() << "LOAD under " << (boost::chrono::duration_cast<milliseconds>(hrclock::now() - t1)).count() << " ms with " << res.size() << " points";

		file_index++;
		qDebug() << "\tOUR BOX IS: " << box.getXMin() << " " << box.getXMax() << " " << box.getYMin() << " " << box.getYMax() << " " << box.getZMin() << " " << box.getZMax();
	}

	void getCloudFromObj(std::set<tuple>& points)
	{
		qDebug() << "getting cloud from obj";

		time_p t1 = hrclock::now();

		boost::filesystem::ifstream file(getFile());
		//boost::filesystem::ifstream tex(NULL);

		if (file.fail()) return;

		float minX = FLT_MAX, maxX = FLT_MIN;
		float minY = FLT_MAX, maxY = FLT_MIN;
		float minZ = FLT_MAX, maxZ = FLT_MIN;

		std::string line;
		float x, y, z;

		while (std::getline(file, line))
		{
			std::istringstream ls(line);
			std::string token;
			ls >> token;

			if (token != "v") continue;
			
			ls >> x >> y >> z;
		//	qDebug() << x << " " << y << z;

			auto res = points.insert(tuple(P(x, y, z), std::vector<T>()));
			std::vector<T>* non_const = const_cast<std::vector<T>*>(&res.first->second);

			if (UV)
			{
				if (false)
				{
					float u, v;
				//	tex >> u >> v;
				//	non_const->push_back(T(0, std::pair<float, float>(u, v)));
				}
				else
					non_const->push_back(T(-1, std::pair<float, float>(0.0f, 0.0f)));
			}

			minX = minX > x ? x : minX;
			maxX = maxX < x ? x : maxX;
			minY = minY > y ? y : minY;
			maxY = maxY < y ? y : maxY;
			minZ = minZ > z ? z : minZ;
			maxZ = maxZ < z ? z : maxZ;

		}
	
		file.close();
	//	if (!tex) tex.close();

		qDebug() << "LOAD under " << (boost::chrono::duration_cast<milliseconds>(hrclock::now() - t1)).count() << " ms with " << points.size() << " points";

		if (!box.isInside(minX, minY, minZ) || !box.isInside(maxX, maxY, maxZ))
		{
			box = Box(minX < box.getXMin() ? minX : box.getXMin()
				, maxX > box.getXMax() ? maxX : box.getXMax()
				, minY < box.getYMin() ? minY : box.getYMin()
				, maxY > box.getYMax() ? maxY : box.getYMax()
				, minZ < box.getZMin() ? minZ : box.getZMin()
				, maxZ > box.getZMax() ? maxZ : box.getZMax());
		}

		qDebug() << "\tOUR BOX IS: " << box.getXMin() << " " << box.getXMax() << " " << box.getYMin() << " " << box.getYMax() << " " << box.getZMin() << " " << box.getZMax();
	}

	void export_to_obj(std::vector<P>& points, std::vector<unsigned int>& indicies)
	{
		qDebug() << "in export" << file_path << " " << file_name;
		char out_file[100];
		memset(out_file, 0, sizeof(out_file));
		strcpy(out_file, "result/");
		qDebug() << out_file;

		char out_filename[50];
		memset(out_filename, 0, sizeof(out_filename));

		postfix = strrchr(file_path, '.');
		qDebug() << postfix;
		qDebug() << strlen(file_name) - strlen(postfix);
		strncpy(out_filename, file_name, strlen(file_name) - strlen(postfix));
		qDebug() << out_filename;
		strcat(out_filename, "_mesh");
		qDebug() << out_filename;

		strcat(out_file, out_filename);
		qDebug() << "OutFile is " << out_file;
		strcat(out_file, ".obj");
		qDebug() << "OutFile is " << out_file;

		boost::filesystem::ofstream file(out_file);
		qDebug() << "ok";
		file << "#" << out_filename << " with PowerCrust by Peter Lukacs\n\n";
		for (auto& poi : points)
			file << "v " << std::setprecision(15) << poi.x() << " " << poi.y() << " " << poi.z() <<"\n";

		file << "#" << points.size() << " vertices\n";

		for (size_t i = 0; i < indicies.size(); i+=3)
			file << "f " << indicies[i] << " " << indicies[i+1] << " " << indicies[i+2] << "\n";
		file << "#" << indicies.size()/3 << " faces\n";
		file.close();
	}

protected:
	void loadCalibrations();
    char* getFile(const int& i = 0);
	char* getTexFile(const int& i);
	char* getPngFile(const int& i);

	void read_header_element(std::istream& in, const unsigned int& id)
	{
		std::string type;
		in >> type;
		if (type == "vertex") in >> vertex_number[id];
	}
	
	bool parseHeader(std::istream& in, const unsigned int& id)
	{
		uintptr_t act = 0;
		
		std::string line;
		while (std::getline(in, line))
		{
			std::istringstream ls(line);
			std::string token;
			ls >> token;
			if (token == "ply" || token == "PLY" || token == "") { act = token == "ply" ? (act << 1) | 0x1 : act; continue; }
			else if (token == "format")   continue; //  read_header_format(ls);
			else if (token == "element") { act = (act << 1) | 0x1; read_header_element(ls, id); continue; }//  read_header_element(ls);
			else if (token == "property") continue; //  read_header_property(ls);
			else if (token == "end_header") { act = (act << 1) | 0x1; break; }
			else return false;
		}
 
		return !(act ^ headerMask); //xor
	}

	void parseContext(std::istream* in, std::istream* tex, std::set<tuple>* points, const unsigned int& id)
	{
		float x,y,z;
		uint16_t r,g,b;

		float minX = FLT_MAX, maxX = FLT_MIN;
		float minY = FLT_MAX, maxY = FLT_MIN;
		float minZ = FLT_MAX, maxZ = FLT_MIN;

		for (int i = 0; i < this->vertex_number[id]; ++i)
		{
			//TODO: any error??? lets think
			//mx.lock();

			(*in) >> x >> y >> z >> r >> g >> b;
			
			auto res = points->insert(tuple(P(x, -y, z), std::vector<T>()));
			std::vector<T>* non_const = const_cast<std::vector<T>*>(&res.first->second);

			if (res.second != true)
			{
				qDebug() << "Already in";
			}

			if (UV)
			{
				if (tex != NULL)
				{
					float u, v;
					(*tex) >> u >> v;
					non_const->push_back(T(id, std::pair<float, float>(u,v)));
				}
				else
					non_const->push_back(T(id, std::pair<float, float>(0.0f, 0.0f)));
			}
			//else
			//	non_const->push_back(T(id, std::pair<float, std::pair<std::pair<uint8_t, uint8_t>, uint8_t> >(r,g,b)));
			

			minX = minX > x ? x : minX;
			maxX = maxX < x ? x : maxX;
			minY = minY > -y ? -y : minY;
			maxY = maxY < -y ? -y : maxY; 
			minZ = minZ > z ? z : minZ;
			maxZ = maxZ < z ? z : maxZ;
			//mx.unlock();
		}
	
		if (!box.isInside(minX, minY, minZ) || !box.isInside(maxX, maxY, maxZ))
		{
			box = Box(minX < box.getXMin() ? minX : box.getXMin()
				    , maxX > box.getXMax() ? maxX : box.getXMax()
					, minY < box.getYMin() ? minY : box.getYMin()
					, maxY > box.getYMax() ? maxY : box.getYMax()
					, minZ < box.getZMin() ? minZ : box.getZMin()
					, maxZ > box.getZMax() ? maxZ : box.getZMax() );
		}
		qDebug() << "\tOUR BOX IS: " << box.getXMin() << " " << box.getXMax() << " " << box.getYMin() << " " << box.getYMax() << " " << box.getZMin() << " " << box.getZMax();

		//in->clear();
		//return true;
	}
	
	point_t_map getCloudRandom();
	point_t_map  getCloudSpecific();
	point_t_map getCloudFromFile();

private:
	cloud_sources cloud_s;
	file_extension ext;
	bool started;
	int run_count;

	Box box;
	//*******************
	// Random
	unsigned int n_point;

	//*******************
	// File
	//boost::thread_group thread_pool;
	//std::vector<boost::thread*> thread_ptrs;
	//boost::mutex mx;

	char counter_prefix[20];	//changable counter
	char* file_name;			//full filename
	char file_path[100];			//filepath
	char* adress_c;				//ip adress
	char* postfix;				//.ply
		
	unsigned int file_index;
	unsigned int file_index_start;
	
	static const int kinectSources = 4;
	const char* ips[kinectSources] = { "192.168.0.75","192.168.0.71","192.168.0.100","127.0.0.1"};
	int vertex_number[kinectSources];

	static const uintptr_t headerMask = 0b111; // az alap információk megvannak-e
	
	//calibration stuff...
	//*******************
};

/******************************************





***************************************/
template <typename P, typename T, bool UV = true>
void C_IO<P,T,UV>::loadCalibrations()
{
	const char cali_prefix[13] = "calibration_";
	char file[50];

	for (int i = 0; i<kinectSources; ++i)
	{
		//todo: file is missing
		memset(file, 0, sizeof(file));
		strcat(file, cali_prefix);
		strcat(file, ips[i]);
		strcat(file, ".txt");

		std::ifstream in(file, std::ios::binary);

		if (in.fail()) continue;

		/*float x,y,z;
		in >> x >> y >> z;
		translate[i] = Coord(x,y,z);

		for (int j=0; j < 3;++j)
		{
		rotate[i][j] = ...;
		}

		file.close();*/
	}
}

template <typename P, typename T, bool UV = true>
char* C_IO<P, T, UV>::getFile(const int& i )
{
	switch (ext)
	{
		case PLY:
		{
			return assemble_ply_filename(i);
		}
		case OBJ:
		{
			return assemble_obj_filename();
		}
		default:
		{
			return "";
		}
	}

	return "";
}

template <typename P, typename T, bool UV = true>
char* C_IO<P, T, UV>::getTexFile(const int& i)
{
	char c_file_index[50];
	itoa(file_index, c_file_index, 10);
	char* file_index_c = c_file_index;

	strncpy(counter_prefix + strlen(counter_prefix) - strlen(c_file_index), c_file_index, strlen(file_index_c));
	strncpy(file_path + (file_name - file_path), counter_prefix, strlen(counter_prefix));
	strncpy(file_path + (adress_c - file_path), ips[i], strlen(ips[i]));

	char tmp[100];
	memset(tmp, 0, sizeof(tmp));
	strncpy(tmp, file_path, adress_c - file_path + strlen(ips[i]));
	strcat(tmp, ".tex");
	strcpy(file_path, tmp);

	return file_path;
}

template <typename P, typename T, bool UV = true>
char* C_IO<P, T, UV>::getPngFile(const int& i)
{
	char c_file_index[50];
	itoa(file_index, c_file_index, 10);
	char* file_index_c = c_file_index;

	strncpy(counter_prefix + strlen(counter_prefix) - strlen(c_file_index), c_file_index, strlen(file_index_c));
	strncpy(file_path + (file_name - file_path), counter_prefix, strlen(counter_prefix));
	strncpy(file_path + (adress_c - file_path), ips[i], strlen(ips[i]));

	char tmp[100];
	memset(tmp, 0, sizeof(tmp));
	strncpy(tmp, file_path, adress_c - file_path + strlen(ips[i]));
	strcat(tmp, ".png");
	strcpy(file_path, tmp);

	return file_path;
}

template <typename P, typename T, bool UV = true>
typename C_IO<P, T, UV>::point_t_map C_IO<P, T, UV>::getCloudFromFile()
{
	std::set<tuple> res;

	switch (ext)
	{
		case PLY:
		{
			getCloudFromPly(res);
			return point_t_map(res.begin(), res.end());
		}
		case OBJ:
		{
			getCloudFromObj(res);
			return point_t_map(res.begin(), res.end());

		}
		default:
		{
			return point_t_map(res.begin(),res.end());
		}
	}

	return point_t_map(res.begin(), res.end());
}

template <typename P, typename T, bool UV = true>
typename C_IO<P, T, UV>::point_t_map C_IO<P, T, UV>::getCloudRandom()
{
	srand(static_cast <unsigned> (time(0)));

	std::set<tuple> res;

	float xmin = box.getXMin();
	float xmax = box.getXMax();
	float ymin = box.getYMin();
	float ymax = box.getYMax();
	float zmin = box.getZMin();
	float zmax = box.getZMax();

	while (res.size() != n_point)
	{
		float x = xmin + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (xmax - xmin)));
		float y = ymin + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (ymax - ymin)));
		float z = zmin + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (zmax - zmin)));

		std::pair<std::set<tuple>::iterator, bool> pos = res.insert(tuple(P(x, y, z), std::vector<T>()));
		//T tmp = T(0, std::pair<float, float>(0.0f, 0.0f));
		//pos.first->second.push_back(tmp);
	}

	return point_t_map(res.begin(),res.end());
}

template <typename P, typename T, bool UV = true>
typename C_IO<P, T, UV>::point_t_map C_IO<P, T, UV>::getCloudSpecific()
{
	point_t_map tmp;

	/*tmp.push_back(Data_3D(2, 3, 3));
	tmp.push_back(Data_3D(5, 4, 2));
	tmp.push_back(Data_3D(9, 6, 7));
	tmp.push_back(Data_3D(4, 7, 9));
	tmp.push_back(Data_3D(8, 1, 5));
	tmp.push_back(Data_3D(7, 2, 6));
	tmp.push_back(Data_3D(9, 4, 1));
	tmp.push_back(Data_3D(8, 4, 2));
	tmp.push_back(Data_3D(9, 7, 8));
	tmp.push_back(Data_3D(6, 3, 1));
	tmp.push_back(Data_3D(3, 4, 5));
	tmp.push_back(Data_3D(1, 6, 8));
	tmp.push_back(Data_3D(9, 5, 3));
	tmp.push_back(Data_3D(2, 1, 3));
	tmp.push_back(Data_3D(8, 7, 6));
	*/
	return tmp;
}

template <typename P, typename T, bool UV = true>
std::vector<std::string> C_IO<P, T, UV>::get_file_name_list()
{
	std::vector<std::string> res;

	for (int i = 0; i < 4; ++i)
	{
		res.push_back(std::string(getPngFile(i)));
	}
	for (auto& it : res)
	{
		qDebug() << it.c_str();
	}
	return res;
}
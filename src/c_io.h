#pragma once
//#include "cgal_types.h"
#include "types.h"
#include "defines.h"
#include <vector>
#include <fstream>
#include <string.h>

#include <cstdlib>
#include <stdio.h>

#include <set>
//#include <boost\thread\mutex.hpp>
//#include <boost\thread\thread.hpp>

#include <boost\filesystem\fstream.hpp>

#include <QDebug>

template <typename Data_3D>
class C_IO
{
	enum cloud_sources
	{
		RANDOM,
		FILE
	};

public:
	C_IO() : box(Box(-0.1f, 0.1f)), n_point(0), file_index(1), cloud_s(FILE)
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
	void setFileName(char* fn)
	{
		strcpy(file_path, fn);
		//file_path = fn;
		file_name=strrchr(file_path,'/');
		file_name++;
	  
		adress_c = strrchr(file_name,'_');
		adress_c++;
		
		postfix = strrchr(adress_c,'.');	//TODO: if there is no ip
		
		memset(counter_prefix, 0, sizeof(counter_prefix));
		strncpy(counter_prefix, file_name, adress_c - file_name -1);
	  
		file_index = file_index_start = atoi(counter_prefix);
		
	    //memset(counter_prefix_c,'_',counter_prefix_c);
		//strncpy (file_path + (file_name - file_path)  ,counter_prefix_c, strlen(counter_prefix_c));
		//strncpy (file_path + (adress_c - file_path)  ,ipadress, strlen(ipadress));
	}

	void resetFileIndex() { file_index = file_index_start; }

	bool checkFileType() 
	{
		bool existOne = false;
		
		for (int i=0; i < kinectSources; ++i)
		{
			std::ifstream file(getFile(i), std::ios::binary);
			if (!file.fail()) existOne = true;
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

	ProcessReturns load_process(std::vector<Data_3D>& points)
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
				if (!nextFileExist())
				{
					//qDebug() << "Next FIle doesnt exist!";
					resetFileIndex();
					/*if (!repeatable)*/ return NO_INPUT_DATA;
				}

				points.clear();
				points = getCloudFromFile();

				return CLOUD_LOAD_SUCCESS;
			}
			default:
				return UNDEFINED_ERROR;
		}
	}

protected:
	void loadCalibrations();
    char* getFile(const int& i);

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

	void parseContext(std::istream* in, std::set<Data_3D>* points, const unsigned int& id)
	{
		float x,y,z;
		int r,g,b;

		float minX = FLT_MAX, maxX = FLT_MIN;
		float minY = FLT_MAX, maxY = FLT_MIN;
		float minZ = FLT_MAX, maxZ = FLT_MIN;

		for (int i = 0; i < this->vertex_number[id]; ++i)
		{
			//TODO: any error??? lets think
			//mx.lock();

			(*in) >> x >> y >> z >> r >> g >> b;
			points->insert(Data_3D(x,-y,z));

			minX = minX > x ? x : minX;
			maxX = maxX < x ? x : maxX;
			minY = minY > -y ? -y : minY;
			maxY = maxY < -y ? -y : maxY; 
			minZ = minZ > z ? z : minZ;
			maxZ = maxZ < z ? z : maxZ;
			//mx.unlock();

			//TODO: eredeti koordrendszer - talán a színekhez??
			// R * (v-t) - inverze:  transp(R)*v + t
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
	
	std::vector<Data_3D> getCloudRandom();
	std::vector<Data_3D>  getCloudSpecific();
	std::vector<Data_3D> getCloudFromFile();

private:
	cloud_sources cloud_s;
	bool started;
	int run_count;

	//*******************
	// Random
	Box box;
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
	
	//calibrations
	//TODO CGAL osra cseréljük le
//	Coord translate[kinectSources];
	//todo: 3x3 matrix, saját típus kéne (mátrix, mátrixszorzással) :/ vagy a coordot lecserélni
//	std::vector<Coord> rotation[kinectSources];
	
	//*******************
};

template <typename Data_3D>
void C_IO<Data_3D>::loadCalibrations()
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

template <typename Data_3D>
char* C_IO<Data_3D>::getFile(const int& i)
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

template <typename Data_3D>
std::vector<Data_3D> C_IO<Data_3D>::getCloudFromFile()
{
	std::set<Data_3D> tmp;

	//boost::filesystem::ifstream files[kinectSources];
	time_p t1 = hrclock::now();

	for (int i = 0; i<kinectSources; ++i)
	{
		boost::filesystem::ifstream file(getFile(i));

		if (file.fail()) continue;

		if (!parseHeader(file, i)) { file.close(); continue; }

		//thread_ptrs.push_back(thread_pool.create_thread( boost::bind(&C_IO::parseContext, this, &files[i], &tmp, i)));
		//if (!parseContext(file,tmp)) { file.close(); continue; }
		parseContext(&file, &tmp, i);
		file.close();
	}

	/*thread_pool.join_all();
	for (auto it : thread_ptrs)
	thread_pool.remove_thread(it);
	thread_ptrs.clear();
	*/

	qDebug() << "LOAD under " << (boost::chrono::duration_cast<milliseconds>(hrclock::now() - t1)).count() << " ms with " << tmp.size() << " points";

	file_index++;
	qDebug() << "\tOUR BOX IS: " << box.getXMin() << " " << box.getXMax() << " " << box.getYMin() << " " << box.getYMax() << " " << box.getZMin() << " " << box.getZMax();

	//box = Box(-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
	//box = Box(minx, maxx, miny, maxy, minz, maxz);

	return std::vector<Data_3D>(tmp.begin(), tmp.end());
}

template <typename Data_3D>
std::vector<Data_3D> C_IO<Data_3D>::getCloudRandom()
{
	srand(static_cast <unsigned> (time(0)));

	std::set<Data_3D> res;

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

		res.insert(Data_3D(x, y, z));
	}

	return std::vector<Data_3D>(res.begin(), res.end());
}

template <typename Data_3D>
std::vector<Data_3D> C_IO<Data_3D>::getCloudSpecific()
{
	std::vector<Data_3D> tmp;

	tmp.push_back(Data_3D(2, 3, 3));
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

	return tmp;
}

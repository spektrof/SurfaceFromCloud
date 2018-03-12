#pragma once
#include "defines.h"

#include <boost/lockfree/queue.hpp>
#include <boost/thread/thread.hpp>

#include <boost/atomic.hpp>
#include <vector>

#include <boost/chrono.hpp>

//TODO: make it more flexible - later
template<class K, typename T>
class TaskSchedular
{
	typedef void (K::*builder)(T*,bool,const int&);
	typedef void (K::*nearest)(T*,bool,const int&);

public:
	TaskSchedular()
	{
		wait_time = boost::chrono::milliseconds(0);
	}
	~TaskSchedular() { }

	bool isEmptyTask() const { return tasks.empty(); }
	bool isEmptyThreadPool() const { return thread_pool.size() == 0; }

	void addTask(T* task)
	{
		while (!tasks.push(task))
		{
			std::cout << "there is no space in task or double push\n";
		}	//do it until the task get into the que
	}

	void addSubscribeWithoutAutomechanism(K* src, builder build)		//try to launch a thread
	{
		if (!isEmptyTask() && thread_pool.size() < maximumThread)
		{
			T* act;
			while (!tasks.pop(act)) {}

			thread_ptrs.push_back(thread_pool.create_thread(boost::bind(build, src, act, false, 0)));	//cref actal unable to read volt
		}
	}

	void addSubscribeShit(K* src, nearest _n)
	{
		if (!isEmptyTask() && thread_pool.size() < maximumThread)
		{
			T* act;
			while (!tasks.pop(act)) {}

			thread_ptrs.push_back(thread_pool.create_thread(boost::bind(_n, src, act, false, 0)));	//cref actal unable to read volt
		}
	}

	void joinAll()
	{
		boost::chrono::high_resolution_clock::time_point t1 = boost::chrono::high_resolution_clock::now();

		thread_pool.join_all();
		for (auto it : thread_ptrs)
			thread_pool.remove_thread(it);
		thread_ptrs.clear();

		boost::chrono::high_resolution_clock::time_point t2 = boost::chrono::high_resolution_clock::now();
		wait_time += (boost::chrono::duration_cast<boost::chrono::milliseconds>(t2 - t1));
	}

	float getJoinTime()
	{
		float tmp = wait_time.count();
		wait_time = boost::chrono::milliseconds(0);
		return tmp;
	}

	void setCloudThreadCapacity(const unsigned int& tc) { n_cloudThread = tc; }
	unsigned int getCloudThreadCapacity() const { return n_cloudThread; }

private:
	boost::lockfree::queue<T*, boost::lockfree::capacity<MAXIMUMTASK> > tasks;
	boost::thread_group thread_pool;
	std::vector<boost::thread*> thread_ptrs;

	const size_t maximumThread = (size_t) MAXIMUMTHREAD; // TODO: REFACTOR
	unsigned int n_cloudThread;

	boost::chrono::milliseconds wait_time;
};
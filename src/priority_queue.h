#pragma once
#include <queue>
#include <QDebug>

#define max_prior(a,b) (((a)>(b))?(a):(b))

namespace power_crust_utils
{
	template<
		class T,
		class Container = std::vector<T>,
		class _T = int,
		class Compare = std::less<typename Container::value_type> >
	class priority_queue_for_pole : public std::priority_queue<T, Container, Compare>
	{
	public:
		typedef typename std::priority_queue<T, Container, Compare>::container_type::const_iterator const_iterator;
		typedef typename std::priority_queue<T, Container, Compare>::container_type::iterator iterator;

		iterator find(const _T&val)
		{
			auto first = this->c.begin();
			auto last = this->c.end();
			while (first != last) {
				if (first->second == val) return first;
				++first;
			}
			return last;
		}

		bool find_higher_prior(const float& val)
		{
			auto first = this->c.begin();
			auto last = this->c.end();
			while (first != last) {
				if (first->first > val) return true;
				++first;
			}
			return false;
		}

		void update(iterator& update_pos)
		{
			T upd_element = *update_pos;

			this->c.erase(update_pos, update_pos + 1);
			std::make_heap(this->c.begin(), this->c.end(), this->comp);

			float in_val = upd_element.second->getInValue();
			float out_val = upd_element.second->getOutValue();
			if (in_val > 0 && out_val > 0)
				upd_element.first = abs(in_val - out_val) - 1.0f;
			else
				upd_element.first = max_prior(in_val,out_val);

		//	qDebug() << "\t" << upd_element.second << " " << in_val << ", " << out_val << " result: " << upd_element.first;

			push(upd_element);
		//	return upd_element;
		}

		int count(const _T&val)
		{
			int c = 0;
			auto first = this->c.begin();
			auto last = this->c.end();
			while (first != last) {
				if (first->second == val) c++;
				++first;
			}
			return c;
		}

		iterator operator [](const unsigned int& ind)
		{
			//assert(this->c.size()-1 >= ind);
			return (this->c.begin() + ind);
		}

		void erase(iterator e_s, iterator e_e)
		{
			this->c.erase(e_s, e_e);
		}

		T erase_one(const int& ind)
		{
			//assert(this->c.size() - 1 >= ind);
			T tmp = *(this->c.begin() + ind);
			iterator it = this->c.begin() + ind;
			this->c.erase(it, it + 1);
			return tmp;
		}

		iterator begin()
		{
			return this->c.begin();
		}

		iterator end()
		{
			return this->c.end();
		}

	};
}
#pragma once
#include <iostream>
#include <cstdio>
#include <queue>

using namespace std;

template<typename T,
		typename Sequence = std::vector<T>,
		typename Compare = std::less<typename Sequence::value_type> >
class _priority_queue : public priority_queue<T, Sequence, Compare> {
	using priority_queue::priority_queue;
public:
	bool remove(const T& value) {
		auto it = std::find(this->c.begin(), this->c.end(), value);
		if (it != this->c.end()) {
			this->c.erase(it);
			std::make_heap(this->c.begin(), this->c.end(), this->comp);
			return true;
		}
		else {
			return false;
		}
	}

	// Call this after update cost
	void re_push(const T& target) {
		remove(target);
		push(target);
	}
};

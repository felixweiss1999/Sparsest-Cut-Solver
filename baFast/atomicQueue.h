#pragma once
#include "TreeDecomposition.h"

//template<typename T>
class atomicQueue
{
public:
	void push(const size_t& value) {
		q.push(value);
	}
	size_t pop() {
		std::lock_guard<std::mutex> lock(m);
		if (q.empty())
			return 0;
		size_t obj = q.front();
		size_t obj2 = obj;
		q.pop();
		return obj2;
	}

private:
	std::queue<size_t> q;
	mutable std::mutex m;
};


#include "BFQuadPool.hh"

void* BFQuadPool::Chunk::alloc()
{
	if (_cur + elem_size > _end)
		return 0;
	void *res = _cur;
	_cur += elem_size;
	return res;
}

BFQuadPool::~BFQuadPool()
{
	while (!_chunks.empty())
	{
		Chunk *chunk = _chunks.front();
		_chunks.pop_front();
		delete chunk;
	}
}

void* BFQuadPool::alloc()
{
	void *res = _chunks.front()->alloc();
	if (!res)
	{
		_chunks.push_front(new Chunk());
		res = _chunks.front()->alloc();
	}
	return res;
}

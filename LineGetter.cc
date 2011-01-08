#include <cstring>
#include "LineGetter.hh"

const char* LineGetter::next(bool except)
{
	while (true)
	{
		_is.getline(_buf, _buf_size);
		if (_is.eof())
		{
			if (except)
				throw UnexpectedEOF();
			return 0;
		}
		if (_is.fail())
			throw ReadFailure();
		_line_nr++;
		if (_buf[0] && _buf[0] != _cmt)
			break;
	}
	return _buf;
}

void LineGetter::setBufferSize(size_t size)
{
	char *new_buf = new char[size];
	size_t n = std::min(size, _buf_size) - 1;
	std::strncpy(new_buf, _buf, n);
	new_buf[n] = '\0';
	delete[] _buf;
	_buf = new_buf;
	_buf_size = size;
}

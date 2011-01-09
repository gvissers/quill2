#include "IndentingStreambuf.hh"

IndentingStreambuf::int_type IndentingStreambuf::overflow(int_type c)
{
	if (traits_type::eq_int_type(c, traits_type::eof()))
		return c;
	if (traits_type::eq(traits_type::to_char_type(c), '\n'))
	{
		_indent_now = true;
	}
	else if (_indent_now)
	{
		for (int i = 0; i < _level; i++)
			_buf->sputn(_indent.data(), _indent.size());
		_indent_now = false;
	}
	_buf->sputc(c);
	return c;
}


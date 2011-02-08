#include "io/IndentFilter.hh"

IndentFilter::int_type IndentFilter::operator()(std::streambuf *sb, int_type c)
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
		{
			std::streamsize n = _indent.size();
			if (sb->sputn(_indent.data(), n) < n)
				return traits_type::eof();
		}
		_indent_now = false;
	}
	return sb->sputc(c);
}


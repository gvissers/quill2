#include "io/JobFilter.hh"

JobFilter::int_type JobFilter::operator()(std::streambuf *sb)
{
	int_type ic;
	char_type c;

	ic = sb->sbumpc();
	if (traits_type::eq_int_type(ic, traits_type::eof()))
		return ic;
	c = traits_type::to_char_type(ic);

	if (traits_type::eq(c, ','))
		return _string_open ? ic : traits_type::to_int_type(' ');

	if (traits_type::eq(c, '#') && !_string_open)
	{
		while (!traits_type::eq(c, '\n'))
		{
			ic = sb->sbumpc();
			if (traits_type::eq_int_type(ic, traits_type::eof()))
				return ic;
			c = traits_type::to_char_type(ic);
		}
	}
	else if (traits_type::eq(c, '\\') && _string_open)
	{
		_escaped = !_escaped;
	}
	else if (traits_type::eq(c, '\n'))
	{
		_string_open = '\0';
		_escaped = false;
	}
	else if (traits_type::eq(c, '\'') || traits_type::eq(c, '"'))
	{
		if (!_string_open)
		{
			_string_open = c;
		}
		else if (c == _string_open && !_escaped)
		{
			_string_open = '\0';
		}
	}
	return ic;
}

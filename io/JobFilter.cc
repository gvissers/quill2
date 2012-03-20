#include "io/JobFilter.hh"

JobFilter::int_type JobFilter::operator()(std::streambuf *sb)
{
	int_type ic;
	char_type c;

	ic = sb->sbumpc();
	if (traits_type::eq_int_type(ic, traits_type::eof()))
		return ic;
	c = traits_type::to_char_type(ic);

	if (_string_open)
	{
		if (traits_type::eq(c, '\n'))
		{
			_string_open = '\0';
			_escaped = false;
		}
		else if (_escaped)
		{
			_escaped = false;
		}
		else if (traits_type::eq(c, '\\'))
		{
			_escaped = true;
		}
		else if (traits_type::eq(c, _string_open))
		{
			_string_open = '\0';
		}
	}
	else
	{
		if (traits_type::eq(c, ','))
			return traits_type::to_int_type(' ');

		if (traits_type::eq(c, '#'))
		{
			while (!traits_type::eq(c, '\n'))
			{
				ic = sb->sbumpc();
				if (traits_type::eq_int_type(ic, traits_type::eof()))
					return ic;
				c = traits_type::to_char_type(ic);
			}
		}
		else if (traits_type::eq(c, '\'') || traits_type::eq(c, '"'))
		{
			_string_open = c;
		}
	}
	return ic;
}

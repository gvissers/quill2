#include "CommentFilter.hh"

CommentFilter::int_type CommentFilter::operator()(std::streambuf *sb)
{
	int_type ic;
	char_type c;

	ic = sb->sbumpc();
	if (traits_type::eq_int_type(ic, traits_type::eof()))
		return ic;
	c = traits_type::to_char_type(ic);

	if (_cmt_chars.find(c) != std::string::npos)
	{
		while (!traits_type::eq(c, '\n'))
		{
			ic = sb->sbumpc();
			if (traits_type::eq_int_type(ic, traits_type::eof()))
				return ic;
			c = traits_type::to_char_type(ic);
		}
	}
	else if (traits_type::eq(c, ','))
	{
                ic = traits_type::to_int_type(' ');
	}

	return ic;
}

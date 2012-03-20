#include "Job.hh"
#include "support.hh"

class ValueGetter
{
public:
	ValueGetter(): _value() {}
	
	JobIStream& get(JobIStream& is);

	const std::string& value() const { return _value; }

private:
	std::string _value;

	JobIStream& getQuoted(char quote, JobIStream& is, bool keepEscape=false);
	JobIStream& getBlock(JobIStream& is);

	void warn(const std::string& msg)
	{
		std::cerr << msg << "\n";
	}
};

JobIStream& ValueGetter::get(JobIStream& is)
{
	char c;

	is >> c;
	if (is.eof())
		throw UnexpectedEOF();
	
	if (c == '"' || c == '\'')
		return getQuoted(c, is);
	if (c == '{')
		return getBlock(is);

	is.unget();
	is >> _value;
	return is;
}

JobIStream& ValueGetter::getQuoted(char quote, JobIStream& is, bool keepEscape)
{
	int c;
	bool escaped = false;
	
	_value.clear();
	while (true)
	{
		c = is.get();
		if (is.eof() || (c == quote && !escaped) || c == '\n')
			break;

		if (keepEscape)
		{
			_value += c;
		}
		else if (escaped)
		{
			switch (c)
			{
				case '\\': _value += '\\'; break;
				case 'a' : _value += '\a'; break;
				case 'b' : _value += '\b'; break;
				case 'f' : _value += '\f'; break;
				case 'n' : _value += '\n'; break;
				case 'r' : _value += '\r'; break;
				case 't' : _value += '\t'; break;
				case 'v' : _value += '\v'; break;
				default  : _value += c;
			}
			escaped = false;
		}
		else if (c == '\\')
		{
			escaped = true;
		}
		else
		{
			_value += c;
		}
	}

	if (c != quote)
		warn("Unterminated string value");

	return is;
}

JobIStream& ValueGetter::getBlock(JobIStream& is)
{
	bool open = true;

	_value.clear();
	while (open)
	{
		int c = is.get();
		if (is.eof())
		{
			is.getLine();
			c = is.get();
			if (is.eof())
				throw UnexpectedEOF();
		}

		switch (c)
		{
			case '\'':
			case '"':
			{
				ValueGetter vg;
				vg.getQuoted(c, is, true);
				_value += char(c) + vg.value() + char(c);
				break;
			}
			case '{':
			{
				ValueGetter vg;
				vg.getBlock(is);
				_value += '{' + vg.value() + '}';
				break;
			}
			case '}':
				open = false;
				break;
			default:
				_value += c;
		}
	}

	return is;
}

JobIStream& operator>>(JobIStream& is, ValueGetter& vg)
{
	return vg.get(is);
}

JobIStream& Job::scan(JobIStream& is)
{
	while (true)
	{
		is >> getline;
		if (is.eof())
			break;

		std::string key;
		is >> key;

		ValueGetter vg;
		is >> vg;

		_params[lower(key)] = vg.value();
	}

	return is;
}

const std::string Job::get(const std::string& key) const
{
	ParamMap::const_iterator it = _params.find(lower(key));
	if (it == _params.end())
		throw InvalidIndex(key);
	return it->second;
}

const std::string Job::get(const std::string& key,
	const std::string& defval) const
{
	try
	{
		return get(key);
	}
	catch (const InvalidIndex&)
	{
		return defval;
	}
}

template <>
const Geometry Job::get<Geometry>(const std::string& key) const
{
	std::istringstream is(get(key));
	JobIStream jis(is);
	Geometry geom;
	jis >> geom;
	return geom;
}
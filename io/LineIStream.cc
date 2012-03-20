#include "io/LineIStream.hh"

void LineIStream::getLine()
{
	if (!_lines.empty())
	{
		clear();
		stream().str(_lines.front());
		_lines.pop_front();
	}
	else
	{
		std::string line = nextLine();
		if (line.empty())
		{
			setstate(eofbit);
		}
		else
		{
			clear();
			stream().str(line);
		}
	}

}

std::string LineIStream::nextLine()
{
	if (_input.eof())
		return std::string();

	// Loop while line is empty
	while (true)
	{
		std::string line;
		std::getline(_input, line);
		if (_input.eof() && line.empty())
			return std::string();
		if (_input.fail())
			throw ReadFailure();
		_line_nr++;
		if (!line.empty())
			return line + '\n';
	}
}


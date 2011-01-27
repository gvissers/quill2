#include "JobIStream.hh"

void JobIStream::getLine()
{
	if (_lines.empty())
	{
		std::string line = nextLine();
		if (line.empty())
		{
			setstate(eofbit);
			return;
		}
		_lines.push_back(line);
	}

	clear();
	str(_lines.front());
	_lines.pop_front();
}

std::string JobIStream::nextLine()
{
	if (_input.eof())
		return std::string();

	// Loop while line is empty
	while (true)
	{
		std::string line;
		std::getline(_input, line);
		if (_input.fail())
			throw ReadFailure();
		if (!line.empty() || _input.eof())
			return line;
	}
}


#ifndef JOB_HH
#define JOB_HH

/*!
 * \file Job.hh
 * \brief Definition of the Job class
 */

#include <map>
#include <util.hh>
#include "io/JobIStream.hh"
#include "Geometry.hh"

/*!
 * \brief Class for Quill jobs
 *
 * Class Job holds the parameters that tell Quill what to compute, such as the
 * type of calculation, the basis set, and the geometry of the system. The
 * parameters are stored as key,value pairs of strings, where the key is
 * case-insensitive. Extraction of values to other types is supported, the
 * conversion from string is done on the fly.
 */
class Job
{
public:
	//! Local typedef for the parameter map
	typedef std::map<std::string, std::string> ParamMap;

	/*!
	 * \brief Constructor
	 *
	 * Create a new and empty job
	 */
	Job(): _params() {}
	/*!
	 * \brief Constructor
	 *
	 * Create a new job and read the parameters from input stream \a is.
	 * \param is The job input stream to read this job from
	 */
	Job(JobIStream& is): _params() { scan(is); }

	//! Clear the job, removing all parameters
	void clear() { _params.clear(); }

	/*!
	 * \brief Read a job
	 *
	 * Read a job from input stream \a is, and store the parameters found.
	 * Parameters currently stored in the job are preserved, but will be
	 * overwritten if a parameter with the same name is found in the input.
	 * \param is The job input stream to read this job from
	 * \return The updated input stream
	 */
	JobIStream& scan(JobIStream& is);

	/*!
	 * \brief Extract a string value
	 *
	 * Extract the value for parameter \a key from the job. If it is not
	 * found, an InvalidIndex exception is thrown.
	 */
	const std::string get(const std::string& key) const;
	/*!
	 * \brief Extract a string value
	 *
	 * Extract the value for parameter \a key from the job. If it is not
	 * found, return \a defval.
	 */
	const std::string get(const std::string& key, const std::string& defval) const;
	/*!
	 * \brief Extract avalue
	 *
	 * Extract the value for parameter \a key from the job and convert it to
	 * type \a T. If the key is not found, an InvalidIndex exception is
	 * thrown. If the conversion cannot be made. a Li::InvalidConversion
	 * exception is thrown.
	 */
	template <typename T>
	const T get(const std::string& key) const
	{
		return Li::fromString<T>(get(key));
	}

private:
	//! The parameters themselves
	ParamMap _params;
};

template <>
const Geometry Job::get<Geometry>(const std::string& key) const;

namespace {

/*!
 * \brief Read a job
 *
 * Read a job from input stream \a is, and store the parameters found in \a job.
 * Parameters currently stored in the job are preserved, but will be overwritten
 * if a parameter with the same name is found in the input.
 * \param is  The input stream to read from
 * \param job The job in which to store the parameters
 * \return The updated input stream
 */
JobIStream& operator>>(JobIStream& is, Job& job)
{
	return job.scan(is);
}

} // namespace

#endif // JOB_HH
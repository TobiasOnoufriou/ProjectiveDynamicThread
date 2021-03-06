#include "simpleTimer.h"
#include <string>
#include <iostream>

//__int64 simpleTimer::start_time = GetTimeMs64();

simpleTimer::simpleTimer()
{
	start_time = GetTimeMs64();
}

simpleTimer::~simpleTimer()
{
}

void simpleTimer::start()
{
	start_time = GetTimeMs64();
}

// Function adjusted to return __int64 value edited by Tobias Onoufriou.
__int64 simpleTimer::stop( std::string some_desriptor )
{
	__int64 end_time = GetTimeMs64();
	return end_time - start_time;
}

__int64 simpleTimer::GetTimeMs64()
{
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	* to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	__int64 ret = li.QuadPart;
	ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
	ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

	return ret;
}
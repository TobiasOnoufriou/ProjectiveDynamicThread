#pragma once

#include <Windows.h>
#include <iostream>

class simpleTimer
{
public:
	simpleTimer( void );
	~simpleTimer( void );

	// methods
	void start( void );
	// function below adjusted to return __int64 rather than just printing to the console edited by Tobias Onoufriou.
	__int64 stop( std::string some_desriptor );
	
private:
	// method
	__int64 GetTimeMs64( void );

	// member
	__int64 start_time;
};
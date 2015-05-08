/***********************************************************************
Timer - Simple class to provide high-resolution timers in a somewhat OS-
independent fashion. This is the windows-compatible version using
performance counters.
Copyright (c) 2001-2005 Oliver Kreylos

This file is part of the Miscellaneous Support Library (Misc).

The Miscellaneous Support Library is free software; you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

The Miscellaneous Support Library is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with the Miscellaneous Support Library; if not, write to the Free
Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA
***********************************************************************/

#include <Misc/Timer.h>

namespace Misc {

/**********************
Methods of class Timer:
**********************/

#ifdef WIN32

Timer::Timer(void)
	{
	/* Test if a performance counter is available, and query its resolution: */
	if(QueryPerformanceFrequency(&ticksPerSecond))
		{
		/* Query the current counter value: */
		QueryPerformanceCounter(&lastMeasured);
		}
	else
		{
		/* Use the regular inaccurate timeGetTime() function instead: */
		ticksPerSecond.HighPart=-1;
		lastMeasured.LowPart=timeGetTime();
		}
	}

void Timer::elapse(void)
	{
	/* Check which timer version is used: */
	LARGE_INTEGER newMeasured;
	if(ticksPerSecond.HighPart>=0)
		{
		/* Query the current counter value: */
		QueryPerformanceCounter(&newMeasured);
		
		/* Calculate the elapsed time: */
		elapsedSeconds=double(newMeasured.QuadPart-lastMeasured.QuadPart)/double(ticksPerSecond.QuadPart);
		}
	else
		{
		/* Query the current time: */
		newMeasured.LowPart=timeGetTime();
		
		/* Calculate the elapsed time: */
		elapsedSeconds=double(newMeasured.LowPart-lastMeasured.LowPart)/1000.0;
		}
	
	/* Restart the timer: */
	lastMeasured=newMeasured;
	}

double Timer::peekTime(void) const
	{
	/* Check which timer version is used: */
	if(ticksPerSecond.HighPart>=0)
		{
		/* Query the current counter value: */
		LARGE_INTEGER newMeasured;
		QueryPerformanceCounter(&newMeasured);
		
		/* Calculate the elapsed time: */
		return double(newMeasured.QuadPart-lastMeasured.QuadPart)/double(ticksPerSecond.QuadPart);
		}
	else
		{
		/* Query the current time: */
		DWORD newMeasured=timeGetTime();
		
		/* Calculate the elapsed time: */
		return double(newMeasured-lastMeasured.LowPart)/1000.0;
		}
	}

#else

Timer::Timer(void)
	:elapsedSeconds(0),elapsedMicrons(0)
	{
	gettimeofday(&lastMeasured,0);
	}

void Timer::elapse(void)
	{
	struct timeval newMeasured;
	gettimeofday(&newMeasured,0);
	
	elapsedMicrons=newMeasured.tv_usec-lastMeasured.tv_usec;
	int secondCarry=0;
	if(elapsedMicrons<0) // Correct for end-of-second wraparound
		{
		elapsedMicrons+=1000000;
		secondCarry=1;
		}
	elapsedSeconds=(newMeasured.tv_sec-secondCarry)-lastMeasured.tv_sec;
	if(elapsedSeconds<0) // Correct for end-of-day wraparound
		elapsedSeconds+=86400;
	lastMeasured=newMeasured;
	}

double Timer::peekTime(void) const
	{
	struct timeval newMeasured;
	gettimeofday(&newMeasured,0);
	
	int passedMicrons=newMeasured.tv_usec-lastMeasured.tv_usec;
	int secondCarry=0;
	if(passedMicrons<0) // Correct for end-of-second wraparound
		{
		passedMicrons+=1000000;
		secondCarry=1;
		}
	int passedSeconds=(newMeasured.tv_sec-secondCarry)-lastMeasured.tv_sec;
	if(passedSeconds<0) // Correct for end-of-day wraparound
		passedSeconds+=86400;
	
	return double(passedSeconds)+double(passedMicrons)/1000000.0;
	}

#endif

}

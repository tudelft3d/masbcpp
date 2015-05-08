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

#ifndef MISC_TIMER_INCLUDED
#define MISC_TIMER_INCLUDED

#ifdef WIN32
#include <math.h>
#include <windows.h>
#else
#include <sys/time.h>
#endif

namespace Misc {

class Timer
	{
	/* Elements: */
	private:
	#ifdef WIN32
	LARGE_INTEGER ticksPerSecond; // Number of performance counter ticks per second
	LARGE_INTEGER lastMeasured; // Performance counter value at last measuring point
	double elapsedSeconds; // Number of seconds in the last timing period
	#else
	struct timeval lastMeasured; // Time value at last measuring point
	int elapsedSeconds,elapsedMicrons; // Number of seconds and microseconds in the last timing period
	#endif
	
	/* Constructors and destructors: */
	public:
	Timer(void); // Creates a timer and initializes with the current resource usage (starts timing)
	
	/* Methods: */
	void elapse(void); // Takes a snapshot of the current timer values
	int getSeconds(void) const // Returns the number of whole seconds in the last timing period
		{
		#ifdef WIN32
		return int(floor(elapsedSeconds));
		#else
		return elapsedSeconds;
		#endif
		}
	int getMicrons(void) const // Returns the number of fractional microseconds in the last timing period
		{
		#ifdef WIN32
		return int(floor((elapsedSeconds-floor(elapsedSeconds))*1000000.0));
		#else
		return elapsedMicrons;
		#endif
		}
	double getTime(void) const // Returns the amount of measured time, in seconds
		{
		#if WIN32
		return elapsedSeconds;
		#else
		return double(elapsedSeconds)+double(elapsedMicrons)/1000000.0;
		#endif
		}
	double peekTime(void) const; // Returns the amount of time passed since the last time elapse() was called
	};

}

#endif

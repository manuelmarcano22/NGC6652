/* $Id: piltimer.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
 *
 * This file is part of the VIMOS pipeline library
 * Copyright (C) 2000-2004 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: cizzo $
 * $Date: 2008-10-21 09:10:13 $
 * $Revision: 1.1.1.1 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_UNISTD
#  include <unistd.h>
#endif

#include <sys/time.h>
#include <time.h>
#include <errno.h>
#include <assert.h>

#include "pilmemory.h"
#include "piltimer.h"


#define TIMER_START  0
#define TIMER_END    1

#define TIME_STRING_MAX  28

#ifdef HAVE_GETTIMEOFDAY
#  define GET_TIME(t) \
          gettimeofday(&t, NULL)
#else
#  define GET_TIME(t) \
          t = time(NULL)
#endif


/** 
 * @defgroup pilTimer Timer
 *
 * The module provides basic timer functionalities.
 */

/**@{*/

struct _PilTimer {
#ifdef HAVE_GETTIMEOFDAY
  struct timeval start;
  struct timeval end;
#else
  time_t start;
  time_t end;
#endif
  unsigned int active;
};


/*
 * Abbreviations for the day of the week and the months used for string
 * conversions.
 */

static char *day[] = {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};

static char *month[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};


/*
 * Conver a time value into a PilTime value. The flag selects either the
 * timer's start or end time.
 */

 static PilTime
_pilTimerEncodeTime(PilTimer *timer, int flag, unsigned long *usec)
{

  PilTime tv;

#ifdef HAVE_GETTIMEOFDAY

  if (flag == 0) {
    tv = timer->start.tv_sec + ((PilTime)timer->start.tv_usec / 1.e6);

    if (usec)
      *usec = timer->start.tv_usec;
  }
  else {
    tv = timer->end.tv_sec + ((PilTime)timer->end.tv_usec / 1.e6);

    if (usec)
      *usec = timer->end.tv_usec;
  }

#else

  if (flag == 0) {
    tv = timer->start;
  }
  else {
    tv = timer->end;
  }

  if (usec)
    *usec = 0;

#endif

  return tv;

}


/**
 * @brief
 *  Create a new timer.
 *
 * @return The newly created timer of @c NULL if the timer could not be
 *   created.
 *
 * The function allocates memory for a new timer and initializes the
 * structure to be an inactive timer. The timer is not started. To
 * start the timer use @b pilTimerStart().
 *
 * @see pilTimerStart()
 */

PilTimer *
newPilTimer(void)
{

  PilTimer *timer = pil_malloc(sizeof *timer);

  if (timer) {
    timer->active = 0;
  }

  return timer;

}


/**
 * @brief
 *  Destroy a timer.
 *
 * @param timer  Timer to destory.
 *
 * @return Nothing.
 *
 * The function deallocates the timer @em timer.
 */

void
deletePilTimer(PilTimer *timer)
{

  if (timer)
    pil_free(timer);

  return;

}


/**
 * @brief
 *  Check whether a timer is active.
 *
 * @param timer  Timer to check for activity.
 *
 * @return The function returns 1 if @em timer is active, or 0 otherwise.
 *
 * The function checks whether the timer @em timer is active or not.
 */

int
pilTimerIsActive(PilTimer *timer)
{

  assert(timer != NULL);

  return timer->active ? 1 : 0;

}

/**
 * @brief
 *  Start a timer.
 *
 * @param timer  Timer to start.
 *
 * @return The start time in seconds since 00:00:00 UTC, January 1, 1970.
 *
 * The function starts the timer @em timer, regardless if it is already
 * active or not. The starting time is returned as a @c PilTime value
 * for further usage within the caller. In addition, if @em microseconds is
 * not a @c NULL pointer, it is set to the microseconds making up the fraction
 * of the returned time.
 */

PilTime
pilTimerStart(PilTimer *timer, unsigned long *microseconds)
{

  PilTime tv;

  assert(timer != NULL);

  GET_TIME(timer->start);
  tv = _pilTimerEncodeTime(timer, TIMER_START, microseconds);

  timer->active = 1;

  return tv;

}


/**
 * @brief
 *  Stop a timer.
 *
 * @param timer  Timer to stop.
 *
 * @return The stop time in seconds since 00:00:00 UTC, January 1, 1970.
 *
 * The function stops the timer @em timer. The stopping time is returned as
 * a @c PilTime value for further usage within the caller. In addition, if
 * @em microseconds is not a @c NULL pointer, it is set to the microseconds
 * making up the fraction of the returned time.
 */

PilTime
pilTimerStop(PilTimer *timer, unsigned long *microseconds)
{

  PilTime tv;

  assert(timer != NULL);

  GET_TIME(timer->end);
  tv = _pilTimerEncodeTime(timer, TIMER_END, microseconds);

  timer->active = 0;

  return tv;

}


/**
 * @brief
 *  Reset a timer.
 *
 * @param timer  Timer to reset.
 *
 * @return The time in seconds since 00:00:00 UTC, January 1, 1970 to which
 *   the timer is reset.
 *
 * The function resets the start time of the timer @em timer, leaving the
 * the timer is active or not. The start time to which @em timer is reset
 * is returned as a @c PilTime value for further usage within the caller. In
 * addition, if @em microseconds is not a @c NULL pointer, it is set to the
 * microseconds making up the fraction of the returned time.
 */

PilTime
pilTimerReset(PilTimer *timer, unsigned long *microseconds)
{

  PilTime tv;

  assert(timer != NULL);

  GET_TIME(timer->start);
  tv = _pilTimerEncodeTime(timer, TIMER_END, microseconds);

  return tv;

}


/**
 * @brief
 *  Get the time elapsed since timer start.
 *
 * @param timer  Timer for which the elapsed time is computed.
 *
 * @return The time in seconds since the timer was started.
 *
 * The function computes for the timer @em timer the time inseconds elapsed
 * since the timer was started. The elapsed time is returned as a @c PilTime
 * value. If @em microseconds is not a @c NULL pointer, it is set to the
 * microseconds making up the fraction of the returned time.
 */

PilTime
pilTimerElapsed(PilTimer *timer, unsigned long *microseconds)
{

  unsigned long usec;
  PilTime tv;

#ifdef HAVE_GETTIMEOFDAY
  struct timeval elapsed;
  struct timeval end;
#else
  time_t elapsed;
#endif


  assert(timer != NULL);

  if (timer->active)
    GET_TIME(elapsed);

#ifdef HAVE_GETTIMEOFDAY

  end.tv_sec = timer->end.tv_sec;
  end.tv_usec = timer->end.tv_usec;

  if (timer->start.tv_usec > timer->end.tv_usec) {
    end.tv_usec += 1000000;
    end.tv_sec--;
  }

  elapsed.tv_usec = end.tv_usec - timer->start.tv_usec;
  elapsed.tv_sec = end.tv_sec - timer->start.tv_sec;

  tv = elapsed.tv_sec + ((PilTime)elapsed.tv_usec / 1.e6);
  usec = elapsed.tv_usec;

#else

  tv = timer->end - timer->start;
  usec = 0;

#endif

  if (tv < 0) {
    tv = 0;

    if (microseconds)
      *microseconds = 0;
  }
  else {
    if (microseconds)
      *microseconds = usec;
  }

  return tv;

}


/**
 * @brief
 *  Convert a time value into a date string.
 *
 * @param tvalue  Time value to convert.
 *
 * @return A string of the form YYYY-MM-DD.
 *
 * The function computes the date given by the time value @em tvalue and
 * converts it into a string with the form YYYY-MM-DD.
 * 
 * @note
 *   The result string is written to a static buffer, which will be
 *   overwritten by subsequent calls to this function.
 */

const char *
pilTimerGetDate(PilTime tvalue)
{

    static char tstring[11];

    struct tm *ts;

    time_t seconds = (time_t)tvalue;


    ts = localtime(&seconds);

    sprintf(tstring, "%4d-%02d-%02d",
            ts->tm_year + 1900,
            ts->tm_mon + 1,
            ts->tm_mday);

    return tstring;
    
}


/**
 * @brief
 *  Convert a time value into a time string.
 *
 * @param tvalue  Time value to convert.
 *
 * @return A string of the form hh:mm:ss.SSS.
 *
 * The function computes the time given by the time value @em tvalue and
 * converts it into a string of the form hh:mm:ss.SSS, where sss indicates
 * the fraction of a second.
 * 
 * @note
 *   The result string is written to a static buffer, which will be
 *   overwritten by subsequent calls to this function.
 */

const char *
pilTimerGetTime(PilTime tvalue)
{

    static char tstring[12];

    struct tm *ts;

    time_t seconds = (time_t)tvalue;
    unsigned long milliseconds =  (tvalue - seconds) * 1000.;

    ts = localtime(&seconds);

    sprintf(tstring, "%02d:%02d:%02d.%03ld",
            ts->tm_hour,
            ts->tm_min,
            ts->tm_sec,
            milliseconds);

    return tstring;
    
}


/**
 * @brief
 *  Convert a time value into a calendar time string.
 *
 * @param tvalue  Time value to convert.
 *
 * @return A formatted string as it is created by the @b ctime() function,
 *   but including the milliseconds.
 *
 * The function converts the time value @em tvalue into a string as it is
 * created by the @b ctime() function. The created string also shows the
 * milliseconds.
 * 
 * @note
 *   The result string is written to a static buffer, which will be
 *   overwritten by subsequent calls to this function.
 */

const char *
pilTimerGetCalendarTime(PilTime tvalue)
{

    static char tstring[29];

    struct tm *ts;

    time_t seconds = (time_t)tvalue;
    unsigned long milliseconds = (tvalue - seconds) * 1000.;


    ts = localtime(&seconds);

    sprintf(tstring, "%3s %3s %2d %2d:%02d:%02d.%03ld %4d",
            day[ts->tm_wday],
            month[ts->tm_mon],
            ts->tm_mday,
            ts->tm_hour,
            ts->tm_min,
            ts->tm_sec,
            milliseconds,
            1900 + ts->tm_year);

    return tstring;
    
}


/**
 * @brief
 *  Convert a time value into an ISO8601 date and time.
 *
 * @param tvalue  Time value to convert.
 *
 * @return A string with date and time formatted as specified by the ISO8601
 *   standard.
 *
 * The function converts the time value @em tvalue into a string as it is
 * specified by the ISO8601 standard, i.e. YYYY-MM-DDThh:mm:ss.SSS with
 * SSS representing the fraction of a second.
 * 
 * @note
 *   The result string is written to a static buffer, which will be
 *   overwritten by subsequent calls to this function.
 */

const char *
pilTimerGetTimeISO8601(PilTime tvalue)
{

    static char tstring[24];

    struct tm *ts;

    time_t seconds = (time_t)tvalue;
    unsigned long milliseconds = (tvalue - seconds) * 1000.;


    ts = localtime(&seconds);

    sprintf(tstring, "%4d-%02d-%02dT%02d:%02d:%02d.%03ld",
            ts->tm_year + 1900,
            ts->tm_mon + 1,
            ts->tm_mday,
            ts->tm_hour,
            ts->tm_min,
            ts->tm_sec,
            milliseconds);

    return tstring;
    
}
/**@}*/

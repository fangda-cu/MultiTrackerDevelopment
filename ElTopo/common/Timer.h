/*
 *  Timer.h
 *  CSim
 *
 *  Created by Fang Da on 10/9/11.
 *  Copyright 2011 Columbia. All rights reserved.
 *
 */

#ifndef TIMER_H__
#define TIMER_H__

#if defined(_WIN32) || defined(__CYGWIN__)
#include <windows.h>
#else
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#endif
#include <string>
#include <map>

namespace CSim
{
    //
    //  Usage:
    //
    //  Get absolute time:
    //      double t = Timer::time();
    //      double t = Timer::time(Timer::CPU);
    //      double t = Timer::time(Timer::REAL);
    //
    //  Profiling:
    //      Timer timer(Timer::CPU);    // or other types
    //      timer.start();
    //      // code
    //      timer.stop();
    //      double t = timer.lasttime();
    //
    //  Multiple Profiling:
    //      Timer timer(Timer::CPU);
    //      timer.start();
    //      // profiled code 1
    //      timer.stop();
    //      double t1 = timer.lasttime();
    //      // irrelevant code
    //      timer.start();
    //      // profiled code 2
    //      timer.stop();
    //      double t2 = timer.lasttime();
    //      double t1andt2 = timer.total();
    //
    
    class Timer
    {
    public:
        enum TimerType
        {
            CPU,
            REAL,
            TIMER_TYPE_COUNT
        };
        
    public:    
        Timer(TimerType type = REAL)
        {
            m_type = type;
#if defined(_WIN32) || defined(__CYGWIN__)
            LARGE_INTEGER freq;
            QueryPerformanceFrequency(&freq);
            m_frequency = (double)freq.QuadPart;
#endif
            reset();
        }
        
        ~Timer() 
        { }
        
        inline void reset()
        {
            m_start = 0;
            m_time = 0;
            m_total = 0;
            m_count = 0;
        }
        
        inline void start()
        {
            m_start = time(m_type);
            m_count++;
        }
        
        inline void stop()
        {
            m_time = time(m_type) - m_start;
            m_total += m_time;
        }
        
        inline double lasttime()
        {
            return m_time;
        }
        
        inline double total()
        {
            return m_total;
        }
        
        inline int count()
        {
            return m_count;
        }
        
        inline static double time(TimerType type = REAL)
        {
#if defined(_WIN32) || defined(__CYGWIN__)
            if (type == CPU)
            {
                LARGE_INTEGER query_ticks;
                QueryPerformanceCounter(&query_ticks);
                return query_ticks.QuadPart / m_frequency;
            } else 
            {
                SYSTEMTIME st;
                GetSystemTime(&st);
                return (double)st.wSecond + (double)st.wMilliseconds * 1e-3;
            }
#else 
#if defined(__MACH__) || defined(__APPLE__)
            if (true)    // temporary placeholder implementation, ms precision only and same for CPU and REAL
            {
                struct timeval tv;
                gettimeofday(&tv, NULL);
                return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
            }
#else
            if (type == CPU)
            {
                timespec ts;
                clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
                return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
            } else
            {
                timespec ts;
                clock_gettime(CLOCK_REALTIME, &ts);
                return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
            }
#endif
#endif
        }
        
    protected:
        TimerType m_type;
        
#if defined(_WIN32) || defined(__CYGWIN__)
        double m_frequency;
#endif
        double m_start;
        double m_time;
        double m_total;
        int m_count;
        
    };

    //
    //  Example usage:
    //
    //      while (...)
    //      {
    //          TimerMan::timer("code").start();
    //          // code
    //          TimerMan::timer("code").stop();
    //      }
    //      ...
    //      std::cout << TimerMan::timer("code").total << std::endl;
    //
    class TimerMan
    {
    public:
        static Timer & timer(const std::string & name);
        
        static void setReport(bool r);
        static void report();
        
        TimerMan();
        ~TimerMan();

        static TimerMan * getSingleton();
        
    protected:
        static TimerMan * s_singleton;
        
        bool m_report;
        std::map<std::string, Timer *> m_timers;
        
    };
}

#endif

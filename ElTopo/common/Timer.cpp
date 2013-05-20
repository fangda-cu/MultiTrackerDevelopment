/*
 *  Timer.cpp
 *  CSim
 *
 *  Created by Fang Da on 10/9/11.
 *  Copyright 2011 Columbia. All rights reserved.
 *
 */

#include "Timer.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <stdio.h>

using namespace CSim;

namespace CSim
{
    TimerMan g_timer_singleton;
}

TimerMan * TimerMan::s_singleton = &g_timer_singleton;

TimerMan::TimerMan() :
    m_report(true),
    m_timers()
{
    
}

TimerMan::~TimerMan()
{
    if (m_report && m_timers.size() > 0)
        report();
}

TimerMan * TimerMan::getSingleton()
{
    return s_singleton;
}

Timer & TimerMan::timer(const std::string & name)
{
    std::map<std::string, Timer *> & timers = getSingleton()->m_timers;
    std::map<std::string, Timer *>::iterator i = timers.find(name);
    if (i == timers.end())
    {
        // timer with that name does not exist, create it
        Timer * t = new Timer();
        timers[name] = t;
        return *t;
    } else
    {
        // timer exists
        return *(i->second);
    }
}

void TimerMan::setReport(bool r)
{
    getSingleton()->m_report = r;
}

void TimerMan::report()
{
    std::map<std::string, Timer *> & timers = getSingleton()->m_timers;
    int total_scaling = 1;
    std::string total_unit = "s";
    int mean_scaling = 1;
    std::string mean_unit = "s";
    
    std::vector<std::string> names;
    for (std::map<std::string, Timer *>::iterator i = timers.begin(); i != timers.end(); i++)
        names.push_back(i->first);
    
    // infer time scaling
    std::vector<int> total_digits;
    std::vector<int> mean_digits;
    for (size_t i = 0; i < names.size(); i++)
    {
        Timer & t = timer(names[i]);
        total_digits.push_back(ceil(log10(t.total())));
        mean_digits.push_back(ceil(log10(t.total() / t.count())));
    }
    
    std::sort(total_digits.begin(), total_digits.end());
    int total_digit_median = total_digits[total_digits.size() / 2];
    if (total_digit_median > 0)
        total_scaling = 1, total_unit = "s";
    else if (total_digit_median > -3)
        total_scaling = 1000, total_unit = "ms";
    else if (total_digit_median > -6)
        total_scaling = 1000000, total_unit = "us";
    else
        total_scaling = 1000000000, total_unit = "ns";

    std::sort(mean_digits.begin(), mean_digits.end());
    int mean_digit_median = mean_digits[mean_digits.size() / 2];
    if (mean_digit_median > 0)
        mean_scaling = 1, mean_unit = "s";
    else if (mean_digit_median > -3)
        mean_scaling = 1000, mean_unit = "ms";
    else if (mean_digit_median > -6)
        mean_scaling = 1000000, mean_unit = "us";
    else
        mean_scaling = 1000000000, mean_unit = "ns";
    
    // infer containment relations from timer names based on delimiter '/'
    std::vector<std::vector<std::string> > exploded;
    for (size_t i = 0; i < names.size(); i++)
    {
        std::vector<std::string> tokens;
        std::string token;
        std::stringstream ss(names[i]);
        while(std::getline(ss, token, '/')) 
            tokens.push_back(token);
        exploded.push_back(tokens);
    }
    
    std::vector<std::pair<int, int> > quotients;
    for (size_t i = 0; i < names.size(); i++)
    {
        bool match = false;
        for (size_t j = 0; j < names.size(); j++)
        {
            if (i == j)
                continue;   // skip comparison to self
            
            if (exploded[i].size() + 1 != exploded[j].size())
                continue;   // not immediate parent-child relations

            size_t k = 0;
            for (k = 0; k < exploded[i].size(); k++)
                if (exploded[i][k] != exploded[j][k])
                    break;
            if (k < exploded[i].size())
                continue;   // mismatched ancestors
            
            // parent-child relation confirmed
            quotients.push_back(std::pair<int, int>(i, j));
            match = true;
        }
        if (match)
            quotients.push_back(std::pair<int, int>(i, -1));
    }
    
    // pre-printing, detecting the table dimensions
    size_t maxnamelen = 4;
    size_t maxcountlen = 5;
    size_t maxtotallen = 10;
    size_t maxmeanlen = 9;
    char buffer[100];
    for (size_t i = 0; i < names.size(); i++)
    {
        Timer & t = timer(names[i]);
        maxnamelen = std::max(maxnamelen, names[i].size());
        snprintf(buffer, 99, "%d", t.count());
        maxcountlen = std::max(maxcountlen, strlen(buffer));
        snprintf(buffer, 99, "%8.6f", t.total() * total_scaling);
        maxtotallen = std::max(maxtotallen, strlen(buffer));
        snprintf(buffer, 99, "%8.6f", t.total() / t.count() * mean_scaling);
        maxmeanlen = std::max(maxmeanlen, strlen(buffer));
    }
  
    for (size_t i = 0; i < quotients.size(); i++)
    {
        if (quotients[i].second < 0)
            maxnamelen = std::max(maxnamelen, names[quotients[i].first].size() + 16);
    }
  
    size_t len = maxnamelen + 3 + maxcountlen + 3 + maxtotallen + 3 + maxmeanlen;
    snprintf(buffer, 99, "%%%lud", maxcountlen);
    std::string count_format(buffer);
    snprintf(buffer, 99, "%%%lu.6f", maxtotallen);
    std::string total_format(buffer);
    snprintf(buffer, 99, "%%%lu.6f", maxmeanlen);
    std::string mean_format(buffer);
    
    // print the table
    std::cout << std::string(len, '=') << std::endl;
    std::cout << std::string((len - 15) / 2, ' ') << "TimerMan Report" << std::endl;
    std::cout << std::string(len, '-') << std::endl;
    std::cout << "Name" << std::string(maxnamelen - 4, ' ') << " | ";
    std::cout << "Count" << std::string(maxcountlen - 5, ' ') << " | ";
    std::cout << "Total (" << total_unit << ")" << std::string(maxtotallen - 8 - total_unit.size(), ' ') << " | ";
    std::cout << "Mean (" << mean_unit << ")" << std::string(maxmeanlen - 7 - mean_unit.size(), ' ') << std::endl;
    for (size_t i = 0; i < names.size(); i++)
    {
        Timer & t = timer(names[i]);
        std::cout << names[i] << std::string(maxnamelen - names[i].size(), ' ') << " | ";
        snprintf(buffer, 99, count_format.c_str(), t.count());
        std::cout << buffer << " | ";
        snprintf(buffer, 99, total_format.c_str(), t.total() * total_scaling);
        std::cout << buffer << " | ";
        snprintf(buffer, 99, mean_format.c_str(), t.total() / t.count() * mean_scaling);
        std::cout << buffer << std::endl;
    }
    std::cout << std::string(len, '-') << std::endl;
    std::cout << "Name" << std::string(maxnamelen - 4, ' ') << " | ";
    std::cout << "Percentage" << std::string(maxcountlen + 3 + maxtotallen - 10, ' ') << std::endl;
    for (size_t i = 0; i < quotients.size(); i++)
    {
        double childtime = 0;
        std::string childname = "";
        if (quotients[i].second < 0)
        {
            double sumtime = 0;
            for (size_t j = 0; j < quotients.size(); j++)
                if (quotients[j].first == (int)quotients[i].first && quotients[j].second >= 0)
                    sumtime += timer(names[quotients[j].second]).total();
            childtime = timer(names[quotients[i].first]).total() - sumtime;
            childname = names[quotients[i].first] + "/unaccounted_for";
        } else
        {
            childtime = timer(names[quotients[i].second]).total();
            childname = names[quotients[i].second];
        }
        snprintf(buffer, 99, "%6.2f%%", 100.0 * childtime / timer(names[quotients[i].first]).total());
        std::cout << childname << std::string(maxnamelen - childname.size(), ' ') << " | ";
        std::cout << buffer << std::endl;
    }
    std::cout << std::string(len, '=') << std::endl;

}

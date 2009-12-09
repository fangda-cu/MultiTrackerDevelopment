// $Id: /wmMisc/tags/wmMisc-6-2/WsplineAttr.hh 579 2006-07-07T06:24:07.165057Z phunter  $
//
// $HeadURL: /wmMisc/tags/wmMisc-6-2/WsplineAttr.hh $
// $LastChangedDate: 2006-07-07T06:24:07.165057Z $
// $LastChangedBy: phunter $
//
// Copyright 2003 Weta Digital Limited
// Jeff Hameluck
//

#ifndef WFOZSPLINEATTR_H_
#define WFOZSPLINEATTR_H_

//#include "WFozCore.h"

#include <inttypes.h>
#include <map>

typedef std::map<float, std::pair<float, int16_t> > cvDataMap;


class WsplineAttr 
{
public: 
	// Methods
	~WsplineAttr() {}

	enum interpType 
    {
        kLinear = 0, 
        kNone = 1, 
        kSpline = 2, 
        kSmooth = 3
    };  // See Maya's MRampAttribute.h
	
	void set(cvDataMap &cvMap);
	float getValue(float param);
	/*
protected: // Methods
public: // Data

protected:  // Data*/
	cvDataMap cvData;
};

#endif // WFOZSPLINEATTR_H_

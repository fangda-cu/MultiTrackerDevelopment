// $Id: /wmMisc/tags/wmMisc-6-2/WsplineAttr.cc 579 2006-07-07T06:24:07.165057Z phunter  $
//
// $HeadURL: /wmMisc/tags/wmMisc-6-2/WsplineAttr.cc $
// $LastChangedDate: 2006-07-07T06:24:07.165057Z $
// $LastChangedBy: phunter $
//
// Copyright � 2003 Weta Digital Limited
// Jeff Hameluck
//

//#include <weta/Wiostream.hh>
#include "SplineAttrEval.hh"


void WsplineAttr::set(cvDataMap &cvMap)
{
	cvData = cvMap;
	cvDataMap::iterator cvIter;
#if 0
	for (cvIter = cvData.begin(); cvIter != cvData.end(); cvIter++)
	{
		fprintf(stderr, "Position:%f \t Value: %f \t Interp: ", cvIter->first, cvIter->second.first);
		switch (cvIter->second.second)
		{
			case LINEAR:
				fprintf(stderr, "LINEAR\n");
				break;
			case SMOOTH:
				fprintf(stderr, "SMOOTH\n");
				break;
			case SPLINE:
				fprintf(stderr, "SPLINE\n");
				break;
			case NONE:
				fprintf(stderr, "NONE\n");
				break;
			default:
				fprintf(stderr, "UNKNOWN\n");
				break;
		}
	}
#endif
}



float WsplineAttr::getValue(float param)
{
	float value = 0;
	
	cvDataMap::const_iterator cvIter = cvData.upper_bound(param);
	
	if (cvIter == cvData.end())  		// at the end
  	{
		if (cvData.size())  			// get the value of the last CV, otherwise will ret 0
	  	{
			--cvIter;
		  	value = cvIter->second.first;
    	}
  	}
	else if (cvIter == cvData.begin())  // before the beginning
	{
		value = cvIter->second.first;
	} 
	else 
	{
    	float k1 = cvIter->second.first;
    	float t1 = cvIter->first;
    	--cvIter;
    	float k0 = cvIter->second.first;
    	float t0 = cvIter->first;
    	
		switch (cvIter->second.second) 
		{
			case WsplineAttr::kNone:
      			value = k0;
      			break;
    		case WsplineAttr::kLinear:
				{
					float u = (param - t0) / (t1 - t0);
					value = k0 + u*(k1 - k0);
      			}
      			break;
    		case WsplineAttr::kSmooth:
      			{
					float u = (param - t0) / (t1 - t0);
					value = k0*(2*u*u*u - 3*u*u + 1) + k1*(3*u*u - 2*u*u*u);
      			}
      			break;
    		case WsplineAttr::kSpline:
      			{
					float scale = 1.0f / (t1 - t0);
					float u = (param - t0) * scale;
					float k0p = 0;
					if (cvIter != cvData.begin()) 
					{
	  					--cvIter;
	  					k0p = .5 * ((k1 - cvIter->second.first) / (scale * (t1 - cvIter->first)));
	  					++cvIter;
					}
					++cvIter;
					++cvIter;
					float k1p = 0;
					if (cvIter != cvData.end()) 
					{
	  					k1p = .5 * ((cvIter->second.first - k0) / (scale * (cvIter->first - t0)));
					}
					
					value = k0*(2*u*u*u - 3*u*u + 1) + k1*(3*u*u - 2*u*u*u) +
	  				k0p*(u*u*u - 2*u*u + u) + k1p*(u*u*u - u*u);
      			}
      			break;
    		default:
      			break;
    	}
  	}

  return value;
}


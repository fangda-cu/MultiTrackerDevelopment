/*
 * ElasticStrandParameters.hh
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#ifndef ELASTICSTRANDPARAMETERS_HH_
#define ELASTICSTRANDPARAMETERS_HH_

namespace strandsim
{

class ElasticStrandParameters
{
public:
    ElasticStrandParameters() :
        m_radius(1.0), m_YoungsModulus(1.0)
    {
        setup();
    }

    ElasticStrandParameters(const ElasticStrandParameters& other) :
        m_radius(other.m_radius), m_YoungsModulus(other.m_YoungsModulus)
    {
        setup();
    }

    void setup()
    {
        m_ks = M_PI * m_radius * m_radius * m_YoungsModulus;
    }

    // Physical parameters. For now these are all constant along the rod
    const Scalar m_radius;
    const Scalar m_YoungsModulus;

    // Computed parameters. Make sure setup() is called each time any of the above is changed.
    Scalar m_ks;
};

}

#endif /* ELASTICSTRANDPARAMETERS_HH_ */

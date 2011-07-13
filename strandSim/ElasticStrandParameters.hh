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
        m_radius(1.0), m_YoungsModulus(1.0), m_shearModulus(1.0) // Is it really a good idea to have defaults?
    {
        setup();
    }

    ElasticStrandParameters(const ElasticStrandParameters& other) :
        m_radius(other.m_radius), m_YoungsModulus(other.m_YoungsModulus), m_shearModulus(other.m_shearModulus)
    {
        setup();
    }

    void setup()
    {
        m_ks = M_PI * m_radius * m_radius * m_YoungsModulus;
        m_kt = 0.5 * M_PI * square(square(m_radius)) * m_shearModulus;
    }

    // Physical parameters. For now these are all constant along the rod
    const Scalar m_radius;
    const Scalar m_YoungsModulus;
    const Scalar m_shearModulus;

    // Computed parameters. Make sure setup() is called each time any of the above is changed.
    Scalar m_ks;
    Scalar m_kt;
};

}

#endif /* ELASTICSTRANDPARAMETERS_HH_ */

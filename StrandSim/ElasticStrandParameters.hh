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
private:
    ElasticStrandParameters();

public:
    ElasticStrandParameters( Scalar radiusA, Scalar radiusB, Scalar YoungsModulus,
            Scalar shearModulus, Scalar density, Scalar baseRotation = 0.0 ) :
        m_radiusA( radiusA ), m_radiusB( radiusB ), m_YoungsModulus( YoungsModulus ),
                m_shearModulus( shearModulus ), m_density( density ), m_baseRotation( baseRotation )
    {
        setup();
    }

    ElasticStrandParameters( const ElasticStrandParameters& other ) :
        m_radiusA( other.m_radiusA ), m_radiusB( other.m_radiusB ),
                m_YoungsModulus( other.m_YoungsModulus ), m_shearModulus( other.m_shearModulus ),
                m_density( other.m_density ), m_baseRotation( other.m_baseRotation )
    {
        setup();
    }

    void setup()
    {
        m_ks = M_PI * m_radiusA * m_radiusB * m_YoungsModulus;
        m_kt = 0.25 * M_PI * m_radiusA * m_radiusB * ( square( m_radiusA ) + square( m_radiusB ) )
                * m_shearModulus;
    }

    // Physical parameters. For now these are all constant along the rod
    //    const Scalar m_radius;
    const Scalar m_radiusA;
    const Scalar m_radiusB;
    const Scalar m_YoungsModulus;
    const Scalar m_shearModulus;
    const Scalar m_density;
    const Scalar m_baseRotation;

    // Computed parameters. Make sure setup() is called each time any of the above is changed.
    Scalar m_ks;
    Scalar m_kt;
};

}

#endif /* ELASTICSTRANDPARAMETERS_HH_ */

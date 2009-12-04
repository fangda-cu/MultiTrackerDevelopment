#include "RodData.hh"

RodData::RodData() : shouldSimulate( true )
{
    rod = NULL; stepper = NULL; rodRenderer = NULL; 
}

RodData::RodData( ElasticRod* i_rod, RodCollisionTimeStepper* i_stepper, RodRenderer* i_rodRenderer ) :
    shouldSimulate( true )
{
    rod = i_rod; stepper = i_stepper; rodRenderer = i_rodRenderer; 
}

RodData::~RodData()
{
    if ( rod != NULL )
        delete rod;
    if ( stepper != NULL )
        delete stepper;
    if ( rodRenderer != NULL )
        delete rodRenderer;
}

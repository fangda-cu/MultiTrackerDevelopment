#ifndef WMFIGRODPARTICLEINPUT_H_
#define WMFIGRODPARTICLEINPUT_H_

/**
  * @file Apps/WmFigaro/WmFigRodParticleInput.h
  * @author Alasdair Coull (acoull@wetafx.co.nz)
  * @date 30/11/10
  */

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MVectorArray.h>
#include <maya/MIntArray.h>
#include <maya/MFnIntArrayData.h>

#include <set>

#include "WmFigRodInputType.hh"

#include "WmFigRodGroup.hh"

/**
  * @brief Class to convert nParticle inputs into rods
  *
  *  
  *
  *
  */ 

class WmFigRodParticleInput : public WmFigRodInputType
{ 
public:
     /**
      * @brief Default constructor
      */
     WmFigRodParticleInput( MObject& i_verticesAttribute, MObject& i_perRodParticleCountAttribute,
                        RodOptions& i_rodOptions,
                        double i_massDamping, BASim::Vec3d& i_gravity, WmFigRodGroup& i_rodGroup,
                        RodTimeStepper::Method i_solverType, std::set< int >& i_simulationSet );

    /**
      * @brief Default destructor
      */
     ~WmFigRodParticleInput();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock );
    virtual void updateRodDataFromInput( MDataBlock& i_dataBlock );
    virtual int numberOfInputs( MDataBlock& i_dataBlock );

private:
    void getRodVertices( MDataBlock& i_dataBlock, MVectorArray& o_vertices, 
                         MIntArray& o_perRodParticleCount );

    MObject& m_verticesAttribute;
    MObject& m_perRodParticleCountAttribute;
    
    WmFigRodGroup& m_rodGroup;

    RodOptions m_rodOptions;
    double m_massDamping;
    BASim::Vec3d m_gravity;
    RodTimeStepper::Method m_solverType;
    std::set< int >& m_simulationSet;
};

#endif // WMFIGRODPARTICLEINPUT_H_

 


#ifndef WMFIGRODBARBINPUT_H_
#define WMFIGRODBARBINPUT_H_

/**
  * @file Apps/WmFigaro/WmFigRodBarbInput.h
  * @author Alasdair Coull (acoull@wetafx.co.nz)
  * @date 20-04-2010
  */

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MVectorArray.h>

#include <set>

#include "WmFigRodInputType.hh"

#include "WmFigRodGroup.hh"

/**
  * @brief Class to convert Barbershop strands into a format suitable for creating rods
  *
  *  
  *
  *
  */ 

class WmFigRodBarbInput : public WmFigRodInputType
{ 
public:
     /**
      * @brief Default constructor
      */
     WmFigRodBarbInput( MObject& i_verticesAttribute, MObject& i_strandRootFramesAttribute, 
                        double i_percentageOfBarbStrands, size_t i_verticesPerRod,
                        bool i_lockFirstEdgeToInput, double i_vertexSpacing,
                        double i_minimumRodLength, RodOptions& i_rodOptions,
                        double i_massDamping, WmFigRodGroup& i_rodGroup,
                        RodTimeStepper::Method i_solverType, std::set< size_t >& i_simulationSet );

    /**
      * @brief Default destructor
      */
     ~WmFigRodBarbInput();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock );
    virtual void updateRodDataFromInput( MDataBlock& i_dataBlock );
    virtual size_t numberOfInputs( MDataBlock& i_dataBlock );

private:
    void getStrandRootFrames( MDataBlock& i_dataBlock, std::vector<MaterialFrame>& o_strandRootFrames );
    void getStrandVertices( MDataBlock& i_dataBlock, MVectorArray& o_vertices, size_t& o_numStrands );

    MObject& m_verticesAttribute;
    MObject& m_strandRootFramesAttribute;
    double m_percentageOfBarbStrands;
    size_t m_verticesPerRod;
    bool m_lockFirstEdgeToInput;
    double m_vertexSpacing;
    double m_minimumRodLength;

    WmFigRodGroup& m_rodGroup;

    RodOptions m_rodOptions;
    double m_massDamping;
    RodTimeStepper::Method m_solverType;
    std::set< size_t >& m_simulationSet;
};

#endif // WMFIGRODNURBSINPUT_H_

 

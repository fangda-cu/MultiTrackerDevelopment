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
                        double i_percentageOfBarbStrands, int i_verticesPerRod,
                        bool i_lockFirstEdgeToInput, double i_vertexSpacing,
                        double i_minimumRodLength, RodOptions& i_rodOptions,
                        double i_massDamping, BASim::Vec3d& i_gravity, WmFigRodGroup& i_rodGroup,
                        RodTimeStepper::Method i_solverType, std::set< int >& i_simulationSet );

    /**
      * @brief Default destructor
      */
     ~WmFigRodBarbInput();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock );
    virtual void updateRodDataFromInput( MDataBlock& i_dataBlock );
    virtual int numberOfInputs( MDataBlock& i_dataBlock );

private:
    void getStrandRootFrames( MDataBlock& i_dataBlock, std::vector<MaterialFrame>& o_strandRootFrames );
    void getStrandVertices( MDataBlock& i_dataBlock, MVectorArray& o_vertices, int& o_numStrands );

    MObject& m_verticesAttribute;
    MObject& m_strandRootFramesAttribute;
    double m_percentageOfBarbStrands;
    int m_verticesPerRod;
    bool m_lockFirstEdgeToInput;
    double m_vertexSpacing;
    double m_minimumRodLength;

    WmFigRodGroup& m_rodGroup;

    RodOptions m_rodOptions;
    double m_massDamping;
    BASim::Vec3d m_gravity;
    RodTimeStepper::Method m_solverType;
    std::set< int >& m_simulationSet;
};

#endif // WMFIGRODNURBSINPUT_H_

 

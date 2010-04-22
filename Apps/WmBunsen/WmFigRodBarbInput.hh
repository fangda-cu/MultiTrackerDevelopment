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

#include "WmFigRodInputType.hh"

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
                        bool i_lockFirstEdgeToInput );

    /**
      * @brief Default destructor
      */
     ~WmFigRodBarbInput();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock, std::vector<RodData*>* i_pRodData  );
    virtual void updateRodDataFromInput( MDataBlock& i_dataBlock, std::vector<RodData*>* i_pRodData );
    virtual size_t numberOfInputs( MDataBlock& i_dataBlock );

private:
    void getStrandRootFrames( MDataBlock& i_dataBlock, std::vector<MaterialFrame>& o_strandRootFrames );
    void getStrandVertices( MDataBlock& i_dataBlock, MVectorArray& o_vertices, size_t& o_numStrands );

    MObject& m_verticesAttribute;
    MObject& m_strandRootFramesAttribute;
    double m_percentageOfBarbStrands;
    size_t m_verticesPerRod;
    bool m_lockFirstEdgeToInput;
};

#endif // WMFIGRODNURBSINPUT_H_

 

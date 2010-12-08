#ifndef WMFIGRODINPUTTYPE_H_
#define WMFIGRODINPUTTYPE_H_

/**
  * @file Apps/WmFigaro/WmFigRodInputType.h
  * @author Alasdair Coull (acoull@wetafx.co.nz)
  * @date 20-04-2010
  */

#include <maya/MDataBlock.h>
#include <maya/MVector.h>

#include "RodData.hh"
#include <vector>

/**
  * @brief Base class for classes which convert Maya types into Rods structures.
  *
  * This class should be derived from when you want to allow a new Maya object type to be used to
  * drive rods. For example, NURBS curves or Barbershop nodes.
  *
  */ 

class WmFigRodInputType
{ 
public:
     /**
      * @brief Default constructor
      */
     WmFigRodInputType();

    /**
      * @brief Default destructor
      */
     virtual ~WmFigRodInputType();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock ) = 0;
    virtual void updateRodDataFromInput( MDataBlock& i_dataBlock );
    virtual int numberOfInputs( MDataBlock& i_dataBlock ) { return 0; };
    void resampleCurve( int i_numVerticesToResample, vector<MVector>& i_curve, 
                        vector<MVector>& o_resampledCurve );


    // may be helpful for the input class to know if sim has been turned off
    void setSimulating(bool);
    bool getSimulating();

protected:

    bool m_simulating;

private:
};

#endif // WMFIGRODINPUTTYPE_H_

 

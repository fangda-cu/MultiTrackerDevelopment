#ifndef WMFIGRODINPUTTYPE_H_
#define WMFIGRODINPUTTYPE_H_

/**
  * @file Apps/WmFigaro/WmFigRodInputType.h
  * @author Alasdair Coull (acoull@wetafx.co.nz)
  * @date 20-04-2010
  */

#include <maya/MDataBlock.h>

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
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock, std::vector<RodData*>* i_pRodData  ) = 0;
    virtual void updateRodDataFromInput( MDataBlock& i_dataBlock, std::vector<RodData*>* i_pRodData );
    virtual size_t numberOfInputs( MDataBlock& i_dataBlock ) { return 0; };

protected:
private:
};

#endif // WMFIGRODINPUTTYPE_H_

 

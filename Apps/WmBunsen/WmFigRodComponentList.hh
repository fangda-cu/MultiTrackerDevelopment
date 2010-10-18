/*
 * WmFigRodComponentList.hh
 *
 *  Created on: Oct 19, 2010
 *      Author: dgould
 */

#ifndef WMFIGRODCOMPONENTLIST_HH
#define WMFIGRODCOMPONENTLIST_HH

#include <maya/MIntArray.h>
#include <maya/MString.h>
#include <map>

class WmFigRodComponentList
{
public:
	bool containsRodVertex( const unsigned int rodId, const unsigned int vertexId );
	void addOrRemoveRodVertex( const unsigned int rodId, const unsigned int vertexId, const bool doAdd=true );
    void removeAllRodVertices();

    void getRodIds( MIntArray &rodIds ) const;
    void getRodVertexIds( const unsigned int rodId, MIntArray &rodVertexIds ) const;

    bool serialise( MString &toString ) const;
    bool unserialise( MString &fromString );

private:
	typedef MIntArray VertexIndices;
    typedef std::map<unsigned int,VertexIndices> RodVertexIndices;

    RodVertexIndices rodVertexIndices;
};

#endif

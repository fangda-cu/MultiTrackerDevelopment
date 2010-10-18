/*
 * WmFigRodComponentList.cc
 *
 *  Created on: Oct 19, 2010
 *      Author: dgould
 */
#include "WmFigRodComponentList.hh"
#include <maya/MStringArray.h>
#include <maya/MGlobal.h>

bool WmFigRodComponentList::containsRodVertex( const unsigned int rodId, const unsigned int vertexId )
{
	unsigned int key = rodId;

	RodVertexIndices::iterator it;
	it = rodVertexIndices.find( key );
	if( it != rodVertexIndices.end() ) { // Found matching rod
		MIntArray &rodVertices = it->second;
		if( vertexId < rodVertices.length() )
			return bool( rodVertices[ vertexId ] );
	}

	return false;
}

void WmFigRodComponentList::addOrRemoveRodVertex( const unsigned int rodId, const unsigned int vertexId, const bool doAdd )
{
	unsigned int key = rodId;

	RodVertexIndices::iterator it;
	it = rodVertexIndices.find( key );
	if( it == rodVertexIndices.end() ) { // Didn't find matching rod
		std::pair<RodVertexIndices::iterator, bool> result;
		result = rodVertexIndices.insert( std::pair<unsigned int,VertexIndices>( key, VertexIndices() ) );
		it = result.first;
	}

	MIntArray &rodVertices = it->second;
	if( vertexId >= rodVertices.length() )  {
		const unsigned int prevLength = rodVertices.length();

		// Resize
		rodVertices.setLength( vertexId+1 );

		// Set all the newly added values to 0 (not selected)
		unsigned int i;
		for( i=prevLength; i < rodVertices.length(); i++ )
			rodVertices[i] = 0; // Not included
	}

	rodVertices[ vertexId ] = int( doAdd );
}


void WmFigRodComponentList::removeAllRodVertices()
{
	unsigned int vi;
	RodVertexIndices::iterator it;
	for( it=rodVertexIndices.begin(); it != rodVertexIndices.end(); it++ ) {
		MIntArray &rodVertices = (*it).second;
		for( vi=0; vi < rodVertices.length(); vi++ )
			rodVertices[vi] = 0; // Not selected
	}
}

void WmFigRodComponentList::getRodIds( MIntArray &rodIds ) const
{
	rodIds.clear();

	unsigned int vi;
	RodVertexIndices::const_iterator it;
	for( it=rodVertexIndices.begin(); it != rodVertexIndices.end(); it++ ) {
		const MIntArray &rodVertices = (*it).second;
		for( vi=0; vi < rodVertices.length(); vi++ ) {
			if( rodVertices[vi] ) {
				rodIds.append( (*it).first );
				break;
			}
		}
	}
}

void WmFigRodComponentList::getRodVertexIds( const unsigned int rodId, MIntArray &rodVertexIds ) const
{
	rodVertexIds.clear();

	unsigned int key = rodId;
	RodVertexIndices::const_iterator it;
	it = rodVertexIndices.find( key );
	if( it != rodVertexIndices.end() ) { // Found matching rod
		const MIntArray &rodVertices = (*it).second;
		unsigned int vi;
		for( vi=0; vi < rodVertices.length(); vi++ ) {
			if( rodVertices[vi] )
				rodVertexIds.append( vi );
		}
	}
}

/*
 * History
 * --------------------------------
 * 1 - First version
 */

const unsigned int WM_FIGRODCOMPONENTLIST_DATAVERSION = 1;

bool WmFigRodComponentList::serialise( MString &toString ) const
{
	toString.clear();

	toString += WM_FIGRODCOMPONENTLIST_DATAVERSION;

	MIntArray rodIds;
	MIntArray rodVertexIds;

	getRodIds( rodIds );

	unsigned int iVertex;
	unsigned int iRod;
	unsigned int rodId, vertexId;
	for( iRod=0; iRod < rodIds.length(); iRod++ ) {
		rodId = rodIds[iRod];

		toString += MString(" rod[") + rodId + "].vtx[";

		getRodVertexIds( rodId, rodVertexIds );
		for( iVertex=0; iVertex < rodVertexIds.length(); iVertex++ ) {
			vertexId = rodVertexIds[ iVertex ];

			if( iVertex > 0 )
				toString += ",";
			toString += vertexId;
		}
		toString += MString("]");
	}

	return true;
}

bool WmFigRodComponentList::unserialise( MString &fromString )
{
	if( !fromString.length() )
		return false;

	MStringArray stringElements;
	fromString.split( ' ', stringElements );

	unsigned int dataVersion;
	dataVersion = stringElements[0].asInt();

	if( dataVersion > WM_FIGRODCOMPONENTLIST_DATAVERSION ) {
		MGlobal::displayWarning( "WmFigRodComponentList data is a more recent than this version supports.");
		return false;
	}

	// @@@ to do...extract the rest

	return true;
}



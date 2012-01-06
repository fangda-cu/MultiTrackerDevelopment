/*
 * ObjWriter.cc
 *
 *  Created on: Dec 30, 2011
 *      Author: andres
 */

#include "ObjWriter.hh"
#include "../Physics/DeformableObjects/DeformableObject.hh"
#include "../Core/TopologicalObject/TopObjIterators.hh"
#include "../Core/TopologicalObject/TopObjProperty.hh"
#include "../Core/Definitions.hh"
#include "../Physics/DeformableObjects/Shells/ElasticShell.hh"

namespace BASim
{

ObjWriter::ObjWriter()
{
    // TODO Auto-generated constructor stub

}

ObjWriter::~ObjWriter()
{
    // TODO Auto-generated destructor stub
}
void ObjWriter::write( const std::string & filename, const ElasticShell& mesh )
{
    std::ofstream of( filename.c_str() );
    ObjWriter::write( of, mesh );
}
void ObjWriter::write( std::ofstream & of, const VertexProperty<Vec3d> & pos,
        DeformableObject & mesh )
{
    VertexProperty<int> indices( &mesh );
    int current = 0;
    for ( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit )
    {
        const Vec3d v = pos[*vit];
        of << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
        indices[*vit] = ++current;
    }

    //write face data
    for ( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
        of << "f";
        //It should be like: "f v"
        for ( FaceVertexIterator fvit = mesh.fv_iter( *fit ); fvit; ++fvit )
        {
            of << " ";
            of << indices[*fvit]; //v
        }
        of << std::endl;
    }
}
void ObjWriter::write( std::ofstream & of, const VertexProperty<Vec3d> & pos,
        const VertexProperty<Vec3d> & vNormals, DeformableObject & mesh )
{
    VertexProperty<int> indices( &mesh );
    int current = 0;
    for ( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit )
    {
        const Vec3d v = pos[*vit];
        of << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
        indices[*vit] = ++current;
    }

    //write normals
    VertexProperty<int> nIndices( &mesh );
    current = 0;
    for ( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit )
    {
        of << "vn " << vNormals[*vit].x() << " " << vNormals[*vit].y() << " " << vNormals[*vit].z()
                << std::endl;
        nIndices[*vit] = ++current;
    }

    //write face data
    for ( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
        of << "f";
        //It should be like: "f v/vt/vn"
        for ( FaceVertexIterator fvit = mesh.fv_iter( *fit ); fvit; ++fvit )
        {
            of << " ";
            of << indices[*fvit]; //v
            of << "//"; // No texture
            of << nIndices[*fvit]; //vn

        }
        of << std::endl;
    }
}
void ObjWriter::write( std::ofstream & of, const ElasticShell & shell )
{

    const DeformableObject &mesh = shell.getDefoObj();
    std::vector<VertexHandle> vhandles;

    FaceProperty<Vec3d> faceNormals( &shell.getDefoObj() );
    VertexProperty<Vec3d> vertexNormals( &shell.getDefoObj() );
    VertexProperty<Scalar> vertThickness( &shell.getDefoObj() );

    DeformableObject defObj;

    shell.getFaceNormals( faceNormals );
    shell.getVertexNormals( vertexNormals );
    shell.getThickness( vertThickness );

    VertexProperty<int> vhIndices( &shell.getDefoObj() );
    VertexProperty<Vec3d> vertPositions( &defObj );
    VertexProperty<int> indices (&defObj);
    int current = 0;
    int currentH = 1; //OBJ vertex numbering starts from 1

    //Push all the vertex handles
    for ( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit )
    {
        //Add top and bottom vertices
        VertexHandle vht = defObj.addVertex();
        VertexHandle vhb = defObj.addVertex();

        vhandles.push_back( vht );
        vhandles.push_back( vhb );
//        //Assign the correct positions to them
//        vertPositions[vht] =
        vhIndices[*vit] = current++;

        vertPositions[vht] = shell.getVertexPosition( *vit )
                + vertexNormals[*vit] * vertThickness[*vit] * 0.5;
        vertPositions[vhb] = shell.getVertexPosition( *vit )
                - vertexNormals[*vit] * vertThickness[*vit] * 0.5;

        indices[vht] = currentH++;
        indices[vhb] = currentH++;

    }

    //Now add all the faces for top and bottom layers
    for ( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
        std::vector<VertexHandle> triT;
        std::vector<VertexHandle> triB;
        for ( FaceVertexIterator fvit = mesh.fv_iter( *fit ); fvit; ++fvit )
        {
            triT.push_back( vhandles[2 * vhIndices[*fvit]] );
            triB.push_back( vhandles[2 * vhIndices[*fvit] + 1] );
        }
        defObj.addFace( triT[2], triT[1], triT[0] );
        defObj.addFace( triB[0], triB[1], triB[2] );

    }

    ObjWriter::write( of, vertPositions, defObj );
    for ( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
        //Now the stich faces located on the boundaries of the mesh
        for ( FaceEdgeIterator feit = mesh.fe_iter( *fit ); feit; ++feit )
        {
            if ( mesh.edgeIncidentFaces( *feit ) == 1 )
            {
                int orient = mesh.getRelativeOrientation( *fit, *feit );
                int from, to;

                if ( orient == 1 )
                {
                    from = vhIndices[mesh.fromVertex( *feit )];
                    to = vhIndices[mesh.toVertex( *feit )];
                }
                else
                {
                    from = vhIndices[mesh.toVertex( *feit )];
                    to = vhIndices[mesh.fromVertex( *feit )];
                }
                int tl = indices[vhandles[2 * from]];
                int bl = indices[vhandles[2 * from + 1]];
                int tr = indices[vhandles[2 * to]];
                int br = indices[vhandles[2 * to + 1]];


                of << "f " << tl << " " << bl << " " << br << " " << tr << std::endl;
//                defObj.addFace( tl, bl, br );
//                defObj.addFace( tl, br, tr );

            }
        }
    }
}

}/* namespace BASim */


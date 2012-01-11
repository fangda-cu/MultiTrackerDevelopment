/*
 * PlyWriter.cc
 *
 *  Created on: Jan 10, 2012
 *      Author: andres
 */

#include "PlyWriter.hh"

#include <iostream>
#include <fstream>
#include "../Physics/DeformableObjects/DeformableObject.hh"
#include "../Core/TopologicalObject/TopObjIterators.hh"
#include "../Core/TopologicalObject/TopObjProperty.hh"
#include "../Core/Definitions.hh"
namespace BASim
{

PlyWriter::PlyWriter()
{
    // TODO Auto-generated constructor stub

}

PlyWriter::~PlyWriter()
{
    // TODO Auto-generated destructor stub
}

void PlyWriter::write(const std::string & filename, const ElasticShell & mesh)
{
//    ply
//    format ascii 1.0
//    comment Created by baRods, a Grinspun Labs simulation tool
//    element vertex 278
//    property float x
//    property float y
//    property float z
//    property float w
//    property float width
//    element face 0
//    property list int int vertex_indices
//    end_header
    std::ofstream of (filename.c_str());
    of << "ply\nformat ascii 1.0\ncomment Created by BASim, a Grinspun Labs simulation tool\n";
    of << "element vertex " << mesh.getDefoObj().nv() << std::endl;
    of << "property float x" << std::endl;
    of << "property float y" << std::endl;
    of << "property float z" << std::endl;
    of << "property float thickness" << std::endl;

    of << "element face " << mesh.getDefoObj().nf() << std::endl;
    of << "property list int int vertex_index" << std::endl;
    of << "end_header" << std::endl;

    DeformableObject & dobj = mesh.getDefoObj();
    VertexProperty<int> indices( &dobj );

    int current = 0;
    for ( VertexIterator vit = dobj.vertices_begin(); vit != dobj.vertices_end(); ++vit )
    {
        const Vec3d v = mesh.getVertexPosition(*vit);
        of << v.x() << " " << v.y() << " " << v.z();
        //Now the thickness
        of << " " << mesh.getThickness(*vit) << std::endl;
        indices[*vit] = current++;
    }
    //write face data
    for ( FaceIterator fit = dobj.faces_begin(); fit != dobj.faces_end(); ++fit )
    {
        //First the property list
        //Each face points to three ints
        of << "3";
        for ( FaceVertexIterator fvit = dobj.fv_iter( *fit ); fvit; ++fvit )
        {
            of << " ";
            of << indices[*fvit];
        }
        of << std::endl;
    }

    of.close();

}

} /* namespace BASim */


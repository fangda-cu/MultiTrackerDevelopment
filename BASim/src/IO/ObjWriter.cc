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
void ObjWriter::write(const std::string & filename, const ElasticShell& mesh){
    std::ofstream of ( filename.c_str() );
    ObjWriter::write(of, mesh);
}
void ObjWriter::write(std::ofstream & of, const VertexProperty<Vec3d> & pos, const VertexProperty<Vec3d> & vNormals, DeformableObject & mesh){
    VertexProperty<int> indices(&mesh);
    int current = 0;
    for (VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit){
        const Vec3d v = pos[*vit];
        of << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
        indices[*vit] = ++current;
    }

    //write normals
    VertexProperty<int> nIndices(&mesh);
    current = 0;
    for ( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit){
        of << "vn " << vNormals[*vit].x() << " " << vNormals[*vit].y() << " " << vNormals[*vit].z() << std::endl;
        nIndices[*vit] = ++current;
    }

    //write face data
    for ( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit){
        of << "f";
        //It should be like: "f v/vt/vn"
        for (FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit){
            of << " ";
            of << indices[*fvit]; //v
            of << "//"; // No texture
            of << nIndices[*fvit];//vn

        }
        of << std::endl;
    }
}
void ObjWriter::write(std::ofstream & of, const ElasticShell & shell){
    VertexProperty<Vec3d> vNormals (&shell.getDefoObj());
    shell.getVertexNormals(vNormals);
    ObjWriter::write(of, shell.getVertexPositions(), vNormals, shell.getDefoObj());
}

}/* namespace BASim */

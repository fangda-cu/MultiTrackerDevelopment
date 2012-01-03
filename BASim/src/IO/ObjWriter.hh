/*
 * ObjWriter.hh
 *
 *  Created on: Dec 30, 2011
 *      Author: andres
 */

#ifndef OBJWRITER_HH_
#define OBJWRITER_HH_

#include "../Physics/DeformableObjects/Shells/ElasticShell.hh";
#include <iostream>
#include <fstream>

namespace BASim
{

class ObjWriter
{
public:
    ObjWriter();
    virtual ~ObjWriter();

    static void write(std::ofstream & of, const ElasticShell& mesh);
    static void write(const std::string & filename, const ElasticShell& mesh);
    static void write(std::ofstream & of, const VertexProperty<Vec3d> & pos, DeformableObject & mesh);
    static void write(std::ofstream & of, const VertexProperty<Vec3d> & pos, const VertexProperty<Vec3d> & normals, DeformableObject & mesh);
    friend class VertexProperty<Vec3d>;
};

} /* namespace BASim */
#endif /* OBJWRITER_HH_ */

/*
 * PlyWriter.hh
 *
 *  Created on: Jan 10, 2012
 *      Author: andres
 */

#ifndef PLYWRITER_HH_
#define PLYWRITER_HH_

#include "../Physics/DeformableObjects/Shells/ElasticShell.hh"

#include <string>

namespace BASim
{

class PlyWriter
{
public:
    PlyWriter();
    virtual ~PlyWriter();

    static void write(const std::string & filename, const ElasticShell& mesh);
};

} /* namespace BASim */
#endif /* PLYWRITER_HH_ */

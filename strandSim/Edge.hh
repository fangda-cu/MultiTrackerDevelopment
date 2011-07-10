/*
 * Edge.hh
 *
 *  Created on: 8/07/2011
 *      Author: jaubry
 */

#ifndef EDGE_HH_
#define EDGE_HH_

#include "Definitions.hh"

namespace strandsim
{

typedef std::pair<const Vec3d&, const Vec3d&> PairOfConstVertices;

template<typename StrandT>
class Edge: public PairOfConstVertices
{
public:
    Edge(const StrandT& strand, const Vec3d& first, const Vec3d& last) :
        m_strand(strand), PairOfConstVertices(first, last)
    {
    }

    Edge(const StrandT strand, const IndexType vtx) :
        m_strand(strand), PairOfConstVertices(strand.getVertex(vtx), strand.getVertex(vtx + 1))
    {
    }

private:
    const StrandT& m_strand;
    // The strand this edge belongs to. Useful to have a back pointer for collision response.
    // We may need a non-const version btw.

};

}

#endif /* EDGE_HH_ */

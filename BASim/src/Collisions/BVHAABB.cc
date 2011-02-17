#include "BVHAABB.hh"

namespace BASim 
{

bool areEdgesEqual( int e0v0, int e0v1, int e1v0, int e1v1 )
{
  return ( (e0v0 == e1v0) && (e0v1 == e1v1) ) || ( (e0v0 == e1v1) && (e0v1 == e1v0) );
}

bool areEdgeEdgeContinuousTimeTheSame( const EdgeEdgeContinuousTimeCollision& col0, const EdgeEdgeContinuousTimeCollision& col1 )
{
  // a == c && b == d
  if( areEdgesEqual( col0.e0_v0, col0.e0_v1, col1.e0_v0, col1.e0_v1 ) && areEdgesEqual( col0.e1_v0, col0.e1_v1, col1.e1_v0, col1.e1_v1 ) ) return true;
  // a == d && c == d
  if( areEdgesEqual( col0.e0_v0, col0.e0_v1, col1.e1_v0, col1.e1_v1 ) && areEdgesEqual( col0.e1_v0, col0.e1_v1, col1.e0_v0, col1.e0_v1 ) ) return true;

  return false;
}

}

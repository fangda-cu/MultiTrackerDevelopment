
#ifndef SAMPLESEEDER_H
#define SAMPLESEEDER_H

#include <vec.h>

namespace ElTopo{
class SurfTrack;
}

class SampleSeeder
{
   
public:
   
   /// Create a set of BCC vertices.
   ///
   static void generate_bcc_points( const ElTopo::Vec3f& domain_low, 
                                    const ElTopo::Vec3f& domain_high, 
                                    float dx, 
                                    std::vector<ElTopo::Vec3f>& xs );

   /// Create a set of points fitting the given surface.
   ///
   static void generate_adaptive_points( const ElTopo::SurfTrack& surface, 
                                         double desired_dx, 
                                         std::vector<ElTopo::Vec3f>& xs );
      
};


#endif


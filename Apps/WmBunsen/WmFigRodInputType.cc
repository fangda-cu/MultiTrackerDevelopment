#include "WmFigRodInputType.hh"
#include <weta/Wnunbsc.hh>

#include <vector>

using namespace std;

WmFigRodInputType::WmFigRodInputType()
{
}

WmFigRodInputType::~WmFigRodInputType()
{
}

void WmFigRodInputType::updateRodDataFromInput( MDataBlock& i_dataBlock )
{
}

void WmFigRodInputType::resampleCurve( size_t i_numVerticesToResample, vector<MVector>& i_curve, 
    vector<MVector>& o_resampledCurve )
{
    //////////////////////////////////////////////////////
    //
    // Build nurbs curves from the input curve. We will
    // use this curve to do a better resampling.
    //
    //////////////////////////////////////////////////////

    Wvec3d positions[ i_curve.size() ];
    for( size_t v = 0; v < i_curve.size(); ++v )
    {
        positions[ v ] = Wvec3d( i_curve[ v ].x, i_curve[ v ].y, i_curve[ v ].z );
    }

    Wnunbsc cubicCurve;
    cubicCurve.centripetalbesselcubicinterp( positions, (int)i_curve.size() );
    
    Wnunbsc sampleCurve;
    cubicCurve.calcpointsperspan( sampleCurve, 3 );
    
    double start = sampleCurve.Knot.K[ sampleCurve.Knot.Deg   - 1 ];
    double end   = sampleCurve.Knot.K[ sampleCurve.Knot.NumCV - 1 ];

    vector< double > paraArray;
    vector< Wvec3d > tangentArray, ptArray;
    
    o_resampledCurve.resize( i_numVerticesToResample );

    for( size_t v = 0; v < i_numVerticesToResample; ++v )
    {
        double t = 0.0;
        float w = float( v ) / float( i_numVerticesToResample - 1 );
        double step = 1.0 / double( i_numVerticesToResample - 1  );
        
        if ( v == 0 )
        {
            t = start;
        }
        else if ( v == i_numVerticesToResample - 1 )
        {
            t = end;
        }
        else
        {
            t = start +  v * step;
        }
            if ( t > end )
        {
            t = end;
        }
        
        Wvec3d pt, tangent;
        sampleCurve.calcpoint( pt, tangent, t );
        paraArray.push_back( t );
        ptArray.push_back( pt );

        tangent = w_normalize( tangent );
        tangentArray.push_back( tangent );
        
        o_resampledCurve[ v ] = MVector( pt[0], pt[1], pt[2] );
    }
}

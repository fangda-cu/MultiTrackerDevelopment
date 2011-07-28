/*
 * wStrandTest.cc
 *
 *  Created on: 18/07/2011
 *      Author: jaubry
 */

#include "wStrandTest.hh"
#include "ElasticStrand.hh"
#include "ElasticStrandStaticStepper.hh"
#include "Forces/ClumpingForce.hh"
#include "Forces/GravitationForce.hh"

#include "../BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "../BASim/src/Physics/ElasticRods/RodUtils.hh"
#include "../BASim/src/Physics/ElasticRods/GroomingTimeStepper.hh"
#include "../BASim/src/Physics/ElasticRods/RodGravity.hh"
#include "../BASim/src/Core/Timer.hh"

using namespace strandsim;

static const int nVertices = 30;
static const int nDOFs = 4 * nVertices - 1;
static const Scalar totalLength = 10.0;
static const Scalar radiusA = 0.5;
static const Scalar radiusB = 1.5;
static const Scalar YoungsModulus = 10000.0;
static const Scalar shearModulus = 1000.0;
static const Scalar density = 1.0;
static const Scalar baseRotation = 0.5;
static const Vec3d gravity( 0.0, 0.0, -981.0 );
static const int nIterations = 1000;

void testStrandSim( const std::vector<Vec3d>& i_vertices )
{
    GravitationForce::setGravity( gravity );

    ElasticStrandParameters params( radiusA, radiusB, YoungsModulus, shearModulus, density,
            baseRotation );
    VecXd dofs( nDOFs );
    for ( int i = 0; i < dofs.size(); i += 4 )
        dofs.segment<3> ( i ) = i_vertices[i / 4];
    ElasticStrand strand( dofs, params );

    strand.addExternalForce( new ClumpingForce );

    ElasticStrandStaticStepper stepper;

    for ( int i = 0; i < nIterations; ++i )
    {
        // std::cout << "\nStrandSim Iteration number " << i << '\n';

        stepper.execute( strand );

        // std::cout << "Press ENTER to continue...\n";
        // std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }
    std::cout << "Vertices: " << strand << '\n';

    std::cout << std::endl;
}

void testBASim( const std::vector<Vec3d>& i_vertices )
{
    using namespace BASim;

    RodOptions rodOptions;
    rodOptions.YoungsModulus = YoungsModulus; /* megapascal */
    rodOptions.ShearModulus = shearModulus; /* megapascal */
    rodOptions.viscosity = 0; /* poise */
    rodOptions.density = density; /* grams per cubic centimeter */
    rodOptions.radiusA = radiusA; /* millimeter */
    rodOptions.radiusB = radiusB; /* millimeter */
    rodOptions.refFrame = BASim::ElasticRod::TimeParallel;
    rodOptions.numVertices = nVertices;

    // Use the rod helper function to build the rod
    ElasticRod* rod = setupRod( rodOptions, i_vertices, i_vertices,
            findNormal<3> ( i_vertices[1] - i_vertices[0] ), Vec3d(), baseRotation );

    // Create a timeStepper to simulate the rod forward in time
    GroomingTimeStepper* stepper = new GroomingTimeStepper( *rod );
    stepper->setDiffEqSolver( BASim::GroomingTimeStepper::STATICS );

    stepper->addExternalForce( new RodGravity( gravity ) );

    // Set the rod's fixed vertices
    RodBoundaryCondition* boundary = stepper->getBoundaryCondition();
    boundary->setDesiredVertexPosition( 0, rod->getVertex( 0 ) );
    boundary->setDesiredVertexPosition( 1, rod->getVertex( 1 ) );
    boundary->setDesiredEdgeAngle( 0, rod->getTheta( 0 ) );

    for ( int i = 0; i < nIterations; ++i )
    {
        // std::cout << "\nBASim iteration number " << i << '\n';

        rod->setIsInRestState( false );
        stepper->execute();

        //  std::cout << "Press ENTER to continue...\n";
        //  std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }
    std::cout << "Vertices: ";
    std::cout << '{';
    for ( int i = 0; i < rod->nv() - 1; i++ )
    {
        const Vec3d& vertex = rod->getVertex( i );
        std::cout << '{' << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << "}, ";
    }
    const Vec3d& vertex = rod->getVertex( rod->nv() - 1 );
    std::cout << '{' << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << '}';
    std::cout << '}';
    std::cout << '\n';

}

int main()
{
    g_log = new TextLog( std::cerr, MsgInfo::kDebug, true );

    static const int NF = Length<BuiltInForcesList>::value;

    std::cout << "Number of built-in forces = " << NF << '\n';

    // Prepare initial rod/strand position
    std::vector<Vec3d> i_vertices;

    // Store arbitrary vertex coordinates
    for ( int i = 0; i < nVertices; i++ )
        i_vertices.push_back(
                Vec3d( i * 1.0 / ( nVertices - 1 ), sin( 2.0 * i * M_PI / ( nVertices ) ), 0.0 ) );

    // Enforce the total length
    Scalar length = 0.0;
    for ( int i = 0; i < nVertices - 1; i++ )
        length += ( i_vertices[i + 1] - i_vertices[i] ).norm();
    for ( int i = 0; i < nVertices; i++ )
        i_vertices[i] *= totalLength / length;

    START_TIMER("StrandSim");
    std::cout << "This is StrandSim\n";
    testStrandSim( i_vertices );
    STOP_TIMER("StrandSim");

    START_TIMER("BASim");
    std::cout << "This is BASim\n";
    testBASim( i_vertices );
    STOP_TIMER("BASim");

    BASim::Timer::report();

    return 0;
}

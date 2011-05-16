/*
 * TestCase6.cc
 *
 *  Created on: 6/04/2011
 *      Author: jaubry
 */

#include "TestCase6.hh"

namespace BASim
{

TestCase6::TestCase6() :
    Problem("Test Case 6", "Reproducing the eponymous Maya scene."), m_tri_meshes(), m_rods(), m_scripting_controllers(),
            m_steppers(), m_pr_stepper(NULL), m_br_stepper(NULL)
{
    addDynamicsProps();
    addRodOptions();
    addRodTimeStepperOptions();

    // Set Defaults for built-in options
    GetStringOpt("integrator") = "implicit";
    // Gravity in CGS
    GetVecOpt("gravity") = Vec3d(0, -981, 0);
    // Timestep
    GetScalarOpt("dt") = 1.0/24.0;
    // Maximum number of implicit solver itations
    GetIntOpt("iterations") = 100;
    // Number of vertices in each rod
    GetIntOpt("nv") = 25;
    // Assume a quasistatic material frame
    GetBoolOpt("quasistatic") = false;
    // Scale the radii of hairs for rendering/contact so we can see them
    GetScalarOpt("radius-scale") = 10.0;

    GetScalarOpt("mass-damping") = 0;
    GetScalarOpt("viscosity") = 0;

    AddOption("hairrootfile", "file with locations (specified as point on sphere) to grow hair from",
            "assets/hairyball/32normals.txt");
    //AddOption("hairrootfile","file with locations (specified as point on sphere) to grow hair from","assets/hairyball/100normals.txt");
    AddOption("curlyhair", "curly if true, straight if false", false);
    AddOption("enablecollisions", "true if collisions enabled, false otherwise", true);
    AddOption("adaptivestepper", "true if adaptive timestepping should be used, failse otherwise", true);
    AddOption("hairlength", "length if the hair", 60.96 / 4.0);
}

TestCase6::~TestCase6()
{
    // TODO Auto-generated destructor stub
}

void TestCase6::Setup()
{
    loadDynamicsProps();

    RodOptions opts;
    getRodOptions(opts);

    std::vector<Vec3d> vertices;
    Vec3d edgy(1.05289, 0, -0.05828);
    vertices.push_back(Vec3d(-8.40179, 0, 7.43399));
    for (int i = 1; i < opts.numVertices; i++)
        vertices.push_back(vertices.back() + edgy);

    ElasticRod* newrod = setupRod(opts, vertices, vertices);
    RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
    boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
    boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
    boundary->setDesiredEdgeAngle(0, 0.0);

    m_rods.push_back(newrod);

    RodTimeStepper* tstep = getRodTimeStepper(*newrod);
    if (!GetBoolOpt("adaptivestepper"))
    {
        m_steppers.push_back(tstep);
    }
    else
    {
        AdaptiveBinaryStepper* adpstep = new AdaptiveBinaryStepper(newrod, tstep);
        m_steppers.push_back(adpstep);
    }

    // Create a sphere in the center of the screen
    ObjParser objparser;
    TriangleMesh* tri_mesh = new TriangleMesh();
    objparser.loadTriangularMesh( "assets/Parallelepiped.obj", *tri_mesh );
    m_tri_meshes.push_back(tri_mesh);





    for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
    for( int i = 0; i < (int) m_tri_meshes.size(); ++i ) m_world->addObject(m_tri_meshes[i]);







    m_br_stepper = new BridsonStepper(m_rods, m_tri_meshes, m_scripting_controllers, m_steppers, GetScalarOpt("dt"));
    //m_br_stepper->setRodLabels(rod_labels);

    m_world->addController(m_br_stepper);

}

void TestCase6::AtEachTimestep()
{
}
void TestCase6::AfterLoad()
{
}
void TestCase6::AfterStep()
{
}

}

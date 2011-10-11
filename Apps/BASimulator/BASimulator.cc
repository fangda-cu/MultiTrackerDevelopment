/**
 * \file BASimulator.cc
 *
 * \author miklos@cs.columbia.edu (based on problem-setup.cc from sar2120@columbia.edu)
 * \date 09/09/2009
 */

#ifdef WETA
#include <weta/Wfigaro/Core/EigenIncludes.hh>
#else
#include <BASim/BASim>
#endif

#ifdef _MSC_VER
#include <direct.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
//#include <sys/time.h>
#include <tclap/CmdLine.h>
#include <typeinfo>
#include <queue>

#include <iostream>
#include <fstream>

#include "BASimulator.hh"

#ifdef WETA
#include <weta/Wfigaro/Render/TriangleMeshRenderer.hh>
#include <weta/Wfigaro/Render/Camera.hh>
#include <weta/Wfigaro/Render/TrackBall.hh>
#include <weta/Wfigaro/Render/Translator.hh>
#include <weta/Wfigaro/Render/Zoomer.hh>
#include <weta/Wfigaro/Render/ViewController.hh>
#include <weta/Wfigaro/Render/RodRenderer.hh>
//#include <weta/Wfigaro/IO/XMLSceneOutputter.hh>
#include <weta/Wfigaro/Core/StatTracker.hh>
#include <weta/Wfigaro/IO/SerializationUtils.hh>
#include <weta/Wfigaro/IO/SolverIO.hh>
#include <weta/Wfigaro/Render/YImage/YImage.hh>
#include <tclap/SwitchArg.h>
#include <tclap/Constraint.h>
#include <tclap/CmdLine.h>
#define BASIM_VERSION "WETA"
#else
#include "BASim/src/Render/TriangleMeshRenderer.hh"
#include "BASim/src/IO/XMLSceneOutputter.hh"
#include "BASim/src/Core/StatTracker.hh"
#include "BASim/src/IO/SerializationUtils.hh"
#endif

// Might not work in windows? :'(
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

std::vector<Problem*> problems;
#include "CreateProblemVector.inl"

Problem* current_problem = NULL;
ViewController controller;

// Renderable objects include RodRenderers and TriangleMeshRenderers.
std::vector<RenderBase*> renderable_objects;
std::vector<RodRenderer*> rod_renderers;
std::vector<TriangleMeshRenderer*> triangle_mesh_renderers;

bool render = true;
bool paused = true;
bool continuous = true;
int window_width = 512;
int window_height = 512;
int max_frames = -1;
Scalar max_time = -1;
bool progress_indicator = false;
bool generate_movie = false;

// Display the simulation time or not
bool g_dsp_sim_tm = true;

// For outputting time-stamped simulation captures,
// for matching current simulation run to output directory.
time_t g_rawtime;
struct tm* g_timeinfo;

// For dumping movies, etc
Scalar fps = 200;
double frame_period;
int current_frame = 0;
int steps_per_frame;
std::string outputdirectory;

int g_problem_idx = -1;

int last_frame_num = -1;

std::string g_resume_file = "";

// Simulation "breakpoints" for debugging purposes. Times at which the simulation should be paused.
std::queue<double>* g_sim_breakpoints;

void cleanup()
{
    Timer::report();
    Timer::deleteAllTimers();

    IntStatTracker::reportIntTrackers();
    IntStatTracker::deleteAllIntTrackers();
    DoubleStatTracker::reportDoubleTrackers();
    DoubleStatTracker::deleteAllDoubleTrackers();

    PairVectorBase::saveToFiles("stats");
    PairVectorBase::clear();

    if (current_problem != NULL)
        current_problem->BaseFinalize();
    for (size_t i = 1; i < problems.size(); ++i)
    {
        assert(problems[i] != NULL);
        delete problems[i];
    }

    for (int i = 0; i < (int) renderable_objects.size(); ++i)
    {
        // RenderBases should only be deleted here
        assert(renderable_objects[i] != NULL);
        if (renderable_objects[i] != NULL)
            delete renderable_objects[i];
        renderable_objects[i] = NULL;
    }

#ifdef HAVE_PETSC
    PetscUtils::finalizePetsc();
#endif // HAVE_PETSC
}

void SetLighting()
{
    // Create a directional white light with a small ambient component
    glEnable( GL_LIGHT0);
    GLfloat white_ambient[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    glLightfv(GL_LIGHT0, GL_AMBIENT, white_ambient);
    GLfloat white_diffuse[] = { 0.55f, 0.55f, 0.55f, 1.0f };
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white_diffuse);
    GLfloat upper_corner[] = { 1.0f, 1.0f, 1.0f, 0.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, upper_corner);

    // Create a much weaker direction light
    glEnable( GL_LIGHT1);
    GLfloat weak_white_diffuse[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    glLightfv(GL_LIGHT1, GL_DIFFUSE, weak_white_diffuse);
    GLfloat negative_z[] = { 0.0f, 0.0f, 1.0f, 0.0f };
    glLightfv(GL_LIGHT1, GL_POSITION, negative_z);

    glShadeModel( GL_FLAT);
}

//void SetMaterial()
//{
//  GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
//  GLfloat mat_shininess[] = {50.0};
//  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
//  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
//  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
//  glEnable(GL_COLOR_MATERIAL);
//}

void InitMenu()
{
    const int centerMenu = glutCreateMenu(menu);
    glutAddMenuEntry("Object (c)", 'c');
    glutAddMenuEntry("World origin (C)", 'C');

    const int showMenu = glutCreateMenu(menu);
    glutAddMenuEntry("Material (1)", '1');
    glutAddMenuEntry("Reference (2)", '2');
    //glutAddMenuEntry("Bishop (3)", '3');
    //glutAddMenuEntry("Show u axis (4)", '4');
    //glutAddMenuEntry("Show v axis (5)", '5');
    //glutAddMenuEntry("Curvature (6)", '6');
    glutAddMenuEntry("Velocity (v)", 'v');
    glutAddMenuEntry("Mode (m)", 'm');
    glutAddMenuEntry("Draw arrows (a)", 'a');
    //glutAddMenuEntry("Force (f)", 'f');
    //glutAddMenuEntry("Camera (p)", 'p');
    glutAddMenuEntry("Scale to radius (r)", 'r');

    glutCreateMenu(menu);
    glutAddMenuEntry("Quit (q)", 'q');
    glutAddMenuEntry("Save screenshot (s)", 's');

    glutAddSubMenu("View center", centerMenu);
    glutAddSubMenu("Display", showMenu);
    //glutAddMenuEntry("Write to file (w)", 'w');
    //glutAddMenuEntry("Print vertex velocities (v)", 'v');

    glutAttachMenu( GLUT_RIGHT_BUTTON);
}

Vec3d calcSimCenter()
{
    Vec3d center = Vec3d::Zero();
    int n = 0;

    for (int i = 0; i < (int) renderable_objects.size(); ++i)
    {
        center += renderable_objects[i]->calculateObjectCenter();
        ++n;
    }

    if (n != 0)
        center /= ((double) n);

    return center;
}

Scalar calcViewRadius(Vec3d& simCenter)
{
    Scalar radius = 0.0;

    for (int i = 0; i < (int) renderable_objects.size(); ++i)
    {
        Vec3d center = renderable_objects[i]->calculateObjectCenter();
        Scalar r = renderable_objects[i]->calculateObjectBoundingRadius(center);
        radius = std::max(radius, r + (center - simCenter).norm());
    }

    if (radius == 0.0)
        radius += 0.1;

    return radius;
}

void centerObject()
{
    Vec3d simCenter = calcSimCenter();
    controller.setCenterMode(ViewController::CENTER_OBJECT);
    controller.setViewCenter(simCenter);

    Scalar radius = calcViewRadius(simCenter);
    controller.setBoundingRadius(radius);
}

void InitCamera()
{
    controller.setViewDirection(Vec3d(0, 0, -2));
    centerObject();
    /*
     if (current_problem->GetBoolOpt("read-camera")) {
     Vec3d& eye = current_problem->GetVecOpt("eye");
     Vec3d& up = current_problem->GetVecOpt("up");
     Vec3d& center = current_problem->GetVecOpt("center");

     SMV::Vec3 e(eye(0), eye(1), eye(2));
     SMV::Vec3 u(up(0), up(1), up(2));
     SMV::Vec3 c(center(0), center(1), center(2));
     SMV::Camera& cam = controller.getCamera();
     cam.setEye(e);
     cam.setUp(u);
     cam.setViewCenter(c);

     } else {
     SMV::Vec3 view_dir(0.0, 0.0, -2.0);
     controller.setViewDirection(view_dir);
     centerObject();
     }*/
}

void setOrthographicProjection()
{
    glMatrixMode( GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    gluOrtho2D(0, window_width, 0, window_height);

    glMatrixMode( GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
}

void renderBitmapString(float x, float y, float z, void *font, std::string s)
{
    glDisable( GL_LIGHTING);
    glColor3f(0.0, 0.0, 0.0);
    glRasterPos3f(x, y, z);
    for (std::string::iterator i = s.begin(); i != s.end(); ++i)
    {
        char c = *i;
        glutBitmapCharacter(font, c);
    }
    glEnable(GL_LIGHTING);
}

void drawSimpleAxis()
{
    glDisable( GL_LIGHTING);

    // Label x axis
    GLvoid *font_style = GLUT_BITMAP_HELVETICA_18;
    glColor3f(1.0, 0.0, 0.0);
    glRasterPos3f(1.0, 0.0, 0.0);
    glutBitmapCharacter(font_style, 'x');

    // Label y axis
    glColor3f(0.0, 1.0, 0.0);
    glRasterPos3f(0.0, 1.0, 0.0);
    glutBitmapCharacter(font_style, 'y');

    // Label z axis
    glColor3f(0.0, 0.0, 1.0);
    glRasterPos3f(0.0, 0.0, 1.0);
    glutBitmapCharacter(font_style, 'z');

    glLineWidth(1.0);

    // Red x axis
    glColor3f(1.0, 0.0, 0.0);
    glBegin( GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
    glEnd();

    // Green y axis
    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
    glEnd();

    // Blue z axis
    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);
    glEnd();

    glEnable(GL_LIGHTING);
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();

    controller.ApplyCamera();

    //controller.setDefault2D();
    //drawSimpleAxis();

    for (int i = 0; i < (int) renderable_objects.size(); ++i)
        renderable_objects[i]->render();

    setOrthographicProjection();
    glColor3f(0.0, 0.0, 0.0);
    //if( g_dsp_sim_tm ) renderBitmapString( 5, window_height-20, 0.0, GLUT_BITMAP_HELVETICA_18, current_problem->ProblemName() );
    if (g_dsp_sim_tm)
        renderBitmapString(5, (float)window_height - 40, 0.0, GLUT_BITMAP_HELVETICA_18, toString(current_problem->getTime()));

    glutSwapBuffers();
    glPopMatrix();
}

void reshape(int w, int h)
{
    window_width = w;
    window_height = h;

    Camera& c = controller.getCamera();
    c.setPerspective(60, 1);
    Scalar radius = controller.getBoundingRadius();
    c.setZClipping(0.01 * radius, 3 * radius);
    c.setViewport(w, h);

    glutPostRedisplay();
}

void serializeScene(const std::string& output_file);

//bool saved30 = false;
//bool saved100 = false;

void idle()
{
    // Some stuff for dumping movies
    if (fps > 0)
    {
        double seconds_per_frame = 1.0 / ((double) fps);
        double steps_per_second = 1.0 / ((double) current_problem->getDt());
        steps_per_frame = (int)std::max(1.0, floor(seconds_per_frame * steps_per_second + 0.5)); // < rounds
        //std::cout << steps_per_frame << std::endl;
    }

    if (max_frames != -1)
    {
        if (current_problem->dynamicsPropsLoaded())
        {
            if (current_problem->getTime() * fps >= max_frames)
            {
                std::cout << "Computed " << max_frames << " frames. Exiting" << std::endl;
                exit(0);
            }
        }
    }

    if (max_time != -1)
    {
        if (current_problem->dynamicsPropsLoaded())
        {
            if (current_problem->getTime() >= max_time)
            {
                std::cout << "Computed up to time " << max_time << ". Exiting" << std::endl;
                exit(0);
            }
        }
    }

    if (g_sim_breakpoints && current_problem->getTime() >= g_sim_breakpoints->front())
    {
        std::cout << "\033[35;1mBASIMULATOR MESSAGE:\033[m Pausing for breakpoint at time " << current_problem->getTime()
                << "." << std::endl;
        paused = true;
        g_sim_breakpoints->pop();
    }

    if (!paused)
    {
        if (!continuous)
            paused = true;
        current_problem->BaseAtEachTimestep();
        //std::cout << "Time: " << current_problem->getTime() << std::endl;
        if (render)
            glutPostRedisplay();
        //if(render&&generate_movie)
        //{
        //  int file_width = 8;
        //  std::stringstream name;
        //  name << std::setfill('0');
        //  name << "screencap/frame_" << std::setw(file_width) << int(current_problem->getTime()/current_problem->getDt()) << ".png";
        //  saveScreen(name.str());
        //  std::cout << "Saving screencap as: " << name.str() << std::endl;
        //}
    }

    if (render && generate_movie)
    {
        int frame = floor(current_problem->getTime() / current_problem->getDt() + 0.5);
        //if( floor(current_problem->getTime()/frame_period) >= current_frame )
    
        if (frame % steps_per_frame == 0 && last_frame_num != frame)
        {
            last_frame_num = frame;
            std::cout << outputdirectory << std::endl;
#ifdef _MSC_VER
            _mkdir(outputdirectory.c_str());
#else
            mkdir(outputdirectory.c_str(), 0755);
#endif

            int file_width = 20;

            //if( generate_movie == 2 || generate_movie == 3 )
            //{
            //std::stringstream framedirname;
            //framedirname << std::setfill('0');
            //framedirname << outputdirectory << "/frame" << std::setw(file_width) << current_frame;
            //
            ////std::cout << framedirname.str() << std::endl;
            //mkdir(framedirname.str().c_str(), 0755);
            //
            //XMLSceneOutputter scenesaver;
            //scenesaver.outputScene( framedirname.str(), "scene.xml", *current_problem );
            //
            //std::cout << "Frame: " << current_frame << "   Time: " << current_problem->getTime() << "   Scenedir: " << framedirname.str() << std::endl;
            //}

            //if( generate_movie == 1 || generate_movie == 3 )
            //{
            glutPostRedisplay();

            std::stringstream name;
            name << std::setfill('0');
            name << outputdirectory << "/frame" << std::setw(file_width) << current_frame << ".png";

            saveScreen(name.str());
            std::cout << "Frame: " << current_frame << "   Time: " << current_problem->getTime() << "   Screencapture: "
                    << name.str() << std::endl;
            //}

            ++current_frame;
        }
    }

    //  int frame = floor( current_problem->getTime()/current_problem->getDt() + 0.5 );
    //
    //  if( frame == 30 && !saved30 )
    //  {
    //    saved30 = true;
    //    std::cout << "Saving binary to testoutput30.bin" << std::endl;
    //    serializeScene("testoutput30.bin");
    //  }
    //
    //  if( frame == 100 && !saved100 )
    //  {
    //    saved100 = true;
    //    std::cout << "Saving binary to testoutput100.bin" << std::endl;
    //    serializeScene("testoutput100.bin");
    //  }

    if (progress_indicator && current_problem->dynamicsPropsLoaded())
        std::cout << "\rtime = " << current_problem->getTime() << std::flush;
}

void scaleMousePos(const int x, const int y, Scalar& xx, Scalar& yy)
{
    int w, h;
    controller.getCamera().getViewport(&w, &h);

    xx = 2 * x / (Scalar) (w - 1) - 1.0;
    yy = 2 * (h - y - 1) / (Scalar) (h - 1) - 1.0;
}

void mouse(int button, int state, int x, int y)
{
    const bool zooming = (button == GLUT_MIDDLE_BUTTON) || ((button == GLUT_LEFT_BUTTON) && (glutGetModifiers()
            & GLUT_ACTIVE_CTRL));
    const bool translating = (button == GLUT_LEFT_BUTTON) && (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
    const bool rotating = (button == GLUT_LEFT_BUTTON) && (glutGetModifiers() == 0);

    Scalar xx, yy;
    scaleMousePos(x, y, xx, yy);
    if (state == GLUT_DOWN)
    {
        if (translating)
            controller.beginTranslationDrag(xx, yy);

        if (zooming)
            controller.beginZoomDrag(xx, yy);

        if (rotating)
            controller.beginRotationDrag(xx, yy);

    }
    else
    {
        controller.endTranslationDrag(xx, yy);
        controller.endRotationDrag(xx, yy);
        controller.endZoomDrag(xx, yy);
    }

    glutPostRedisplay();
}

void motion(int x, int y)
{
    Scalar xx, yy;
    scaleMousePos(x, y, xx, yy);
    controller.updateDrag(xx, yy);
    glutPostRedisplay();
}

void keyboard(unsigned char key, int, int)
{
    menu(key);
}

void menu(int id)
{
    switch (id)
    {
    case 'q':
        exit(0);
        break;

    case ' ':
        paused = !paused;
        break;

    case 'i':
        continuous = !continuous;
        break;

    case '1':
    {
        for (int i = 0; i < (int) rod_renderers.size(); ++i)
        {
            bool& b = rod_renderers[i]->drawMaterial();
            b = !b;
        }
        glutPostRedisplay();
        break;
    }

    case '2':
    {
        for (int i = 0; i < (int) rod_renderers.size(); ++i)
        {
            bool& b = rod_renderers[i]->drawReference();
            b = !b;
        }
        glutPostRedisplay();
        break;
    }

    case 'v':
    {
        for (int i = 0; i < (int) rod_renderers.size(); ++i)
        {
            bool& bv = rod_renderers[i]->drawVelocity();
            bool& br = rod_renderers[i]->drawResponse();

            if (bv)
            {
                bv = false;
                br = true;
            }
            else if (br)
                br = false;
            else
                bv = true;
        }
        glutPostRedisplay();
        break;
    }

    case 'm':
    {
       
       for(unsigned int i = 0; i < renderable_objects.size(); ++i)
          renderable_objects[i]->cycleMode();

       //for (int i = 0; i < (int) rod_renderers.size(); ++i)
       // {
       //     RodRenderer::DrawMode mode = rod_renderers[i]->getMode();
       //     mode = (RodRenderer::DrawMode) ((mode + 1) % RodRenderer::NONE);
       //     rod_renderers[i]->setMode(mode);
       // }
       // for (int i = 0; i < (int) triangle_mesh_renderers.size(); ++i)
       // {
       //     TriangleMeshRenderer::DrawMode mode = triangle_mesh_renderers[i]->getMode();
       //     mode = (TriangleMeshRenderer::DrawMode) ((mode + 1) % TriangleMeshRenderer::NONE);
       //     triangle_mesh_renderers[i]->setMode(mode);
       // }
        glutPostRedisplay();
        break;
    }

    case 'r':
    {
        for (int i = 0; i < (int) rod_renderers.size(); ++i)
        {
            bool& b = rod_renderers[i]->scaleToRadius();
            b = !b;
        }
        glutPostRedisplay();
        break;
    }

    case 'a':
    {
        for (int i = 0; i < (int) rod_renderers.size(); ++i)
        {
            bool& b = rod_renderers[i]->drawArrows();
            b = !b;
        }
        glutPostRedisplay();
        break;
    }

    case 's':
        saveScreen("screen.png");
        break;

    case 'd':
        generate_movie = !generate_movie;
        std::cout << "Dump movie: " << generate_movie << std::endl;
        break;

    case 'c':
        centerObject();
        controller.setCenterMode(ViewController::CENTER_OBJECT);
        glutPostRedisplay();
        break;

    case 'C':
        controller.setCenterMode(ViewController::CENTER_WORLD_ORIGIN);
        glutPostRedisplay();
        break;

    case 'f':
        current_problem->BaseAtEachTimestep();
        glutPostRedisplay();
        break;

    }
    /*
     switch(id) {
     case '2':
     {
     bool& showBishop = rod_renderer->GetOptions().show_bishop;
     showBishop = !showBishop;
     glutPostRedisplay();
     break;
     }
     case '4':
     {
     bool& showU = rod_renderer->GetOptions().show_u;
     showU = !showU;
     glutPostRedisplay();
     break;
     }
     case '5':
     {
     bool& showV = rod_renderer->GetOptions().show_v;
     showV = !showV;
     glutPostRedisplay();
     break;
     }
     case '6':
     {
     bool& showKappa = rod_renderer->GetOptions().show_kappa;
     showKappa = !showKappa;
     glutPostRedisplay();
     break;
     }
     case 'f':
     {
     bool& showForce = rod_renderer->GetOptions().show_force;
     showForce = !showForce;
     glutPostRedisplay();
     break;
     }
     case 'p':
     {
     SMV::Camera& c = controller.getCamera();
     std::cout << "eye: " << c.getEye() << std::endl;
     std::cout << "up: " << c.getUp() << std::endl;
     std::cout << "center: " << c.getViewCenter() << std::endl;
     break;
     }
     }
     case 'v':
     {
     std::vector<Rod*>& rods = current_problem->GetRods();
     for (uint i = 0; i < rods.size(); ++i) {
     Rod& rod = *rods[i];
     for (uint j = 0; j < rod.nv(); ++j) {
     std::cout << rod.velocity(j).norm() << std::endl;
     }
     }
     break;
     }

     case 'w':
     std::vector<Rod*>& rods = current_problem->GetRods();
     for (uint i = 0; i < rods.size(); ++i) {
     std::ostringstream filename;
     filename << "rod_" << i << ".txt";

     std::ofstream out_file;
     out_file.open(filename.str().c_str());
     if (!out_file.is_open()) {
     std::cout << "Failed to write rod to file " << filename << std::endl;

     } else {
     out_file << *rods[i];
     out_file.close();
     }
     }
     break;
     }
     */
}

// Note that serilizeScene is slightly buggy currently.
void serializeScene(const std::string& output_file)
{
    std::ofstream of(output_file.c_str(), std::ios::binary);
    if (!of.is_open())
    {
        std::cerr << "\033[31;1mWARNING IN BASIMULATOR:\033[m Failed to open file for scene serialization: " << output_file
                << ". No output saved." << std::endl;
        return;
    }

    //std::cout << "BASimulator::serializeScene()" << std::endl;

    int BASIMMAGIC = 0xF0F0F0F0;
    of.write((char*) &BASIMMAGIC, sizeof(int));
    of.write((char*) &BASIMMAGIC, sizeof(int));

    // Save the current problem as an integer
    serializeVal(of, g_problem_idx);

    of.write((char*) &BASIMMAGIC, sizeof(int));
    of.write((char*) &BASIMMAGIC, sizeof(int));

    // TODO: save info about dumping frames
    // TODO: save info about stat tracking

    current_problem->serialize(of);

    of.close();
}

// Note that resumeSerializedScene is slightly buggy currently.
void resumeSerializedScene(const std::string& binary_file_name)
{
    std::ifstream ifs(binary_file_name.c_str(), std::ios::binary);
    if (!ifs.is_open())
    {
        std::cerr << "\033[31;1mERROR IN BASIMULATOR:\033[m Failed to resume from serialized scene: " << binary_file_name
                << ". Exiting." << std::endl;
        exit(1);
    }

    int byteeater;
    ifs.read((char*) &byteeater, sizeof(int));
    ifs.read((char*) &byteeater, sizeof(int));

    // Load the current problem
    loadVal(ifs, g_problem_idx);
    assert(g_problem_idx >= 0);
    assert(g_problem_idx < (int) problems.size());
    std::cout << "Resuming problem: " << g_problem_idx << std::endl;

    ifs.read((char*) &byteeater, sizeof(int));
    ifs.read((char*) &byteeater, sizeof(int));

    problems[g_problem_idx]->resumeFromfile(ifs);

    current_problem = problems[g_problem_idx];

    ifs.close();
}

/// Copy the current OpenGL render buffer and save it to a PNG-format file.
bool saveScreen(const std::string& filename)
{
    //int frame = floor( current_problem->getTime()/current_problem->getDt() + 0.5 );
    //int file_width = 8;
    //std::stringstream name;
    //name << std::setfill('0');
    //name << "frame" << std::setw(file_width) << frame << ".bin";

    //std::cout << "Serializing scene to file: " << name.str() << std::endl;

    //std::cout << "Breakpoints: ";
    //std::queue<double> queuecopy = *g_sim_breakpoints;
    //while( !queuecopy.empty() )
    //{
    //  std::cout << queuecopy.front() << " ";
    //  queuecopy.pop();
    //}
    //std::cout << std::endl;

    //serializeScene(name.str());

    //return true;

#ifdef HAVE_PNG
    int w, h;
    controller.getCamera().getViewport(&w, &h);

    YImage image;
    image.resize(w, h);

    glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, image.data());
    image.flip();
    image.setAllAlpha(255);

    return image.save(filename.c_str());
#else
    std::cerr << "Not compiled with PNG support, can't save images." << std::endl;
    return false;
#endif
}

void initializeOpenGL(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
    glutInitWindowPosition(500, 500);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow(argv[0]);
    glEnable( GL_DEPTH_TEST);
    glDepthFunc( GL_LESS);
    //glClearColor(161.0/255.0,161.0/255.0,161.0/255.0,0.0);
    //glClearColor(0.0/255.0,0.0/255.0,0.0/255.0,0.0);
    glClearColor(1.0, 1.0, 1.0, 0.0);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);

    InitMenu();
    InitCamera();

    SetLighting();
    //SetMaterial();
}

void setOptions()
{
    current_problem->AddOption("render", "display output in OpenGL", render);
    current_problem->AddOption("paused", "start simulation paused", paused);
    current_problem->AddOption("continuous", "pause after every step", continuous);
    current_problem->AddOption("window-width", "width of the window", 512);
    current_problem->AddOption("window-height", "height of the window", 512);
    current_problem->AddOption("fps", "frames per second for output movie", fps);
    current_problem->AddOption("max-frames", "number of frames to compute", max_frames);
    current_problem->AddOption("max-time", "maximum (simulation) time", max_time);
    current_problem->AddOption("progress-indicator", "prints out time", progress_indicator);
}

void getOptions()
{
    render = current_problem->GetBoolOpt("render");
    paused = render && current_problem->GetBoolOpt("paused");
    continuous = render && current_problem->GetBoolOpt("continuous");

    window_width = current_problem->GetIntOpt("window-width");
    window_height = current_problem->GetIntOpt("window-height");

    fps = current_problem->GetScalarOpt("fps");
    max_frames = current_problem->GetIntOpt("max-frames");
    max_time = current_problem->GetScalarOpt("max-time");
    progress_indicator = current_problem->GetBoolOpt("progress-indicator");
}

void RunProblem(int argc, char** argv)
{
    if (render)
    {
        initializeOpenGL(argc, argv);
        World& world = current_problem->getWorld();
        World::Objects& objects = world.getObjects();

        //for (World::Objects::const_iterator it = objects.begin(); it != objects.end(); ++it)
        //{
        //    assert((*it) != NULL);

        //    ElasticRod* rod = dynamic_cast<ElasticRod*> (*it);
        //    TriangleMesh* tripobj = dynamic_cast<TriangleMesh*> (*it);

        //    if (rod)
        //    {
        //        rod_renderers.push_back(new RodRenderer(*rod));
        //        renderable_objects.push_back(rod_renderers.back());
        //    }
        //    else if (tripobj)
        //    {
        //        triangle_mesh_renderers.push_back(new TriangleMeshRenderer(*tripobj));
        //        renderable_objects.push_back(triangle_mesh_renderers.back());
        //    }
        //    else
        //    {
        //        std::cout << "Unknown object encountered" << std::endl;
        //    }
        //}

        World::Renderers& rndrs = world.getRenderers();
        for (World::Renderers::const_iterator renderitr = rndrs.begin(); renderitr != rndrs.end(); ++renderitr)
        {
            renderable_objects.push_back(*renderitr);
        }

        centerObject();

        // Extract simulation time breakpoints from the problem
        g_sim_breakpoints = &current_problem->getBreakpoints();
        if (g_sim_breakpoints->size() == 0 || g_sim_breakpoints->back() != std::numeric_limits<double>::infinity())
            g_sim_breakpoints->push(std::numeric_limits<double>::infinity());
         
        glutMainLoop();
    }
    else
    {
        paused = false;
        continuous = true;
        for (;;)
            idle();
    }
}

void PrintProblemTypes()
{
    for (size_t i = 1; i < problems.size(); ++i)
    {
        Problem& prob = *problems[i];
        std::cout << i << ": " << prob.ProblemName() << std::endl << "Description: " << prob.ProblemDescription() << std::endl
                << std::endl;
    }
}

void CreateOptionsFile(int idx, const std::string& filename)
{
    std::ofstream file;

    file.open(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Failed to open file " << filename << std::endl;
        return;
    }

    problems[idx]->PrintOptions(file);

    file.close();

    std::cout << "Generated options file for " << problems[idx]->ProblemName() << ": " << filename << std::endl;
}

template<typename T>
class ProblemConstraint: public TCLAP::Constraint<T>
{
public:

    ProblemConstraint(const T& lower, const T& upper) :
        m_lower(lower), m_upper(upper)
    {
        m_typeDesc = toString(lower) + "-" + toString(upper);
    }

    virtual std::string description() const
    {
        return m_typeDesc;
    }
    virtual std::string shortID() const
    {
        return m_typeDesc;
    }

    virtual bool check(const T& value) const
    {
        return ((value >= m_lower) && (value <= m_upper));
    }

protected:

    T m_lower;
    T m_upper;
    std::string m_typeDesc;
};

int parseCommandLine(int argc, char** argv)
{
    ProblemConstraint<int> allowedProblems(1, (int) (problems.size() - 1));

    try
    {
        std::string version(BASIM_VERSION);
        TCLAP::CmdLine cmd("Send inquiries to basim.support@gmail.com", ' ', version);

        // One and only one of the following arguments is required

        const std::string p1("p");
        const std::string p2("print");
        const std::string p3("print stuff");

        TCLAP::SwitchArg print(p1, p2, p3);
        // TCLAP::SwitchArg print("p", "print", "Print available problems");
        TCLAP::ValueArg<int> run("r", "run", "Run a problem", false, 1, &allowedProblems);
        TCLAP::ValueArg<int> opts("o", "options", "Print a problem's options", false, 1, &allowedProblems);
        TCLAP::ValueArg<int> gen("g", "generate", "Generate a default options file for a problem", false, 1, &allowedProblems);
        TCLAP::ValueArg<std::string> cntu("c", "continue", "Continue a simulation from a saved binary", false, "", "string");

        // Arguments that modify the behavior of r
        TCLAP::ValueArg<std::string> file("f", "file", "Options file for a problem", false, "", "string", cmd);
        TCLAP::ValueArg<std::string> solver("s", "solver", "File describing options for solver", false, "", "string", cmd);

        // Require one and only one of print, run, options, or generate modes
        std::vector<TCLAP::Arg*> xorlist;
        xorlist.push_back(&print);
        xorlist.push_back(&run);
        xorlist.push_back(&opts);
        xorlist.push_back(&gen);
        xorlist.push_back(&cntu);
        cmd.xorAdd(xorlist);

        cmd.parse(argc, argv);

        if (solver.isSet())
        {
            if (readSolverFile(solver.getValue()) != 0)
                return -1;
        }

        if (print.getValue())
        {
            PrintProblemTypes();
            return -1;
        }

        if (opts.isSet())
        {
            int idx = opts.getValue();
            std::cout << "# Options for " << problems[idx]->ProblemName() << std::endl;
            problems[idx]->PrintOptions(std::cout);
            return -1;
        }

        if (gen.isSet())
        {
            if (!file.isSet())
            {
                TCLAP::ArgException e("Must also specify -" + file.getFlag(), "-" + gen.getFlag());
                throw(e);
            }

            CreateOptionsFile(gen.getValue(), file.getValue());
            return -1;
        }

        if (run.isSet())
        {
            int idx = run.getValue();
            current_problem = problems[idx];
            g_problem_idx = idx;
            setOptions();

            if (file.isSet())
            {
                if (current_problem->LoadOptions(file.getValue()) == -1)
                    return -1;
            }
            else
            {
                std::cout << "No options file specified. Using default options." << std::endl;
            }
            if (current_problem->LoadOptions(argc, argv) == -1)
                return -1;
            getOptions();
            return idx;
        }

        if (cntu.isSet())
        {
            g_resume_file = cntu.getValue();
            return 0;
        }

        std::cerr << cmd.getProgramName() << ": missing operand" << std::endl << "Try `" << cmd.getProgramName() << " --help'"
                << " for more information" << std::endl;

    } catch (TCLAP::ArgException& e)
    {
        std::cerr << "ERROR: " << e.argId() << std::endl << "       " << e.error() << std::endl;
    }

    return -1;
}

void printCommandLineSplashScreen()
{
    assert(g_problem_idx >= 0);
    assert(g_problem_idx < (int) problems.size());

    std::cout << "\033[32;1m";
    std::cout << "----------------------------------" << std::endl;
    std::cout << "  ____    _    ____  _" << std::endl;
    std::cout << " | __ )  / \\  / ___|(_)_ __ ___" << std::endl;
    std::cout << " |  _ \\ / _ \\ \\___ \\| | '_ ` _ \\" << std::endl;
    std::cout << " | |_) / ___ \\ ___) | | | | | | |" << std::endl;
    std::cout << " |____/_/   \\_\\____/|_|_| |_| |_|" << std::endl;
    std::cout << "----------------------------------" << std::endl;
    std::cout << "\033[m";
    std::cout << std::endl;

#ifdef DEBUG
    std::cout << " Build mode: DEBUG" << std::endl;
#else
    std::cout << " Build mode: RELEASE" << std::endl;
#endif
    std::cout << std::setfill('0');
    std::cout << " Timestamp: " << (1900 + g_timeinfo->tm_year) << "/" << std::setw(2) << (1 + g_timeinfo->tm_mon) << "/";
    std::cout << (g_timeinfo->tm_mday) << "  " << std::setw(2) << (g_timeinfo->tm_hour) << ":" << std::setw(2)
            << (g_timeinfo->tm_min);
    std::cout << ":" << std::setw(2) << (g_timeinfo->tm_sec) << std::endl;

    std::cout << " Simulation: " << current_problem->ProblemName() << std::endl;

    std::cout << " Linear Solver: " << SolverUtils::instance()->getSolverName() << std::endl;

#ifdef HAVE_OPENMP
    std::cout << " OpenMP: Enabled" << std::endl;
    std::cout << " OpenMP Max Threads: " << omp_get_max_threads() << std::endl;
#else
    std::cout << " OpenMP: Disabled" << std::endl;
#endif

    std::cout << std::endl;
}

std::string generateOutputDirName()
{
    assert(g_timeinfo != NULL);

    std::stringstream datestream;
    datestream.fill('0');
    datestream << std::setw(4) << (1900 + g_timeinfo->tm_year) << "_" << std::setw(2) << (1 + g_timeinfo->tm_mon) << "_";
    datestream << std::setw(2) << (g_timeinfo->tm_mday) << "_" << std::setw(2) << (g_timeinfo->tm_hour) << "_" << std::setw(2)
            << (g_timeinfo->tm_min) << "_";
    datestream << std::setw(2) << (g_timeinfo->tm_sec) << "_" << "simulation_capture";

    return datestream.str();
}

int main(int argc, char** argv)
{
#ifdef HAVE_PETSC
    PetscUtils::initializePetsc(&argc, &argv);
#endif // HAVE_PETSC
    // Generate a directory name with the date and time
    time(&g_rawtime);
    g_timeinfo = localtime(&g_rawtime);
    outputdirectory = generateOutputDirName();
    std::cout << "Output dir: " << outputdirectory << std::endl;
    CreateProblemVector();
    atexit(cleanup);
    if (parseCommandLine(argc, argv) < 0)
        return -1;

    // TODO: Load from serialized output goes here
    if (g_resume_file == std::string(""))
        current_problem->BaseSetup(argc, argv);
    else
        resumeSerializedScene(g_resume_file);

    printCommandLineSplashScreen();

    if(g_log == NULL)
      g_log = new TextLog(std::cerr, MsgInfo::kDebug, true);

    RunProblem(argc, argv);
    return 0;
}

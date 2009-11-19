/**
 * \file BASimulator.cc
 *
 * \author miklos@cs.columbia.edu (based on problem-setup.cc from sar2120@columbia.edu)
 * \date 09/09/2009
 */

#include <BASim/BASim>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <tclap/CmdLine.h>

#include "BASimulator.hh"

std::vector<Problem*> problems;
#include "CreateProblemVector.inl"

Problem* current_problem = NULL;
ViewController controller;
RodRenderer* rod_renderer = NULL;
bool render = true;
bool paused = true;
bool continuous = true;
int window_width = 512;
int window_height = 512;
Scalar fps = 30;
int max_frames = -1;
Scalar max_time = -1;
bool progress_indicator = false;

void cleanup()
{
  if (current_problem != NULL) current_problem->BaseFinalize();
  for (size_t i = 1; i < problems.size(); ++i) {
    delete problems[i];
  }
}

void SetLighting()
{
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2);
  glShadeModel(GL_SMOOTH);

  GLfloat light_position1[] = {-1.0, 1.0, 1.0, 1.0};
  GLfloat light_position2[] = {1.0, 0.5, 0.2, 1.0};
  GLfloat light_position3[] = {-0.1, 1.0, -1.0, 1.0};

  //white lights
  GLfloat bright_light[] = {0.7, 0.7, 0.7, 1};
  GLfloat mid_light[] = {0.5, 0.5, 0.5, 1.0};
  GLfloat low_light[] = {0.2, 0.2, 0.2, 1.0};
  GLfloat no_light[] = {0, 0, 0, 1};

  //GLfloat lmodel_ambient[] = {0.1, 0.1, 0.1, 1.0};

  glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, bright_light);
  glLightfv(GL_LIGHT0, GL_SPECULAR, bright_light);

  glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, low_light);
  glLightfv(GL_LIGHT1, GL_SPECULAR, low_light);

  glLightfv(GL_LIGHT2, GL_POSITION, light_position3);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, low_light);
  glLightfv(GL_LIGHT2, GL_SPECULAR, low_light);

  //glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
}

void SetMaterial()
{
  GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat mat_shininess[] = {50.0};
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
  //glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  //glEnable(GL_COLOR_MATERIAL);
}

void InitMenu()
{
  const int centerMenu = glutCreateMenu(menu);
  glutAddMenuEntry("Object (c)", 'c');
  glutAddMenuEntry("World origin (C)", 'C');

  const int showMenu = glutCreateMenu(menu);
  glutAddMenuEntry("Reference (1)", '1');
  glutAddMenuEntry("Material (2)", '2');
  //glutAddMenuEntry("Bishop (3)", '3');
  //glutAddMenuEntry("Show u axis (4)", '4');
  //glutAddMenuEntry("Show v axis (5)", '5');
  //glutAddMenuEntry("Curvature (6)", '6');
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

  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

Vec3d calcSimCenter()
{
  Vec3d center = Vec3d::Zero();
  size_t i = 0;

  const World::Objects& objects = current_problem->getWorld().getObjects();
  for (; i < objects.size(); ++i) {
    if (dynamic_cast<ElasticRod*>(objects[i]) != NULL)
      center += calculateObjectCenter(dynamic_cast<ElasticRod&>(*objects[i]));
  }
  center /= i;

  return center;
}

Scalar calcViewRadius(Vec3d& simCenter)
{
  Scalar radius = 0.0;
  size_t i = 0;

  const World::Objects& objects = current_problem->getWorld().getObjects();
  for (; i < objects.size(); ++i) {
    if (dynamic_cast<ElasticRod*>(objects[i]) != NULL) {
      ElasticRod& rod = dynamic_cast<ElasticRod&>(*objects[i]);
      Vec3d center = calculateObjectCenter(rod);
      Scalar r = calculateObjectBoundingRadius(rod, center);
      radius = std::max(radius, r + (center - simCenter).norm());
    }
  }

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

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();

  controller.ApplyCamera();

  rod_renderer->render();

  glutSwapBuffers();
  glPopMatrix();
}

void reshape(int w, int h)
{
  Camera& c = controller.getCamera();
  c.setPerspective(60, 1);
  Scalar radius = controller.getBoundingRadius();
  c.setZClipping(0.1 * radius, 3 * radius);
  c.setViewport(w, h);
}

void idle()
{
  if (max_frames != -1) {
    if (current_problem->dynamicsPropsLoaded()) {
      if (current_problem->getTime() * fps >= max_frames) {
        std::cout << "Computed " << max_frames << " frames. Exiting"
                  << std::endl;
        exit(0);
      }
    }
  }

  if (max_time != -1) {
    if (current_problem->dynamicsPropsLoaded()) {
      if (current_problem->getTime() >= max_time) {
        std::cout << "Computed up to time " << max_time << ". Exiting"
                  << std::endl;
        exit(0);
      }
    }
  }

  if (!paused) {
    if (!continuous) paused = true;
    current_problem->BaseAtEachTimestep();
    if (render) glutPostRedisplay();
  }

  if (progress_indicator && current_problem->dynamicsPropsLoaded()) {
    std::cout << "\rtime = " << current_problem->getTime() << std::flush;
  }
}

void scaleMousePos(const int x, const int y, Scalar* xx, Scalar* yy)
{
  int w, h;
  controller.getCamera().getViewport(&w, &h);
  assert(xx && yy);
  *xx = 2 *           x / (Scalar) (w - 1) - 1.0;
  *yy = 2 * (h - y - 1) / (Scalar) (h - 1) - 1.0;
}

void mouse(int button, int state, int x, int y)
{
  const bool zooming     = (button == GLUT_MIDDLE_BUTTON) ||
    ((button == GLUT_LEFT_BUTTON) && (glutGetModifiers() & GLUT_ACTIVE_CTRL));
  const bool translating = (button == GLUT_LEFT_BUTTON)
    && (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
  const bool rotating    = (button == GLUT_LEFT_BUTTON)
    && (glutGetModifiers() == 0);

  Scalar xx, yy;
  scaleMousePos(x, y, &xx, &yy);
  if (state == GLUT_DOWN) {
    if (translating)
      controller.beginTranslationDrag(xx, yy);

    if (zooming)
      controller.beginZoomDrag(xx, yy);

    if (rotating)
      controller.beginRotationDrag(xx, yy);

  } else {
    controller.endTranslationDrag(xx, yy);
    controller.endRotationDrag(xx, yy);
    controller.endZoomDrag(xx, yy);
  }

  glutPostRedisplay();
}

void motion(int x, int y)
{
  Scalar xx, yy;
  scaleMousePos(x, y, &xx, &yy);
  controller.updateDrag(xx, yy);
  glutPostRedisplay();
}

void keyboard(unsigned char key, int, int)
{
  menu(key);
}

void menu(int id)
{
  switch(id) {
  case 'q':
    exit(0);
    break;

  case ' ':
    paused = !paused;
    break;

  case 'i':
    continuous = !continuous;
    break;

  case '1': {
    bool& b = rod_renderer->drawMaterial();
    b = !b;
    glutPostRedisplay();
    break;
  }

  case '2': {
    bool& b = rod_renderer->drawReference();
    b = !b;
    glutPostRedisplay();
    break;
  }

  case 'm': {
    RodRenderer::DrawMode mode = rod_renderer->getMode();
    mode = (RodRenderer::DrawMode) ((mode + 1) % RodRenderer::NONE);
    rod_renderer->setMode(mode);
    glutPostRedisplay();
    break;
  }

  case 'r': {
    bool& b = rod_renderer->scaleToRadius();
    b = !b;
    glutPostRedisplay();
    break;
  }

  case 'a': {
    bool& b = rod_renderer->drawArrows();
    b = !b;
    glutPostRedisplay();
    break;
  }

  case 's':
    saveScreen("screen.png");
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

/// Copy the current OpenGL render buffer and save it to a PNG-format file.
bool saveScreen(const std::string& filename)
{
#ifndef NO_LIBPNG
  int w, h;
  controller.getCamera().getViewport(&w, &h);

  YImage image;
  image.resize(w, h);

  glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, image.data());
  image.flip();

  return image.save(filename.c_str());
#else
  std::cerr << "Not compiled with libPNG support, can't save images.\n"
            << "Recompile with NO_LIBPNG undefined to enable saving images."
            << std::endl;
#endif
}

void initializeOpenGL(int argc, char** argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
  glutInitWindowPosition(500, 500);
  glutInitWindowSize(window_width, window_height);
  glutCreateWindow(argv[0]);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glClearColor(0,0,0,0);

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);

  InitMenu();
  InitCamera();

  SetLighting();
  SetMaterial();
}

void setOptions()
{
  current_problem->AddOption("render", "display output in OpenGL", render);
  current_problem->AddOption("paused", "start simulation paused", paused);
  current_problem->AddOption("continuous", "pause after every step",
                             continuous);
  current_problem->AddOption("window-width", "width of the window", 512);
  current_problem->AddOption("window-height", "height of the window", 512);
  current_problem->AddOption("fps", "frames per second", fps);
  current_problem->AddOption("max-frames", "number of frames to compute",
                             max_frames);
  current_problem->AddOption("max-time", "maximum (simulation) time",
                             max_time);
  current_problem->AddOption("progress-indicator", "prints out time",
                             progress_indicator);
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
  atexit(cleanup);

  current_problem->BaseSetup(argc, argv);

  if (render) {
    initializeOpenGL(argc, argv);
    World& world = current_problem->getWorld();
    World::Objects& objects = world.getObjects();
    World::Objects::iterator it;
    for (it = objects.begin(); it != objects.end(); ++it) {
      ElasticRod* rod = dynamic_cast<ElasticRod*>(*it);
      if (rod == NULL) continue;
      rod_renderer = new RodRenderer(*rod);
      break;
    }

    glutMainLoop();

  } else {
    paused = false;
    continuous = true;
    while(1) idle();
  }
}

void PrintProblemTypes()
{
  for (size_t i = 1; i < problems.size(); ++i) {
    Problem& prob = *problems[i];
    std::cout << i << ": " << prob.ProblemName() << std::endl
              << "Description: " << prob.ProblemDescription() << std::endl
              << std::endl;
  }
}

void CreateOptionsFile(int idx, const std::string& filename)
{
  std::ofstream file;

  file.open(filename.c_str());
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename << std::endl;
    return;
  }

  problems[idx]->PrintOptions(file);

  file.close();

  std::cout << "Generated options file for "
            << problems[idx]->ProblemName() << ": "
            << filename << std::endl;
}

template <typename T>
class ProblemConstraint : public TCLAP::Constraint<T>
{
public:

  ProblemConstraint(const T& lower, const T& upper)
    : m_lower(lower)
    , m_upper(upper)
  {
    m_typeDesc = toString(lower) + "-" + toString(upper);
  }

  virtual std::string description() const { return m_typeDesc; }
  virtual std::string shortID() const { return m_typeDesc; }

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
  ProblemConstraint<int> allowedProblems(1, problems.size() - 1);

  try {
    std::string version(BASIM_VERSION);
    TCLAP::CmdLine cmd("Send inquiries to basim.support@gmail.com", ' ',
                       version);

    TCLAP::SwitchArg
      print("p", "print", "Print available problems", cmd);

    TCLAP::ValueArg<int>
      run("r", "run", "Run a problem", false, 1, &allowedProblems, cmd);

    TCLAP::ValueArg<int>
      opts("o", "options", "Print a problem's options", false, 1,
           &allowedProblems, cmd);

    TCLAP::ValueArg<std::string>
      file("f", "file", "Options file for a problem", false, "",
           "string", cmd);

    TCLAP::ValueArg<int>
      gen("g", "generate","Generate a default options file for a problem",
          false, 1, &allowedProblems, cmd);

    cmd.parse(argc, argv);

    if (print.getValue()) {
      PrintProblemTypes();
      return -1;
    }

    if (opts.isSet()) {
      int idx = opts.getValue();
      std::cout << "Options for " << problems[idx]->ProblemName() << std::endl;
      problems[idx]->PrintOptions(std::cout);
      return -1;
    }

    if (gen.isSet()) {
      if (!file.isSet()) {
        TCLAP::ArgException e("Must also specify -" + file.getFlag(),
                              "-" + gen.getFlag());
        throw(e);
      }

      CreateOptionsFile(gen.getValue(), file.getValue());
      return -1;
    }

    if (run.isSet()) {
      int idx = run.getValue();
      current_problem = problems[idx];
      setOptions();

      if (file.isSet()) {
        if (current_problem->LoadOptions(file.getValue()) == -1) return -1;
      } else {
        std::cout << "No options file specified. Using default options."
                  << std::endl;
      }
      if (current_problem->LoadOptions(argc, argv) == -1) return -1;
      getOptions();
      return idx;
    }

    std::cerr << cmd.getProgramName() << ": missing operand" << std::endl
              << "Try `" << cmd.getProgramName() << " --help'"
              << " for more information" << std::endl;

  } catch (TCLAP::ArgException& e) {
    std::cerr << "ERROR: " << e.argId() << std::endl
              << "       " << e.error() << std::endl;
  }

  return -1;
}

int main(int argc, char** argv)
{
  CreateProblemVector();
  if (parseCommandLine(argc, argv) < 0) return -1;
  RunProblem(argc, argv);
  return 0;
}

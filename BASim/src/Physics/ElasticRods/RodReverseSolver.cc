
#include "RodReverseSolver.hh"
#include "RodStretchingForce.hh"
#include "RodTwistingForce.hh"
#include "RodTwistingForceSym.hh"
#include "RodBendingForce.hh"
#include "RodBendingForceSym.hh"
#include "RodAnisoForce.hh"

// NOTE
// time step should be assigned in the rod instance if viscous forces are included.
//  ===> WRONG. I'm not going to update undeformed starins of viscous froce

using namespace std;

namespace BASim {
  
RodReverseSolver::RodReverseSolver(ElasticRod *rod, RodTimeStepper *stepper) : m_rod(rod), m_stepper(stepper)
{
  size = (rod->nv() - 2) * 4;
  
  x = VecXd(size);
  dx = VecXd(size);
  r = VecXd(size);
  
  m_A = SolverUtils::instance()->createBandMatrix(size, size, 10, 10);
  m_solver = SolverUtils::instance()->createLinearSolver(m_A);
}

RodReverseSolver::~RodReverseSolver()
{
  if( m_A != NULL ) 
  {
    delete m_A;
    m_A = NULL;
  }
    
  if( m_solver != NULL )
  {
    delete m_solver;
    m_solver = NULL;
  }
}

bool RodReverseSolver::execute() 
{
  std::cout << "Reverse Solver / execute\n";
  
  if (m_rod->nv() < 3) return true;
  
  x.setZero();
  dx.setZero();
  r.setZero();
  m_A->setZero();
  
  // Gather unknowns from force classes
  // size = n : (number of vertices - 2) * 4
  // UNKNOWNS : 
  //   undeformed material curvature, kappa bar [n * 2]   (from second vertex ~ second-to-last one)
  //   undeformed twist, m bar [n]   (from second vertex ~ second-to-last one)
  //   undeforemd edge lengths, e bar [n]  (from second edge ~)
  // kappa_0[2], m_0, |e|_0, kappa_1[2], .........
  // Each force class should know index
  
  // set x
  
  std::vector<RodForce*>& forces = m_rod->getForces();
  std::vector<RodForce*>::iterator fIt;
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
    if ((*fIt)->viscous()) continue;
    
//    if (dynamic_cast<RodStretchingForce*> (*fIt) != NULL ) {
      //(*fIt)->updateUndeformedConfiguration(ldata);
    //}
    if (dynamic_cast<RodBendingForceSym*> (*fIt) != NULL ) {
      RodBendingForceSym* f = dynamic_cast<RodBendingForceSym*> (*fIt);
      
      ElasticRod::vertex_iter vit = m_rod->vertices_begin();
      ++vit;
      
      for(int i=1; i<m_rod->nv()-1; i++, ++vit) {
        Vec2d k = f->getKappaBar(*vit);
        x[(i-1)*4 + 0] = k(0);
        x[(i-1)*4 + 1] = k(1);
      }
    }
    if (dynamic_cast<RodTwistingForceSym*> (*fIt) != NULL ) {
      RodTwistingForceSym* f = dynamic_cast<RodTwistingForceSym*> (*fIt);
      
      ElasticRod::vertex_iter vit = m_rod->vertices_begin();
      ++vit;
      
      // undeformed twist is always set to zero
      // safe to delete
      for(int i=1; i<m_rod->nv()-1; i++, ++vit) {
        x[(i-1)*4 + 2] = f->getUndeformedTwist(*vit);
      }
    }
  }
  
  for(int i=1; i<m_rod->ne(); i++) {
    x[(i-1) * 4 + 3] = m_rod->getEdgeLength(i);
  }
  
    
  // Newton iteration
  const int max_iterations = 50;
  const double stol = 0.00000001;
  int curr_it = 0;
  
  while(curr_it < max_iterations)
  {
//    std::cout << curr_it << " \033[31;1m ITERATION :\033[m " << "\n";
    
    //std::cout << "\033[31;1mCURRENT UNDEFORMED STRAIN:\033[m : " << x << "\n";
    
    VecXd f(m_rod->ndof());
    f.setZero();
    // Fill r = force through force classes
    m_stepper->evaluatePDot(f);

//    std::cout << "entire force \n" << f << "\n";
    r.setZero();
    r = - f.segment(7, size); // b as in Ax=b
    
    //std::cout << "\033[31;1mFORCE:\033[m  " << r.norm() << "\n" << r << "\n";
    
    if ( r.norm() < stol ) {
      return true;
    }
        
    // dx = 0
    dx.setZero();
    m_A->setZero();
    
    // Fill m_A = Jacobian w.r.t. x through force classes
    m_rod->computeReverseJacobian(*m_A);
    m_A->finalize(); // need this?
    
    //m_A->print();
        
    // Run linear solver
    // m_solver->solve
    int status = m_solver->solve(dx, r);
    if( status < 0 )
    {
//      successfull_solve = false;
      std::cerr << "\033[31;1mWARNING IN REVERSE SOLVER:\033[m Problem during linear solve detected. " << std::endl;
    }
        
    //std::cout << "dx  " << dx.norm() << "\n" << dx << "\n";
    
    //VecXd test_Adx (size);
    //test_Adx.setZero();
    
    //m_A->multiply(test_Adx, 1.0, dx);
    //std::cout << "A.dx " << test_Adx <<"\n";
    //std::cout << "A.dx - (-F) " << test_Adx - r <<" " << (test_Adx - r).norm() << "\n";

    // Update x
    x = x + dx;
    
    // Update undeformed configurations in force classes 
    // 
    m_rod->updateReverseUndeformedStrain(x);
    
    // check dx.norm()
    if ( dx.norm() < stol) {
        //std::cout << "converged stol" << std::endl;
      return true;
    }
    
    curr_it++;
  }
  
  
  std::cout << "\033[31;1mWARNING IN REVERSE SOLVER:\033[m Max iterations reached\n";
  return false;
}

  
  


} // namespace BASim

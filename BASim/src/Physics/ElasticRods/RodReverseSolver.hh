// REVERSE HAIRDO

#ifndef ROD_REVERSE_SOLVER_HH
#define ROD_REVERSE_SOLVER_HH

#include "ElasticRod.hh"
#include "RodTimeStepper.hh";
#include "../../Math/MatrixBase.hh"
#include "../../Math/LinearSolverBase.hh"
#include "../../Math/SolverUtils.hh"


namespace BASim {

  class RodReverseSolver //: public RodForceT<EdgeStencil>
  {
    public:

      RodReverseSolver(ElasticRod *rod, RodTimeStepper *stepper);
      ~RodReverseSolver();
  
      bool RodReverseSolver::execute();
      
    private:
      ElasticRod* m_rod;
      RodTimeStepper* m_stepper;
      
      int size; // size of unknown vector - 
      
      VecXd x;  // Unknown vector - material properties
      VecXd dx;  // delta x
      VecXd r;  // residual - force
      
      MatrixBase* m_A;
      LinearSolverBase* m_solver;
      
  };

} // namespace BASim

#endif // ROD_REVERSE_SOLVER_HH
 

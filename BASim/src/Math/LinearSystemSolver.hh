#include "LinearSolverBase.hh"

namespace BASim
{

// Linear solver and matrix that will work together
class LinearSystemSolver
{
public:

    LinearSystemSolver(MatrixBase* lhs, LinearSolverBase* solver) :
        m_sys_size(lhs->rows()), m_lhs(lhs), m_solver(solver)
    {
        assert(lhs->rows() == lhs->cols());
        assert(lhs->rows() == m_sys_size);
    }

    ~LinearSystemSolver()
    {
        if (m_lhs != NULL)
        {
            delete m_lhs;
            m_lhs = NULL;
        }
        if (m_solver != NULL)
        {
            delete m_solver;
            m_solver = NULL;
        }
    }

    int m_sys_size;
    MatrixBase* m_lhs;
    LinearSolverBase* m_solver;
};

// Collection of linear system solvers of different sizes
//  NOT THREAD SAFE
//  ASSUMES YOU WANT A BAND MATRIX OF SIZE 10 :)
class LinearSystemSolverCollection
{
public:

    ~LinearSystemSolverCollection()
    {
        std::map<int, LinearSystemSolver*>::iterator it = m_solver_map.begin();
        for (; it != m_solver_map.end(); ++it)
        {
            assert(it->second != NULL);
            delete it->second;
            it->second = NULL;
        }
    }

    LinearSystemSolver* getLinearSystemSolver(int size)
    {
        assert(size > 0);

        // Attempt to locate a solver of the requested size
        std::map<int, LinearSystemSolver*>::iterator it = m_solver_map.find(size);
        // If a solver of the size exists, return it
        if (it != m_solver_map.end())
        {
            return it->second;
        }
        // Otherwise create a new solver
        int band = 10;
        MatrixBase* lhs = SolverUtils::instance()->createBandMatrix(size, size, band, band);
        LinearSolverBase* solver = SolverUtils::instance()->createLinearSolver(lhs);
        m_solver_map.insert(std::pair<int, LinearSystemSolver*>(size, new LinearSystemSolver(lhs, solver)));

        it = m_solver_map.find(size);
        return it->second;
    }

private:
    std::map<int, LinearSystemSolver*> m_solver_map;
};

}

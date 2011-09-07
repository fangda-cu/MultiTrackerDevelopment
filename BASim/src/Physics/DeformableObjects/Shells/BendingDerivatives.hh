#ifndef _BENDING_DERIVATIVES_H_
#define _BENDING_DERIVATIVES_H_

#include <iostream>

#include "BASim/src/Core/Definitions.hh"
namespace BASim {
// given a dimension, creates a dim x dim sized matrix, where *each entry is a
// 3x3 matrix
// i.e. the dimension is the number of free vertices in the optimization

class EnergyHessian
{
public:
    EnergyHessian(void)
    {
        m_Dim = -1;
        m_Hessian = NULL;
    }

    EnergyHessian(int dim)
    {
      m_Hessian = new Mat3d*[dim];

        for (int i = 0;i < dim;i++)
            m_Hessian[i] = new Mat3d[dim];

        m_Dim = dim;

        ClearHessian();
    }

    ~EnergyHessian()
    {
        if (m_Hessian) {
            for (int i = 0;i < m_Dim;i++)
                if (m_Hessian[i]) {
                    delete [] m_Hessian[i];
                    m_Hessian[i] = NULL;
                }
        }

        delete [] m_Hessian;

        m_Hessian = NULL;
        m_Dim = -1;
    }

    void ClearHessian()
    {
        for (int i = 0;i < m_Dim;i++)
            for (int j = 0;j < m_Dim;j++)
                for (int k = 0;k < 3;k++)
                    for (int l = 0;l < 3;l++)
                        m_Hessian[i][j](k, l) = 0;
    }

    void Add(int index1, int index2, Mat3d value)
    {
        if (index1 > m_Dim || index2 > m_Dim) {
            std::cerr << "Error! Hessian matrix too small" << std::endl;
            assert(0);
        }

        m_Hessian[index1][index2] = m_Hessian[index1][index2] + value;
    }

    Mat3d** m_Hessian;
    int m_Dim;
};

}

#endif // _BENDING_DERIVATIVES_H_

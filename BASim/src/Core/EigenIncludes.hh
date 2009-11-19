/**
 * \file EigenIncludes.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/15/2009
 */

#ifndef EIGENINCLUDES_HH
#define EIGENINCLUDES_HH

#define EIGEN_VECTOR_IO Eigen::IOFormat(8, Eigen::Raw, ",", ",", "", "", "{", "}")
#define EIGEN_MATRIX_IO Eigen::IOFormat(8, Eigen::Raw, ",", ",", "{", "}", "{", "}")
#define EIGEN_SPACES_ONLY_IO Eigen::IOFormat(8, Eigen::Raw, " ", " ", "", "", "", "")
#define EIGEN_DEFAULT_IO_FORMAT EIGEN_MATRIX_IO

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#endif // EIGENINCLUDES_HH

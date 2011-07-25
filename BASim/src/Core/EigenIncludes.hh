/**
 * \file EigenIncludes.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/15/2009
 */

#ifndef EIGENINCLUDES_HH
#define EIGENINCLUDES_HH

#define EIGEN_RAW 0
#define EIGEN_VECTOR_IO Eigen::IOFormat(8, EIGEN_RAW, ",", ",", "", "", "{", "}")
#define EIGEN_MATRIX_IO Eigen::IOFormat(8, EIGEN_RAW, ",", ",", "{", "}", "{", "}")
#define EIGEN_SPACES_ONLY_IO Eigen::IOFormat(8, EIGEN_RAW, " ", " ", "", "", "", "")
#define EIGEN_DEFAULT_IO_FORMAT EIGEN_MATRIX_IO

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO 1

#ifdef WETA
#include <weta/Wfigaro/Eigen/LU>
#include <weta/Wfigaro/Eigen/Core>
#include <weta/Wfigaro/Eigen/Geometry>
#include <weta/Wfigaro/Eigen/StdVector>
#include <weta/Wfigaro/Eigen/Dense>
#else
#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/Dense>
#endif

#endif // EIGENINCLUDES_HH

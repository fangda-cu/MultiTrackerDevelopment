This document describes how to set up the necessary external libraries
and environment variables to compile a fairly optimized version of
BASim.

=== Standard libraries needed ===

You will need OpenGL and GLUT, which should be fairly standard. You
may also want to install the PNG library headers. On Ubuntu, you can
install the libpng-dev package from the repositories.


=== Eigen ===

This is a linear algebra library we use for small, dense matrix and
vector computations used in local force and Jacobian computations. It
is a purely header library, so it is easy to install (no compilation
necessary). We have tested the code with version 2.0.12 of Eigen. You
can obtain the library from:

http://eigen.tuxfamily.org/index.php?title=Main_Page

Download the library and decompress it. This will create a folder
called eigen. After this, define the following environment variable:

export Eigen_INC_DIR=/path/to/eigen


=== TCLAP ===

This is also a purely header library that handles parsing command-line
options. If you are using Ubuntu, then you can install it from from
the repositories by typing

sudo apt-get install libtclap-dev

If you are not using Ubuntu, or this doesn't work, then you can
download the source from:

http://tclap.sourceforge.net/

We have tested the code with version 1.2.0 of TCLAP. After downloading
it, decompress it, which creates a directory called tclap-1.2.0. Define
the following environment variable:

export TCLAP_INC_DIR=/path/to/tclap-1.2.0/include


=== LAPACK ===

LAPACK is the Linear Algebra PACKage. It is optional, but if you have
LAPACK, it will speed up the linear solve needed for simulating rods
dramatically. It can be installed from the Ubuntu repositories by
running

sudo apt-get install liblapack-dev

You *may* need to define the environment variable LAPACK_LIB to point
to the library (usually found in /usr/lib in Ubuntu):

export LAPACK_LIB=/usr/lib


=== PETSc ===

*** INSTALLING PETSC IS OPTIONAL. YOU MAY SKIP THIS STEP. ***

This is a linear algebra library we use for large, sparse linear
solves used in the global time stepping. This can be a bit of a pain
to install correctly. You'll need a fortran compiler (like gfortran)
in addition to g++. You can download release 3.0.0-p9 from

http://www.mcs.anl.gov/petsc/petsc-as/

Decompress the archive, which creates a directory called
petsc-3.0.0-p9. Switch to this directory and run the following command
(note that it spans multiple lines)

./configure --doCleanup --with-clanguage=C++ --with-errorchecking=no \
--with-log=no --with-info=no --with-shared --with-dynamic \
--with-debugging=no --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 \
--FOPTFLAGS=-O3 --with-mpi=0 --download-f-blas-lapack=ifneeded

After this is finished, you will be instructed to define the
environment variables PETSC_DIR and PETSC_ARCH. Make sure to set these
first before you continue. Next, you can compile the library by typing

make


=== BASim ===

You will need cmake (http://www.cmake.org/) in order to compile the
library. Go into the BASim directory. First, create a directory called
build, and switch to it. From there, type

ccmake ..

This will bring up a console configuration utility for cmake. Press
'c' to configure the project. If you don't get any errors, check to
make sure all the settings shown are correct.

By default, LAPACK is USE_LAPACK is set to OFF and USE_MKL is set to
ON. If you followed the instructions above and have LAPACK (but don't
have MKL), then switch both of these settings (USE_LAPACK=ON and
USE_MKL=OFF).

To compile with optimizations, you should change the CMAKE_BUILD_TYPE
from Debug to Release.

Once you've checked all the settings, you need to press 'c' again, and
then press 'g' to generate the actual makefiles. This will return you
back to the command line. To compile, type

make all

After the build is complete, you can test the program out by typing

Apps/BASimulator/BASimulator -r 1

Press the space bar to start the simulation. To specify command-line
arguments, you need to put two dashes (--) after the command above and
then write the arguments, such as

Apps/BASimulator/BASimulator -r 1 -- shape-radius 1

To get a list of command-line options that are available for a
specific problem, type

Apps/BASimulator/BASimulator -o 1


=== Documentation ===

If you have Doxygen installed, you can automatically generate html 
documentation after running cmake by entering the build/doc directory,
and executing the command

doxygen Doxyfile

After this command successfully executes, html documentation will be
available in build/doc/html. If you have the graphviz package 
installed, Doxygen will also generate plots of class hierarchies.
To enable this, edit the above Doxyfile and ensure that
HAVE_DOT is set to yes.

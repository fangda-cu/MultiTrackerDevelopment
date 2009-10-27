=== Standard libraries needed ===

You will need OpenGL and GLUT, which should be fairly standard. You
may also want to install the PNG library headers. On Ubuntu, you can
install the libpng-dev package from the repositories.


=== Eigen ===

This is a linear algebra library we use for small, dense matrix and
vector computations using in local force and Jacobian computations. It
is a purely header library, so it is easy to install (no compilation
necessary). You can obtain the library from:

http://eigen.tuxfamily.org/index.php?title=Main_Page

Download version 2.0.6 of the library and decompress it. This will
create a folder called eigen2. After this, define the following
environment variable:

export Eigen_INC_DIR=/path/to/eigen2

=== TCLAP ===

This is a simple library that handles parsing command-line options. If
you are using Ubuntu, then you can simply install it from from the
repositories by typing

sudo apt-get install libtclap-dev

If you are not using Ubuntu, or this doesn't work, then you can
download the source from:

http://tclap.sourceforge.net/

Download version 1.2.0 of the library. This is also a purely header
library, so just decompress it, which creates a directory called
tclap-1.2.0. After this, define the following environment variable:

export TCLAP_INC_DIR=/path/to/tclap-1.2.0/include

=== PETSc ===

This is a linear algebra library we use for large, sparse linear
solves used in the global time stepping. This can be a bit of a pain
to install correctly. You'll need a fortran compiler (like gfortran)
in addition to g++. You can download release 3.0.0-p8 from

http://www.mcs.anl.gov/petsc/petsc-as/


Decompress the archive, which creates a directory called
petsc-3.0.0-p8. Switch to this directory and run the following command
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
make sure all the settings shown are correct. You can change the
CMAKE_BUILD_TYPE from Debug to Release if you want. Once you've
checked these settings, you need to press 'c' again, and then press
'g' to generate the actual makefiles. This will return you back to the
command line. To compile, type

make

After the build is complete, you can test the program out by typing

Apps/BASimulator/BASimulator -r 1

Press the space bar to start the simulation.

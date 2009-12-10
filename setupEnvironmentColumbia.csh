export MAYA_LOCATION=/data/shared/autodesk/maya2008-x64/
export PATH=/data/shared/autodesk/maya2008-x64/bin/:$PATH
export Eigen_INC_DIR=/data/shared/libraries/eigen2
export TCLAP_INC_DIR=/data/shared/libraries/tclap-1.2.0/include
export PETSC_ARCH=linux-gnu-cxx-opt
export PETSC_DIR=/data/shared/libraries/petsc-3.0.0-p8

export MAYA_PLUG_IN_PATH=`pwd`/Apps/WmBunsen/:${MAYA_PLUG_IN_PATH}
export MAYA_SCRIPT_PATH=`pwd`/../Apps/WmBunsen/:${MAYA_SCRIPT_PATH}
export XBMLANGPATH=`pwd`/../Apps/WmBunsen/icons/%B:${XBMLANGPATH}

#export MKLROOT=/data/shared/libraries/intel/cmkl/10.2.2.025/
export MKLDIR=/data/shared/libraries/intel/cmkl/10.2.2.025/
export MKLLIB=/data/shared/libraries/intel/cmkl/10.2.2.025/lib/em64t/
export LD_LIBRARY_PATH="${MKLROOT}/lib/64:${MKLROOT}/lib/32":${LD_LIBRARY_PATH}
                                                                      

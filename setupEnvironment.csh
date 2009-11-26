setenv Eigen_INC_DIR /local1/apps/eigen2
setenv TCLAP_INC_DIR /local1/apps/tclap-1.2.0/include
setenv PETSC_ARCH linux-gnu-cxx-opt
setenv PETSC_DIR /local1/apps/petsc-3.0.0-p8

setenv MAYA_PLUG_IN_PATH `pwd`/Apps/WmBunsen/:${MAYA_PLUG_IN_PATH}
setenv MAYA_SCRIPT_PATH `pwd`/../Apps/WmBunsen/:${MAYA_SCRIPT_PATH}
setenv XBMLANGPATH `pwd`/../Apps/WmBunsen/icons/%B:${XBMLANGPATH}

# weta specific...
#need intel-11.1.056_64
setenv PATH /vol/apps/intel_64/compiler11.1.056/bin/intel64:${PATH}
setenv LM_LICENSE_FILE /vol/apps/intel_64/compiler11.1.056/licenses:28518@linux-license.wetafx.co.nz:${LM_LICENSE_FILE}
setenv LD_LIBRARY_PATH /vol/apps_master/apps.Linux64/intel_64/compiler11.1.056/lib/intel64/:${LD_LIBRARY_PATH}

#need mkl-10.2.1.017
setenv MKLROOT /local1/apps/mkl-10.2.1.017/
#setenv INCLUDE "${MKLROOT}/include"
setenv LD_LIBRARY_PATH "${MKLROOT}/lib/64:${MKLROOT}/lib/32":${LD_LIBRARY_PATH}


#need cuda-2.3_64
#setenv OCU_LOCATION /vol/weta/src/amcadams/opencurrent

#Need to fix the license path as the VPN connection messes up sometimes. Get Jason in systems to
# change it on the laptop next time I'm up at Manuka...
#setenv LM_LICENSE_FILE 28518@linux-license.wetafx.co.nz:$LM_LICENSE_FILE



                                                                                      

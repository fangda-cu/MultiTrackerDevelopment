
# See /vol/bob/make/bob_head.mk for detailed comments.

# stop Bob picking up old stuff from /vol/weta/include as there are antique versions of tclap
ignoreold=1

##############################################################################
include /vol/bob/make/head.mk
##############################################################################
IPATHS 	+= -I$(MKLROOT)/include 


###{ Name of executable to be built.
###}
PRG1_CF		:= -DWETA -DHAVE_PNG -pthread -DHAVE_LAPACK -DWETA #-openmp -DHAVE_OPENMP

#Nasty but we need it to build BASimulator which uses a file in the BASim library which references a file back in BASim...
PRG1_CF += -I../../../ -I. #-I/vol/bob/src/code/Figaro/work/acoull-breannanCollisions/
PRG1_CF += -Itclap-1.2.0/include/

PRG1		:= basimulator


PRG1_H		:= \
BASimulator.hh \
Problems/ProblemBase.hh \
Problems/BentTwisting.hh \
Problems/CollisionTest.hh \
Problems/CollisionTestTwo.hh \
Problems/CollisionTestWeta.hh \
Problems/RegressionTests.hh \
Problems/HairyBall.hh \
Problems/HairyBallWeta.hh \
Problems/SHO.hh \
Problems/ViscosityTests.hh \
Problems/TestCase6.hh \
tclap-1.2.0/tclap/include/SwitchArg.h

#PRG1_HDIR	:=./Problems

###}
###{ Names of source files for this executable.
PRG1_S		:= \
		BASimulator.cc \
		./Problems/BentTwisting.cc \
		./Problems/CollisionTest.cc \
		./Problems/CollisionTestTwo.cc \
		./Problems/CollisionTestWeta.cc \
		./Problems/ProblemBase.cc \
		./Problems/RegressionTests.cc \
		./Problems/HairyBall.cc \
		./Problems/HairyBallWeta.cc \
		./Problems/SHO.cc \
		./Problems/ViscosityTests.cc  \
		./Problems/TestCase6.cc

###}
###{ Linker specifications for libraries to use when
###  linking this executable.
PRG1_L		:= \
	-lWfigaro_dbg \
    -lglut -lGLU \
    -lpng

# add in the MKL libraries

ifeq ($(LINUX_COMPILER) ,I)
PRG1_L        += -Wl,--start-group \
                  $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a \
                  $(MKLROOT)/lib/em64t/libmkl_intel_thread.a \
                  $(MKLROOT)/lib/em64t/libmkl_core.a \
                 -Wl,--end-group -openmp -lpthread
endif

ifeq ($(LINUX_COMPILER) ,GCC)
PRG1_L        += $(MKLROOT)/lib/em64t/libmkl_solver_ilp64_sequential.a \
				-Wl,--start-group $(MKLROOT)/lib/em64t/libmkl_intel_ilp64.a \
				$(MKLROOT)/lib/em64t/libmkl_sequential.a \
				$(MKLROOT)/lib/em64t/libmkl_core.a \
				-Wl,--end-group -lpthread
endif

#
##############################################################################
include /vol/bob/make/tail.mk
##############################################################################

# Stuff for Emacs:
# Local Variables:
# mode:makefile-gmake
# indent-tabs-mode:t
# tab-width:8
# End:

# Stuff for Vim [by default, this needs to be within the last 5 lines]:
# ex:set filetype=make tabstop=8 softtabstop=8 shiftwidth=8 noexpandtab:

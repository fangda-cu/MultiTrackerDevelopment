##############################################################################
include /vol/bob/make/head.mk
##############################################################################

BUILD := "$(shell perl -e '$$pwd=$$ENV{PWD}; $$time=localtime(); print "\\\"$$pwd $$time\\\"";')"

WITHKITS +=  ext/boost/tags/boost-1_41_0-1
IPATHS 	+= -I$(MKLROOT)/include

######

PRG1       := wStrandTest
PRG1_S     := wStrandTest.cc
PRG1_H     := wStrandTest.hh

PRG1_CF    := -DBUILD=${BUILD} -DWETA

PRG1_L     := -lwStrandSim_dbg -Wl,--start-group \
                  $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a \
                  $(MKLROOT)/lib/em64t/libmkl_intel_thread.a \
                  $(MKLROOT)/lib/em64t/libmkl_core.a \
                 -Wl,--end-group \
                 -openmp \
                 -lWfigaro # testing and comparison
                 
######

LIB1       := wStrandSim

LIB1_CF    += -DBUILD=${BUILD} -DWETA -UNDEBUG
LIB1_CF    += -w1

LIB1_HDIR := .

LIB1_H += BandMatrix.hh
LIB1_H += Definitions.hh
LIB1_H += ElasticStrand.hh
LIB1_H += ElasticStrandParameters.hh
LIB1_H += ElasticStrandStaticStepper.hh
LIB1_H += ElasticStrandUtils.hh
LIB1_H += Forces/BendingForce.hh
LIB1_H += Forces/ForceBase.hh
LIB1_H += Forces/GravitationForce.hh
LIB1_H += Forces/StretchingForce.hh
LIB1_H += Forces/TwistingForce.hh
LIB1_H += LinearSolver.hh
LIB1_H += StepperBase.hh
LIB1_H += StrandBase.hh
LIB1_H += StrandGeometry.hh
LIB1_H += StrandStepManager.hh

LIB1_SDIR := .

LIB1_S += ElasticStrand.cc
LIB1_S += ElasticStrandStaticStepper.cc
LIB1_S += Forces/BendingForce.cc
LIB1_S += Forces/ForceBase.cc
LIB1_S += Forces/GravitationForce.cc
LIB1_S += Forces/StretchingForce.cc
LIB1_S += Forces/TwistingForce.cc
LIB1_S += LinearSolver.cc
LIB1_S += StrandBase.cc
LIB1_S += StrandGeometry.cc
LIB1_S += StrandStepManager.cc



##############################################################################
include /vol/bob/make/tail.mk
##############################################################################

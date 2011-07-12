##############################################################################
include /vol/bob/make/head.mk
##############################################################################

BUILD := "$(shell perl -e '$$pwd=$$ENV{PWD}; $$time=localtime(); print "\\\"$$pwd $$time\\\"";')"

WITHKITS +=  ext/boost/tags/boost-1_41_0-1

LIB1       := wStrandSim

LIB1_CF    += -DBUILD=${BUILD} -openmp -DWETA
LIB1_CF    += -w1

LIB1_H += ElasticStrandUtils.hh
LIB1_H += ElasticStrandParameters.hh
LIB1_H += StrandStepManager.hh
LIB1_H += BandMatrix.hh
LIB1_H += StepperBase.hh
LIB1_H += ElasticStrandStaticStepper.hh
LIB1_H += StrandBase.hh
LIB1_H += ElasticStrand.hh
LIB1_H += Definitions.hh
LIB1_H += Edge.hh
LIB1_H += Force.hh
LIB1_H += StretchingForce.hh
LIB1_H += BendingForce.hh
LIB1_H += TwistingForce.hh
LIB1_H += ForceAccumulator.hh

LIB1_S += StrandStepManager.cc
LIB1_S += StrandBase.cc
LIB1_S += ElasticStrandStaticStepper.cc
LIB1_S += ElasticStrand.cc
LIB1_S += Force.cc
LIB1_S += StretchingForce.cc
LIB1_S += BendingForce.cc
LIB1_S += TwistingForce.cc
LIB1_S += ForceAccumulator.cc




##############################################################################
include /vol/bob/make/tail.mk
##############################################################################

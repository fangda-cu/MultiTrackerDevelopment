USEMAYA         := 1

BUILD := "$(shell perl -e '$$pwd = $$ENV{PWD}; $$time=localtime(); print "\\\"$$pwd $$time\\\"";')"

##############################################################################
include /vol/bob/make/head.mk
##############################################################################

WITHKITS += ext/boost/tags/boost-1_41_0-1

VERSION         := $(KVERMAJ)
PLG_INSD        :=#
SCR_INSD        :=#
GLSL_INSD       :=#
ifdef HAVEMAYA
ifndef INSTALL_AS_DEFAULT
PLG_INSD        := $(wmayaplg)/wmFigaro-$(VERSION)
SCR_INSD        := $(wmayascr)/wmFigaro-$(VERSION)
GLSL_INSD       := $(wetaglsl)/wmFigaro-$(VERSION)
endif
endif

IPATHS 	+= -I$(MKLROOT)/include -I. -I./utils/

PLG1            := wmFigaro.so
#PLG1_CF	        += -Werror -DBUILD=${BUILD}
#PLG1_CF	        += -DBUILD=${BUILD} -openmp -pthread -DHAVE_MKL -DHAVE_LAPACK -DWETA
PLG1_CF	        += -DBUILD=${BUILD} -pthread -DHAVE_LAPACK -DWETA -openmp -DHAVE_OPENMP
PLG1_CF	        += -no-ipo

#PhysBAM libraries
#PHYSBAM_LOCATION := /vol/weta/src/amcadams/PhysBAM-svn
#PHYSBAM_LOCATION := /vol/bob/src/code/PhysBAM/work/girving-trunk
#PHYSBAM_LIBS 	 := -lPhysBAM_Advection_Equations -lPhysBAM_Arrays -lPhysBAM_Articulated_Rigid_Bodies -lPhysBAM_Boundaries -lPhysBAM_Collisions_And_Interactions \
#		    -lPhysBAM_Conservation_Law_Solvers -lPhysBAM_Constitutive_Models -lPhysBAM_OpenGL -lPhysBAM_Data_Structures -lPhysBAM_Optimization \
#		    -lPhysBAM_Deformable_Objects -lPhysBAM_Ordinary_Differential_Equations -lPhysBAM_Euler_Equations -lPhysBAM_Parallel_Computation \
#		    -lPhysBAM_Forces_And_Torques -lPhysBAM_Particles -lPhysBAM_Fourier_Transforms -lPhysBAM_Poisson_Equations \
#		    -lPhysBAM_Fracture -lPhysBAM_Random_Numbers -lPhysBAM_Geometry -lPhysBAM_Read_Write -lPhysBAM_Grid_Based_Fields \
#		    -lPhysBAM_Grids -lPhysBAM_Rigid_Bodies -lPhysBAM_Heat_Flows -lPhysBAM_Solids_And_Fluids -lPhysBAM_Incompressible_Flows -lPhysBAM_Solids_Evolution \
#		    -lPhysBAM_Arrays -lPhysBAM_Interpolation -lPhysBAM_Standard_Tests -lPhysBAM_Lagrange_Equations -lPhysBAM_Turbulence \
#		    -lPhysBAM_Level_Sets -lPhysBAM_Utilities -lPhysBAM_Matrices_And_Vectors -lPhysBAM_Meshing -lPhysBAM_Motion \
#		    -lboost_iostreams -lboost_filesystem -lpng -ljpeg -lglut
#PHYSBAM_RPATH := -Wl,-rpath,$(PHYSBAM_LOCATION)/build/release/Public_Library -Wl,-rpath,$(PHYSBAM_LOCATION)/External_Libraries/boost/lib 
#PHYSBAM_LIB_PATH := -L$(PHYSBAM_LOCATION)/build/release/Public_Library \
 #                   -L/vol/bob/built/linux64/L40/ext/boost/tags/boost-1_41_0-1/lib

#PHYSBAM_INCLUDE  := -I$(PHYSBAM_LOCATION)/Public_Library

PLG1_CF += -Wno-overloaded-virtual -Wno-unused -pipe -fPIC -g -funroll-loops # $(PHYSBAM_INCLUDE)


#ifeq "$(ABI)" "LI40"

#PLG1_L += -L$(MKLROOT)/lib/em64t/ $(MKLROOT)/lib/em64t/libmkl_solver_lp64.a -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread
                   
#endif
#ifeq "$(ABI)" "L40"
# L40 doesn't support -openmp so this won't build and you're out of luck...  
#endif 

PLG1_MAYA       := 1
PLG1_INSD       := $(PLG_INSD)

# Turn of IPO as it takes ages and may not do much at all.
#ABIOPTS_o   := $(filter-out -ipo-separate, $(ABIOPTS_o))

PLG1_I          := \
Beaker.hh \
RodData.hh \
SplineAttrEval.hh \
WmFigaroCmd.hh \
WmBunsenHelpers.hh \
WmBunsenNode.hh \
WmBunsenCollisionMeshNode.hh \
WmFigRodNode.hh \
WmFigConnectionNode.hh \
WmFigSelectionContext.hh \
WmFigSelectionContextCommand.hh \
WmFigSelectionToolCommand.hh \
WmFigControllerContext.hh \
WmFigControllerContextCommand.hh \
WmFigControllerToolCommand.hh \
WmFigRodInputType.hh \
WmFigRodNurbsInput.hh \
WmFigRodBarbInput.hh \
WmFigRodFileIO.hh \
WmFigRodFileInput.hh \
WmFigRodGroup.hh \
WmFigaroRodShape/WmFigaroRodShape.hh \
WmFigaroRodShape/WmFigaroRodShapeIterator.hh \
WmFigaroRodShape/WmFigaroRodShapeUI.hh \
WmFigSelectionDisplayNode.hh \
WmFigUtils.hh \
WmFigRodComponentList.hh  \
WmFigRodParticleInput.hh \
WmFigMeshController.hh \
Sweeney/WmSweeneyNode.hh \
Sweeney/WmSweeneyRodManager.hh \
Sweeney/WmSweeneySubsetNode.hh \
Sweeney/WmSweeneyVolumetricNode.hh \
Sweeney/WmSweeneyCmd.hh \
Sweeney/WmSweeneyUtils.hh \
Sweeney/Tools/AddRod/WmSwAddRodContext.hh \
Sweeney/Tools/AddRod/WmSwAddRodContextCommand.h \
Sweeney/Tools/AddRod/WmSwAddRodToolCommand.hh \
#VolumetricCollisions/VolumetricCollisions.hh \
#VolumetricCollisions/VolumetricCollisionsCPU.hh \
#WmBunsenCollisionMeshNode.hh \
#constraints/WmFigConstraintNode.hh \
#WmFigCombContext.hh \
#WmFigCombContextCommand.hh \
#WmFigCombToolCommand.hh \

PLG1_S          := \
Beaker.cc \
RodData.cc \
SplineAttrEval.cc \
WmFigaroCmd.cc \
WmBunsenNode.cc \
WmBunsenPluginMain.cc \
WmBunsenCollisionMeshNode.cc \
WmFigRodNode.cc \
WmFigConnectionNode.cc \
WmFigSelectionContext.cc \
WmFigSelectionContextCommand.cc \
WmFigSelectionToolCommand.cc \
WmFigControllerContext.cc \
WmFigControllerContextCommand.cc \
WmFigControllerToolCommand.cc \
WmFigRodInputType.cc \
WmFigRodNurbsInput.cc \
WmFigRodBarbInput.cc \
WmFigRodFileIO.cc \
WmFigRodFileInput.cc \
WmFigRodGroup.cc \
WmFigaroRodShape/WmFigaroRodShape.cc \
WmFigaroRodShape/WmFigaroRodShapeIterator.cc \
WmFigaroRodShape/WmFigaroRodShapeUI.cc \
WmFigSelectionDisplayNode.cc \
utils/WmFigUtils.cc \
WmFigRodComponentList.cc \
WmFigRodParticleInput.cc \
WmFigMeshController.cc \
Sweeney/WmSweeneyNode.cc \
Sweeney/WmSweeneySubsetNode.cc \
Sweeney/WmSweeneyVolumetricNode.cc \
Sweeney/WmSweeneyRodManager.cc \
Sweeney/WmSweeneyCmd.cc \
Sweeney/WmSweeneyUtils.cc \
Sweeney/Tools/AddRod/WmSwAddRodContext.cc \
Sweeney/Tools/AddRod/WmSwAddRodContextCommand.cc \
Sweeney/Tools/AddRod/WmSwAddRodToolCommand.cc \
#WmBunsenCollisionMeshNode.cc \
#VolumetricCollisions/VolumetricCollisionsCPU.cc \
#constraints/WmFigConstraintNode.cc \
#WmFigCombContext.cc \
#WmFigCombContextCommand.cc \
#WmFigCombToolCommand.cc \

PLG1_L          += \
	-lWfigaro_dbg \
	-lWexcept_dbg \
	-lWnurbs_dbg \
	-lWmath_dbg \
	-lOpenMaya \
	-lOpenMayaAnim \
	-lOpenMayaUI \
    -lglut

PLG1_L        += -Wl,--start-group \
                  $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a \
                  $(MKLROOT)/lib/em64t/libmkl_intel_thread.a \
                  $(MKLROOT)/lib/em64t/libmkl_core.a \
                 -Wl,--end-group -openmp -lpthread

# PhysBAM libraries
#PLG1_L += $(PHYSBAM_RPATH) $(PHYSBAM_LIB_PATH) $(PHYSBAM_LIBS)


SCR1            := AEwmFigaroNodeTemplate.mel
SCR1_MAYA       := 1
SCR1_INSD       := $(SCR_INSD)

SCR2            := AEwmFigRodNodeTemplate.mel
SCR2_MAYA       := 1
SCR2_INSD       := $(SCR_INSD)

SCR3            := AEwmFigCollisionNodeTemplate.mel
SCR3_MAYA       := 1
SCR3_INSD       := $(SCR_INSD)

SCR4            := wmFigSelectionContextProperties.mel
SCR4_MAYA       := 1
SCR4_INSD       := $(SCR_INSD)

SCR5            := wmFigSelectionContextValues.mel
SCR5_MAYA       := 1
SCR5_INSD       := $(SCR_INSD)

SCR5            := wmFigControllerContextProperties.mel
SCR5_MAYA       := 1
SCR5_INSD       := $(SCR_INSD)

SCR6            := wmFigControllerContextValues.mel
SCR6_MAYA       := 1
SCR6_INSD       := $(SCR_INSD)

SCR7            := WmFigaro.mel
SCR7_MAYA       := 1
SCR7_INSD       := $(SCR_INSD)

SCR8            := wmFigCombContextProperties.mel
SCR8_MAYA       := 1
SCR8_INSD       := $(SCR_INSD)

SCR9            := wmFigCombContextValues.mel
SCR9_MAYA       := 1
SCR9_INSD       := $(SCR_INSD)

SCR10            := WmFigaroRodShape/AEwmFigaroRodShapeTemplate.mel
SCR10_MAYA       := 1
SCR10_INSD       := $(SCR_INSD)

SCR11            := constraints/AEwmFigConstraintNodeTemplate.mel
SCR11_MAYA       := 1
SCR11_INSD       := $(SCR_INSD)

SCR12            := Sweeney/wmSweeney.mel
SCR12_MAYA       := 1
SCR12_INSD       := $(SCR_INSD)

SCR13            := Sweeney/AEwmSweeneyNodeTemplate.mel
SCR13_MAYA       := 1
SCR13_INSD       := $(SCR_INSD)

SCR14            := Sweeney/AEwmSweeneySubsetNodeTemplate.mel
SCR14_MAYA       := 1
SCR14_INSD       := $(SCR_INSD)

##############################################################################
include /vol/bob/make/tail.mk
##############################################################################

# Stuff for Emacs:
# Local Variables:
# mode:makefile
# indent-tabs-mode:t
# tab-width:8
# End:

# Stuff for Vim [by default, this needs to be within the last 5 lines]:
# ex:set filetype=make tabstop=8 softtabstop=8 shiftwidth=8 noexpandtab:



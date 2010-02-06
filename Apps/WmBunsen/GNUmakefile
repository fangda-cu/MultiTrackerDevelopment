WITHKITS += Wmaya-2B1370

USEMAYA         := 1

BUILD := "$(shell perl -e '$$pwd = $$ENV{PWD}; $$time=localtime(); print "\\\"$$pwd $$time\\\"";')"

##############################################################################
include /vol/bob/make/head.mk
##############################################################################

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

IPATHS 	+= -I$(MKLROOT)/include 

PLG1            := wmFigaro.so
#PLG1_CF	        += -Werror -DBUILD=${BUILD}
PLG1_CF	        += -DBUILD=${BUILD}

ifeq "$(ABI)" "LI40"
PLG1_LCF        += $MKLROOT/libmkl_solver_lp64.a \
                   -Wl,--start-group \
                   $MKLROOT/libmkl_intel_lp64.a \
                   $MKROOT/libmkl_intel_thread.a \
                   $MKLROOT/libmkl_core.a \
                   -Wl,--end-group -openmp -lpthread
endif
ifeq "$(ABI)" "L40"
PLG1_LCF        +=  $MKLROOT/libmkl_solver_lp64.a \
                    -Wl,--start-group \
                    $MKLPATH/libmkl_intel_lp64.a \
                    $MKLPATH/libmkl_intel_thread.a \
                    $MKLPATH/libmkl_core.a \
                    -Wl,--end-group -liomp5 -lpthread
endif

PLG1_MAYA       := 1
PLG1_INSD       := $(PLG_INSD)

# Turn of IPO as it takes ages and may not do much at all.
ABIOPTS_o   := $(filter-out -ipo-separate, $(ABIOPTS_o))

PLG1_I          := \
Beaker.hh \
RodData.hh \
SplineAttrEval.hh \
WmBunsenCmd.hh \
WmBunsenCollisionMeshNode.hh \
WmBunsenHelpers.hh \
WmBunsenNode.hh \
WmBunsenRodNode.hh

PLG1_S          := \
Beaker.cc \
RodData.cc \
SplineAttrEval.cc \
WmBunsenCmd.cc \
WmBunsenCollisionMeshNode.cc \
WmBunsenNode.cc \
WmBunsenPluginMain.cc \
WmBunsenRodNode.cc

        
PLG1_L          := \
	-lWfigaro_dbg \
	-lWexcept_dbg \
	-lOpenMaya \
	-lOpenMayaAnim \
	-lOpenMayaUI \
    -lglut
	
SCR1            := AEwmBunsenNodeTemplate.mel
SCR1_MAYA       := 1
SCR1_INSD       := $(SCR_INSD)


SCR2            := AEwmBunsenRodNodeTemplate.mel
SCR2_MAYA       := 1
SCR2_INSD       := $(SCR_INSD)


SCR3            := WmBunsen.mel
SCR3_MAYA       := 1
SCR3_INSD       := $(SCR_INSD)


# i somehow broke this .. butg do we need it? -SHANE -
#ifdef HAVEMAYA
#package_help:
#        @$(ECHO)\
#"\nAdditional notes for wmFur:\n\
# If you are in a working copy (or 'wip' is set), then wmFur will be installed\n\
# below these directories:\n\
#\t\t$(wetaos)/$(MAYANEED)/plug-ins/wetaWIP/wmFur-"'$$(VERSION)'"\n\
#\t\t$(weta)/$(MAYANEED)/mel/wetaWIP/wmFur-"'$$(VERSION)'"\n\
# If you are in a tag (and 'wip' is not set) and INSTALL_AS_DEFAULT is not set,\n\
# then wmFur will be installed below these directories:\n\
#\t\t$(wetaos)/$(MAYANEED)/plug-ins/weta/wmFur-"'$$(VERSION)'"\n\
#\t\t$(weta)/$(MAYANEED)/mel/weta/wmFur-"'$$(VERSION)'"\n\
# If you are in a tag (and 'wip' is not set) and INSTALL_AS_DEFAULT is set,\n\
# then wmFur will be installed below these directories:\n\
#\t\t$(wetaos)/$(MAYANEED)/plug-ins/weta\n\
#\t\t$(weta)/$(MAYANEED)/mel/weta\n\
#Note: 'wip' is currently $(if $(WIP),set,not set).\n\
#Note: Versioned paths would currently have: wmFur-$(VERSION)"
#endif

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



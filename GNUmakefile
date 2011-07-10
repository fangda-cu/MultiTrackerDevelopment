##############################################################################
include /vol/bob/make/head.mk
##############################################################################

LINUX_COMPILER := I
WITHKITS += dev/Figaro/work/jaubry-work

lib_h:
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C BASim USEMAYA=1 $@
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C strandSim USEMAYA=1 $@

lib_d lib_o lib: lib_h
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C BASim USEMAYA=1 $@
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C strandSim USEMAYA=1 $@
		
lib_fresh lib_zap lib_clean lib_wipe:
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C BASim USEMAYA=1 $@		
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C strandSim USEMAYA=1 $@

all:	lib
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/BASimulator
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen/icons

dbg:	lib_d
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/BASimulator  $@
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen  $@
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen/icons

opt:	lib_o
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/BASimulator  $@
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen  $@
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen/icons $@

fresh:	lib_fresh
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/BASimulator  fresh
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen  fresh
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen/icons fresh

zap:	lib_zap
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/BASimulator zap
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen zap
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen/icons zap

clean:	lib_clean
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/BASimulator  clean
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen  clean
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen/icons clean

wipe:	lib_wipe
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/BASimulator  wipe
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen  wipe
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen/icons wipe

init: lib_h

all: dbg opt


override allow_inst	:=#
ifdef release
override allow_inst	:= 1
endif
ifdef WIP
override allow_inst	:= 1
endif
ifdef allow_inst
inst:		lib
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C Apps/WmBunsen inst
else
inst:		inst_help
endif

maya_help:
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C maya plg_help scr_help package_help

icons_help:
		$(MAKE) LINUX_COMPILER=$(LINUX_COMPILER) -C maya/icons scr_help package_help

inst_help:
		@$(ECHO)\
"\nAdditonal notes for installation:\
\n If you want to install a work-in-procress version, please use:\
\n\tmake wip=1 $(MAKECMDGOALS)\
\n or:\
\n\tenv wip=1 make $(MAKECMDGOALS)\
\n If you want to install a release version, please use:\
\n\tmake release=1 $(MAKECMDGOALS)\
\n or:\
\n\tenv release=1 make $(MAKECMDGOALS)\
\n If you want to install a release version as the default version, use:\
\n\tmake INSTALL_AS_DEFAULT=1 release=1 $(MAKECMDGOALS)\
\n or:\
\n\tenv INSTALL_AS_DEFAULT=1 release=1 make $(MAKECMDGOALS)\
\n Note: You can only do a release install from a tag.\n"

help:		maya_help prman_help inst_help

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

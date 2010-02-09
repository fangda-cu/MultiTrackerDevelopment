
USEMAYA		:= 1

BUILD := "$(shell perl -e '$$pwd = $$ENV{PWD}; $$time=localtime(); print "\\\"$$pwd $$time\\\"";')"

##############################################################################
include /vol/bob/make/head.mk
##############################################################################

VERSION		:= $(KVERMAJ)
SCR_INSD	:=#
ifdef HAVEMAYA
ifndef INSTALL_AS_DEFAULT
SCR_INSD	:= $(wmayaicon)/wmFigaro-$(VERSION)
endif
endif

SCR1		:= out_wmFigaroNode.xpm
SCR1_MAYAIC	:= 1
SCR1_INSD	:= $(SCR_INSD)

SCR2		:= out_wmFigRodNode.xpm
SCR2_MAYAIC	:= 1
SCR2_INSD	:= $(SCR_INSD)

SCR3		:= out_wmFigCollisionNode.xpm
SCR3_MAYAIC	:= 1
SCR3_INSD	:= $(SCR_INSD)

ifdef HAVEMAYA
package_help:
		@$(ECHO)\
"\nAdditional notes for wmBarberShop icons:\n\
 If you are in a working copy (or 'wip' is set), then the icons will be installed\n\
 below this directory:\n\
\t\t$(weta)/$(MAYANEED)/bitmaps/wmBarberShop-"'$$(VERSION)'"\n\
 If you are in a tag (and 'wip' is not set) and INSTALL_AS_DEFAULT is not set,\n\
 then the icons will be installed below this directory:\n\
\t\t$(weta)/$(MAYANEED)/bitmaps/wmBarberShop-"'$$(VERSION)'"\n\
 If you are in a tag (and 'wip' is not set) and INSTALL_AS_DEFAULT is set,\n\
 then the icons will be installed below this directory:\n\
\t\t$(weta)/$(MAYANEED)/bitmaps\n\
Note: 'wip' is currently $(if $(WIP),set,not set).\n\
Note: Versioned paths would currently have: wmBarberShop-$(VERSION)"
endif

##############################################################################
include /vol/bob/make/tail.mk
##############################################################################


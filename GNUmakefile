 # $Id: GNUmakefile 22 2009-12-22 12:36:46Z schaelic $
 # --------------------------------------------------------------
 # GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
 # --------------------------------------------------------------
 
 name := esr2015
 G4TARGET := $(name)
 G4EXLIB := true
 
 .PHONY: all
 all: lib bin
 
 include $(G4INSTALL)/config/architecture.gmk
 
 #Add ROOT options for compilation
# CPPFLAGS += `root-config --cflags`
# LDFLAGS  += `root-config --libs`
 CPPFLAGS  += $(shell $(ROOTSYS)/bin/root-config --cflags)
 EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --glibs)
 
 include $(G4INSTALL)/config/binmake.gmk
 

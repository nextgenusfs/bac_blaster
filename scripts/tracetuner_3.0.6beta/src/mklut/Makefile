# *************************************************************************
#  This file is part of TraceTuner, the DNA sequencing quality value,
#  base calling and trace processing software.
# 
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received (LICENSE.txt) a copy of the GNU General Public
#  License along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# ************************************************************************/

#
# Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
#
# $Id: Makefile,v 1.8 2008/11/27 13:03:32 gdenisov Exp $
#

include ../include.mk

.PHONY: all lut 
 
all:  lut 

lut:		$(RELDIR)/lut

INCTOOLSDIR = ../../..
INCJTCTTDIR = $(INCTOOLSDIR)/jtc_tt
INCDIR      = ../compute_qv
INCTRAINDIR = ../mktrain
CURDIR      = .
LIBS        = -lm -lz
INSTALLDIR  = /home/gdenisov/build
IOLIBDIR    = $(INSTALLDIR)/lib
INCIOLIBDIR = $(INSTALLDIR)/include/io_lib
INCCHKDIR   = ../mkchk

LUTOBJS		= $(OBJDIR)/lut.o $(OBJDIR)/select.o $(OBJDIR)/func_name.o \
                  $(OBJDIR)/get_thresholds.o $(OBJDIR)/check_data.o
TTLIB		= $(LIBDIR)/libtt.a
CFLAGS     += -I$(INCDIR)
CFLAGS     += -I$(INCTRAINDIR)
CFLAGS     += -I$(CURDIR)
CFLAGS     += -I$(INCIOLIBDIR) -I$(INCJTCTTDIR) -I$(CURDIR)
CFLAGS     += -I$(INCTOOLSDIR)/io_lib
CFLAGS     += -I$(INCTOOLSDIR)/io_lib/utils
CFLAGS     += -I$(INCCHKDIR)

$(OBJDIR)/%.o: $(INCCHKDIR)/%.c
	mkdir -p $(OBJDIR)
	$(COMPILE.c) $< -o $@

LINK.c = g++

$(RELDIR)/lut: $(LUTOBJS) $(TTLIB) 
	@mkdir -p $(RELDIR)
	$(LINK.c) $(LUTOBJS) -o $@ $(LIBS) $(TTLIB) $(IOLIB)

$(OBJDIR)/lut.o:	lut.c lut.h get_thresholds.h select.h func_name.h params.h \
			$(INCDIR)/Btk_atod.h $(INCDIR)/Btk_qv.h $(INCTRAINDIR)/train.h

$(DIRS):
	mkdir -p $@

clean:
	@/bin/rm -f $(LUTOBJS)
	@/bin/rm -f $(RELDIR)/lut

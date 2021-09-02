#    This file is part of sndsys, a digital signal processing system.
#    Copyright (c) 2004-07 Volker Schatz (noise at volkerschatz dot com).
#
#    sndsys is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    sndsys is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with sndsys; if not, see <http://www.gnu.org/licenses/>.

ifneq "$(shell which icc)" ""
  CC	  =   icc ${ICCFLAGS}
  CC_LINK =   icc
else
  CC	  =   gcc ${GCCFLAGS}
  CC_LINK =   gcc
endif

LIBS	=   -L/usr/X11R6/lib -lm -lX11 -lpthread -lasound
GCCFLAGS =  -O3 -c -std=gnu99 -DUNIX -DLINUX -funsigned-char -funroll-loops\
	    -fomit-frame-pointer -Wstrict-prototypes -Wmissing-prototypes
ICCFLAGS =  -O3 -c -std=c99 -DUNIX -DLINUX -funsigned-char \
	    -Wall -w1 -wd1572 -wd47 -we140 -we147
DEBUGOPT =  -g -ggdb 

SYSOBJ	=   sndwt sndmodel sndscore sndasync sndeff sndcron sndfilt sndmisc \
	    sndsink sndsource sndsys

TARGETS = $(DISTTARGETS)

DISTTARGETS = frame vrec piano distort model model2
DISTFILES = piano.xbm pianostring.asc README.SNDSYS LICENCE melody2.abc objectdoc.pl grepobj.pl sndsys.h sndtut.tex

BACKUPTGZ = sndsys$(shell date +%Y%m%d).tgz
DISTTRUNC = ../distrib
DISTDIR = sndsys_distrib$(shell date +%Y%m%d)

$(TARGETS): %: %.o $(addsuffix .o,$(SYSOBJ))
	     $(CC_LINK) $(LDFLAGS) $(LIBS) -o $@ $^

all:	    $(TARGETS)

clean:
	    rm -f $(TARGETS) $(addsuffix .o, $(TARGETS)) $(addsuffix .wav, $(TARGETS)) $(addsuffix .o, $(SYSOBJ))

%.o:	    %.c sndsys.h Makefile
	    $(CC) $(CFLAGS) $(DEBUGOPT) $<

%.wav:	    %
	    nice -19 ./$<

model2.wav:   melody2.abc

melody4.wav:	melody4.xyz

instrtest.o:  instruments/*/*.c

tut:	    sndtut.dvi

sndtut.dvi: sndtut.tex
	    latex sndtut.tex
	    latex sndtut.tex
	    latex sndtut.tex

pdftut:	    sndtut.pdf

sndtut.pdf: sndtut.tex
	    pdflatex sndtut.tex
	    pdflatex sndtut.tex
	    pdflatex sndtut.tex

backup:
	    cd .. ; tar czf $(BACKUPTGZ) \
	    sndsys/*.h sndsys/*.c sndsys/*.pl sndsys/README.SNDSYS \
	    sndsys/sndtut.tex \
	    sndsys/Makefile sndsys/oldstuff sndsys/instruments/*/*.c \
	    sndsys/piano.xbm sndsys/*.abc sndsys/*.xyz sndsys/*.asc sndsys/*.dat

dist:
	    rm -rf $(DISTTRUNC)/$(DISTDIR)
	    mkdir -p $(DISTTRUNC)/$(DISTDIR)
	    cp $(addsuffix .c, $(SYSOBJ)) $(addsuffix .c, $(DISTTARGETS)) \
	        $(DISTFILES) $(DISTTRUNC)/$(DISTDIR)
	    sed -e 's/^TARGETS[[:space:]]*=.*$$/TARGETS = $$(DISTTARGETS)/' \
		-e 's///' < Makefile > $(DISTTRUNC)/$(DISTDIR)/Makefile
	    cd $(DISTTRUNC) && tar czf $(DISTDIR).tgz $(DISTDIR)


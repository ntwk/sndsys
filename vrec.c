/*
    This file is part of sndsys, a digital signal processing system.
    Copyright (c) 2004-07 Volker Schatz (noise at volkerschatz dot com).

    sndsys is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    sndsys is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with sndsys; if not, see <http://www.gnu.org/licenses/>.
*/

#ifdef NO_OSS
#error vrec needs Open Sound System!  Make sure OSS is available and comment "#define NO_OSS" in sndsys.h.
#endif

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <float.h>
#include <ctype.h>
#include "sndsys.h"


char *command;
sndobj *top;
static void (*vrec_oldinthandler)(int);


void usagemessage();
double readnumarg(int argc, char **argv, int index);
void vrec_oninterrupt(int signr);
void vrec_atexit(void);


int main(int argc, char *argv[])
{
  sndobj *pickup;
  char *outfile;
  double duration, thresh, ampl;
  int linein, count, dashfile;

  command= argv[0];
  duration= (double)(LONG_MAX-1)/(double)SAMPRATE;
  thresh= 0.01;
  ampl= 1.0;
  outfile= "vrec.wav";
  dashfile= 0;
  linein= 0;
  for( count= 1; count< argc; ++count )
  {
    if( !strcmp(argv[count], "-h") || !strcmp(argv[count], "--help") )
      usagemessage();
    else if( !strcmp(argv[count], "-l") )
      linein= 1;
    else if( !strcmp(argv[count], "-t") ) {
      duration= readnumarg(argc, argv, count);
      ++count;
    }
    else if( !strcmp(argv[count], "-m") ) {
      thresh= readnumarg(argc, argv, count);
      ++count;
    }
    else if( !strcmp(argv[count], "-a") ) {
      ampl= readnumarg(argc, argv, count);
      ++count;
    }
    else if( !strcmp(argv[count], "--") )
      dashfile= 1;
    else if( *argv[count] == '-' && !dashfile ) {
      fprintf(stderr, "%s: Unknown option %s.  Use -- to give output files starting with '-'.\n", command, argv[count]);
      usagemessage();
    }
    else
      outfile= argv[count];
  }
  pickup= mike(linein, 1.0);
  top= writeout( outfile, linear(0, ampl, voiceact( thresh, 0.3, 
	  rms(c(20), pickup), timeout(duration, NULL, pickup) ) ) );
  vrec_oldinthandler= signal( SIGINT, vrec_oninterrupt );
  if( vrec_oldinthandler==SIG_ERR || vrec_oldinthandler==SIG_DFL || 
              vrec_oldinthandler==SIG_IGN )
    vrec_oldinthandler= NULL;
  sndexecute( 0.0, duration, top );
  top= NULL;  // to prevent vrec_atexit trying to call writeout's exit() again
}


double readnumarg(int argc, char **argv, int index)
{
  double value;

  if( ++index >= argc ) {
    fprintf(stderr, "%s: Command-line switch %s requires an argument.\n", command, argv[index-1]);
    usagemessage();
  }
  if( !isdigit(*argv[index]) ) {
    fprintf(stderr, "%s: Command-line switch %s requires a numerical argument.\n", command, argv[index-1]);
    usagemessage();
  }
  return strtod(argv[index], NULL);
}


void usagemessage()
{
  fprintf(stderr, "usage: %s [-t <max. seconds>] [-m <threshold>] "
      "[-a amplification] [-l] [--] [<output wav file>]\n"
      "Voice-activated recording, %d Hz stereo.  Records when the RMS of the "
      "signal\nis above the given threshold, with gaps of 0.3 seconds.  The "
      "default threshold\nis 0.01, which is pretty high for spoken words "
      "without a preamp.  The\namplification is applied after the thresholding."
    "  Recording ends at the given\ntimeout or on interrupt signal (Ctrl-C).  "
      "The default output file name is\nvrec.wav.  -l selects line input "
      "instead of microphone.\n",
      command, SAMPRATE );
  exit(1);
}


void vrec_oninterrupt(int signr)
{
  top->exit(top);
  if( vrec_oldinthandler )
    vrec_oldinthandler(signr);
  else
    exit(-1);
}

void vrec_atexit(void)
{
  if( top )
    top->exit(top);
}


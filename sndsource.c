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

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <glob.h>
#include <sys/ioctl.h>
#include <unistd.h>
#ifndef NO_OSS
#include <linux/soundcard.h>
#endif

#include "sndsys.h"


/*=============================================================================
    harmonic [tab] <freq> <ampl>

Generates a chord of sine oscillations.  [tab] points to a table which contains
repeatedly frequency multiplier, amplitude multiplier and phase divided by 2pi
and is terminated by an entry with zero frequency (or END).

See also: sine, rect, saw
=============================================================================*/

struct harmonic {
  float     *tab, *phase, *sine;    // phase[..] is phase divided by 2pi
  long      sinelen;
};

float calcharmonic(sndobj *me, sndobj *, int);
void skipharmonic(sndobj *me, sndobj *, int);

sndobj *harmonic( float *tab, sndobj *freq, sndobj *ampl )
{
  sndobj *p;
  struct harmonic *d;
  float *har;
  char *sharedname;
  int nhar;
  
  p= newsndo( calcharmonic, "harmonic", "harmonic", 1, 2, freq, ampl );
  p->skip= skipharmonic;
  p->private[0]= d= (struct harmonic *)malloc(sizeof(struct harmonic));
  d->tab= tab;
  for( har= d->tab, nhar= 0; *har!=END && *har; har += 3, ++nhar );
  *har= 0.0;
  p->private[1]= d->phase= (float*)malloc(nhar*sizeof(float));
  for( har= d->tab, nhar= 0; *har; har += 3, ++nhar )
    d->phase[nhar]= har[2] - floorf(har[2]);
  d->sinelen= SAMPRATE;
  sharedname= news(char, NAMESTRLEN);
  snprintf(sharedname, NAMESTRLEN, "sinetab%ld", d->sinelen );
  if( (d->sine= findshared(sharedname, 1))!=NULL )
    free(sharedname);
  else {
    d->sine= sinelookup( d->sinelen );
    registershared( d->sine, sharedname );
  }
  p->private[2]= d->sine;
  return p;
}


float calcharmonic(sndobj *me, sndobj *caller, int innr)
{
  struct harmonic *d;
  float *har, *phase;
  float basefreq, baseampl, outval;
  
  d= (struct harmonic *)me->private[0];
  basefreq= INCALC(0);
  baseampl= INCALC(1);
  phase= d->phase;
  outval= 0.0;
  for( har= d->tab; *har; har += 3, ++phase )
  {
    outval += baseampl * har[1] *
		lininterpol( d->sine, d->sinelen, *phase*d->sinelen );
    *phase += basefreq * har[0] / (float)SAMPRATE;
    if( *phase > 1.0 )
      *phase -= 1.0;
  }
  OUTPUT(0, outval);
  CALCRETURN;
}


void skipharmonic(sndobj *me, sndobj *caller, int innr)
{
  struct harmonic *d;
  float *har, *phase;
  float basefreq;
  
  d= (struct harmonic *)me->private[0];
  basefreq= INCALC(0);
  INSKIP(1);
  phase= d->phase;
  for( har= d->tab; *har; har += 3, ++phase )
  {
    *phase += basefreq * har[0] / (float)SAMPRATE;
    if( *phase > 1.0 )
      *phase -= 1.0;
  }
}


/*=============================================================================
    sine (phase/2pi) <freq> <ampl>

Generates a sine oscillation.  The float argument phase_2pi determines the
initial phase divided by 2 pi (360 degrees).  The following argument, the first
input, is the frequency, the second the amplitude.

See also: harmonic, rect, saw
=============================================================================*/

sndobj *sine( float phase_2pi, sndobj *freq, sndobj *ampl )
{
  sndobj *p;
  float *tab;
  
  tab= (float*)malloc(4*sizeof(float));
  tab[0]= tab[1]= 1.0;
  tab[2]= phase_2pi;
  tab[3]= 0.0;
  p= harmonic( tab, freq, ampl );
  p->name= "sine";
  p->private[PLENTY-1]= tab;
  return p;
}


/*=============================================================================
    rect (phase) <freq> <ampl> <ratio>

Creates a rectangular waveform with frequency <freq> which alternates between
-<ampl> and <ampl>.  (phase) is starting phase / 2pi.  Phase=0 is defined as
the point where the output signal becomes high.  <ratio> should be between -1
and 1.  -1 means constantly low except one value per period, 1 the opposite.  0
is symmetric.  Values outside this range give constant output.

See also: sine, saw
=============================================================================*/

float calcrect(sndobj *me, sndobj *, int);

sndobj *rect( float phase, sndobj *freq, sndobj *ampl, sndobj *ratio )
{
  sndobj *p;
  double *ph;
  
  p= newsndo( calcrect, "rect", "rect", 1, 3, freq, ampl, ratio );
  p->private[0]= ph= (double*)malloc(sizeof(double));
  *ph= fmod(phase, 1.0);
  if( *ph< 0.0 )
    *ph += 1.0;
  return p;
}

float calcrect(sndobj *me, sndobj *caller, int innr)
{
  double *ph;
  double bound;
  float ratio, f, lastf, amp, wlen;
  
  ph= (double*)me->private[0];
  lastf= INPUT(0);
  f= INCALC(0);
  amp= INCALC(1);
  ratio= INCALC(2);
  if( !f )
    CALCRETURN;
  if( lastf )
    wlen= (float)SAMPRATE/lastf;
  else
    wlen= (float)SAMPRATE/f;
  bound= (1.0 + 0.5*(1.0+ratio)*(wlen-2.0))/wlen;
  if( ratio > 1.0 )
    OUTPUT(0, amp);
  else if( ratio < -1.0 )
    OUTPUT(0, -amp);
  else if( *ph<=bound )
    OUTPUT(0, amp);
  else
    OUTPUT(0, -amp);
  *ph += f/SAMPRATE;
  while( *ph>1.0 )
    *ph -= 1.0;
  CALCRETURN;
}


/*=============================================================================
    saw (phase) <freq> <ampl> <ratio>

Creates a sawtooth waveform with frequency <freq> and amplitude <ampl>.
(phase) is starting phase / 2pi.  Phase=0 is defined as the zero crossing on
the rising part of the sawtooth.  <ratio> should be between -1 and 1.  1 means
constant rise with immediate fall, -1 the opposite.  0 is symmetric.

See also: sine, rect
=============================================================================*/

float calcsaw(sndobj *me, sndobj *, int);

sndobj *saw( float phase, sndobj *freq, sndobj *ampl, sndobj *ratio )
{
  sndobj *p;
  double *ph;
  
  p= newsndo( calcsaw, "saw", "saw", 1, 3, freq, ampl, ratio );
  p->private[0]= ph= (double*)malloc(sizeof(double));
  *ph= fmod(phase, 1.0);
  if( *ph< 0.0 )
    *ph += 1.0;
  return p;
}

float calcsaw(sndobj *me, sndobj *caller, int innr)
{
  double *ph;
  double bound;
  float ratio, f, amp, outval;
  
  ph= (double*)me->private[0];
  f= INCALC(0);
  amp= INCALC(1);
  ratio= INCALC(2);
  if( !f )
    CALCRETURN;
  if( ratio > 1.0 )
    ratio= 1.0;
  else if( ratio < -1.0 )
    ratio= -1.0;
  bound= (1.0+ratio)/4.0;
  if( *ph<=bound )
    outval= *ph/bound;
  else if( *ph>1.0-bound )
    outval= -(1.0-*ph)/bound;
  else
    outval= (0.5-*ph)/(0.5-bound);
  OUTPUT(0, amp*outval);
  *ph += f/SAMPRATE;
  while( *ph>1.0 )
    *ph -= 1.0;
  CALCRETURN;
}


/*=============================================================================
    file "name" <freq> ...

Sample or waveform in .wav or ASCII or XFig file.  The ASCII file can be a
single column (successive sampled values with implied sampling rate SAMPRATE,
extension .asc) or two columns (time in seconds and value, extension .dat).
The two columns of .dat files are taken to be time in seconds and value.  They
should be in ascending order in time.  The highest time wraps around to the
lowest and should therefore have the same value.  The XFig file must contain a
polyline object drawn from left to right which determines the wave form.  The
optional (may be NULL) input <freq> serves to change the frequency; if a
normalisation of <freq> is required, the resampling factor (see below) can be
set.  For .asc, .dat and .fig files containing one period of a waveform, 
<freq>==1 means 1 Hz, but if <freq> is not given (NULL), the .asc file values
are output one at a time (at the sampling rate), and the time coordinate of
.dat files will be assumed to be in seconds.  The file name may be prepended by
special characters, which correspond to additional arguments:

*   resampling factor   double  >1 for higher pitch, <1 for lower.
+   offset              double  [0,1].  Start offset in file/part of file.
				Phase/2pi for waveforms.
only for .wav:
<   time range      2 x double  Use only part of file, starting at first 
				argument, ending at second.  Negative values
				are relative to end of file.  0 as second value
				means end of file.
#   value range     2 x int     Same for sampling value indices.
$   stop		int	Only play once.  Stop sndexecute by setting
				stoptime so many samples after end of playback.

See also: oneshot, rapidfire, multiwave, writeout
=============================================================================*/

typedef enum { wavfile, ascfile, datfile, figfile } filetype;
struct file {
  const char    *name;
  filetype      type;
  void          *buf, *pos;
  double        factor, fpos;
  long          length;
  int		stop, stopdelay;
};

float calcwavfile(sndobj *me, sndobj *, int);
void skipwavfile(sndobj *me, sndobj *, int );
float *parseascfile( FILE *handle, char *sharedname );
float calcascfile(sndobj *me, sndobj *, int);
void skipascfile(sndobj *me, sndobj *, int );
int datcompare(const void *a, const void *b);
float calcdatfile(sndobj *me, sndobj *, int);
void skipdatfile(sndobj *me, sndobj *, int );

sndobj *file(const char *name, sndobj *freq, ...)
{
  sndobj *p;
  struct file *d;
  FILE *handle;
  wavheader *head;
  double *dblptr;
  long *lngptr;
  int *intptr;
  const char *filename, *ext, *argstr;
  char *txtbuf, *line, *eol, *endconv, *sharedname;
  unsigned char *cptr;
  signed short *sptr;
  float *fptr;
  va_list varargs;
  filetype type;
  double ffrom, fto;
  long flen, lcount;
  int commentline, ifrom, ito, fileerr;
  
  //  determine file type
  for( filename= name; *filename=='+' || *filename=='*' || *filename=='<' ||
       *filename=='#' || *filename=='$'; ++filename );
  for( ext= filename; *ext; ++ext );
  if( (ext -= 4) >= filename )
    if( !strcmp(ext, ".wav") )
      type= wavfile;
    else if( !strcmp(ext, ".asc") )
      type= ascfile;
    else if( !strcmp(ext, ".dat") )
      type= datfile;
    else if( !strcmp(ext, ".fig") || !strcmp(ext, ".fig.bak") )
      type= figfile;
    else {
      fprintf(stderr, "Warning: File `%s' has unknown extension. Returning 0.\n", filename );
      return c(0.0);
    }
  else {
    fprintf(stderr, "Warning: File `%s' has unknown extension. Returning 0.\n", filename );
    return c(0.0);
  }

  //  open file and allocate common structures
  handle= fopen( filename, "r" );
  if( !handle ) {
    fprintf(stderr, "Warning: Could not open file `%s'. Returning 0.\n", filename );
    return c(0.0);
  }
  d= (struct file *)malloc(sizeof(struct file));
  d->name= filename;
  d->type= type;
  d->factor= 1.0;
  d->fpos= 0.0;
  d->stop= 0;
  if( !freq )
    p= newsndo(NULL, name, "file", 1, 0);
  else
    p= newsndo(NULL, name, "file", 1, 1, freq);
  p->private[0]= d;

  //  read additional arguments
  ffrom= fto= 0.0;
  ifrom= ito= 0;
  va_start(varargs, freq);
  for( argstr= name; argstr< filename; ++argstr )
    if( *argstr == '*' ) {
      d->factor= va_arg(varargs, double);
      if( !finite(d->factor) || d->factor <= 0.0 )
	d->factor= 1.0;
    }
    else if( *argstr == '+' ) {
      d->fpos= va_arg(varargs, double);
      if( !finite(d->fpos) )
	d->fpos= 0.0;
      else if( d->fpos < 0.0 )
	d->fpos= 1.0 + fmod(d->fpos, 1.0);
      else if( d->fpos > 1.0 )
	d->fpos= fmod(d->fpos, 1.0);
    }
    else if( *argstr == '<' ) {
      ffrom= va_arg(varargs, double);
      fto= va_arg(varargs, double);
      if( !finite(ffrom) || !finite(fto) )
	ffrom= fto= 0.0;
    }
    else if( *argstr == '#' ) {
      ifrom= va_arg(varargs, int);
      ito= va_arg(varargs, int);
    }
    else if( *argstr == '$' ) {
      d->stop= 1;
      d->stopdelay= va_arg(varargs, int);
    }
  va_end(varargs);

  //  read file and set up structures
  fileerr= 0;
  switch( type )
  {
    case wavfile:
      p->calc= calcwavfile;
      p->skip= skipwavfile;
      sharedname= news(char, 2*NAMESTRLEN);
      snprintf( sharedname, 2*NAMESTRLEN, "%s, *%g +%g <(%g,%g) #(%d,%d)", filename, d->factor, d->fpos, ffrom, fto, ifrom, ito );
      if( (d->buf= findshared(sharedname, 1))!=NULL )
      {
	free(sharedname);
	p->private[1]= d->buf;
	dblptr= (double*)d->buf;
	d->factor *= *dblptr++;
	lngptr= (long*)dblptr;
	d->length= *lngptr++;
	intptr= (int*)lngptr;
	p->nch= *intptr++;
	d->buf= intptr;
      }
      else
      {
	head= (wavheader *)malloc(sizeof(wavheader));
	fread(head, sizeof(wavheader), 1L, handle);
	if( ferror(handle) )
	  { fileerr= 1; break; }
	if( head->RIFF!='FFIR' || head->WAVE!='EVAW' || head->fmt!=' tmf' ) {
	  fprintf(stderr, "Warning: Wrong .wav header in `%s'. Returning 0.\n", d->name );
	  fileerr= 1;
	  break;
	}
	if( head->audiofmt!=1 || head->bitspersample==0 || head->bitspersample>16
	    || head->blockalign!=head->nchannels*((head->bitspersample+7)/8) ) {
	  fprintf(stderr, "Warning: Shocking parameters in `%s': Format tag %hu,"
	      " %hu bits per sample, %hu channels, block alignment %hu. "
	      "Returning 0.\n", d->name, head->audiofmt, head->bitspersample, 
	      head->nchannels, head->blockalign );
	  fileerr= 1;
	  break;
	}
	while( head->data!='atad' && !feof(handle) && !ferror(handle) ) {
	  fseek( handle, head->restsize, SEEK_CUR );
	  fread( &head->data, sizeof(unsigned long), 2L, handle );
	}
	if( feof(handle) || ferror(handle) ) {
	  if( feof(handle) )
	    fprintf(stderr, "Warning: No data chunk in `%s'. Returning 0.\n", 
		      d->name);
	  else if( ferror(handle) )
	    fprintf(stderr, "Warning: Error searching for data chunk in `%s'. "
		      "Returning 0.\n", d->name);
	  fileerr= 1;
	  break;
	}
	
	if( ffrom || fto ) {
	  ifrom= (int)round(ffrom*head->samplerate);
	  ito= (int)round(fto*head->samplerate);
	}
	if( ifrom || ito ) {
	  if( ifrom < 0 )
	    ifrom += head->restsize/head->blockalign;
	  if( ito <= 0 )
	    ito += head->restsize/head->blockalign;
	  if( ito > head->restsize/head->blockalign ) {
	    fprintf(stderr, "Warning: Range endpoint for file `%s' too large, replaced by end of file.\n", d->name );
	    ito= head->restsize/head->blockalign;
	  }
	  if( ifrom>=ito ) {
	    fprintf(stderr, "Warning: Ignoring invalid range for file `%s'.\n", d->name );
	    d->length= head->restsize/head->blockalign;
	  }
	  else {
	    d->length= (long)(ito-ifrom);
	    fseek( handle, (long)(ifrom*head->blockalign), SEEK_CUR );
	  }
	}
	else
	  d->length= head->restsize/head->blockalign;
	
	p->nch= head->nchannels;
	d->factor *= (double)head->samplerate/(double)SAMPRATE;
	p->private[1]= d->buf= malloc(d->length*p->nch*2L + sizeof(long) + sizeof(double) + sizeof(int));
	dblptr= (double*)d->buf;
	*dblptr++ = (double)head->samplerate/(double)SAMPRATE;
	lngptr= (long*)dblptr;
	*lngptr++ = d->length;
	intptr= (int*)lngptr;
	*intptr++ = head->nchannels;
	d->buf= intptr;
	fread(d->buf, 1L, d->length*head->blockalign, handle);
	if( ferror(handle) ) {
	  fprintf(stderr, "Warning: Error reading data in `%s'. Returning 0.\n", d->name );
	  fileerr= 1;
	  break;
	}
	if( head->bitspersample<=8 )
	  for( cptr= (unsigned char *)d->buf + d->length*p->nch,
		  sptr= (signed short *)d->buf + d->length*p->nch;
	      cptr >= (unsigned char *)d->buf;  --cptr, --sptr )
	    *sptr= (*cptr-128) << 8;
	
	free(head);
	registershared( (char*)d->buf-sizeof(long)-sizeof(double)-sizeof(int), sharedname);
      }
      d->pos= d->buf;
      d->fpos *= d->length;
      break;
    
    case ascfile:
      p->calc= calcascfile;
      p->skip= skipascfile;
      sharedname= news(char, 2*NAMESTRLEN);
      snprintf(sharedname, 2*NAMESTRLEN, "%s, *%g +%g", filename, d->factor, d->fpos);
      p->private[1]= d->buf= parseascfile(handle, sharedname);
      d->length= *(long*)d->buf;
      fileerr= !d->buf;
      if( fileerr ) {
	fprintf(stderr, "Warning: Error parsing file `%s'.\n", d->name );
	break;
      }
      d->buf= d->pos= ((long*)d->buf+1);
      d->fpos *= d->length;
      if( ifrom || ito || ffrom || fto )
	fprintf( stderr, "Warning: Range for .asc file not yet implemented.\n" );
      if( p->nin )
	d->factor *= d->length/(double)SAMPRATE;
      break;

    case datfile:
      p->calc= calcdatfile;
      p->skip= skipdatfile;
      sharedname= news(char, 2*NAMESTRLEN);
      snprintf(sharedname, 2*NAMESTRLEN, "%s, *%g +%g", filename, d->factor, d->fpos);
      if( (d->buf= findshared(sharedname, 1))!=NULL )
      {
	free(sharedname);
	p->private[1]= d->buf;
	lngptr= (long*)d->buf;
	d->length= *lngptr++;
	d->buf= lngptr;
      }
      else
      {
	fseek( handle, 0L, SEEK_END );
	flen= ftell(handle);
	fseek( handle, 0L, SEEK_SET );
	txtbuf= (char *)malloc(flen+2L);
	fread(txtbuf, 1L, flen, handle);
	txtbuf[flen]= txtbuf[flen+1]= 0;
	if( ferror(handle) ) {
	  fprintf(stderr, "Warning: Error reading from `%s'. Returning 0.\n", d->name );
	  fileerr= 1;
	  free(txtbuf);
	  break;
	}
	flen= 0;
	for( line= txtbuf; *line; ++line )
	  if( *line == '\n' )
	    ++flen;
	p->private[1]= d->buf= malloc(flen*(sizeof(double)+sizeof(float)) + sizeof(long));
	d->buf= (long*)d->buf + 1;
	fptr= (float*)d->buf;
	for( line= txtbuf; *line; line= eol+1 )
	{
	  for( eol= line; *eol && isspace(*eol) && *eol!='\n'; ++eol );
	  commentline= (*eol==';' || *eol=='#' || *eol=='\n' || !*eol);
	  while( *eol && *eol!='\n' )
	    ++eol;
	  if( commentline )
	    continue;
	  *eol= 0;
	  dblptr= (double*)fptr;
	  *dblptr++ = strtod( line, &endconv );
	  fptr= (float *)dblptr;
	  if( endconv>line )
	    *fptr++= strtof( (line= endconv), &endconv );
	  if( !finite( ((double*)(fptr-1))[-1] ) || !finite(fptr[-1]) || 
		  errno==ERANGE || endconv==line ) {
	    fprintf(stderr, "Warning: Error parsing file `%s'.\n", d->name );
	    fileerr= 1;
	    break;
	  }
	}
	free(txtbuf);
	if( fileerr )
	  break;
	if( ifrom || ito || ffrom || fto )
	  fprintf( stderr, "Warning: Range for .dat file not yet implemented.\n" );
	d->length= ((long*)d->buf)[-1]= ((char*)fptr-(char*)d->buf)/(sizeof(double)+sizeof(float));
	registershared((long*)d->buf-1, sharedname);
      }
      d->pos= (char*)d->buf + (sizeof(double) + sizeof(float));
      ffrom= *(double*)d->buf;
      fto= *(double*)((char*)d->buf+(d->length-1)*(sizeof(double)+sizeof(float)));
      for( fptr= (float*)d->buf, lcount= d->length; lcount> 0; --lcount ) {
	*(double*)fptr -= ffrom;
	fptr= (float*)((char*)fptr + sizeof(double)+sizeof(float));
      }
      d->fpos *= fto-ffrom;
      if( !p->nin )
	d->factor /= (double)SAMPRATE;
      else      //  divide by total time interval
	d->factor *= (fto-ffrom) / (double)SAMPRATE;
      break;

    case figfile:
      p->calc= calcdatfile;
      p->skip= skipdatfile;
      fseek( handle, 0L, SEEK_END );
      flen= ftell(handle);
      fseek( handle, 0L, SEEK_SET );
      txtbuf= (char *)malloc(flen+2L);
      fread(txtbuf, 1L, flen, handle);
      txtbuf[flen]= txtbuf[flen+1]= 0;
      if( ferror(handle) ) {
	fprintf(stderr, "Warning: Error reading from `%s'. Returning 0.\n", d->name );
	fileerr= 1;
	free(txtbuf);
	break;
      }
      flen= 0;
      for( line= txtbuf; *line && (*line!='2'|| line[1]!=' ' || line[2]!='1' || line[3]!=' '); ++line )
	while( *line && *line!='\n' )
	  ++line;
      if( !*line ) {
	fprintf(stderr, "Warning: No polyline object found in `%s'. Returning 0.\n", d->name );
	fileerr= 1;
	free(txtbuf);
	break;
      }
      for( eol= line; *eol && *eol!='\n'; ++eol );
      if( eol-line < 31 ) {
	fprintf(stderr, "Warning: Wrong format in polyline object in `%s': too few number fields in main line. Returning 0.\n", d->name );
	fileerr= 1;
	free(txtbuf);
	break;
      }
      while( isspace(*eol) )    --eol;
      while( isdigit(*eol) )    --eol;
      while( isspace(*eol) )    --eol;
      while( isdigit(*eol) )    --eol;
      while( isspace(*eol) )    --eol;
      while( isdigit(*eol) )    --eol;
      ++eol;
      commentline= 0;
      if( *eol!='0' || eol[1]!=' ' )
	++commentline;          // additional line for arrow specification
      while( isdigit(*eol) )    ++eol;
      while( isspace(*eol) )    ++eol;
      if( *eol!='0' || eol[1]!=' ' )
	++commentline;          // additional line for arrow specification
      while( isdigit(*eol) )    ++eol;
      flen= strtol(eol, NULL, 10);
      if( flen<2 ) {
	fprintf(stderr, "Warning: file `%s': number of points in polyline object unreadable or <2. Returning 0.\n", d->name );
	fileerr= 1;
	free(txtbuf);
	break;
      }
      d->length= flen;
      p->private[1]= d->buf= malloc(flen*(sizeof(double)+sizeof(float)));
      line= eol;
      for( ; commentline>=0; --commentline )
	for( ++line; *line && *line!='\n'; ++line );
      for( fptr= (float*)d->buf; flen> 0 && *line; --flen ) {
	*((double*)fptr) = (double)strtol(line, &line, 10);
	fptr= (float*)((double*)fptr + 1);
	*fptr++ = (float)strtol(line, &line, 10);
      }
      if( flen>0 ) {
	fprintf(stderr, "Warning: Encountered premature end of file `%s' while reading polyline points. Returning 0.\n", d->name );
	fileerr= 1;
	free(txtbuf);
	break;
      }
      free(txtbuf);
      if( fileerr )
	break;
      if( ifrom || ito || ffrom || fto )
	fprintf( stderr, "Warning: Range for .fig file not yet implemented.\n" );
      d->pos= (char*)d->buf + (sizeof(double) + sizeof(float));
      ffrom= *(double*)d->buf;
      fto= *(double*)((char*)d->buf+(d->length-1)*(sizeof(float)+sizeof(double)));
      if( ffrom < fto ) {
	for( fptr= (float*)d->buf, lcount= d->length; lcount> 0; --lcount ) {
	  *(double*)fptr -= ffrom;
	  fptr= (float*)((char*)fptr + sizeof(double)+sizeof(float));
	}
	d->fpos *= fto-ffrom;
	d->factor *= (fto-ffrom)/(double)SAMPRATE;
      }
      else {        //  assume reverse order
	for( fptr= (float*)d->buf, lcount= d->length; lcount> 0; --lcount ) {
	  *(double*)fptr = ffrom - *(double*)fptr;
	  fptr= (float*)((char*)fptr + sizeof(double)+sizeof(float));
	}
	d->fpos *= ffrom-fto;
	d->factor *= (ffrom-fto)/(double)SAMPRATE;
      }
      // normalise y values (-> in [-1,1]):
      ffrom= DBL_MAX;
      fto= -DBL_MAX;
      for( fptr= (float*)((double*)d->buf+1), lcount= d->length; lcount> 0; --lcount ) {
	if( *fptr>fto )
	  fto= *fptr;
	if( *fptr<ffrom )
	  ffrom= *fptr;
	fptr= (float*)((char*)fptr + sizeof(double)+sizeof(float));
      }
      if( ffrom<fto )
	for( fptr= (float*)((double*)d->buf+1), lcount= d->length; lcount> 0; --lcount ) {
	  *fptr= -1.0 + 2.0 * (*fptr-fto)/(ffrom-fto);
			// scale and invert - xfig has (0,0) in upper left corner
	  fptr= (float*)((char*)fptr + sizeof(double)+sizeof(float));
	}
      else
	for( fptr= (float*)((double*)d->buf+1), lcount= d->length; lcount> 0; --lcount ) {
	  *fptr= 0.0;
	  fptr= (float*)((char*)fptr + sizeof(double)+sizeof(float));
	}
      break;
  }
  
  if( fileerr )
  {
    errno= 0;
    return c(0.0);
  }
  fclose(handle);
  return p;
}


float calcwavfile(sndobj *me, sndobj *caller, int innr)
{
  struct file *d;
  signed short *read, *nextread;
  float val, nextval;
  int ch;
  
  d= (struct file *)me->private[0];
  if( d->stop< 0 ) {
    if( --d->stopdelay < 0 )
      stoptime= -DBL_MAX;
    CALCRETURN;
  }
  read= d->pos;
  read += (long)floor(d->fpos) * me->nch;
  if( read >= (short *)d->buf+d->length*me->nch )
    if( d->stop> 0 ) {
      d->stop= -1;
      FORCH
	OUTPUT(ch, 0.0f);
      if( --d->stopdelay < 0 )
	stoptime= -DBL_MAX;
      CALCRETURN;
    }
    else
      do {
	read -= d->length*me->nch;
      } while(read >= (short *)d->buf+d->length*me->nch);
  d->fpos -= floor(d->fpos);
  if( !d->fpos )
    for( ch= 0; ch< me->nch; ++ch )
      OUTPUT(ch, (((float)*read++)+32768.0)/65535.0*2.0 - 1.0 );
  else {
    nextread= read + me->nch;
    if( nextread >= (signed short *)d->buf+d->length*me->nch )
      nextread -= d->length*me->nch;
    for( ch= 0; ch< me->nch; ++ch ) {
      val= (((float)*read++)+32768.0)/65535.0*2.0 - 1.0;
      nextval= (((float)*nextread++)+32768.0)/65535.0*2.0 - 1.0;
      OUTPUT(ch, val + d->fpos*(nextval-val) );
    }
  }
  d->pos= read - me->nch;
  if( !me->nin )
    d->fpos += d->factor;
  else
    d->fpos += d->factor * INCALC(0);
  CALCRETURN;
}


void skipwavfile(sndobj *me, sndobj *caller, int innr)
{
  struct file *d;

  d= (struct file *)me->private[0];
  if( !me->nin )
    d->fpos += d->factor;
  else
    d->fpos += d->factor * INCALC(0);
  if( d->fpos > (double)d->length )
    d->fpos -= (double)d->length;
}


float *parseascfile( FILE *handle, char *sharedname )
{
  float *buf, *fptr;
  long *length;
  char *txtbuf, *line, *eol, *endconv;
  long flen;
  int commentline;
  
  if( (buf= findshared(sharedname, 1))!=NULL ) {
    free(sharedname);
    return buf;
  }
  fseek( handle, 0L, SEEK_END );
  flen= ftell(handle);
  fseek( handle, 0L, SEEK_SET );
  txtbuf= (char *)malloc(flen+2L);
  fread(txtbuf, 1L, flen, handle);
  txtbuf[flen]= txtbuf[flen+1]= 0;
  if( ferror(handle) ) {
    fprintf(stderr, "Warning: Error reading from .asc file. Returning 0.\n" );
    free(txtbuf);
    return NULL;
  }
  flen= 1;
  for( line= txtbuf; *line; ++line )
    if( *line == '\n' )
      ++flen;
  buf= malloc(flen*sizeof(float) + sizeof(long));
  length= (long*)buf;
  buf= (float*)(length+1);
  fptr= buf;
  for( line= txtbuf; *line; line= eol+1 )
  {
    for( eol= line; *eol && isspace(*eol) && *eol!='\n'; ++eol );
    commentline= (*eol==';' || *eol=='#' || *eol=='\n' || !*eol);
    while( *eol && *eol!='\n' )
      ++eol;
    if( commentline )
      continue;
    *eol= 0;
    *fptr++= strtof( line, &endconv );
    if( !finite(fptr[-1]) || errno==ERANGE || endconv==line ) {
      fprintf(stderr, "Warning: Error parsing .asc file.\n" );
      free(length);
      free(txtbuf);
      return NULL;
    }
  }
  free(txtbuf);
  *length= fptr-buf;
  registershared( length, sharedname );
  return (float*)length;
}

float calcascfile(sndobj *me, sndobj *caller, int innr)
{
  struct file *d;
  float *next;
  
  d= (struct file *)me->private[0];
  d->pos= (float*)d->pos + (long)floor(d->fpos);
  while( (float*)d->pos >= (float*)d->buf + d->length )
    d->pos= (float*)d->pos - d->length;
  d->fpos -= floor(d->fpos);
  if( !d->fpos )
    OUTPUT(0, *(float*)d->pos);
  else {
    next = (float*)d->pos + 1;
    if( next >= (float*)d->buf + d->length )
      next -= d->length;
    OUTPUT(0, *(float*)d->pos + d->fpos*(*next-*(float*)d->pos));
  }
  if( !me->nin )
    d->fpos += d->factor;
  else
    d->fpos += d->factor * INCALC(0);
  CALCRETURN;
}

void skipascfile(sndobj *me, sndobj *caller, int innr)
{
  struct file *d;

  d= (struct file *)me->private[0];
  if( !me->nin )
    d->fpos += d->factor;
  else
    d->fpos += d->factor * INCALC(0);
  if( d->fpos > (double)d->length )
    d->fpos -= (double)d->length;
}


int datcompare(const void *a, const void *b)
{
  double diff;
  
  diff= *(double*)a - *(double*)b;
  if( diff < 0.0 )      return -1;
  else if( diff > 0.0 ) return 1;
  else                  return 0;
}

/* 
.dat file: two columns, time (in seconds) and value; need not be equidistant
but should be in scending order; highest time will be interpreted as equal to 0
for wraparound and should have the same value.  .fig file: The polyline should
be drawn from left to right (right to left will be reversed).  */

float calcdatfile(sndobj *me, sndobj *caller, int innr)
{
  struct file *d;
  char *read, *previous;
  
  d= (struct file *)me->private[0];
  read= (char *)d->pos;
  while( *(double*)read < d->fpos ) {
    read += sizeof(double)+sizeof(float);
    if( read >= ((char *)d->buf) + d->length*(sizeof(double)+sizeof(float)) ) {
      d->fpos -= *(double*)(read - (sizeof(double)+sizeof(float)));
      read= (char*)d->buf + (sizeof(double)+sizeof(float));
    }
  }
  previous= read-(sizeof(double)+sizeof(float));
  OUTPUT(0, *(float*)(previous+sizeof(double)) +
	    (*(float*)(read+sizeof(double)) - *(float*)(previous+sizeof(double))) *
	    (d->fpos - *(double*)previous)/(*(double*)read - *(double*)previous) );
  d->pos= read;
  if( !me->nin )
    d->fpos += d->factor;
  else
    d->fpos += d->factor * INCALC(0);
  CALCRETURN;
}

void skipdatfile(sndobj *me, sndobj *caller, int innr)
{
  struct file *d;

  d= (struct file *)me->private[0];
  if( !me->nin )
    d->fpos += d->factor;
  else
    d->fpos += d->factor * INCALC(0);
  if( d->fpos > *(double*)((char*)d->buf + (d->length-1)*(sizeof(double)+sizeof(float))) ) {
    d->fpos -= *(double*)((char*)d->buf + (d->length-1)*(sizeof(double)+sizeof(float)));
    d->pos= d->buf;
  }
}


/*=============================================================================
    oneshot "filename" <speed> <trigger>

Plays a file (see above) once every time it is triggered by <trigger>!=0.
<speed> may be NULL as <freq> in file.  However, the behaviour of .asc and .dat
files does not change as in the file object when <speed> is non-NULL (a single
period of a frequency makes little sense).  Rather, <speed> always acts as a
simple resampling factor.  It is read whenever oneshot is triggered and remains
constant for that time.  Only for .fig files is speed a frequency on which one
period with the defined waveform is output.  The value of the <trigger> channel
gives the volume of the sample.  At the next trigger, the previous playback is
cut off.

See also: file, rapidfire
=============================================================================*/

struct oneshot {
  double    basefactor;
  float	    volume;
  double    basecount;
  long      count;
};

float calconeshot(sndobj *me, sndobj *, int);
void skiponeshot(sndobj *me, sndobj *, int);

sndobj *oneshot( const char *filename, sndobj *speed, sndobj *trigger )
{
  sndobj *p, *fileobj;
  struct oneshot *d;
  struct file *filed;
  double *last;

  if( *filename=='*' || *filename=='+' || *filename=='<' || *filename=='#' ) {
    printf("Warning: filename modifiers not allowed with oneshot.  Discarding modifiers.\n");
    while( *filename=='*' || *filename=='+' || *filename=='<' || *filename=='#' )
      ++filename;
  }
  fileobj= file(filename, NULL);
  if( strcmp("file", fileobj->type ) )
    return fileobj;             // error in file()
  if( speed )
    p= newsndo( calconeshot, filename, "oneshot", fileobj->nch, 3, fileobj, trigger, speed);
  else
    p= newsndo( calconeshot, filename, "oneshot", fileobj->nch, 2, fileobj, trigger);
  p->skip= skiponeshot;
  p->private[0]= d= new(struct oneshot);
  filed= (struct file *)fileobj->private[0];
  d->basefactor= filed->factor;
  if( filed->type==figfile )
    d->basecount= SAMPRATE;
  else if( filed->type==datfile ) {
    last= (double*)((char*)filed->buf + (filed->length-1)*(sizeof(double)+sizeof(float)));
    d->basecount= (*last - *(double*)filed->buf)*SAMPRATE;
  }
  else
    d->basecount= filed->length-1;
  d->volume= 1.0;
  d->count= 0;
  return p;
}

float calconeshot(sndobj *me, sndobj *caller, int innr)
{
  struct oneshot *d;
  struct file *filed;
  float speed;
  int ch;

  d= (struct oneshot *)me->private[0];
  filed= (struct file *)me->in[0]->private[0];
  if( INCALC(1) ) {     // triggered
    d->volume= INPUT(1);
    filed= (struct file *)me->in[0]->private[0];
    if( filed->type==datfile || filed->type==figfile ) {
      filed->pos= (char*)filed->buf + sizeof(double)+sizeof(float);
      filed->fpos= 0.0;
    }
    else {
      filed->pos= filed->buf;
      filed->fpos= 0.0;
    }
    if( me->nin > 2 ) {
      speed= INCALC(2);
      if( !speed || !finite(speed) )
	speed= 1.0;
      filed->factor= d->basefactor * speed;
      d->count= (long)floor(d->basecount/speed)+1;
    }
    else {
      filed->factor= d->basefactor;
      d->count= (long)floor(d->basecount)+1;
    }
  }
  else if( me->nin> 2 )
    INSKIP(2);
  if( d->count ) {
    INCALC(0);
    FORCH
      OUTPUT(ch, d->volume*INPUTC(0, ch));
    --d->count;
  }
  else {    // no INSKIP(0) necessary since no one else has access to our file
    FORCH
      OUTPUT(ch, 0.0);
  }
  CALCRETURN;
}

void skiponeshot(sndobj *me, sndobj *caller, int innr)
{
  struct oneshot *d;
  struct file *filed;
  float speed;

  d= (struct oneshot *)me->private[0];
  if( INCALC(1) ) {     // triggered
    d->volume= INPUT(1);
    filed= (struct file *)me->in[0]->private[0];
    if( filed->type==datfile || filed->type==figfile ) {
      filed->pos= (char*)filed->buf + sizeof(double)+sizeof(float);
      filed->fpos= 0.0;
    }
    else {
      filed->pos= filed->buf;
      filed->fpos= 0.0;
    }
    if( me->nin > 2 ) {
      speed= INCALC(2);
      if( !speed || !finite(speed) )
	speed= 1.0;
      filed->factor= d->basefactor * speed;
      d->count= (long)floor(d->basecount/speed)+1;
    }
    else {
      filed->factor= d->basefactor;
      d->count= (long)floor(d->basecount)+1;
    }
  }
  else if( me->nin> 2 )
    INSKIP(2);
  if( d->count ) {
    INSKIP(0);
    --d->count;
  }
}



/*=============================================================================
    rapidfire (nmax) "filename" <speed> <trigger>

This sndobj is similar to oneshot in that it plays a file when receiving a
non-zero value on its <trigger> input.  Unlike oneshot, it does not interrupt a
playback when the next trigger arrives but ignores new triggers until the
current playback is finished.  Besides, rapidfire can play several incarnations
of the file in parallel, adding them up.  The parameter (nmax) determines how
many.  So only if more than (nmax) triggers happen within the time it takes to
play the file once will any of them be ignored.  <speed> has the same effect as
in oneshot: It is a resampling factor, changing the speed of the playback from
one value per sampling interval for .wav and .asc files, and relative to the
timing given in .dat and .fig files.  It can be passed as NULL.  The value of
<trigger>, when non-zero, determines the volume of the file.

See also: file, oneshot
=============================================================================*/

struct onefire {
  float  volume;
  long	 count;
  struct file filed;
};
struct rapidfire {
  double    basefactor, basecount;
  int	    nmax, n;
  struct onefire *par;
};

float calcrapidfire(sndobj *me, sndobj *, int);
void skiprapidfire(sndobj *me, sndobj *, int);

sndobj *rapidfire( int nmax, const char *filename, sndobj *speed, sndobj *trigger )
{
  sndobj *p, *fileobj;
  struct rapidfire *d;
  struct file *origfiled;
  double *last;
  int count;

  if( *filename=='*' || *filename=='+' || *filename=='<' || *filename=='#' ) {
    printf("Warning: filename modifiers not allowed with rapidfire.  Discarding modifiers.\n");
    while( *filename=='*' || *filename=='+' || *filename=='<' || *filename=='#' )
      ++filename;
  }
  fileobj= file(filename, NULL);
  if( strcmp("file", fileobj->type ) )
    return fileobj;             // error in file()
  if( speed )
    p= newsndo( calcrapidfire, filename, "rapidfire", fileobj->nch, 3, fileobj, trigger, speed);
  else
    p= newsndo( calcrapidfire, filename, "rapidfire", fileobj->nch, 2, fileobj, trigger);
  p->skip= skiprapidfire;
  p->private[0]= d= new(struct rapidfire);
  origfiled= (struct file *)fileobj->private[0];
  d->basefactor= origfiled->factor;
  if( origfiled->type==figfile )
    d->basecount= SAMPRATE;
  else if( origfiled->type==datfile ) {
    last= (double*)((char*)origfiled->buf + (origfiled->length-1)*(sizeof(double)+sizeof(float)));
    d->basecount= (long)ceil((*last - *(double*)origfiled->buf)*SAMPRATE);
  }
  else
    d->basecount= origfiled->length-1;
  d->nmax= nmax;
  d->n= 0;
  p->private[1]= d->par= news(struct onefire, nmax);
  for( count= 0; count< nmax; ++count ) {
    d->par[count].volume= 1.0;
    d->par[count].count= 0;
    d->par[count].filed= *origfiled;
  }
  return p;
}

float calcrapidfire(sndobj *me, sndobj *caller, int innr)
{
  struct rapidfire *d;
  struct file *origfiled;
  float speed;
  int count, ch;

  d= (struct rapidfire *)me->private[0];
  if( INCALC(1) && d->n < d->nmax )      // triggered and able to play one more
  {
    for( count= 0; count< d->nmax && d->par[count].count; ++count );
    ++d->n;
    d->par[count].volume= INPUT(1);
    if( d->par[count].filed.type==datfile || 
	d->par[count].filed.type==figfile )
    {
      d->par[count].filed.pos= (char*)d->par[count].filed.buf +
				      sizeof(double)+sizeof(float);
      d->par[count].filed.fpos= 0.0;
    }
    else {
      d->par[count].filed.pos= d->par[count].filed.buf;
      d->par[count].filed.fpos= 0.0;
    }
    if( me->nin > 2 ) {
      speed= fabsf(INCALC(2));
      if( !speed || !finite(speed) )
	speed= 1.0;
      d->par[count].filed.factor= d->basefactor * speed;
      d->par[count].count= (long)floor(d->basecount/speed)+1;
    }
    else {
      d->par[count].filed.factor= d->basefactor;
      d->par[count].count= (long)floor(d->basecount)+1;
    }
  }
  else if( me->nin> 2 )
    INSKIP(2);
  FORCH
    OUTPUT(ch, 0.0);
  if( d->n > 0 )
  {
    origfiled= (struct file *)me->in[0]->private[0];
    for( count= 0; count< d->nmax; ++count )
      if( d->par[count].count )
      {
	me->in[0]->private[0]= &d->par[count].filed;
	INCALC(0);
	FORCH
	  me->ch[ch] += d->par[count].volume*INPUTC(0, ch);
	if( !--d->par[count].count )
	  --d->n;
      }
    me->in[0]->private[0]= origfiled;
    		// just to ensure deallocation on program termination
  }
  CALCRETURN;
}

void skiprapidfire(sndobj *me, sndobj *caller, int innr)
{
  struct rapidfire *d;
  float speed;
  int count, ch;

  d= (struct rapidfire *)me->private[0];
  if( INCALC(1) && d->n < d->nmax ) {
    for( count= 0; count< d->nmax && d->par[count].count; ++count );
    ++d->n;
    d->par[count].volume= INPUT(1);
    if( d->par[count].filed.type==datfile || 
	d->par[count].filed.type==figfile )
    {
      d->par[count].filed.pos= (char*)d->par[count].filed.buf +
				      sizeof(double)+sizeof(float);
      d->par[count].filed.fpos= 0.0;
    }
    else {
      d->par[count].filed.pos= d->par[count].filed.buf;
      d->par[count].filed.fpos= 0.0;
    }
    if( me->nin > 2 ) {
      speed= INCALC(2);
      if( !speed || !finite(speed) )
	speed= 1.0;
      d->par[count].filed.factor= d->basefactor * speed;
      d->par[count].count= (long)floor(d->basecount/speed)+1;
    }
    else {
      d->par[count].filed.factor= d->basefactor;
      d->par[count].count= (long)floor(d->basecount)+1;
    }
  }
  else if( me->nin> 2 )
    INSKIP(2);
  if( d->n > 0 )
    for( count= 0; count< d->nmax; ++count )
      if( d->par[count].count )
      {
	if( !--d->par[count].count )
	  --d->n;
      }
}


/*=============================================================================
    multiwave (dimension) [indextypes] [ascfiles] [positions] <freq> <par1>...

This objects interpolates linearly between different wave forms given in .asc
files.  The different wave forms are thought of as occupying a point in an
n-dimensional parameter space.  They should form a rectangular or cubic grid,
in which interpolation is performed in one direction after another.
(dimension) gives the dimension of the parameter space.  The array [indextypes]
contains 2*(dimension) integer values.  The first of the pair gives the type of
the interpolation: MWTYPE_LINEAR, linear interpolation; MWTYPE_WRAP, linear
interpolation with wrap around at -1 and 1 (the corresponding input is also
wrapped automatically); or MWTYPE_FREQ, where the input determining the
position in this direction is the frequency.  This type value can be OR-ed with
the flag MWFLAG_CONTINUE.  If it is set, and if the type is MWTYPE_LINEAR or
MWTYPE_FREQ, the interpolation is continued for coordinate values which lie
outside the interval of the waveform files.  Unless the type is MWTYPE_FREQ,
the second value of the pair gives the number of the relevant input after
<freq>.  It starts at 0.  MWTYPE_FREQ may be used only once.  The order of the
entries in [indextypes] determines the order in which the interpolation is
conducted: interpolation with respect to coordinates given first in [indextype]
is conducted last.  The following two arrays contain the names of the waveform
.asc files and their positions in the n-dimensional space.  [positions]
contains (dimension) values for each file name in [ascfiles].  The order of the
coordinates is as defined by [indextypes].  [ascfiles] has to be terminated by
an empty string or a NULL pointer.  <freq>, the frequency input, has to be
followed by (dimension) (if MWTYPE_FREQ is not used) or (dimension)-1 (if
MWTYPE_FREQ is used) additional parameter inputs.

See also: file
=============================================================================*/

struct multi_node {
  float     index;
  int       type, flags, inputnr, ngrandchildren, wavenr;
  struct multi_node *firstchild;
};
struct multiwave {
  double    fpos;
  int       dim, freqused, nfiles, nsuccess;
  int       *indextypes;
  int       rootchildren;
  struct multi_node *tree;
  float     **wave;
  long      *size;
  float     *cache;
};

int buildsubtree( struct multi_node **write, struct multiwave *d, 
    const int *indextypes, const float *positions, int depth );
int cmpnodeindex( const void *p1, const void *p2 );
float calcmultiwave(sndobj *me, sndobj *, int);
float evaluatesubtree( sndobj *me, struct multi_node *subtree, int children, int level );
float evaluatemultiwave( sndobj *me, int wavenr );
void skipmultiwave(sndobj *me, sndobj *, int);
void exitmultiwave(sndobj *me);

#define MULTIWAVE_MAXFILES  1000

sndobj *multiwave( int dimension, const int *indextypes, const char **ascfiles, const float *positions, sndobj *freq, ... )
{
  static int checkinputnr[PLENTY];
  sndobj *p;
  struct multiwave *d;
  struct multi_node *write;
  FILE *file;
  va_list varargs;
  float *parsedfile;
  char *sharedname;
  int naddinputs, count, fileind;
  
  if( dimension> PLENTY ) {
    fprintf(stderr, "Warning: multiwave: dimension (%d) too large (max= %d).  Returning 0.\n", dimension, PLENTY );
    return c(0);
  }
  p= newsndo(calcmultiwave, "multiwave", "multiwave", 1, 1, freq);
  p->skip= skipmultiwave;
  p->exit= exitmultiwave;
  p->private[0]= d= new(struct multiwave);
  for( count= 0; count< PLENTY; ++count )
    checkinputnr[count]= 0;
  naddinputs= 0;
  d->freqused= 0;
  for( count= 0; count< dimension; ++count ) {
    if( indextypes[2*count]!=MWTYPE_FREQ )
      ++naddinputs;
    else {
      if( d->freqused )
	fprintf(stderr, "Warning: multiwave: MWTYPE_FREQ should be used only once.\n");
      ++d->freqused;
    }
    if( indextypes[2*count] != MWTYPE_FREQ ) {
      if( indextypes[2*count+1] < dimension-d->freqused )
	++checkinputnr[indextypes[2*count+1]];
      else
	fprintf(stderr, "Warning: multiwave: input number given in indextypes array element %d exceeds dimension (taking into account MWTYPE_FREQ).  It will be set to 0.\n", count );
    }
  }
  for( count= 0; count< dimension-d->freqused; ++count )
    if( !checkinputnr[count] )
      fprintf(stderr, "Warning: multiwave: coordinate %d (excluding freq) not used.\n", count );
    else if( checkinputnr[count]> 1 )
      fprintf(stderr, "Warning: multiwave: coordinate %d (excluding freq) used more than once.\n", count );
  for( ; count< dimension; ++count )
    if( checkinputnr[count] ) {
      fprintf(stderr, "Warning: multiwave: some input numbers are too large and will be set to 0.\n" );
      break;
    }
  va_start(varargs, freq);
  for( count= 0; count< naddinputs; ++count )
    addinput(p, va_arg(varargs, sndobj *));
  va_end(varargs);
  d->fpos= 0.0;
  d->dim= dimension;
  for( d->nfiles= 0; ascfiles[d->nfiles] && *ascfiles[d->nfiles]; ++d->nfiles );
  if( d->nfiles> MULTIWAVE_MAXFILES ) {
    fprintf(stderr, "Warning: multiwave: Too many files given (%d, max= %d). Change limit in source if you need more. Returning 0.\n", d->nfiles, MULTIWAVE_MAXFILES );
    return c(0);
  }
  p->private[1]= malloc( (sizeof(float*)+sizeof(long)+sizeof(float)) * d->nfiles );
  d->wave= (float **)p->private[1];
  d->size= (long *)(d->wave+d->nfiles);
  d->cache= (float*)(d->size+d->nfiles);
  d->nsuccess= 0;
  for( fileind= 0; fileind < d->nfiles; ++fileind )
  {
    file= fopen( ascfiles[fileind], "r" );
    if( !file ) {
      fprintf(stderr, "Warning: multiwave: Could not open file `%s'. Ignoring this file.\n", ascfiles[fileind] );
      d->wave[fileind]= NULL;
      d->size[fileind]= 0;
      continue;
    }
    sharedname= news(char, NAMESTRLEN);
    strncpy(sharedname, ascfiles[fileind], NAMESTRLEN );
    parsedfile= parseascfile(file, sharedname);
    if( !parsedfile ) {
      fprintf(stderr, "Warning: multiwave: Error reading/parsing file `%s'. Ignoring this file.\n", ascfiles[fileind] );
      d->wave[fileind]= NULL;
      d->size[fileind]= 0;
      continue;
    }
    fclose(file);
    d->wave[fileind]= (float *)((long*)parsedfile + 1);
    d->size[fileind]= *(long*)parsedfile;
    ++d->nsuccess;
  }
  if( !d->nsuccess ) {
    fprintf(stderr, "Warning: multiwave: No files successfully read. Returning 0.\n" );
    return c(0);
  }
//fprintf(stderr, "d->nfiles= %d, d->success= %d\n", d->nfiles, d->nsuccess );
  if( d->freqused> 1 )
    d->freqused= 1;
  write= d->tree= news( struct multi_node, d->nfiles*d->dim + 1 );
  d->rootchildren= buildsubtree( &write, d, indextypes, positions, 0 );
  return p;
}


int buildsubtree( struct multi_node **write, struct multiwave *d, const int *indextypes, const float *positions, int depth )
{
  static float indexvalues[PLENTY];
  struct multi_node *localwrite, *reread;
  int prevdepth, fileind, nsiblings, sibling;
  
//printf( "buildsubtree: depth %d.\n", depth );
  localwrite= *write;
  for( fileind= 0; fileind< d->nfiles; ++fileind ) {
    // Ignore files that could not be read:
    if( !d->wave[fileind] ) {
      continue;
    }
    // Ignore files which are part of a different subtree because coordinates
    // higher in the hierarchy differ:
    for( prevdepth= 0; prevdepth< depth; ++prevdepth )
      if( positions[fileind*d->dim+prevdepth] != indexvalues[prevdepth] ) {
	break;
      }
    if( prevdepth< depth )
      continue;
    // Ignore files which have (up to this depth) the same coordinates as ones 
    // we already had:
    for( reread= *write; reread< localwrite; ++reread )
      if( reread->index == positions[fileind*d->dim+depth] ) {
	break;
      }
    // New value for the right coordinate => new node entry
    if( reread == localwrite ) {
      localwrite->index= positions[fileind*d->dim+depth];
//printf( "adding %dth coordinate %g\n", depth,  localwrite->index );
      if( indextypes[2*depth]==MWTYPE_FREQ ) {
	localwrite->type= MWTYPE_LINEAR;
	localwrite->inputnr= 0;
      }
      else {
	localwrite->type= indextypes[2*depth] & MWTYPE_MASK;
	localwrite->inputnr= indextypes[2*depth+1] + 1;
      }
      localwrite->flags= indextypes[2*depth] & MWFLAGS_MASK;
      if( localwrite->type==MWTYPE_WRAP ) {
	while( localwrite->index >= 1.0 )
	  localwrite->index -= 2.0;
	while( localwrite->index < -1.0 )
	  localwrite->index += 2.0;
      }
      localwrite->wavenr= fileind;
      localwrite->ngrandchildren= 0;
      localwrite->firstchild= NULL;
      ++localwrite;
    }
  }
  nsiblings= localwrite - *write;
  *write= localwrite;
  if( !nsiblings ) {
    fprintf(stderr, "Fatal internal error in buildsubtree (called by multiwave): no siblings found! Depth= %d, position %d in tree array.\n", depth, (int)(localwrite - d->tree) );
    exit(-1);
  }
  localwrite -= nsiblings;
  qsort( localwrite, nsiblings, sizeof(struct multi_node), cmpnodeindex );
  if( depth+1 < d->dim )
    for( sibling= 0; sibling< nsiblings; ++sibling ) {
      localwrite->firstchild= *write;
      indexvalues[depth]= localwrite->index;
      localwrite->ngrandchildren= buildsubtree( write, d, indextypes, positions, depth+1 );
      ++localwrite;
    }
  return nsiblings;
}

int cmpnodeindex( const void *p1, const void *p2 )
{
  const struct multi_node *n1= p1, *n2= p2;
  
  if( n1->index < n2->index )
    return -1;
  else if( n1->index > n2->index )
    return 1;
  else
    return 0;
}


float calcmultiwave(sndobj *me, sndobj *caller, int innr)
{
  struct multiwave *d;
  int inputnr;
  
  d= (struct multiwave *)me->private[0];
  for( inputnr= 0; inputnr< me->nin; ++inputnr )
    INCALC(inputnr);
  evaluatemultiwave( me, -1 );  //  clear cache
  OUTPUT(0, evaluatesubtree(me, d->tree, d->rootchildren, d->dim ) );
  d->fpos += INPUT(0)/(float)SAMPRATE;
  while( d->fpos>= 1.0 )
    d->fpos -= 1.0;
  CALCRETURN;
}

/*
Recursively evaluates the multiwave tree structure.  At the lowest level, 
evaluatemultiwave is called, which uses lininterpol to determine the value of
individual wave forms.  level is not the depth used for constructing the tree,
but rather dimension-depth.  This saves one argument, the dimension.
*/
float evaluatesubtree( sndobj *me, struct multi_node *subtree, int children, int level )
{
  float index, leftval, rightval;
  int sibling, leftind, rightind, wraparound;
  
  if( children == 1 ) {
    if( level == 1 )
      return evaluatemultiwave( me, subtree->wavenr );
    else
      return evaluatesubtree( me, subtree->firstchild, subtree->ngrandchildren, level-1 );
  }
  index= INPUT( subtree->inputnr );
  if( subtree->type==MWTYPE_WRAP ) {
    while( index >= 1.0 )
      index -= 2.0;
    while( index < -1.0 )
      index += 2.0;
  }
  if( (index<= subtree->index || index>= subtree[children-1].index) &&
    (subtree->type==MWTYPE_WRAP || (subtree->flags&MWFLAG_CONTINUE)==0) ) {
    if( subtree->type==MWTYPE_WRAP ) {
      wraparound= 1;
      leftind= children-1;
      rightind= 0;
      if( index<= subtree->index )
	index += 2.0;
    }
    else {
      if( index<= subtree->index )
	leftind= 0;
      else
	leftind= children-1;
      if( level == 1 )
	return evaluatemultiwave( me, subtree[leftind].wavenr );
      else
	return evaluatesubtree( me, subtree[leftind].firstchild, subtree[leftind].ngrandchildren, level-1 );
    }
  }
  else {
    wraparound= 0;
    for( sibling= 1; sibling< children-1 && index> subtree[sibling].index; ++sibling );
    leftind= sibling-1;
    rightind= sibling;
  }
  if( level == 1 ) {
    leftval= evaluatemultiwave( me, subtree[leftind].wavenr );
    rightval= evaluatemultiwave( me, subtree[rightind].wavenr );
  }
  else {
    leftval= evaluatesubtree( me, subtree[leftind].firstchild, subtree[leftind].ngrandchildren, level-1 );
    rightval= evaluatesubtree( me, subtree[rightind].firstchild, subtree[rightind].ngrandchildren, level-1 );
  }
  if( wraparound )
    return leftval + (rightval-leftval)*(index-subtree[leftind].index)/
		(1.0-subtree[leftind].index+subtree[rightind].index-(-1.0));
  else
    return leftval + (rightval-leftval)*(index-subtree[leftind].index)/
			    (subtree[rightind].index-subtree[leftind].index);
}

float evaluatemultiwave( sndobj *me, int wavenr )
{
  struct multiwave *d;
  int count;
  
  d= (struct multiwave *)me->private[0];
  if( wavenr< 0 ) {
    for( count= 0; count< d->nfiles; ++count )
      d->cache[count]= MAGIC;
    return 0.0;
  }
  else if( wavenr>= d->nfiles )
    return 0.0;
  else if( d->cache[wavenr]!=MAGIC )
    return d->cache[wavenr];
  else {
//printf("evaluatemultiwave called for wave nr %d, pos %g\n", wavenr, d->fpos );
    d->cache[wavenr]= lininterpol( d->wave[wavenr], d->size[wavenr], d->fpos*d->size[wavenr] );
    return d->cache[wavenr];
  }
}


void skipmultiwave(sndobj *me, sndobj *caller, int innr)
{
  struct multiwave *d;
  int inputnr;
  
  d= (struct multiwave *)me->private[0];
  for( inputnr= 0; inputnr< me->nin; ++inputnr )
    INSKIP(inputnr);
  d->fpos += INPUT(0)/(float)SAMPRATE;
  while( d->fpos>= 1.0 )
    d->fpos -= 1.0;
  return;
}


void exitmultiwave(sndobj *me)
{
  struct multiwave *d;
  int count;
  
  d= (struct multiwave *)me->private[0];
  for( count= 0; count< d->nfiles; ++count )
    freeshared( (long*)d->wave[count] - 1 );
}


/*=============================================================================
    flatrand (max)

Evenly distributed random numbers between -(max) and (max).

See also: gaussrand, fbm
=============================================================================*/

float calcflatrand(sndobj *me, sndobj *, int);

sndobj *flatrand(float max)
{
  sndobj *p;
  
  p= newsndo(calcflatrand, "flatrand", "flatrand", 1, 0);
  p->private[0]= makerng();
  p->private[1]= new(float);
  *(float*)p->private[1]= max;
  return p;
}

float calcflatrand(sndobj *me, sndobj *caller, int innr)
{
  OUTPUT(0, *(float*)me->private[1] * 
	    flatrng( (rngdesc*)me->private[0] ) );
  CALCRETURN;
}


/*=============================================================================
    gaussrand (sigma)

Gaussian-distributed random numbers with mean deviation (sigma).

See also: flatrand, fbm
=============================================================================*/

float calcgaussrand(sndobj *me, sndobj *, int);

sndobj *gaussrand(float sigma)
{
  sndobj *p;
  
  p= newsndo(calcgaussrand, "gaussrand", "gaussrand", 1, 0);
  p->private[0]= makerng();
  p->private[1]= new(float);
  *(float*)p->private[1]= sigma;
  return p;
}

float calcgaussrand(sndobj *me, sndobj *caller, int innr)
{
  OUTPUT(0, *(float*)me->private[1] * 
	    gaussrng( (rngdesc*)me->private[0] ) );
  CALCRETURN;
}


/*=============================================================================
    fbm <lacunarity> <H> <cutoff> <reset>

This is a fairly realistic noise generator based on a fractional brownian
motion algorithm by Dietmar Saupe (or at least presented by him in "The 
Science of Fractal Images",
Springer 1988).  It works by successively interpolating linearly and adding a
random contribution.  <lacunarity> is the factor by which the grid spacings of
successively added structures differ.  It must be between 0 and 0.9.  A larger
value gives a fuller sound.  <H> determines the fractal dimension, and hence
the frequency spectrum, of the output.  Its range is from 0 to somewhat below
1, so that sqrt(1.0 - pow(lacunarity, 2-2*H)) is still > 0 numerically.  A
reasonable limit is 0.9999 for lac.=0.5 on my system and gets smaller for
rising lacunarity.  Small H values give nearly white noise, large values favour
low frequencies.  <cutoff> is a highpass cutoff frequency.  This is necessary
for large H, to stop the fractal random walk from wandering too far.  It must
be larger than 1 Hz.  The maximum lacunarity and the minimal cutoff can be
changed in the source code.  <reset> causes the generator to reset if non-zero
and may be passed as NULL.  Using this input is recommended if and when any of
the other inputs change sharply (e.g. as a step function).  However, resetting
fbm results in a small click.  Examples:

fbm 0.5  0.4      1200    -- air compressor 
fbm 0.5  0.8       150    -- rumbling of a waterfall
fbm 0.08 0.53      500    -- nasty rain heard from indoors
fbm 0.8  0->0.75   150    -- gravel sliding from the back of a truck

See also: flatrand, gaussrand, walk, slide, linnoise
=============================================================================*/

struct fbm {
  float     *currval, *prevval, *sigmatotsqr;
  double    *gridpoint, *gridspacing;
  rngdesc   *rng;
  long      currmaxdepth, count;
  int       firstcall;
};

float calcfbm(sndobj *me, sndobj *, int);

#define FBM_MINCUTOFF   1
#define FBM_MAXLACU     0.9

sndobj *fbm(sndobj *lacunarity, sndobj *H, sndobj *cutoff, sndobj *reset)
{
  sndobj *p;
  struct fbm *d;
  long maxdepthp1;

  p= newsndo( calcfbm, "fbm", "fbm", 1, 3, lacunarity, H, cutoff );
  if( reset )
    addinput(p, reset);
  //  Note: I don't define a skip function for fbm.  As a consequence, if the
  //  inputs change while fbm is being skipped, it will take time to adjust
  //  after the skip is over - a different result than without the skip.  This
  //  is probably tolerable.
  p->private[0]= d= new(struct fbm);
  maxdepthp1= (long) ceil(- log(0.5*(double)SAMPRATE/(double)FBM_MINCUTOFF)/log((double)FBM_MAXLACU) ) + 1L;
  p->private[1]= d->currval= news(float, maxdepthp1);
  p->private[2]= d->prevval= news(float, maxdepthp1);
  p->private[3]= d->sigmatotsqr= news(float, maxdepthp1);
  p->private[4]= d->gridpoint= news(double, maxdepthp1);
  p->private[5]= d->gridspacing= news(double, maxdepthp1);
  p->private[6]= d->rng= makerng();
  d->firstcall= 1;
  return p;
}

float calcfbm(sndobj *me, sndobj *caller, int innr)
{
  struct fbm *d;
  double lacunarity, H, cutoff;
  float sigma;
  long depth;

  d= (struct fbm *)me->private[0];
  lacunarity= INCALC(0);
  H= INCALC(1);
  cutoff= INCALC(2);
  if( lacunarity < 0.0 )
    lacunarity= 0.0;
  else if( lacunarity > FBM_MAXLACU )
    lacunarity= FBM_MAXLACU;
  if( H < 0.0 )
    H= 0.0;
  else if( H > 1.0 )
    H= 1.0;
  if( cutoff< FBM_MINCUTOFF )
    cutoff= FBM_MINCUTOFF;

  if( (me->nin>3 && INCALC(3)) || d->firstcall )    /* ! left to right ! */
  {
    depth= 0;
    sigma= M_SQRT1_2 * pow(lacunarity, H) * sqrt( 1.0 - pow(lacunarity, 2.0-2.0*H) );
    d->currval[0]= sigma*gaussrng(d->rng);
    d->gridpoint[0]= 0.0;
    d->currmaxdepth= 0;
    d->count= 0;
    d->firstcall= 0;
  }
  else {
    if( d->count >= 100000L ) {
      for( depth= 0; depth <= d->currmaxdepth; ++depth )
	d->gridpoint[depth] -= d->count;
      d->count= 0;
    }
    depth= d->currmaxdepth;
  }

  while( depth> 0 && d->gridpoint[depth] + d->gridspacing[depth-1]*lacunarity
			     > d->gridpoint[depth-1] )
    --depth;

  if( depth==0 )
  {
    d->gridspacing[0]= (double)SAMPRATE/(2.0*cutoff);
    d->gridpoint[0] += d->gridspacing[0];
    d->prevval[0]= d->currval[0];
    sigma= M_SQRT1_2 * pow(lacunarity, H) * sqrt( 1.0 - pow(lacunarity, 2.0-2.0*H) );
    d->currval[0]= sigma*gaussrng(d->rng);
    d->sigmatotsqr[0]= sigma*sigma;
    ++depth;
  }
  
  for( ; d->gridspacing[depth-1] > 1.0; ++depth )
  {
    d->gridspacing[depth]= d->gridspacing[depth-1] * lacunarity;
    sigma= M_SQRT1_2 * pow(lacunarity, (double)(depth+1) * H) * 
		    sqrt( 1.0 - pow(lacunarity, 2.0-2.0*H));
    if( d->gridspacing[depth] < 1.0 ) {
      sigma *= (d->gridspacing[depth]-lacunarity)/(1.0-lacunarity);
      d->gridspacing[depth]= 1.0;
    }
    d->sigmatotsqr[depth]= d->sigmatotsqr[depth-1] + sigma*sigma;
    if( depth > d->currmaxdepth )
      d->gridpoint[depth]= (double)d->count;
    else
      d->gridpoint[depth] += d->gridspacing[depth];
    d->prevval[depth]= d->currval[depth];
    if( d->gridpoint[depth]>=d->gridpoint[depth-1] )
      d->currval[depth]= d->currval[depth-1] + sigma*gaussrng(d->rng);
    else if( d->gridpoint[depth]<=d->gridpoint[depth-1]-d->gridspacing[depth-1] )
      d->currval[depth]= d->prevval[depth-1] + sigma*gaussrng(d->rng);
    else
      d->currval[depth]= sigma*gaussrng(d->rng) + d->prevval[depth-1] + 
      (d->gridpoint[depth] - (d->gridpoint[depth-1]-d->gridspacing[depth-1])) *
      (d->currval[depth-1] - d->prevval[depth-1]) / d->gridspacing[depth-1];
  }
  d->currmaxdepth= depth-1;

  /*  divide by 3.5 standard deviations to normalise:  */
  OUTPUT(0, d->currval[d->currmaxdepth]/(3.5*sqrt(d->sigmatotsqr[d->currmaxdepth])));
  ++d->count;
  if( !d->sigmatotsqr[d->currmaxdepth] )
    printf ("zero sigma_tot at %ld\n", d->count );
  CALCRETURN;
}


/*=============================================================================
  linnoise <slope> <nastyness>

This is a simple noise generator creating a piecewise linear waveform.  <slope>
controls the frequency, <nastyness> (typically between 0 and 1) the sharpness
of the noise.  For instance, linnoise 300 1 generates "exploding bomb" noise.
Applying a low pass to the output gives further possibilities.

See also: flatrand, gaussrand, fbm
=============================================================================*/

struct linnoise {
  double    inc, prevval;
  int       samcount;
  rngdesc   *rngstate;
};

float calclinnoise(sndobj *me, sndobj *, int);
void skiplinnoise(sndobj *me, sndobj *, int);

sndobj *linnoise( sndobj *slope, sndobj *nastyness )
{
  sndobj *p;
  struct linnoise *d;

  p= newsndo( calclinnoise, "linnoise", "linnoise", 1, 2, slope, nastyness );
  p->skip= skiplinnoise;
  p->private[0]= d= new(struct linnoise);
  d->prevval= 0.0;
  d->samcount=0;
  d->rngstate= makerng();
  return p;
}

float calclinnoise(sndobj *me, sndobj *caller, int innr)
{
  struct linnoise *d;
  double nextstep;
  float slope, nastyness;

  d= (struct linnoise *)me->private[0];
  slope= INCALC(0);
  nastyness= INCALC(1);
  if( d->samcount<=0 )
  {
    if( slope < 1.0/HORIZON )
      d->inc= (double)HORIZON/SAMPRATE;
    else if( slope > SAMPRATE )
      d->inc= 1.0;
    else
      d->inc= slope/SAMPRATE;
    nextstep= 0.5*nastyness*gaussrng(d->rngstate);
    if( d->prevval + d->inc > 1.0 )
      nextstep= -fabs(nextstep);
    else if( d->prevval - d->inc < -1.0 )
      nextstep= fabs(nextstep);
    if( d->prevval + nextstep > 1.0 )
      nextstep= 1.0-d->prevval;
    else if( d->prevval + nextstep < -1.0 )
      nextstep= -1.0+-d->prevval;
    if( nextstep<0.0 )
      d->inc= -d->inc;
    d->samcount= (int)trunc(nextstep/d->inc);
    if( d->samcount<=0 )
      d->samcount= 1;
  }
  d->prevval += d->inc;
  --d->samcount;
  OUTPUT(0, (float)d->prevval);
  CALCRETURN;
}

void skiplinnoise(sndobj *me, sndobj *caller, int innr)
{
  struct linnoise *d;

  d= (struct linnoise *)me->private[0];
  d->samcount= 0;
  INSKIP(0);
  INSKIP(1);
}


/*=============================================================================
  walk (steprms)

Gaussian random walk wrapping around at 1 and -1.  The result is in [-1,1).
(steprms) is the RMS size of a step from one sample to the next (useful range
0.01 .. 0.1).

See also: fwalk, fbm
=============================================================================*/

struct walk {
  float     rms, outval;
  rngdesc   *rand;
};

float calcwalk(sndobj *me, sndobj *, int);

sndobj *walk(float steprms)
{
  sndobj *p;
  struct walk *d;
  
  p= newsndo(calcwalk, "walk", "walk", 1, 0);
  p->private[0]= d= new(struct walk);
  d->rms= steprms;
  d->outval= 0.0;
  p->private[1]= d->rand= makerng();
  return p;
}

float calcwalk(sndobj *me, sndobj *caller, int innr)
{
  struct walk *d;
  
  d= (struct walk *)me->private[0];
  d->outval += d->rms*gaussrng(d->rand);
  if( d->outval >= 1.0 )
    d->outval -= 2.0;
  else if( d->outval < -1.0 )
    d->outval += 2.0;
  OUTPUT(0, d->outval);
  CALCRETURN;
}


/*=============================================================================
  fwalk <speed>

Gaussian random walk wrapping around at 1 and -1.  The result is in [-1,1).
The input determines what distance the walk covers per second.

See also: walk, fbm
=============================================================================*/

float calcfwalk(sndobj *me, sndobj *, int);

sndobj *fwalk(sndobj *speed)
{
  sndobj *p;
  
  p= newsndo(calcfwalk, "fwalk", "fwalk", 1, 1, speed);
  p->private[0]= makerng();
  p->ch[0]= 0.0;
  return p;
}

float calcfwalk(sndobj *me, sndobj *caller, int innr)
{
  rngdesc *randgen;
  
  randgen= (rngdesc *)me->private[0];
  me->ch[0] += INCALC(0)/SAMPRATE*gaussrng(randgen);
  while( me->ch[0] >= 1.0 )
    me->ch[0] -= 2.0;
  while( me->ch[0] < -1.0 )
    me->ch[0] += 2.0;
  CALCRETURN;
}


/*=============================================================================
    slide <minfreq> <maxfreq> <lacunarity> <exponent> <meandev>

Generates more or less regular positive pulses similar to the lurches of a
block just sliding on a vibrating sloped surface.  <minfreq> and <maxfreq> give
the range of frequencies for series of pulses; |<lacunarity>| is the factor
between the frequencies (or intervals) of different series.  The minimum and
maximum frequencies are limited to SLIDE_MINMINF (1 Hz) and SLIDE_MAXMAXF
(SAMPRATE/2); lacunarity is limited to SLIDE_MAXLACU (0.9).  These constants
can be changed in the source code.  The height of the pulses is related with
the average interval between them by a power law with <exponent>: height is
proportional to interval^<exponent>.  The height of the pulses with frequency
<minfreq> is always 1.  So for positive <exponent>, smaller pulses are more
frequent.  Both the interval and the height are perturbed by a factor (1 +
random), where the random contribution is a gaussian-distributed variable with
mean deviation <meandev>.  If two pulses of different series coincide, they are
not added up; rather, only the higher pulse is output.

Changes in the inputs take effect only when events happen, so an average delay
of 1/(previous minimum frequency) is to be expected.

This object is as versatile as fbm.  It can sound like a violin's bow, or like
crackling fire or other noises.

See also: fbm, interleave
=============================================================================*/

struct slide {
  rngdesc *rng;
  float	  currminf, currmaxf;
  double  maxinter, lacunarity, exponent, meandev;
  int     *countdown;
  int     maxlevels, levels, initialised;
  int	  updatecountdown;
};

float calcslide(sndobj *me, sndobj *caller, int innr);

#define SLIDE_MAXLACU	0.9
#define SLIDE_MINMINF	1.0
#define SLIDE_MAXMAXF	((double)(SAMPRATE/2))
#define SLIDE_MAXUPDATEINTER	0.001		// in seconds

sndobj *slide( sndobj *minfreq, sndobj *maxfreq, sndobj *lacunarity, sndobj *exponent, sndobj *meandev )
{
  sndobj *p;
  struct slide *d;

  p= newsndo(calcslide, "slide", "slide", 1, 5, lacunarity, exponent, meandev, minfreq, maxfreq);
  p->skip= dontskip;
  p->private[0]= d= new(struct slide);
  p->private[1]= d->rng= makerng();
  d->maxlevels= (int)floor(-log(SLIDE_MAXMAXF/SLIDE_MINMINF)/log(SLIDE_MAXLACU)) + 1;
  p->private[2]= d->countdown= news(int, d->maxlevels);
  d->levels= 0;
  d->initialised= 0;
  return p;
}

void updateslideparms(sndobj *me, struct slide *d, int eventlevel);

float calcslide(sndobj *me, sndobj *caller, int innr)
{
  struct slide *d;
  int count, eventlevel, updated;

  d= (struct slide *)me->private[0];
  if( !d->initialised ) {
    d->initialised= 1;
    updateslideparms(me, d, -1);
    for( count= 0, eventlevel= -1; count < d->levels; --count )
      if( d->countdown[count] <= 0 ) {
	eventlevel= count;
	break;
      }
    if( eventlevel < 0 )
      OUTPUT(0, 0.0);
    else {
      double shrinkfact, interval, ampl;

      shrinkfact= pow(d->lacunarity, eventlevel);
      interval= d->maxinter * shrinkfact;
      ampl= pow(shrinkfact, d->exponent);
      OUTPUT(0, ampl * (1.0 + d->meandev*gaussrng(d->rng)) );
      d->countdown[eventlevel]= (int)floor(interval * (1.0 + d->meandev*gaussrng(d->rng)) + 0.5);
      for( count= eventlevel+1, interval *= d->lacunarity; count< d->levels; 
	  			++count, interval *= d->lacunarity )
	if( d->countdown[count] <= 0 )
	  d->countdown[count]= (int)floor(interval * (1.0 + d->meandev*gaussrng(d->rng)) + 0.5);
      
    }
    CALCRETURN;
  }
  if( --d->updatecountdown <= 0 ) {
    updateslideparms(me, d, -1);
    updated= 1;
  }
  else
    updated= 0;
  for( count= d->levels-1, eventlevel= -1; count >= 0; --count )
    if( --d->countdown[count] <= 0 )
      eventlevel= count;
  if( d->levels>=0 && eventlevel < 0 ) {
    OUTPUT(0, 0.0);
    INSKIP(0);
    INSKIP(1);
    INSKIP(2);
    INSKIP(3);
    INSKIP(4);
  }
  else if( !updated )
    updateslideparms(me, d, eventlevel);
  CALCRETURN;
}

// Calls all inputs' calc function and updates the parameters of this sndobj
// accordingly.  If eventlevel is >= 0, outputs the pulse corresponding to the
// event level and the new parameters; otherwise 0 is output.  The event
// counters are updated by setting those which are due anew and rescaling the
// others' countdowns according to the difference in their old and new
// intervals.
void updateslideparms(sndobj *me, struct slide *d, int eventlevel)
{
  double minf, maxf, lac, expo, mdev, range, shrinkfact, ampl, newmaxint, interval, oldinter;
  int count, newlevels, minlevels, parchange;

  d->updatecountdown= SLIDE_MAXUPDATEINTER*SAMPRATE;
  minf= INCALC(3);
  if( minf < SLIDE_MINMINF )
    minf= SLIDE_MINMINF;
  maxf= INCALC(4);
  if( maxf > SLIDE_MAXMAXF )
    maxf= SLIDE_MAXMAXF;
  if( maxf < minf ) {   // no output
    d->levels= 0;
    d->currmaxf= maxf;
    d->currminf= minf;
    INSKIP(0);
    INSKIP(1);
    INSKIP(2);
    OUTPUT(0, 0.0);
    return;
  }
  // read remaining parameters and update event counters
  lac= fabsf(INCALC(0));
  if( lac > SLIDE_MAXLACU )
    lac= SLIDE_MAXLACU;
  expo= INCALC(1);
  mdev= INCALC(2);
  newmaxint= (double)SAMPRATE/minf;
  range= log(maxf/minf);
  shrinkfact= pow(lac, eventlevel);
  if( eventlevel >= 0 ) {
    interval= newmaxint * shrinkfact;
    ampl= pow(shrinkfact, expo);
    OUTPUT(0, ampl * (1.0 + mdev*gaussrng(d->rng)) );
    d->countdown[eventlevel]= (int)floor(interval * (1.0 + mdev*gaussrng(d->rng)) + 0.5);
  }
  else
    OUTPUT(0, 0.0);
  if( lac == 0.0 )	 newlevels= 1;
  else			 newlevels= (int)floor(-range/log(lac)) + 1;
  minlevels= newlevels < d->levels? newlevels : d->levels;
  parchange= d->currmaxf!=maxf || d->currminf!=minf || d->lacunarity!=lac ||
    		d->exponent!=expo || d->meandev!=mdev;
  for( count= 0, interval= newmaxint, oldinter= d->maxinter; count< minlevels; 
      ++count, interval *= lac, oldinter *= d->lacunarity )
    if( d->countdown[count] <= 0 )
      d->countdown[count]= (int)floor(interval * (1.0 + mdev*gaussrng(d->rng)) + 0.5);
    else if( parchange )
      d->countdown[count]= (int)floor(interval/oldinter * d->countdown[count] + 0.5);
  for( ; count < newlevels; ++count, interval *= lac )
    d->countdown[count]= (int)floor(interval * 0.5 * (1.0 + flatrng(d->rng))
				      * (1.0 + mdev*gaussrng(d->rng)) + 0.5);
  d->levels= newlevels;
  d->maxinter= newmaxint;
  d->currmaxf= maxf;
  d->currminf= minf;
  d->lacunarity= lac;
  d->exponent= expo;
  d->meandev= mdev;
}


/*=============================================================================
    mike (linein) (level)

Microphone or line input through open sound system.  (linein) is 0 for
microphone or !=0 for line input.  (level) determines the recording volume
(0..1; <0: leave unchanged).

Unavailable if NO_OSS has been defined in sndsys.h.

See also: oss, alsa
=============================================================================*/

#ifndef NO_OSS

struct mike {
  int           dsp, bufind, bufsize;
  signed short  *buf;
};

float calcmike(sndobj *me, sndobj *, int);

#define MIKE_BUFSIZE    1024

sndobj *mike( int linein, float level )
{
  sndobj *p;
  struct mike *d;
  int ioarg, indevice;
  
  p= newsndo( calcmike, "mike", "mike", 2, 0 );
  p->skip= dontskip;
  p->private[0]= d= new(struct mike);
  d->bufsize= 2*MIKE_BUFSIZE;
  p->private[1]= d->buf= news(signed short, d->bufsize);
  d->bufind= d->bufsize;
  d->dsp= dspdesc(0);
  if( d->dsp < 0 ) {
    fprintf(stderr, "Warning: mike: DSP could not be opened.  Returning 0.\n");
    return c(0);
  }
  if( !linein )
    indevice= SOUND_MIXER_MIC;
  else
    indevice= SOUND_MIXER_LINE;
  ioarg= 1 << indevice;
  ioctl( d->dsp, SOUND_MIXER_WRITE_RECSRC, &ioarg );
  if( level >= 0.0 ) {
    ioarg= (int)round(level*100);
    if( ioarg > 100 )
      ioarg= 100;
    ioarg *= 257;
    ioctl( d->dsp, MIXER_WRITE(indevice), &ioarg );
    ioctl( d->dsp, SOUND_MIXER_WRITE_IGAIN, &ioarg );
  }
  return p;
}

float calcmike(sndobj *me, sndobj *caller, int innr)
{
  struct mike *d;
  
  d= (struct mike *)me->private[0];
  if( d->bufind >= d->bufsize ) {
    read( d->dsp, d->buf, sizeof(short)*d->bufsize );
    d->bufind= 0;
  }
  OUTPUT(0, (((float)d->buf[d->bufind++])+32768.0)/65535.0*2.0 - 1.0 );
  OUTPUT(1, (((float)d->buf[d->bufind++])+32768.0)/65535.0*2.0 - 1.0 );
  CALCRETURN;
}

#endif



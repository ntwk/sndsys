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
#include <float.h>
#include <sys/ioctl.h>
#include <unistd.h>
#ifndef NO_OSS
#include <linux/soundcard.h>
#endif
#ifndef NO_ALSA
#include <alsa/asoundlib.h>
#include <sys/select.h>
#include <sys/time.h>
#endif

#include "sndsys.h"


signed short toshort( float fval );

signed short toshort( float fval )
{
  if( isnan(fval) )
    return 0;
  if( fval > 1.0 )
    return 32767;
  else if( fval < -1.0 )
    return -32768;
  else
    return (signed short)floor(0.5 + (fval+1.0)*65535.0/2.0 - 32768.0);
}


/*=============================================================================
    writeout "filename" <signal>

Write to .wav file as 16 bit data.

See also: file
=============================================================================*/

struct writeout {
  const char    *name;
  FILE          *handle;
  wavheader     *head;
  long          length;
  float         min, max;
};

float calcwriteout(sndobj *me, sndobj *, int);
void exitwriteout(sndobj *me);

wavheader writeouttemplhead= { 'FFIR', 0xFFFFFFFFL, 'EVAW',
' tmf', 16L, 1, 1, (long)SAMPRATE, (long)SAMPRATE*2, 2, 16, 'atad', 0xFFFFFFFFL };
//{ 0x46464952L, 0L, 0x45564157L, 0x20746d66L,
//16L, 1, 1, (long)SAMPRATE, (long)SAMPRATE*2, 2, 16, 0x61746164, 0L };

sndobj *writeout(const char *filename, sndobj *in)
{
  sndobj *p;
  struct writeout *d;
  
  p= newsndo(calcwriteout, filename, "writeout", in->nch, 1, in);
  p->skip= dontskip;
  p->private[0]= d= (struct writeout *)malloc(sizeof(struct writeout));
  p->private[1]= d->head= (wavheader*)malloc(sizeof(wavheader));
  d->name= filename;
  d->length= 0L;
  d->min= FLT_MAX;
  d->max= -FLT_MAX;
  d->handle= fopen(filename, "w");
  if( !d->handle )
    fprintf( stderr, "Warning: Error opening file `%s' for writing.\n", filename );
  p->exit= exitwriteout;
  memcpy( d->head, &writeouttemplhead, sizeof(wavheader) );
  d->head->nchannels *= in->nch;
  d->head->blockalign *= in->nch;
  d->head->byterate *= in->nch;
  fwrite( d->head, sizeof(wavheader), 1L, d->handle );
  return p;
}

float calcwriteout(sndobj *me, sndobj *caller, int innr)
{
  struct writeout *d;
  signed short values[PLENTY];
  float fval;
  int ch;
  
  d= (struct writeout *)me->private[0];
  INCALC(0);
  for( ch= 0; ch < me->in[0]->nch; ++ch ) {
    fval= INPUTC(0, ch);
    OUTPUT(ch, fval);
    values[ch]= toshort( fval );
    if( isnan(fval) )
      continue;
    if( fval>d->max )
      d->max= fval;
    if( fval<d->min )
      d->min= fval;
  }
  fwrite( values, (long)d->head->blockalign, 1L, d->handle );
  ++d->length;
  CALCRETURN;
}

void exitwriteout(sndobj *me)
{
  struct writeout *d;

  d= (struct writeout *)me->private[0];
  d->head->restsize= d->length*d->head->blockalign;
  d->head->totalsizem8= d->head->restsize+36L;
  fseek( d->handle, 0L, SEEK_SET );
  fwrite( d->head, sizeof(wavheader), 1L, d->handle );
  fclose( d->handle );
  printf( "Writeout to `%s' finished.  Min= %g, max= %g.\n", d->name, 
	    (double)d->min, (double)d->max );
}


/*=============================================================================
    print *stream* <signal>

Output to stream, as text.  Channels are separated by two tab characters.

See also: tprint, wprint
=============================================================================*/

float calcprint(sndobj *me, sndobj *, int);

sndobj *print( FILE *stream, sndobj *in )
{
  sndobj *p;
  char *name;
  
  if( stream==stdout )
    name= "stdout";
  else if( stream==stderr )
    name= "stderr";
  else
    name= "unknown stream";

  p= newsndo( calcprint, name, "print", in->nch, 1, in );
  p->skip= dontskip;
  p->private[0]= new(FILE*);
  *(FILE**)p->private[0]= stream;
  return p;
}


float calcprint(sndobj *me, sndobj *caller, int innr)
{
  FILE *stream;
  float val;
  int ch;
  
  stream= *(FILE**)me->private[0];
  if( me->nch )
  {
    val= INCALC(0);
    OUTPUT(0, val);
    fprintf( stream, "%g", (double)val );
    for( ch= 1; ch< me->nch; ++ch ) {
      val= INPUTC(0, ch);
      OUTPUT(ch, val);
      fprintf( stream, "\t\t%g", val );
    }
  }
  fprintf(stream, "\n");
  CALCRETURN;
}


/*=============================================================================
    tprint *stream* <signal>

Output to stream, as text.  The first column of the output is the tim ein
seconds.  It is followed by all the channels of <signal>.  Columns are
separated by two tab characters.

See also: print
=============================================================================*/

struct tprint {
  FILE *outstream;
  double time;
};

float calctprint(sndobj *me, sndobj *, int);

sndobj *tprint( FILE *stream, sndobj *in )
{
  sndobj *p;
  struct tprint *d;
  char *name;
  
  if( stream==stdout )
    name= "stdout";
  else if( stream==stderr )
    name= "stderr";
  else
    name= "unknown stream";

  p= newsndo( calctprint, name, "tprint", in->nch, 1, in );
  p->skip= dontskip;
  p->private[0]= d= (struct tprint *)new(struct tprint);
  d->outstream= stream;
  d->time= starttime;
  return p;
}


float calctprint(sndobj *me, sndobj *caller, int innr)
{
  struct tprint *d;
  float val;
  int ch;
  
  d= (struct tprint *)me->private[0];
  fprintf( d->outstream, "%.10g", d->time );
  INCALC(0);
  for( ch= 0; ch< me->nch; ++ch ) {
    val= INPUTC(0, ch);
    OUTPUT(ch, val);
    fprintf( d->outstream, "\t\t%g", (double)val );
  }
  fprintf(d->outstream, "\n");
  d->time += 1.0/SAMPRATE;
  CALCRETURN;
}


/*=============================================================================
   sndstat (min) (max) (step) <signal>

Statistics for signal input.  Values between (min) and (max) are counted in
bins of size (step).  Contributions from all channels are collected.  
=============================================================================*/

struct sndstat {
  long      *bin;
  float     min, max, minval, maxval, mininrange, maxinrange;
  double    sum, sqrsum, step;
  long      total, inf, nan, nbins;
};

float calcsndstat(sndobj *me, sndobj *, int);
void exitsndstat(sndobj *me);

sndobj *sndstat( float min, float max, float binsize, sndobj *signal )
{
  sndobj *p;
  struct sndstat *d;
  
  p= newsndo(calcsndstat, signal->name, "sndstat", signal->nch, 1, signal );
  p->skip= dontskip;
  p->exit= exitsndstat;
  p->private[0]= d= new(struct sndstat);
  if( min>max ) {
    d->min= max;
    d->max= min;
  }
  else {
    d->min= min;
    d->max= max;
  }
  d->step= fabsf(binsize);
  d->nbins= (long)ceil((d->max-d->min)/d->step);
  if( !d->nbins )
    d->nbins= 1;
  else if( d->nbins<0 )
    d->nbins= -d->nbins;
  p->private[1]= d->bin= calloc(d->nbins, sizeof(long));
  d->minval= d->mininrange= FLT_MAX;
  d->maxval= d->maxinrange= -FLT_MAX;
  d->sum= 0.0;
  d->sqrsum= 0.0;
  d->total= 0;
  d->inf= d->nan= 0;
  return p;
}

float calcsndstat(sndobj *me, sndobj *caller, int innr)
{
  struct sndstat *d;
  float val;
  long index;
  int ch;
  
  d= (struct sndstat *)me->private[0];
  INCALC(0);
  for( ch= 0; ch < me->nch; ++ch )
  {
    val= INPUTC(0, ch);
    OUTPUT(ch, val);
    if( !finite(val) ) {
      if( isnan(val) )
	++d->nan;
      else
	++d->inf;
      continue;
    }
    d->sum += val;
    d->sqrsum += val*val;
    if( val> d->maxval )
      d->maxval= val;
    if( val< d->minval )
      d->minval= val;
    if( d->step!=0.0 && val >= d->min && val <= d->max )
    {
      if( val> d->maxinrange )
	d->maxinrange= val;
      if( val< d->mininrange )
	d->mininrange= val;
      index= (long)floor((val-d->min)/d->step);
      if( index>=d->nbins )
	--index;
      ++d->bin[index++];
    }
  }
  d->total += me->nch;
  CALCRETURN;
}

void exitsndstat(sndobj *me)
{
  struct sndstat *d= (struct sndstat *)me->private[0];
  float limit;
  long i, inrange;

  printf( "Statistics for object `%s', type `%s', %d channel(s):\n", me->in[0]->name, me->in[0]->type, me->nch );
  limit= d->min;
  inrange= 0;
  if( d->step!=0.0 )
  {
    for( i= 0; i< d->nbins-1; ++i ) {
      printf("  Bin from % g to % g:\t%6ld (%5.2f %%)\n", limit, limit+d->step, d->bin[i], 100.0*(double)d->bin[i]/(double)d->total );
      inrange += d->bin[i];
      limit += d->step;
    }
    printf("  Bin from % g to % g:\t%6ld (%5.2f %%)\n", limit, d->max, d->bin[i], 100.0*(double)d->bin[i]/(double)d->total );
    inrange += d->bin[i];
    printf( "  Remaining values: %ld (%.2f %%) outside range, %ld (%.2f %%) "
	  "infinite, %ld (%.2f %%) NaN.\n", d->total-inrange-d->nan-d->inf, 
	  100.0*(double)(d->total-inrange-d->nan-d->inf)/(double)d->total,
	  d->inf, 100.0*(double)d->inf/(double)d->total, d->nan, 
	  100.0*(double)d->nan/(double)d->total );
  }
  else
    printf("  Illegal values: %ld (%.2f %%) infinite, %ld (%.2f %%) NaN.\n", 
	d->inf, 100.0*(double)d->inf/(double)d->total, d->nan,
	          100.0*(double)d->nan/(double)d->total );
  printf( "  %ld values in total (", d->total );
  if( me->nch>1 )
    printf( "%ld per channel, ", d->total/me->nch );
  printf("%g seconds), ", (double)d->total/me->nch/(double)SAMPRATE);
  d->sum /= d->total;
  d->sqrsum /= d->total;
    printf( "average %g, rms %g, variance %g.\n", d->sum, sqrt(d->sqrsum), sqrt(d->sqrsum-d->sum*d->sum) );
  if( d->step!=0.0 )
  printf("  Extremal values within range: %g -- %g\n", d->mininrange, d->maxinrange );
  printf("  Extremal values overall: %g -- %g\n", d->minval, d->maxval );
}


/*=============================================================================
    beancounter <signal>

This is a debug object which counts how often it has been calc()ed and how
often it has been skip()ped.  At program termination, the counts are printed 
out.

See also: skipwatch
=============================================================================*/

struct beancounter {
  int calccount, skipcount;
};

float calcbeancounter(sndobj *me, sndobj *, int);
void skipbeancounter(sndobj *me, sndobj *, int);
void exitbeancounter(sndobj *me);

sndobj *beancounter( sndobj *signal )
{
  sndobj *p;
  struct beancounter *d;

  p= newsndo(calcbeancounter, signal->name, "beancounter", signal->nch, 1, signal);
  p->skip= skipbeancounter;
  p->exit= exitbeancounter;
  p->private[0]= d= new(struct beancounter);
  d->calccount= d->skipcount= 0;
  return p;
}

float calcbeancounter(sndobj *me, sndobj *caller, int innr)
{
  struct beancounter *d;
  int ch;

  d= (struct beancounter *)me->private[0];
  ++d->calccount;
  OUTPUT(0, INCALC(0));
  for( ch= 1; ch< me->nch; ++ch )
    OUTPUT(ch, INPUTC(0, ch));
  CALCRETURN;
}

void skipbeancounter(sndobj *me, sndobj *caller, int innr)
{
  struct beancounter *d;

  d= (struct beancounter *)me->private[0];
  ++d->skipcount;
  INSKIP(0);
}

void exitbeancounter(sndobj *me)
{
  struct beancounter *d;
  double calcsecs, skipsecs;

  d= (struct beancounter *)me->private[0];
  calcsecs= (double)d->calccount/SAMPRATE;
  skipsecs= (double)d->skipcount/SAMPRATE;
  printf("Beancounter `%s': Processed %d samples (%g seconds), of which %d "
		"(%g seconds) were calculated and %d (%g seconds) skipped.\n", 
		me->name, d->calccount+d->skipcount, calcsecs+skipsecs,
		d->calccount, calcsecs, d->skipcount, skipsecs );
}


/*=============================================================================
    oss (volume) <signal>

Speaker/headphones output via open sound system.  Valid volumes are between 0
and 1; <0 means leave unchanged.

Unavailable if NO_OSS has been defined in sndsys.h.

See also: alsa, mike
=============================================================================*/

#ifndef NO_OSS

struct oss {
    int           dsp, bufind, bufsize;
    signed short  *buf;
};

float calcoss(sndobj *me, sndobj *, int);
void exitoss(sndobj *me);

#define OSS_OUTBUFSIZE  1024

sndobj *oss( float volume, sndobj *signal )
{
  sndobj *p;
  struct oss *d;
  int ioarg;
  
  p= newsndo( calcoss, signal->name, "oss", signal->nch, 1, signal );
  p->skip= dontskip;
  p->private[0]= d= new(struct oss);
  d->bufsize= 2*OSS_OUTBUFSIZE;
  p->private[1]= d->buf= news(signed short, d->bufsize);
  p->exit= exitoss;
  d->bufind= 0;
  d->dsp= dspdesc(0);
  if( d->dsp < 0 ) {
    fprintf(stderr, "Warning: oss: DSP could not be opened.  Returning signal.\n");
    return signal;
  }
  if( volume >= 0.0 ) {
    ioarg= (int)round(volume*100);
    if( ioarg > 100 )
      ioarg= 100;
    ioarg *= 257;
    ioctl( d->dsp, SOUND_MIXER_WRITE_PCM, &ioarg );
  }
  return p;
}

float calcoss(sndobj *me, sndobj *caller, int innr)
{
  struct oss *d;
  int ch;
  
  d= (struct oss *)me->private[0];
  INCALC(0);
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  if( d->bufind >= d->bufsize ) {
    write( d->dsp, d->buf, sizeof(short)*d->bufsize );
    d->bufind= 0;
  }
  if( me->nch > 1 ) {
    d->buf[d->bufind++]= toshort(INPUTC(0, 0));
    d->buf[d->bufind++]= toshort(INPUTC(0, 1));
  }
  else {
    d->buf[d->bufind++]= toshort(INPUT(0));
    d->buf[d->bufind++]= toshort(INPUT(0));
  }
  CALCRETURN;
}

void exitoss(sndobj *me)
{
  struct oss *d;
  
  d= (struct oss *)me->private[0];
  write( d->dsp, d->buf, sizeof(short)*d->bufind );
}

#endif


/*=============================================================================
    alsa "device" <signal>

Outputs the <signal> using ALSA.  "device" is the ALSA device to be used for
playback.  Use of a plughw:... device is recommended, as this sndobj does not
do any conversions.  Sample data are passed to ALSA as float values.  Sample
format conversion and resampling (as far as necessary) is left to ALSA.  The
number of channels requested from ALSA is the number of channels in <signal>.
As far as sndsys is concerned, the data are passed through unchanged.

This sndobj is unavailable if NO_ALSA is #defined in sndsys.h.

See also: oss, mike
=============================================================================*/

#ifndef NO_ALSA

struct alsa {
  snd_pcm_t		*pcm;
  snd_pcm_hw_params_t	*pars;
};

float calcalsa(sndobj *me, sndobj *, int);
void exitalsa(sndobj *me);

sndobj *alsa(const char *device, sndobj *signal)
{
  sndobj *p;
  struct alsa *d;
  unsigned srate;
  int err;

  p= newsndo(calcalsa, signal->name, "alsa", signal->nch, 1, signal);
  p->skip= dontskip;
  p->exit= exitalsa;
  p->private[0]= d= new(struct alsa);
  d->pcm= NULL;
  err= snd_pcm_hw_params_malloc(&d->pars);
  if( err<0 ) {
    fprintf(stderr, "Warning: alsa: Error allocating PCM parameter structure: %s.  Returning signal.\n", snd_strerror(err));
    return signal;
  }
  p->private[1]= d->pars;
  err= snd_pcm_open(&d->pcm, device, SND_PCM_STREAM_PLAYBACK, SND_PCM_NONBLOCK);
  if( err<0 ) {
    fprintf(stderr, "Warning: alsa: Error opening device `%s': %s.  Returning signal.\n", device, snd_strerror(err));
    return signal;
  }
  err= snd_pcm_hw_params_any(d->pcm, d->pars);
  if( err<0 ) {
    fprintf(stderr, "Warning: alsa: Error initialising PCM parameter structure: %s.  Returning signal.\n", snd_strerror(err));
    return signal;
  }
  err= snd_pcm_hw_params_set_access(d->pcm, d->pars, SND_PCM_ACCESS_RW_INTERLEAVED);
  if( err<0 ) {
    fprintf(stderr, "Warning: alsa: Could not set interleaved access: %s.  Returning signal.\n", snd_strerror(err));
    return signal;
  }
  err= snd_pcm_hw_params_set_channels(d->pcm, d->pars, p->nch);
  if( err<0 ) {
    fprintf(stderr, "Warning: alsa: Could not set number of channels to %d: %s.  Returning signal.\n", p->nch, snd_strerror(err));
    return signal;
  }
  err= snd_pcm_hw_params_set_format(d->pcm, d->pars, SND_PCM_FORMAT_FLOAT);
  if( err<0 ) {
    fprintf(stderr, "Warning: alsa: Could not set sample format to float: %s.  Consider using a plughw:... device.  Returning signal.\n", snd_strerror(err));
    return signal;
  }
  srate= SAMPRATE;
  err= snd_pcm_hw_params_set_rate_near(d->pcm, d->pars, &srate, NULL);
  if( err<0 ) {
    fprintf(stderr, "Warning: alsa: Could not set sampling rate to %u: %s.  Returning signal.\n", (unsigned)SAMPRATE, snd_strerror(err));
    return signal;
  }
  if( srate!=SAMPRATE )
    fprintf(stderr, "Warning: alsa: Requested sampling rate %u Hz, got %u Hz.\n", (unsigned)SAMPRATE, srate);
  err= snd_pcm_hw_params(d->pcm, d->pars);
  if( err<0 ) {
    fprintf(stderr, "Warning: alsa: Could not activate hardware parameters: %s.  Returning signal.\n", snd_strerror(err));
    return signal;
  }
  return p;
}

float calcalsa(sndobj *me, sndobj *caller, int innr)
{
  struct alsa *d;
  struct timeval waittime;
  int ch, err;

  d= (struct alsa *)me->private[0];
  INCALC(0);
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  while( 13 ) {
    err= snd_pcm_writei(d->pcm, me->inch[0], 1);
    if( err==1 )	// 1 sample transmitted
      break;
    if( err == -EPIPE ) {
      err= snd_pcm_prepare(d->pcm);
      if( err< 0 )
	break;		// give up on this sample value
    }
    else if( err == -ESTRPIPE ) {
      while( (err= snd_pcm_resume(d->pcm)) == -EAGAIN ) {
	waittime.tv_sec= 0;
	waittime.tv_usec= 5*1000000/SAMPRATE;
	select(0, NULL, NULL, NULL, &waittime);
      }
      if( err< 0 ) {
	err= snd_pcm_prepare(d->pcm);
	if( err< 0 )
	  break;		// give up on this sample value
      }
    }
    waittime.tv_sec= 0;
    waittime.tv_usec= 5*1000000/SAMPRATE;
    select(0, NULL, NULL, NULL, &waittime);
  }
  CALCRETURN;
}

void exitalsa(sndobj *me)
{
  struct alsa *d;
  struct timeval waittime;
  int loop;

  d= (struct alsa *)me->private[0];
  if( d->pcm ) {
    snd_pcm_drain(d->pcm);
    // wait for a maximum of 2 seconds in 20th second steps
    for( loop= 0; snd_pcm_state(d->pcm)==SND_PCM_STATE_RUNNING && loop< 40; ++loop ) {
      waittime.tv_sec= 0;
      waittime.tv_usec= 50000;
      select(0, NULL, NULL, NULL, &waittime);
    }
    snd_pcm_close(d->pcm);
  }
}

#endif


/*=============================================================================
    meantime <signal>

Computes the location in time of a signal.  This is done by interpreting its
modulus as a probability density function with respect to which the temporal
average is computed.  The signal is passed through unchanged; the time average
is output when processing has terminated.
=============================================================================*/

struct meantime {
  double time;
  double timesum[PLENTY], norm[PLENTY];
};

float calcmeantime(sndobj *me, sndobj *, int);
void exitmeantime(sndobj *me);

sndobj *meantime( sndobj *signal )
{
  sndobj *p;
  struct meantime *d;
  int ch;

  p= newsndo( calcmeantime, signal->name, "meantime", signal->nch, 1, signal );
  p->exit= exitmeantime;
  p->private[0]= d= new(struct meantime);
  d->time= 0.0;
  for( ch= 0; ch< p->nch; ++ch )
    d->timesum[ch]= d->norm[ch]= 0.0;
  return p;
}

float calcmeantime(sndobj *me, sndobj *caller, int innr)
{
  struct meantime *d;
  double ampl;
  int ch;

  d= (struct meantime *)me->private[0];
  INCALC(0);
  FORCH {
    ampl= (double)fabsf(INPUTC(0, ch));
    d->timesum[ch] += ampl * d->time;
    d->norm[ch] += ampl;
    OUTPUT(ch, INPUTC(0, ch));
  }
  d->time += 1.0/(double)SAMPRATE;
  CALCRETURN;
}

void exitmeantime(sndobj *me)
{
  struct meantime *d;
  double chsum, chnorm;
  int sample, ch;

  d= (struct meantime *)me->private[0];
  chsum= chnorm= 0.0;
  FORCH {
    chsum += d->timesum[ch];
    chnorm += d->norm[ch];
  }
  if( chnorm == 0.0 ) {
    printf("Temporal location of signal from `%s', type `%s' unknown - all channels permanently zero\n", me->in[0]->name, me->in[0]->type );
    return;
  }
  sample= (int)floor(chsum/chnorm*SAMPRATE + 0.5);
  printf("Temporal location of signal from `%s', type `%s':\n  Channel average:"
      " %g seconds (sample %d)\n  By channel: ",
      me->in[0]->name, me->in[0]->type, chsum/chnorm, sample );
  for( ch= 0; ch< me->nch-1; ++ch )
    if( d->norm[ch] > 0.0 ) {
      sample= (int)floor(d->timesum[ch]/d->norm[ch]*SAMPRATE + 0.5);
      printf("%g s (%d), ", d->timesum[ch]/d->norm[ch], sample );
    }
    else
      printf("-- (zero), ");
  if( d->norm[ch] > 0.0 ) {
    sample= (int)floor(d->timesum[ch]/d->norm[ch]*SAMPRATE + 0.5);
    printf("%g s (%d)\n", d->timesum[ch]/d->norm[ch], sample );
  }
  else
    printf("-- (zero)\n");
  return;
}


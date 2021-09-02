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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "sndsys.h"

/*
The objects in this file are asyncronous.  The fanout mechanis does not work
with them, so they must not share inputs with other objects.  Use with care!  */

/*=============================================================================
    linresample <factor> <signal>

Resampling by linear interpolation of sample data.  Factor > 1 means coarser
sampling, ie lower sampling rate or higher frequency.

See also: wlinresample
=============================================================================*/

struct linresample {
  double    fpos;
  float     *buf;
};

float calclinresample(sndobj *me, sndobj *, int);

sndobj *linresample(sndobj *factor, sndobj *in)
{
  sndobj *p;
  struct linresample *d;

  p= newsndo(calclinresample, in->name, "linresample", in->nch, 2, in, factor);
  p->private[0]= d= (struct linresample *)malloc(sizeof(struct linresample));
  p->private[1]= d->buf= (float*)malloc(2*p->nch*sizeof(float));
  d->fpos= 2.0;
  return p;
}


float calclinresample(sndobj *me, sndobj *caller, int innr)
{
  struct linresample *d;
  float *write;
  int ch, i;

  d= (struct linresample *)me->private[0];
  if( d->fpos >= 2.0 )
  {
    i= (int)floor(d->fpos) - 2;
    d->fpos -= floor(d->fpos);
    for( ; i>0; --i )
      INCALC(0);
    write= d->buf;
    *write++ = INCALC(0);
    OUTPUT(0, write[-1]);
    for( ch= 1; ch< me->nch; ++ch ) {
      *write++ = INPUTC(0, ch);
      OUTPUT(ch, INPUTC(0, ch));
    }
    *write++ = INCALC(0);
    for( ch= 1; ch< me->nch; ++ch )
      *write++ = INPUTC(0, ch);
  }
  else if( d->fpos >= 1.0 )
  {
    write= d->buf;
    INCALC(0);
    for( ch= 0; ch < me->nch; ++ch ) {
      *write = write[me->nch];
      write[me->nch]= INPUTC(0, ch);
      ++write;
    }
    d->fpos -= 1.0;
  }
  write= d->buf;
  for( ch= 0; ch< me->nch; ++ch, ++write )
    OUTPUT(ch, *write + d->fpos*(write[me->nch]-*write));
  d->fpos += (double)INCALC(1);
  CALCRETURN;
}


/*=============================================================================
    head (interval) <signal>

Stops evaluating its input after a given time and returns zero afterwards.  For
initialisation of delay loop buffers.

See also: repeat
=============================================================================*/

float calchead(sndobj *me, sndobj *, int);

sndobj *head( double interval, sndobj *in )
{
  sndobj *p;
  
  p= newsndo(calchead, "head", "head", in->nch, 1, in );
  p->private[0]= new(long);
  *(long*)p->private[0]= (long)round(interval*SAMPRATE);
  return p;
}


float calchead(sndobj *me, sndobj *caller, int innr)
{
  int ch;

  if( *(long*)me->private[0] >= 0 ) {
    if( *(long*)me->private[0] > 0 ) {
      INCALC(0);
      FORCH
	OUTPUT(ch, INPUTC(0, ch));
    }
    else
      FORCH
	OUTPUT(ch, 0.0);
    -- *(long*)me->private[0];
  }
  CALCRETURN;
}


/*=============================================================================
    repeat (interval) <freq> <signal>

Stops evaluating its input after a given time and repeats the sequence of
sample values afterwards.  For expensively created waveforms.  If freq is
non-null, the interval is assumed to contain one wavelength which will be
linearly interpolated.  freq is evaluated syncronously.  The original signal
should not be interpolated, or aliasing effects will be heard (at least if the
desired frequency is close to but not equal to the one from which the data
were interpolated).

See also: head
=============================================================================*/

struct repeat {
  long      size, pos;
  int       freqvalid;
  double    fpos, interval;
};

float calcrepeat(sndobj *me, sndobj *, int);

sndobj *repeat( double interval, sndobj *freq, sndobj *in )
{
  sndobj *p;
  struct repeat *d;
  
  p= newsndo(calcrepeat, "repeat", "repeat", in->nch, 1, in );
  buf( 0.0, interval, 0.0, in->nch, p, 0 );
  p->private[0]= d= new(struct repeat);
  d->size= (long)round(interval*SAMPRATE);
  d->pos= 0;
  d->fpos= 0.0;
  d->interval= interval;
  if( freq ) {
    addinput(p, freq);
    d->freqvalid= 1;
  }
  else
    d->freqvalid= 0;
  return p;
}


float calcrepeat(sndobj *me, sndobj *caller, int innr)
{
  struct repeat *d;
  int ch;

  d= (struct repeat *)me->private[0];
  BUFINCALC( 0, d->pos, d->fpos );
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  if( d->freqvalid ) {
    d->fpos += d->interval * INCALC(1);
    if( d->fpos > 1.0 ) {
      d->pos += (long)floor(d->fpos);
      d->fpos -= floor(d->fpos);
    }
  }
  else
    ++d->pos;
  while( d->pos >= d->size )
    d->pos -= d->size;
  CALCRETURN;
}



/*=============================================================================
    convolve (klength) (advance) <kernel> <signal>

Convolution of a <kernel> with the <signal>.  The length of the interval to be
convoluted is (klength) samples.  After each convolution, <advance> new samples
are read from <kernel>.  (Before the first convolution, <klength> are read.)
The <signal> input is processed syncronously.  For large kernels, this is very,
very, very expensive computationally.

See also: filter
=============================================================================*/

struct convolve {
  int	length, advance, kernelch, klen1, klen2, getkptr;
  float	*kbuf1, *kbuf2;
};

float calcconvolve( sndobj *me, sndobj *, int );

sndobj *convolve( int klength, int advance, sndobj *kernel, sndobj *signal )
{
  sndobj *p;
  struct convolve *d;
  
  if( klength < 1 || advance < 0 ) {
    fprintf(stderr, "convolve `%s': Warning: kernel length (%d) has to be >= 1, advancement (%d) has to be >= 0.  Returning signal.\n", signal->name, klength, advance);
    return signal;
  }
  p= newsndo( calcconvolve, signal->name, "convolve", signal->nch, 2, kernel, signal );
  p->skip= dontskip;
  buf( 0.0, (double)klength/(double)SAMPRATE, 0.0, kernel->nch, p, 0 );
  buf( (double)klength/(double)SAMPRATE, 0.0, 0.0, signal->nch, p, 1 );
  p->private[0]= d= new(struct convolve);
  d->length= klength;
  d->advance= advance;
  d->getkptr= 1;
  return p;
}

float calcconvolve( sndobj *me, sndobj *caller, int innr )
{
  struct convolve *d;
  float *sigread, *sigptr1, *sigptr2;
  int ch, ch0, count, siglen1, siglen2, sigch;
  int section1len, section3len;

  d= (struct convolve *)me->private[0];
  if( d->getkptr ) {
    BUFINCALC( 0, d->length-1, 0.0 );
    BUFPTR(0, 0, &d->kbuf1, &d->klen1, &d->kernelch);
    BUFPTR(0, d->klen1, &d->kbuf2, &d->klen2, &d->kernelch);
    if( d->klen1+d->klen2 < d->length ) {
      fprintf(stderr, "convolve `%s': Error: cannot access entire convolution kernel buffer\n", me->name);
      fprintf(stderr, "got %d + %d, want %d\n", d->klen1, d->klen2, d->length);
      exit(-1);
    }
    if( d->klen1 > d->length )
      d->klen1= d->length;
    if( d->klen1+d->klen2 > d->length )
      d->klen2= d->length - d->klen1;
    if( d->advance )
      BUFINCR(0, d->advance);
    else
      d->getkptr= 0;
  }
  BUFINCALC( 1, 0, 0.0 );
  FORCH
    me->ch[ch]= 0.0;
  BUFPTR( 1, -d->length+1, &sigptr1, &siglen1, &sigch );
  BUFPTR( 1, -d->length+1+siglen1, &sigptr2, &siglen2, &sigch );
  if( siglen1+siglen2 < d->length ) {
    fprintf(stderr, "convolve `%s': Error: cannot access entire signal buffer\n", me->name);
    exit(-1);
  }
  if( siglen1 > d->length )
    siglen1= d->length;
  if( siglen1+siglen2 > d->length )
    siglen2= d->length - siglen1;
  section1len= siglen2 > d->klen1 ? d->klen1 : siglen2;
  for( count= 0; count< section1len; ++count )
    for( ch= 0, ch0= 0; ch< me->nch; ++ch, ++ch0 ) {
      if( ch0 >= d->kernelch )
	ch0= 0;
      me->ch[ch] += d->kbuf1[d->kernelch*count + ch0] * 
			sigptr2[sigch*(siglen2-count) + ch];
    }
  if( siglen2 > d->klen1 ) {
    for( ; count< siglen2; ++count )
      for( ch= 0, ch0= 0; ch< me->nch; ++ch, ++ch0 ) {
	if( ch0 >= d->kernelch )
	  ch0= 0;
	me->ch[ch] += d->kbuf2[d->kernelch*(count-d->klen1) + ch0] * 
			  sigptr2[sigch*(siglen2-count) + ch];
      }
    section3len= siglen1;
  }
  else {
    for( ; count< d->klen1; ++count )
      for( ch= 0, ch0= 0; ch< me->nch; ++ch, ++ch0 ) {
	if( ch0 >= d->kernelch )
	  ch0= 0;
	me->ch[ch] += d->kbuf1[d->kernelch*count + ch0] * 
			  sigptr1[sigch*(d->length-1-count) + ch];
      }
    section3len= d->klen2;
  }
  for( count= section3len; count> 0; --count )
    for( ch= 0, ch0= 0; ch< me->nch; ++ch, ++ch0 ) {
      if( ch0 >= d->kernelch )
	ch0= 0;
      me->ch[ch] += d->kbuf2[d->kernelch*(d->klen2-count) + ch0] * 
			sigptr1[sigch*(count-1) + ch];
    }

  BUFINCR( 1, 1 );
  CALCRETURN;
}



/*=============================================================================
    interleave (inigap) <gap> <signal>

Outputs the values of <signal> interleaved by the number of zero samples given
by <gap>.  Both inputs are evaluated asyncronously (w.r.t. the output), but are
syncronised relative to each other.  Both are evaluated once at the end of each
gap.  (inigap) gives the size of the initial gap before the first evaluation of
the inputs, in samples.  <gap> is rounded to the nearest integer; .5 rounds up
and negative values are equivalent to 0.

This object is intended for creating a series of cues.

See also: at, slide
=============================================================================*/

float calcinterleave(sndobj *me, sndobj *, int);

sndobj *interleave( int inigap, sndobj *gap, sndobj *signal )
{
  sndobj *p;
  p= newsndo( calcinterleave, "interleave", "interleave", signal->nch, 2, gap, signal);
  p->skip= dontskip;
  p->private[0]= new(signed long);
  *(signed long*)p->private[0]= (inigap > 0 ? inigap : 0);
  return p;
}

float calcinterleave(sndobj *me, sndobj *caller, int innr)
{
  signed long *gapcount;
  int ch;
  
  gapcount= (signed long *)me->private[0];
  if( *gapcount < 0 ) {
    FORCH
      OUTPUT(ch, 0.0);
    *gapcount= -*gapcount-1;
  }
  else if( !*gapcount ) {
    INCALC(1);
    FORCH
      OUTPUT(ch, INPUTC(1, ch));
    *gapcount= (signed long)floorf(INCALC(0) + 0.5f);
    if( *gapcount<0 )
      *gapcount= 0;
    else
      *gapcount= -*gapcount;
  }
  else 
    --*gapcount;
  CALCRETURN;
}


/*=============================================================================
    powersave <gate> <in1> <in2>

Mixing of two inputs depending on <gate>.  1 gives input1, 0 input 2.  Unlike
mix2, powersave it does not evaluate the unneeded input if gate is 0 or 1.
Besides, in2 may be passed as NULL if only gating is desired.  Then zero is
output if gate is 0.  

For syncronous objects, this can be done with mix2 since skipping was
introduced.  The powersave object offers a way to switch off asyncronous
objects which typically have "dontskip" as their skip function.
=============================================================================*/

float calcpowersave(sndobj *me, sndobj *, int);

sndobj *powersave( sndobj *gate, sndobj *in1, sndobj *in2 )
{
  if( in2 )
    return newsndo( calcpowersave, "powersave", "powersave", 
	    in1->nch> in2->nch? in1->nch: in2->nch, 3, gate, in1, in2 );
  else
    return newsndo( calcpowersave, "powersave", "powersave", 
		    in1->nch, 2, gate, in1 );
}


float calcpowersave(sndobj *me, sndobj *caller, int innr)
{
  float gateval;
  int ch;
  
  gateval= INCALC(0);
  if( gateval>= 1.0 ) {
    INCALC(1);
    FORCH
      OUTPUT(ch, INPUTC(1, ch));
  }
  else if( gateval <= 0.0 ) {
    if( me->nin==3 ) {
      INCALC(2);
      FORCH
	OUTPUT(ch, INPUTC(2, ch));
    }
    else
      FORCH
	OUTPUT(ch, 0.0);
  }
  else {
    INCALC(1);
    if( me->nin==3 ) {
      INCALC(2);
      FORCH
	OUTPUT(ch, INPUTC(2, ch) + gateval*(INPUTC(1, ch)-INPUTC(2, ch)));
    }
    else
      FORCH
	OUTPUT(ch, gateval*INPUTC(1, ch));
  }
  CALCRETURN;
}


/*=============================================================================
    reshape (maxtime) (mode) <shapesig> <signal>

This object reshapes the half-waves of a <signal> according to its other
(first) input, <shapesig>.  This is asyncronous w. r. t. <shapesig> since its
half-waves do not necessarily have the same length as those of <signal>.
<mode> is the mode of operation: for RESHAPE_SLOPE, a half-wave is a time
interval in which the signal slopes upwards or downwards.  For the other two
modes, it is an interval between zero-crossings.  RESHAPE_KEEPAMPL retains the
amplitude of the half-waves, which leads to kinks at the zero-crossings.
RESHAPE_KEEPANGLE retains the aspect ratio of <shapesig>'s half-waves, which
may lead to high peaks in the output.  RESHAPE_KEEPRMS scales the new
half-waves so that the RMS (root of mean square) is the same.  <maxtime> is the
maximal duration of a half-wave, i. e. half the period of the minimal
frequency.  This is needed for determining buffer sizes.  If it is too small,
this will lead to artifacts; beware of slow "floating" of an approximately
constant (or zero) signal.  Only the first channel of the signal is processed.

See also: repeathw, avghw
=============================================================================*/

struct reshape {
  int	mode, maxind, sigcount, shapecount;
  double scale, offset;
  float sigpeak, shapepeak;
  double shapeind, shapeincr;
};

float calcreshape(sndobj *me, sndobj *, int);

sndobj *reshape( double maxtime, int mode, sndobj *shapesig, sndobj *signal )
{
  sndobj *p;
  struct reshape *d;
  
  p= newsndo( calcreshape, "reshape", "reshape", 1, 2, shapesig, signal );
  p->skip= dontskip;
  buf( 1.0/(double)SAMPRATE, maxtime, 0.0, 1, p, 0);
  buf( 1.0/(double)SAMPRATE, maxtime, 0.0, 1, p, 1);
  p->private[0]= d= new(struct reshape);
  d->maxind= (int)floor(SAMPRATE*maxtime);
  d->mode= mode;
  d->sigcount= 0;
  d->scale= d->offset= 0.0;
  d->sigpeak= d->shapepeak= 0.0;
  d->shapeind= d->shapeincr= 0.0;
  return p;
}

float calcreshape(sndobj *me, sndobj *caller, int innr)
{
  static char debugmsg[1024];
  static int indtrunc= 0;

  struct reshape *d;

  d= (struct reshape *)me->private[0];
  if( d->sigcount <= 0 ) {
    debugmsg[0]= debugmsg[1023]= 0;
    if( d->mode == RESHAPE_SLOPE )
    {
      int searchind, safety;
      float sigstart, prev, sigend, shapestart, shapeend;
      
      BUFINCALC(1, 0, 0.0);
      if( INPUT(1) == d->sigpeak ) {
        OUTPUT(0, d->sigpeak);
        BUFINCR(1, 1);
        CALCRETURN;
      }
      sigstart= prev= d->sigpeak;
      if( INPUT(1) > prev )	// rising half-wave
        for( searchind= 1; searchind< d->maxind; ++searchind ) {
          prev= INPUT(1);
  	  BUFINCALC(1, searchind, 0.0);
  	  if( INPUT(1) <= prev )
  	    break;
        }
      else			// falling half-wave
        for( searchind= 1; searchind< d->maxind; ++searchind ) {
          prev= INPUT(1);
          BUFINCALC(1, searchind, 0.0);
  	  if( INPUT(1) >= prev )
  	    break;
        }
      sigend= prev;
      d->sigpeak= prev;
      d->sigcount= searchind;
      if( d->shapecount > 0 )
      // last BUFINCR may not have taken place due to roundoff error
        BUFINCR(0, d->shapecount);
      BUFINCALC(0, 0, 0.0);
      shapestart= prev= d->shapepeak;
      for( safety= (int)SAMPRATE*10; INPUT(0)==prev && safety> 0; --safety )
      {				    // guard against infinite loop
        BUFINCR(0, 1);
        BUFINCALC(0, 0, 0.0);
      }
      if( INPUT(0) > prev )	// rising half-wave
        for( searchind= 1; searchind < d->maxind; ++searchind ) {
          prev= INPUT(0);
  	  BUFINCALC(0, searchind, 0.0);
	  if( INPUT(0) <= prev )
  	    break;
        }
      else			// falling half-wave
        for( searchind= 1; searchind < d->maxind; ++searchind ) {
          prev= INPUT(0);
          BUFINCALC(0, searchind, 0.0);
  	  if( INPUT(0) >= prev )
  	    break;
        }
      d->shapecount= searchind;
      shapeend= prev;
      d->shapepeak= prev;
      d->shapeincr= (double)searchind / (double)d->sigcount;
      if( shapeend==shapestart )
	d->scale= 1.0;
      else
	d->scale= (double)(sigend-sigstart)/(double)(shapeend-shapestart);
      d->offset= sigstart - d->scale * shapestart;
      d->shapeind= 0.0;
    }
    else
    {
      int searchind, safety;
      float previous, maxsig, maxshape;
      double fsigcount, fshapecount;
      
      BUFINCALC(1, -1, 0.0);
      previous= INPUT(1);
      BUFINCALC(1, 0, 0.0);
      if( INPUT(1) == 0.0 ) {
        OUTPUT(0, 0.0);
        BUFINCR(1, 1);
        CALCRETURN;
      }
      fsigcount= (double)INPUT(1) / (INPUT(1) - previous);
      previous= INPUT(1);
      maxsig= INPUT(1);
      if( INPUT(1) > 0 )	// upper half-wave
        for( searchind= 1; searchind< d->maxind; ++searchind ) {
  	  BUFINCALC(1, searchind, 0.0);
  	  if( INPUT(1) <= 0 )
  	    break;
	  if( INPUT(1) > maxsig )
	    maxsig= INPUT(1);
	  previous= INPUT(1);
	  fsigcount += 1.0;
        }
      else			// lower half-wave
        for( searchind= 1; searchind< d->maxind; ++searchind ) {
          BUFINCALC(1, searchind, 0.0);
  	  if( INPUT(1) >= 0 )
  	    break;
	  if( INPUT(1) < maxsig )
	    maxsig= INPUT(1);
	  previous= INPUT(1);
	  fsigcount += 1.0;
        }
      d->sigcount= searchind;
      fsigcount += (double)previous / (previous - INPUT(1));
      if( d->shapecount > 0 )
      // last BUFINCR may not have taken place due to roundoff error
        BUFINCR(0, d->shapecount);
      BUFINCALC(0, -1, 0.0);
      previous= INPUT(0);
      BUFINCALC(0, 0, 0.0);
      for( safety= (int)SAMPRATE*10; INPUT(0)==0.0 && safety> 0; --safety )
      {				    // guard against infinite loop
        BUFINCR(0, 1);
        BUFINCALC(0, 0, 0.0);
	previous= 0.0f;
      }
      d->shapeind= (double)previous / (previous - INPUT(0));
      fshapecount= (double)INPUT(0) / (INPUT(0) - previous);
      previous= INPUT(0);
      maxshape= INPUT(0);
      if( INPUT(0) > 0 )	// upper half-wave
        for( searchind= 1; searchind < d->maxind; ++searchind ) {
  	  BUFINCALC(0, searchind, 0.0);
	  if( INPUT(0) <= 0 )
  	    break;
  	  if( INPUT(0) > maxshape )
	    maxshape= INPUT(0);
	  previous= INPUT(0);
	  fshapecount += 1.0;
        }
      else			// lower half-wave
	for( searchind= 1; searchind < d->maxind; ++searchind ) {
          BUFINCALC(0, searchind, 0.0);
  	  if( INPUT(0) >= 0 )
  	    break;
	  if( INPUT(0) < maxshape )
	    maxshape= INPUT(0);
	  previous= INPUT(0);
	  fshapecount += 1.0;
        }
      d->shapecount= searchind;
      fshapecount += (double)previous / (previous - INPUT(0));
      d->shapeincr= fshapecount / fsigcount;
      if( d->mode == RESHAPE_KEEPRMS )
      {
        int count;
	double sigsqrsum, shapesqrsum;
	
	sigsqrsum= 0.0;
	for( count= 0; count< d->sigcount; ++count ) {
	  BUFINCALC(1, count, 0.0);
	  sigsqrsum += INPUT(1)*INPUT(1);
	}
	shapesqrsum= 0.0;
	for( count= 0; count< d->shapecount; ++count ) {
	  BUFINCALC(0, count, 0.0);
	  shapesqrsum += INPUT(0)*INPUT(0);
	}
	if( !shapesqrsum )
	  d->scale= 1;
	else
	  d->scale= (float)sqrt(sigsqrsum/shapesqrsum*d->shapeincr);
	// = sqrt((sigsqrsum/d->sigcount)/(shapesqrsum/searchind))
      }
      else if( d->mode  == RESHAPE_KEEPANGLE )
        d->scale= 1.0/d->shapeincr;
      else	// mode == RESHAPE_KEEPAMPL
	if( !maxshape )
	  d->scale= 1;
	else
	  d->scale= maxsig/maxshape;
      // d->offset= 0.0; from initialisation
//printf("sigcount= %d, shapecount= %d, maxsig= %g, maxshape= %g\n", d->sigcount, searchind, maxsig, maxshape);
    }
  }
  BUFINCALC(0, -1, d->shapeind);
  OUTPUT(0, d->offset + d->scale*INPUT(0));
  d->shapeind += d->shapeincr;
  if( d->shapeind >= 1.0 )
    if( d->shapecount >= (int)trunc(d->shapeind) ) {
      int increment= (int)trunc(d->shapeind);
      BUFINCR(0, increment);
      d->shapecount -= increment;
      d->shapeind -= (double)increment;
    }
    else {
      BUFINCR(0, d->shapecount);
      d->shapecount= 0;
      d->shapeind= 1.0;
    }
  --d->sigcount;
  BUFINCR(1, 1);
  CALCRETURN;
}


/*=============================================================================
    repeathw (maxlength) <repetitions> <grouplength> <signal>

This sndobj repeats a signal's half-waves (the parts between zero crossings)
for a number of times.  <repetitions> gives the number of repetitions.
<grouplength> is usually zero.  Otherwise it gives the length in time of a
group of half-waves which is to be repeated.  Before every repetition is
started, the number of half-waves closest to the requested time interval is
marked for repetition.  (maxlength) gives the maximum value for <grouplength>.
Changes in <repetitions> are noticed after each repetition, changes in
<grouplength> after each full set of repetitions.  The half-waves of the output
alternate in sign, as for a normal oscillation.

This object is syncronous w. r. t. its first two (parameter) inputs.  It works
only on the first signal channel.

For large values of <grouplength> (100ths of a second), this gives a slightly
distorted granular expansion.  Small values result in distortion for
<repetitions> ~ 2, strange bubbly sounds for larger <repetitions> ~ 5.  Still
larger values for <repetitions> give a melody with pitches corresponding to the
length of half-waves.  This is definitely a sndobj which warrants a lot of
experimentation.

See also: reshape, avghw
=============================================================================*/

struct repeathw {
  int	maxind;
  float	glength, sign;
  int	reps, gsamples;
  int	samcount, repcount;
};

float calcrepeathw(sndobj *me, sndobj *, int);

sndobj *repeathw(double maxlength, sndobj *repetitions, sndobj *grouplength, sndobj *signal)
{
  sndobj *p;
  struct repeathw *d;
  
  p= newsndo(calcrepeathw, "repeathw", "repeathw", 1, 3, repetitions, grouplength, signal);
  p->skip= dontskip;
  buf( 0.0, maxlength, 0.0, 1, p, 2);
  p->private[0]= d= new(struct repeathw);
  d->maxind= (int)floor(SAMPRATE*maxlength);
  d->glength= 0.0;
  d->sign= 1.0;
  d->gsamples= d->reps= 0;
  d->samcount= d->repcount= INT_MAX;
  return p;
}

float calcrepeathw(sndobj *me, sndobj *caller, int innr)
{
  struct repeathw *d;
  
  d= (struct repeathw *)me->private[0];
  if( d->samcount >= d->gsamples )
  {
    d->samcount= 0;
    d->reps= (int)floorf(INCALC(0)+0.5);
    if( d->repcount >= d->reps )
    {
      int searchind, afterlastzero, sign;
      
      d->repcount= 1;
      BUFINCR(2, d->gsamples);
      d->glength= INCALC(1);	// <=0 is OK, handled gracefully below
      BUFINCALC(2, 0, 0.0);
      afterlastzero= 0;
      if( INPUT(2)> 0 )
        sign= 1;
      else if( INPUT(2)< 0 )
        sign= -1;
      else
        sign= 0;
      for( searchind= 1; searchind< d->maxind; ++searchind ) {
	BUFINCALC(2, searchind, 0.0);
	if( (sign>= 0 && INPUT(2)<= 0.0) || (sign<= 0 && INPUT(2)>= 0.0) ) {
	  if( (double)searchind/(double)SAMPRATE > d->glength )
	    break;
	  afterlastzero= searchind;
          if( INPUT(2)> 0 )
            sign= 1;
          else if( INPUT(2)< 0 )
            sign= -1;
          else
            sign= 0;
	}
      }
      if( !afterlastzero || (searchind < d->maxind &&
	    (double)searchind/(double)SAMPRATE - d->glength < 
	    d->glength - (double)afterlastzero/(double)SAMPRATE) )
        d->gsamples= searchind;
      else
        d->gsamples= afterlastzero;
    }
    else {
      INSKIP(1);
      ++d->repcount;
      d->sign= -d->sign;
    }
  }
  else {
    INSKIP(0);
    INSKIP(1);
  }
  BUFINCALC(2, d->samcount, 0.0);
  OUTPUT(0, d->sign * INPUT(2));
  ++d->samcount;
  CALCRETURN;
}


/*=============================================================================
    voiceact (threshold) (gap) <cmp> <signal>

This object only returns samples which exceed a certain threshold and ignores
the others.  (threshold) is the noise level which must be reached or exceeded
by <cmp>.  (gap) is the length of the period of silence which is let through to
separate above-threshold stretches.  When the modulus of <cmp> falls below the
threshold, up to <gap> seconds of the <signal> are passed through, before
<signal> is ignored until <cmp> exceeds the threshold again.  The gap between
recordings is restricted to be at most twice as long as the preceding recorded
stretch - this avoids prolongation of the gap by brief crackling noise.

<cmp> and <signal> are evaluated syncronously w.r.t. each other (<cmp> will in
many applications be the RMS of <signal>).  If <cmp> has multiple channels, it
is sufficient for one to exceed the threshold to enable output.
=============================================================================*/

#define VOICEACT_TIMEOUT	(10*60*SAMPRATE)

struct voiceact {
  float noiselvl;
  int   gap, gapcount, passcount;
};

float calcvoiceact(sndobj *me, sndobj *, int);

sndobj *voiceact( float threshold, double gap, sndobj *cmp, sndobj *signal )
{
  sndobj *p;
  struct voiceact *d;

  p= newsndo(calcvoiceact, "voiceact", "voiceact", signal->nch, 2, cmp, signal);
  p->skip= dontskip;
  p->private[0]= d= new(struct voiceact);
  d->noiselvl= fabsf(threshold);
  d->gap= (int)ceil(gap*(double)SAMPRATE);
  d->gapcount= 0;
  d->passcount= d->gap;
  return p;
}


float calcvoiceact(sndobj *me, sndobj *caller, int innr)
{
  struct voiceact *d;
  int ch, accept, skipped;

  d= (struct voiceact*)me->private[0];
  if( d->gapcount< d->gap && d->gapcount < d->passcount ) {
    INCALC(0);
    accept= 0;
    for( ch= 0; ch< me->in[0]->nch; ++ch )
      if( fabsf(INPUTC(0, ch)) >= d->noiselvl ) {
	accept= 1;
        break;
      }
    if( accept ) {
      d->gapcount= 0;
      d->passcount= 0;
    }
    else
      ++d->gapcount;
    INCALC(1);
  }
  else {
    skipped= 0;
    while( 13 ) {
      INCALC(0);
      accept= 0;
      for( ch= 0; ch< me->in[0]->nch; ++ch )
	if( fabsf(INPUTC(0, ch)) >= d->noiselvl ) {
	  accept= 1;
	  d->gapcount= 0;
	  break;
	}
      ++skipped;
      if( accept || skipped >= VOICEACT_TIMEOUT )
        break;
      INSKIP(1);
    }
    INCALC(1);
    d->passcount += 2;
  }
  FORCH
    OUTPUT(ch, INPUTC(1, ch));
  CALCRETURN;
}



/*=============================================================================
    wlinresample (bottom) <factor> <signal>

Resamples a wavelet-transformed signal by <factor>.  A <factor> > 0 shifts
frequencies upwards, a <factor> < 1 downwards.  Each series of wavelet
coefficients is interpreted as a separate time series which is interpolated
linearly between existing values.

See also: wtany, linresample
=============================================================================*/

// from sndwt.c:
int wtdepth( double bottom );
int pow2div(int num);

struct wlinresample {
  int     depth, setlen, srcsetcount, destsetcount, depth1, initialised;
  float   *last, *next;
  double  srcpos;
};

float calcwlinresample(sndobj *me, sndobj *, int);

sndobj *wlinresample( double bottom, sndobj *factor, sndobj *signal )
{
  sndobj *p;
  struct wlinresample *d;

  p= newsndo( calcwlinresample, signal->name, "wlinresample", signal->nch, 2, factor, signal );
  p->skip= dontskip;
  p->private[0]= d= new(struct wlinresample);
  d->depth= wtdepth(bottom);
  d->depth1= d->depth+1;
  d->setlen= 1 << d->depth;
  buf( 0.0, 2.0*(double)d->setlen/(double)SAMPRATE, 0.0, p->nch, p, 1 );
  // might as well set jitter to 0 since we're asyncronous anyway
  d->srcsetcount= 0;
  d->destsetcount= 1;
  p->private[1]= d->last= news(float, (d->depth+1)*p->nch);
  p->private[2]= d->next= news(float, (d->depth+1)*p->nch);
  d->initialised= 0;
  d->srcpos= 1.0;
  return p;
}

float calcwlinresample(sndobj *me, sndobj *caller, int innr)
{
  struct wlinresample *d;
  float nextwgt;
  int scale, depth, ch;

  d= (struct wlinresample *)me->private[0];
  if( !d->initialised )
  {
    d->initialised= 1;
    for( depth= 0, scale= 1; depth<= d->depth; ++depth, scale <<= 1 ) {
      BUFINCALC(1, scale-1, 0.0);
      FORCH {
	d->last[ch*d->depth1 + depth]= INPUTC(1, ch);
	d->next[ch*d->depth1 + depth]= INPUTC(1, ch);
      }
    }
  }
  if( d->srcpos >= d->srcsetcount+1 ) {
    d->srcsetcount= (int)floor(d->srcpos);
    depth= pow2div(d->srcsetcount);
    if( depth< d->depth ) {
      BUFINCALC(1, (1 << (depth+1)) + d->srcsetcount-1, 0.0);
      FORCH {
	d->last[ch*d->depth1 + depth]= d->next[ch*d->depth1 + depth];
	d->next[ch*d->depth1 + depth]= INPUTC(1, ch);
      }
      if( depth==d->depth-1 ) {
	BUFINCALC(1, 2*d->setlen-1, 0.0);
	FORCH {
	  d->last[ch*d->depth1 + d->depth]= d->next[ch*d->depth1 + d->depth];
	  d->next[ch*d->depth1 + d->depth]= INPUTC(1, ch);
	}
      }
    }
    else {
      d->srcsetcount= 0;
      d->srcpos -= d->setlen;
      BUFINCR(1, d->setlen);
    }
  }
  depth= pow2div(d->destsetcount);
  if( depth==d->depth )
    scale= d->setlen;
  else
    scale= 1 << (depth+1);
  nextwgt= (((d->srcsetcount+scale/2)&(scale-1)) + d->srcpos - floor(d->srcpos))/(double)scale;
  FORCH
    OUTPUT(ch, d->last[ch*d->depth1+depth] + nextwgt *
		(d->next[ch*d->depth1+depth] - d->last[ch*d->depth1+depth]));
  if( ++d->destsetcount > d->setlen )
    d->destsetcount= 1;
  d->srcpos += INCALC(0);
  CALCRETURN;
}


/*=============================================================================
    wasyncinterleave (bottom) (ninputs) <input1> ...

The asyncronous version of a reverse of wuninterleave.  Unlike winterleave,
this object takes multiple inputs, which provide the time series of the
individual wavelet data.  All of them are evaluated asyncronously.  The first
represents the highest frequency band and is evaluated every other time, the
second the next lower band, evaluated every four output values, and so on.
(ninputs) gives the number of input signals.  If lower than the number of bands
computed from the (bottom) frequency (floor(log2(SAMPRATE/(bottom))+0.5)+1),
the remaining wavelet coefficients are set to 0.  The number of output channels
is the minimum of the number of the channels of all inputs, and every output
channel is generated from the corresponding channel of all inputs.  In contrast
to winterleave, wasyncinterleave adjusts for differences in normalisations
between wavelet data of different scales.

See also: winterleave, wuninterleave, wtany
=============================================================================*/

struct wasyncinterleave {
  int setcount, setlen, depth;
};

float calcwasyncinterleave(sndobj *me, sndobj *, int);

sndobj *wasyncinterleave(double bottom, int ninputs, ... )
{
  sndobj *p;
  struct wasyncinterleave *d;
  va_list in;
  int i, depth, maxinputs;

  depth= wtdepth(bottom);
  p= newsndo(calcwasyncinterleave, "wasyncinterleave", "wasyncinterleave", 1, 0);
  p->skip= dontskip;
  if( ninputs > PLENTY ) {
    fprintf(stderr, "Warning: add: Number of inputs (%d) limited to %d.\n", ninputs, (int)PLENTY );
    ninputs= PLENTY;
  }
  va_start(in, ninputs);
  maxinputs= depth+1 > PLENTY ? PLENTY : depth+1;
  if( ninputs > maxinputs ) {
    if( depth+1 <= PLENTY )
      fprintf(stderr, "Warning: wasyncinterleave: %d surplus inputs given - ignored.\n", ninputs-maxinputs );
    else
      fprintf(stderr, "Warning: wasyncinterleave: number of inputs limited to %d.  Last %d inputs ignored.\n", (int)PLENTY, ninputs-maxinputs );
    ninputs= maxinputs;
  }
  p->nch= depth+1;
  for( i= 0; i< ninputs; ++i ) {
    addinput(p, va_arg(in, sndobj*));
    if( p->in[i]->nch < p->nch )
      p->nch= p->in[i]->nch;
  }
  p->private[0]= d= new(struct wasyncinterleave);
  d->depth= depth;
  d->setlen= 1<<depth;
  d->setcount= 1;
  return p;
}

float calcwasyncinterleave(sndobj *me, sndobj *caller, int innr)
{
  struct wasyncinterleave *d;
  float coeff;
  int ch, currdepth;

  d= (struct wasyncinterleave *)me->private[0];
  currdepth= pow2div(d->setcount);
  if( currdepth >= me->nin )
    FORCH
      OUTPUT(ch, 0.0);
  else {
    coeff= (float)(1 << (currdepth/2));
    if( (currdepth&1) != 0 )
      coeff *= M_SQRT2;
    INCALC(currdepth);
    FORCH
      OUTPUT(ch, coeff * INPUTC(currdepth, ch));
  }
  if( ++d->setcount > d->setlen )
    d->setcount= 1;
  CALCRETURN;
}


/*=============================================================================
    resync (id) (mode) <signal> <factor>

This object allows in some cases to re-synchronise asynchronous objects or 
object trees.  To that end, resync objects with the same (id) have to be
inserted both into the asynchronous oject's input and output data streams.
(mode) has to be RESYNC_OUT for the single resync object which processes the
embedded object's output.  The others are normally of mode RESYNC_IN.  If (mode)
is RESYNC_TRIG, care is taken that non-zero triggers are not lost.  The modes
RESYNC_TIME and RESYNC_FREQ can take an additional argument <factor>.  This is
assumed to be the resampling factor of a linresample or similar object just
upstream of the RESYNC_OUT object.  The values passing through RESYNC_TIME
objects are multiplied, those passing through RESYNC_FREQ objects divided by
<factor> to compensate the effect of resampling on time and frequency inputs if
desired.  Note that <factor> is ignored for the other modes, and in particular
that RESYNC_OUT objects perform no resampling.

Method: resync objects are registered in linked lists, one for every (id).
When the RESYNC_OUT object's calc function is called, it calls the calc
functions of the inputs of all other resync objects with the same (id).  This
ensures that inputs are synchronous with the output whatever the object(s) in
between may do.  RESYNC_IN objects just copy their input values, losing or
duplicating some if their calc function is called less or more frequently than
that of the RESYNC_OUT object.  RESYNC_TRIG objects never lose a non-zero value
unless they are called so infrequently and there are so many non-zero values
that this is impossible.  RESYNC_TIME and RESYNC_FREQ objects adjust their
output values as described above.

Limitations: input resync objects should not be placed upstream of buffered
inputs or inputs which are frequently read ahead unless their value changes
slowly.  As their calculations are triggered by the RESYNC_OUT object, they
will just repeat the last value if they are being read ahead.  resync objects
are never skipped because the RESYNC_OUT object cannot know which input values
will be needed.
=============================================================================*/

struct resync {
  int		id, mode;
  struct resync	*nextid, *nextobj;
  sndobj	*obj;
  int		dirty;
};

static struct resync *resyncreg= NULL;

float calcresync(sndobj *me, sndobj *caller, int innr);
void incalcresync(sndobj *me);
void exitresync(sndobj *me);

sndobj *resync(int id, int mode, sndobj *signal, sndobj *factor)
{
  sndobj *p;
  struct resync *d, *search, **link;

  d= new(struct resync);
  d->id= id;
  d->mode= mode;
  if( resyncreg ) {
    for( link= &resyncreg, search= resyncreg; search->id!=id && search->nextid; 
	link= &search->nextid, search= search->nextid );
    if( mode==RESYNC_OUT && search->id==id && search->mode==RESYNC_OUT ) {
      fprintf(stderr, "Error: resync object with id %d and mode RESYNC_OUT already exists! Aborting.\n", id );
      exit(1);
    }
    if( search->id==id ) {	// insert in nextobj chain for this ID
      if( mode==RESYNC_OUT ) {	// insert RESYNC_OUT mode objects at the top
	d->nextid= search->nextid;
	search->nextid= NULL;
	d->nextobj= search;
	*link= d;
      }
      else {
	d->nextobj= search->nextobj;
	search->nextobj= d;
	d->nextid= NULL;
      }
    }
    else {	// first resync object with this ID - insert in nextid chain
      search->nextid= d;
      d->nextobj= d->nextid= NULL;
    }
  }
  else {
    resyncreg= d;
    d->nextobj= d->nextid= NULL;
  }
  if( !factor || (mode!=RESYNC_TIME && mode!=RESYNC_FREQ) ) {
    if( factor )
      fprintf(stderr, "resync: Warning: factor input can only be used with modes RESYNC_TIME and RESYNC_FREQ. Ignoring factor.\n");
    p= newsndo( calcresync, signal->name, "resync", signal->nch, 1, signal );
  }
  else
    p= newsndo( calcresync, signal->name, "resync", signal->nch, 2, signal, factor );
  d->obj= p;
  d->dirty= 0;
  p->private[0]= d;
  p->skip= dontskip;
  p->exit= exitresync;
  return p;
}


float calcresync(sndobj *me, sndobj *caller, int innr)
{
  struct resync *d;
  int ch;

  d= (struct resync *)me->private[0];
  if( d->mode==RESYNC_OUT ) {
    struct resync *inobj;

    for( inobj= d->nextobj; inobj; inobj= inobj->nextobj )
      incalcresync(inobj->obj);
    OUTPUT(0, INCALC(0));
    FORCH1
      OUTPUT(ch, INPUTC(0, ch));
  }
  else if( d->mode==RESYNC_TRIG )
    d->dirty= 0;
  CALCRETURN;
}


void incalcresync(sndobj *me)
{
  struct resync *d;
  int ch, ch1, non0;

  d= (struct resync *)me->private[0];
  INCALC(0);
  if( me->nin>1 ) {
    INCALC(1);
    if( d->mode==RESYNC_TIME )
      FORCH
	OUTPUT(ch, INPUTC(0, ch)*INPUT(1));
    else
      FORCH
	OUTPUT(ch, INPUTC(0, ch)/INPUT(1));
  }
  else if( d->mode==RESYNC_TRIG ) {
    if( !d->dirty ) {	// don't overwrite non-0 trigger which hasn't been read
      non0= 0;
      FORCH {
	OUTPUT(ch, INPUTC(0, ch));
	if( INPUTC(0, ch) )
	  non0= 1;
      }
      d->dirty= non0;
    }
  }
  else {
    FORCH
      OUTPUT(ch, INPUTC(0, ch));
  }
  return;
}


// Because some resync objects may have been created by mistake and may be
// removed by prune(), we need an exit function which neatly removes an object
// from the resync registry.
void exitresync(sndobj *me)
{
  struct resync *d, *search, **link;

  d= (struct resync *)me->private[0];
  if( !resyncreg )
    return;
  for( link= &resyncreg, search= resyncreg; 
      search->id!=d->id && search->nextid; 
      link= &search->nextid, search= search->nextid );
  if( search->id != d->id )
    return;
  if( search->obj == me ) {
    if( search->nextobj ) {
      search->nextobj->nextid= search->nextid;
      *link= search->nextobj;
    }
    else	// only object with this ID
      *link= search->nextid;
  }
  else {
    for( link= &search->nextobj, search= search->nextobj; 
	search && search->obj != me; 
	link= &search->nextobj, search= search->nextobj );
    if( !search )
      return;
    *link= search->nextobj;
  }
}


void printresyncreg()
{
  struct resync *idscan, *objscan;

  printf("resync object registry:\n");
  for( idscan= resyncreg; idscan; idscan= idscan->nextid ) {
    printf("id %d: ", idscan->id);
    for( objscan= idscan; objscan; objscan= objscan->nextobj )
      printf("`%s' (mode %d),  ", objscan->obj->name, objscan->mode);
    printf("\n");
  }
  printf("end of resync object registry\n");
}



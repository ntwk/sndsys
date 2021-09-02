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
#include <ctype.h>
#include <float.h>

#include "sndsys.h"


/*=============================================================================
    avg <freq> <fade> <signal>

Moving average filter which averages signal over period of <freq>.  When
applied to a step between two constant values, a ramp with the length of the
period will result.  <fade> determines mix between original and averaged
signal; 0 means averaged only, 1 original only, 2 means original minus averaged
(for highpass).

See also: gaussavg, rms, lowpass, shelve, highpass
=============================================================================*/

struct average {	/* also used by lee object - copy changes there */
  double    addsum[PLENTY], subsum[PLENTY];
  int       nch, maxind, subcountm1, nonorm, firstcall;
  float     freq;
};
struct avg {
  char              namestr[NAMESTRLEN];
  struct average    *a;
};

float calcavg(sndobj *me, sndobj *, int);
float average( sndobj *me, int innr, struct average *a, float newfreq, float *dest );

sndobj *avg( sndobj *freq, sndobj *fade, sndobj *signal )
{
  sndobj *p;
  struct avg *d;
  
  d= new(struct avg);
  snprintf( d->namestr, NAMESTRLEN, "[%s] Hz ~ [%s], [%s] %%", freq->name, signal->name, fade->name );
  p= newsndo( calcavg, d->namestr, "avg", signal->nch, 3, freq, fade, signal );
  p->skip= dontskip;
  p->private[0]= d;
  p->private[1]= d->a= new(struct average);
  d->a->nch= p->nch;
  d->a->nonorm= 0;
  d->a->firstcall= 1;
  buf( HORIZON/2.0, HORIZON/2.0 + 1.0/(double)SAMPRATE, 0.0, p->nch, p, 2 );
  return p;
}

float calcavg(sndobj *me, sndobj *caller, int innr )
{
  struct avg *d;
  float fade, origval;
  int ch;
  
  d= (struct avg *)me->private[0];
  average( me, 2, d->a, INCALC(0), me->ch );
  if( (fade= INCALC(1))!=0 ) {
    BUFINCALC(2, 0, 0.0);
    if( fade < 1 )
      for( ch= 0; ch< me->nch; ++ch ) {
	origval= INPUTC(2, ch);
	me->ch[ch] += fade*(origval-me->ch[ch]);
      }
    else {
      fade -= 1.0;
      for( ch= 0; ch< me->nch; ++ch ) {
	origval= INPUTC(2, ch);
	me->ch[ch]= origval - fade*me->ch[ch];
      }
    }
  }
  BUFINCR( 2, 1 );
  CALCRETURN;
}

// gnuplot command for frequency response (f is cutoff frequency):
// plot [0:10][0:1] f(x,f)= f * 1/(2*pi*exp(x)) * sqrt((cos(2*pi*exp(x)/f)-1)**2+sin(2*pi*exp(x)/f)**2), f(x,50)
float average( sndobj *me, int innr, struct average *a, float newfreq, float *dest )
{
  float retval;
  int i, ch, newmaxind;

  if( a->firstcall ) {
    if( !finite(newfreq) )
      return 0.0;
    if( newfreq <= 1.0/(float)HORIZON )
      newfreq= 1.0/(float)HORIZON;
    a->maxind= (int)round(0.5*(double)SAMPRATE/newfreq);
    if( !a->maxind )
      a->maxind= 1;     // otherwise subcountm1 becomes negative
    a->freq= newfreq;
    BUFINCALC(innr, 0, 0.0);
    for( ch= 0; ch< a->nch; ++ch ) {
      a->addsum[ch]= INPUTC(innr, ch);
      a->subsum[ch]= 0.0;
    }
    for( i= a->maxind; i> 0; --i ) {
      BUFINCALC(innr, i, 0.0);
      for( ch= 0; ch< a->nch; ++ch )
	a->addsum[ch] += INPUTC(innr, ch);
      BUFINCALC(innr, -i, 0.0);
      for( ch= 0; ch< a->nch; ++ch )
	a->subsum[ch] += INPUTC(innr, ch);
    }
    a->subcountm1= a->maxind - 1;
		    // number of values contained in subsum minus 1
    a->firstcall= 0;
  }
  else {
    if( newfreq!=a->freq && finite(newfreq) ) {
      if( newfreq <= 1.0/(float)HORIZON )
	newfreq= 1.0/HORIZON;
      newmaxind= (int)round(0.5*(double)SAMPRATE/newfreq);
      a->freq= newfreq;
    }
    else
      newmaxind= a->maxind;

    if( newmaxind > a->maxind ) {
      BUFINCALC(innr, a->maxind, 0.0);
      for( ch= 0; ch < a->nch; ++ch )
	a->addsum[ch] += INPUTC(innr, ch);
      BUFINCALC(innr, a->maxind+1, 0.0);
      for( ch= 0; ch < a->nch; ++ch )
	a->addsum[ch] += INPUTC(innr, ch);
      ++a->maxind;
    }
    else if( newmaxind < a->maxind ) {
      if( a->subcountm1<2 ) {
	if( !a->subcountm1 )
	  BUFINCALC(innr, -a->maxind, 0.0);
	for( ch= 0; ch < a->nch; ++ch ) {
	  a->subsum[ch]= a->addsum[ch];
	  a->addsum[ch]= 0.0;
	  if( !a->subcountm1 )
	    a->subsum[ch] -= INPUTC(innr, ch);
	}
	a->subcountm1= 2*a->maxind - 2;
      }
      else {
	BUFINCALC(innr, -a->maxind-1, 0.0);
	for( ch= 0; ch < a->nch; ++ch )
	  a->subsum[ch] -= INPUTC(innr, ch);
	BUFINCALC(innr, -a->maxind, 0.0);
	for( ch= 0; ch < a->nch; ++ch )
	  a->subsum[ch] -= INPUTC(innr, ch);
	a->subcountm1 -= 2;
      }
      --a->maxind;
    }
    else {
      if( !a->subcountm1 ) {
	BUFINCALC(innr, a->maxind, 0.0);
	for( ch= 0; ch < a->nch; ++ch ) {
	  a->subsum[ch]= a->addsum[ch] + INPUTC(innr, ch);
	  a->addsum[ch]= 0.0;
	}
	a->subcountm1= 2*a->maxind;
      }
      else {
	BUFINCALC(innr, -a->maxind-1, 0.0);
	for( ch= 0; ch < a->nch; ++ch )
	  a->subsum[ch] -= INPUTC(innr, ch);
	BUFINCALC(innr, a->maxind, 0.0);
	for( ch= 0; ch < a->nch; ++ch )
	  a->addsum[ch] += INPUTC(innr, ch);
	--a->subcountm1;
      }
    }
  }
  if( !a->nonorm ) {
    retval= (a->subsum[0]+a->addsum[0])/(2*a->maxind+1);
    if( dest ) {
      dest[0]= retval;
	for( ch= 1; ch < a->nch; ++ch )
	dest[ch] = (a->subsum[ch]+a->addsum[ch])/(2*a->maxind+1);
    }
  }
  else {
    retval= a->subsum[0]+a->addsum[0];
    if( dest ) {
      dest[0]= retval;
	for( ch= 1; ch < a->nch; ++ch )
	dest[ch] = a->subsum[ch]+a->addsum[ch];
    }
  }
  return retval;
}


/*=============================================================================
    gaussavg (quality) <freq> <fade> <signal>

Convolves signal with a gaussian with a mean deviation of 1/2 period of <freq>.
(quality) is the number of moving averages used to model the gaussian (5 to
10 are useful values).  <fade> determines mix between original and averaged
signal; 0 means averaged only, 1 original only, 2 means original minus averaged
(for highpass).

See also: avg, lowpass, shelve
=============================================================================*/

struct gaussavg {
  int               quality;
  float             *freqmul;
  struct average    *a;
};

float calcgaussavg(sndobj *me, sndobj *, int);

#define GAUSSAVG_MAXQUAL    100

sndobj *gaussavg( int quality, sndobj *freq, sndobj *fade, sndobj *signal )
{
  sndobj *p;
  struct gaussavg *d;
  double gaussy, ystep;
  int i;
  
  p= newsndo( calcgaussavg, signal->name, "gaussavg", signal->nch, 3, freq, fade, signal );
  p->skip= dontskip;
  p->private[0]= d= new(struct gaussavg);
  if( quality > GAUSSAVG_MAXQUAL ) {
    printf( "Warning: gaussavg: Limiting quality (%d) to %d.\n", quality, (int)GAUSSAVG_MAXQUAL );
    quality= GAUSSAVG_MAXQUAL;
  }
  d->quality= quality;
  p->private[1]= d->freqmul= news(float, quality);
  p->private[2]= d->a= news(struct average, quality);
  gaussy= ystep= 1.0/(quality+1);
  for( i= 0; i< quality; ++i ) {
    d->a[i].nch= p->nch;
    d->a[i].nonorm= 1;
    d->a[i].firstcall= 1;
    d->freqmul[i]= 1.0/sqrt(-2.0*log(gaussy));
//printf("Gaussian freq multiplier # %d: %g\n", i, d->freqmul[i] );
    gaussy += ystep;
  }
  buf( HORIZON/2.0, HORIZON/2.0 + 1.0/(double)SAMPRATE, 0.0, p->nch, p, 2 );
  return p;
}

float calcgaussavg(sndobj *me, sndobj *caller, int innr )
{
  struct gaussavg *d;
  float avgresult[PLENTY];
  float freq, fade, origval;
  int norm, i, ch;
  
  d= (struct gaussavg *)me->private[0];
  freq= INCALC(0);
  FORCH
    me->ch[ch]= 0.0;
  norm= 0;
  for( i= 0; i< d->quality; ++i ) {
    average( me, 2, d->a+i, freq*d->freqmul[i], avgresult );
    FORCH
      me->ch[ch] += avgresult[ch];
    norm += 2*d->a[i].maxind+1;
  }
  FORCH
    me->ch[ch] /= norm;
  if( (fade= INCALC(1))!=0 ) {
    BUFINCALC(2, 0, 0.0);
    if( fade < 1 )
      for( ch= 0; ch< me->nch; ++ch ) {
	origval= INPUTC(2, ch);
	me->ch[ch] += fade*(origval-me->ch[ch]);
      }
    else {
      fade -= 1.0;
      for( ch= 0; ch< me->nch; ++ch ) {
	origval= INPUTC(2, ch);
	me->ch[ch]= origval - fade*me->ch[ch];
      }
    }
  }
  BUFINCR( 2, 1 );
  CALCRETURN;
}


/*=============================================================================
    rms <freq> <signal>

Root of mean square, for all channels separately.  Applies a moving average
filter to the square of the signal before taking the square root.  The cutoff
frequency should change smoothly; a change takes as long as half the difference
in wave periods to take effect.

See also: fastrms, avg, gaussavg
=============================================================================*/
/*
Note: Faster reaction to changes in <freq> could be achieved by using the same
procedure as on the first call of the calc function, but that would be
computationally expensive.
*/

struct rms {
  char      namestr[NAMESTRLEN];
  double    addsum[PLENTY], subsum[PLENTY];
  int       maxind, subcountm1, firstcall;
};

float calcrms(sndobj *me, sndobj *, int);

sndobj *rms( sndobj *freq, sndobj *signal )
{
  sndobj *p;
  struct rms *d;
  
  d= (struct rms *)malloc(sizeof(struct rms));
  snprintf(d->namestr, NAMESTRLEN, "[%s] Hz ~^2 [%s]", freq->name, signal->name );
  p= newsndo( calcrms, d->namestr, "rms", signal->nch, 2, freq, signal );
  p->skip= dontskip;
  p->private[0]= d;
  d->firstcall= 1;
  buf( HORIZON/2.0, HORIZON/2.0 + 1.0/(double)SAMPRATE, 0.0, p->nch, p, 1 );
  return p;
}


float calcrms(sndobj *me, sndobj *caller, int innr )
{
  struct rms *d;
  double val, sqravg;
  float freq;
  int newmaxind, ch, i;

  d= (struct rms *)me->private[0];
  freq= INCALC(0);
  if( freq <= 1.0/(float)HORIZON )
    freq= 1.0/HORIZON;
  newmaxind= (int)round(0.5*(double)SAMPRATE/freq);
  if( d->firstcall ) {
    d->maxind= newmaxind;
    BUFINCALC(1, 0, 0.0);
    for( ch= 0; ch< me->nch; ++ch ) {
      val= INPUTC(1, ch);
      d->addsum[ch]= val*val;
      d->subsum[ch]= 0.0;
    }
    for( i= d->maxind; i> 0; --i ) {
      BUFINCALC(1, i, 0.0);
      for( ch= 0; ch< me->nch; ++ch ) {
	val= INPUTC(1, ch);
	d->addsum[ch] += val*val;
      }
      BUFINCALC(1, -i, 0.0);
      for( ch= 0; ch< me->nch; ++ch ) {
	val= INPUTC(1, ch);
	d->subsum[ch] += val*val;
      }
    }
    d->subcountm1= d->maxind - 1;
		    // number of values contained in subsum minus 1
    d->firstcall= 0;
  }
  else if( newmaxind > d->maxind ) {
    BUFINCALC(1, d->maxind, 0.0);
    for( ch= 0; ch < me->nch; ++ch ) {
      val= INPUTC(1, ch);
      d->addsum[ch] += val*val;
    }
    BUFINCALC(1, d->maxind+1, 0.0);
    for( ch= 0; ch < me->nch; ++ch ) {
      val= INPUTC(1, ch);
      d->addsum[ch] += val*val;
    }
    ++d->maxind;
  }
  else if( newmaxind < d->maxind ) {
    if( d->subcountm1<2 ) {
      if( !d->subcountm1 )
	BUFINCALC(1, -d->maxind, 0.0);
      for( ch= 0; ch < me->nch; ++ch ) {
	d->subsum[ch]= d->addsum[ch];
	d->addsum[ch]= 0.0;
	if( !d->subcountm1 ) {
	  val= INPUTC(1, ch);
	  d->subsum[ch] -= val*val;
	}
      }
      d->subcountm1= 2*d->maxind - 2;
    }
    else {
      BUFINCALC(1, -d->maxind-1, 0.0 );
      for( ch= 0; ch < me->nch; ++ch ) {
	val= INPUTC(1, ch);
	d->subsum[ch] -= val*val;
      }
      BUFINCALC(1, -d->maxind, 0.0 );
      for( ch= 0; ch < me->nch; ++ch ) {
	val= INPUTC(1, ch);
	d->subsum[ch] -= val*val;
      }
      d->subcountm1 -= 2;
    }
    --d->maxind;
  }
  else {
    if( !d->subcountm1 ) {
      BUFINCALC(1, d->maxind, 0.0 );
      for( ch= 0; ch < me->nch; ++ch ) {
	val= INPUTC(1, ch);
	d->subsum[ch]= d->addsum[ch] + val*val;
	d->addsum[ch]= 0.0;
      }
      d->subcountm1= 2*d->maxind;
    }
    else {
      BUFINCALC(1, -d->maxind-1, 0.0 );
      for( ch= 0; ch < me->nch; ++ch ) {
	val= INPUTC(1, ch);
	d->subsum[ch] -= val*val;
      }
      BUFINCALC(1, d->maxind, 0.0 );
      for( ch= 0; ch < me->nch; ++ch ) {
	val= INPUTC(1, ch);
	d->addsum[ch] += val*val;
      }
      --d->subcountm1;
    }
  }
  BUFINCR(1, 1);
  ch= 0;
  for( ch= 0; ch < me->nch; ++ch ) {
    sqravg= (d->subsum[ch]+d->addsum[ch])/(2*d->maxind+1);
    if( sqravg < 0 )   // possible due to rounding errors
      sqravg= 0.0;
    OUTPUT(ch, sqrt( sqravg ) );
  }
  CALCRETURN;
}


/*=============================================================================
    bands [tab] <fmult> <signal>

Split signal into frequency bands.  <tab> contains boundary frequencies from
top to bottom, followed by END.  <fmult> is a variable frequency multiplier.
The output channels contain the bands from top to bottom in frequency.  All
bands from the first signal channel precede all from the next.  Not all bands
for all channels are output if the total exceeds PLENTY.  Successive freqencies
should differ by at least a factor of 2.  Separating the bands is done by
moving average filters with a width corresponding to a period of the boundary
frequencies.

See also: graphiceq, avg
=============================================================================*/

float canonicalbands[]= { 7000., 5000., 2500., 300., 60., END };

struct bands {
  float     *tab;
  struct average    *avgs;
  int       avgnum, innch;
};

float calcbands(sndobj *me, sndobj *, int);

sndobj *bands( float *tab, sndobj *fmult, sndobj *signal )
{
  sndobj *p;
  struct bands *d;
  float *read;
  int i;
  
  p= newsndo( calcbands, signal->name, "bands", 1, 2, fmult, signal );
  p->private[0]= d= new(struct bands);
  buf( HORIZON/2.0, HORIZON/2.0 + 1.0/SAMPRATE, 0.0, signal->nch, p, 1 );
  for( read= tab, d->avgnum= 0; *read!=END; ++read )
    ++d->avgnum;
  if( d->avgnum>PLENTY ) {
    fprintf(stderr, "Warning: bands: too many bands (%d, max. %d).\n", (int)d->avgnum, (int)PLENTY );
    d->avgnum= PLENTY;
  }
  d->tab= tab;
  d->innch= signal->nch;
  p->nch= (d->avgnum+1)*d->innch;
  if( p->nch>PLENTY )
    p->nch= PLENTY;
  if( !d->avgnum )
    d->avgs= NULL;
  else {
    p->private[1]= d->avgs= news(struct average, d->avgnum);
    for( i= 0; i< d->avgnum; ++i ) {
      d->avgs[i].nch= d->innch;
      d->avgs[i].nonorm= 0;
      d->avgs[i].firstcall= 1;
    }
  }
  return p;
}

float calcbands(sndobj *me, sndobj *caller, int innr )
{
  struct bands *d;
  float avgvalues[PLENTY], prevvalues[PLENTY];
  float fmult;
  int ch, outch, band;
  
  d= (struct bands *)me->private[0];
  fmult= INCALC(0);
  BUFINCALC(1, 0, 0.0);
  for( ch= 0; ch< d->innch; ++ch )
    prevvalues[ch]= INPUTC(1, ch);
  for( band= 0; band< d->avgnum; ++band )
  {
    average( me, 1,  d->avgs+band, d->tab[band]*fmult, avgvalues );
    for( outch= band, ch= 0; outch< me->nch; outch += d->avgnum+1, ++ch )
      me->ch[outch]= prevvalues[ch] - avgvalues[ch];
    for( ch= 0; ch< d->innch; ++ch )
      prevvalues[ch]= avgvalues[ch];
  }
  for( outch= d->avgnum, ch= 0; outch< me->nch; outch += d->avgnum+1, ++ch )
    me->ch[outch]= prevvalues[ch];
  BUFINCR(1, 1);
  CALCRETURN;
}


/*=============================================================================
    graphiceq [tab] <fmult> <signal>

Graphic equaliser (equalizer to Americans ;)).  <tab> must contain floats,
which are alternatingly interpreted as gain/attenuation factors and frequencies
limiting the bands.  The first value is a gain factor, the frequencies must be
given from top to bottom, and the end is signalled by the magic value END after
the last gain value.  Successive freqencies should be at least a factor 2
apart.  The second argument is a variable multiplier for the frequencies.
The frequency bands are separated by a bank of moving average filters, as in
the bands object.

See also: bands, avg
=============================================================================*/

struct graphiceq {
  float     *tab;
  struct average    *avgs;
  int       avgnum;
};

float calcgraphiceq(sndobj *me, sndobj *, int);

sndobj *graphiceq( float *tab, sndobj *fmult, sndobj *signal )
{
  sndobj *p;
  struct graphiceq *d;
  float *read;
  int tabvals, i;
  
  p= newsndo( calcgraphiceq, signal->name, "graphiceq", signal->nch, 2, fmult, signal );
  p->private[0]= d= new(struct graphiceq);
  buf( HORIZON/2.0, HORIZON/2.0 + 1.0/(double)SAMPRATE, 0.0, signal->nch, p, 1 );
  for( read= tab, tabvals= 0; *read!=END; ++read )
    ++tabvals;
  d->tab= tab;
  d->avgnum= tabvals/2;
  if( !d->avgnum )
    d->avgs= NULL;
  else {
    p->private[1]= d->avgs= news(struct average, d->avgnum);
    for( i= 0; i< d->avgnum; ++i ) {
      d->avgs[i].nch= p->nch;
      d->avgs[i].nonorm= 0;
      d->avgs[i].firstcall= 1;
    }
  }
  return p;
}

float calcgraphiceq(sndobj *me, sndobj *caller, int innr )
{
  struct graphiceq *d;
  float avgvalues[PLENTY], prevvalues[PLENTY];
  float fmult;
  int ch, i;
  
  d= (struct graphiceq *)me->private[0];
  fmult= INCALC(0);
  BUFINCALC(1, 0, 0.0);
  for( ch= 0; ch< me->nch; ++ch ) {
    prevvalues[ch]= INPUTC(1, ch);
    me->ch[ch]= 0.0;
  }
  for( i= 0; i< d->avgnum; ++i )
  {
    average( me, 1, d->avgs+i, d->tab[2*i+1]*fmult, avgvalues );
    for( ch= 0; ch< me->nch; ++ch ) {
      me->ch[ch] += (prevvalues[ch]-avgvalues[ch])*d->tab[2*i];
      prevvalues[ch]= avgvalues[ch];
    }
  }
  if( d->tab[2*d->avgnum]!=END )
    for( ch= 0; ch< me->nch; ++ch )
      me->ch[ch] += prevvalues[ch]*d->tab[2*d->avgnum];
  else
    for( ch= 0; ch< me->nch; ++ch )
      me->ch[ch] += prevvalues[ch];
  BUFINCR(1, 1);
  CALCRETURN;
}


/*=============================================================================
    filter [forwardfeed] [backfeed] <signal>

The tables [forwardfeed] and [backfeed] contain delays and coefficients
alternatingly.  The delays are given in sample values, not seconds and should
be integers (they will be truncated).  The forwardfeed delays may be negative
by up to one second.  The end of the tables is signalled by END.  [backfeed]
may be NULL.  Coefficients but not delays may be changed between calls of
calcfilter.  The feedback coefficients are not negated, as is often done (the
number added to the output is +coeff.*y(n-delay)).

See also: freqfir, convolve, lowpass, highpass, shelve and pretty much anything
else in sndfilt.c
=============================================================================*/

struct filter {
  int       nforw, nback;
  float     *forwco, *backco;
  int       *forwdel, *backdel;
};

float calcfilter(sndobj *me, sndobj *, int);

sndobj *filter( float *forwardfeed, float *backfeed, sndobj *signal )
{
  sndobj *p;
  struct filter *d;
  float *read;
  float maxxdelay, minxdelay, maxydelay;
  int i;
  
  if( !forwardfeed || *forwardfeed==END )
    return c(0.0);
  p= newsndo( calcfilter, signal->name, "filter", signal->nch, 1, signal );
  p->private[0]= d= new(struct filter);
  maxxdelay= -FLT_MAX;
  minxdelay= FLT_MAX;
  for( read= forwardfeed; *read!=END; read+=2 ) {
    if( *read>maxxdelay )   maxxdelay= *read;
    if( *read<minxdelay )   minxdelay= *read;
  }
  if( minxdelay < -HORIZON*SAMPRATE )
    minxdelay= -HORIZON*SAMPRATE;
  else if( minxdelay > HORIZON*SAMPRATE )
    minxdelay= HORIZON*SAMPRATE;
  d->nforw= (read-forwardfeed)/2;
  d->forwco= forwardfeed + 1;
  p->private[1]= d->forwdel= news(int, d->nforw);
  for( i= 0, read= forwardfeed; i< d->nforw; ++i, read+=2 )
    d->forwdel[i]= (int)*read;
  buf( maxxdelay/SAMPRATE, -minxdelay/SAMPRATE, 0.0, p->nch, p, 0 );
  if( !backfeed || *backfeed==END )
    d->nback= 0;
  else {
    maxydelay= 0.0;
    for( read= backfeed; *read!=END; read += 2 )
      if( *read > maxydelay )
	maxydelay= *read;
    d->nback= (read-backfeed)/2;
    d->backco= backfeed+1;
    p->private[2]= d->backdel= news(int, d->nback);
    for( i= 0, read= backfeed; i< d->nback; ++i, read+=2 )
      if( *read <= 0 ) {
	fprintf(stderr, "Warning: filter: backfeed delay %ld is <=0 (%g) and will be set to 1.\n", (long)(read-backfeed)/2, *read );
	d->backdel[i]= 1;
      }
      else d->backdel[i]= (int)*read;
    addinput(p, p);
    buf( maxydelay/SAMPRATE, 0.0, 0.0, p->nch, p, 1 );
  }
  return p;
}

float calcfilter(sndobj *me, sndobj *caller, int innr )
{
  struct filter *d;
  float result[PLENTY];
  int i, ch;
  
  d= (struct filter *)me->private[0];
  for( ch= 0; ch< me->nch; ++ch )
    result[ch]= 0.0;
  for( i= 0; i< d->nforw; ++i ) {
    BUFINCALC(0, -d->forwdel[i], 0.0);
    for( ch= 0; ch< me->nch; ++ch )
      result[ch] += INPUTC(0, ch) * d->forwco[2*i];
  }
  for( i= 0; i< d->nback; ++i ) {
    BUFINCALC(1, -d->backdel[i], 0.0);
    for( ch= 0; ch< me->nch; ++ch )
      result[ch] += INPUTC(1, ch) * d->backco[2*i];
  }
  for( ch= 0; ch< me->nch; ++ch )
    OUTPUT(ch, result[ch]);
  BUFINCR(0, 1);
  if( d->nback )
    BUFINCR(1, 1);
  CALCRETURN;
}


/*=============================================================================
    freqfir [response] (order) <signal>

FIR filter with given frequency response.  [response] has to contain
frequency-factor (0..1) pairs which are to be linearly interpolated.  If the
first and last points are not at zero and half the sampling rate, respectively,
the response outside the given range is taken to be constant.  (order)
determines the quality of the result.  It has to be even and should be at least
20 for a smallish number of points in [response] (minimum is 10).  The
corresponding FIR filter coefficients will be determined and a filter object
returned.

See also: filter
=============================================================================*/

sndobj *freqfir(float *response, int order, sndobj *signal)
{
  sndobj *p;
  float *sampled, *read, *forwfeed, *write;
  float f, fstep, coeff, fhalforderm1;
  int i, j, halforder;
  
  if( order<10 )
    order= 10;
  if( (order&1)!=0 )
    ++order;
  // First sample the response envelope at equidistant points
  for( read= response; *read!=END; read += 2 );
  if( read-response < 4 ) {
    if( read==response )
      return signal;
    else
      return linear(0.0, read[1], signal);
  }
  halforder= order/2;
  sampled= news(float, halforder);
  fstep= (float)SAMPRATE/(order-1);
  read= response;
  write= sampled;
  for( i= 0, f= 0.0; f< read[0] &&  i< halforder; ++i, f+=fstep )
    *write++= read[1];
  for( ; i< halforder; ++i, f+=fstep ) {
    while( read[2]!=END && read[2]< f )
      read += 2;
    if( read[2]==END )
      *write++ = read[1];
    else
      *write++ = read[1] + (read[3]-read[1])/(read[2]-read[0])*(f-read[0]);
  }

  // Now "Fourier transform" into FIR filter forwardfeed coefficients
  forwfeed= news(float, 2*order+1);
  fhalforderm1= 0.5 * (order - 1);
  f= 2.0*M_PI / (float)order;
  for( j= 0; j < halforder; ++j )
  {
    coeff= sampled[0] * 0.5;
    for( i= 1; i < halforder; ++i )
      coeff += sampled[i] * cos(f * (fhalforderm1 - j) * i);
    // indices -> delay 1/2 sample independent of order
    forwfeed[2*j]= j-halforder+1;
    forwfeed[2*(order-1-j)]= halforder-j;
    forwfeed[2*(order-1-j)+1]= forwfeed[2*j+1]= 2.0 * coeff / (float)order;
  }
  forwfeed[2*order]= END;
  free(sampled);
  
  // Last: return filter object
  p= filter( forwfeed, NULL, signal );
  for( i= 0; i< PLENTY && p->private[i]; ++i );
  if( i<PLENTY )
    p->private[i]= forwfeed;
  p->type= "freqfir";
  return p;
}


/*=============================================================================
    formant <frequency> <bandwidth> <signal>

Band-pass filter.
Difference equation:
y(n) = x(n) - r*x(n-2) + 2*r*cos(2*pi*frequency/srate)*y(n-1) - r*r*y(n-2)
r = exp( -pi * bandwidth / srate )
Resonance gain for this equation: 0.5201306 + 14033.69386 / bw [Hz]
The result is normalised by dividing by resonance gain, but this can be
numerically imprecise for narrow resonances (small bandwidth).

See also: peak, peaknotch
=============================================================================*/

struct formant {
  double    xcoeffs[2], ycoeffs[2];
  float     currfreq, currbw;
  double    gain, radius;
  double    x[2*PLENTY], y[2*PLENTY];
};

float calcformant(sndobj *me, sndobj *caller, int innr);

sndobj *formant( sndobj *freq, sndobj *bandwidth, sndobj *signal )
{
  sndobj *p;
  struct formant *d;
  int ch;

  p= newsndo( calcformant, "formant", "formant", signal->nch, 3, signal, freq, bandwidth );
  p->private[0]= d= new(struct formant);
  d->currfreq= 0.0;
  d->currbw= 0.0;
  d->gain= 1.0;
  for( ch= 0; ch< p->nch; ++ch )
    d->x[ch]= d->x[ch+PLENTY]= d->y[ch]= d->y[ch+PLENTY]= 0.0;
  d->xcoeffs[0]= d->xcoeffs[1]= 1.0;
  d->ycoeffs[0]= 2.0;
  d->ycoeffs[1]= 1.0;
  return p;
}

float calcformant(sndobj *me, sndobj *caller, int innr)
{
  struct formant *d;
  float newfreq, newbw;
  double phase, y;
  int ch;

  d= (struct formant *)me->private[0];
  newfreq= INCALC(1);
  newbw= INCALC(2);
  if( newbw != d->currbw )
  {
    if( newbw )
      d->xcoeffs[0]= 1.0 / (0.5201306 + 14033.69386 / (double)newbw);
    else
      d->xcoeffs[0]= 0.0;
    d->radius= exp( - (double)newbw * M_PI / SAMPRATE );
    phase= (double)newfreq*2.0*M_PI/SAMPRATE;
    d->xcoeffs[1]= -d->radius * d->xcoeffs[0];
    d->ycoeffs[0]= 2.0*d->radius*cos(phase);
    d->ycoeffs[1]= -d->radius*d->radius;
    d->currfreq= newfreq;
    d->currbw= newbw;
//    fprintf(stderr, "new bandwidth %g, freq %g, coeffs %g, %g, %g\n", 
//	  newbw, newfreq, d->coeffs[3], d->coeffs[5+1], d->coeffs[5+3]);
  }
  else if( newfreq!=d->currfreq )
  {
    phase= (double)newfreq*2.0*M_PI/SAMPRATE;
    d->ycoeffs[0]= 2.0*d->radius*cos(phase);
    d->currfreq= newfreq;
  }
  INCALC(0);
  FORCH {
    y= d->xcoeffs[0]*INPUTC(0, ch) + d->xcoeffs[1]*d->x[ch+PLENTY] +
      d->ycoeffs[0]*d->y[ch] + d->ycoeffs[1]*d->y[ch+PLENTY];
    OUTPUT(ch, y);
    d->x[ch+PLENTY]= d->x[ch];
    d->x[ch]= INPUTC(0, ch);
    d->y[ch+PLENTY]= d->y[ch];
    d->y[ch]= y;
  }
  CALCRETURN;
}


/*=============================================================================
    shelve <freq> <ratio> <signal>

Shelve filter.  <ratio> is the factor in amplitude between very high and very
low frequencies.  For <ratio> > 1, this object is a shelve highpass, for
<ratio> < 1, a lowpass.  <freq> is the cutoff frequency.

See also: lowpass, highpass
=============================================================================*/

struct shelve {
  float     coeffs[12];
  float     freq, ratio;
};

float calcshelve(sndobj *me, sndobj *caller, int innr);

sndobj *shelve( sndobj *freq, sndobj *ratio, sndobj *signal )
{
  sndobj *p;
  struct shelve *d;
  float *forwfeed, *backfeed;
  
  d= new(struct shelve);
  forwfeed= d->coeffs;
  forwfeed[0]= 0.0;
  forwfeed[1]= 1.0;
  forwfeed[2]= 1.0;
  forwfeed[4]= 2.0;
  forwfeed[3]= forwfeed[5]= 0.0;
  forwfeed[6]= END;
  backfeed= d->coeffs+7;
  backfeed[0]= 1.0;
  backfeed[2]= 2.0;
  backfeed[1]= backfeed[3]= 0.0;
  backfeed[4]= END;
  d->freq= MAGIC;
  d->ratio= MAGIC;
  p= filter( forwfeed, backfeed, signal );
  p->calc= calcshelve;
  addinput(p, freq);
  addinput(p, ratio);
  p->private[PLENTY-1]= d;
  p->type= "shelve";
  return p;
}

float calcshelve(sndobj *me, sndobj *caller, int innr)
{
  struct shelve *d;
  double newratio;
  float newfreq;
  
  d= (struct shelve *)me->private[PLENTY-1];
  newfreq= INCALC(2);   // first two inputs are forwardfeed and backfeed
  if( newfreq < 0 )
    newfreq= -newfreq;
  if( newfreq > SAMPRATE )
    newfreq= SAMPRATE;
  newratio= INCALC(3);
  if( newratio < 0 )
    newratio= -newratio;
  if( newfreq != d->freq || (float)newratio != d->ratio )
  {
    double tanf, tanf2, F, F2, tmp, gamman, gammad, gamma2, siggam2, gam2p1;
    double ta0, ta1, ta2, tb0, tb1, tb2, aa1, ab1, b0, scale;

    d->freq= newfreq;
    d->ratio= newratio;
    
    // *** invert ratio here - high-f-to-low-f -> low-f-to-high-f ***
    if( newratio )
      newratio= 1.0/newratio;
    else
      newratio= DBL_MAX;
    
    tanf = tan(M_PI*(newfreq/SAMPRATE-0.25));
    tanf2 = tanf*tanf;
    if( newratio > 0.5 && newratio < 2.0 )
      F = sqrt(newratio);
    else if( newratio > 1.0 )
      F = newratio / M_SQRT2;
    else
      F = newratio * M_SQRT2;
  
    F2 = F * F;
    tmp = newratio*newratio - F2;
    if( fabs(tmp) <= DBL_MIN ) gammad = 1.0;
    else gammad = pow( (F2 - 1.0)/tmp, 0.25);
    gamman = sqrt(newratio) * gammad;
  
    gamma2 = gamman * gamman;
    gam2p1 = 1.0 + gamma2;
    siggam2 = M_SQRT2 * gamman;
    ta0 = gam2p1 + siggam2;
    ta1 = -2.0 * (1.0 - gamma2);
    ta2 = gam2p1 - siggam2;
  
    aa1 = tanf * ta1;
    d->coeffs[1] = ta0 + aa1 + tanf2 * ta2;
    d->coeffs[3] = 2.0 * tanf * (ta0 + ta2) + (1.0 + tanf2) * ta1;
    d->coeffs[5] = tanf2 * ta0 + aa1 + ta2;
  
    gamma2 = gammad * gammad;
    gam2p1 = 1.0 + gamma2;
    siggam2 = M_SQRT2 * gammad;
    tb0 = gam2p1 + siggam2;
    tb1 = -2.0 * (1.0 - gamma2);
    tb2 = gam2p1 - siggam2;
  
    ab1 = tanf * tb1;
    b0 = tb0 + ab1 + tanf2 * tb2;
    d->coeffs[8] = - (2.0 * tanf * (tb0 + tb2) + (1.0 + tanf2) * tb1);
    d->coeffs[10] = - (tanf2 * tb0 + ab1 + tb2);
  
    scale = 1.0/b0;
    d->coeffs[8] *= scale;
    d->coeffs[10] *= scale;
    if( newratio > 1.0 )
      scale /= newratio;
    d->coeffs[1] *= scale;
    d->coeffs[3] *= scale;
    d->coeffs[5] *= scale;
  }
  return calcfilter(me, caller, innr);
}


/*=============================================================================
    _allpass (delay) (coeff) <signal>

Schroeder-like allpass.
y(n) = (coeff)*x(n-1) + x(n-delay-1) - (coeff)*y(n-delay)

This sndobj's name starts with an underscore to indicate it is primarily for
internal use.  It does not have inputs; its parameters remain constant.

See also: _foallpass, filter
=============================================================================*/

sndobj *_allpass( int delay, float coeff, sndobj *signal )
{
  sndobj *p;
  float *coeffs, *forw, *back;
  
  if( delay<=0 )
    delay= 1;
  coeffs= news(float, 8);
  forw= coeffs;
  forw[0]= 1.0;
  forw[1]= coeff;
  forw[2]= delay+1;
  forw[3]= 1.0;
  forw[4]= END;
  back= coeffs+5;
  back[0]= delay;
  back[1]= -coeff;
  back[2]= END;
  p= filter( forw, back, signal );
  p->type= "_allpass";
  p->private[PLENTY-1]= coeffs;
  return p;
}

/*=============================================================================
    _foallpass (coeff) <signal>

First-order allpass filter.
y(n) = (coeff)*x(n) - x(n-1) + (coeff)*y(n-1)

This sndobj's name starts with an underscore to indicate it is primarily for
internal use.  It does not have inputs; its parameters remain constant.

See also: _allpass, filter
=============================================================================*/

sndobj *_foallpass( float coeff, sndobj *signal )
{
  sndobj *p;
  float *coeffs, *forw, *back;
  
  coeffs= news(float, 8);
  forw= coeffs;
  forw[0]= 0.0;
  forw[1]= coeff;
  forw[2]= 1.0;
  forw[3]= -1.0;
  forw[4]= END;
  back= coeffs+5;
  back[0]= 1.0;
  back[1]= coeff;
  back[2]= END;
  p= filter( forw, back, signal );
  p->type= "_foallpass";
  p->private[PLENTY-1]= coeffs;
  return p;
}

/*=============================================================================
    _comb (delay) (backfeed) <signal>

Feedback comb filter of order (delay), with backfeed coefficient (backfeed).  A
comb filter amplifies frequencies which are equally spaced in frequency, at a
distance of SAMPRATE/(delay).  The height (or depth) of the peaks is determined
by (backfeed): The peak gain is 1/(1-|(backfeed)|).  Its difference equation is:
y(n)= x(n) + (backfeed) y(n-(delay)) .

This sndobj's name starts with an underscore to indicate it is primarily for
internal use.  It does not have inputs; its parameters remain constant.

See also: _fcomb, filter; nrev
=============================================================================*/

sndobj *_comb( int delay, float backfeed, sndobj *signal )
{
  sndobj *p;
  float *coeffs, *forw, *back;

  if( delay<=0 )
    delay= 1;
  coeffs= news(float, 6);
  forw= coeffs;
  forw[0]= 0.0;
  forw[1]= 1.0;
  forw[2]= END;
  back= coeffs+3;
  back[0]= delay;
  back[1]= backfeed;
  back[2]= END;
  p= filter( forw, back, signal );
  p->type= "_comb";
  p->private[PLENTY-1]= coeffs;
  return p;
}

/*=============================================================================
    _fcomb (delay) (bf1) (bf2) <signal>

Comb filter with double backfeed.  Difference equation:
y(n)= x(n) + (bf1) * y(n-(delay)) + (bf2) * y(n-(delay)-1)

This sndobj's name starts with an underscore to indicate it is primarily for
internal use.  It does not have inputs; its parameters remain constant.

See also: _comb, filter
=============================================================================*/

sndobj *_fcomb( int delay, float bf1, float bf2, sndobj *signal )
{
  sndobj *p;
  float *coeffs, *forw, *back;
  
  if( delay<=0 )
    delay= 1;
  coeffs= news(float, 8);
  forw= coeffs;
  forw[0]= 0.0;
  forw[1]= 1.0;
  forw[2]= END;
  back= coeffs+3;
  back[0]= delay;
  back[1]= bf1;
  back[2]= delay+1;
  back[3]= bf2;
  back[4]= END;
  p= filter( forw, back, signal );
  p->type= "_fcomb";
  p->private[PLENTY-1]= coeffs;
  return p;
}

/*=============================================================================
    _12_0pole (iszero) (a0) (ab1) (ab2) <signal>

One-zero, Two-zero, one-pole or two-pole filter.  If <iszero> is !=0, a zero
filter is created, other wise a pole filter.  The other parameters are the
(constant) filter coefficients.  To obtain a one-zero/pole filter, set <ab2>=0.
<ab1> and (if applicable) <ab2> are either forwardfeed or backfeed coefficients
depending on the filter type.

The difference equations are:
one-zero  y(n) = a0 x(n) + a1 x(n-1)
two-zero  y(n) = a0 x(n) + a1 x(n-1) + a2 x(n-2)
one-pole  y(n) = a0 x(n) - b1 y(n-1)
two-pole  y(n) = a0 x(n) - b1 y(n-1) - b2 y(n-2)

This sndobj's name starts with an underscore to indicate it is primarily for
internal use.  It does not have inputs; its parameters remain constant.

See also: filter
=============================================================================*/

sndobj *_12_0pole( int iszero, float a0, float ab1, float ab2, sndobj *signal )
{
  sndobj *p;
  float *coeffs, *forw, *back;
  char *type;
  
  coeffs= news(float, 8);
  forw= coeffs;
  forw[0]= 0;
  forw[1]= a0;
  if( iszero ) {
    forw[2]= 1;
    forw[3]= ab1;
    if( ab2 ) {
      forw[4]= 2;
      forw[5]= ab2;
      forw[6]= END;
      type= "_two-zero";
    }
    else {
      forw[4]= END;
      type= "_one-zero";
    }
    back= NULL;
  }
  else {
    forw[2]= END;
    back= coeffs+3;
    back[0]= 1;
    back[1]= -ab1;
    if( ab2 ) {
      back[2]= 2;
      back[3]= -ab2;
      back[4]= END;
      type= "_two-pole";
    }
    else {
      back[2]= END;
      type= "_one-pole";
    }
  }
  p= filter( forw, back, signal );
  p->type= type;
  p->private[PLENTY-1]= coeffs;
  return p;
}


/*=============================================================================
    lowpass <cutoff> <signal>

Simple fast first-order IIR filter favouring low frequencies.  Difference
equation:
y(n)= c x(n) + (1-c) y(n-1)

See also: shelve, avg, gaussavg, highpass
=============================================================================*/

struct fofilter {
  float  cx, cx1, cy1;
  float  x1[PLENTY], y1[PLENTY];
  float  currparam;
};

float calclowpass(sndobj *me, sndobj *, int);

sndobj *lowpass( sndobj *cutoff, sndobj *signal )
{
  sndobj *p;
  struct fofilter *d;
  int i;
  
  p= newsndo(calclowpass, signal->name, "lowpass", signal->nch, 2, cutoff, signal );
  p->private[0]= d= new(struct fofilter);
  d->cx= d->cx1= d->cy1= 0.0;
  for( i= 0; i< p->nch; ++i )
    d->y1[i]= 0.0;
  d->currparam= MAGIC;
  return p;
}

float calclowpass(sndobj *me, sndobj *caller, int innr )
{
  struct fofilter *d;
  float newcut, val;
  int ch;
  
  d= (struct fofilter *)me->private[0];
  newcut= INCALC(0);
  if( newcut!=d->currparam ) {
    if( newcut>(float)SAMPRATE/2.0 )
      newcut= (float)SAMPRATE/2.0;
    val= 2.0 - cos(newcut*2.0*M_PI/(float)SAMPRATE);
    d->cx= 1.0 - (val - sqrt(val*val - 1.0));
    d->currparam= newcut;
  }
// In effect, cy1= 1-cx and cx1= 0
  d->y1[0]= (INCALC(1) - d->y1[0])*d->cx + d->y1[0]; 
  OUTPUT(0, d->y1[0]);
  for( ch= 1; ch< me->nch; ++ch ) {
    d->y1[ch]= (INPUTC(1, ch) - d->y1[ch])*d->cx + d->y1[ch]; 
    OUTPUT(ch, d->y1[ch]);
  }
  CALCRETURN;
}


/*=============================================================================
    highpass <cutoff> <signal>

Simple fast first-order IIR-based filter favouring high frequencies.

See also: lowpass, shelve, avg
=============================================================================*/

float calchighpass(sndobj *me, sndobj *, int);

sndobj *highpass( sndobj *cutoff, sndobj *signal )
{
  sndobj *p;
  struct fofilter *d;
  int i;
  
  p= newsndo(calchighpass, signal->name, "highpass", signal->nch, 2, cutoff, signal );
  p->private[0]= d= new(struct fofilter);
  d->cx= d->cx1= d->cy1= 0.0;
  for( i= 0; i< p->nch; ++i )
    d->y1[i]= 0.0;
  d->currparam= MAGIC;
  return p;
}

float calchighpass(sndobj *me, sndobj *caller, int innr )
{
  struct fofilter *d;
  float newcut, val;
  int ch;
  
  d= (struct fofilter *)me->private[0];
  newcut= INCALC(0);
  if( newcut!=d->currparam ) {
    if( newcut>(float)SAMPRATE/2.0 )
      newcut= (float)SAMPRATE/2.0;
    val= 2.0 - cos(newcut*2.0*M_PI/(float)SAMPRATE);
    d->cy1= val - sqrt(val*val - 1.0);
    d->currparam= newcut;
  }
// In effect, cx= cy1-1 < 0 and cx1= 0
// Not y is output, but y+x, which cancels low frequencies because cx<0
  d->y1[0]= d->cy1 * (d->y1[0] + INCALC(1));
  OUTPUT(0, d->y1[0]);
  d->y1[0] -= INPUT(1);
  for( ch= 1; ch< me->nch; ++ch ) {
    d->y1[ch]= d->cy1 * (d->y1[ch] + INPUTC(1, ch));
    OUTPUT(ch, d->y1[ch]);
    d->y1[ch] -= INPUTC(1, ch);
  }
  CALCRETURN;
}


/*=============================================================================
    fastrms <cutoff> <signal>

Square root of (fast) lowpass for squared signal.

See also: rms, lowpass
=============================================================================*/

float calcfastrms(sndobj *me, sndobj *, int);

sndobj *fastrms( sndobj *cutoff, sndobj *signal )
{
  sndobj *p;
  struct fofilter *d;
int i;

  p= newsndo(calcfastrms, signal->name, "fastrms", signal->nch, 2, cutoff, signal );
  p->private[0]= d= new(struct fofilter);
  d->cx= d->cx1= d->cy1= 0.0;
  for( i= 0; i< p->nch; ++i )
    d->y1[i]= 0.0;
  d->currparam= MAGIC;
  return p;
}

float calcfastrms(sndobj *me, sndobj *caller, int innr )
{
  struct fofilter *d;
  float newcut, val;
  int ch;
  
  d= (struct fofilter *)me->private[0];
  newcut= INCALC(0);
  if( newcut!=d->currparam ) {
    if( newcut>(float)SAMPRATE/2.0 )
      newcut= (float)SAMPRATE/2.0;
    val= 2.0 - cos(newcut*2.0*M_PI/(float)SAMPRATE);
    d->cx= 1.0 - (val - sqrt(val*val - 1.0));
    d->currparam= newcut;
  }
// In effect, cy1= 1-cx and cx1= 0
  val= INCALC(1);
  d->y1[0]= (val*val - d->y1[0])*d->cx + d->y1[0]; 
  OUTPUT(0, sqrt(d->y1[0]));
  for( ch= 1; ch< me->nch; ++ch ) {
    d->y1[ch]= (INPUTC(1, ch)*INPUTC(1, ch) - d->y1[ch])*d->cx + d->y1[ch]; 
    OUTPUT(ch, sqrt(d->y1[ch]));
  }
  CALCRETURN;
}


/*=============================================================================
    median (notmedianmode) <freq> <signal>

Outputs the median of an interval of samples.  The width of the interval is
given by the wave length corresponding to the frequency freq.  If
(notmedianmode) is !=0, not the median is computed but the average of the
highest and lowest values.

See also: avg, gaussavg, lowpass
=============================================================================*/

struct median {
  float     *sorted[PLENTY];
  long      size, halfsize, maxsize;
  float     currfreq;
  int       notmedianmode;
};

int cmpfloat(const void *a, const void *b);
float calcmedian(sndobj *me, sndobj *, int);

sndobj *median( int notmedianmode, sndobj *freq, sndobj *signal )
{
  sndobj *p;
  struct median *d;
  int ch;
  
  p= newsndo( calcmedian, signal->name, "median", signal->nch, 2, freq, signal );
  p->private[0]= d= new(struct median);
  d->size= 0;
  d->maxsize= (long)(HORIZON*SAMPRATE);
  if( (d->maxsize&1)==0 )
    --d->maxsize;
  d->currfreq= MAGIC;
  d->notmedianmode= notmedianmode;
  p->private[1]= d->sorted[0]= (float*)calloc(p->nch*d->maxsize, sizeof(float));
  for( ch= 1; ch < p->nch; ++ch )
    d->sorted[ch]= d->sorted[0] + d->maxsize;
  buf( HORIZON/2.0, HORIZON/2.0 + 1.0/SAMPRATE, 0.0, p->nch, p, 1 );
  return p;
}

int cmpfloat(const void *a, const void *b)
{
  if( *(float*)a < *(float*)b )
    return -1;
  else
    return *(float*)a > *(float*)b;
}

float calcmedian(sndobj *me, sndobj *caller, int innr )
{
  struct median *d;
  float newfreq;
  long newsize, newpos, oldpos;
  int i, ch;
  
  d= (struct median *)me->private[0];
  newfreq= INCALC(0);
  if( d->currfreq==MAGIC ) {      // first call
    if( !finite(newfreq) )
      CALCRETURN;
    d->currfreq= newfreq;
    d->size= (long)round((float)SAMPRATE/newfreq);
    if( d->size > d->maxsize )
      d->size= d->maxsize;
    if( (d->size&1)==0 )
      ++d->size;
    d->halfsize= d->size/2;
    for( i= d->halfsize; i>=0; --i ) {
      BUFINCALC(1, i, 0.0);
      FORCH
	d->sorted[ch][d->halfsize+i]= INPUTC(1, ch);
    }
    FORCH
      qsort( d->sorted[ch]+d->halfsize, d->halfsize+1, sizeof(float), cmpfloat );
  }
  else {
    if( finite(newfreq) && newfreq != d->currfreq ) {
      d->currfreq= newfreq;
      newsize= (long)round((float)SAMPRATE/newfreq);
      if( (newsize&1)==0 )
	++newsize;
      if( newsize > d->maxsize )
	newsize= d->maxsize;
      else if( !newsize )
	newsize= 1;
    }
    else
      newsize= d->size;
    if( newsize > d->size ) {
      BUFINCALC(1, d->halfsize, 0.0);
      FORCH {
	newpos= binfindf( d->sorted[ch], d->size, INPUTC(1, ch) );
	memmove( d->sorted[ch]+newpos+1, d->sorted[ch]+newpos,
		    (d->size-newpos)*sizeof(float) );
	d->sorted[ch][newpos]= INPUTC(1, ch);
      }
      ++d->size;
      ++d->halfsize;
      BUFINCALC(1, d->halfsize, 0.0);
      FORCH {
	newpos= binfindf( d->sorted[ch], d->size, INPUTC(1, ch) );
	memmove( d->sorted[ch]+newpos+1, d->sorted[ch]+newpos,
		    (d->size-newpos)*sizeof(float) );
	d->sorted[ch][newpos]= INPUTC(1, ch);
      }
      ++d->size;
    }
    else if( newsize < d->size ) {
      BUFINCALC(1, -d->halfsize-1, 0.0 );
      FORCH {
	oldpos= binfindf( d->sorted[ch], d->size, INPUTC(1, ch) );
	memmove( d->sorted[ch]+oldpos, d->sorted[ch]+oldpos+1,
		    (d->size-oldpos-1)*sizeof(float) );
      }
      --d->size;
      BUFINCALC(1, -d->halfsize, 0.0);
      FORCH {
	oldpos= binfindf( d->sorted[ch], d->size, INPUTC(1, ch) );
	memmove( d->sorted[ch]+oldpos, d->sorted[ch]+oldpos+1,
		    (d->size-oldpos-1)*sizeof(float) );
      }
      --d->halfsize;
      --d->size;
    }
    else {
      BUFINCALC(1, -d->halfsize-1, 0.0 );
      FORCH {
	oldpos= binfindf( d->sorted[ch], d->size, INPUTC(1, ch) );
//      if( oldpos==d->size )   --oldpos;// unnecessary unless there's a bug somewhere else
	memmove( d->sorted[ch]+oldpos, d->sorted[ch]+oldpos+1,
		    (d->size-oldpos-1)*sizeof(float) );
      }
      BUFINCALC(1, d->halfsize, 0.0);
      FORCH {
	newpos= binfindf( d->sorted[ch], d->size-1, INPUTC(1, ch) );
	memmove( d->sorted[ch]+newpos+1, d->sorted[ch]+newpos,
		    (d->size-1-newpos)*sizeof(float) );
	d->sorted[ch][newpos]= INPUTC(1, ch);
      }
    }
  }

  BUFINCR(1, 1);
  if( !d->notmedianmode )
    FORCH
      OUTPUT(ch, d->sorted[ch][d->halfsize]);
  else
    FORCH
      OUTPUT(ch, (d->sorted[ch][0]+d->sorted[ch][d->size-1])/2.0);
  CALCRETURN;
}



/*=============================================================================
    dcblock <in>

DC blocker - eliminates very low frequencies and baseline drift.

See also: highpass, shelve
=============================================================================*/

float calcdcblock(sndobj *me, sndobj *, int);

sndobj *dcblock( sndobj *in )
{
  sndobj *p;
  float *decay;
  
  p= newsndo( calcdcblock, in->name, "dcblock", in->nch, 1, in);
  p->private[0]= calloc(p->nch, sizeof(float));
  p->private[1]= decay= new(float);
  // determine decay constant from its value for 44.1k sampling rate:
  *decay= exp(log(0.995)*44100.0/SAMPRATE);
  return p;
}

float calcdcblock(sndobj *me, sndobj *caller, int innr )
{
  float *prevx, *decay;
  float newx;
  int ch;
  
  prevx= (float*)me->private[0];
  decay= (float*)me->private[1];
  INCALC(0);
  FORCH {
    newx= INPUTC(0, ch);
    OUTPUT(ch, *decay * me->ch[ch] + newx - prevx[ch]);
    prevx[ch]= newx;
  }
  CALCRETURN;
}


/*=============================================================================
    peak <sharpness> <freq> <signal>

Peak resonance filter.  sharpness determines the narrowness of the peak and has
to be in [0, 1] (useful range: above 0.9).  Good for wahwah.

See also: formant, peaknotch
=============================================================================*/

struct peak {
  float     coeffs[10];
  float     currR, currf;
};

float calcpeak(sndobj *me, sndobj *caller, int innr);

sndobj *peak( sndobj *sharpness, sndobj *freq, sndobj *signal )
{
  sndobj *p;
  struct peak *d;
  float *forw, *back;
  
  d= new(struct peak);
  d->currR= -1.0;
  d->currf= -1.0;
  forw= d->coeffs;
  forw[0]= 0.0;
  forw[1]= 1.0;
  forw[2]= 2.0;
  forw[3]= -1.0;
  forw[4]= END;
  back= d->coeffs+5;
  back[0]= 1.0;
  back[1]= 0.0;
  back[2]= 2.0;
  back[3]= 0.0;
  back[4]= END;
  p= filter( forw, back, signal );
  p->type= "peak";
  p->calc= calcpeak;
  p->private[PLENTY-1]= d;
  addinput(p, sharpness);
  addinput(p, freq);
  return p;
}


float calcpeak(sndobj *me, sndobj *caller, int innr)
{
  struct peak *d;
  float *forw, *back;
  float newR, newfreq, old1pR2;
  
  d= (struct peak *)me->private[PLENTY-1];
  newR= INCALC(2);  // first two inputs are forwardfeed and backfeed
  newfreq= INCALC(3);
  if( newR != d->currR ) {
    if( newR > 1.0 )
      newR= 1.0;
    else if( newR < 0.0 )
      newR= 0.0;
  }
  if( newR != d->currR ) {
    forw= d->coeffs;
    back= d->coeffs+5;
    if( newfreq != d->currf ) {
      back[3]= -newR*newR;
      back[1]= cos(2*M_PI*newfreq/SAMPRATE)*(1-back[3]);
      d->currf= newfreq;
    }
    else {
      old1pR2= 1.0-back[3];
      back[3]= -newR*newR;
      back[1] *= (1.0-back[3])/old1pR2;
    }
    forw[1]= (1.0 + back[3])/2.0;
    forw[3]= -forw[1];
    d->currR= newR;
  }
  else if( newfreq != d->currf ) {
    back= d->coeffs+5;
    back[1]= cos(2*M_PI*newfreq/SAMPRATE)*(1-back[3]);
    d->currf= newfreq;
  }
  return calcfilter( me, caller, innr );
}


/*=============================================================================
    peaknotch <freq> <gain> <bandwidth> <signal>

Peak/notch filter.  <freq> gives the centre frequency, <gain> the gain factor,
and <bandwidth> the bandwidth of the filter.

See also: peak, formant
=============================================================================*/

struct peaknotch {
  int   baseinp;
  float currfreq, currgain, currbw;
  float back[5], forw[7];
};

float calcpeaknotch(sndobj *me, sndobj *caller, int innr);

sndobj *peaknotch(sndobj *freq, sndobj *gain, sndobj *bandwidth, sndobj *signal)
{
  sndobj *p;
  struct peaknotch *d;

  d= new(struct peaknotch);
  d->forw[0]= 0;
  d->forw[1]= 1;
  d->forw[2]= 1;
  d->forw[3]= 0;
  d->forw[4]= 2;
  d->forw[5]= 0;
  d->forw[6]= END;
  d->back[0]= 1;
  d->back[1]= 0;
  d->back[2]= 2;
  d->back[3]= 0;
  d->back[4]= END;
  d->currfreq= -1;
  d->currgain= -1;
  d->currbw= -1;
  p= filter(d->forw, d->back, signal);
  d->baseinp= p->nin;
  p->name= p->type= "peaknotch";
  p->calc= calcpeaknotch;
  p->private[PLENTY-1]= d;
  addinput(p, freq);
  addinput(p, gain);
  addinput(p, bandwidth);
  return p;
}

float calcpeaknotch(sndobj *me, sndobj *caller, int innr)
{
  struct peaknotch *d;
  float freq, gain, bw;

  d= (struct peaknotch *)me->private[PLENTY-1];
  freq= INCALC(d->baseinp);
  if( freq < 0 )
    freq= 0;
  else if( freq > SAMPRATE/2 )
    freq= SAMPRATE/2;
  gain= INCALC(d->baseinp+1);
  if( gain < 0 )
    gain= 0;
  bw= INCALC(d->baseinp+2);
  if( bw < 0 )
    bw= 0;
  else if( bw > SAMPRATE/2 )
    bw= SAMPRATE/2;
  if( freq != d->currfreq || gain != d->currgain || bw != d->currbw )
  {
    double relf, relbw, tanquot;

//    printf("freq %g %g  gain %g %g  bw %g %g\n", (double)freq, (double)d->currfreq, (double)gain, (double)d->currgain, (double)bw, (double)d->currbw);
    d->currfreq= freq;
    d->currgain= gain;
    d->currbw= bw;
    relf= 2.0 * M_PI * (double)freq / (double)SAMPRATE;
    relbw= 2.0 * M_PI * (double)bw / (double)SAMPRATE;
    tanquot= (1.0 - tan(relbw/2)) / (1.0 + tan(relbw/2));
    d->forw[1]= 0.5 * (1.0 + gain + tanquot - tanquot * gain);
    d->forw[3]= 0.5 * -2 * cos(relf) * (1.0 + tanquot);
    d->forw[5]= 0.5 * (1 - gain + tanquot + gain * tanquot);
    d->back[1]= 2 * cos(relf) / (1 + tan(relbw/2));
    d->back[3]= -tanquot;
//    printf("forw: %g %g %g  backw: %g %g \n", (double)d->forw[1], (double)d->forw[3], (double)d->forw[5], (double)d->back[1], (double)d->back[3] );
  }
  return calcfilter(me, caller, innr);
}


/*=============================================================================
    boost <freq> <gain> <bandwidth> <signal>

Another peak (or notch) filter centred around frequency <freq>, with gain <gain>
and bandwidth <bandwidth>.  Guaranteed to diverge under any circumstances.

See also: formant, peak, peaknotch
=============================================================================*/

struct boost {
  int   baseinp;
  float currfreq, currgain, currbw;
  float forw[7], back[7];
};

float calcboost(sndobj *me, sndobj *caller, int innr);

sndobj *boost(sndobj *freq, sndobj *gain, sndobj *bandwidth, sndobj *signal)
{
  sndobj *p;
  struct boost *d;

  d= new(struct boost);
  d->forw[0]= 0;
  d->forw[1]= 1;
  d->forw[2]= 1;
  d->forw[3]= 0;
  d->forw[4]= 2;
  d->forw[5]= 0;
  d->forw[6]= END;
  d->back[0]= 1;
  d->back[1]= 0;
  d->back[2]= 2;
  d->back[3]= 0;
  d->back[4]= END;
  d->currfreq= -1;
  d->currgain= -1;
  d->currbw= -1;
  p= filter( d->forw, d->back, signal );
  d->baseinp= p->nin;
  p->name= p->type= "boost";
  p->calc= calcboost;
  p->private[PLENTY-1]= d;
  addinput(p, freq);
  addinput(p, gain);
  addinput(p, bandwidth);
  return p;
}

float calcboost(sndobj *me, sndobj *caller, int innr)
{
  struct boost *d;
  float freq, gain, bw;

  d= (struct boost *)me->private[PLENTY-1];
  freq= INCALC(d->baseinp);
  if( freq < 0 )
    freq= 0;
  else if( freq > SAMPRATE/2 )
    freq= SAMPRATE/2;
  gain= INCALC(d->baseinp+1);
  if( gain < 0 )
    gain= 0;
  bw= INCALC(d->baseinp+2);
  if( bw < 0 )
    bw= 0;
  else if( bw > SAMPRATE/2 )
    bw= SAMPRATE/2;
  if( freq != d->currfreq || gain != d->currgain || bw != d->currbw )
  {
    double tanf, tanf_q, tanf2, norm;

    d->currfreq= freq;
    d->currgain= gain;
    d->currbw= bw;
    tanf= tan( M_PI*(double)freq/(double)SAMPRATE );
    tanf2= tanf * tanf;
    tanf_q= tanf * bw / (double)SAMPRATE;
    norm= 1.0 / (1.0 + tanf_q + tanf2);
    printf("%g %g %g %g \n", tanf, tanf2, tanf_q, norm );
    d->forw[1]= norm * (1.0 + gain * tanf_q + tanf2);
    d->forw[3]= norm * 2.0 * (tanf2 - 1.0);
    d->forw[5]= norm * (1.0 - gain * tanf_q + tanf2);
    d->back[1]= norm * 2.0 * (tanf2 - 1.0);
    d->back[3]= norm * (1.0 - tanf_q + tanf2);
    printf("boost: coeffs: forw %g %g %g  back %g %g \n", (double)d->forw[1],
(double)d->forw[3], (double)d->forw[5], (double)d->back[1], (double)d->back[3]);
  }
  return calcfilter(me, caller, innr);
}


/*=============================================================================
    rmpeak (peakminheight) <signal>

This object was created to remove artificial 1-2 samples wide peaks from a
badly copied music CD track.  It has nothing to do with the peak object.
=============================================================================*/

struct rmpeak {
  float     minheight;
  float     lastout[PLENTY];
};

float calcrmpeak(sndobj *me, sndobj *, int);

sndobj *rmpeak( float peakminheight, sndobj *signal )
{
  sndobj *p;
  struct rmpeak *d;
  int count;
  
  p= newsndo( calcrmpeak, "rmpeak", "rmpeak", signal->nch, 1, signal );
  buf( 0, 3.5/SAMPRATE, 0.0, signal->nch, p, 0 );
  p->private[0]= d= new(struct rmpeak);
  d->minheight= fabsf(peakminheight);
  for( count= 0; count< PLENTY; ++count )
    d->lastout[count]= 0.0;
  return p;
}


float calcrmpeak(sndobj *me, sndobj *caller, int innr )
{
  static float diffs_[PLENTY], diffs[PLENTY], diffs1[PLENTY], diffs2[PLENTY];
  struct rmpeak *d;
  int ch;
  
  d= (struct rmpeak *)me->private[0];
  FORCH
    diffs_[ch]= -d->lastout[ch];
  BUFINCALC(0, 1, 0.0);
  FORCH {
    diffs[ch]= INPUTC(0, ch);
    diffs1[ch]= -INPUTC(0, ch);
  }
  BUFINCALC(0, 2, 0.0);
  FORCH {
    diffs1[ch] += INPUTC(0, ch);
    diffs2[ch]= -INPUTC(0, ch);
  }
  BUFINCALC(0, 3, 0.0);
  FORCH
    diffs2[ch] += INPUTC(0, ch);
  BUFINCALC(0, 0, 0.0);
  FORCH {
    diffs_[ch] += INPUTC(0, ch);
    diffs[ch] -= INPUTC(0, ch);
    if( fabsf(diffs_[ch]) < d->minheight )
      OUTPUT(ch, INPUTC(0, ch));
    else if( (diffs_[ch]>0 && diffs[ch]>0 && diffs1[ch]>0) || 
		(diffs_[ch]<0 && diffs[ch]<0 && diffs1[ch]<0) )
      OUTPUT(ch, INPUTC(0, ch));
    else if( fabsf(diffs_[ch]+diffs[ch]) < 2*d->minheight )
      OUTPUT(ch, INPUTC(0, ch)-(diffs_[ch]-diffs[ch])/2.0);
    else if( fabsf(diffs_[ch]+diffs[ch]+diffs1[ch]) < 3*d->minheight )
      OUTPUT(ch, INPUTC(0, ch)-diffs_[ch]*2.0/3.0+(diffs[ch]+diffs1[ch])/3.0);
    else if( fabsf(diffs_[ch]+diffs[ch]+diffs1[ch]+diffs2[ch]) < 4*d->minheight )
      OUTPUT(ch, INPUTC(0, ch)-diffs_[ch]*3.0/4.0+(diffs[ch]+diffs1[ch]+diffs2[ch])/4.0);
    else  {
//      printf("calcrmpeak: giving up, lastval: %g, current val: %g, first "
//                "diff: %g, 2: %g, 3: %g, 4: %g\n", 
//        d->lastout[ch], INPUTC(0, ch), diffs_[ch], diffs_[ch]+diffs[ch], 
//        diffs_[ch]+diffs[ch]+diffs1[ch], 
//        (diffs_[ch]+diffs[ch]+diffs1[ch]+diffs2[ch]));
//      if( fabsf(d->lastout[ch]) < fabsf(INPUTC(0, ch)) )
//        OUTPUT(ch, d->lastout[ch]);
//      else
//        OUTPUT(ch, INPUTC(0, ch));
      OUTPUT(ch, d->lastout[ch]);
    }
    d->lastout[ch]= me->ch[ch];
  }
  BUFINCR(0, 1);
  CALCRETURN;
}



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
#include <math.h>
#include "sndsys.h"


/*=============================================================================
    karplusfilt (loss) (blend) <signal>

Filter used in loop of Karplus-Strong string model.  It outputs the mean of two
successive sample values, possibly with the following modifications: (loss) is
a decay factor by which the result is multiplied.  If (blend) is < 1, the
result is negated with the probability 1 - (blend).  The latter is useful for
producing drum sounds.

See also: karplusbypass, stringfeedback
=============================================================================*/

struct karplusfilt {
  float loss, blend, prevval;
  rngdesc *rng;
};

float calckarplusfilt(sndobj *me, sndobj *caller, int innr);

sndobj *karplusfilt(float loss, float blend, sndobj *signal)
{
  sndobj *p;
  struct karplusfilt *d;
  
  p= newsndo(calckarplusfilt, "karplusfilt", "karplusfilt", 1, 1, signal);
  p->private[0]= d= new(struct karplusfilt);
  d->loss= 0.5f * loss;
  d->blend= blend * 2 - 1.0f;
  d->prevval= 0.0f;
  p->private[1]= d->rng= makerng();
  return p;
}

float calckarplusfilt(sndobj *me, sndobj *caller, int innr)
{
  struct karplusfilt *d;
  
  d= (struct karplusfilt *)me->private[0];
  INCALC(0);
  if( flatrng(d->rng) <= d->blend )
    OUTPUT(0, d->loss * (INPUT(0) + d->prevval));
  else
    OUTPUT(0, - d->loss * (INPUT(0) + d->prevval));
  d->prevval= INPUT(0);
  CALCRETURN;
}


/*=============================================================================
    karplusbypass (len) (ratio) (exponent) <signal>

Nonlinear filter for modification of the Karplus-Strong string model.  The
signal value is raised to the power (exponent) (keeping the sign); the part of
it given by (ratio) is used right away while the rest is passed through a delay
line of length (len).  This is a nonlinear operation for (exponent)!=1 because
the (ratio) is applied to the power of the sample value.

See also: karplusfilt, stringfeedback
=============================================================================*/

struct karplusbypass {
  float *pipe, *rw, *top;
  rngdesc *rng;
  int	pipelen;
  float	ratio;
  double exponent, exponentm1, inverseexp;
};

float calckarplusbypass(sndobj *me, sndobj *caller, int innr);

sndobj *karplusbypass( double len, float ratio, double exponent, sndobj *signal)
{
  sndobj *p;
  struct karplusbypass *d;
  
  p= newsndo(calckarplusbypass, "karplusbypass", "karplusbypass", 1, 1,
				    signal);
  p->private[0]= d= new(struct karplusbypass);
  d->pipelen= (int)(len * SAMPRATE);
  p->private[1]= d->pipe= news(float, d->pipelen);
  d->top= d->pipe + d->pipelen;
  for( d->rw= d->pipe; d->rw< d->top; ++d->rw )
    *d->rw= 0.0;
  d->rw= d->pipe;
  d->ratio= ratio;
  d->exponent= exponent;
  d->exponentm1= exponent - 1.0;
  d->inverseexp= 1.0/exponent;
  d->rng= makerng();
  return p;
}

float calckarplusbypass(sndobj *me, sndobj *caller, int innr)
{
  struct karplusbypass *d;
  float *read;
  double sqrout;
  float pownew;
  
  d= (struct karplusbypass *)me->private[0];
  pownew= INCALC(0);
  pownew *= (float)pow(fabsf(pownew), d->exponentm1);
  read= d->rw + (int)floor( (0.5+0.5*flatrng(d->rng))*d->pipelen );
  if( read >= d->top )
    read -= d->pipelen;
  sqrout= *d->rw + d->ratio * pownew;
  if( sqrout >= 0.0 )
    OUTPUT(0, pow(sqrout, d->inverseexp));
  else
    OUTPUT(0, -pow(-sqrout, d->inverseexp));
  *d->rw++ = (1.0f - d->ratio) * pownew;
  if( d->rw >= d->top )
    d->rw= d->pipe;
  CALCRETURN;
}


/*=============================================================================
    stringfeedback (superinit) (si_loopscale) (si_initscale) <cue> <delay> <feed> <init>

This object combines feedback (for instance of a Karplus-Strong string loop)
with an initialisation signal.  A nonzero value from the <cue> input indicates
that the string was plucked and the delay loop should be initialised.  Values
from <init> are output exclusively for one period of the loop, given by the
value of <delay> at the time <cue> is non-zero.  After that, <init> values are
mixed with the feedback <feed>.  The initialisation data is scaled by
(si_initscale) and the feedback data by (si_loopscale) before they are added
up.  If (superinit) is negative, this continues forever.  Otherwise, after
(superinit) times the delay length, only the feedback is passed through.  A
negative value from <cue> mutes the string - 0 is output.

See also: karplusfilt, karplusbypass
=============================================================================*/

struct stringfeedback {
  int initcount, addinitcount;
  float cueval, addinitfact, loopscale, initscale;
};

float calcstringfeedback(sndobj *me, sndobj *caller, int innr);

sndobj *stringfeedback(float superinit, float si_loopscale, float si_initscale,
    			sndobj *cue, sndobj *delay, sndobj *feed, sndobj *init)
{
  sndobj *p;
  struct stringfeedback *d;
  
  p= newsndo(calcstringfeedback, "stringfeedback", "stringfeedback", 1, 4,
	    cue, delay, feed, init );
  p->private[0]= d= new(struct stringfeedback);
  d->initcount= -1;
  d->addinitcount= 0;
  d->cueval= 0.0;
  d->addinitfact= superinit;
  d->loopscale= si_loopscale;
  d->initscale= si_initscale;
  if( superinit < 0 )
    d->addinitcount= 1.0;
  return p;
}

float calcstringfeedback(sndobj *me, sndobj *caller, int innr)
{
  struct stringfeedback *d= (struct stringfeedback *)me->private[0];

  if( INCALC(0)!=0.0 ) {
    d->cueval= INPUT(0);
    INCALC(1);
    d->initcount= (int)((float)SAMPRATE*INPUT(1));
    if( d->addinitfact >= 0 )
      d->addinitcount= (int)((float)SAMPRATE*d->addinitfact*INPUT(1));
  }
  else
    INSKIP(1);
  if( d->initcount < 0 ) {
    INSKIP(2);
    INSKIP(3);
    OUTPUT(0, 0.0);
  }
  else if( d->initcount > 0 ) {
    INSKIP(2);
    OUTPUT(0, d->cueval*INCALC(3));
    --d->initcount;
  }
  else if( d->addinitcount > 0 ) {
    OUTPUT(0, d->loopscale*INCALC(2)+d->initscale*d->cueval*INCALC(3));
    if( d->addinitfact >= 0 )
      --d->addinitcount;
  }
  else {
    OUTPUT(0, INCALC(2));
    INSKIP(3);
  }
  CALCRETURN;
}




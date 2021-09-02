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
#include <math.h>

#include "sndsys.h"


/*=============================================================================
    ramp (offset) [list]

Outputs a function of time that is piecewise either linear or exponential.
(offset) is a time interval that will be subtracted from the internal time
variable, so in effect added to the times in [list].  [list] contains 3*n+2
double values plus a concluding END.  The values are time-value pairs with mode
values in between.  mode can be RAMP_L (linear ramp), RAMP_E (exponential
decay/rise) or RAMP_C (piecewise constant, abrupt step when next time is
reached).

See also: saw, smoothstep, avg
=============================================================================*/

struct ramp {
  double    time;
  double    *list;
  double    outval, inc;
  int       mode, countdown, firstcall, skipped, lastramp;
};

float calcramp(sndobj *me, sndobj *, int);
void skipramp(sndobj *me, sndobj *, int);

#define RAMP_L_RECALC       441000
#define RAMP_E_RECALC       44100
#define RAMP_C_RECALC       INT_MAX

sndobj *ramp( double offset, double *list )
{
  sndobj *p;
  struct ramp *d;
  
  if( *list==END ) {
    fprintf(stderr, "Warning: empty list for ramp. Returning c(0).\n");
    return c(0);
  }
  p= newsndo( calcramp, "ramp", "ramp", 1, 0);
  p->skip= skipramp;
  p->private[0]= d= new(struct ramp);
  d->time= -offset;
  d->list= list;
  d->outval= 0.0;
  d->inc= 0.0;
  d->firstcall= 1;
  d->skipped= 0;
  d->lastramp= 0;   // important for firstcall!
  return p;
}

float calcramp(sndobj *me, sndobj *caller, int innr)
{
  struct ramp *d;
  double *list;
  
  d= (struct ramp *)me->private[0];
  list= d->list;
  if( !d->lastramp && (d->time >= list[3] || d->firstcall || d->skipped) ) {
    d->skipped= 0;
    if( d->firstcall ) {
      d->time += starttime;
      d->firstcall= 0;
    }
    while( list[2]!=END && d->time >= list[3] )
      list += 3;
    d->list= list;
    if( list[2]!=END ) {
      d->mode= (int)list[2];
      if( d->mode==RAMP_E && (!list[1] || !list[4] || 
	  (list[1]>0 && list[4]<0) || (list[1]<0 && list[4]>0)) )
	d->mode= RAMP_L;        // no sign change possible for exponential ramp
    }
    else {                      // end of list reached
      d->mode= RAMP_C;
      d->lastramp= 1;
    }
    switch( d->mode ) {
      case RAMP_E:d->inc= exp(log(list[4]/list[1])/(list[3]-list[0])/SAMPRATE);
		  d->outval= list[1]*exp(log(list[4]/list[1])*(d->time-list[0])/(list[3]-list[0]));
		  d->countdown= RAMP_E_RECALC;
		  break;
      case RAMP_C:d->outval= list[1];
		  d->countdown= RAMP_C_RECALC;
		  break;
      case RAMP_L:
	  default:d->inc= (list[4]-list[1])/(list[3]-list[0])/SAMPRATE;
		  d->outval= list[1]+(list[4]-list[1])*(d->time-list[0])/(list[3]-list[0]);
		  d->countdown= RAMP_L_RECALC;
		  break;
    }
  }
  if( !d->countdown )
    switch( d->mode ) {
      case RAMP_E:d->inc= exp(log(list[4]/d->outval)/(list[3]-d->time)/SAMPRATE);
		  d->countdown= RAMP_E_RECALC;
		  break;
      case RAMP_C:d->countdown= RAMP_C_RECALC;
		  break;
      case RAMP_L:
	  default:d->inc= (list[4]-d->outval)/(list[3]-d->time)/SAMPRATE;
		  d->countdown= RAMP_L_RECALC;
		  break;
    }
  OUTPUT(0, d->outval);
  switch( d->mode ) {
    case RAMP_E:d->outval *= d->inc;
		break;
    case RAMP_C:break;
    case RAMP_L:
	default:d->outval += d->inc;
		break;
  }
  d->time += 1.0/SAMPRATE;
  --d->countdown;
  CALCRETURN;
}

void skipramp(sndobj *me, sndobj *caller, int innr)
{
  struct ramp *d;
  
  d= (struct ramp *)me->private[0];
  d->skipped= 1;
  d->time += 1.0/SAMPRATE;
  if( d->countdown > 0 )
    --d->countdown;
}


/*=============================================================================
    at (offset) [list]

Outputs one sample wide pulses at the times given in [list], otherwise zero.
[list] contains times and corresponding pulse heights in ascending temporal
order and is concluded by the value END.  (offset) allows to shift the pulses
in time without replacing all times in [list].  (offset) is subtracted from the
time variable internal to this sndobj, so in effect added to the times in
[list].  A positive (offset) delays the events in [list].

See also: cron, interleave, slide
=============================================================================*/

struct at {
  double    time;
  double    *list;
  int       firstcall;
};

float calcat(sndobj *me, sndobj *, int);
void skipat(sndobj *me, sndobj *, int);

sndobj *at( double offset, double *list )
{
  sndobj *p;
  struct at *d;
  
  if( *list==END || list[1]==END ) {
    printf("Warning: at: Empty list. Returning c(0).\n");
    return c(0);
  }
  p= newsndo(calcat, "at", "at", 1, 0);
  p->skip= skipat;
  p->private[0]= d= new(struct at);
  d->time= -offset;
  d->list= list;
  d->firstcall= 1;
  return p;
}

float calcat(sndobj *me, sndobj *caller, int innr)
{
  struct at *d;
  
  d= (struct at *)me->private[0];
  if( d->firstcall ) {
    d->firstcall= 0;
    d->time += starttime;
    if( *d->list!=END && d->list[1]!=END && d->time >= *d->list ) {
      while( *d->list!=END && d->list[1]!=END && d->time >= *d->list )
	d->list += 2;
      if( d->list[-2] > d->time - 1.0/SAMPRATE )
	OUTPUT(0, d->list[-1]);
    }
    else
      OUTPUT(0, 0.0);
  }
  else if( *d->list!=END && d->list[1]!=END && d->time >= *d->list ) {
    while( *d->list!=END && d->list[1]!=END && d->time >= *d->list )
      d->list += 2;
    OUTPUT(0, d->list[-1]);
  }
  else
    OUTPUT(0, 0.0);
  d->time += 1.0/SAMPRATE;
  CALCRETURN;
}

void skipat(sndobj *me, sndobj *caller, int innr)
{
  struct at *d;
  
  d= (struct at *)me->private[0];
  if( d->firstcall ) {
    d->firstcall= 0;
    d->time += starttime;
  }
  d->time += 1.0/SAMPRATE;
}


/*=============================================================================
    cron (offset) (interval) (value) [exceptions]

Outputs a (value) at regular intervals, starting at the earliest time which
differs from (offset) by a multiple of (interval).  If [exceptions] is
non-NULL, it has to contain pairs of times and values followed by the value
END.  The values will then be substituted for the regular (value) at the time
closest to the one given.  The times in [exceptions] are understood to be
relative to (offset) (ie defined for (offset)==0) so that all events can be
shifted in time by changing (offset).

See also: at, interleave
=============================================================================*/

struct cron {
  double    inttime, halfpoint, interval;
  double    *exc;
  float     value;
  int       firstcall;
};

float calccron(sndobj *me, sndobj *, int);

sndobj *cron( double offset, double interval, float value, double *exceptions )
{
  sndobj *p;
  struct cron *d;
  
  p= newsndo(calccron, "cron", "cron", 1, 0);
  p->skip= dontskip;
  p->private[0]= d= new(struct cron);
  d->inttime= -offset;
  d->halfpoint= -offset;
  d->interval= interval;
  if( exceptions && *exceptions!=END && exceptions[1]!=END )
    d->exc= exceptions;
  else
    d->exc= NULL;
  d->value= value;
  d->firstcall= 1;
  return p;
}

float calccron(sndobj *me, sndobj *caller, int innr)
{
  struct cron *d;
  
  d= (struct cron *)me->private[0];
  if( d->firstcall ) {
    d->firstcall= 0;
    d->inttime += starttime;
    while( d->inttime > d->interval )
      d->inttime -= d->interval;
    while( d->inttime < 1.0/SAMPRATE )
      d->inttime += d->interval;
    if( d->exc ) {
      d->halfpoint += starttime - d->inttime + 0.5*d->interval;
      while( *d->exc!=END && d->exc[1]!=END && d->halfpoint >= *d->exc )
	d->exc += 2;
    }
  }
  if( d->inttime >= d->interval ) {
    d->inttime -= d->interval;
    if( d->exc ) {
      d->halfpoint += d->interval;
      if( *d->exc!=END && d->exc[1]!=END && d->halfpoint >= *d->exc ) {
	OUTPUT(0, d->exc[1]);
	do {  d->exc += 2; }
	while( *d->exc!=END && d->exc[1]!=END && d->halfpoint >= *d->exc );
      }
      else OUTPUT(0, d->value);
    }
    else OUTPUT(0, d->value);
  }
  else OUTPUT(0, 0.0);
  d->inttime += 1.0/SAMPRATE;
  CALCRETURN;
}


/*=============================================================================
    exttrigger [data]

External trigger (eg for user input).  Every time *data is !=0, the value is
output and it is reset to zero.  Zero is output otherwise.

See also: extc
=============================================================================*/

float calcexttrigger(sndobj *me, sndobj *, int);

sndobj *exttrigger( float *data )
{
  sndobj *p;
  
  p= newsndo( calcexttrigger, "exttrigger", "exttrigger", 1, 0 );
  p->private[0]= new(float*);
  *(float**)p->private[0]= data;
  return p;
}

float calcexttrigger(sndobj *me, sndobj *caller, int innr)
{
  float *d;
  
  d= *(float **)me->private[0];
  if( *d )
  {
    OUTPUT(0, *d);
    *d= 0.0;
  }
  else
    OUTPUT(0, 0.0);
  CALCRETURN;
}


/*=============================================================================
    expodecay <decaylength> <trigger>

Converts a trigger input to exponential decay.  Each time <trigger> is !=0, the
output is set to its value, which subsequently decays exponentially.
<decaylength> gives the time after which the output signal should have decayed
to 1/e of its value.  It may change continuously.

See also: adsr, idexp, halfgauss
=============================================================================*/

struct expodecay {
  double    currval, factor;
  float     currentlength;
};

float calcexpodecay(sndobj *me, sndobj *, int);

sndobj *expodecay( sndobj *decaylength, sndobj *trigger )
{
  sndobj *p;
  struct expodecay *d;
  
  p= newsndo( calcexpodecay, "expodecay", "expodecay", 1, 2, decaylength, trigger );
  p->skip= dontskip;
  p->private[0]= d= new( struct expodecay );
  d->currentlength= -1.0;
  d->currval= 0.0;
  return p;
}

float calcexpodecay(sndobj *me, sndobj *caller, int innr)
{
  struct expodecay *d;
  float decaylength, newval;
  
  d= (struct expodecay *)me->private[0];
  decaylength= fabsf(INCALC(0));
  if( decaylength != d->currentlength ) {
    if( decaylength > 0 )
      d->factor= exp( -1.0/((double)decaylength*SAMPRATE) );
    else
      d->factor= 0.0;
    d->currentlength= decaylength;
  }
  if( 0!=(newval= INCALC(1)) ) {
    d->currval= newval;
    OUTPUT(0, newval);
  }
  else {
    d->currval *= d->factor;
    OUTPUT(0, d->currval);
  }
  CALCRETURN;
}


/*=============================================================================
    idexp <decaylength> <trigger>

Outputs a curve proportional to x*exp(-x/decaylength) at each trigger.  The
curve's maximum is normalised to the value of the trigger input.  To avoid
sudden changes in the output because of changing normalisation, decaylength is
evaluated only when a trigger occurs.

This envelope has it maximum <decaylength> after the trigger; its rising slope
reaches half the maximal value at approximately 0.234 times this value.

See also: adsr, expodecay, halfgauss
=============================================================================*/

struct idexp {
  double    currexpval, factor, time, norm;
};

float calcidexp(sndobj *me, sndobj *, int);

sndobj *idexp( sndobj *decaylength, sndobj *trigger )
{
  sndobj *p;
  struct idexp *d;
  
  p= newsndo( calcidexp, "idexp", "idexp", 1, 2, decaylength, trigger );
  p->skip= dontskip;
  p->private[0]= d= new( struct idexp );
  d->currexpval= 0.0;
  d->factor= 0.0;
  d->norm= 1.0;
  return p;
}

float calcidexp(sndobj *me, sndobj *caller, int innr)
{
  struct idexp *d;
  float trigger, decaylength;

  d= (struct idexp *)me->private[0];
  if( 0!=(trigger= INCALC(1)) ) {
    d->currexpval= trigger;
    d->time= 0.0;
    decaylength= fabsf(INCALC(0));
    if( !decaylength ) {
      d->factor= 0.0;
      d->norm= 1.0;
    }
    else {
      d->factor= exp( -1.0/(decaylength*SAMPRATE) );
      d->norm= M_E/decaylength;
    }
    OUTPUT(0, 0.0);
  }
  else {
    INSKIP(0);
    d->currexpval *= d->factor;
    d->time += 1.0/(double)SAMPRATE;
    OUTPUT(0, d->currexpval*d->time*d->norm);
  }
  CALCRETURN;
}


/*=============================================================================
    halfgauss <halflength> <trigger>

Envelope in the shape of half a Gaussian.  Whenever <trigger> is non-zero,
a new curve normalised to the trigger value is started.  It reaches half its
amplitude after <halflength> seconds.  <halflength> is evaluated at the time
of the trigger.

Compared to idexp with the same length parameter, this envelope shaper stays
loud for longer, but decays much more quickly after ~1.5*<halflength>.

See also: adsr, expodecay, idexp
=============================================================================*/

struct halfgauss {
  double time, trigger, expcoeff;
};

float calchalfgauss(sndobj *me, sndobj *, int);

sndobj *halfgauss( sndobj *halflength, sndobj *trigger )
{
  sndobj *p;
  struct halfgauss *d;

  p= newsndo( calchalfgauss, "halfgauss", "halfgauss", 1, 2, halflength, trigger );
  p->skip= dontskip;
  p->private[0]= d= new(struct halfgauss);
  d->time= d->trigger= d->expcoeff= 0.0;
  return p;
}

float calchalfgauss(sndobj *me, sndobj *caller, int innr)
{
  struct halfgauss *d;
  double halflength;
  float trigger;

  d= (struct halfgauss *)me->private[0];
  if( 0!= (trigger= INCALC(1)) ) {
    d->trigger= trigger;
    d->time= 0;
    halflength= INCALC(0);
    d->expcoeff= -M_LN2 / (halflength*halflength);
    OUTPUT(0, trigger);
  }
  else {
    INSKIP(0);
    d->time += 1.0/SAMPRATE;
    OUTPUT(0, d->trigger * exp( d->expcoeff*d->time*d->time ) );
  }
  CALCRETURN;
}


/*=============================================================================
    adsr (attack) (decay) (sustain) (release) (sustainfrac) <trigger>

Attack, Decay, Sustain, Release envelope shaper.  Whenever the trigger input is
non-zero, a new envelope is started.  Its value gives the volume.  The double
parameters (attack), (decay), (sustain) and (release) give the length of the
respective phases.  (sustainfrac) is the fraction of the total volume which is
retained during the sustain phase.  If the argument (sustain) is negative, the
end of the sustain phase is signaled by a negative value on the trigger input.

See also: expodecay, idexp, halfgauss
=============================================================================*/

struct adsr {
  double    attack, adecay, adsustain, adsrelease;
  double    time, volume;
  double    sustainfrac;
  int	    variablesustain;
};

float calcadsr(sndobj *me, sndobj *, int);

sndobj *adsr( double attack, double decay, double sustain, double release, 
		float sustainfrac, sndobj *trigger )
{
  sndobj *p;
  struct adsr *d;
  
  p= newsndo( calcadsr, "adsr", "adsr", 1, 1, trigger );
  p->skip= dontskip;
  p->private[0]= d= new(struct adsr);
  d->attack= attack;
  d->adecay= attack + decay;
  if( sustain >= 0 ) {
    d->adsustain= d->adecay + sustain;
    d->adsrelease= d->adsustain + release;
    d->variablesustain= 0;
  }
  else {
    d->adsustain= 0.0;
    d->adsrelease= release;
    d->variablesustain= 1;
  }
  d->sustainfrac= sustainfrac;
  d->time= 0.0;
  return p;
}

float calcadsr(sndobj *me, sndobj *caller, int innr)
{
  struct adsr *d;
  float trig;
  
  d= (struct adsr *)me->private[0];
  trig= INCALC(0);
  if( trig != 0.0 )
    if( d->variablesustain )
      if( trig < 0.0 )
      	d->adsustain= d->time;
      else {
      	d->adsustain= 0.0;
        d->time= 0.0;
        d->volume= trig;
      }
    else {
      d->time= 0.0;
      d->volume= trig;
    }
  if( d->time < d->attack ) {
    OUTPUT(0, d->volume * d->time / d->attack);
    if( !finite(me->ch[0]) )
      OUTPUT(0, d->volume);
  }
  else if( d->time < d->adecay ) {
    OUTPUT(0, d->volume * (1.0 - (d->time - d->attack) *
			    (1.0 - d->sustainfrac)/(d->adecay - d->attack)) );
    if( !finite(me->ch[0]) )
      OUTPUT(0, (1.0 - d->sustainfrac)*d->volume);
  }
  else
    if( d->variablesustain )
    {
      if( d->adsustain > 0.0 ) {
        if( d->time < d->adsustain + d->adsrelease ) {
	  OUTPUT(0, d->sustainfrac*d->volume * (1.0 - 
		    (d->time - d->adsustain) / d->adsrelease) );
	  if( !finite(me->ch[0]) )
	    OUTPUT(0, 0.0);
	}
	else  OUTPUT(0, 0.0);
      }
      else
	OUTPUT(0, d->sustainfrac*d->volume);
    }
    else {
      if( d->time < d->adsustain )
        OUTPUT(0, d->sustainfrac*d->volume);
      else if( d->time < d->adsrelease ) {
        OUTPUT(0, d->sustainfrac*d->volume * (1.0 - (d->time - d->adsustain) /
			 (d->adsrelease - d->adsustain)));
	if( !finite(me->ch[0]) )
	  OUTPUT(0, 0.0);
      }
      else
        OUTPUT(0, 0.0);
    }

  d->time += 1.0/(double)SAMPRATE;
  CALCRETURN;
}



/*=============================================================================
    flipflop <trigger> <signal>

Does what a flip-flop does: Whenever the (first channel of the) trigger is !=0,
the input is stored and output until the next trigger.  The number of output
channels is the same as for the input.  Unlike in real flip-flops, the stored
datum _is_ real-valued.

See also: holdtrigger, holdnon0
=============================================================================*/

float calcflipflop(sndobj *me, sndobj *caller, int innr);

sndobj *flipflop( sndobj *trigger, sndobj *input )
{
  sndobj *p;
  
  p= newsndo( calcflipflop, "flipflop", "flipflop", input->nch, 2, trigger, input );
  p->skip= dontskip;
  return p;
}


float calcflipflop(sndobj *me, sndobj *caller, int innr)
{
  int ch;
  
  if( INCALC(0) ) {
    INCALC(1);
    FORCH
      OUTPUT(ch, INPUTC(1, ch));
  }
  else
    INSKIP(1);
  
  CALCRETURN;
}


/*=============================================================================
    holdtrigger <trigger> <duration>

Whenever <trigger> is !=0, the value is held for <duration> seconds, but at
least for one sample.

See also: flipflop, holdnon0
=============================================================================*/

struct holdtrigger {
  int	countdown;
};

float calcholdtrigger(sndobj *me, sndobj *, int);

sndobj *holdtrigger( sndobj *trigger, sndobj *duration )
{
  sndobj *p;
  struct holdtrigger *d;
  
  p= newsndo(calcholdtrigger, trigger->name, "holdtrigger", 1, 2, trigger, 
		    duration);
  p->skip= dontskip;
  p->private[0]= d= new(struct holdtrigger);
  d->countdown= 0;
  return p;
}

float calcholdtrigger(sndobj *me, sndobj *caller, int innr)
{
  struct holdtrigger *d;
  
  d= (struct holdtrigger *)me->private[0];
  if( INCALC(0) )
  {
    OUTPUT(0, INPUT(0));
    d->countdown= (int)round((double)SAMPRATE*(double)INCALC(1));
  }
  else if( d->countdown )
  {
    INSKIP(1);
    if( --d->countdown <= 0 )
      OUTPUT(0, 0);
  }
  else
    INSKIP(1);
  CALCRETURN;
}


/*=============================================================================
  holdnon0 (mode) <signal>

Holds the last non-zero input value.  The initial output value is 1.  If mode
is !=0, the output value is 1 rather than the last input value.

See also: holdtrigger, flipflop
=============================================================================*/

float calcholdnon0(sndobj *me, sndobj *, int);

sndobj *holdnon0(int mode, sndobj *signal)
{
  sndobj *p;
  int ch;
  
  p= newsndo( calcholdnon0, "holdnon0", "holdnon0", signal->nch, 1, signal);
  p->private[0]= new(int);
  *(int*)p->private[0]= mode;
  for( ch= 0; ch< p->nch; ++ch )
    p->ch[ch]= 1.0;
  return p;
}


float calcholdnon0(sndobj *me, sndobj *caller, int innr)
{
  int ch;
  
  INCALC(0);
  if( *(int*)me->private[0] ) {
    FORCH
      if( INPUTC(0, ch) )
        OUTPUT(ch, 1.0);
  }
  else
    FORCH
      if( INPUTC(0, ch) )
        OUTPUT(ch, INPUTC(0, ch));
  CALCRETURN;
}


/*=============================================================================
    timeout (delay) <trigger> <signal>

This object exits the program (delay) seconds after a non-zero value is
received from <trigger>.  <trigger> may be NULL, in which case the timeout
starts when evaluation is started.  <signal> is passed through unchanged and
may have any number of channels, while only the first channel of <trigger> is
watched.

Note that downstream delays or worse, asynchronous objects, will lead to
discrepancies between the given delay and the actual interval of sample values
output.  This object can only keep track of time by the number of samples which
are passed through itself.
=============================================================================*/
  
struct timeout {
  long countdown;
  int counting;
};

float calctimeout(sndobj *me, sndobj *caller, int innr);

sndobj *timeout( double delay, sndobj *trigger, sndobj *signal )
{
  sndobj *p;
  struct timeout *d;

  p= newsndo( calctimeout, "timeout", "timeout", signal->nch, 1, signal );
  p->skip= dontskip;
  if( trigger )
    addinput(p, trigger);
  p->private[0]= d= new(struct timeout);
  d->countdown= (long)floor(0.5 + delay*(double)SAMPRATE);
  d->counting= !trigger;
  return p;
}

float calctimeout(sndobj *me, sndobj *caller, int innr)
{
  struct timeout *d;
  int ch;

  d= (struct timeout *)me->private[0];
  if( me->nin > 1 )
    if( d->counting )
      INSKIP(1);
    else if( INCALC(1)!=0.0 )
      d->counting= 1;
  if( d->counting )
    if( --d->countdown < 0 )
      exit(0);
  INCALC(0);
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  CALCRETURN;
}


/*=============================================================================
    scheduler (mode) (nclients) <trigger>

This object is intended for distributing successive <trigger>s to (nclients)
objects, which will usually be of the same type.  Every time <trigger> is
non-zero, its value will be output on a different output channel.  The total
number of output channels equals <nclients>.  For now, the only (mode) is
SCHED_RR, round-robin scheduling, which causes each of the clients to be used
in turn.
=============================================================================*/

float calcscheduler(sndobj *me, sndobj *caller, int innr);

sndobj *scheduler(int mode, int nclients, sndobj *trigger)
{
  sndobj *p;

  if( nclients < 2 )
    return trigger;
  p= newsndo( calcscheduler, "scheduler", "scheduler", nclients, 1, trigger );
  p->private[0]= new(int);
  *(int*)p->private[0]= nclients-1;
  return p;
}

float calcscheduler(sndobj *me, sndobj *caller, int innr)
{
  int *outind;
  float trigval;

  outind= me->private[0];
  trigval= INCALC(0);
  OUTPUT(*outind, 0.0f);
  if( trigval != 0.0 )
  {
    if( ++*outind >= me->nch )
      *outind= 0;
    OUTPUT(*outind, trigval);
  }
  CALCRETURN;
}



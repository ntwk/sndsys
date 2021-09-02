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
#include <float.h>
#include <ctype.h>
#include <string.h>
#include <sys/time.h>
#ifndef NO_PTHREADS
#include <pthread.h>  // doesn't declare pthread_mutexattr_settype ?
extern int pthread_mutexattr_settype (pthread_mutexattr_t *__attr, int __kind);
#endif

#include "sndsys.h"

/*=============================================================================
    linear (offset) (slope) <signal>

Linear function of input.

See also: add, mul, arithmetic
=============================================================================*/

struct linear {
  char  namestr[NAMESTRLEN];
  float offset, slope;
};

float calclinear(sndobj *me, sndobj *, int);

sndobj *linear(float offset, float slope, sndobj *in)
{
  sndobj *p;
  struct linear *d;
  
  if( slope==0.0 ) {
    fprintf(stderr, "Warning: linear: zero slope, returning constant object.\n");
    return c(offset);
  }
  if( offset==0.0f && slope==1.0f )
    return in;
  d= (struct linear *)malloc(sizeof(struct linear));
  snprintf(d->namestr, NAMESTRLEN, "%g + %g x", offset, slope );
  p= newsndo(calclinear, d->namestr, "linear", in->nch, 1, in);
  p->private[0]= d;
  d->offset= offset;
  d->slope= slope;
  return p;
}

float calclinear(sndobj *me, sndobj *caller, int innr)
{
  struct linear *d;
  int ch;

  d= (struct linear *)me->private[0];
  OUTPUT(0, d->offset + d->slope*INCALC(0) );
  for( ch= 1; ch < me->nch; ++ch )
    OUTPUT(ch, d->offset + d->slope*INPUTC(0, ch) );
  CALCRETURN;
}


/*=============================================================================
    c (value)
    cc (number of channels) (value)

c - constant.  cc - multi-channel constant.

See also: extc
=============================================================================*/

float calcc(sndobj *me, sndobj *, int);

sndobj *c(float val)
{
  sndobj *p;
  char *valstr;
  
  valstr= malloc(NAMESTRLEN);
  snprintf(valstr, NAMESTRLEN, "%g", (double)val );
  p= newsndo(calcc, valstr, "c", 1, 0);
  p->private[0]= valstr;
  p->ch[0]= val;
  return p;
}

sndobj *cc(int nch, float val)
{
  sndobj *p;
  char *valstr;
  int ch;
  
  valstr= malloc(NAMESTRLEN);
  snprintf(valstr, NAMESTRLEN, "%g", (double)val );
  p= newsndo(calcc, valstr, "cc", nch, 0);
  p->private[0]= valstr;
  for( ch= 0; ch< nch; ++ch )
    p->ch[ch]= val;
  return p;
}

float calcc(sndobj *me, sndobj *caller, int innr)
{
  CALCRETURN;
}


/*=============================================================================
    extc (number of channels) [data pointer] [flag pointer]

External constant, for tuning of instruments or similar things depending on
user input.  data must contain nch floating-point numbers.  They are read when
creating this object and every time the flag is >0.  While the object is
reading the data, it is set to -1 (the writing function should wait during that
time), and to 0 after it has finished.  This object is not thread-safe in the
sense that the data may only be modified by the same thread which runs this
object.

See also: c, cc
=============================================================================*/

struct extc {
  float *data;
  int   *flag;
};

float calcextc(sndobj *me, sndobj *, int);

sndobj *extc( int nch, float *data, int *flag )
{
  sndobj *p;
  struct extc *d;
  int ch;
  
  p= newsndo( calcextc, "extc", "extc", nch, 0 );
  p->private[0]= d= new(struct extc);
  d->data= data;
  d->flag= flag;
  for( ch= 0; ch< nch; ++ch )
    p->ch[ch]= data[ch];
  return p;
}

float calcextc(sndobj *me, sndobj *caller, int innr)
{
  struct extc *d;
  int ch;
  
  d= (struct extc *)me->private[0];
  if( *d->flag > 0 )
  {
    *d->flag= -1;
    FORCH
      OUTPUT(ch, d->data[ch]);
    *d->flag= 0;
  }
  CALCRETURN;
}


/*=============================================================================
    add (n) <signal> ...

Adds up the number of inputs given by its first argument.  The number of output
channels is the maximum of the input channel numbers: inputs with fewer
channels wrap.  Additional inputs can be added with addinput() after
construction.  In that case the number of output channels may have to be
adjusted by hand.  As many sndobj's have the same number of channels as their
inputs, this object should be passed as an argument only after all inputs have
been added.  The cases n= 0 and 1 are not optimised in the obvious way (to
allow adding more inputs later), and the default number of output channels (for
n=0) is 1.

See also: linear, mul, arithmetic
=============================================================================*/

float calcadd(sndobj *me, sndobj *, int);

sndobj *add(int n, ...)
{
  sndobj *p;
  va_list in;
  int i;

  p= newsndo(calcadd, "add", "add", 1, 0);
  if(n>PLENTY) {
    fprintf(stderr, "Warning: add: Number of inputs (%d) limited to %d.\n", n, (int)PLENTY );
    n=PLENTY;
  }
  va_start(in, n);
  if( n ) {
    addinput(p, va_arg(in, sndobj*));
    p->nch= p->in[0]->nch;
  }
  for( i= 1; i< n; ++i ) {
    addinput(p, va_arg(in, sndobj*));
    if( p->in[i]->nch > p->nch )
      p->nch= p->in[i]->nch;
  }
  va_end(in);
  return p;
}


float calcadd(sndobj *me, sndobj *caller, int innr)
{
  float sum;
  int i, ch;
  
  sum= 0.0;
  for( i= 0; i< me->nin; ++i )
    sum += INCALC(i);
  OUTPUT(0, sum);
  for( ch= 1; ch< me->nch; ++ch ) {
    sum= 0.0;
    for( i= 0; i< me->nin; ++i )
      sum += INPUTC(i, ch % me->in[i]->nch);
    OUTPUT(ch, sum);
  }
  CALCRETURN;
}



/*=============================================================================
    mul (n) <signal> ...

Multiplies the number of inputs given by its first argument.  The number of
output channels is the maximum of the input channel numbers: inputs with fewer
channels wrap.  Additional inputs can be added with addinput() after
construction, in which case the number of output channels may have to be
adjusted by hand.  This object should be passed as an argument to downstream
objects only after all inputs have been added, as those downstream sndobj's may
set their number of channels acccording to it.  In the case of no inputs, 1 is
returned, and the number of output channels is 1.  The inputs are evaluated
from left to right.  If all channels of the product become zero at some point,
the remaining inputs are skipped.  Therefore it makes sense to pass the input
most likely to be zero as the first argument.

See also: linear, add, arithmetic
=============================================================================*/

float calcmul(sndobj *me, sndobj *, int);

sndobj *mul(int n, ...)
{
  sndobj *p;
  va_list in;
  int i;

  p= newsndo(calcmul, "mul", "mul", 1, 0);
  if( n<0 )
    n= 0;
  else if(n>PLENTY) {
    fprintf(stderr, "Warning: mul: Number of inputs (%d) limited to %d.\n", n, (int)PLENTY );
    n=PLENTY;
  }
  va_start(in, n);
  if( n ) {
    addinput(p, va_arg(in, sndobj*));
    p->nch= p->in[0]->nch;
  }
  for( i= 1; i< n; ++i ) {
    addinput(p, va_arg(in, sndobj*));
    if( p->in[i]->nch > p->nch )
      p->nch= p->in[i]->nch;
  }
  va_end(in);
  return p;
}


float calcmul(sndobj *me, sndobj *caller, int innr)
{
  int i, ch, zeroflag;
  
  if( me->nin==0 ) {
    OUTPUT(0, 1.0);
    CALCRETURN;
  }
  INCALC(0);
  zeroflag= 1;
  FORCH {
    me->ch[ch]= INPUTC(0, ch % me->in[0]->nch);
    if( me->ch[ch] )
      zeroflag= 0;
  }
  for( i= 1; i< me->nin && !zeroflag; ++i ) {
    INCALC(i);
    zeroflag= 1;
    FORCH {
      me->ch[ch] *= INPUTC(i, ch % me->in[i]->nch);
      if( me->ch[ch] )
        zeroflag= 0;
    }
  }
//  if( i< me->nin )
//    printf("mul: skipping inputs %d to %d.\n", i, me->nin-1);
  for( ; i< me->nin; ++i )
    INSKIP(i);
  CALCRETURN;
}



/*=============================================================================
    arithmetic "operation string" <argument1> <argument2>

The operation string may be: "+", "-", "_", "*", "/", "^", "|", "=", "\" (add,
subtract, absolute difference of abs. values, multiply, divide, to power,
absolute value of <argument1>, negative of <argument1>, inverse of <argument1>.
Remember to write "\" as "\\" in C code).  If a binary operator is given, the
argument with more channels determines the number of output channels; channels
of the other argument are wrapped.  For "^", the sign of the base is kept no
matter what the exponent.  If a unary operand is given, the second argument may
be passed as NULL.  For the inversion operation "\\" the second operand, if
given, is interpreted as an upper limit for the result.  Its channels wrap
around.

Usage tips: "-" and "_" are useful for automatically comparing sound streams.
"^"  with an exponent <1 can be used for soft clipping (cf. softclip).

See also: add, mul, linear
=============================================================================*/

struct arithmetic {
  char  op;
  char  namestr[NAMESTRLEN];
};

float calcarithmetic(sndobj *me, sndobj *, int);

sndobj *arithmetic(const char *op, sndobj *arg1, sndobj *arg2 )
{
  sndobj *p;
  struct arithmetic *d;
  
  d= (struct arithmetic *)malloc( sizeof(struct arithmetic) );
  d->op= *op;
  if( d->op!='+' && d->op!='-' && d->op!='_' && d->op!='*' && d->op!='/' && d->op!='^' && d->op!='|' && d->op!='=' && d->op!='\\' ) {
    fprintf(stderr, "Warning: Illegal operator argument (`%c') to arithmetic. Substituting `+'.\n", d->op );
    d->op= '+';
  }
  p= newsndo( calcarithmetic, d->namestr, "arithmetic", arg1->nch, 1, arg1 );
  p->private[0]= d;
  if( d->op=='+' || d->op=='-' || d->op=='_' || d->op=='*' || d->op=='/' || d->op=='^' ) {
    addinput(p, arg2);
    if( arg2->nch>arg1->nch )
      p->nch= arg2->nch;
    snprintf(d->namestr, NAMESTRLEN, "[%s] %c [%s]", arg1->name, d->op, arg2->name);
  }
  else {
    p->nin= 1;
    if( d->op=='|' )
      snprintf(d->namestr, NAMESTRLEN, "|%s|", arg1->name );
    else if( d->op=='=' )
      snprintf(d->namestr, NAMESTRLEN, "- %s", arg1->name );
    else {
      if( arg2 )
	addinput(p, arg2);
      snprintf(d->namestr, NAMESTRLEN, "1/ %s", arg1->name );
    }
  }
  return p;
}

float calcarithmetic(sndobj *me, sndobj *caller, int innr)
{
  struct arithmetic *d;
  double base, power;
  int ch, ch1, ch2, sign;
  
  d= (struct arithmetic *)me->private[0];
  INCALC(0);
  if( d->op!='|' && d->op!='=' && d->op!='\\' )
    INCALC(1);
  switch( d->op )
  {
    case '+':for( ch= 0, ch1= 0, ch2= 0; ch< me->nch; ++ch, ++ch1, ++ch2 ) {
		if( ch1>= me->in[0]->nch )
		  ch1= 0;
		if( ch2>= me->in[1]->nch )
		  ch2= 0;
		OUTPUT(ch, INPUTC(0, ch1) + INPUTC(1, ch2));
	     }
	    break;
    case '-':for( ch= 0, ch1= 0, ch2= 0; ch< me->nch; ++ch, ++ch1, ++ch2 ) {
		if( ch1>= me->in[0]->nch )
		  ch1= 0;
		if( ch2>= me->in[1]->nch )
		  ch2= 0;
		OUTPUT(ch, INPUTC(0, ch1) - INPUTC(1, ch2));
	     }
	    break;
    case '_':for( ch= 0, ch1= 0, ch2= 0; ch< me->nch; ++ch, ++ch1, ++ch2 ) {
		if( ch1>= me->in[0]->nch )
		  ch1= 0;
		if( ch2>= me->in[1]->nch )
		  ch2= 0;
		OUTPUT(ch, fabsf(fabsf(INPUTC(0, ch1)) - fabsf(INPUTC(1, ch2))));
	     }
	    break;
    case '*':for( ch= 0, ch1= 0, ch2= 0; ch< me->nch; ++ch, ++ch1, ++ch2 ) {
		if( ch1>= me->in[0]->nch )
		  ch1= 0;
		if( ch2>= me->in[1]->nch )
		  ch2= 0;
		OUTPUT(ch, INPUTC(0, ch1) * INPUTC(1, ch2));
	     }
	    break;
    case '/':for( ch= 0, ch1= 0, ch2= 0; ch< me->nch; ++ch, ++ch1, ++ch2 ) {
		if( ch1>= me->in[0]->nch )
		  ch1= 0;
		if( ch2>= me->in[1]->nch )
		  ch2= 0;
		if( INPUTC(1, ch2) )
		  OUTPUT(ch, INPUTC(0, ch1) / INPUTC(1, ch2));
		else if( INPUTC(0, ch1)>0.0 )
		  OUTPUT(ch, FLT_MAX);
		else if( INPUTC(0, ch1)<0.0 )
		  OUTPUT(ch, -FLT_MAX);
		else;
	    // leave previous value as limit of the quotient of two 0-series
	     }
	    break;
    case '^':for( ch= 0, ch1= 0, ch2= 0; ch< me->nch; ++ch, ++ch1, ++ch2 ) {
		if( ch1>= me->in[0]->nch )
		  ch1= 0;
		if( ch2>= me->in[1]->nch )
		  ch2= 0;
		base= INPUTC(0, ch1);
		sign= base>=0.0? 1:-1;
		base= fabs(base);
		if( !base )
		  OUTPUT(ch, INPUTC(1, ch2)>=0.0? 0.0: FLT_MAX);
		else {
		  power= exp(log(base) * INPUTC(1, ch2));
		  if( finite(power) )
		    OUTPUT(ch, sign*power );
		  else
		    OUTPUT(ch, sign*FLT_MAX);
		}
	      }
	    break;
    case '|':for( ch= 0; ch< me->nch; ++ch )
		OUTPUT(ch, fabsf( INPUTC(0, ch) ) );
	    break;
    case '=':for( ch= 0; ch< me->nch; ++ch )
		OUTPUT(ch, -INPUTC(0, ch));
	    break;
    case '\\':for( ch= 0; ch< me->nch; ++ch )
	       if( INPUTC(0, ch) )
		 OUTPUT(ch, 1.0f/INPUTC(0, ch));
	       else if( me->ch[ch]>=0 )
		 OUTPUT(ch, FLT_MAX);
	       else
		 OUTPUT(ch, -FLT_MAX);
	     if( me->nin > 1 ) {
	       int ch2;
	       INCALC(1);
	       for( ch= 0, ch2= 0; ch< me->nch; ++ch, ++ch2 ) {
		 if( ch2 >= me->in[1]->nch )
		   ch2= 0;
		 if( fabsf(me->ch[ch]) > fabsf(INPUTC(1, ch2)) )
		   if( me->ch[ch] > 0 )
		     me->ch[ch]= fabsf(INPUTC(1, ch2));
		   else
		     me->ch[ch]= -fabsf(INPUTC(1, ch2));
	       }
	     }
    default:break;
  }
  CALCRETURN;
}


/*=============================================================================
    ch (channel) <signal>
    c0 <signal>
    c1 <signal>
    c2 <signal>
    c3 <signal>

Channel selection.  ch extracts the channel which is given in its first
argument, the others the one which is part of their name.  If the channel does
not exist in the input, zero is returned.

See also: switchboard, chsel
=============================================================================*/

float calcch(sndobj *me, sndobj *, int);

sndobj *ch(int ch, sndobj *in)
{
  sndobj *p;
  char *typestr;
  
  typestr= malloc(NAMESTRLEN);
  snprintf(typestr, NAMESTRLEN, "c%d", ch );
  p= newsndo(calcch, in->name, typestr, 1, 1, in);
  p->private[0]= malloc(sizeof(int));
  p->private[1]= typestr;
  if( ch >= in->nch )
    fprintf( stderr, "Warning: channel for ch (%d) outside range (0..%d).\n", ch, in->nch );
  *(int*)p->private[0]= ch;
  return p;
}


sndobj *c0(sndobj *in)  {   return ch(0, in);   }
sndobj *c1(sndobj *in)  {   return ch(1, in);   }
sndobj *c2(sndobj *in)  {   return ch(2, in);   }
sndobj *c3(sndobj *in)  {   return ch(3, in);   }


float calcch(sndobj *me, sndobj *caller, int innr)
{
  INCALC(0);
  OUTPUT(0, INPUTC(0, *(int*)me->private[0]));
  CALCRETURN;
}


/*=============================================================================
     chsum <signal>

chsum forms the sum of all channels of the input signal and outputs a
one-channel signal.

See also: chavg, ch, mix
=============================================================================*/

float calcchsum(sndobj *me, sndobj *, int);

sndobj *chsum(sndobj *signal)
{
  return newsndo( calcchsum, "chavg", "chavg", 1, 1, signal );
}

float calcchsum(sndobj *me, sndobj *caller, int innr)
{
  int ch;

  INCALC(0);
  me->ch[0]= 0;
  for( ch= 0; ch < me->in[0]->nch; ++ch )
    me->ch[0] += INPUTC(0, ch);
  CALCRETURN;
}


/*=============================================================================
     chavg <signal>

chavg forms the average of all channels of the input signal and outputs a
one-channel signal.

See also: chsum, ch, mix
=============================================================================*/

float calcchavg(sndobj *me, sndobj *, int);

sndobj *chavg(sndobj *signal)
{
  return newsndo( calcchavg, "chavg", "chavg", 1, 1, signal );
}

float calcchavg(sndobj *me, sndobj *caller, int innr)
{
  int ch;

  INCALC(0);
  me->ch[0]= 0;
  for( ch= 0; ch < me->in[0]->nch; ++ch )
    me->ch[0] += INPUTC(0, ch);
  me->ch[0] /= me->in[0]->nch;
  CALCRETURN;
}


/*=============================================================================
    chsel (nch) <channels> <signal>

A routing object.  chsel has (nch) output channels.  Which channels of <signal>
go to which output channels is determined by the respective channel of
<channels>.  <signal> and <channels> have to have the same number of channels;
otherwise the surplus channels of whichever has more are discarded.  Output
channels which receive no contribution are silent; signal channels assigned to
the same output channel are added up.  Channel numbers start at 0 and are
rounded to the nearest integer away from zero.  Out-of-range channel numbers
are set to the minimum or maximum.

See also: switchboard, ch
=============================================================================*/

struct chsel {
  int nchin;
};

float calcchsel(sndobj *me, sndobj *, int);

sndobj *chsel(int nch, sndobj *channels, sndobj *signal)
{
  sndobj *p;
  struct chsel *d;

  if( nch < 1 ) {
    fprintf(stderr, "Warning: chsel: number of output channels has to be at least 1.  Setting to 1.\n");
    nch= 1;
  }
  else if( nch > PLENTY ) {
    fprintf(stderr, "Warning: chsel: number of output channels can be at most %d.  Setting to %d.\n", PLENTY, PLENTY);
    nch= PLENTY;
  }
  p= newsndo(calcchsel, signal->name, "chsel", nch, 2, channels, signal);
  p->private[0]= d= new(struct chsel);
  d->nchin= channels->nch > signal->nch? signal->nch : channels->nch;
  return p;
}

float calcchsel(sndobj *me, sndobj *caller, int innr)
{
  struct chsel *d;
  float outch;
  int ch, inch;

  d= (struct chsel *)me->private[0];
  INCALC(0);
  INCALC(1);
  FORCH
    OUTPUT(ch, 0);
  for( inch= 0; inch< d->nchin; ++inch ) {
    outch= INPUTC(0, inch);
    if( outch < 0 || !finitef(outch) )
      ch= 0;
    else if( outch > me->nch-1 )
      ch= me->nch - 1;
    else
      ch= (int)truncf(outch+0.5f);
    me->ch[ch] += INPUTC(1, inch);
  }
  CALCRETURN;
}


/*=============================================================================
    switchboard  "specification string"  <input> ...

Switching matrix, eg for creating a stereo signal from two mono signals.
Argument string gives the inputs to be connected to all outputs in the form
"<input>[.<channel>] <input2>[.<channel2>] ...".  Number of inputs is largest
input number plus 1.  Input channels default to 0.  Input = -1 means return 0
for this channel.

See also: chsel, ch
=============================================================================*/

float calcswitchboard(sndobj *me, sndobj *, int);

sndobj *switchboard( char *sboard, ... )
{
  sndobj *p;
  int *tab, *write;
  char *parse;
  va_list varg;
  int nin, in, nch, ch;
  
  p= newsndo( calcswitchboard, sboard, "switchboard", 1, 0 );
  p->private[0]= tab= (int*)malloc(PLENTY*2*sizeof(int));
  write= tab;
  nin= 0;
  nch= 0;
  parse= sboard;
  while( isspace(*parse) )
    ++parse;
  while( *parse )
  {
    if( !isdigit(*parse) && *parse!='-' ) {
      fprintf(stderr, "Warning: switchboard `%s': Illegal character `%c' at position %d.\n", sboard, *parse, (int)(parse-sboard) );
      break;
    }
    if( nch==PLENTY ) {
      fprintf(stderr, "Warning: switchboard `%s': Too many channels ( > %d).\n", sboard, PLENTY );
      break;
    }
    in= (int)strtol( parse, &parse, 10 );
    if( in<-1 || in>=PLENTY ) {
      fprintf(stderr, "Warning: switchboard `%s': Input number %d out of range in channel %d.\n", sboard, in, nch );
      in= -1;
    }
    if( *parse=='.' ) {
      ch= (int)strtol( parse+1, &parse, 10 );
      if( ch<-1 || ch>=PLENTY ) {
	fprintf(stderr, "Warning: switchboard `%s': Input channel %d out of range in channel %d.\n", sboard, ch, nch );
	ch= 0;
      }
    }
    else
      ch= 0;
    if( in+1>nin )
      nin= in+1;
    *write++ = in;
    *write++ = ch;
    ++nch;
    while( *parse==',' || *parse==';' || isspace(*parse) )
      ++parse;
  }
  p->nin= nin;
  p->nch= nch;
  va_start(varg, sboard);
  for( in= 0; in< nin; ++in ) {
    p->in[in]= va_arg(varg, sndobj*);
    p->inch[in]= p->in[in]->ch;
    for( ch= 0; ch< p->nch; ++ch )
      if( tab[2*ch]==in )
	break;
    if( ch==p->nch )
      fprintf(stderr, "Warning: switchboard `%s': Input %d not used.\n", sboard, in );
  }
  va_end(varg);
  for( ch= 0; ch< p->nch; ++ch )
    if( tab[2*ch+1]>=p->in[tab[2*ch]]->nch )
      fprintf(stderr, "Warning: switchboard `%s': Input channel (%d) >= input channel count (%d) in channel %d.\n", sboard, tab[2*ch+1], p->in[tab[2*ch]]->nch, ch );
  return p;
}


float calcswitchboard(sndobj *me, sndobj *caller, int innr)
{
  int *tab;
  int i, ch;
  
  for( i= 0; i< me->nin; ++i )
    INCALC(i);
  tab= (int*)me->private[0];
  for( ch= 0; ch < me->nch; ++ch, tab += 2 )
    if( *tab < 0 )
      OUTPUT(ch, 0.0);
    else
      OUTPUT(ch, INPUTC(*tab, tab[1]));
  CALCRETURN;
}


/*=============================================================================
    mix "specification string" <input> ...

Mixing with constant coefficients.  The inputs are preceded by a string giving
the coefficients, inputs and channels to use.  Its syntax is as follows:

mixing-string = ( channel-mode-string | stereo-mode-string )
channel-mode-string = output-channel-spec { ";" output-channel-spec } [";"]
output-channel-spec = channel-term-spec { "+" channel-term-spec }
channel-term-spec = inputnr "." channelnr [ "*" coefficient ]
stereo-mode-string = input-term-spec { "+" input-term-spec } [";"]
input-term-spec = inputnr [ "*" coefficient ]

inputnr and channelnr are integers, coefficients are floating-point numbers.
In "stereo" mode, the number of output channels is the largest channel number
of all inputs.  The channel number of narrower inputs is wrapped around so that
mono and stereo inputs can be mixed together to give stereo output.

    matrix (nch) [coeffs] ...

The low-level interface of mix.  Its arguments are the number of output
channels and a float array containing the matrix coefficients row by row,
concluded by END.  Each column corresponds to one channel of one input, in the
order they are given as arguments.  Note that this is dangerous to use since
you have to know exactly how many channels each input has, otherwise
coefficients will not match the intended channels.

See also: switchboard, chavg, chsum
=============================================================================*/

struct mix {
  int   stereomode, nterms;
  int   *inp;
  float *coeff;
};

float calcmix(sndobj *me, sndobj *, int);

sndobj *mix( char *mixstr, ... )
{
  sndobj *p;
  struct mix *d;
  float *cwrite;
  int *iwrite;
  char *name, *parse, *endconv;
  va_list varg;
  float co;
  int t, nin, in, inch, ch, nplus, nsemi;
  
  if( strlen(mixstr) < NAMESTRLEN )
    name= mixstr;
  else {
    name= news(char, NAMESTRLEN);
    strncpy(name, mixstr, NAMESTRLEN);
    name[NAMESTRLEN-1]= 0;
    name[NAMESTRLEN-2]= '.';
    name[NAMESTRLEN-3]= '.';
  }
  p= newsndo( calcmix, name, "mix", 1, 0 );
  if( name != mixstr )
    p->private[1]= name;
  nplus= nsemi= 0;
  for( parse= mixstr; *parse; ++parse )
    if( *parse==';' )
      ++nsemi;
    else if( *parse=='+' )
      ++nplus;
  while( --parse > mixstr )
    if( !isspace(*parse) ) {
      if( *parse!=';' )
	++nsemi;
      break;
    }
  if( nsemi+nplus == 0 ) {
    fprintf(stderr, "Warning: mix: Empty mixing string. Returning c(0).\n");
    return c(0);
  }
  if( nsemi > PLENTY ) {
    fprintf(stderr, "Warning: mix `%s': Too many output channels (%d, max=%d).  Discarding excess.\n", p->name, nsemi, PLENTY);
    p->nch= PLENTY;
  }
  else p->nch= nsemi;
  p->private[0]= d= new(struct mix);
  if( nsemi == 1 )
    d->stereomode= -1;
  else
    d->stereomode= 0;
  d->nterms= nsemi+nplus;
  p->private[2]= d->inp= iwrite= news(int, 3*d->nterms);
  p->private[3]= d->coeff= cwrite= news(float, d->nterms);
  parse= mixstr;
  nin= 0;
  for( ch= 0; ch< p->nch; ++ch )
  {
    do {
      while( *parse && (isspace(*parse) || *parse=='+') )
	++parse;
      if( !isdigit(*parse) ) {
	fprintf(stderr, "Warning: mix `%s': Illegal character `%c' at "
			"position %d.  Input number expected.\n", 
			p->name, *parse, (int)(parse-mixstr) );
	break;
      }
      in= (int)strtol( parse, &endconv, 10 );
      if( endconv==parse ) {
	fprintf(stderr, "Warning: mix `%s': Could not read input number "
	    "at position %d. Sibstituting 0.\n", p->name, (int)(parse-mixstr) );
	in= 0;
      }
      else if( in>=PLENTY ) {
	fprintf(stderr, "Warning: mix `%s': Input number %d out of range "
			"in channel %d.\n", p->name, in, ch );
	in= 0;
      }
      parse= endconv;
      if( *parse=='.' ) {
	inch= (int)strtol( ++parse, &endconv, 10 );
	if( d->stereomode > 0 ) {
	  fprintf( stderr, "Warning: mix `%s': Mixing of single channel "
		"(position %d) and multichannel (previously) not allowed.\n", 
		p->name, (int)(parse-mixstr));
	}
	else {
	  d->stereomode= 0;
	  if( endconv==parse ) {
	    fprintf(stderr, "Warning: mix `%s': Could not read input "
    "number at position %d. Substituting 0.\n", p->name, (int)(parse-mixstr));
	    inch= 0;
	  }
	  else if( inch>=PLENTY || inch < 0 ) {
	    fprintf(stderr, "Warning: mix `%s': Input channel %d out of "
			    "range in channel %d.\n", p->name, inch, ch );
	    inch= 0;
	  }
	}
	parse= endconv;
      }
      else {
	if( d->stereomode == 0 ) {
	  fprintf( stderr, "Warning: mix `%s': Missing channel number at "
		"position %d. It must be given if multiple output channels are"
		" defined or if it has been specified in other terms. "
		"Substituting 0.\n", p->name, (int)(parse-mixstr) );
	  inch= 0;
	}
	d->stereomode= 1;
	inch= 0;
      }
      while( *parse && isspace(*parse) )
	++parse;
      if( *parse=='*' ) {
	co= strtof( ++parse, &endconv );
	if( !finite(co) || parse==endconv ) {
	  fprintf( stderr, "Warning: mix `%s': Can't read coefficient at "
	    "position %d. Substituting 1.0.\n", p->name, (int)(parse-mixstr) );
	  co= 1.0;
	}
	else if( co == 0.0 )
	  fprintf( stderr, "Warning: mix `%s': Zero coefficient at "
		"postition %d.\n", p->name, (int)(parse-mixstr) );
	parse= endconv;
      }
      else
	co= 1.0;
      if( in+1>nin )
	nin= in+1;
      *iwrite++ = in;
      *iwrite++ = inch;
      *iwrite++ = ch;
      *cwrite++ = co;
      while( *parse && isspace(*parse) )
	++parse;
    }
    while( *parse=='+' );
    if( *parse && *parse != ';' ) {
      fprintf(stderr, "Warning: mix `%s': Illegal character `%c' at "
"position %d.  Expecting `;' or `+'.\n", p->name, *parse, (int)(parse-mixstr) );
    }
    ++parse;
  }
  va_start(varg, mixstr);
  for( in= 0; in< nin; ++in ) {
    addinput(p, va_arg(varg, sndobj*));
    if( d->stereomode && p->in[in]->nch > p->nch )
      p->nch= p->in[in]->nch;
    for( t= 0; t< d->nterms; ++t )
      if( d->inp[3*t]==in )
	break;
    if( t==d->nterms )
      fprintf(stderr, "Warning: mix `%s': Input %d not used.\n", p->name, in );
  }
  va_end(varg);
  for( t= 0; t< d->nterms; ++t )
    if( d->inp[3*t+1] >= p->in[d->inp[3*t]]->nch )
      fprintf(stderr, "Warning: mix `%s': Input channel (%d) >= "
	    "input channel count (%d) of input %d in term %d.\n", p->name, 
	    d->inp[3*t+1], p->in[d->inp[3*t]]->nch, d->inp[3*t], t );
  return p;
}


sndobj *matrix( int nch, float *coeff, ... )
{
  sndobj *p;
  struct mix *d;
  va_list varg;
  float *readco, *cwrite;
  int *iwrite;
  int t, nin, ch, in, ich;

  p= newsndo( calcmix, "matrix", "mix", nch, 0);
  p->private[0]= d= new(struct mix);
  d->stereomode= 0;
  for( readco= coeff; *readco!=END; ++readco );
  d->nterms= readco-coeff;
  nin= d->nterms / p->nch;
  if( nin*p->nch != d->nterms ) {
    fprintf(stderr, "Warning: matrix: Number of coefficients no multiple of output channels.\n");
    d->nterms= nin*p->nch;
  }
  p->private[2]= d->inp= iwrite= news(int, 3*d->nterms);
  p->private[3]= d->coeff= cwrite= news(float, d->nterms);
  va_start(varg, coeff);
  addinput(p, va_arg(varg, sndobj*));
  in= 0;
  ich= 0;
  for( t= 0, ch= 0; t< d->nterms; ++t, ++ch, ++ich ) {
    if( ich >= p->in[in]->nch ) {
      addinput(p, va_arg(varg, sndobj*));
      ++in;
      ich= 0;
    }
    if( ch >= p->nch )
      ch= 0;
    *iwrite++ = in;
    *iwrite++ = ich;
    *iwrite++ = ch;
    *cwrite++ = *readco++;
  }
  va_end(varg);
  return p;
}


float calcmix(sndobj *me, sndobj *caller, int innr)
{
  struct mix *d;
  float *readco;
  int *readterm;
  int i, ch, ich;
  
  for( i= 0; i< me->nin; ++i )
    INCALC(i);
  FORCH
    OUTPUT(ch, 0.0);
  d= (struct mix *)me->private[0];
  readterm= d->inp;
  readco= d->coeff;
  if( d->stereomode )
    for( i= 0; i< d->nterms; ++i ) {
      for( ch= ich= 0; ch< me->nch; ++ch ) {
	me->ch[ch] += *readco * INPUTC(readterm[0], ich);
	if( ++ich >= me->in[readterm[0]]->nch )
	  ich= 0;
      }
      readterm += 3;
      ++readco;
    }
  else
    for( i= 0; i< d->nterms; ++i ) {
      me->ch[readterm[2]] += *readco * INPUTC(readterm[0], readterm[1]);
      readterm += 3;
      ++readco;
    }
  CALCRETURN;
}


/*=============================================================================
    mix2 <gate> <input1> <input2>

Mixing of two inputs depending on <gate>.  1 gives <input1>, 0 <input2>.  The
output has the maximum number of channels of <input1> and <input2>.  The input
with the fewer channels wraps around.

See also: mix, chavg, chsum
=============================================================================*/

float calcmix2(sndobj *me, sndobj *caller, int innr);

sndobj *mix2( sndobj *gate, sndobj *input1, sndobj *input2 )
{
  return newsndo( calcmix2, "mix2", "mix2", 
		    input1->nch> input2->nch? input1->nch: input2->nch, 
		    3, gate, input1, input2 );
}

float calcmix2(sndobj *me, sndobj *caller, int innr)
{
  float gate;
  int ch, ch1, ch2;
  
  gate= INCALC(0);
  if( gate<=0.0 ) {
    INSKIP(1);
    INCALC(2);
    ch2= 0;
    FORCH {
      if( ch2 >= me->in[2]->nch )
	ch2= 0;
      OUTPUT(ch, INPUTC(2, ch2));
      ++ch2;
    }
    CALCRETURN;
  }
  else if( gate>=1.0 ) {
    INCALC(1);
    INSKIP(2);
    ch1= 0;
    FORCH {
      if( ch1 >= me->in[1]->nch )
	ch1= 0;
      OUTPUT(ch, INPUTC(1, ch1));
      ++ch1;
    }
    CALCRETURN;
  }
  else {
    INCALC(1);
    INCALC(2);
    ch1= ch2= 0;
    FORCH {
      if( ch1 >= me->in[1]->nch )
	ch1= 0;
      if( ch2 >= me->in[2]->nch )
	ch2= 0;
      OUTPUT(ch, INPUTC(2, ch2) + gate*(INPUTC(1, ch1)-INPUTC(2, ch2)));
      ++ch1;
      ++ch2;
    }
    CALCRETURN;
  }
}


/*=============================================================================
    cdelay (time interval) <signal>
    sdelay (number of samples) <signal>

Constant delay, positive or negative.  cdelay takes its delay time in secs,
maximal value +/-10 sec (sanity test).  sdelay takes the number of samples as
its argument.

See also: delay
=============================================================================*/

struct cdelay {
  long	dsamples;
  int	buffered;
};

float calccdelay(sndobj *me, sndobj *, int);

#define CDELAY_MAX      10.0
sndobj *cdelay(double time, sndobj *in)
{
  sndobj *p;
  struct cdelay *d;
  char *namestr;

  namestr= malloc(NAMESTRLEN);
  snprintf(namestr, NAMESTRLEN, "%g sec", time );
  p= newsndo(calccdelay, namestr, "cdelay", in->nch, 1, in);
  p->private[0]= d= malloc(sizeof(struct cdelay));
  p->private[1]= namestr;
  if( time<-CDELAY_MAX ) {
    fprintf(stderr, "Warning: cdelay: Delay (%g sec) too negative. Changed to maximum: %g sec.\n", time, -CDELAY_MAX );
    time= -CDELAY_MAX;
    snprintf(namestr, NAMESTRLEN, "%g sec", -CDELAY_MAX );
  }
  else if( time>CDELAY_MAX ) {
    fprintf(stderr, "Warning: cdelay: Delay (%g sec) too large. Changed to maximum: %g sec.\n", time, CDELAY_MAX );
    time= CDELAY_MAX;
    snprintf(namestr, NAMESTRLEN, "%g sec", CDELAY_MAX );
  }
  d->buffered= 1;
  if( time > HORIZON )
    buf( time, 0.0, 0.0, p->nch, p, 0 );
  else if( -time > HORIZON )
    buf( 0.0, -time, 0.0, p->nch, p, 0 );
  else
    d->buffered= 0;
  d->dsamples= (long)round(time*SAMPRATE);
  return p;
}


sndobj *sdelay(long samples, sndobj *in)
{
  sndobj *p;
  struct cdelay *d;
  char *namestr;

  if( !samples ) 
    return in;
  namestr= malloc(NAMESTRLEN);
  snprintf(namestr, NAMESTRLEN, "%ld samples", samples );
  p= newsndo(calccdelay, namestr, "sdelay", in->nch, 1, in);
  p->private[0]= d= malloc(sizeof(struct cdelay));
  p->private[1]= namestr;
  if( samples<(long)(-CDELAY_MAX*SAMPRATE) ) {
    fprintf(stderr, "Warning: sdelay: Delay (%ld samples) too negative. Changed to maximum: %ld samples.\n", samples, (long)(-CDELAY_MAX*SAMPRATE) );
    samples= (long)(-CDELAY_MAX*SAMPRATE);
    snprintf(namestr, NAMESTRLEN, "%ld samples", samples );
  }
  else if( samples>(long)(CDELAY_MAX*SAMPRATE) ) {
    fprintf(stderr, "Warning: sdelay: Delay (%ld samples) too large. Changed to maximum: %ld samples.\n", samples, (long)(CDELAY_MAX*SAMPRATE) );
    samples= (long)(CDELAY_MAX*SAMPRATE);
    snprintf(namestr, NAMESTRLEN, "%ld samples", samples );
  }
  d->buffered= 1;
  if( samples > (long)(HORIZON*SAMPRATE) )
    buf( (double)samples/(double)SAMPRATE, 0.0, 0.0, p->nch, p, 0 );
  else if( -samples > (long)(HORIZON*SAMPRATE) )
    buf( 0.0, (double)samples/(double)SAMPRATE, 0.0, p->nch, p, 0 );
  else
    d->buffered= 0;
  d->dsamples= samples;
  return p;
}


float calccdelay(sndobj *me, sndobj *caller, int innr)
{
  struct cdelay *d;
  long i;
  int ch;

  d= (struct cdelay *)me->private[0];
  if( !d->buffered ) {
  //  printf("calccdelay: %ld to go \n", *(long*)me->private[0] );
    if( d->dsamples )
    {
      if( d->dsamples < 0 )
      {
	for( i= -d->dsamples; i>0; --i )
	  INSKIP(0);
	d->dsamples= 0;
      }
      else {
	--d->dsamples;
	CALCRETURN;
      }
    }
    INCALC(0);
  }
  else {
    BUFINCALC(0, -d->dsamples, 0.0);
    BUFINCR(0, 1);
  }
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  CALCRETURN;
}


/*=============================================================================
    delay (maximal delay) <delay input> <signal>

Delay second input by delay time given by first.  Delay time is in secs,
maximal value is an argument (has to be <= 30 secs, sanity check).  Delay may
be negative up to HORIZON.  May be used for phase modulation.  If maximal delay
passed is negative, the maximal delay is its modulus, and the actual delay
value is rounded to a whole number of samples.  Otherwise, linear interpolation
is performed.

See also: cdelay, sdelay, rdelay
=============================================================================*/

struct delay {
  char  namestr[NAMESTRLEN];
  float maxdelay;
  int	nointerpol;
};

float calcdelay(sndobj *me, sndobj *, int);
void skipdelay(sndobj *me, sndobj *, int);

#define DELAY_MAXMAX    30.0

sndobj *delay( float maxdelay, sndobj *delay, sndobj *in )
{
  sndobj *p;
  struct delay *d;
  int nointer;
  
  nointer= maxdelay < 0;
  maxdelay= fabsf(maxdelay);
  if( maxdelay > DELAY_MAXMAX ) {
    maxdelay= DELAY_MAXMAX;
    fprintf( stderr, "Warning: Argument maxdelay to delay (%g) exceeds safety limit (%g).  Argument adjusted.\n", maxdelay, (double)DELAY_MAXMAX );
  }
  d= (struct delay *)malloc(sizeof(struct delay));
  snprintf(d->namestr, NAMESTRLEN, "%s << %s", in->name, delay->name );
  d->maxdelay= maxdelay;
  d->nointerpol= nointer;
  p= newsndo( calcdelay, d->namestr, "delay", in->nch, 2, delay, in );
  p->skip= skipdelay;
  p->private[0]= d;
  buf( maxdelay, HORIZON, 0.0, p->nch, p, 1 );
  return p;
}


float calcdelay(sndobj *me, sndobj *caller, int innr)
{
  struct delay *d;
  float delay;
  int ch, index;
  
  d= (struct delay *)me->private[0];
  delay= INCALC(0);
  if( delay>d->maxdelay )
    delay= d->maxdelay;
  else if( delay<-HORIZON )
    delay= -HORIZON;
  delay *= -SAMPRATE;
  if( !d->nointerpol ) {
    index= (int)floorf(delay);
    delay -= index;
  }
  else {
    index= (int)floorf(delay+0.5f);
    delay= 0.0;
  }
  BUFINCALC(1, index, delay);
  for( ch= 0; ch< me->nch; ++ch )
    OUTPUT(ch, INPUTC(1, ch));
  BUFINCR(1, 1);
  CALCRETURN;
}


void skipdelay(sndobj *me, sndobj *caller, int innr)
{
  INSKIP(0);
  BUFINCR(1, 1);
}



/*=============================================================================
    rdelay (maximal delay) <delay input> <reset> <signal>

Delay object which can be reset by a non-zero value on the <reset> input.  The
delay line is cleared, and only zeros are output for the delay which was
current when the reset happened.  Otherwise the same as the delay object.

See also: delay
=============================================================================*/

struct rdelay {
  char  namestr[NAMESTRLEN];
  float maxdelay, deadtime;
};

float calcrdelay(sndobj *me, sndobj *, int);
void skiprdelay(sndobj *me, sndobj *, int);

#define RDELAY_MAXMAX    30.0

sndobj *rdelay( float maxdelay, sndobj *delay, sndobj *reset, sndobj *in )
{
  sndobj *p;
  struct rdelay *d;
  
  maxdelay= fabsf(maxdelay);
  if( maxdelay > RDELAY_MAXMAX ) {
    maxdelay= RDELAY_MAXMAX;
    fprintf( stderr, "Warning: Argument maxdelay to rdelay (%g) exceeds safety limit (%g).  Argument adjusted.\n", maxdelay, (double)RDELAY_MAXMAX );
  }
  d= (struct rdelay *)malloc(sizeof(struct rdelay));
  snprintf(d->namestr, NAMESTRLEN, "%s << %s", in->name, delay->name );
  d->maxdelay= maxdelay;
  d->deadtime= 0.0f;
  p= newsndo( calcrdelay, d->namestr, "rdelay", in->nch, 3, delay, reset, in );
  p->skip= skiprdelay;
  p->private[0]= d;
  buf( maxdelay, HORIZON, 0.0, p->nch, p, 2 );
  return p;
}

float calcrdelay(sndobj *me, sndobj *caller, int innr)
{
  struct rdelay *d;
  float delay;
  int ch, index;

  d= (struct rdelay *)me->private[0];
  delay= INCALC(0);
  if( delay>d->maxdelay )
    delay= d->maxdelay;
  else if( delay<-HORIZON )
    delay= -HORIZON;
  if( INCALC(1) ) {
    d->deadtime= delay - 1.0/SAMPRATE;
    FORCH
      OUTPUT(ch, 0.0);
  }
  else if( d->deadtime >= 0.5/SAMPRATE ) {
    d->deadtime -= 1.0/SAMPRATE;
    FORCH
      OUTPUT(ch, 0.0);
  }
  else
  {
    delay *= -SAMPRATE;
    index= (int)floorf(delay);
    delay -= index;
    BUFINCALC(2, index, delay);
    FORCH
      OUTPUT(ch, INPUTC(2, ch));
  }
  BUFINCR(2, 1);
  CALCRETURN;
}

void skiprdelay(sndobj *me, sndobj *caller, int innr)
{
  struct rdelay *d;

  d= (struct rdelay *)me->private[0];
  if( d->deadtime >= 0.5/SAMPRATE )
    d->deadtime -= 1.0/SAMPRATE;
  INSKIP(0);
  INSKIP(1);
  BUFINCR(2, 1);
}


/*=============================================================================
    limiter (requested maximal value) <signal>

Limits output amplitude to max by rescaling all samples starting from the
offending one.  This is a kludge intended for developing sounds only.  It can
save you rerunning a program when the resulting sound is out of range.
=============================================================================*/

struct limiter {
  float max, sigmin, sigmax, scale;
};

float calclimiter(sndobj *me, sndobj *, int);
void exitlimiter(sndobj *me);

sndobj *limiter( float max, sndobj *signal )
{
  sndobj *p;
  struct limiter *d;

  p= newsndo( calclimiter, signal->name, "limiter", signal->nch, 1, signal );
  p->exit= exitlimiter;
  p->private[0]= d= new(struct limiter);
  d->max= max;
  d->sigmax= -FLT_MAX;
  d->sigmin= FLT_MAX;
  d->scale= 1.0;
  return p;
}

float calclimiter(sndobj *me, sndobj *caller, int innr)
{
  struct limiter *d;
  float val;
  int ch;
  
  d= (struct limiter *)me->private[0];
  INCALC(0);
  for( ch= 0; ch < me->nch; ++ch ) {
    val= INPUTC(0, ch);
    if( val*d->scale > d->max )
      d->scale= d->max/val;
    else if( val*d->scale < -d->max )
      d->scale= -d->max/val;
    if( val > d->sigmax )
      d->sigmax= val;
    if( val < d->sigmin )
      d->sigmin= val;
  }
  for( ch= 0; ch < me->nch; ++ch )
    OUTPUT(ch, d->scale*INPUTC(0, ch));
  CALCRETURN;
}

void exitlimiter(sndobj *me)
{
  struct limiter *d;

  d= (struct limiter *)me->private[0];
  printf( "Limiter of `%s', type `%s'.  Original min= %g, max= %g.  "
	    "Final scale: %g\n", me->in[0]->name, me->in[0]->type, 
			d->sigmin, d->sigmax, d->scale );
}



/*=============================================================================
    comptime <subtree>

Times the amount of time spent calculating in a given subtree.  For this to
work, the subtree must have no fanout connections to any other part of the main
tree, and no other process should be running which takes any significant amount
of CPU time.  For more reliable results, time a constant along with the subtree
you are interested in and compare.
=============================================================================*/

struct comptime {
  long  secs, usecs, samples, skipsamples;
};

float calccomptime(sndobj *me, sndobj *, int);
void skipcomptime(sndobj *me, sndobj *, int);
void exitcomptime(sndobj *me);

sndobj *comptime(sndobj *subtree)
{
  sndobj *p;
  struct comptime *d;
  
  p= newsndo( calccomptime, subtree->name, "comptime", subtree->nch, 1, subtree );
  p->skip= skipcomptime;
  p->exit= exitcomptime;
  p->private[0]= d= new(struct comptime);
  d->secs= d->usecs= d->samples= d->skipsamples= 0;
  return p;
}

float calccomptime(sndobj *me, sndobj *caller, int innr)
{
  struct comptime *d;
  struct timeval before, after;
  int ch;
  
  d= (struct comptime *)me->private[0];
  gettimeofday(&before, NULL);
  INCALC(0);
  gettimeofday(&after, NULL);
  d->secs += after.tv_sec - before.tv_sec;
  d->usecs += after.tv_usec - before.tv_usec;
  if( d->usecs< 0 ) {
    --d->secs;
    d->usecs += 1000000;
  }
  else if( d->usecs > 1000000 ) {
    ++d->secs;
    d->usecs -= 1000000;
  }
  ++d->samples;
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  CALCRETURN;
}

void skipcomptime(sndobj *me, sndobj *caller, int innr)
{
  struct comptime *d;
  struct timeval before, after;
  
  d= (struct comptime *)me->private[0];
  gettimeofday(&before, NULL);
  INSKIP(0);
  gettimeofday(&after, NULL);
  d->secs += after.tv_sec - before.tv_sec;
  d->usecs += after.tv_usec - before.tv_usec;
  if( d->usecs< 0 ) {
    --d->secs;
    d->usecs += 1000000;
  }
  else if( d->usecs > 1000000 ) {
    ++d->secs;
    d->usecs -= 1000000;
  }
  ++d->samples;
  ++d->skipsamples;
  return;
}

void exitcomptime(sndobj *me)
{
  struct comptime *d;
  long sampersec;

  d= (struct comptime *)me->private[0];
  sampersec= (long)round((double)d->samples/((double)d->secs+1e-6*(double)d->usecs));
  printf("comptime: Time taken by object `%s', type `%s', and its subtree:\n",
	    me->in[0]->name, me->in[0]->type );
  if( d->secs>3600 ) {
    printf(" %ld hours", d->secs/3600);
    d->secs %= 3600;
    if( d->secs || d->usecs ) {
      printf(" %ld minutes", d->secs/60 );
      d->secs %= 60;
      if( d->secs || d->usecs )
	printf(" %ld.%03ld seconds", d->secs, d->usecs );
    }
  }
  else if( d->secs>60 ) {
    printf(" %ld minutes", d->secs/60 );
    d->secs %= 60;
    if( d->secs || d->usecs )
      printf(" %ld.%03ld seconds", d->secs, d->usecs );
  }
  else
    printf(" %ld.%03ld seconds", d->secs, (d->usecs+500)/1000 );
  printf(" for %ld samples (%g seconds), %ld skipped; %ld samples per second.\n", d->samples, (double)d->samples/(double)SAMPRATE, d->skipsamples, sampersec );
  
}


/*=============================================================================
    diff <signal>

Differentiates a signal by subtracting a sample from the previous one.

See also: integrate
=============================================================================*/

float calcdiff(sndobj *me, sndobj *, int);

sndobj *diff(sndobj *in)
{
  sndobj *p;
  float *prev;
  int i;
  
  p= newsndo( calcdiff, in->name, "diff", in->nch, 1, in);
  p->private[0]= prev= news(float, PLENTY);
  for( i= 0; i< PLENTY; ++i )
    prev[i]= 0.0;
  return p;
}

float calcdiff(sndobj *me, sndobj *caller, int innr)
{
  float *prev;
  int ch;
  
  prev= (float*)me->private[0];
  INCALC(0);
  FORCH {
    OUTPUT(ch, INPUTC(0, ch) - prev[ch]);
//printf("diff: channel %d, previous %g, current %g\n", ch, prev[ch], INPUTC(0, ch));
    prev[ch]= INPUTC(0, ch);
  }
  CALCRETURN;
}


/*=============================================================================
    integrate <signal>

Integrates a signal by adding up sample values.

See also: diff
=============================================================================*/

float calcintegrate(sndobj *me, sndobj *, int);

sndobj *integrate(sndobj *in)
{
  sndobj *p;
  double *prev;
  int i;
  
  p= newsndo( calcintegrate, in->name, "integrate", in->nch, 1, in);
  p->private[0]= prev= news(double, PLENTY);
  for( i= 0; i< PLENTY; ++i )
    prev[i]= 0.0;
  return p;
}

float calcintegrate(sndobj *me, sndobj *caller, int innr)
{
  double *prev;
  int ch;
  
  prev= (double*)me->private[0];
  INCALC(0);
  FORCH {
    prev[ch] += INPUTC(0, ch);
    OUTPUT(ch, prev[ch]);
  }
  CALCRETURN;
}


/*=============================================================================
    bifurcation (interval) <Verhulst parameter input>

This object can illustrate the period doubling of the Verhulst process x ->
c*x*(1-x).  The output changes only every (interval) seconds.  After each
interval, the cycle (if any) of the Verhulst process is computed by performing
1000 iterations.  In a further 200 iterations, the element of the cycle which
is closest to the previously output value is searched for.  The next element of
the cycle becomes the next value to be output.  This allows to hear each cycle,
even if the parameter c changes slowly.

The parameter c is given by the input and should be between 2.5 and
4.  The result is always between 0 and 1.  Period doubling occurs around c= 3,
sqrt(6)=3.449, 3.5, ... up to 3.57, where chaotic behaviour starts.  At
c=3.8284, there is a stable 3-cycle, which doubles infinitely many times until
c=3.8495.  By ramping the parameter input up through these ranges, one can
actually hear the period doubling behaviour.
=============================================================================*/

struct bifurcation {
  int       interval, countdown;
  float     oldparam;
};

float calcbifurcation(sndobj *me, sndobj *, int);

sndobj *bifurcation( float interval, sndobj *parameter )
{
  sndobj *p;
  struct bifurcation *d;
  
  p= newsndo(calcbifurcation, "bifurcation", "bifurcation", 1, 1, parameter);
  p->skip= dontskip;
  p->private[0]= d= new(struct bifurcation);
  d->interval= (int)round(interval*SAMPRATE);
  d->countdown= -1;
  d->oldparam= 0;
  return p;
}

#define BIFURC_PREPERIOD    1000
#define BIFURC_SEARCHCLOSEST 200

float calcbifurcation(sndobj *me, sndobj *caller, int innr)
{
  struct bifurcation *d;
  float x, param, closest;
  int count;
  
  d= (struct bifurcation *)me->private[0];
  param= INCALC(0);
  if( d->countdown <= 0 )
  {
    if( param> 4.0 )
      param= 4.0;
    else if( param< 0.0 )
      param= 0.0;
    if( param==d->oldparam )
      closest= me->ch[0];
    else {
      x= 0.5;
      for( count= BIFURC_PREPERIOD; count> 0; --count )
	x= x*(1-x)*param;
      closest= 2.0;
      for( count= BIFURC_SEARCHCLOSEST; count> 0; --count ) {
	x= x*(1-x)*param;
	if( fabsf(x-me->ch[0]) < fabsf(closest-me->ch[0]) )
	  closest= x;
      }
    }
    OUTPUT(0, closest*(1-closest)*param);
    d->countdown= d->interval;
  }
  else
    --d->countdown;
  CALCRETURN;
}


/*=============================================================================
    smoothstep (mode) (minimum step) <steptime> <signal>

Smooths steps in the input.  Depending on (mode), the result is linear
(SMST_LINEAR) or has the shape of half a sine oscillation (between -pi/2 and
pi/2; SMST_SINE).  The mode can be OR combined with flags which restrict the
type of steps to smooth: SMST_UP and SMST_DOWN cause only steps which go up or
down to be smoothed, SMST_ABSUP and SMST_ABSDOWN only steps which reduce the
modulus of the input by <minstep>.  Multiple flags may be OR combined; if none
are given, all steps qualify.  <minstep> is the minimum step height to be
recognised as a step, <steptime> is the time for the smooth ramping of each
step.  While one ramp is in progress, occurring steps are not recognised, so
the steps have to be at least <steptime> apart.  This limitation can be worked
around by using several smoothstep objects in a row, which handle steps in
turn.  Only the first channel of the signal <in> is used.  The second output
channel contains the difference between the current ramp and the starting or
finishing value, whichever is closer, or zero if no ramp is in progress.

See also: smoothstart, avg, lowpass
=============================================================================*/

struct smoothstep {
  int       mode, flags, init_done, countdown;
  float     minstep, curr, inc, start, goal;
  float     *sinelut;
};

float calcsmoothstep(sndobj *me, sndobj *, int);

#define SMST_SINESIZE   SAMPRATE

sndobj *smoothstep( int mode, float minstep, sndobj *steptime, sndobj *in )
{
  sndobj *p;
  struct smoothstep *d;
  
  p= newsndo( calcsmoothstep, "smoothstep", "smoothstep", 2, 2, steptime, in );
  p->private[0]= d= new(struct smoothstep);
  d->mode= mode & SMST_MODEMASK;
  d->flags= mode & SMST_FLAGMASK;
  d->init_done= 0;
  d->countdown= 0;
  d->minstep= minstep;
  if( d->mode==SMST_SINE )
    p->private[2]= d->sinelut= sinelookup(SMST_SINESIZE);
  else
    d->sinelut= NULL;
  buf( 0.0, HORIZON, 0.0, 1, p, 1 );
  return p;
}


float calcsmoothstep(sndobj *me, sndobj *caller, int innr)
{
  struct smoothstep *d;
  float pretime, oldin, sine;
  int preread, havestep;

  d= (struct smoothstep *)me->private[0];
  pretime= INCALC(0)/2.0;
  if( !d->countdown )
  {
    if( pretime > 0.5*HORIZON )
      pretime= 0.5*HORIZON;
    else if( pretime < 0.0 )
      pretime= 0.0;
    preread= (int)round(SAMPRATE*pretime);
    BUFINCALC(1, preread-1, 0.0);
    oldin= INPUTC(1, 0);
    BUFINCALC(1, preread, 0.0);
    if( d->flags == 0 )
      havestep= fabsf(INPUTC(1, 0) - oldin) >= d->minstep;
    else {
      havestep= 0;
      if( (d->flags & SMST_UP) != 0 )
	havestep= havestep || INPUTC(1, 0) - oldin >= d->minstep;
      if( (d->flags & SMST_DOWN) != 0 )
	havestep= havestep || INPUTC(1, 0) - oldin <= -d->minstep;
      if( (d->flags & SMST_ABSUP) != 0 )
	havestep= havestep || fabsf(INPUTC(1, 0)) - fabsf(oldin) >= d->minstep;
      if( (d->flags & SMST_ABSDOWN) != 0 )
	havestep= havestep || fabsf(INPUTC(1, 0)) - fabsf(oldin) <= -d->minstep;
    }
    if( havestep )
    {
      d->countdown= 2*preread;
      BUFINCALC(1, 0, 0.0);
      d->start= INPUTC(1, 0);
      BUFINCALC(1, 2*preread, 0.0);
      d->goal= INPUTC(1, 0);
      if( d->mode==SMST_LINEAR ) {
	d->curr= d->start;
	d->inc= (INPUTC(1, 0)-d->curr)/(float)(2*preread);
      }
      else {
	d->curr= -0.25 * SMST_SINESIZE;
	d->inc= (0.5 * SMST_SINESIZE)/(float)(2*preread);
      }
    }
  }
  if( d->countdown )
  {
    if( d->mode==SMST_LINEAR )
      OUTPUT(0, d->curr);
    else {
      sine= 1.0 + lininterpol(d->sinelut, SMST_SINESIZE, d->curr);
      OUTPUT(0, d->start + (d->goal-d->start)*sine/2.0);
    }
    OUTPUT(1, fabsf(me->ch[0]-d->start) < fabsf(me->ch[0]-d->goal) ?
		fabsf(me->ch[0]-d->start) : fabsf(me->ch[0]-d->goal));
    d->curr += d->inc;
    --d->countdown;
  }
  else {
    BUFINCALC(1, 0, 0.0);
    OUTPUT(0, INPUTC(1, 0));
    OUTPUT(1, 0.0);
  }
  BUFINCR(1, 1);
  CALCRETURN;
}



/*=============================================================================
    conditional <input1> <input2> <greaterequal> <lessthan>

Compares the first channel of the first two inputs and outputs the values of
<greaterequal> if the first is larger than or equal to the second, or of
<lessthan> otherwise.  The number of output channels is the maximum of the
channels of <greaterequal> and <lessthan>.

See also: clip, limiter
=============================================================================*/

float calcconditional(sndobj *me, sndobj *, int);

sndobj *conditional( sndobj *in1, sndobj *in2, sndobj *greaterequal, sndobj *lessthan )
{
  return newsndo(calcconditional, "conditional", "conditional", 
	greaterequal->nch> lessthan->nch ? greaterequal->nch: lessthan->nch, 4,
	in1, in2, greaterequal, lessthan );
}

float calcconditional(sndobj *me, sndobj *caller, int innr)
{
  if( INCALC(0) >= INCALC(1) ) {
    INCALC(2);
    INSKIP(3);
    caller->inch[innr]= me->inch[2];
    return me->inch[2][0];
  }
  else {
    INSKIP(2);
    INCALC(3);
    caller->inch[innr]= me->inch[3];
    return me->inch[3][0];
  }
}


/*=============================================================================
    smoothstart (number of samples to modify) <signal>

Smooths the discontinuity which may result when a sound is suddenly switched on
or when a sample whose first value is not zero is played.  <samples> zero
samples followed by a non-zero sample are replaced by a third-order polynomial
which starts at zero with a horizontal derivative and matches the first
non-zero sample and the following one (in effect matching the derivative as
well).  Works only on one channel.

See also: smoothstep
=============================================================================*/

struct smoothstart {
  int	samples, countdown, zeros;
  float	lincoeff, quadcoeff, cubecoeff, nonzeroval;
};

float calcsmoothstart(sndobj *me, sndobj *caller, int innr);

sndobj *smoothstart(int samples, sndobj *signal)
{
  sndobj *p;
  struct smoothstart *d;
  
  p= newsndo(calcsmoothstart, "smoothstart", "smoothstart", 1, 1, signal);
  p->private[0]= d= new(struct smoothstart);
  if( samples<=1 ) {
    fprintf(stderr, "smoothstart: Number of samples (%d) given <=1.  Using default (5).\n", samples);
    samples= 5;
  }
  else if( (double)samples> (double)HORIZON*SAMPRATE ) {
    fprintf(stderr, "smoothstart: Number of samples (%d) given too large.  ", samples);
    samples= (int)round((double)HORIZON*SAMPRATE);
    fprintf(stderr, "Using maximum (%d).\n", samples);
  }
  d->samples= samples;
  d->countdown= -1;
  d->zeros= 0;
  d->nonzeroval= 0.0;
  return p;
}

float calcsmoothstart(sndobj *me, sndobj *caller, int innr)
{
  struct smoothstart *d;
  
  d= (struct smoothstart *)me->private[0];
  if( d->countdown>= 0 )
  {
    double x= d->samples-d->countdown;
    OUTPUT(0, x*(d->lincoeff + x*(d->quadcoeff + d->cubecoeff*x)));
    --d->countdown;
  }
  else if( d->zeros> 0 )
  {
    OUTPUT(0, 0.0);
    if( d->zeros==d->samples )
      if( INCALC(0) )
      {
	double y1= INPUT(0);
	
	INCALC(0);
	d->countdown= d->samples-1;
        d->zeros= 0;
	d->nonzeroval= 0.0;
	d->cubecoeff= (double)INPUT(0)/d->samples/(d->samples+1) - 
			y1/d->samples/(d->samples-1);
	d->quadcoeff= (double)INPUT(0)/d->samples/(d->samples+1) - 
			d->cubecoeff*(d->samples-1);
	d->lincoeff= d->quadcoeff - d->cubecoeff;
      }
      else;
    else
      --d->zeros;
  }
  else if( d->nonzeroval ) {
    OUTPUT(0, d->nonzeroval);
    d->nonzeroval= 0.0;
  }
  else {
    OUTPUT(0, INCALC(0));
    if( !INPUT(0) ) {
      ++d->zeros;
      for( ; d->zeros< d->samples; ++d->zeros )
        if( INCALC(0) ) {
          d->nonzeroval= INPUT(0);
	  break;
        }
    }
  }
  CALCRETURN;
}


/*=============================================================================
    sndthread (fifolen) <subtree>

This sndobj creates a separate thread for computing the data of its input.
Intended for making full use of multi-CPU or multi-core systems.  Choose one or
several sub-trees of snd objects which take an approximately equal share of CPU
time and separate all except one of them from the top level with sndthread.
These subtrees must be true subhierarchies; no data must flow from one of its
sndobj's to the main hierarchy except through one sndthread object.  This
object can't be skipped, since it computes its output concurrently with the
objects above it and cannot know whether its output will be needed.

The first argument gives the number of samples to be calculated in advance by
the thread taking care of the subtree.  The larger this number, the better
variations in computation time per sample can be equalised, but the more
unneeded samples are computed at the end.

This object is unavailable if NO_PTHREADS is defined in sndsys.h.
=============================================================================*/

#ifndef NO_PTHREADS

struct sndthread {
  int		    detach_done;
  pthread_t	    thread;
  pthread_mutex_t   mutex;
  pthread_cond_t    writedone;
  pthread_cond_t    readdone;
  int		    fifolen, watermark, datasize;
  float		    *data, *write;
};

float calcsndthread(sndobj *me, sndobj *caller, int innr);
void exitsndthread(sndobj *me);
void *threadsndthread(void *voidme);

sndobj *sndthread(int fifolen, sndobj *signal)
{
  sndobj *p;
  struct sndthread *d;
  pthread_mutexattr_t mattr;
  int threaderr, count;
  
  p= newsndo(calcsndthread, "thread", "thread", signal->nch, 1, signal);
  p->skip= dontskip;
  p->exit= exitsndthread;
  p->private[0]= d= new(struct sndthread);
  p->private[1]= d->data= news(float, p->nch*fifolen);
  d->write= d->data;
  d->fifolen= fifolen;
  d->watermark= 0;
  d->datasize= fifolen*p->nch;
  d->detach_done= 0;
  threaderr= 0;
  threaderr |= pthread_mutexattr_init(&mattr);
  threaderr |= pthread_mutexattr_settype(&mattr, PTHREAD_MUTEX_ERRORCHECK_NP);
  threaderr |= pthread_mutex_init(&d->mutex, &mattr);
  threaderr |= pthread_mutexattr_destroy(&mattr);
  threaderr |= pthread_cond_init(&d->writedone, NULL);
  threaderr |= pthread_cond_init(&d->readdone, NULL);
  if( threaderr ) {
    fprintf(stderr, "sndthread: Thread data initialisation error. Aborting.\n");
    exit(1);
  }
  return p;
}

float calcsndthread(sndobj *me, sndobj *caller, int innr)
{
  struct sndthread *d;
  float *read;
  int ch;

  d= (struct sndthread *)me->private[0];
  if( !d->detach_done )
  {
    pthread_attr_t attr;
    
    if( pthread_mutex_lock(&d->mutex) ) {
      fprintf(stderr, "calcsndthread for `%s\': Error when first locking mutex. Aborting.\n", me->name);
      exit(1);
    }
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_attr_setinheritsched(&attr, PTHREAD_INHERIT_SCHED);
    if( pthread_create(&d->thread, &attr, threadsndthread, me) ) {
      fprintf(stderr, "calcsndthread for `%s\': Error creating thread. Aborting.\n", me->name);
      exit(1);
    }
    pthread_attr_destroy(&attr);
    d->detach_done= 1;
  }
  pthread_mutex_lock(&d->mutex);
  // Has to be while not if, see below
  while( !d->watermark )
    pthread_cond_wait(&d->writedone, &d->mutex);
  read= d->write - me->nch*d->watermark;
  if( read < d->data )
    read += d->datasize;
  FORCH
    OUTPUT(ch, read[ch]);
  --d->watermark;
  pthread_mutex_unlock(&d->mutex);
  pthread_cond_signal(&d->readdone);
  CALCRETURN;
}

void exitsndthread(sndobj *me)
{
  struct sndthread *d;
  
  d= (struct sndthread *)me->private[0];
  pthread_cancel(d->thread);
  pthread_cond_destroy(&d->readdone);
  pthread_cond_destroy(&d->writedone);
  pthread_mutex_destroy(&d->mutex);
}

void *threadsndthread(void *voidme)
{
  sndobj *me;
  struct sndthread *d;
  int ch;
  
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
  me= (sndobj *)voidme;
  d= (struct sndthread *)me->private[0];
  while( 13 )
  {
    INCALC(0);
    pthread_mutex_lock(&d->mutex);
    // For reasons I don't understand, this while cannot be replaced by an if.
    // Otherwise it happens (with nice -19) that watermark is still ==fifolen
    // even though it must just have been decremented by the reading thread.
    // (I have since learnt that these "spurious wakeups" are explicitly allowed
    // by POSIX.)
    while( d->watermark == d->fifolen )
      pthread_cond_wait(&d->readdone, &d->mutex);
    FORCH
      d->write[ch]= INPUTC(0, ch);
    d->write += me->nch;
    if( d->write >= d->data+d->datasize )
      d->write -= d->datasize;
    ++d->watermark;
    pthread_mutex_unlock(&d->mutex);
    pthread_cond_broadcast(&d->writedone);
  }
}

#endif


/*=============================================================================
    skipwatch "filename" <signal>

A debugging object which allows to watch which values of a signal would
normally be skipped.  This sound object uses a writeout object and a further
auxiliary object to write all values to the file, transformed in the following
way: Values for which calc was called are linearly mapped to the interval [0,
1] (from [-1, 1]); values which would have been skipped are mapped to [-1, 0].

See also: beancounter
=============================================================================*/

struct skipwatchinternal {
  int isaskip;
};

sndobj *skipwatchinternal(sndobj *signal);

float calcskipwatch(sndobj *me, sndobj *caller, int innr);
void skipskipwatch(sndobj *me, sndobj *caller, int innr);

sndobj *skipwatch(const char *filename, sndobj *signal)
{
  sndobj *p, *pp;

  pp= skipwatchinternal(signal);
  p= newsndo(calcskipwatch, "skipwatch", "skipwatch", signal->nch, 1, 
      						writeout(filename, pp) );
  p->skip= skipskipwatch;
  p->private[0]= new(struct skipwatchinternal *);
  *(struct skipwatchinternal **)p->private[0]= 
    		(struct skipwatchinternal *)pp->private[0];
  return p;
}

float calcskipwatch(sndobj *me, sndobj *caller, int innr)
{
  int ch;

  (*(struct skipwatchinternal **)me->private[0])->isaskip= 0;
  INCALC(0);
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  CALCRETURN;
}

void skipskipwatch(sndobj *me, sndobj *caller, int innr)
{
  int ch;

  (*(struct skipwatchinternal **)me->private[0])->isaskip= 1;
  INCALC(0);
}

float calcskipwatchinternal(sndobj *me, sndobj *caller, int innr);

sndobj *skipwatchinternal(sndobj *signal)
{
  sndobj *p;

  p= newsndo(calcskipwatchinternal, "skipwatchinternal", "skipwatchinternal", signal->nch, 1, signal);
  p->skip= dontskip;
  p->private[0]= new(struct skipwatchinternal);
  ((struct skipwatchinternal *)p->private[0])->isaskip= 0;
  return p;
}

float calcskipwatchinternal(sndobj *me, sndobj *caller, int innr)
{
  int ch;

  INCALC(0);
  if( ! ((struct skipwatchinternal *)me->private[0])->isaskip )
    FORCH
      if( INPUTC(0, ch) > 1.0 )
	OUTPUT(ch, 1.0);
      else if( INPUTC(0, ch) < -1.0 )
	OUTPUT(ch, 0.0);
      else
	OUTPUT(ch, 0.5 + 0.5*INPUTC(0, ch));
  else
    FORCH
      if( INPUTC(0, ch) > 1.0 )
	OUTPUT(ch, 0.0);
      else if( INPUTC(0, ch) < -1.0 )
	OUTPUT(ch, -1.0);
      else
	OUTPUT(ch, -0.5 + 0.5*INPUTC(0, ch));
  CALCRETURN;
}



/*=============================================================================
    panpot (method) <right> <signal>

Creates stereo signals from the channels of its last input <signal>.
<right> gives the signal's perceived position: Only on right channel for 1,
only on left channel for -1, centred for 0.  For other values, the secondary
(softer) channel is delayed with respect to the primary channel depending on
(method).  If (method) is 0, the perceived direction is created solely by
varying the intensity of the sound (intensity stereophony), if 1, only a delay
is applied (time-of-arrival stereophony).  Note that pure time-of-arrival
stereophony is no good for extreme left or right panning.  Besides, time-of-
arrival stereophony can lead to comb filter effects when the stereo signal is
later converted to mono, and to beat effects when the position <right>
changes.

<right> should change slowly.  As panpot doubles the number of channels, 
<signal> can have at most PLENTY/2 channels.  The even output channels (0, 2,
...) contain the left signals, the odd channels the right signals.

See also: pansur, mix2
=============================================================================*/

struct panpot {
  float method;
};

float calcpanpot(sndobj *me, sndobj *caller, int innr);
void skippanpot(sndobj *me, sndobj *caller, int innr);

sndobj *panpot( float method, sndobj *right, sndobj *signal )
{
  sndobj *p;
  struct panpot *d;
  int channels;

  if( signal->nch > PLENTY/2 ) {
    fprintf(stderr, "panpot: Warning: input signal has too many channels.  Ignoring surplus channels.\n");
    channels= PLENTY;
  }
  else
    channels= signal->nch*2;
  p= newsndo(calcpanpot, "panpot", "panpot", channels, 2, right, signal);
  p->skip= skippanpot;
  p->private[0]= d= new(struct panpot);
  if( method< 0 )
    d->method= 0;
  else if( method> 1 )
    d->method= 1;
  else
    d->method= method;
  buf(1.0, 0, 0, 1, p, 1);
  return p;
}

float calcpanpot(sndobj *me, sndobj *caller, int innr)
{
  struct panpot *d;
  float rightin, righttoa, fdelay, leftampl, rightampl;
  int idelay, ch;

  d= (struct panpot *)me->private[0];
  rightin= INCALC(0);
  rightampl= sqrtf(0.5f + 0.5f * rightin * (1.0f - d->method));
  leftampl= sqrtf(0.5f - 0.5f * rightin * (1.0f - d->method));
  righttoa= 0.5f + 0.5f * rightin * d->method;
  fdelay= 0.00073f * 10 * fabsf(logf(righttoa/(1.0f-righttoa))) / (float)M_LN10;
  if( !finite(fdelay) || fdelay > 0.0015f ) {
    idelay= 0.0015*SAMPRATE;
    fdelay= 0;
  }
  else {
    idelay= ceilf( fdelay*SAMPRATE );
    fdelay= SAMPRATE*fdelay - idelay;
  }
  BUFINCALC(1, 0, 0.0);
  if( rightin >= 0 )
    for( ch= 1; ch< me->nch; ch += 2 )
      OUTPUT(ch, INPUTC(1, ch/2) * rightampl);
  else
    for( ch= 0; ch< me->nch; ch += 2 )
      OUTPUT(ch, INPUTC(1, ch/2) * leftampl);
  BUFINCALC(1, -idelay, -fdelay);
  if( rightin >= 0 )
    for( ch= 0; ch< me->nch; ch += 2 )
      OUTPUT(ch, INPUTC(1, ch/2) * leftampl);
  else
    for( ch= 1; ch< me->nch; ch += 2 )
      OUTPUT(ch, INPUTC(1, ch/2) * rightampl);
  BUFINCR(1, 1);
  CALCRETURN;
}

void skippanpot(sndobj *me, sndobj *caller, int innr)
{
  INSKIP(0);
  BUFINCR(1, 1);
}


/*=============================================================================
    pansur [speakerpos] <sigpos> <sigrad> <signal>

Surround panning object.  The first channel of the input <signal> is panned to
a set of surround speakers, which correspond to output channels.  [speakerpos]
gives the angles of the speakers with respect to the listening position,
<sigpos> gives the angular position of the signal.  Angles are given in units
of 2pi, so for instance 0.25 is 90 degrees.  This object is indifferent about
whether you define these angles as clockwise or anticlockwise; you merely have
to be consistent.  The [speakerpos] array has to be concluded by the magic
value END after the last speaker angle.  Alternatively, the macro
POS_REGULAR(n) can be passed for [speakerpos], where 1<=n<=128.  Then pansur
will compute the speaker positions at the corners of an n-sided polygon,
ordered with increasing angles.

The <sigrad> input determines how localised the sound is.  It can be
interpreted as the radius of the signal w.r.to the listening position.  If
<sigrad> is 0, the signal will be split up equally among the speakers; if it is
1, only the speaker(s) which is closest to <sigpos> will receive any
contribution.

The panning is done so that the total intensity (squared amplitude) of the
signal is preserved.  For <sigrad> < 1, the output intensities are proportional
to (cos(speakerpos[i]-<sigpos>))**(<sigrad>/(1-<sigrad>)).  Multiple speakers
having the same position are handled gracefully.

See also: panpot, mix
=============================================================================*/

struct pansur {
  float *speakerpos;
};

float calcpansur(sndobj *me, sndobj *caller, int innr);

sndobj *pansur(float *speakerpos, sndobj *sigpos, sndobj *sigrad, sndobj *signal)
{
  sndobj *p;
  struct pansur *d;
  int nch;

  p= newsndo(calcpansur, signal->name, "pansur", 1, 3, sigpos, sigrad, signal);
  p->private[0]= d= new(struct pansur);
  fprintf(stderr, "speakerpos= %p\n", speakerpos);
  if( ((unsigned long)speakerpos & ~(POS_REG_MASK<<POS_REG_SHIFT)) == POS_REG_COMM ) {
    int ind;

    fprintf(stderr, "found regular pos pseudo-pointer\n");
    p->nch= ((unsigned long)speakerpos >> POS_REG_SHIFT) & POS_REG_MASK;
    if( !p->nch )
      p->nch= POS_REG_MASK+1;
    if( p->nch > PLENTY ) {
      fprintf(stderr, "pansur: too many output channels (%d; maximum is %d).  Ignoring excess.\n", p->nch, PLENTY);
      p->nch= PLENTY;
    }
    p->private[PLENTY-1]= d->speakerpos= news(float, p->nch+1);
    for( ind= 0; ind <= p->nch/2; ++ind )
      d->speakerpos[ind]= (float)ind/(float)p->nch;
    for( ; ind < p->nch; ++ind )
      d->speakerpos[ind]= (float)ind/(float)p->nch - 1.0f;
    d->speakerpos[ind]= END;	// don't really need this; just in case
  }
  else {
    float *scanpos;

    fprintf(stderr, "found NO regular pos pseudo-pointer\n");
    for( scanpos= speakerpos; *scanpos!=END; ++scanpos );
    p->nch= scanpos-speakerpos;
    if( p->nch > PLENTY ) {
      fprintf(stderr, "pansur: too many output channels (%d; maximum is %d).  Ignoring excess.\n", p->nch, PLENTY);
      p->nch= PLENTY;
      speakerpos[PLENTY]= END;
    }
    d->speakerpos= speakerpos;
  }
  return p;
}

float calcpansur(sndobj *me, sndobj *caller, int innr)
{
  struct pansur *d;
  double sp_levels[PLENTY], pos, rad, exponent;
  int ch;

  d= (struct pansur *)me->private[0];
  pos= INCALC(0);
  rad= INCALC(1);
  if( rad<1 ) {
    double norm;

    if( rad <= 0 )	exponent= 0;
    else		exponent= rad/(1.0-rad);
    norm= 0;
    FORCH {
      sp_levels[ch]= pow(0.5+0.5*cos(2.0*M_PI*(pos-d->speakerpos[ch])), exponent);
      norm += sp_levels[ch];
    }
    FORCH
      sp_levels[ch]= sqrt(sp_levels[ch]/norm);
  }
  else {	// search for closest speaker for exponent = infinity
    double thisdist, mindist, level;
    int closest[PLENTY], closeind;

    mindist= 1;
    closeind= 0;
    FORCH {
      thisdist= fabs(fmod(d->speakerpos[ch]-pos, 1.0));
      if( thisdist > 0.5 )
	thisdist= 1.0-thisdist;
      if( thisdist <= mindist ) {
	if( thisdist < mindist )
	  closeind= 0;
	closest[closeind++]= ch;
	mindist= thisdist;
      }
    }
    level= sqrt(1.0/(double)closeind);
    FORCH
      sp_levels[ch]= 0;
    for( --closeind; closeind>= 0; --closeind )
      sp_levels[closest[closeind]]= level;
  }
  INCALC(2);
  FORCH
    OUTPUT(ch, sp_levels[ch]*INPUT(2));
  CALCRETURN;
}

 
/*=============================================================================
    invA <freq>

returns the square root of the inverse of the A-weighting factor for a specific
frequency.  The intended use is to compensate for frequency-dependent human
loudness perception by multiplying with this factor.  Due to the square root
inherent in the calculation, this sndobj's output has to be multiplied with the
waveform (amplitude), even though the A-weighting is defined with respect to
the intensity (squared amplitude).
=============================================================================*/

float calcinvA(sndobj *me, sndobj *caller, int innr);

sndobj *invA(sndobj *freq)
{
  return newsndo(calcinvA, freq->name, "invA", freq->nch, 1, freq);
}

float calcinvA(sndobj *me, sndobj *caller, int innr)
{
  double freq, freq2, freq8, tmp1, tmp2, invafact;
  int ch= 0;

  freq= (double)INCALC(0);
  while( 13 )
  {
    freq2= freq*freq;
    freq8= freq2*freq2;
    freq8= freq8*freq8;
    tmp1= freq2+20.598997*20.598997;
    tmp2= freq2+12194.22*12194.22;
    invafact= sqrt( (freq2+107.65265*107.65265)* (freq2+737.86223*737.86223) * 
        		tmp1*tmp1 * tmp2*tmp2 / (1.562339*2.242881e16*freq8) );
    OUTPUT(ch, (float)invafact);
    if( ++ch >= me->nch )
      break;
    freq= INPUTC(0, ch);
  }
  CALCRETURN;
}


/*=============================================================================
    sndassert "cmp" "msg" <signal>

This object aborts the program if <signal> does not compare to zero as
specified by "cmp".  The message "msg" is then output before aborting.  "cmp"
can be one of "<", ">", "=", "<=", ">=", "==" and "!=".  The input <signal> is
compared to constant zero, so "<" causes the program to terminate if the signal
contains a value which is >=0.  "=" and "==" are interchangeable.  Be aware of
numerical issues when comparing floating-point values (as samples are within
sndsys) for equality:  Two signals which are equal theoretically need not yield
numerical zero when subtracted from each other.  Sine waves will compare !=0
unless sampled exactly at their zero, which tends not to happen except for low
frequencies or for frequencies dividing the sampling rate.
=============================================================================*/

struct sndassert {
  int	gt, neq;
  const char *msg;
  int	samplecount;
};

float calcsndassert(sndobj *me, sndobj *caller, int innr);

sndobj *sndassert(const char *cmp, const char *msg, sndobj *signal)
{
  sndobj *p;
  struct sndassert *d;
  int gt, neq;

  if( !*cmp || (cmp[1] && cmp[1]!='=') ) {
    fprintf(stderr, "Warning: sndassert: Illegal assert comparison `%s'.  Returning signal.\n", cmp);
    return signal;
  }
  neq= !cmp[1];
  if( *cmp=='=' ) {
    gt= 0;
    neq= 0;
  }
  else if( *cmp=='<' )
    gt= -1;
  else if( *cmp=='>' )
    gt= 1;
  else if( *cmp=='!' && cmp[1] ) {
    gt= 0;
    neq= 1;
  }
  else {
    fprintf(stderr, "Warning: sndassert: Illegal assert comparison `%s'.  Returning signal.\n", cmp);
    return signal;
  }
  p= newsndo(calcsndassert, signal->name, "sndassert", signal->nch, 1, signal);
  p->skip= dontskip;
  p->private[0]= d= new(struct sndassert);
  d->gt= gt;
  d->neq= neq;
  d->msg= msg;
  d->samplecount= 0;
  return p;
}

float calcsndassert(sndobj *me, sndobj *caller, int innr)
{
  struct sndassert *d;
  int ch, ok;

  d= (struct sndassert *)me->private[0];
  INCALC(0);
  ok= 1;
  if( d->gt> 0 )
    if( d->neq )
      FORCH
	ok= ok && (INPUTC(0, ch) > 0.0f);
    else
      FORCH
	ok= ok && (INPUTC(0, ch) >= 0.0f);
  else if( d->gt< 0 )
    if( d->neq )
      FORCH
	ok= ok && (INPUTC(0, ch) < 0.0f);
    else
      FORCH
	ok= ok && (INPUTC(0, ch) <= 0.0f);
  else
    if( d->neq )
      FORCH
	ok= ok && (INPUTC(0, ch) != 0.0f);
    else
      FORCH
	ok= ok && (INPUTC(0, ch) == 0.0f);
  if( !ok ) {
    fprintf(stderr, "sndassert `%s' on sndobj `%s', type %s failed\n"
	"at sample %d, time %g", me->name, me->in[0]->name, me->in[0]->type,
	d->samplecount, (double)d->samplecount/(double)SAMPRATE+starttime);
    if( starttime )
      fprintf(stderr, " (start time %g)", starttime);
    fprintf(stderr, ".\n%s\nAborting.\n", d->msg);
    exit(1);
  }
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  ++d->samplecount;
  CALCRETURN;
}





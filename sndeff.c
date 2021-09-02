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
    freeverb <roomsize> <damping> <stereomix> <revmix> <signal>

Freeverb Version 3, by Jezar.  Very sopisticated reverberator simulating sound
within a room.  Its two components are eight parallel comb filters with lowpass
filters in their feedback line, which simulate low-frequency reverberation and
the damping of high frequencies on the walls of the room; and four sequential
not-quite-allpass filters which provide high-frequency reverberation.

This object always outputs a stereo signal; its signal input may be either mono
or stereo.  The output usually has to be scaled since the effective gain
depends on the nature of the signal and on the parameters in a non-trivial way.

All the parameter inputs (all except the signal) should be in the range [0,1].
<roomsize> stands for the size of the simulated room.  This is somewhat
misnamed, as it controls not the length of the delay lines and the
reverberation, but the feedback coefficient of the allpasses.  A more accurate
interpretation is as the low-frequency reflectivity of the walls.  <damping>
describes the damping by the walls of the room, especially of high frequencies.
<stereomix> determines how much of the opposite channel's reverberation is
mixed into each channel - 0 means only contributions from the same side, 1 only
from the opposite side.  <revmix> gives the mixing of the reverberation with
the original signal: 0 only original, 1 only reverberation.

See also: nrev, josrev, fdnrev, boxrev, convrev; _lbcf, _fvap
=============================================================================*/

#define FREEVERB_COMBS	    8
#define FREEVERB_ALLPASSES  4

#define FREEVERB_SCALEDAMP  0.4f
#define FREEVERB_SCALEROOM  0.28f
#define FREEVERB_OFFSETROOM 0.7f
#define FREEVERB_GAIN	    0.1f

struct freeverb {
  float	rmsize, stereo, damp, revmix;
  float stereostraight, stereocross, orig;
  sndobj *combs[2*FREEVERB_COMBS], *allpasses[2*FREEVERB_ALLPASSES];
};
struct _lbcf {
  float	feedback, damp1, damp2, filterstore;
  int	bufsize, index;
  float *buf;
};

float calcfreeverb(sndobj *me, sndobj *, int);

sndobj *freeverb( sndobj *roomsize, sndobj *damping, sndobj *stereomix,
		    sndobj *revmix, sndobj *signal )
{
  const int stereo_spread = 23;
  static const int comb_tuning[FREEVERB_COMBS] = {1116, 1188, 1277, 1356, 1422, 1491, 1557, 1617};
  static const int allpass_tuning[FREEVERB_ALLPASSES] = {556, 441, 341, 225};

  sndobj *p, *left, *right, *rightsig;
  struct freeverb *d;
  float comb_feedback, comb_lowcoeff;
  int i;
  
  p= newsndo(calcfreeverb, "freeverb", "freeverb", 2, 5,
	roomsize, damping, stereomix, revmix, signal /*further inputs below*/);
  p->private[0]= d= new(struct freeverb);
  d->rmsize= 0.5;
  d->stereo= 1.;
  d->damp= 0.5;
  d->revmix= 0.1;
  d->orig= FREEVERB_GAIN*(1.0 - d->revmix);
  d->stereostraight= FREEVERB_GAIN * d->revmix * (0.5*d->stereo + 0.5);
  d->stereocross= FREEVERB_GAIN * d->revmix * (-0.5*d->stereo + 0.5);
  left= add(0);	    // implies nch=1
  right= add(0);
  rightsig= signal->nch > 1 ? c1(signal) : signal;
  comb_feedback= d->rmsize*FREEVERB_SCALEROOM + FREEVERB_OFFSETROOM;
  comb_lowcoeff= d->damp*FREEVERB_SCALEDAMP;
  for( i= 0; i< FREEVERB_COMBS; ++i )  {
    d->combs[i]= _lbcf(comb_tuning[i], comb_feedback, comb_lowcoeff, signal);
    addinput(left, d->combs[i]);
    d->combs[FREEVERB_COMBS+i]= 
  _lbcf(comb_tuning[i] + stereo_spread, comb_feedback, comb_lowcoeff, rightsig);
    addinput(right, d->combs[FREEVERB_COMBS+i]);
  }
  for( i= 0; i< FREEVERB_ALLPASSES; ++i ) {
    left= _fvap(allpass_tuning[i], 0.5, left);
    right= _fvap(allpass_tuning[i] + stereo_spread, 0.5, right);
  }
  addinput(p, left);
  addinput(p, right);
  return p;
}

float calcfreeverb(sndobj *me, sndobj *caller, int innr)
{
  struct freeverb *d;
  float newroom, newdamp, newstereo, newmix;
  int i;
  
  d= (struct freeverb *)me->private[0];
  newroom= INCALC(0);
   if( newroom != d->rmsize )
   {
     float feedback;
     
     feedback= FREEVERB_SCALEROOM * newroom + FREEVERB_OFFSETROOM;
     for( i= 0; i< 2*FREEVERB_COMBS; ++i )
       ((struct _lbcf *)d->combs[i]->private[0])->feedback= feedback;
     d->rmsize= newroom;
  }
  newdamp= INCALC(1);
  if( newdamp != d->damp )
  {
    float damp1, damp2;
    
    damp1= FREEVERB_SCALEDAMP * newdamp;
    damp2= 1.0f - damp1;
    for( i= 0; i< 2*FREEVERB_COMBS; ++i ) {
      ((struct _lbcf *)d->combs[i]->private[0])->damp1= damp1;
      ((struct _lbcf *)d->combs[i]->private[0])->damp2= damp2;
    }
    d->damp= newdamp;
  }
  newstereo= INCALC(2);
  newmix= INCALC(3);
  if( newstereo != d->stereo || newmix != d->revmix ) {
    d->stereo= newstereo;
    d->revmix= newmix;
    d->orig= FREEVERB_GAIN*(1.0 - d->revmix);
    d->stereostraight= FREEVERB_GAIN * d->revmix * (0.5*d->stereo + 0.5);
    d->stereocross= FREEVERB_GAIN * d->revmix * (-0.5*d->stereo + 0.5);
  }
  INCALC(4);
  INCALC(5);
  INCALC(6);
  if( me->in[4]->nch > 1 ) {
    OUTPUT(0, d->stereostraight*INPUT(5) + d->stereocross * INPUT(6)
		+ d->orig * INPUTC(4, 0) );
    OUTPUT(1, d->stereostraight*INPUT(6) + d->stereocross * INPUT(5)
		+ d->orig * INPUTC(4, 1) );
  }
  else {
    OUTPUT(0, d->stereostraight*INPUT(5) + d->stereocross * INPUT(6)
		+ d->orig * INPUT(4) );
    OUTPUT(1, d->stereostraight*INPUT(6) + d->stereocross * INPUT(5)
		+ d->orig * INPUT(4) );
  }
  CALCRETURN;
}


/*=============================================================================
    _lbcf (delay) (feedback) (lowcoeff) <signal>

Lowpass backfeed comb filter used internally by freeverb.  (feedback) is the
comb filter backfeed coefficient, (delay) the length of its delay line in
samples.  The backfeed loop of the comb filter contains a lowpass filter with
the difference equation y(n) = (1-(lowcoeff)) x(n) + (lowcoeff) y(n-1).  A
higher (lowcoeff) causes a lower cutoff frequency.

The cutoff frequency of the feedback lowpass is:
SAMPRATE/2/pi*arccos((x**2+1)/2/x - sqrt(((x**2+1)/2/x)**2 -1))
where x = (lowcoeff).

This sndobj's name starts with an underscore to indicate it is primarily for
internal use.  It does not have inputs, but its private data can be adjusted
between calls of calc().  This sndobj is one of the two components used by the
freeverb reverberator.  It works only on one channel.

See also: freeverb, _fvap, _comb
=============================================================================*/

// see freeverb for struct _lbcf

float calc_lbcf(sndobj *me, sndobj *, int);

sndobj *_lbcf( int delay, float feedback, float lowcoeff, sndobj *signal )
{
  sndobj *p;
  struct _lbcf *d;
  int i;
  
  p= newsndo(calc_lbcf, "_lbcf", "_lbcf", 1, 1, signal);
//  p->skip= dontskip;
  p->private[0]= d= new(struct _lbcf);
  p->private[1]= d->buf= news(float, delay);
  d->bufsize= delay;
  for( i= 0; i< delay; ++i )
    d->buf[i]= 0.0;
  d->index= 0;
  d->filterstore= 0.0;
  d->feedback= feedback;
  d->damp1= lowcoeff;
  d->damp2= 1.0f - d->damp1;
  return p;
}

float calc_lbcf(sndobj *me, sndobj *caller, int innr)
{
  struct _lbcf *d;
  
  d= (struct _lbcf *)me->private[0];
  OUTPUT(0, d->buf[d->index]);
  d->filterstore= d->damp2 * me->ch[0] + d->damp1 * d->filterstore;
  d->buf[d->index]= INCALC(0) + d->feedback * d->filterstore;
  if( ++d->index >= d->bufsize )
    d->index= 0;
  CALCRETURN;
}


/*=============================================================================
    _fvap (delay) (feedback) <signal>

Not-quite-allpass filter used internally by freeverb.  Its amplitude response
is that of a comb filter turned upside down, with narrow spikes extending to
zero.  Its phase response resembes that of an allpass but has additional jumps
between the wraparound points from -pi to pi.  Its difference equation is:
y(n)= -x(n) + x(n-(delay)) + (feedback) y(n-(delay))

This sndobj's name starts with an underscore to indicate it is primarily for
internal use.  It does not have inputs, but its private data can be adjusted
between calls of calc().  This sndobj is one of the two components used by the
freeverb reverberator.  It works only on one channel.

See also: freeverb, _lbcf, _allpass
=============================================================================*/

struct _fvap {
  float	feedback;
  int	bufsize, index;
  float *buf;
};

float calc_fvap(sndobj *me, sndobj *, int);

sndobj *_fvap( int bufsize, float feedback, sndobj *signal )
{
  sndobj *p;
  struct _fvap *d;
  int i;
  
  p= newsndo(calc_fvap, "_fvap", "_fvap", 1,
		1, signal);
  p->private[0]= d= new(struct _fvap);
  p->private[1]= d->buf= news(float, bufsize);
  for( i= 0; i< bufsize; ++i )
    d->buf[i]= 0.0;
  d->bufsize= bufsize;
  d->index= 0;
  d->feedback= feedback;
  return p;
}

float calc_fvap(sndobj *me, sndobj *caller, int innr)
{
  struct _fvap *d;
  
  d= (struct _fvap *)me->private[0];
  OUTPUT(0, - INCALC(0) + d->buf[d->index]);
  d->buf[d->index]= d->feedback * d->buf[d->index] + (1.0f-d->feedback)*INPUT(0);
  if( ++d->index >= d->bufsize )
    d->index= 0;
  CALCRETURN;
}


/*=============================================================================
    nrev (length) (lowcoeff) (feedback) <signal>

Mike McNabb's nrev reverberator.  This one has to be scaled down (~0.1) and
added to the original signal.  (length) determines the reverberation time (~ 1
..  10).  (feedback) (~ 0..1) determines how much the reverberation is
superimposed over itself.  (lowcoeff) sets the lowpass coefficient (0..1, NOT a
frequency).

This reverberation tends to sound a bit like an echo.  It simulates large
rooms with bare walls.

See also: freeverb, josrev, fdnrev, boxrev, convrev
=============================================================================*/

sndobj *nrev( float length, float lowpass, float feedback, sndobj *signal ) 
{
  /* Mike McNabb's nrev from Mus10 days (ca. 1978) */
  #define NREV_BASE_DLY_LEN 14
  static int base_dly_len[NREV_BASE_DLY_LEN] = {1433, 1601, 1867, 2053, 2251, 2399, 347, 113, 37, 59, 43, 37, 29, 19};
  #define NREV_COMBS 6      /*  ALSO CHANGE BELOW!  */
  static float comb_factors[NREV_COMBS] = {0.822, 0.802, 0.773, 0.753, 0.753, 0.733};
  #define NREV_ALLPASSES    4
  
  sndobj *rev, *combs[NREV_COMBS];
  int dly_len[NREV_BASE_DLY_LEN];
  int i, delayind, len;

  length *= (float)SAMPRATE / 25641.0;
  for (i = 0; i < NREV_BASE_DLY_LEN; ++i )
    dly_len[i] = getprime((int)(length * base_dly_len[i]));
  delayind= 0;
  for( i = 0; i < NREV_COMBS; i++ )
    combs[i] = _comb( dly_len[delayind++], comb_factors[i] * feedback, signal);
  // CHANGE THIS TOO IF NREV_COMBS IS CHANGED!
  rev= add( NREV_COMBS, combs[0], combs[1], combs[2], 
			combs[3], combs[4], combs[5] );
  rev->name= "nrev-combs";
  for( i = 0; i < NREV_ALLPASSES-1; ++i )
    rev= _allpass( dly_len[delayind++], 0.7, rev );
  rev= _12_0pole( 0, lowpass, lowpass - 1.0, 0.0, rev );
//  for( i = 0; i < channels; i++ )  normally for all channels separately...
  if( delayind < NREV_BASE_DLY_LEN )
    len = dly_len[delayind++]; 
  else
    len = getprime((int)(length * (20.0 + 40.0*random()/(float)RAND_MAX)));
  rev= _allpass( len, 0.7, rev );
  return rev;
}


/*=============================================================================
    josrev <lowpassf> <in>

Good general-purpose reverberator based on suggestions in J. O. Smith's online
books.  Has to be scaled and added to the original signal.

See also: freeverb, nrev, fdnrev, boxrev, convrev
=============================================================================*/

float josrevforwcombs[]= { 0, 1, 4799, 0.742/4, 4999, 0.733/4, 5399, 0.715, 5801, 0.697, END };

sndobj *josrev( sndobj *lowpassf, sndobj *in )
{
  return filter( josrevforwcombs, NULL, lowpass( lowpassf,
		  _allpass( 113, 0.7, _allpass( 337, 0.7, 
		  _allpass( 1051, 0.7, _allpass( 3781, 0.7, 
		  _allpass( 10013, 0.7, in)))))));
}


/*=============================================================================
    fdnrev [delays] (t60low) (t60high) <signal>

Single-channel feedback delay network reverberator with a Householder
reflection matrix for backfeed channel mixing.  [delays] contains the delay
lengths in samples, followed by a 0.  They have to be 5 or larger and will be
substituted by the next larger prime to obtain mutually prime delay line
lengths.  Useful values are around 300 and larger; for smaller delays this
reverberator is susceptible to ringing.  (t60low) and (t60high) are
reverberation times for low/high frequencies in seconds.

Due to the low passes in the backfeed loop, this reverberator tends to make the sound smoother but lower.

See also: freeverb, nrev, josrev, boxrev, convrev
=============================================================================*/

float calchhreflect(sndobj *me, sndobj *, int);

#define FDNREV_MAX      20
sndobj *fdnrev( int *delays, double t60low, double t60high, sndobj *signal )
{
  sndobj *hhr, **del, **bf, *result;
  int *readdel;
  double deltime, g, a1, normcutoff;
  int ndel, i, dellength;

  if( signal->nch > 1 ) {
    fprintf(stderr, "Warning: fdnrev works only on single-channel inputs.  Processing first channel only.\n");
    signal= c0(signal);
  }
  if( t60low < 1.0/SAMPRATE || t60high < 1.0/SAMPRATE ) {
    fprintf(stderr, "Warning: fdnrev: Reverberation time too low (minimum %g). Returning signal.\n", 1.0/SAMPRATE);
    return signal;
  }
  for( readdel= delays; *readdel; ++readdel );
  ndel= readdel-delays;
  if( !ndel ) {
    fprintf(stderr, "Warning: fdnrev: No delays given. Returning signal.\n");
    return signal;
  }
  if( ndel > FDNREV_MAX ) {
    fprintf(stderr, "Warning: fdnrev: Number of delays (%d) exceeds sane "
		"maximum (%d). Setting to max.\n", ndel, (int)FDNREV_MAX );
    ndel= FDNREV_MAX;
  }
  del= news(sndobj*, ndel);
  bf= news(sndobj*, ndel);
  hhr= newsndo( calchhreflect, "fdnrev", "hhreflect", ndel, 0 );
  signal= ch(0, signal);   // single-channel only
  readdel= delays;
  for( i= 0; i< ndel; ++i ) {
    bf[i]= add(1, signal );
    if( *readdel<5 ) {
      fprintf(stderr, "Warning: fdnrev: cannot have delays < 5 samples. Setting to 5.\n");
      *readdel= 5;
    }
    dellength= getprime(*readdel++);
    deltime= (double)dellength/(double)SAMPRATE;
    g= pow( 10.0, -3.0*deltime/(t60low + 2.0/(double)dellength*(t60high-t60low)) );
    normcutoff= 2.0*M_PI/(double)dellength;
    a1= (2.0 + sin(normcutoff)*sin(normcutoff) - sin(normcutoff)
	*sqrt(8.0 + sin(normcutoff)*sin(normcutoff))) / (2.0*cos(normcutoff));
    del[i]= sdelay(dellength, bf[i]);
    addinput( hhr, _12_0pole(0, g*(1.0-a1), -a1, 0.0, del[i]) );
  }
  for( i= 0; i< ndel; ++i )
    addinput( bf[i], ch(i, hhr) );
  result= add( 1, del[0] );
  for( i= 1; i< ndel; ++i )
    addinput( result, del[i] );
  free(del);
  free(bf);
  return result;
}

// Householder reflection for single-channel inputs
float calchhreflect(sndobj *me, sndobj *caller, int innr )
{
  float sum;
  int i;
  
  sum= 0.0;
  for( i= 0; i< me->nin; ++i ) {
    OUTPUT(i, INCALC(i));
    sum += INPUT(0);
  }
  sum *= 2.0/(float)me->nin;
  for( i= 0; i< me->nin; ++i )
    me->ch[i] -= sum;
  CALCRETURN;
}


/*=============================================================================
    boxrev [del_refl] (cutoff) (random) <signal>

Reverberation simulating the ways sound can travel in a box-shaped room.  The
sound source is assumed to be located in the middle of one wall of the room,
and the listener at the same position opposite.  The array [del_refl] contains
pairs of delays and reflection coefficients, followed by the value END.  The
delays give the time sound takes to traverse a certain dimension of the room
(3.3 ms per metre; values around 0.1 sec give decent results); the reflection
coefficients determine the attenuation as a result of reflection on one of the
walls limiting the room in this direction.  There can be an arbitrary number of
delay/reflection coefficient pairs, so this is in fact the simulation of a
"hyper"box.  This object is implemented as a FIR filter.  <cutoff> determines
when the series of filter coefficients should be truncated.  If > 0, it gives a
maximal delay in seconds (corresponding to the order of the FIR filter); if <
0, the negated minimum amplitude (filter coefficient) to be retained.  <random>
is the RMS ratio of a random contribution added to the delays.

See also: freeverb, nrev, josrev, fdnrev, convrev
=============================================================================*/

#define BOXREV_MAXREFL	1.0
#define BOXREV_MINMAXC	256
#define BOXREV_MAXMAXC	(5*SAMPRATE)

int boxrev_addcoeff( float **coeffs, int *maxcoeffs, float delay, float coefficient );

sndobj *boxrev(float *del_refl, float cutoff, float random, sndobj *signal)
{
  sndobj *p;
  rngdesc *rng;
  float delay[PLENTY], refl[PLENTY], mode[PLENTY], depthdelay[PLENTY];
  float dlysqr_[PLENTY+1], atten_[PLENTY+1];
  float *dlysqr= dlysqr_+1, *atten= atten_+1;
  float *coeffs;
  float basedelay, thisdelay;
  int maxcoeffs, count, depth, maxdepth;

  if( *del_refl == END ) {
    fprintf(stderr, "boxrev: Warning: No delays / reflection coefficients given.  Returning unchanged signal.\n");
    return signal;
  }
  if( cutoff == 0.0f ) {
    fprintf(stderr, "boxrev: Warning: Cutoff = 0.  Returning unchanged signal.\n");
    return signal;
  }
  for( count= 0; count< PLENTY && del_refl[2*count]!=END; ++count ) {
    delay[count]= fabsf(del_refl[2*count]);
    if( cutoff < 0 && fabsf(del_refl[2*count+1]) > BOXREV_MAXREFL ) {
      fprintf(stderr, "boxrev: Warning: reflection index %d (counting from 0) is larger than the allowed maximum (%g) for an amplitude cutoff.  Replaced by maximum.\n", count, (double)BOXREV_MAXREFL);
      refl[count]= del_refl[2*count+1] > 0 ? BOXREV_MAXREFL : -BOXREV_MAXREFL;
    }
    else
      refl[count]= del_refl[2*count+1];
  }
  if( del_refl[2*count] != END ) {
    fprintf(stderr, "boxrev: Warning: Can use a maximum dimension of %d.  Some delay-value pairs are ignored.\n", PLENTY);
    maxdepth= PLENTY-1;
  }
  else
    maxdepth= count-1;
  if( cutoff < 0 )
    cutoff /= delay[0];	// for easier comparison below
  else
    cutoff += delay[0];	// for easier comparison below
  basedelay= delay[0];	// sound has to travel across the room before first reaching the listener
  delay[0] *= 2;
  refl[0] *= refl[0];	// to facilitate calculation of total attenuation
  maxcoeffs= BOXREV_MINMAXC;
  coeffs= malloc( (2*maxcoeffs+1) * sizeof(float) );
  *coeffs= END;
  rng= makerng();
  mode[0]= 1;
  depthdelay[0]= basedelay;
  dlysqr[-1]= 0.0f;
  dlysqr[0]= basedelay*basedelay;
  atten[-1]= atten[0]= 1.0f;
  depth= 0;
  thisdelay= basedelay;
  while( depth >= 0 )
  {
    while( depth< maxdepth ) {
      ++depth;
      mode[depth]= 0;
      depthdelay[depth]= 0;
      dlysqr[depth]= dlysqr[depth-1];
      atten[depth]= atten[depth-1];
    }

    boxrev_addcoeff( &coeffs, &maxcoeffs, thisdelay - basedelay,
			atten[depth]*basedelay/thisdelay );

    while( depth >= 0 ) {
      if( depth==0 )
	mode[depth] += 2;   // sound has to travel twice the length of the room before it reaches listener again
      else
	++mode[depth];
      depthdelay[depth] += delay[depth] * (1.0f + random*(float)gaussrng(rng));
      dlysqr[depth]= dlysqr[depth-1] + depthdelay[depth]*depthdelay[depth];
      thisdelay= sqrtf(dlysqr[depth]);
      atten[depth] *= refl[depth];
      if( cutoff < 0 )	// cut off according to amplitude
	if( fabsf(atten[depth])/sqrtf(dlysqr[depth]) > -cutoff )
	  break;	// still over the threshold -> continue at this depth
	else;
      else	// cut off according to delay
	if( thisdelay < cutoff )
	  break;	// still below maximum delay -> continue at this depth
      --depth;
    }
  }
  free(rng);
  if( *coeffs == END )
    return signal;
  p= filter(coeffs, NULL, signal);
  p->name= p->type= "boxrev";
  p->private[PLENTY-1]= coeffs;
  return p;
}

int boxrev_addcoeff( float **coeffs, int *maxcoeffs, float delay, float coefficient )
{
  float *search, *newcoeffs;
  int newsize;

//  printf("addcoeff: %g\t%g\n", (double)delay, (double)coefficient );
  if( coefficient == 0.0f )
    return 0;
  delay= floorf(delay*(float)SAMPRATE + 0.5f);
  for( search= *coeffs; *search != END; search += 2 )
    if( *search == delay ) {
      search[1] += coefficient;
      return 0;
    }
  if( search - *coeffs == 2 * *maxcoeffs )
  {
    if( *maxcoeffs == BOXREV_MAXMAXC )
      return -1;
    newsize= 2 * *maxcoeffs;
    if( newsize > BOXREV_MAXMAXC )
      newsize= BOXREV_MAXMAXC;
    newcoeffs= malloc( (2*newsize+1) * sizeof(float) );
    if( !newcoeffs )
      return 1;
    memcpy(newcoeffs, *coeffs, 2 * *maxcoeffs * sizeof(float));
    // no need to copy END  - we are going to overwrite it
    search += newcoeffs - *coeffs;
    free(*coeffs);
    *coeffs= newcoeffs;
    *maxcoeffs= newsize;
  }
  *search++ = delay;
  *search++ = coefficient;
  *search++ = END;
  return 0;
}


/*=============================================================================
    convrev (revtime) (refltime) (directtime) (direct2noise) (refl2noise) <mixratio> <signal>
    revkernel (revtime) (refltime) (directtime) (direct2noise) (refl2noise) (stereo)

Reverberation by convolution with a kernel.  The kernel is generated by the
sndobj revkernel.  Its arguments are the reverberation time, the total duration
of all reflections and the duration of the direct reflections, the relative
volume of the direct reflections w.r.to the ambient noise, and the relative
volume of the indirect reflections.  The reflections are created using the 
"slide" sound object, the ambient noise with shelve-filtered white noise.
The convrev object has as additional arguments two inputs, the mixing ratio
between dry signal and reverberation (1 - only signal, 0 - only reverb),
and the signal.  If the signal has an even number of channels, it is taken to
be a collection of stereo signals.  Then the odd and even channels are
convoluted with kernels generated from different pseudo-random number series.

See also: freeverb, nrev, josrev, fdnrev, boxrev
=============================================================================*/

sndobj *convrev( double revtime, double refltime, double directtime, double direct2noise, double refl2noise, sndobj *mixratio, sndobj *signal )
{
  return mix2(mixratio, signal, convolve(SAMPRATE*revtime, 0,
		lowpass(c(6000), revkernel(revtime, refltime, directtime,
		direct2noise, refl2noise, !(signal->nch&1))), signal) );
}

sndobj *revkernel( double revtime, double refltime, double directtime, double direct2noise, double refl2noise, int stereo )
{
  sndobj *p, *purenoise, *noise, *refl;
  double *rampdata;
  double norm;

  if( refltime < directtime )
    refltime= directtime;
  if( revtime < refltime )
    revtime= refltime;
  if( direct2noise > refl2noise )
    norm= 1.0/(1.0 + direct2noise);
  else
    norm= 1.0/(1.0 + refl2noise);
  rampdata= news(double, 27);
  // amplitude shaping of reflections:
  rampdata[0]= 0;
  rampdata[1]= norm * direct2noise;
  rampdata[2]= RAMP_C;
  rampdata[3]= directtime;
  rampdata[4]= norm * refl2noise;
  rampdata[5]= RAMP_L;
  rampdata[6]= refltime;
  rampdata[7]= 0;
  rampdata[8]= END;
  // amplitude shaping of noise:
  rampdata[9]= 0;
  rampdata[10]= 0;
  rampdata[11]= RAMP_L;
  rampdata[12]= directtime;
  rampdata[13]= norm * 37;
  rampdata[14]= RAMP_L;
  rampdata[15]= revtime;
  rampdata[16]= 0;
  rampdata[17]= END;
  // common amplitude shaping:
  rampdata[18]= 0;
  rampdata[19]= 1;
  rampdata[20]= RAMP_C;
  rampdata[21]= directtime;
  rampdata[22]= 1;
  rampdata[23]= RAMP_E;
  rampdata[24]= revtime;
  rampdata[25]= 0;
  rampdata[26]= END;
  if( stereo ) {
    refl= mul(2, ramp(0, rampdata), mix("0.0*1 + 1.0*-1; 2.0*1 + 3.0*-1",
		slide(c(200), c(600), c(0.5), c(0.7), c(0.4)), 
		slide(c(200), c(600), c(0.5), c(0.7), c(0.4)),
		slide(c(200), c(600), c(0.5), c(0.7), c(0.4)), 
		slide(c(200), c(600), c(0.5), c(0.7), c(0.4)) ) );
    purenoise= switchboard( "0 1", flatrand(0.05), flatrand(0.05) );
  }
  else {
    refl= mul(2, ramp(0, rampdata), arithmetic("-", 
		slide(c(200), c(600), c(0.5), c(0.7), c(0.4)), 
		slide(c(200), c(600), c(0.5), c(0.7), c(0.4)) ));
    purenoise= flatrand(0.05);
  }
  noise= mul(2, ramp(0, rampdata+9), 
    mix2(c(0.66), shelve( c(50), c(10), shelve( c(3000), c(0.15), purenoise ) ),
      		shelve( c(10000), c(10), purenoise ) ) );
  p= mul(2, ramp(0, rampdata+18), add(2, noise, refl));
  p->private[PLENTY-1]= rampdata;
  return p;
}



/*=============================================================================
    clip (quench) <max> <signal>

Clip all channels of second input by modulus of corresponding channel of first.
If the second input has fewer channels, its channel number wraps around.
<quench> is the factor by which the parts of the signal which exceed the limit
are multiplied (0.0=strict clipping).  If quench is <0, the parts of the signal
which are within limits are replaced by the limit (but keep their sign).

See also: softclip, crossoverdist, lee, conditional, limiter
=============================================================================*/

float calcclip(sndobj *me, sndobj *, int);

sndobj *clip( float quench, sndobj *max, sndobj *signal )
{
  sndobj *p;
  float *savequench;
  char *namestr;
  
  namestr= (char *)malloc(NAMESTRLEN*sizeof(char));
  snprintf(namestr, NAMESTRLEN, "q=%g", quench );
  p= newsndo(calcclip, namestr, "clip", signal->nch, 2, max, signal );
  p->private[0]= savequench= new(float);
  *savequench= quench;
  p->private[1]= namestr;
  return p;
}

float calcclip(sndobj *me, sndobj *caller, int innr)
{
  float quench, limit, sigmod, sigsign;
  int ch, lch;
  
  quench= *(float*)me->private[0];
  INCALC(0);
  INCALC(1);
  for( ch= 0, lch= 0; ch< me->nch; ++ch, ++lch ) {
    if( lch>=me->in[0]->nch )
      lch -= me->in[0]->nch;
    limit= fabsf(INPUTC(0, lch));
    sigmod= INPUTC(1, ch);
//if( ch==1 )
//  printf("channel 1 clip: limit= %g, signal= %g\n", limit, sigmod );
    if( sigmod<0.0 ) {
      sigmod= -sigmod;
      sigsign= -1.0;
    }
    else sigsign= 1.0;
    if( quench >= 0.0 )
      if( sigmod > limit )
	OUTPUT(ch, sigsign*(limit + (sigmod-limit)*quench));
      else
	OUTPUT(ch, INPUTC(1, ch));
    else
      if( sigmod < limit )
	OUTPUT(ch, sigsign*limit);
      else
	OUTPUT(ch, INPUTC(1, ch));
  }
  CALCRETURN;
}


/*=============================================================================
    softclip <max> <signal>

Clip all channels of second input by modulus of corresponding channel of first.
If the first input has fewer channels, its channel number wraps around.
Unlike clip, this clipping is soft, ie the waveform is rounded off towards the
limit by replacing each in-range sample value y by y*(1.5 - 0.5*(y/<max>)^2).

See also: clip, crossoverdist, lee, limiter, arithmetic
=============================================================================*/

float calcsoftclip(sndobj *me, sndobj *, int);

sndobj *softclip( sndobj *max, sndobj *signal )
{
  return newsndo( calcsoftclip, "softclip", "softclip", signal->nch, 2, max, signal );
}

float calcsoftclip(sndobj *me, sndobj *caller, int innr)
{
  float limit, normval;
  int ch, lch;
  
  INCALC(0);
  INCALC(1);
  for( ch= 0, lch= 0; ch< me->nch; ++ch, ++lch ) {
    if( lch>=me->in[0]->nch )
      lch -= me->in[0]->nch;
    limit= fabsf(INPUTC(0, lch));
    if( !limit ) {
      OUTPUT(ch, 0.0);
      continue;
    }
    if( fabsf(INPUTC(1, ch)) > limit ) {
      if( INPUTC(1, ch) > 0 )
	OUTPUT(ch, limit);
      else
	OUTPUT(ch, -limit);
    }
    else {
      normval= INPUTC(1, ch)/limit;
      OUTPUT(ch, INPUTC(1, ch) * (1.5f - 0.5f*normval*normval));
    }
  }
  CALCRETURN;
}


/*=============================================================================
    lee <freq> <noiselvl> <signal>

Lee's statistical filter.  Intended for noise removal on silent stretches but
also suitable for distortion.  If the RMS of the <signal> exceeds the noise
level <noiselvl>, the <signal> is mixed with a signal average with a prefactor
of (1 - (noise level)^2 / RMS of signal^2) for the average.  <freq> gives the
cutoff frequency of both average and RMS.

See also: clip, softclip, crossoverdist
=============================================================================*/

/* from the avg object: */
struct average {
  double    addsum[PLENTY], subsum[PLENTY];
  int       nch, maxind, subcountm1, nonorm, firstcall;
  float     freq;
};
struct lee {
  struct average *data, *datasqr;
};

float calclee(sndobj *me, sndobj *, int);

sndobj *lee( sndobj *freq, sndobj *noiselvl, sndobj *signal )
{
  sndobj *p;
  struct lee *d;
  
  d= new(struct lee);
  p= newsndo( calclee, "lee", "lee", signal->nch, 4, freq, noiselvl, signal, 
						    mul(2, signal, signal) );
  p->private[0]= d;
  p->private[1]= d->data= new(struct average);
  p->private[2]= d->datasqr= new(struct average);
  d->data->nch= p->nch;
  d->data->nonorm= 0;
  d->data->firstcall= 1;
  d->datasqr->nch= p->nch;
  d->datasqr->nonorm= 0;
  d->datasqr->firstcall= 1;
  buf( HORIZON/2.0, HORIZON/2.0 + 1.0/(double)SAMPRATE, 0.0, p->nch, p, 2 );
  buf( HORIZON/2.0, HORIZON/2.0 + 1.0/(double)SAMPRATE, 0.0, p->nch, p, 3 );
  return p;
}

/* from the avg object: */
float average( sndobj *me, int innr, struct average *a, float newfreq, float *dest );

float calclee(sndobj *me, sndobj *caller, int innr )
{
  static float mean[PLENTY], meansqr[PLENTY];
  
  struct lee *d;
  float noiselvl, fade;
  int ch;
  
  d= (struct lee *)me->private[0];
  average( me, 2, d->data, INCALC(0), mean );
  average( me, 3, d->datasqr, INPUT(0), meansqr );
  INCALC(0);
  noiselvl= INCALC(1);
  BUFINCALC(2, 0, 0.0);
  FORCH
    if( meansqr[ch] <= noiselvl*noiselvl )
      OUTPUT(ch, mean[ch]);
    else {
      fade= (meansqr[ch] - noiselvl*noiselvl)/meansqr[ch];
      OUTPUT(ch, fade*INPUTC(2, ch) + (1.0-fade)*mean[ch]);
    }
  BUFINCR(2, 1);
  BUFINCR(3, 1);
  CALCRETURN;
}


/*=============================================================================
    crossoverdist <crossover amplitude> <smoothing> <signal>

Based on crossover distortion by Steve Harris, which I didn't manage to port
one-to-one (different conventions?).  Apparently "a simulation of the
distortion that happens in class B and AB power amps when the signal crosses
0".

The crossover amplitude should be much smaller than the signal.  The smoothing
value means: 0 no smoothing; 1 really smooth signal; <0 and >>1 really nasty.

The output is = +/- (|input| - ampl) if |input| > ampl; 
+/- (1 - smooth) * (|input| - ampl) otherwise (+/- = sign(input)).
The smoothing value input is only evaluated if needed.

See also: clip, softclip, lee
=============================================================================*/

float calccrossoverdist(sndobj *me, sndobj *, int);

sndobj *crossoverdist( sndobj *ampl, sndobj *smooth, sndobj *signal )
{
  return newsndo(calccrossoverdist, "crossoverdist", "crossoverdist", 
		    signal->nch, 3, ampl, smooth, signal );
}

float calccrossoverdist(sndobj *me, sndobj *caller, int innr)
{
  float ampl, smooth, sigmod;
  int ch, smoothneeded;

  ampl= fabsf(INCALC(0));
  smoothneeded= 0;
  INCALC(2);
  FORCH {
    sigmod= fabsf(INPUTC(2, ch));
    if( sigmod > ampl )
      if( INPUTC(2, ch) >= 0 )
        OUTPUT(ch, sigmod - ampl);
      else
        OUTPUT(ch, -(sigmod - ampl));
    else {
      if( !smoothneeded ) {
        smooth= 1.0f - INCALC(1);
	smoothneeded= 1;
      }
      if( INPUTC(2, ch) >= 0 )
        OUTPUT(ch, smooth * (sigmod - ampl));
      else
        OUTPUT(ch, -smooth * (sigmod - ampl));
    }
  }
  if( !smoothneeded )
    INSKIP(1);
  CALCRETURN;
}



/*=============================================================================
    flanger (coeff) (backfeed) (maxdelay) <del> <signal>

Sum of signal and signal delayed by <del> (maximum= (maxdelay)), scaled by
(coeff) < 1.  If (backfeed)!=0, the delayed but not scaled signal is fed back
to the input, scaled by (backfeed).  This leads to resonances at frequencies
which are integral multiples of 1/del, equidistant on the frequency axis.

See also: phaser
=============================================================================*/

sndobj *flanger( float coeff, float backfeed, float maxdelay, sndobj *del, sndobj *signal )
{
  sndobj *bf, *delayed;
  
  if( backfeed ) {
    bf= add(1, signal);
    delayed= delay(maxdelay, del, bf);
    addinput( bf, linear( 0.0, backfeed, delayed));
    return add(2, bf, linear(0.0, coeff, delayed));
  }
  else
    return add(2, signal, linear(0.0, coeff, delay(maxdelay, del, signal)));
}


/*=============================================================================
    phaser (nfreq) (directgain) <basefreq> <freqfactor> <signal>

Sum of signal scaled by (directgain) and signal passed through (nfreq) allpasses
with "break" (=notch) frequencies at <basefreq> * <freqfactor>^n.  Similar to
flanger, but the notches are uniformly spaced on the logarithmic frequency
axis.

See also: flanger
=============================================================================*/

struct phaser {
  sndobj    **allpasses;
  int       nfreq;
  float     currbasef, currfactor;
};

float calcphaser(sndobj *me, sndobj *, int);

#define PHASER_MAXNF    100
sndobj *phaser( int nfreq, float directgain, sndobj *basefreq, sndobj *freqfactor, sndobj *signal )
{
  sndobj *p, *ap;
  struct phaser *d;
  int i;

  if( nfreq<1 )
    nfreq= 1;
  if( nfreq > PHASER_MAXNF ) {
    fprintf(stderr, "Warning: phaser: number of frequencies (%d) too large. Reduced to maximum (%d).\n", nfreq, (int)PHASER_MAXNF );
    nfreq= PHASER_MAXNF;
  }
  p= newsndo( calcphaser, signal->name, "phaser", signal->nch, 2, basefreq, freqfactor );
  p->private[0]= d= new(struct phaser);
  p->private[1]= d->allpasses= news( sndobj *, nfreq );
  d->nfreq= nfreq;
  d->currbasef= 0.0;
  d->currfactor= 1.0;
  ap= signal;
  for( i= 0; i< nfreq; ++i )
    d->allpasses[i]= ap= _foallpass( 1.0, ap );
  addinput( p, add(2, ap, linear(0.0, directgain, signal)) );
  return p;
}

float calcphaser(sndobj *me, sndobj *caller, int innr)
{
  struct phaser *d;
  float *forw, *back;
  float newbasef, newfactor;
  double fnorm, tanval;
  int i, ch;
  
  d= (struct phaser *)me->private[0];
  newbasef= INCALC(0);
  newfactor= INCALC(1);
  if( newbasef != d->currbasef || newfactor != d->currfactor ) {
    if( newbasef==d->currbasef )
      i= 1;
    else
      i= 0;
    d->currbasef= newbasef;
    d->currfactor= newfactor;
    fnorm= M_PI*newbasef/SAMPRATE;
    for( ; i< d->nfreq; ++i ) {
      forw= (float*)d->allpasses[i]->private[PLENTY-1];
      back= forw+5;
      tanval= tan( fnorm );
      forw[1]= back[1]= (float)((1.0-tanval)/(1.0+tanval));
      fnorm *= d->currfactor;
    }
  }
  OUTPUT(0, INCALC(2));
  for( ch= 1; ch < me->nch; ++ch )
    OUTPUT(ch, INPUTC(2, ch));
  CALCRETURN;
}


/*=============================================================================
    leslie <freq> <strength> <signal>

Classical Hammond organ effect caused by phase modulation with a sine.  Useful
range for freq: 1..10 Hz; for strength: 0.0001..0.001.
=============================================================================*/

sndobj *leslie( sndobj *freq, sndobj *strength, sndobj *signal )
{
  return delay(1.0, sine(0.0, freq, strength), signal );
}


/*=============================================================================
    chirp (chirptime) <signal>

Effect simulating transmission on a dispersive medium exhibiting bending
vibrations, on which high frequencies arrive first.  chirptime is the delay
between the fastest and the slowest frequency.  For chirptime < 0, low
frequencies arrive first.
=============================================================================*/

sndobj *chirp( float chirptime, sndobj *signal )
{
  sndobj *p;
  float *chirpfilter;
  double cdecay, delaylen, baselen;
  int pulsecount, inversechirp;
  
  if( chirptime < 0 ) {
    inversechirp= 1;
    chirptime= -chirptime;
  }
  else
    inversechirp= 0;
  if( chirptime <= 1/22.05 )
    chirptime= 1/22.05;
  if( inversechirp ) {
    cdecay= -chirptime*SAMPRATE/3/M_LN10;
    baselen= 2000;
  }
  else {
    // tau= (tmax-lambdamax)*log(lambdamax/lambdamin)
    cdecay= (chirptime-1/22.05)*SAMPRATE/3/M_LN10;
    baselen= 2;
  }
  for( pulsecount= 1, delaylen= baselen; delaylen< (long)((float)SAMPRATE*chirptime); ++pulsecount )
    delaylen += baselen*exp((double)delaylen/cdecay);
  chirpfilter= news(float, 2*pulsecount+1);
  chirpfilter[0]= 0.0;
  chirpfilter[1]= 1.0/(float)pulsecount;
  for( pulsecount= 2, delaylen= baselen; delaylen< (long)((float)SAMPRATE*chirptime); ) {
    chirpfilter[pulsecount++]= delaylen;
    chirpfilter[pulsecount++]= chirpfilter[1];
    delaylen += baselen*exp((double)delaylen/cdecay);
  }
  chirpfilter[pulsecount]= END;
  p= filter( chirpfilter, NULL, signal );
  p->private[PLENTY-1]= chirpfilter;
  return p;
}


/*=============================================================================
    spike (nsam) <maxamp> <signal>

Multiplies (nsam) successive samples.  If (nsam)<0, the signal value's sign is
kept, otherwise the modulus is taken.  <maxamp> gives the value which remains
unchanged and thereby determines the scaling of the output.

See also: arithmetic
=============================================================================*/

float calcspike(sndobj *me, sndobj *, int);

sndobj *spike( int nsam, sndobj *maxamp, sndobj *signal )
{
  sndobj *p;
  
  if( abs(nsam)<2 ) {
    fprintf(stderr, "Warning: spike multiplies at least two samples (argument: %d). Returning original signal or c(1).\n", nsam );
    if( nsam )
      return signal;
    else
      return c(1);
  }
  p= newsndo( calcspike, signal->name, "spike", signal->nch, 2, maxamp, signal );
  p->private[0]= new(int);
  *(int*)p->private[0]= nsam;
  buf( (float)abs(nsam)/SAMPRATE, 0.0, 0.0, signal->nch, p, 1 );
  return p;
}

float calcspike(sndobj *me, sndobj *caller, int innr)
{
  float norm;
  int i, ch;
  
  FORCH
    me->ch[ch]= 1.0;
  for( i= abs(*(int*)me->private[0])-1; i>= 0; --i ) {
    BUFINCALC(1, -i, 0.0);
    FORCH
      me->ch[ch] *= INPUTC(1, ch);
  }
  norm= INCALC(0);
  for( i= abs(*(int*)me->private[0])-1; i> 0; --i )
    norm *= INPUT(0);
  if( *(int*)me->private[0] > 0 )
    FORCH
      me->ch[ch]= fabsf(me->ch[ch])/norm;
  else
    FORCH
      if( INPUTC(1, ch) >= 0.0 )
	me->ch[ch]= fabsf(me->ch[ch])/norm;
      else
	me->ch[ch]= -fabsf(me->ch[ch])/norm;
  BUFINCR(1, 1);
  CALCRETURN;
}


/*=============================================================================
    fbmdistort <lacunarity> <H> <cutoff> <source>

Similar to fractional Brownian motion noise generator, but the midpoint
displacement at each level is given by the corresponding channel of the input
<source>.  This object is experimental and did not quite work as I expected.

See also: fbm
=============================================================================*/
 
struct fbmdistort {
  float     *currval, *prevval, *totalscale;
  double    *gridpoint, *gridspacing;
  long      currmaxdepth, count;
  int       firstcall;
};

float calcfbmdistort(sndobj *me, sndobj *, int);

#define FBMD_MINCUTOFF  10
#define FBMD_MAXLACU    0.9

sndobj *fbmdistort(sndobj *lacunarity, sndobj *H, sndobj *cutoff, sndobj *source, sndobj *reset)
{
  sndobj *p;
  struct fbmdistort *d;
  long maxdepthp1;

  p= newsndo( calcfbmdistort, "fbm", "fbm", 1, 4, lacunarity, H, cutoff, source );
  if( reset )
    addinput(p, reset);
  p->private[0]= d= new(struct fbmdistort);
  maxdepthp1= (long) ceil(- log(0.5*(double)SAMPRATE/(double)FBMD_MINCUTOFF)/log((double)FBMD_MAXLACU) ) + 1L;
  p->private[1]= d->currval= news(float, maxdepthp1);
  p->private[2]= d->prevval= news(float, maxdepthp1);
  p->private[3]= d->totalscale= news(float, maxdepthp1);
  p->private[4]= d->gridpoint= news(double, maxdepthp1);
  p->private[5]= d->gridspacing= news(double, maxdepthp1);
  d->firstcall= 1;
  return p;
}

float calcfbmdistort(sndobj *me, sndobj *caller, int innr)
{
  struct fbmdistort *d;
  double lacunarity, H, cutoff;
  float scale;
  long depth;

  d= (struct fbmdistort *)me->private[0];
  lacunarity= INCALC(0);
  H= INCALC(1);
  cutoff= INCALC(2);
  INCALC(3);
  if( lacunarity < 0.0 )
    lacunarity= 0.0;
  else if( lacunarity > FBMD_MAXLACU )
    lacunarity= FBMD_MAXLACU;
  if( H < 0.0 )
    H= 0.0;
  else if( H > 1.0 )
    H= 1.0;
  if( cutoff< FBMD_MINCUTOFF )
    cutoff= FBMD_MINCUTOFF;

  if( (me->nin>4 && INCALC(4)) || d->firstcall )    /* ! left to right ! */
  {
    depth= 0;
    scale= M_SQRT1_2 * pow(lacunarity, H) * sqrt( 1.0 - pow(lacunarity, 2.0-2.0*H) );
    d->currval[0]= scale*INCALC(3);
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
    scale= M_SQRT1_2 * pow(lacunarity, H) * sqrt( 1.0 - pow(lacunarity, 2.0-2.0*H) );
    d->currval[0]= scale*INPUT(3);
    d->totalscale[0]= scale;
    ++depth;
  }
  
  for( ; d->gridspacing[depth-1] > 1.0; ++depth )
  {
    d->gridspacing[depth]= d->gridspacing[depth-1] * lacunarity;
    scale= M_SQRT1_2 * pow(lacunarity, (double)(depth+1) * H) * 
		    sqrt( 1.0 - pow(lacunarity, 2.0-2.0*H));
    if( d->gridspacing[depth] < 1.0 ) {
      scale *= (d->gridspacing[depth]-lacunarity)/(1.0-lacunarity);
      d->gridspacing[depth]= 1.0;
    }
    d->totalscale[depth]= d->totalscale[depth-1] + scale;
    if( depth > d->currmaxdepth )
      d->gridpoint[depth]= (double)d->count;
    else
      d->gridpoint[depth] += d->gridspacing[depth];
    d->prevval[depth]= d->currval[depth];
    if( d->gridpoint[depth]>=d->gridpoint[depth-1] )
      if( depth< PLENTY )
        d->currval[depth]= d->currval[depth-1] + scale*INPUTC(3, depth);
      else
        d->currval[depth]= d->currval[depth-1] + scale*INPUTC(3, PLENTY-1);
  }

  /*  divide by 3.5 standard deviations to normalise:  */
  OUTPUT(0, d->currval[d->currmaxdepth]/(3.5*sqrt(d->totalscale[d->currmaxdepth])));
  ++d->count;
  if( !d->totalscale[d->currmaxdepth] )
    printf ("zero sigma_tot at %ld\n", d->count );
  CALCRETURN;
}


/*=============================================================================
    avghw (nmax) <n> <signal>

Averages the shape of a number of successive half-waves of the signal.  A
half-wave is the time interval between two zero-crossings.  This snd object
only works on one channel.  The first input, <n>, gives the number of
half-waves to average.  The parameter <nmax> is the maximal value <n> can take.

See also: reshape, repeathw
=============================================================================*/

struct avghw {
  int	    first, last, ncurr, nmax, samcount;
  float	    *hwmax, scale;
  double    *hwstart, *ind, *incr;
};

float calcavghw(sndobj *me, sndobj *, int);

sndobj *avghw(int nmax, sndobj *n, sndobj *signal)
{
  sndobj *p;
  struct avghw *d;
  int count;
  
  p= newsndo(calcavghw, "avghw", "avghw", 1, 2, n, signal);
  p->skip= dontskip;
  buf(HORIZON, HORIZON, HORIZON, 1, p, 1);
  p->private[0]= d= new(struct avghw);
  d->nmax= nmax;
  p->private[1]= d->hwstart= news(double, nmax+1);
  p->private[2]= d->hwmax= news(float, nmax+1);
  p->private[3]= d->ind= news(double, 2*nmax);
  d->incr= d->ind + nmax;
  d->ncurr= 1;
  d->first= 0;
  d->last= nmax/2;
  for( count= 0; count< nmax; ++count ) {
    d->ind[count]= 0.0;
    d->incr[count]= 0.0;
  }
  for( count= 0; count<= nmax/2+1; ++count )
    d->hwstart[count]= 0.0;
  for( ; count<= nmax; ++count )
    d->hwstart[count]= MAGIC;
  for( count= 0; count<= nmax; ++count )
    d->hwmax[count]= 1.0;
  d->samcount= 0;
  return p;
}

float calcavghw(sndobj *me, sndobj *caller, int innr)
{
  struct avghw *d;
  double index;
  float outval;
  int hw;
  
  d= (struct avghw *)me->private[0];
  if( !d->samcount )
  {
    float max, previous;
    int shift, searchind;
    
    shift= (int)ceil(d->hwstart[d->nmax/2+1]) - (int)ceil(d->hwstart[d->nmax/2]);
    BUFINCR(1, shift);
    for( hw= 0; hw <= d->last && d->hwstart[hw+1] != MAGIC; ++hw )  {
      d->hwstart[hw]= d->hwstart[hw+1] - shift;
      d->hwmax[hw]= d->hwmax[hw+1];
    }
    d->ncurr= (int)floorf(INCALC(0)+0.5);
    if( d->ncurr> d->nmax )
      d->ncurr= d->nmax;
    d->first= d->nmax/2 - d->ncurr/2;
    d->last= d->first + d->ncurr - 1;
    // Now search for boundaries of more half-waves as necessary.  Search for
    // non-zero value first - only necessary at the beginning of processing or
    // if the previous search hit the HORIZON.
    for( searchind= (int)ceil(d->hwstart[hw]) - shift; 
		searchind< (int)(SAMPRATE*HORIZON); ++searchind ) {
      BUFINCALC(1, searchind, 0.0);
      if( INPUT(1)!=0.0 )
        break;
    }
    previous= INPUT(1);
    for( ; hw <= d->last+1; ++hw ) {
      if( previous> 0.0 ) {
        max= INPUT(1);
	for( ; searchind< (int)(SAMPRATE*HORIZON); ++searchind ) {
	  BUFINCALC(1, searchind, 0.0);
	  if( INPUT(1) < 0.0 )
	    break;
	  previous= INPUT(1);
	  if( INPUT(1) > max )
	    max= INPUT(1);
	}
      }
      else {
        max= INPUT(1);
	for( ; searchind< (int)(SAMPRATE*HORIZON); ++searchind ) {
	  BUFINCALC(1, searchind, 0.0);
	  if( INPUT(1) > 0.0 )
	    break;
	  previous= INPUT(1);
	  if( INPUT(1) < max )
	    max= INPUT(1);
	}
      }
      d->hwstart[hw]= searchind - (double)INPUT(1) / ((double)INPUT(1) - previous);
      d->hwmax[hw-1]= max? max: 1.0;
      previous= INPUT(1);
    }
    d->samcount= (int)ceil(d->hwstart[d->nmax/2+1]) - (int)ceil(d->hwstart[d->nmax/2]);
    for( hw= d->first; hw<= d->last; ++hw ) {
      d->incr[hw]= (d->hwstart[hw+1] - d->hwstart[hw]) /
		    (d->hwstart[d->nmax/2+1] - d->hwstart[d->nmax/2]);
      d->ind[hw]= d->hwstart[hw] + (ceil(d->hwstart[d->nmax/2])-d->hwstart[d->nmax/2]) * d->incr[hw];
    }
    d->scale= d->hwmax[d->nmax/2] / (d->last-d->first+1);
  }
  else
    INSKIP(0);
  outval= 0;
  for( hw= d->first; hw<= d->last; ++hw ) {
    index= d->ind[hw];
    BUFINCALC(1, (int)floor(index), index-floor(index));
    outval += INPUT(1) / d->hwmax[hw];
    d->ind[hw] += d->incr[hw];
  }
  OUTPUT(0, d->scale * outval);
  --d->samcount;
  CALCRETURN;
}



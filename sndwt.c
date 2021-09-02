/*
    This file is part of sndsys, a digital signal processing system.
    Copyright (c) 2005-07 Volker Schatz (noise at volkerschatz dot com).

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

/*
  This source file contains objects performing wavelet transforms and working
  on wavelet-transformed signals.  The wavelet transform divides the signal up
  into several frequency bands comprising an octave (a factor of two) on the
  frequency scale each.  In addition, there is a band containing the
  frequencies remaining below the lowest band.  The boundary between the lowest
  wavelet band and the remainder is computed from the parameter <bottom> which
  all wavelet objects need.  <bottom> is approximately the frequency bounding
  the lowest wavelet band from the remainder.  The real bounding frequency is
  the one closest to <bottom> on the logarithmic frequency axis.  The number of
  wavelet bands for a given bottom frequency is computed by wtdepth(), which
  should be used by all objects dealing with wavelet signals.  After taking the
  modulus of <bottom> and limiting it to at least MINBOTTOM, wtdepth()
  calculates floor(log(SAMPRATE/bottom)/M_LN2+0.5).  MINBOTTOM can be changed
  below and is currently set to 20 Hz, which puts the top end of the remainder
  band between 14.14 and 28.3 Hz.

  An important feature of wavelet transforms is that (unlike the Fourier
  transform) it retains temporal information of the frequency-space data.
  There are twice as many data corresponding to a given band as corresponding
  to the next lower band.  In sndsys, the wavelet signals are organised in data
  sets 2^depth in length.  Half of them represent the top band, a quarter the
  next-highest and so on.  The location of the data within the data set is such
  that every wavelet datum is as close as possible to the centre of the
  interval of the source data from which it was computed.  This is achieved
  with the following representation: the first, third and so on value
  represents the highest frequency band or smallest scale; the second, sixth,
  tenth value and so on the next larger scale; down to the last value in the
  data set which represents the largest scale.  Here is a graphical
  representation of a data set for depth=4:

  S S S S S S S S S S S S S S S S
  |/  |/  |/  |/  |/  |/  |/  |/  
  - + - + - + - + - + - + - + - +
    | /     | /     | /     | /  
    - +     - +     - +     - +
       \ ,----'        \ ,----'
        - +             - +
           `---. ,--------'
	        - + ----------> +

  "S" are the original sample values, "-" are the "detail" coefficients, and
  "+" are the smoothed values resulting from a convolution with the scaling
  function.  Every "+" is computed from the same original values as the "-" to
  its left, only it represents the smooth part as opposed to the detail at its
  scale.  The lowest "+" or "-" in each column is the wavelet coefficient
  output in lieu of the input sample value in the top row.  There is only one
  "+" or smoothed value which is output, the last value in the set which
  represents the lowest frequency band.  The others just serve as intermediate
  stages in the computation.  The one smooth output value represents the
  "remainder" band.  For Haar wavelets, each "-" (and the "+" to its right) is
  computed as the difference between the two smoothed values at the next
  smaller scale with which it is connected.  For other wavelets, it results
  from convolving several (usually more than two) smoothed values with the
  scaling filter (which is a free parameter of wtany, and chosen for a specific
  wavelet in the other transforms).  For the predefined wavelets except wtdaub,
  the filter is centred between the two data on which the haar transform is
  computed, which makes for good temporal coincidence between the wavelets and
  the original signal.

  Manipulating wavelet-transformed signals tends to create distorted, metallic
  sounds.  This is because if you do "the same thing" to the series of wavelet
  coefficients representing different bands, you in effect apply a filter with
  a different time constant to each.  Because the high-frequency values are
  more frequent, their associated time scale is smaller.  This is similar to
  wave propagation in metals, in which the speed of sound is larger for high
  frequencies.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "sndsys.h"

#define WTMINBOTTOM	20.0

//======================================================================
//  Auxiliary functions

int wtdepth( double bottom );
int floorlog2(int power);
int pow2div(int num);
float bitinverse(unsigned num);

int wtdepth( double bottom )
{
  bottom= fabs(bottom);
  if( bottom < WTMINBOTTOM )
    bottom= WTMINBOTTOM;
  return (int)floor(log((double)SAMPRATE/bottom)/M_LN2 + 0.5);
}


// Floor value of 2's logarithm of an integer; equivalent to index of most
// significant non-zero bit
int floorlog2(int power)
{
  static int nibblelog[]= { 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 };
  
  int thelog= 0;
  
//  printf("floorlog2 called with argument %d\n", power);
  if( power <= 0 )
    return 0;
  if( power >= 0x10000 ) {
    thelog += 16;
    power >>= 16;
  }
  if( power >= 0x100 ) {
    thelog += 8;
    power >>= 8;
  }
  if( power >= 0x10 ) {
    thelog += 4;
    power >>= 4;
  }
//  printf("returning %d\n", thelog + nibblelog[power] );
  return thelog + nibblelog[power];
}


// Exponent of largest power of 2 which divides a number; or power of 2 in a
// number's prime decomposition; in other words, the number of zeros at its end
// when written in binary notation.
int pow2div(int num)
{
  int result;
  
  if( !num )
    return -1;
  for(result= 0; (num & 1)==0; ++result, num >>= 1 );
  return result;
}

// Converts a binary number to a fraction between 0 and 1 by reversing the
// order of the bits and treating them as bits after the decimal point.
float bitinverse(unsigned num)
{
  float result, thisbit;

  result= 0.0f;
  for( thisbit= 0.5f; num; thisbit /= 2.0f, num >>= 1 )
    if( (num&1) != 0 )
      result += thisbit;
  return result;
}



/*=============================================================================
    wthaarnn (bottom) <signal>

Fast wavelet transform with Haar wavelets.  (bottom) is the minimum frequency
(in Hertz) to be represented, corresponding to the upper limit of the lowest
frequency band.  The wavelet coefficients generated by this object are
not normalised (hence wthaar"nn").  This allows a faster computation (no
multiplications with 1/sqrt(2)) but means they are not compatible with
normalised Haar wavelets which can be computed with wthaar.

See also: iwthaarnn, wthaar, wtany
=============================================================================*/

struct wthaarnn {
  int	depth, setlen, setcount;
  float *output;
};

float calcwthaarnn(sndobj *me, sndobj *, int);

sndobj *wthaarnn(double bottom, sndobj *signal)
{
  sndobj *p;
  struct wthaarnn *d;
  int depth;
  
  if( (depth= wtdepth(bottom)) == 0 )
    return signal;
  p= newsndo( calcwthaarnn, "wthaarnn", "wthaarnn", signal->nch, 1, signal );
  p->skip= dontskip;
  p->private[0]= d= new(struct wthaarnn);
  d->depth= depth;
  d->setlen= 1 << d->depth;
  d->setcount= d->setlen;
  buf( 0.0, (double)d->setlen/(double)SAMPRATE, 0.0, signal->nch, p, 0 );
  p->private[1]= d->output= news(float, p->nch*d->setlen);
  return p;
}

float calcwthaarnn(sndobj *me, sndobj *caller, int innr)
{
  struct wthaarnn *d;
  int ch;
  
  d= (struct wthaarnn *)me->private[0];
  if( d->setcount >= d->setlen )
  {
    float tmp;
    int count, offset;
    
    for( count= 0; count< d->setlen; count += 2 )
    {
      BUFINCALC(0, count, 0.0);
      FORCH
        d->output[ch*d->setlen+count]= d->output[ch*d->setlen+count+1]= 
	    0.5f * INPUTC(0, ch);
      BUFINCALC(0, count+1, 0.0);
      FORCH {
        d->output[ch*d->setlen+count] -= 0.5f * INPUTC(0, ch);
	d->output[ch*d->setlen+count+1] += 0.5f * INPUTC(0, ch);
      }
    }
    for( offset= 2; offset < d->setlen; offset <<= 1 )
      for( count= offset-1; count< d->setlen; count += 2*offset )
        FORCH {
	  d->output[ch*d->setlen + count] *= 0.5f;
	  d->output[ch*d->setlen + count + offset] *= 0.5f;
	  tmp= d->output[ch*d->setlen + count + offset];
	  d->output[ch*d->setlen + count + offset] += 
			d->output[ch*d->setlen + count];
	  d->output[ch*d->setlen + count] -= tmp;
	}
    d->setcount= 0;
  }
  FORCH
    OUTPUT(ch, d->output[d->setlen*ch+d->setcount]);
  ++d->setcount;
  BUFINCR(0, 1);
  CALCRETURN;
}


/*=============================================================================
    iwthaarnn

Inverse wavelet transform corresponding to wthaarnn.  <bottom> has to be the
same minimum frequency passed to wthaarnn.  This reverse transformation expects
the wavelet data normalisation to be off, see wthaarnn.

See also: wthaarnn, wthaar, wtany
=============================================================================*/

struct iwthaarnn {
  int	depth, setlen, setcount;
  float *output;
};

float calciwthaarnn(sndobj *me, sndobj *, int);

sndobj *iwthaarnn(double bottom, sndobj *signal)
{
  sndobj *p;
  struct iwthaarnn *d;
  int depth;
  
  if( (depth= wtdepth(bottom)) == 0 )
    return signal;
  p= newsndo( calciwthaarnn, "iwthaarnn", "iwthaarnn", signal->nch, 1, signal );
  p->skip= dontskip;
  p->private[0]= d= new(struct iwthaarnn);
  d->depth= depth;
  d->setlen= 1 << d->depth;
  d->setcount= d->setlen;
  buf( 0.0, (double)d->setlen/(double)SAMPRATE, 0.0, signal->nch, p, 0 );
  p->private[1]= d->output= news(float, p->nch*d->setlen);
  return p;
}

float calciwthaarnn(sndobj *me, sndobj *caller, int innr)
{
  struct iwthaarnn *d;
  int ch;
  
  d= (struct iwthaarnn *)me->private[0];
  if( d->setcount >= d->setlen )
  {
    int count, scale;
    
    memset(d->output, 0, me->nch*d->setlen*sizeof(float));
    BUFINCALC(0, d->setlen-1, 0.0);
    FORCH
      d->output[(ch+1)*d->setlen - 1]= INPUTC(0, ch);
    for( scale= d->setlen/2; scale >= 1; scale >>= 1 )
      for( count= scale-1; count< d->setlen; count += 2*scale ) {
        BUFINCALC(0, count, 0.0);
        FORCH {
	  d->output[ch*d->setlen + count]= INPUTC(0, ch) + 
				d->output[ch*d->setlen + count + scale];
	  d->output[ch*d->setlen + count + scale] -= INPUTC(0, ch);
	}
      }
    d->setcount= 0;
  }
  FORCH
    OUTPUT(ch, d->output[d->setlen*ch+d->setcount]);
  ++d->setcount;
  BUFINCR(0, 1);
  CALCRETURN;
}


/*=============================================================================
    wtany (bottom) [scalingfunc] (llap) <signal>
    iwtany (bottom) [scalingfunc] (llap) <signal>

Wavelet transform with any wavelet with compact support which is based on
quadrature mirror filters.  This object takes as its first parameter the
requested frequency bounding the lowest wavelet band against the "remainder"
band.  The actual minimum frequency is the power-of-two quotient of the
sampling rate closest to (bottom) on the logarithmic frequency axis.  The
requested minimum frequency is limited to WTMINBOTTOM which is currently
defined as 20 Hz and can be changed in the source.

wtany takes a pointer to a float array as the second argument [scalingfunc].
The array must contain an even number of values followed by the special value
END.  These values are the coefficients of the scaling (smoothing) filter.
(llap) is the index in [scalingfunc] which is to be aligned with the current
sample.  This parameter affects the phase behaviour of the transform.

The wavelet transform is implemented as a moving filter bank, so that there are
no boundary effects between data sets.  The wavelet and scaling filters are
aligned w. r. to the current position, which for appropriate (llap) yields a
good temporal coincidence between PCM and wavelet-transformed signal.  There
remains the problem of boundary effects at the start and end of the signal.
The first value transmitted is partly represented by wavelet data which should
precede it.  To take this into account (the wavelet transform is not
orthonormal without these before-the-start data), the wavelet data are delayed
w. r. to the PCM signal by a certain number of samples.  This delay has to be
taken into account when doing something time-dependent to the
wavelet-transformed signal.  For that purpose, the object wdelay is provided,
which should _always_ be used in case the wavelet data delay is changed in a
future version of the wavelet transform implementation.

All calculations in the wavelet transform objects are done in single precision.
This is accurate to 1e-7 or so, which is quite sufficient for sound processing.

In the wavelet-transformed signals, the wavelet coefficients are organised in
data sets 2^depth in length, where depth = round(log2(SAMPRATE/bottom)).  Half
of them represent the top band, a quarter the next-highest and so on.  The
location of the data within the data set is such that every wavelet datum is as
close as possible to the centre of the interval of the source data from which
it was computed.  This is achieved with the following representation: the
first, third and so on value represents the highest frequency band or smallest
scale; the second, sixth, tenth value and so on the next larger scale; down to
the last value in the data set which represents the largest scale.  Here is a
graphical representation of a data set for depth=4:

  S S S S S S S S S S S S S S S S
  |/  |/  |/  |/  |/  |/  |/  |/  
  - + - + - + - + - + - + - + - +
    | /     | /     | /     | /  
    - +     - +     - +     - +
       \ ,----'        \ ,----'
        - +             - +
           `---. ,--------'
	        - + ----------> +

"S" are the original sample values, "-" are the "detail" coefficients, and "+"
are the smoothed values resulting from a convolution with the scaling function.
Every "+" is computed from the same original values as the "-" to its left,
only it represents the smooth part as opposed to the detail at its scale.  The
lowest "+" or "-" in each column is the wavelet coefficient output in place of,
which here means wdelay() samples after, the input sample value in the top row.
There is only one "+" or smoothed value which is output, the last value in the
set which represents the lowest frequency band, below SAMPRATE/2**depth
(approximately (bottom)).

See also: anything in sndwt.c
=============================================================================*/

sndobj *wtany_( const char *name, int ncoeffs, float *coeffs, int llap, double bottom, sndobj *signal );
sndobj *iwtany_( const char *name, int ncoeffs, float *coeffs, int llap, double bottom, sndobj *signal );


sndobj *wtany( double bottom, float *scalingfunc, int llap, sndobj *signal )
{
  int ncoeffs;

  for( ncoeffs= 0; scalingfunc[ncoeffs]!=END; ++ncoeffs);
  if( ncoeffs< 2 ) {
    fprintf(stderr, "wtany: Error: scaling function has to have at least 2 coefficients. Aborting.\n");
    exit(1);
  }
  if( (ncoeffs&1)!=0 ) {
    fprintf(stderr, "wtany: Warning: scaling function cannot have odd number of coefficients.  Discarding last.\n");
    ncoeffs &= ~1;
  }
  if( llap < 0 || llap > ncoeffs-2 ) {
    fprintf(stderr, "wtany: Warning: llap has to be in [0; #coeffs-2].\n");
    if( llap < 0 )
      llap= 0;
    if( llap > ncoeffs-2 )
      llap= ncoeffs-2;
  }
  return wtany_( "wtany", ncoeffs, scalingfunc, llap, bottom, signal );
}

sndobj *iwtany( double bottom, float *scalingfunc, int llap, sndobj *signal )
{
  int ncoeffs;

  for( ncoeffs= 0; scalingfunc[ncoeffs]!=END; ++ncoeffs);
  if( ncoeffs< 2 ) {
    fprintf(stderr, "iwtany: Error: scaling function has to have at least 2 coefficients. Aborting.\n");
    exit(1);
  }
  if( (ncoeffs&1)!=0 ) {
    fprintf(stderr, "iwtany: Warning: scaling function cannot have odd number of coefficients.  Discarding last.\n");
    ncoeffs &= ~1;
  }
  if( llap < 0 || llap > ncoeffs-2 ) {
    fprintf(stderr, "wtany: Warning: llap has to be in [0; #coeffs-2].\n");
    if( llap < 0 )
      llap= 0;
    if( llap > ncoeffs-2 )
      llap= ncoeffs-2;
  }
  return iwtany_( "iwtany", ncoeffs, scalingfunc, llap, bottom, signal );
}


/*===========================================================================*/

struct wtany {
  int	initialised, depth, setlen, setcount, ncoeffs, smchsize, maxahead;
  int	llap, rlap;
  float *smooth, *coeffs, *detcoeffs, *smoothpos;
  int	*readahead;
};

float calcwtany(sndobj *me, sndobj *, int);

// Wavelet transform with arbitrary coefficients.
// For the calculation of wavelet coefficients of all but the smallest scale,
// the result of a series of smoothing filters is needed.  The intermediate
// smoothed signals at every scale are saved in the array d->smooth for reuse.
// They are stored in the same format as the wavelet coefficients (see above,
// just substitute "+" for all "-").  Let us say D is the depth of the wavelet
// transform, so that d->setlen= 2^D. Then the number of samples we need to read
// ahead to compute one set of smoothed values is (counted from the start of a
// data set): 2^(D-1)-1 + 2^(D-2) + rlap*2^(D-1) + 2^(D-3) + rlap*2^(D-2) - ...
// + 1 + rlap*2 + 1 + rlap = 2^(D-1) - 1 + 2^(D-1) + rlap*(2^D - 1) =
// (rlap + 1)*(2^D - 1).  The largest-scale smoothed value
// has the offset 2^(D-1)-1 from the start of the set.  For each smoothed value
// at depth d+1 (scale 2^(d+1)), a level-d value is located 2^(d-1) to the
// right.  Taking into account proper centering of the convolution filter, the
// last needed value is rlap*2^d to the right of that.  Hence we get a sum over
// -2^(d-1) + rlap*2^d, for d from D-1 down to 1.  The smallest scale smoothing
// uses sample values up to rlap+1 to its right, so we add that.
// The position of the rightmost needed smoothed value at a certain scale is
// computed by truncating the sum after the relevant rlap*2^d, which gives 
// (rlap + 1)*(2^D - 2^(d+1)) + 2^d - 1.  These formulas are not quite
// everything, however.  We compute a smoothed value of a given depth at the
// time at which a wavelet datum of the same depth occurs in the input.  In
// order to have the smooth data available in time, the last of each depth has
// to be computed at the time the first wavelet datum of that depth in the data
// set is in the input.  Therefore we have to subtract 2^d - 1 from the result
// above to get the readahead relative to the current position.

sndobj *wtany_( const char *name, int ncoeffs, float *coeffs, int llap, double bottom, sndobj *signal )
{
  sndobj *p;
  struct wtany *d;
  float sign;
  int depth, count, offset;
  
  if( (depth= wtdepth(bottom)) == 0 )
    return signal;
  p= newsndo( calcwtany, name, name, signal->nch, 1, signal );
  p->skip= dontskip;
  p->private[0]= d= new(struct wtany);
  d->initialised= 0;
  d->depth= depth;
  d->setlen= 1 << d->depth;
  d->setcount= 1;
  d->ncoeffs= ncoeffs;
  d->coeffs= coeffs;
  d->llap= llap;
  d->rlap= ncoeffs - llap - 2;
  d->smchsize= d->ncoeffs*d->setlen;
  d->maxahead= (d->rlap + 1) * (d->setlen - 1);
  if( ncoeffs > 2 ) {
    p->in[0]= sdelay( (d->rlap+1) * d->setlen, signal );
    p->inch[0]= p->in[0]->ch;
  }
  buf( (double)d->llap/SAMPRATE,
       (double)d->maxahead/(double)SAMPRATE, (double)d->setlen/SAMPRATE,
           signal->nch, p, 0 );
  p->private[1]= d->smooth= news( float, p->nch*d->smchsize );
  memset( d->smooth, 0x01, sizeof(float)*p->nch*d->smchsize );
  d->smoothpos= d->smooth;
  p->private[2]= d->detcoeffs= news(float, d->ncoeffs);
  for( count= 0, sign= 1.0; count< d->ncoeffs; ++count, sign= -sign ) {
    d->detcoeffs[count]= sign * d->coeffs[d->ncoeffs-1-count];
  }
  p->private[3]= d->readahead= news(int, d->depth);
  for(count= 0, offset= 2; count< d->depth; ++count, offset<<=1 )
    d->readahead[count]= (d->rlap + 1) * (d->setlen - offset);
  return p;
}

float calcwtany(sndobj *me, sndobj *caller, int innr)
{
  struct wtany *d;
  float *write, *read;
  int sample, c_index, ch, offset;

  d= (struct wtany *)me->private[0];
  if( (d->setcount&1) != 0 )
  {
    FORCH
      OUTPUT(ch, 0.0);
    for( sample= 0; sample < d->ncoeffs; ++sample ) {
      BUFINCALC(0, sample - d->llap, 0.0);
      FORCH
        me->ch[ch] += d->detcoeffs[sample] * INPUTC(0, ch);
    }
    if( d->initialised ) 
    {
      write= d->smoothpos + (d->rlap+1)*(d->setlen-2);
      if( write >= d->smooth + d->smchsize )
	write -= d->smchsize;
      FORCH
	write[ch*d->smchsize] = 0.0;
      sample= d->maxahead - d->ncoeffs + 1;
      for( c_index= 0; c_index< d->ncoeffs; ++c_index, ++sample ) {
	BUFINCALC(0, sample, 0.0);
	FORCH
	  write[ch*d->smchsize] += INPUTC(0, ch)*d->coeffs[c_index];
      }
    }
    else
    {
      // Initialisation: compute all smoothed values we need ahead of time;
      // later we will add new values as we go along.
      d->initialised= 1;
      // Smallest-scale smoothing.  We do it sample-wise, not result-value-wise,
      // to save some overhead in BUFINCALC.
      for( sample= 0; sample < (d->rlap+2)*d->setlen/2; ++sample )
      {
	BUFINCALC(0, sample, 0.0);
	write= d->smooth + ((sample - d->llap)&~1);
	if( write < d->smooth )
	  write += d->smchsize;
	else if( write> d->smooth+d->smchsize )
	  write -= d->smchsize;
	//  The complicated bit manipulation in the initialisation of c_index
	//  ensures that samples alternately are multiplied with only even and
	//  only odd-indexed wavelet coefficients.  These are just every second
	//  entry in the columns of the transformation matrix, those entries
	//  which give the smoothed results.
	for( c_index= d->ncoeffs-1-(((d->rlap+1)&1)^(sample&1)); c_index > 0; c_index -= 2 ) {
	  if( write > d->smoothpos + d->maxahead - d->rlap - 1 )
	    break;
	  FORCH
	    write[ch*d->smchsize] += d->coeffs[c_index] * INPUTC(0, ch);
	  write += 2;
	}
	if( c_index == 0 )
	  FORCH
	    write[ch*d->smchsize]= d->coeffs[c_index] * INPUTC(0, ch);
      }
      // Smoothing at larger scales, from the bottom up.  offset is the scale
      // of the smoothing which is done already, half the scale of the
      // smoothing which is currently computed.  With the conventions of above
      // formulas, offset = 2^(d-1).
      for( offset= 2; offset< d->setlen; offset <<= 1 )
	for( write= d->smooth + offset - 1; write < d->smooth + 
	    (d->rlap+1)*(d->setlen-2*offset) + offset - 1; write += 2*offset )
	{
	  FORCH
	    write[ch*d->smchsize]= 0.0;
	  read= write - offset/2 + d->ncoeffs/2*offset;
	  for( c_index= d->ncoeffs-1; c_index >= 0; --c_index, read -= offset )
	  {
	    if( read< d->smoothpos )
	      read += d->smchsize;
	    FORCH
	      write[ch*d->smchsize] += d->coeffs[c_index] *
				       read[ch*d->smchsize];
	  }
	}
    }
  }
  else if( d->setcount == d->setlen )
  {
    // Return the one smooth wavelet coefficient corresponding to the largest
    // scale:
    FORCH
      OUTPUT(ch, d->smoothpos[ch*d->smchsize - d->setlen/2]);
    d->setcount= 0;
    if( d->smoothpos >= d->smooth + d->smchsize - 1 )
      d->smoothpos= d->smooth - 1;
  }
  else
  {
    offset= 1 << pow2div(d->setcount); 
    //  First compute one more smoothed value at our scale:
    write= d->smoothpos + (d->rlap+1)*(d->setlen-2*offset);
    if( write >= d->smooth + d->smchsize )
      write -= d->smchsize;
    FORCH {
      *write = 0.0;
      read= write - offset/2 - d->llap*offset;
      if( read < d->smooth + ch*d->smchsize )
	read += d->smchsize;
      for( c_index= 0; c_index< d->ncoeffs; ++c_index ) {
	*write += d->coeffs[c_index] * *read;
	read += offset;
	if( read >= d->smooth + (ch+1)*d->smchsize )
	  read -= d->smchsize;
      }
      write += d->smchsize;
    }
    //  Now compute and output wavelet coefficient:
    read= d->smoothpos - offset/2 - d->llap*offset;
    FORCH {
      if( read < d->smooth + ch*d->smchsize )
	read += d->smchsize;
      me->ch[ch] = 0.0;
      for( c_index= 0; c_index< d->ncoeffs; ++c_index ) {
	me->ch[ch] += d->detcoeffs[c_index] * *read;
	read += offset;
	if( read >= d->smooth + (ch+1)*d->smchsize )
	  read -= d->smchsize;
      }
      read += d->smchsize - d->ncoeffs*offset;
    }
  }
  BUFINCR(0, 1);
  ++d->setcount;
  ++d->smoothpos;
  CALCRETURN;
}


/*===========================================================================*/

struct iwtany {
  int	initialised, depth, setlen, setcount, ncoeffs, ncoeffs2, smchsize;
  int	llap, rlap;
  float *smooth, **smch, **smpos, *coeffs, *detcoeffs;
  int	*readahead, *rwoff;
};

float calciwtany(sndobj *me, sndobj *, int);

sndobj *iwtany_( const char *name, int ncoeffs, float *coeffs, int llap, double bottom, sndobj *signal )
{
  sndobj *p;
  struct iwtany *d;
  float sign;
  int depth, count, offset;
  
  if( (depth= wtdepth(bottom)) == 0 )
    return signal;
  p= newsndo( calciwtany, name, name, signal->nch, 1, signal );
  p->skip= dontskip;
  p->private[0]= d= new(struct iwtany);
  d->initialised= 0;
  d->depth= depth;
  d->setlen= 1 << d->depth;
  d->setcount= 0;
  d->ncoeffs= ncoeffs;
  d->ncoeffs2= ncoeffs/2;
  d->coeffs= coeffs;
  d->llap= llap;
  d->rlap= ncoeffs - llap - 2;
  d->smchsize= 2*d->ncoeffs + 4;
  p->private[1]= d->smooth= news( float, p->nch*d->smchsize*d->depth );
  // clear smooth data storage - this is expected in initialisation of calc
  memset( d->smooth, 0, sizeof(float)*p->nch*d->smchsize*d->depth );
  p->private[2]= d->smch= news(float*, 2*d->depth);
  d->smpos= d->smch + d->depth;
  for( count= 0; count< d->depth; ++count )
    d->smpos[count]= d->smch[count]= d->smooth + count*p->nch*d->smchsize;
  p->private[3]= d->detcoeffs= news(float, d->ncoeffs);
  for( count= 0, sign= 1.0; count< d->ncoeffs; ++count, sign= -sign )
    d->detcoeffs[count]= sign * d->coeffs[d->ncoeffs-1-count];
  p->private[4]= d->readahead= news(int, 2*d->depth);
  d->rwoff= d->readahead + d->depth;
  // readahead[i] is the sample to be read at the time a value of depth i would 
  // be read from input if it were not for readahead.  It is an offset relative
  // to the wavelet datum _following_ the current depth-i value.  How far this
  // extends into the future is determined by the requirement that the
  // calculation should complete all values needed by depth i-1 (or by the
  // output in the case of i=0) until the next depth-i value occurs.
  d->readahead[0]= d->llap | 1;
  d->rwoff[0]= (d->llap + 1) & 1;
//  printf("depth 0: readahead %d, rwoff %d\n", d->readahead[0], d->rwoff[0]);
  for( count= 1, offset= 2; count< d->depth; ++count, offset <<= 1 ) { 
    // The further readaheads are computed with a tough recursion formula.
    // Note that the spread of depth-i values is 2^(i+1)=2*offset.  The -1 at
    // the beginning makes the result relative to the following input value as
    // described above and facilitates the recursion.  Adding 3*offset/2 gets us
    // the next but one depth-(i-1) value, which is the last one before the
    // next depth-i value.  Adding the depth-(i-1) readahead gets us the
    // rightmost depth-i value which is needed for that one.  offset*(llap+1)
    // is the depth-i readahead for computing the corresponding smoothed value.
    // The logical & at the end serves to align the result.
    d->readahead[count]= -1 + ((3*(offset/2) + d->readahead[count-1]
	  			+ offset*(d->llap+1)) & ~(2*offset-1));
    // Here we compute additionally the readahead needed for only the next
    // depth-(i-1) value.  If it differs from the previously computed
    // readahead, we have to leave a 1-value gap between the reading and
    // writing position at this depth.  There must be a closed formula for
    // this, but it's hard to find.
    if( d->readahead[count] == -1 + ((offset/2 + d->readahead[count-1]
	    			+ offset*(d->llap+1)) & ~(2*offset-1)) )
      d->rwoff[count]= 0;
    else
      d->rwoff[count]= 1;
//    printf("depth %d: readahead %d, rwoff %d\n", count, d->readahead[count], d->rwoff[count]);
  }
  buf( 0.0, (double)(d->readahead[d->depth-1]+d->setlen/2)/(double)SAMPRATE,
      (double)d->setlen/SAMPRATE, signal->nch, p, 0 );
  if( ncoeffs > 2 ) {
    return sdelay( -(d->rlap+1)*d->setlen, p );
  }
  else
    return p;
}

//  The calculation of the reverse wavelet transform is implemented in a quite
//  different way from the transform to wavelets.  The inverse transform
//  necessarily works from the top level (coarsest in time resolution)
//  downward.  To keep storage requirements to a minimum, every fully computed
//  smoothed value at a given level is added (with the right prefactors) to
//  every value one level down to which it contributes and can then be thrown
//  away.  This method is similar to the smallest-scale smoothing in the
//  initialisation of the transform to wavelets in calcwtany().  As values are
//  used as soon as their calculation is completed, only d->ncoeffs values at
//  each level have to be stored while terms are added to them.  Actually we
//  round that amount up to the next multiple of 4 to make the algorithm
//  slightly simpler.

float calciwtany(sndobj *me, sndobj *caller, int innr)
{
  struct iwtany *d;
  float *write, *coeffs;
  int sample, c_index, ch, offset, depth;

  d= (struct iwtany *)me->private[0];
  if( !d->initialised )
  {
    d->initialised= 1;
    // In this first part of the initilisation we compute the smoothed data
    // at the largest but one scale.  This is done by reading the
    // largest-scale smooth and detail data and adding their product with the
    // coefficients to the corresponding elements of d->smooth.
    for( sample= d->setlen/2-1; sample < d->setlen/2 + d->readahead[d->depth-1]; sample += d->setlen/2 )
    {
      c_index= d->llap - ((sample >> (d->depth-1)) & ~1);
      write= d->smch[d->depth-1];
      if( c_index < 0 ) {
	write -= c_index;
	c_index= 0;
      }
      BUFINCALC(0, sample, 0.0);
      if( (sample&(d->setlen/2))==0 )
	coeffs= d->detcoeffs;
      else
	coeffs= d->coeffs;
      for(; c_index< d->ncoeffs; ++c_index, ++write ) {
	FORCH
	  write[ch*d->smchsize] += coeffs[c_index]*INPUTC(0, ch);
      }
    }
    // Now successively compute the smaller-scale smoothed data
    for( depth= d->depth-2, offset= d->setlen/2; depth >= 0; --depth, offset /= 2 ) {
      for( sample= offset/2-1; sample < offset/2 + d->readahead[depth]; sample += offset )
      {
	c_index= d->llap - ((sample >> depth) & ~1);
	write= d->smch[depth];
	if( c_index < 0 ) {
	  write -= c_index;
	  c_index= 0;
	}
	BUFINCALC(0, sample, 0.0);
	for( ; c_index< d->ncoeffs; ++c_index, ++write ) {
	  FORCH
	    write[ch*d->smchsize] += d->detcoeffs[c_index]*INPUTC(0, ch) +
		  d->coeffs[c_index]*d->smpos[depth+1][ch*d->smchsize];
	}
	if( ++d->smpos[depth+1] >= d->smch[depth+1] + d->smchsize )
	  d->smpos[depth+1]= d->smch[depth+1];
      }
    }
  }

  FORCH
    OUTPUT(ch, d->smpos[0][ch*d->smchsize]);
  if( ++d->smpos[0] >= d->smch[0]+d->smchsize )
    d->smpos[0]= d->smch[0];
  ++d->setcount;
  BUFINCR(0, 1);

  if( d->setcount == d->setlen )
    d->setcount= 0;
  else if( d->setcount == d->setlen/2 )
  {
    sample= d->readahead[d->depth-1];
    BUFINCALC(0, sample, 0.0);
    write= d->smpos[d->depth-1];
    if( d->rwoff[d->depth-1] && ++write >= d->smch[d->depth-1]+d->smchsize )
      write -= d->smchsize;
    for( c_index= 0; c_index< d->ncoeffs-2; ++c_index ) {
      FORCH
	write[ch*d->smchsize] += d->detcoeffs[c_index]*INPUTC(0, ch);
      if( ++write >= d->smch[d->depth-1] + d->smchsize )
	write= d->smch[d->depth-1];
    }
    FORCH
      write[ch*d->smchsize]= d->detcoeffs[c_index]*INPUTC(0, ch);
    ++c_index;
    if( ++write >= d->smch[d->depth-1] + d->smchsize )
      write= d->smch[d->depth-1];
    FORCH
      write[ch*d->smchsize]= d->detcoeffs[c_index]*INPUTC(0, ch);
    BUFINCALC(0, sample+d->setlen/2, 0.0);
    write= d->smpos[d->depth-1];
    if( d->rwoff[d->depth-1] && ++write >= d->smch[d->depth-1]+d->smchsize )
      write -= d->smchsize;
    for( c_index= 0; c_index< d->ncoeffs; ++c_index ) {
      FORCH
	write[ch*d->smchsize] += d->coeffs[c_index]*INPUTC(0, ch);
      if( ++write >= d->smch[d->depth-1] + d->smchsize )
	write= d->smch[d->depth-1];
    }
  }
  else
  {
    depth= pow2div(d->setcount);
    sample= d->readahead[depth];
    BUFINCALC(0, sample, 0.0);
    write= d->smpos[depth];
    if( d->rwoff[depth] && ++write >= d->smch[depth] + d->smchsize )
      write -= d->smchsize;
    for( c_index= 0; c_index< d->ncoeffs-2; ++c_index ) {
      FORCH
	write[ch*d->smchsize] += d->detcoeffs[c_index]*INPUTC(0, ch) +
		       d->coeffs[c_index]*d->smpos[depth+1][ch*d->smchsize];
      if( ++write >= d->smch[depth] + d->smchsize )
	write= d->smch[depth];
    }
    FORCH
      write[ch*d->smchsize]= d->detcoeffs[c_index]*INPUTC(0, ch) +
		     d->coeffs[c_index]*d->smpos[depth+1][ch*d->smchsize];
    ++c_index;
    if( ++write >= d->smch[depth] + d->smchsize )
      write= d->smch[depth];
    FORCH
      write[ch*d->smchsize]= d->detcoeffs[c_index]*INPUTC(0, ch) +
		     d->coeffs[c_index]*d->smpos[depth+1][ch*d->smchsize];
    if( ++d->smpos[depth+1] >= d->smch[depth+1] + d->smchsize )
      d->smpos[depth+1]= d->smch[depth+1];
  }

  CALCRETURN;
}



/*=============================================================================
    wthaar (bottom) <signal>
    iwthaar (bottom) <signal>

Haar wavelet transform, properly normalised.  Computed with the general-purpose
wavelet transform wtany.  As a consequence, this is slower than the
non-normalised (i)wthaarnn.  (bottom) is the boundary frequency delimiting the
lowest frequency band, and <signal> the signal to be transformed.  See (i)wtany
for details concerning the transformation.

See also: wthaarnn, iwthaarnn, wtany
=============================================================================*/

typedef struct {
  float *scaling;
  char  *name;
}
wttype;

float haarscaling[]= { 0.70710678118654752440, 0.70710678118654752440 };

sndobj *wthaar(double bottom, sndobj *signal)
{
  return wtany_("wthaar", 2, haarscaling, 0, bottom, signal);
}

sndobj *iwthaar(double bottom, sndobj *signal)
{
  return iwtany_("iwthaar", 2, haarscaling, 0, bottom, signal);
}


/*=============================================================================
    wtdaub (n) (bottom) <signal>
    iwtdaub (n) (bottom) <signal>

Daubechies wavelet transformation.  The order (n) of the transform has to be
even and between 2 and 20.  For n=2, this transform is identical to the
(normalised) Haar wavelet transform.  The scaling and wavelet filters are
aligned left, so this is not a constant-phase transform.  (bottom) is the
boundary frequency delimiting the lowest frequency band, and <signal> the
signal to be transformed.  See (i)wtany for details concerning the
transformation.

See also: wtdaubla, wtdaubbl, wtany
=============================================================================*/

static float daub4scaling[]= { 0.4829629131445341, 0.8365163037378079, 
			    0.2241438680420134, -0.1294095225512604 };

static float daub6scaling[]= { 0.3326705529500826, 0.8068915093110926,
  	0.4598775021184916, -0.1350110200102546, -0.0854412738820267,
	0.0352262918857095 };

static float daub8scaling[]= {
   0.2303778133074431, 0.7148465705484058, 0.6308807679358788,
   -0.0279837694166834, -0.1870348117179132, 0.0308413818353661,
   0.0328830116666778, -0.0105974017850021 };

static float daub10scaling[]= {
   0.1601023979741930, 0.6038292697971898, 0.7243085284377729,
   0.1384281459013204, -0.2422948870663824, -0.0322448695846381,
   0.0775714938400459, -0.0062414902127983, -0.0125807519990820,
   0.0033357252854738 };

static float daub12scaling[]= {
  0.1115407433501094, 0.4946238903984530, 0.7511339080210954,
  0.3152503517091980, -0.2262646939654399, -0.1297668675672624,
  0.0975016055873224, 0.0275228655303053, -0.0315820393174862,
  0.0005538422011614, 0.0047772575109455, -0.0010773010853085 };

static float daub14scaling[]= {
  0.0778520540850081, 0.3965393194819136, 0.7291320908462368,
  0.4697822874052154, -0.1439060039285293, -0.2240361849938538,
  0.0713092192668312, 0.0806126091510820, -0.0380299369350125,
  -0.0165745416306664, 0.0125509985560993, 0.0004295779729214,
  -0.0018016407040474, 0.0003537137999745 };

static float daub16scaling[]= {
  0.0544158422431049, 0.3128715909143031, 0.6756307362972904,
  0.5853546836541907, -0.0158291052563816, -0.2840155429615702,
  0.0004724845739124, 0.1287474266204837, -0.0173693010018083,
  -0.0440882539307952, 0.0139810279173995, 0.0087460940474061,
  -0.0048703529934518, -0.0003917403733770, 0.0006754494064506,
  -0.0001174767841248 };

static float daub18scaling[]= {
  0.0380779473638791, 0.2438346746125939, 0.6048231236901156,
  0.6572880780512955, 0.1331973858249927, -0.2932737832791761,
  -0.0968407832229524, 0.1485407493381306, 0.0307256814793395,
  -0.0676328290613302, 0.0002509471148340, 0.0223616621236805,
  -0.0047232047577520, -0.0042815036824636, 0.0018476468830564,
  0.0002303857635232, -0.0002519631889427, 0.0000393473203163 };

static float daub20scaling[]= {
  0.0266700579005546, 0.1881768000776863, 0.5272011889317202,
  0.6884590394536250, 0.2811723436606485, -0.2498464243272283,
  -0.1959462743773399, 0.1273693403357890, 0.0930573646035802,
  -0.0713941471663697, -0.0294575368218480, 0.0332126740593703,
  0.0036065535669880, -0.0107331754833036, 0.0013953517470692,
  0.0019924052951930, -0.0006858566949566, -0.0001164668551285,
  0.0000935886703202, -0.0000132642028945 };

wttype daub[]= {
  haarscaling, "iwthaar", daub4scaling, "iwtdaub4", daub6scaling, "iwtdaub6",
  daub8scaling, "iwtdaub8", daub10scaling, "iwtdaub10", daub12scaling,
  "iwtdaub12", daub14scaling, "iwtdaub14", daub16scaling, "iwtdaub16",
  daub18scaling, "iwtdaub18", daub20scaling, "iwtdaub20" };

#define DAUBMAX		20

int wtdaubcheckn( int n, const char *name );

int wtdaubcheckn( int n, const char *name )
{
  const char *ntoolarge= "%s: Warning: Wavelet filter order given (%d) is larger than the maximum available (%d).  Substituting maximum.\n";
  const char *nodd= "%s: Warning: Wavelet filter order has to be even.  Replacing by next largest even number.\n";
  const char *n0= "%s: Warning: You gave wavelet filter order <= 0.  You kidding?  Substituting n=2, Haar wavelets.\n";

  if( n > DAUBMAX ) {
    fprintf(stderr, ntoolarge, name, n, DAUBMAX);
    n= DAUBMAX;
  }
  if( n <= 0 ) {
    fprintf(stderr, n0, name);
    n= 2;
  }
  if( 0 != (n&1) ) {
    fprintf(stderr, nodd, name);
    ++n;
  }
  return n;
}

sndobj *wtdaub( int n, double bottom, sndobj *signal )
{
  n= wtdaubcheckn(n, "wtdaub");
  return wtany_(daub[n/2-1].name+1, n, daub[n/2-1].scaling, 0, bottom, signal);
}

sndobj *iwtdaub( int n, double bottom, sndobj *signal )
{
  n= wtdaubcheckn(n, "iwtdaub");
  return iwtany_(daub[n/2-1].name, n, daub[n/2-1].scaling, 0, bottom, signal);
}


/*=============================================================================
    wtdaubla (n) (bottom) <signal>
    iwtdaubla (n) (bottom) <signal>

Least asymmetric Daubechies wavelet transformation.  These transforms are a
good compromise between good phase behaviour, separation of frequency bands and
lack of artifacts.  The order of the coefficients has been reversed as
recommended by Percival and Walden.  The order (n) of the transform has to be
even and between 2 and 20.  For n<=6, this transform has coefficients identical
to the normal Daubechies wavelet transform.  The scaling and wavelet filters
are centred, which makes this transform close to zero-phase.  (bottom) is the
boundary frequency delimiting the lowest frequency band, and <signal> the
signal to be transformed.  See (i)wtany for details concerning the
transformation.

See also: wtdaub, wtdaubbl, wtany
=============================================================================*/

static float daubla8scaling[]= {
  -0.0757657147893407, -0.0296355276459541, 0.4976186676324578,
  0.8037387518052163, 0.2978577956055422, -0.0992195435769354,
  -0.0126039672622612, 0.0322231006040713 };

static float daubla10scaling[]= {
  0.0195388827353869, -0.0211018340249298, -0.1753280899081075,
  0.0166021057644243, 0.6339789634569490, 0.7234076904038076,
  0.1993975339769955, -0.0391342493025834, 0.0295194909260734,
  0.0273330683451645 };

static float daubla12scaling[]= {
  0.0154041093273377, 0.0034907120843304, -0.1179901111484105,
  -0.0483117425859981, 0.4910559419276396, 0.7876411410287941,
  0.3379294217282401, -0.0726375227866000, -0.0210602925126954,
  0.0447249017707482, 0.0017677118643983, -0.0078007083247650 };

static float daubla14scaling[]= {
  0.0102681767084968, 0.0040102448717033, -0.1078082377036168,
  -0.1400472404427030, 0.2886296317509833, 0.7677643170045710,
  0.5361019170907720, 0.0174412550871099, -0.0495528349370410,
  0.0678926935015971, 0.0305155131659062, -0.0126363034031526,
  -0.0010473848889657, 0.0026818145681164 };

static float daubla16scaling[]= {
  -0.0033824159513594, -0.0005421323316355, 0.0316950878103452,
  0.0076074873252848, -0.1432942383510542, -0.0612733590679088,
  0.4813596512592012, 0.7771857516997478, 0.3644418948359564,
  -0.0519458381078751, -0.0272190299168137, 0.0491371796734768,
  0.0038087520140601, -0.0149522583367926, -0.0003029205145516,
  0.0018899503329007 };

static float daubla18scaling[]= {
  0.0010694900326538, -0.0004731544985879, -0.0102640640276849,
  0.0088592674935117, 0.0620777893027638, -0.0182337707798257,
  -0.1915508312964873, 0.0352724880359345, 0.6173384491413523,
  0.7178970827642257, 0.2387609146074182, -0.0545689584305765,
  0.0005834627463312, 0.0302248788579895, -0.0115282102079848,
  -0.0132719677815332, 0.0006197808890549, 0.0014009155255716 };

static float daubla20scaling[]= {
  0.0007701598091030, 0.0000956326707837, -0.0086412992759401,
  -0.0014653825833465, 0.0459272392237649, 0.0116098939129724,
  -0.1594942788575307, -0.0708805358108615, 0.4716906668426588,
  0.7695100370143388, 0.3838267612253823, -0.0355367403054689,
  -0.0319900568281631, 0.0499949720791560, 0.0057649120455518,
  -0.0203549398039460, -0.0008043589345370, 0.0045931735836703,
  0.0000570360843390, -0.0004593294205481 };

wttype daubla[]= {
  haarscaling, "iwthaar", daub4scaling, "iwtdaub4", daub6scaling, "iwtdaub6",
  daubla8scaling, "iwtdaubla8", daubla10scaling, "iwtdaubla10",
  daubla12scaling, "iwtdaubla12", daubla14scaling, "iwtdaubla14",
  daubla16scaling, "iwtdaubla16", daubla18scaling, "iwtdaubla18",
  daubla20scaling, "iwtdaubla20" };

sndobj *wtdaubla( int n, double bottom, sndobj *signal )
{
  n= wtdaubcheckn(n, "wtdaubla");
  return wtany_(daubla[n/2-1].name+1, n, daubla[n/2-1].scaling, n/2-1, bottom,
      			signal);
}

sndobj *iwtdaubla( int n, double bottom, sndobj *signal )
{
  n= wtdaubcheckn(n, "iwtdaubla");
  return iwtany_(daubla[n/2-1].name, n, daubla[n/2-1].scaling, n/2-1, bottom,
      			signal);
}


/*=============================================================================
    wtdaubbl (n) (bottom) <signal>
    iwtdaubbl (n) (bottom) <signal>

Best localised Daubechies wavelet transformation.  These transforms are
optimised with regard to linearity of phase, which makes them close to zero
phase when shifted appropriately, but may cause artifacts if the wavelet
data are manipulated.  The order of the coefficients has been reversed
as recommended by Percival and Walden.  The order (n) of the transform has to
be even and between 2 and 20.  For n<=12 and n=16, this transform is identical
to the least asymmetric Daubechies wavelet transform.  The scaling and wavelet
filters are centred, which is rather close to the optimum phase behaviour.
(bottom) is the boundary frequency delimiting the lowest frequency band, and
<signal> the signal to be transformed.  See (i)wtany for details concerning the
transformation.

See also: wtdaub, wtdaubla, wtany
=============================================================================*/

float daubbl14scaling[]= {
  0.0120154192834842, 0.0172133762994439, -0.0649080035533744,
  -0.0641312898189170, 0.3602184608985549, 0.7819215932965554,
  0.4836109156937821, -0.0568044768822707, -0.1010109208664125,
  0.0447423494687405, 0.0204642075778225, -0.0181266051311065,
  -0.0032832978473081, 0.0022918339541009 };

float daubbl18scaling[]= {
  0.0002594576266544, -0.0006273974067728, -0.0019161070047557,
  0.0059845525181721, 0.0040676562965785, -0.0295361433733604,
  -0.0002189514157348, 0.0856124017265279, -0.0211480310688774,
  -0.1432929759396520, 0.2337782900224977, 0.7374707619933686,
  0.5926551374433956, 0.0805670008868546, -0.1143343069619310,
  -0.0348460237698368, 0.0139636362487191, 0.0057746045512475 };

float daubbl20scaling[]= {
  0.0008625782242896, 0.0007154205305517, -0.0070567640909701,
  0.0005956827305406, 0.0496861265075979, 0.0262403647054251,
  -0.1215521061578162, -0.0150192395413644, 0.5137098728334054,
  0.7669548365010849, 0.3402160135110789, -0.0878787107378667,
  -0.0670899071680668, 0.0338423550064691, -0.0008687519578684,
  -0.0230054612862905, -0.0011404297773324, 0.0050716491945793,
  0.0003401492622332, -0.0004101159165852 };

wttype daubbl[]= {
  haarscaling, "iwthaar", daub4scaling, "iwtdaub4", daub6scaling, "iwtdaub6",
  daubla8scaling, "iwtdaubla8", daubla10scaling, "iwtdaubla10",
  daubla12scaling, "iwtdaubla12", daubbl14scaling, "iwtdaubbl14",
  daubla16scaling, "iwtdaubla16", daubbl18scaling, "iwtdaubbl18",
  daubbl20scaling, "iwtdaubbl20" };

sndobj *wtdaubbl( int n, double bottom, sndobj *signal )
{
  n= wtdaubcheckn(n, "wtdaubbl");
  return wtany_(daubbl[n/2-1].name+1, n, daubbl[n/2-1].scaling, n/2-1, bottom,
      			signal);
}

sndobj *iwtdaubbl( int n, double bottom, sndobj *signal )
{
  n= wtdaubcheckn(n, "iwtdaubbl");
  return iwtany_(daubbl[n/2-1].name, n, daubbl[n/2-1].scaling, n/2-1, bottom,
      			signal);
}


/*=============================================================================
    wtcoif (n) (bottom) <signal>
    iwtcoif (n) (bottom) <signal>

Coifman wavelet transformation.  These wavelets are derived by requiring a
zero-moment condition for the wavelet as well as scaling filters.  They tend to
create artifacts if wavelet data are manipulated.  The order (n) of the
transform has to be a multiple of 6 and between 6 and 30.  The scaling and
wavelet filters are centred.  The order of the coefficients has been reversed
as recommended by Percival and Walden.  (bottom) is the boundary frequency
delimiting the lowest frequency band, and <signal> the signal to be
transformed.  See (i)wtany for details concerning the transformation.

See also: wtdaub, wtany
=============================================================================*/

float coif6scaling[]= {
  -0.0156557285289848, -0.0727326213410511, 0.3848648565381134,
  0.8525720416423900, 0.3378976709511590, -0.0727322757411889 };

float coif12scaling[]= {
  -0.0007205494453679, -0.0018232088707116, 0.0056114348194211,
  0.0236801719464464, -0.0594344186467388, -0.0764885990786692,
  0.4170051844236707, 0.8127236354493977, 0.3861100668229939,
  -0.0673725547222826, -0.0414649367819558, 0.0163873364635998 };

float coif18scaling[]= {
  -0.0000345997728362, -0.0000709833031381, 0.0004662169601129,
  0.0011175187708906, -0.0025745176887502, -0.0090079761366615,
  0.0158805448636158, 0.0345550275730615, -0.0823019271068856,
  -0.0717998216193117, 0.4284834763776168, 0.7937772226256169,
  0.4051769024096150, -0.0611233900026726, -0.0657719112818552,
  0.0234526961418362, 0.0077825964273254, -0.0037935128644910 };

float coif24scaling[]= {
  -0.0000017849850031, -0.0000032596802369, 0.0000312298758654,
  0.0000623390344610, -0.0002599745524878, -0.0005890207562444,
  0.0012665619292991, 0.0037514361572790, -0.0056582866866115,
  -0.0152117315279485, 0.0250822618448678, 0.0393344271233433,
  -0.0962204420340021, -0.0666274742634348, 0.4343860564915321,
  0.7822389309206135, 0.4153084070304910, -0.0560773133167630,
  -0.0812666996808907, 0.0266823001560570, 0.0160689439647787,
  -0.0073461663276432, -0.0016294920126020, 0.0008923136685824 };

float coif30scaling[]= {
  -0.0000000951765727, -0.0000001674428858, 0.0000020637618516,
  0.0000037346551755, -0.0000213150268122, -0.0000413404322768,
  0.0001405411497166, 0.0003022595818445, -0.0006381313431115,
  -0.0016628637021860, 0.0024333732129107, 0.0067641854487565,
  -0.0091642311634348, -0.0197617789446276, 0.0326835742705106,
  0.0412892087544753, -0.1055742087143175, -0.0620359639693546,
  0.4379916262173834, 0.7742896037334738, 0.4215662067346898,
  -0.0520431631816557, -0.0919200105692549, 0.0281680289738655,
  0.0234081567882734, -0.0101311175209033, -0.0041593587818186,
  0.0021782363583355, 0.0003585896879330, -0.0002120808398259 };

wttype coif[]= {
  coif6scaling, "iwtcoif6", coif12scaling, "iwtcoif12", coif18scaling,
  "iwtcoif18", coif24scaling, "iwtcoif24", coif30scaling, "iwtcoif30" };

#define COIFMAX		30

int wtcoifcheckn( int n, const char *name );

int wtcoifcheckn( int n, const char *name )
{
  const char *ntoolarge= "%s: Warning: Coiflet filter order given (%d) is larger than the maximum available (%d).  Substituting maximum.\n";
  const char *n6= "%s: Warning: Coiflet filter order has to be multiple of six.  Replacing by next largest multiple of six.\n";
  const char *n0= "%s: Warning: You gave coiflet filter order <= 0.  You kidding?  Substituting n=6.\n";

  if( n > COIFMAX ) {
    fprintf(stderr, ntoolarge, name, n, COIFMAX);
    n= COIFMAX;
  }
  if( n <= 0 ) {
    fprintf(stderr, n0, name);
    n= 6;
  }
  if( 6*(n/6) != n ) {
    fprintf(stderr, n6, name);
    n= 6*((n+5)/6);
  }
  return n;
}

sndobj *wtcoif( int n, double bottom, sndobj *signal )
{
  n= wtcoifcheckn(n, "wtcoif");
  return wtany_(coif[n/6-1].name+1, n, coif[n/6-1].scaling, n/2-1, bottom,
      			signal);
}

sndobj *iwtcoif( int n, double bottom, sndobj *signal )
{
  n= wtcoifcheckn(n, "iwtcoif");
  return iwtany_(coif[n/6-1].name, n, coif[n/6-1].scaling, n/2-1, bottom,
      			signal);
}



/*=============================================================================
    wdelay (bottom) (n) (llap) <signal>

Delay of wavelet-transformed signal relative to original PCM signal for the
given maximum scale and number of wavelet filter coefficients.  This delay
should always be used with time-dependent signals controlling operations on
wavelet signals.  When a signal is synthesised with different wavelets than
it was analysed, this delay should be used first (upstream), as it is an
actual positive delay.

(bottom) is the bottom requency of the lowest wavelet band, (n) is the number
of coefficients of the wavelet scaling function, (llap) is the number of
positions the wavelet and scaling filters are shifted to the right, and (signal)
the input signal.  (llap) is zero for "ordinary" Daubechies wavelet transforms
and Haar wavelets; for the others it is =n/2-1, which is automatically computed
when the constant WDELAY_CENTRE is passed.

See also: wtany, iwdelay
=============================================================================*/

sndobj *wdelay( double bottom, int n, int llap, sndobj *signal )
{
  int depth, setlen;

  if( (n&1)!=0 ) {
    fprintf(stderr, "wdelay: Warning: scaling function cannot have odd number of coefficients.  Discarding last.\n");
    n &= ~1;
  }
  if( llap == WDELAY_CENTRE )
    llap= n/2-1;
  depth= wtdepth(bottom);
  if( depth==0 || n<=2 )
    return signal;
  setlen= 1 << depth;
  return sdelay( (n-llap-1)*setlen, signal );
}


/*=============================================================================
    iwdelay (bottom) (n) (llap) <signal>

Delay of original PCM signal relative to wavelet-transformed signal for the
given maximum scale and number of wavelet filter coefficients.  This delay
is negative.  It is useful when using a wavelet-transformed signal to control
a different PCM signal.  When used to translate between wavelet signals from
different wavelet transforms, wtdelay should be used first, upstream of this
delay, which is negative and will discard initial data.

(bottom) is the bottom requency of the lowest wavelet band, (n) is the number
of coefficients of the wavelet scaling function, (llap) is the number of
positions the wavelet and scaling filters are shifted to the right, and <signal>
the input signal.  (llap) is zero for "ordinary" Daubechies wavelet transforms
and Haar wavelets; for the others it is =n/2-1, which is automatically computed
when the constant WDELAY_CENTRE is passed.

See also: wtany, wdelay
=============================================================================*/

sndobj *iwdelay( double bottom, int n, int llap, sndobj *signal )
{
  int depth, setlen;

  if( (n&1)!=0 ) {
    fprintf(stderr, "iwdelay: Warning: scaling function cannot have odd number of coefficients.  Discarding last.\n");
    n &= ~1;
  }
  if( llap == WDELAY_CENTRE )
    llap= n/2-1;
  depth= wtdepth(bottom);
  if( depth==0 || n<=2 )
    return signal;
  setlen= 1 << depth;
  return sdelay( -(n-llap-1)*setlen, signal );
}


/*=============================================================================
    wprint *stream* (bottom) <signal>

Print out wavelet data as text.  The values of different scales (frequency
bands) are printed in different columns, separated by pairs of tab characters.
The first column (the scale 2/(sampling rate)) is always present, the following
columns successively rarer, since every larger scale has a time resolution
worse by a factor of 2.  Values of different channels are printed in
immediately successive lines.

See also: wtany, print
=============================================================================*/

struct wprint {
  int	depth, setlen, setcount;
  FILE	*out;
};

float calcwprint(sndobj *me, sndobj *, int);

sndobj *wprint(FILE *out, double bottom, sndobj *signal)
{
  sndobj *p;
  struct wprint *d;

  p= newsndo( calcwprint, "wprint", "wprint", signal->nch, 1, signal );
  p->skip= dontskip;
  p->private[0]= d= new(struct wprint);
  d->depth= wtdepth(bottom);
  d->setlen= 1 << d->depth;
  d->setcount= d->setlen;
  d->out= out;
  buf( 0.0, (double)d->setlen/(double)SAMPRATE, 0.0, signal->nch, p, 0 );
  return p;
}

float calcwprint(sndobj *me, sndobj *caller, int innr)
{
  struct wprint *d;
  int ch;
  
  d= (struct wprint *)me->private[0];
  if( d->setcount >= d->setlen )
  {
    int count, count2, bits;
    
    FORCH {
      BUFINCALC(0, 0, 0.0);
      fprintf( d->out, "%g", (double)INPUTC(0, ch));
      for( bits= 2; bits<= d->setlen; bits <<= 1 ) {
        BUFINCALC(0, bits-1, 0.0);
	fprintf( d->out, "\t\t%g", (double)INPUTC(0, ch));
      }
      fprintf(d->out, "\n");
    }
    for( count= 2; count < d->setlen; count += 2 )
      FORCH {
        BUFINCALC(0, count, 0.0);
        fprintf( d->out, "%g", (double)INPUTC(0, ch));
        count2= pow2div(count);
        for( bits= 2; count2> 1; bits <<= 1, --count2 ) {
	  BUFINCALC(0, count+bits-1, 0.0);
	  fprintf( d->out, "\t\t%g", (double)INPUTC(0, ch));
	}
        fprintf(d->out, "\n");
      }
    d->setcount= 0;
  }
  ++d->setcount;
  BUFINCALC(0, 0, 0.0);
  FORCH
    OUTPUT(ch, INPUTC(0, ch));
  BUFINCR(0, 1);
  CALCRETURN;
}


/*=============================================================================
    wtrunc (bottom) <freq> <crossover> <wtsignal>

Truncates the wavelet series.  Depending on the sign of <freq>, the wavelet
scales corresponding to frequencies above (for freq > 0) or below <freq> are
suppressed.  <crossover> == 0 yields strict truncation - all wavelet
coefficients above/below <freq> are set to 0.  Otherwise, the coefficients are
filtered with a hyperbolic tangent, ie multiplied with the fator 0.5 * (1 +
tanh( waveletfreq-<freq>/|<crossover>|)).

For Haar wavelets, this does not result in a highpass or lowpass as one might
expect but creates an interesting effect.

See also: wselect, wthresh, highpass, lowpass, shelve
=============================================================================*/

struct wtrunc {
  int	depth, setlen, setcount;
};

float calcwtrunc(sndobj *me, sndobj *, int);
void skipwtrunc(sndobj *me, sndobj *, int);

sndobj *wtrunc(double bottom, sndobj *freq, sndobj *crossover, sndobj *wtsignal)
{
  sndobj *p;
  struct wtrunc *d;

  p= newsndo( calcwtrunc, "wtrunc", "wtrunc", wtsignal->nch, 3, freq, crossover, wtsignal );
  p->private[0]= d= new(struct wtrunc);
  d->depth= wtdepth(bottom);
  d->setlen= 1 << d->depth;
  d->setcount= 1;
  return p;
}

float calcwtrunc(sndobj *me, sndobj *caller, int innr)
{
  struct wtrunc *d;
  float freq, crossover, datafreq;
  int bits, ch;
  
  d= (struct wtrunc *)me->private[0];
  freq= fabsf(INCALC(0));
  crossover= fabsf(INCALC(1));
  INCALC(2);
  datafreq= (float)SAMPRATE/2.0f;
  for( bits= d->setcount; bits && (bits&1)==0; bits >>= 1 )
    datafreq /= 2.0f;
  if( !crossover )
    if( datafreq==freq )
      FORCH
        OUTPUT(ch, 0.5f*INPUTC(2, ch));
    else if( datafreq > freq )
      FORCH
        OUTPUT(ch, INPUT(0) > 0 ? INPUTC(2, ch): 0.0);
    else
      FORCH
        OUTPUT(ch, INPUT(0) > 0 ? 0.0 : INPUTC(2, ch));
  else
    if( INPUT(0) > 0 )
      FORCH
        OUTPUT(ch, 0.5f*(1.0f+tanhf( (datafreq-freq)/crossover ))*INPUTC(2, ch) );
    else
      FORCH
        OUTPUT(ch, 0.5f*(1.0f-tanhf( (datafreq-freq)/crossover ))*INPUTC(2, ch) );
  if( ++d->setcount > d->setlen )
    d->setlen= 1;
  CALCRETURN;
}

void skipwtrunc(sndobj *me, sndobj *caller, int innr)
{
  struct wtrunc *d;

  d= (struct wtrunc *)me->private[0];
  if( ++d->setcount > d->setlen )
        d->setlen= 1;
  INSKIP(0);
  INSKIP(1);
  INSKIP(2);
  return;
}


/*=============================================================================
    wthresh <threshold> <wtsignal>

Applies a threshold to the wavelet coefficients.  The coefficients in the
wavelet-transformed signal <wtsignal> are compared to <threshold> and set to 0
if they are below the threshold in terms of absolute values.  (The modulus of
both <wtsignal> and <threshold> is taken.)  The channel index of <threshold> is
incremented along with that of <wtsignal> and wraps if necessary.

See also: wselect, wtrunc
=============================================================================*/

float calcwthresh(sndobj *me, sndobj *, int);

sndobj *wthresh( sndobj *threshold, sndobj *wtsignal )
{
  return newsndo( calcwthresh, "wthresh", "wthresh", wtsignal->nch, 2, threshold, wtsignal );
}

float calcwthresh(sndobj *me, sndobj *caller, int innr)
{
  float threshold, sig;
  int thch, ch;
  
  thch= 0;
  INCALC(0);
  INCALC(1);
  FORCH
  {
    threshold= fabsf(INPUTC(0, thch));
    sig= fabsf(INPUTC(1, ch));
    if( sig >= threshold )
      OUTPUT(ch, INPUTC(1, ch));
    else
      OUTPUT(ch, 0.0f);
    if( ++thch >= me->in[0]->nch )
      thch= 0;
  }
  CALCRETURN;
}



/*=============================================================================
    wshift (bottom) (mode) (nshifts) <signal>

Shifts the wavelet data up or down in frequency.  This can be used to transpose
sounds up or down by a number of octaves, but does not work well for sines
which are far from the sampling rate divided by a power-of-two.  As there are a
different amount of wavelet data at different scales, The data have to be
reduced or extended somehow.  (mode) determines how this is done.
WSHIFT_SELECT causes the first half/fourth/... of the smaller-scale data to be
used for shifting down in frequency by one/two/... .  For upshifting, the
larger-scale data are repeated the appropriate number of times.  WSHIFT_SPREAD
averages over all halfs/fourths/... of the smaller-scale data for downshifting
and interpolates between data sets for upshifting.  WSHIFT_SUCC averages /
interpolates between successive data, which creates a metallic distortion.  The
other two modes create artifacts at dataset boundaries, but these are no longer
audible when shifted wavelets are added to the original signal (as by
wspread and wlinear).  (nshifts) determines by how many bands the data are
shifted; >0 means upward in frequency.  The differences in normalisation (by a
factor sqrt(2) between adjacent bands) is taken into account.

See also: wspread, wlinear
=============================================================================*/

struct wshift {
  int    mode, nshifts, depth, depth1, setlen, setcount;
  double sqrt2fact;
};

float calcwshiftup(sndobj *me, sndobj *, int);
float calcwshiftdown(sndobj *me, sndobj *, int);
void skipwshift(sndobj *me, sndobj *, int);

sndobj *wshift( double bottom, int mode, int nshifts, sndobj *signal )
{
  sndobj *p;
  struct wshift *d;

  if( !nshifts )
    return signal;
  p= newsndo( calcwshiftup, "wshift", "wshift", signal->nch, 1, signal );
  p->skip= skipwshift;
  p->private[0]= d= new(struct wshift);
  d->mode= mode;
  d->nshifts=  nshifts;
  // For shifting downwards in frequency, sqrt2fact also contains the
  // normalisation of the average.
  if( nshifts < 0 && mode == WSHIFT_SELECT )
    d->sqrt2fact= pow( 2.0, 0.5*abs(nshifts) );
  else
    d->sqrt2fact= pow( 2.0, -0.5*abs(nshifts) );
  d->depth= wtdepth(bottom);
  d->setlen= 1 << d->depth;
  buf( 0.0, 2.0*(double)d->setlen/(double)SAMPRATE, (double)d->setlen/(double)SAMPRATE, p->nch, p, 0 );
  d->setcount= 1;
  if( nshifts< 0 )
    p->calc= calcwshiftdown;
  return p;
}

void dowshift( sndobj *me, int mode, int nshifts, int setcount,
    		int totaldepth, int signalin, float *dest );

float calcwshiftup(sndobj *me, sndobj *caller, int innr)
{
  struct wshift *d;
  double coeff;
  int ch, depth, srcdepth;

  d= (struct wshift *)me->private[0];
  depth= pow2div(d->setcount);
  srcdepth= depth + d->nshifts;
  if( srcdepth< d->depth )
    coeff= d->sqrt2fact;
  else
    coeff= d->sqrt2fact * M_SQRT2;
  dowshift(me, d->mode, d->nshifts, d->setcount, d->depth, 0, me->ch);
  FORCH
    me->ch[ch] *= coeff;
  if( ++d->setcount > d->setlen ) {
    d->setcount= 1;
    BUFINCR(0, d->setlen);
  }
  CALCRETURN;
}

float calcwshiftdown(sndobj *me, sndobj *caller, int innr)
{
  struct wshift *d;
  double coeff;
  int ch;

  d= (struct wshift *)me->private[0];
  if( d->setcount != d->setlen )
    coeff= d->sqrt2fact;
  else
    coeff= d->sqrt2fact * (d->mode == WSHIFT_SELECT ? M_SQRT1_2 : M_SQRT2);
  dowshift(me, d->mode, d->nshifts, d->setcount, d->depth, 0, me->ch);
  FORCH
    me->ch[ch] *= coeff;
  if( ++d->setcount > d->setlen ) {
    d->setcount= 1;
    BUFINCR(0, d->setlen );
  }
  CALCRETURN;
}

void skipwshift(sndobj *me, sndobj *caller, int innr)
{
  struct wshift *d;

  d= (struct wshift *)me->private[0];
  if( ++d->setcount > d->setlen ) {
    d->setcount= 1;
    BUFINCR(0, d->setlen );
  }
  return;
}

// Calculate shifted wavelet coefficients, without any renormalisation.
// For shifts downward in frequency, the required linear combination of
// wavelet data of the source band is computed without renormalisation.
// For upward shifts and mode!=WSHIFT_SELECT, the sum is computed.  The average
// has to be obtained by a division in the calling function.
// This function assumes the current read position on the relevant (buffered)
// signal input <signalin> to be at the start of the data set.
// Arguments: me - caller; mode - shift mode; nshifts - number of bands to 
// shift by; setcount - 1-based data set index; totaldepth - depth of wavelet
// transform; signalin - index of signal input to shift; dest - float array to
// write result to, needs to provide space for all channels.

void dowshift( sndobj *me, int mode, int nshifts, int setcount,
    		int totaldepth, int signalin, float *dest )
{
  float nextwgt;
  int srcdepth, currdepth, setlen, navg, ch, srcoffset, sample, count;

  if( nshifts == 0 ) {
    BUFINCALC(signalin, setcount-1, 0.0);
    FORCH
      dest[ch]= INPUTC(signalin, ch);
    return;
  }
  navg= 1 << abs(nshifts);
  currdepth= pow2div(setcount);
  srcdepth= currdepth + nshifts;
  setlen= 1 << totaldepth;
  if( srcdepth > totaldepth || srcdepth < 0 ) {
    FORCH
      dest[ch]= 0.0f;
    return;
  }
  if( nshifts > 0 )  // upward shift
  {
    if( mode == WSHIFT_SUCC ) {
      if( srcdepth == totaldepth ) {
	srcoffset= setlen;
	sample= setlen - 1;
      }
      else {
	srcoffset= 1 << (srcdepth+1);
	sample= (setcount&~(srcoffset-1)) + srcoffset/2 - 1;
      }
      nextwgt= (float)((setcount >> (currdepth+1))&(navg-1))/(float)navg;
    }
    else {
      sample= ((setcount << nshifts) - 1) & (setlen-1);
      srcoffset= setlen;
      if( mode != WSHIFT_SELECT )
	nextwgt= (float)(setcount - (1<<currdepth))/(float)setlen;
      else
	nextwgt= 0.0;
    }
    BUFINCALC(signalin, sample, 0.0);
    FORCH
      dest[ch]= INPUTC(signalin, ch);
    if( nextwgt ) {
      BUFINCALC(signalin, sample + srcoffset, 0.0);
      FORCH
	dest[ch] += nextwgt * (INPUTC(signalin, ch) - dest[ch]);
    }
  }
  else    // downward shift
  {
    nshifts= -nshifts;
    if( mode == WSHIFT_SUCC ) {  // add successive data
      sample= setcount-1 - (((1 << nshifts) - 1) << srcdepth);
      srcoffset= 1 << (srcdepth+1);
    }
    else {
      sample= (setcount >> nshifts) - 1;
      if( currdepth != totaldepth )
	srcoffset= setlen >> nshifts;
      else
	srcoffset= setlen >> (nshifts-1);
    }
    BUFINCALC(signalin, sample, 0.0);
    FORCH
      dest[ch]= INPUTC(signalin, ch);
    if( mode == WSHIFT_SELECT )
      return;
    if( currdepth == totaldepth )
      count= navg/2-1;
    else
      count= navg-1;
    for( sample += srcoffset; count> 0; --count, sample += srcoffset ) {
      BUFINCALC(signalin, sample, 0.0);
      FORCH
	dest[ch] += INPUTC(signalin, ch);
    }
  }
}



/*=============================================================================
    wspread (bottom) (mode) <spreadup> <spreaddown> <signal>

This object spreads the value of each wavelet datum to wavelets of different
scales.  Its value multiplied with <spreadup> is added to the wavelet with the
next smaller scale, multiplied with <spreadup>^2 to the next but one and so on,
and analogously with <spreaddown> and the larger-scale wavelets.  For spreading
to lower frequencies (larger scales), the original wavlet datum is averaged
over the scale of the one it is added to.  <spreadup> and <spreaddown> are 
adjusted by factors which account for differences in normalisation between
wavelets.  (mode) can take the same values as for wshift.

This creates an interesting, polyphonous effect similar to spectral spread.
Some types of wavelets cause artifacts.  For (mode)= WSHIFT_SUCC, an interesting
distortion is created.

See also: wshift, wlinear
=============================================================================*/

struct wspread {
  int	mode, depth, depth1, setlen, setcount;
};

float calcwspread(sndobj *me, sndobj *, int);
void skipwspread(sndobj *me, sndobj *, int);

sndobj *wspread( double bottom, int mode, sndobj *spreadup, sndobj *spreaddown, sndobj *signal )
{
  sndobj *p;
  struct wspread *d;

  p= newsndo( calcwspread, "wspread", "wspread", signal->nch, 3, spreadup, spreaddown, signal );
  p->skip= skipwspread;
  p->private[0]= d= new(struct wspread);
  d->mode= mode;
  d->depth= wtdepth(bottom);
  d->depth1= d->depth+1;
  d->setlen= 1 << d->depth;
  buf( 0.0, 2.0*(double)d->setlen/(double)SAMPRATE, (double)d->setlen/(double)SAMPRATE, p->nch, p, 2 );
  d->setcount= 1;
  return p;
}

float calcwspread(sndobj *me, sndobj *caller, int innr)
{
  struct wspread *d;
  float buf[PLENTY];
  float upcoeff, downcoeff, coeff;
  int ch, currdepth, nshifts;

  d= (struct wspread *)me->private[0];
  upcoeff= fabsf(INCALC(0));
  downcoeff= fabsf(INCALC(1));
  upcoeff *= M_SQRT1_2;
  downcoeff *= d->mode==WSHIFT_SELECT? M_SQRT2 : M_SQRT1_2;
  currdepth= pow2div(d->setcount);

  BUFINCALC(2, d->setcount-1, 0.0);
  FORCH
    me->ch[ch]= INPUTC(2, ch);
  for( nshifts= 1, coeff= upcoeff; currdepth+nshifts < d->depth; ++nshifts, coeff *= upcoeff )
  {
    dowshift(me, d->mode, nshifts, d->setcount, d->depth, 2, buf);
    FORCH
      me->ch[ch] += coeff * buf[ch];
  }
  if( currdepth+nshifts == d->depth ) {
    coeff *= M_SQRT2;
    dowshift(me, d->mode, nshifts, d->setcount, d->depth, 2, buf);
    FORCH
      me->ch[ch] += coeff * buf[ch];
  }
  coeff= downcoeff;
  if( currdepth==d->depth )
    coeff *= (d->mode == WSHIFT_SELECT ? M_SQRT1_2 : M_SQRT2);
  for( nshifts= 1; nshifts <= currdepth; ++nshifts, coeff *= downcoeff )
  {
    dowshift(me, d->mode, -nshifts, d->setcount, d->depth, 2, buf);
    FORCH
      me->ch[ch] += coeff * buf[ch];
  }
  if( ++d->setcount > d->setlen ) {
    d->setcount= 1;
    BUFINCR(2, d->setlen);
  }
  CALCRETURN;
}

void skipwspread(sndobj *me, sndobj *caller, int innr)
{
  struct wspread *d;

  INSKIP(0);
  INSKIP(1);
  d= (struct wspread *)me->private[0];
  if( ++d->setcount > d->setlen ) {
    d->setcount= 1;
    BUFINCR(2, d->setlen );
  }
  return;
}


/*=============================================================================
    wlinear (bottom) (mode) (dimension) [matrix] <signal>

Applies a linear transformation to a wavelet-transformed signal.  <matrix>
should point to an array of (dim) by (dim)+1 coefficients, stored in row-major
format (ie first all values of the first row and so on).  Each row of (dim)+1
values describes the linear combination of wavelet data from the input signal
to obtain the output wavelet data.  To form the linear combination, wavelet
data are shifted according to (mode), which is explained in the documentation
of wshift.  The last value in a row is a constant which is added.  If (dim) is
smaller than the depth of the wavelet transform, only the (dim) smallest scales
are transformed.  The first value in each row or column refers to the highest
frequency band.

Alternatively, when (dim) is negative, [matrix] can contain a one-dimensional
array of size 2*|(dim)|-1.  These values give the coefficients of data from
higher and lower frequency bands which are mixed in with each wavelet datum.
The first value corresponds to the highest frequency band, |(dim)|-1 bands up
from the current datum; the value with index |(dim)|-1 is the prefactor of the
original datum; and the following ones correspond to the lower frequency bands.

See also: wshift, wspread
=============================================================================*/

struct wlinear {
  int	dim, mode, setlen, depth, setcount;
  float *matrix, *c;
};

float calcwlinear(sndobj *me, sndobj *, int);
void skipwlinear(sndobj *me, sndobj *, int);

sndobj *wlinear( double bottom, int mode, int dim, float *matrix, sndobj *signal )
{
  sndobj *p;
  struct wlinear *d;
  double downcoeffbase;
  int srcind, destind;

  if( !dim )
    return signal;
  p= newsndo(calcwlinear, "wlinear", "wlinear", signal->nch, 1, signal);
  p->private[0]= d= new(struct wlinear);
  d->depth= wtdepth(bottom);
  d->setlen= 1 << d->depth;
  d->setcount= 1;
  if( dim > 0 ) {
    d->dim= dim;
    if( d->dim > d->depth+1 ) {
      fprintf(stderr, "wlinear: Warning: matrix dimension (%d) exceeds number of frequency bands (%d).  Ignoring surplus entries.\n", d->dim, d->depth+1);
      d->dim= d->depth+1;
    }
  }
  else
    d->dim= d->depth+1;
  d->mode= mode;
  downcoeffbase= mode==WSHIFT_SELECT? 2.0 : 0.5;
  p->private[1]= d->matrix= news(float, d->dim*(d->dim+1));
  d->c= d->matrix + d->dim*d->dim;
  for( destind= 0; destind< d->dim; ++destind )
    for( srcind= 0; srcind< d->dim; ++srcind ) {
      if( dim > 0 )
	d->matrix[destind*d->dim+srcind]= matrix[destind*(dim+1) + srcind];
      else if( abs(srcind-destind) < -dim )
	d->matrix[destind*d->dim+srcind]= matrix[-dim-1+srcind-destind];
      else
	d->matrix[destind*d->dim+srcind]= 0.0f;
      if( destind > srcind )
	if( destind == d->depth )
	  d->matrix[destind*d->dim+srcind] *= pow(downcoeffbase, 0.5*(destind-srcind-1));
	else
	  d->matrix[destind*d->dim+srcind] *= pow(downcoeffbase, 0.5*(destind-srcind));
      else if( destind < srcind )
	if( srcind == d->depth )
	  d->matrix[destind*d->dim+srcind] *= pow(0.5, 0.5*(srcind-destind-1));
	else
	  d->matrix[destind*d->dim+srcind] *= pow(0.5, 0.5*(srcind-destind));
    }
  if( dim > 0 )
    for( destind= 0; destind< d->dim; ++destind )
      d->c[destind]= matrix[destind*(dim+1) + dim];
  buf( 0.0, 2.0*(double)d->setlen/(double)SAMPRATE, (double)d->setlen/(double)SAMPRATE, p->nch, p, 0 );
  return p;
}

float calcwlinear(sndobj *me, sndobj *caller, int innr)
{
  struct wlinear *d;
  float buf[PLENTY];
  int ch, currdepth, srcdepth;

  d= (struct wlinear *)me->private[0];
  currdepth= pow2div(d->setcount);
  if( currdepth < d->dim )
  {
    FORCH
      me->ch[ch]= d->c[currdepth];
    for( srcdepth= 0; srcdepth< d->dim; ++srcdepth )
      if( d->matrix[currdepth*d->dim+srcdepth] != 0.0f )
      {
	dowshift(me, d->mode, srcdepth-currdepth, d->setcount, d->depth, 0, buf);
	FORCH
	  me->ch[ch] += d->matrix[currdepth*d->dim+srcdepth] * buf[ch];
      }
  }
  else {
    BUFINCALC(0, d->setcount-1, 0.0);
    FORCH
      me->ch[ch]= INPUTC(0, ch);
  }
  if( ++d->setcount > d->setlen ) {
    d->setcount= 1;
    BUFINCR(0, d->setlen);
  }
  CALCRETURN;
}

void skipwlinear(sndobj *me, sndobj *caller, int innr)
{
  struct wspread *d;

  d= (struct wspread *)me->private[0];
  if( ++d->setcount > d->setlen ) {
    d->setcount= 1;
    BUFINCR(0, d->setlen );
  }
  return;
}



/*=============================================================================
    wuninterleave (bottom) <signal>

This object extracts the wavelet coefficients from a wavelet signal and outputs
them in separate channels.  <signal> may have multiple channels if the number
of wavelet coefficients is not too large.  (The number of output channels is
limited to PLENTY.)  The number of output channels for each input channel is
depth+1, where depth is approximately the two's logarithm of the sampling rate
divided by <bottom>, or the two's logarithm of the wavelet data set length.
(See the top of sndwt.c for a more thorough description of how the length of
the data set is computed.)  The smallest-scale wavelet data (highest frequency
band) are output on the first channel.  The output data have the sampling rate
as their time resolution while the wavelet data are coarser in time, so the
wavelet data are interpolated linearly.

See also: winterleave, wasyncinterleave, wtany
=============================================================================*/

struct wuninterleave {
  int setcount, setlen, depth, initialised;
  float last[PLENTY], next[PLENTY];
};

float calcwuninterleave(sndobj *me, sndobj *, int);
void skipwuninterleave(sndobj *me, sndobj *, int);

sndobj *wuninterleave( double bottom, sndobj *signal )
{
  sndobj *p;
  struct wuninterleave *d;
  int nch, depth;

  d= new(struct wuninterleave);
  depth= wtdepth(bottom);
  d->depth= depth;
  d->setlen= (1 << depth);
  d->setcount= 1;
  d->initialised= 0;
  nch= (depth+1)*signal->nch > PLENTY ? PLENTY : (depth+1)*signal->nch;
  p= newsndo(calcwuninterleave, "wuninterleave", "wuninterleave", nch, 1, signal);
  p->skip= skipwuninterleave;
  p->private[0]= d;
  buf( 0.0, 1.5*(double)d->setlen/(double)SAMPRATE, 0.0, (nch+depth)/(depth+1), p, 0 );
  return p;
}

float calcwuninterleave(sndobj *me, sndobj *caller, int innr)
{
  struct wuninterleave *d;
  float nextwgt;
  int ch, scale, currdepth, count;

  d= (struct wuninterleave *)me->private[0];
  if( !d->initialised )
  {
    d->initialised= 1;
    for( count= 0, scale= 1, ch= 0; count< me->nch; ++count ) {
      BUFINCALC( 0, scale-1, 0.0 );
      d->next[count]= INPUTC(0, ch);  // first value in this band
      d->last[count]= INPUTC(0, ch);
      if( scale==d->setlen ) {
        scale= 1;
        ++ch;
      }
      else
        scale *= 2;
    }
  }
  currdepth= pow2div( d->setcount );
  if( currdepth < d->depth ) {
    BUFINCALC(0, 1<<(currdepth+1), 0.0);
    for( ch= 0, count= currdepth; count< me->nch; ++ch, count += d->depth+1 ) {
      d->last[count]= d->next[count];
      d->next[count]= INPUTC(0, ch);
    }
    if( currdepth==d->depth-1 ) {
      BUFINCALC(0, d->setlen/2+d->setlen, 0.0);
      for( ch= 0, count= currdepth+1; count< me->nch;
	  			++ch, count += d->depth+1 ) {
	d->last[count]= d->next[count];
	d->next[count]= INPUTC(0, ch);
      }
    }
  }
  for( ch= 0, scale= 2; ch< me->nch; ++ch )
    if( scale> d->setlen ) {
      scale= d->setlen;
      nextwgt= (float)((d->setcount+scale/2)&(scale-1))/(float)scale;
      OUTPUT(ch, d->last[ch] + nextwgt*(d->next[ch] - d->last[ch]));
      scale= 2;
    }
    else {
      nextwgt= (float)((d->setcount+scale/2)&(scale-1))/(float)scale;
      OUTPUT(ch, d->last[ch] + nextwgt*(d->next[ch] - d->last[ch]));
      scale *= 2;
    }
  if( ++d->setcount > d->setlen )
    d->setcount= 1;
  BUFINCR(0, 1);
  CALCRETURN;
}

void skipwuninterleave(sndobj *me, sndobj *caller, int innr)
{
  struct wuninterleave *d;

  d= (struct wuninterleave *)me->private[0];
  if( ++d->setcount > d->setlen )
    d->setcount= 1;
  BUFINCR(0, 1);
  return;
}


/*=============================================================================
    winterleave (bottom) <signal>

The reverse of wuninterleave.  The wavelet data in <signal> are supposed to be
organised in channels (smallest scale in first channel) and are output in the
usual wavelet data channel.  If <signal> contains several full sets of wavelet
data (with floor(log2(SAMPRATE/<bottom>)+0.5)+1 channels each), the output has
multiple channels.  If not even one full set is available, a warning is printed
and zero returned.  The time resolution of the input wavelet data is reduced by
subsampling, not averaging.

See also: wuninterleave, wasyncinterleave, wtany
=============================================================================*/

struct winterleave {
  int setcount, setlen, depth;
};

float calcwinterleave(sndobj *me, sndobj *, int);
void skipwinterleave(sndobj *me, sndobj *, int);

sndobj *winterleave(double bottom, sndobj *signal)
{
  sndobj *p;
  struct winterleave *d;
  int depth, nch;

  depth= wtdepth(bottom);
  nch= signal->nch/(depth+1);
  if( !nch ) {
    fprintf(stderr, "winterleave: warning: no complete set of data in input signal.  Did you give the right bottom frequency?  Returning 0.\n");
    return c(0);
  }
  p= newsndo(calcwinterleave, "winterleave", "winterleave", nch, 1, signal);
  p->private[0]= d= new(struct winterleave);
  d->depth= depth;
  d->setlen= 1<<depth;
  d->setcount= 1;
  return p;
}

float calcwinterleave(sndobj *me, sndobj *caller, int innr)
{
  struct winterleave *d;
  int ch, currdepth;

  d= (struct winterleave *)me->private[0];
  currdepth= pow2div(d->setcount);
  INCALC(0);
  for( ch= 0; ch< me->nch; ++ch, currdepth+=d->depth+1 )
    OUTPUT(ch, INPUTC(0, currdepth));
  if( ++d->setcount > d->setlen )
    d->setcount= 1;
  CALCRETURN;
}

void skipwinterleave(sndobj *me, sndobj *caller, int innr)
{
  struct winterleave *d;

  INSKIP(0);
  d= (struct winterleave *)me->private[0];
  if( ++d->setcount > d->setlen )
      d->setcount= 1;
  return;
}


/*=============================================================================
    wselect (bottom) <mode> <mask> <signal>

Allows to zero specific wavelets.  <mode> and <mask> are both interpreted as
integer values; they are rounded if necessary.  Both are read only at the start
of each wavelet data set, that is at the start of each interval of length
1/(bottom).  The output should be muted when they change, or clicks may be
heard.

<mode> specifies the mode of operation and the meaning of <mask>:  for
WSELECT_INDEX, bit 0 of <mask> refers to the highest wavelet band; for the
other modes, it corresponds to the one largest in magnitude.  WSELECT_SORTMAX,
WSELECT_SORTMEAN and WSELECT_SORTMIN differ in how the magnitude is determined:
the maximum, mean or minimum over a dataset is taken.  The determination of
which wavelet bands to keep is applied independently to each channel of
<signal>.

If <mode> is ORed with WSELECT_SAMERMS, changes in RMS amplitude due to the
selection of wavelet bands are compensated.  This causes *uncorrelated* signals
to have the same RMS amplitude after wselect as before.

See also: wtrunc, wthresh, wfreq
=============================================================================*/

struct wselect {
  int		setcount, setlen, depth, depth1, samerms;
  unsigned	realmask[PLENTY];
  double	scale[PLENTY];
};
struct wselectsort {
  float value;
  int depth;
};

float calcwselect(sndobj *me, sndobj *, int);

sndobj *wselect( double bottom, sndobj *mode, sndobj *mask, sndobj *signal)
{
  sndobj *p;
  struct wselect *d;
  int ch;

  p= newsndo( calcwselect, "wselect", "wselect", signal->nch, 3, signal, mode, mask );
  p->skip= dontskip;
  p->private[0]= d= new(struct wselect);
  d->depth= wtdepth(bottom);
  d->depth1= d->depth+1;
  d->setlen= 1<<d->depth;
  buf( 0.0, (double)d->setlen/(double)SAMPRATE, (double)d->setlen/(double)SAMPRATE, p->nch, p, 0 );
  d->setcount= 1;
  d->samerms= 0;
  return p;
}

int cmpwselectsort(const void *, const void *);

#define WSELECT_MINMODE	WSELECT_INDEX
#define WSELECT_MAXMODE WSELECT_SORTRMS

float calcwselect(sndobj *me, sndobj *caller, int innr)
{
  struct wselect *d;
  unsigned mask;
  int ch, depth, depthbit, mode;

  d= (struct wselect *)me->private[0];
  if( d->setcount == 1 )
  {
    mode= (int)floorf(INCALC(1) + 0.5f);
    d->samerms= mode & WSELECT_SAMERMS;
    mode= mode & ~WSELECT_SAMERMS;
    if( mode < WSELECT_MINMODE || mode > WSELECT_MAXMODE )
      mode= WSELECT_INDEX;
    mask= (unsigned)(int)floorf(INCALC(2) + 0.5f);
    if( d->samerms && (mask&(d->setlen*2-1))==0 )
      d->samerms= 0;
    if( mode != WSELECT_INDEX )
    {
      struct wselectsort *sortarray;
      int offset, sample;

      sortarray= news(struct wselectsort, me->nch*d->depth1);
      for( depth= 0, offset= 2; depth<= d->depth; ++depth, offset <<= 1 ) {
	FORCH
	  sortarray[ch*d->depth1+depth].depth= depth;
	if( mode != WSELECT_SORTMIN )
	  FORCH
	    sortarray[ch*d->depth1+depth].value= 0.0;
	else
	  FORCH
	    sortarray[ch*d->depth1+depth].value= FLT_MAX;
	for( sample= offset/2-1; sample < d->setlen; sample += offset ) {
	  BUFINCALC(0, sample, 0.0);
	  if( mode==WSELECT_SORTMEAN )
	    FORCH
	      sortarray[ch*d->depth1+depth].value += INPUTC(0, ch);
	  else if( mode==WSELECT_SORTRMS )
	    FORCH
	      sortarray[ch*d->depth1+depth].value += INPUTC(0, ch)*INPUTC(0, ch);
	  else if( mode==WSELECT_SORTMAX )
	    FORCH
	      if( fabsf(INPUTC(0, ch)) > sortarray[ch*d->depth1+depth].value )
		sortarray[ch*d->depth1+depth].value= fabsf(INPUTC(0, ch));
	      else;
	  else
	    FORCH
	      if( fabsf(INPUTC(0, ch)) < sortarray[ch*d->depth1+depth].value )
		sortarray[ch*d->depth1+depth].value= fabsf(INPUTC(0, ch));
	}
	if( (mode==WSELECT_SORTMEAN || mode==WSELECT_SORTRMS) 
	    && depth< d->depth-1 )
	  FORCH
	    sortarray[ch*d->depth1+depth].value /= d->setlen/offset;
	if( mode == WSELECT_SORTRMS )
	  FORCH
	    sortarray[ch*d->depth1+depth].value=
				  sqrtf(sortarray[ch*d->depth1+depth].value);
      }
      FORCH {
	qsort( sortarray+ch*d->depth1, d->depth1, sizeof(struct wselectsort),
			    cmpwselectsort );
	d->realmask[ch]= 0;
	for( depth= 0, depthbit= 1; depth<= d->depth; ++depth, depthbit <<= 1 )
	  if( (mask & depthbit) != 0 )
	    d->realmask[ch] |= 1 << sortarray[ch*d->depth1+depth].depth;
      }
      free(sortarray);
      if( d->samerms )
	FORCH
	  if( (d->realmask[ch] & d->setlen) != 0 ) {
	    d->scale[ch]= 1.0f / sqrtf(1.0f/(float)d->setlen
			  + bitinverse(d->realmask[ch] & (d->setlen-1)));
	  }
	  else
	    d->scale[ch]= 1.0f / 
	      		  sqrtf(bitinverse(d->realmask[ch] & (d->setlen-1)));
    }
    else {
      FORCH
	d->realmask[ch]= mask;
      if( d->samerms ) {
	float scale;
	
	if( (mask & d->setlen) != 0 ) {
	  scale= 1.0f / sqrtf(1.0f/(float)d->setlen + 
	      				bitinverse(mask & (d->setlen-1)));
	}
	else
	  scale= 1.0f / sqrtf(bitinverse(mask & (d->setlen-1)));
	FORCH
	  d->scale[ch]= scale;
      }
    }
  }
  else {
    INSKIP(1);
    INSKIP(2);
  }
  depth= pow2div(d->setcount);
  depthbit= 1<<depth;
  BUFINCALC(0, d->setcount-1, 0.0);
  FORCH
    if( (d->realmask[ch] & depthbit) != 0 )
      if( d->samerms )
	OUTPUT(ch, d->scale[ch]*INPUTC(0, ch));
      else
	OUTPUT(ch, INPUTC(0, ch));
    else
      OUTPUT(ch, 0.0);
  if( ++d->setcount > d->setlen ) {
    d->setcount= 1;
    BUFINCR(0, d->setlen);
  }
  CALCRETURN;
}

int cmpwselectsort(const void *a, const void *b)
{
  struct wselectsort *c, *d;

  c= (struct wselectsort *)a;
  d= (struct wselectsort *)b;
  if( c->value > d->value )
    return -1;
  else if( c->value < d->value )
    return 1;
  else
    return 0;
}


/*=============================================================================
    wfreq (bottom) (mode) <freq> <mask>

Calculates a wavelet band mask (for wselect in mode WSELECT_INDEX) and a
resampling factor from a frequency and an arbitrarily shifted mask.  Intended
for making sounds synthesised with wavelets independent of the sampling rate.
(bottom) is, as always, the bottom frequency of the transformation.  Here, it
may be passed as zero and will then be ignored.

The low byte of (mode) may be WFREQ_CONST or WFREQ_NOISE, depending on whether
the input to wselect is more or less constant or more or less random (fbm can
tend towards either case, depending on H).  The relevant difference between the
two modes is that wavelet synthesis of a constant yields the frequency which is
the upper limit of the band (which is then also the centre frequency), while a
synthesis of random data yields the whole band, whose centre frequency is in
the middle.  If the wavelet spectrum to be synthesised contains both contant
and random parts, WFREQ_CONST is probably a better choice as the human hearing
tends to pick out the harmonic part.  

<freq> gives the requested frequency.  <mask> is the wavelet band mask for
wselect (with mode==WSELECT_INDEX).  It may be shifted left or right
arbitrarily.  Depending on flags with which (mode) can be ORed, the frequency
is the requested centre, top or bottom frequency of the output of wselect:
WFREQ_CENTER (default, =0), WFREQ_TOP, WFREQ_BOTTOM.  The centre frequency is
defined on a logarithmic scale, ie as the geometric mean of the top and bottom
frequencies.

The output has two channels.  The first gives the mask shifted such that the
centre frequency of the output of wselect will be below, but within an octave
of the requested.  The second channel contains the resampling factor in [1,2)
with which that output has to be resampled to achieve exactly the requested
centre frequency.

See also: wselect
=============================================================================*/

struct wfreq {
  int	 mode;
  unsigned maskmask;
};

float calcwfreq(sndobj *me, sndobj *caller, int innr);

sndobj *wfreq(double bottom, int mode, sndobj *freq, sndobj *mask)
{
  sndobj *p;
  struct wfreq *d;

  p= newsndo(calcwfreq, "wfreq", "wfreq", 2, 2, freq, mask);
  p->private[0]= d= new(struct wfreq);
  d->mode= mode;
  if( bottom > 0.0 )
    d->maskmask= (1 << (wtdepth(bottom)+1)) - 1;
  else
    d->maskmask= 0;
  return p;
}

float calcwfreq(sndobj *me, sndobj *caller, int innr)
{
  struct wfreq *d;
  double reqfreq, reqoct, diffoct;
  unsigned mask, shifter;
  int octaves, shiftoct;

  d= (struct wfreq *)me->private[0];
  reqfreq= fabs(INCALC(0));
  if( reqfreq< 1.0 )
    reqfreq= 1.0;
  mask= (unsigned)(int)INCALC(1);
  if( !mask ) {
    OUTPUT(0, 0.0f);
    OUTPUT(1, 1.0f);
    CALCRETURN;
  }
  while( (mask&1) == 0 )
    mask >>= 1;
  for( octaves= 1, shifter= mask>>1; shifter; shifter >>= 1 )
    ++octaves;
  reqoct= log(SAMPRATE/2/reqfreq)/M_LN2;
  if( (d->mode & WFREQ_ALIGNMASK) == WFREQ_BOTTOM )
    diffoct= reqoct - octaves;
  else if( (d->mode & WFREQ_ALIGNMASK) == WFREQ_TOP )
    diffoct= reqoct;
  else
    diffoct= reqoct - 0.5*octaves;
  if( (d->mode & WFREQ_MODEMASK)==WFREQ_CONST )
    if( (d->mode & WFREQ_ALIGNMASK) == WFREQ_BOTTOM )
      diffoct += 1.0;
    else if( (d->mode & WFREQ_ALIGNMASK) == WFREQ_TOP )
      ;
    else
      diffoct += 0.5;
  shiftoct= (int)ceil(diffoct);
  if( shiftoct > 0 )
    mask <<= shiftoct;
  else if( shiftoct < 0 )
    mask >>= -shiftoct;
  if( d->maskmask )
    mask &= d->maskmask;
  OUTPUT(0, mask);
  OUTPUT(1, exp2(shiftoct-diffoct));	// remaining freq difference -> resample
  CALCRETURN;
}




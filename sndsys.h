/*
    This file is part of sndsys, a digital signal processing system.
    Copyright (c) 2004-2007 Volker Schatz (noise at volkerschatz dot com).

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

#ifndef __sndsys_h
#define __sndsys_h

#define __USE_MISC

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

/******************************************************************************
    Conditional compilation
******************************************************************************/

//#define NO_X
//#define NO_PTHREADS
//#define NO_OSS
//#define NO_ALSA

/******************************************************************************
    Macros
******************************************************************************/

#define PLENTY          20
#define NAMESTRLEN      80
#define SAMPRATE        44100
#define HORIZON         1.0   // maximal time out of sync in either direction for non-buffering objects

#define INCALC(n)       (me->in[n]->calc(me->in[n], me, (n)))
#define INSKIP(n)       (me->in[n]->skip(me->in[n], me, (n)))
#define INPUT(n)        (me->inch[n][0])
#define INPUTC(n,c)     (me->inch[n][c])
#define BUFINCALC(n,ind,frac)   (intercalcbuf(me->in[n], me, (n), (ind), (frac)))
#define BUFINCR(n,incr)         (incrbuf(me->in[n], me, (n), (incr)))
#define BUFPTR(n,ind,pptr,pnsamp,pnch)	(getbufptr(me->in[n], me, (n), (ind), (pptr), (pnsamp), (pnch)))
#define OUTPUT(c,val)   (me->ch[c]= (val))
#define CALCRETURN      return me->ch[0]

#define FORCH           for( ch= 0; ch< me->nch; ++ch )
#define FORCH1          for( ch= 1; ch< me->nch; ++ch )

#define new(type)       ((type *)malloc(sizeof(type)))
#define news(type,n)    ((type *)malloc(sizeof(type)*(unsigned long)(n)))
#define END             ((float)-0x1234560000000000ll)
#define MAGIC           ((float)-0x1234560000000000ll)
#define SEMITONE        1.05946309435929526455      //  2^(1/12)


/******************************************************************************
    Central stuff
******************************************************************************/

typedef struct sndo {
  float         (*calc)(struct sndo *, struct sndo *, int);
  void          (*skip)(struct sndo *, struct sndo *, int);
  struct sndo   *in[PLENTY];
  float         *inch[PLENTY];
  float         ch[PLENTY];
  int           nin, nch;
  const char    *name, *type;
  void          *private[PLENTY];
  void          (*exit)(struct sndo *me);
  struct sndo   *nextinchain;
}
sndobj;

extern double starttime, stoptime;

void skipinputs(sndobj *me, sndobj *, int);
void dontskip(sndobj *me, sndobj *caller, int innr);

sndobj *newsndo(float (*calc)(sndobj *, sndobj *, int), const char *name, 
		const char *type, int nch, int nin, ...);
int addinput(sndobj *obj, sndobj *input);
#define DUMP_IN     1
#define DUMP_INNAME 2
#define DUMP_INCH   4
#define DUMP_INALL  (DUMP_IN|DUMP_INNAME|DUMP_INCH)
#define DUMP_OUTCH  8
#define DUMP_PRIV   16
#define DUMP_ALL    (DUMP_INALL|DUMP_OUTCH|DUMP_PRIV)
void dumpsndo(sndobj *p, int dumpwhat);
void disableskip();

void closeloops(void);
void prune(sndobj *top);
void makefanouts(void);

void killsndo(sndobj *victim);
void killemall(void);
void sndexecute(double t_start, double t_end, sndobj *top);
void dumpchain(void);
void dumptree(sndobj *top);

sndobj *loop( const char *name, int id, int nch );
sndobj *defloop( const char *name, int id, sndobj *source );
float calcloop( sndobj *me, sndobj *, int);

void buf( float backtime, float forwtime, float jitter, int nch, sndobj *client, int innr );
float intercalcbuf( sndobj *me, sndobj *caller, int innr, int index, float frac );
void incrbuf( sndobj *me, sndobj *caller, int innr, int incr );
void getbufptr( sndobj *me, sndobj *caller, int innr, int index, float **ptr, int *nsamples, int *nch );

void registershared( void *buf, char *name );
void *findshared( const char *name, int registerme );
int freeshared( void *buf );

float *sinelookup( long size );

typedef struct {
  int	inext, inextp;
  int 	state[56];
}
knuthrngstate;

#define RNG_ARRAYSIZE  256
typedef struct {
  struct random_data re_buf;		// reentrancy data
  char  state[RNG_ARRAYSIZE];
}
rngdesc;
rngdesc *makerng( void );
float flatrng( rngdesc *rng );
float gaussrng( rngdesc *rng );

float lininterpol( float *buf, long size, float fpos );

int getprime( int lowlimit );

long binfindf( float *sorted, long size, float val );

int dspdesc(int restore);


/******************************************************************************
    Objects
******************************************************************************/

// structs used in multiple .c files:

typedef struct wavhead {
  unsigned int      RIFF;
  unsigned int      totalsizem8;
  unsigned int      WAVE;
  unsigned int      fmt;
  unsigned int      fmtsize;
  unsigned short    audiofmt;
  unsigned short    nchannels;
  unsigned int      samplerate;
  unsigned int      byterate;
  unsigned short    blockalign;
  unsigned short    bitspersample;
  unsigned int      data;
  unsigned int      restsize;
}
wavheader;


/*-----------------------------------------------------------------------------
    sndmisc.c
-----------------------------------------------------------------------------*/

sndobj *linear(float offset, float slope, sndobj *in);
sndobj *c(float val);
sndobj *cc(int nch, float val);
sndobj *extc( int nch, float *data, int *flag );
sndobj *add(int n, ...);
sndobj *mul(int n, ...);
sndobj *arithmetic(const char *op, sndobj *arg1, sndobj *arg2 );
sndobj *ch(int ch, sndobj *in);
sndobj *c0(sndobj *in);
sndobj *c1(sndobj *in);
sndobj *c2(sndobj *in);
sndobj *c3(sndobj *in);
sndobj *chsum(sndobj *signal);
sndobj *chavg(sndobj *signal);
sndobj *chsel(int nch, sndobj *channels, sndobj *signal);
sndobj *switchboard( char *sboard, ... );
sndobj *mix( char *mixstr, ... );
sndobj *matrix( int nch, float *coeff, ... );
sndobj *mix2( sndobj *gate, sndobj *input1, sndobj *input2 );
sndobj *cdelay(double time, sndobj *in);
sndobj *sdelay(long samples, sndobj *in);
sndobj *delay( float maxdelay, sndobj *delay, sndobj *in );
sndobj *rdelay( float maxdelay, sndobj *delay, sndobj *reset, sndobj *in );
sndobj *limiter( float max, sndobj *signal );
sndobj *comptime(sndobj *subtree);
sndobj *diff(sndobj *in);
sndobj *integrate(sndobj *in);
sndobj *bifurcation( float interval, sndobj *parameter );
#define SMST_MODEMASK	0xFF
#define SMST_LINEAR     0
#define SMST_SINE       1
#define SMST_FLAGMASK	0xFFFF00
#define SMST_UP		0x100
#define SMST_DOWN	0x200
#define SMST_ABSUP	0x400
#define SMST_ABSDOWN	0x800
sndobj *smoothstep( int mode, float minstep, sndobj *steptime, sndobj *in );
sndobj *conditional( sndobj *in1, sndobj *in2, sndobj *greaterequal, sndobj *lessthan );
sndobj *smoothstart( int samples, sndobj *signal );
#ifndef NO_PTHREADS
sndobj *sndthread( int fifolen, sndobj *signal );
#endif
sndobj *skipwatch(const char *filename, sndobj *signal);
sndobj *panpot( float method, sndobj *right, sndobj *signal );
#define POS_REG_COMM	0x1234501
#define POS_REG_MASK	0x7F
#define POS_REG_SHIFT	1
#define POS_REGULAR(n)	((float*)(POS_REG_COMM|((n&POS_REG_MASK)<<POS_REG_SHIFT)))
sndobj *pansur(float *speakerpos, sndobj *sigpos, sndobj *sigrad, sndobj *signal);
sndobj *invA(sndobj *freq);
sndobj *sndassert(const char *cmp, const char *msg, sndobj *signal);


/*-----------------------------------------------------------------------------
    sndasync.c
-----------------------------------------------------------------------------*/

sndobj *linresample(sndobj *factor, sndobj *in);
sndobj *head( double interval, sndobj *in );
sndobj *repeat( double interval, sndobj *freq, sndobj *in );
sndobj *convolve( int klength, int advance, sndobj *kernel, sndobj *signal );
sndobj *interleave( int inigap, sndobj *gap, sndobj *signal );
sndobj *powersave( sndobj *gate, sndobj *in1, sndobj *in2 );
#define RESHAPE_SLOPE	    0
#define RESHAPE_KEEPAMPL    1
#define RESHAPE_KEEPANGLE   2
#define RESHAPE_KEEPRMS	    3
sndobj *reshape( double maxtime, int mode, sndobj *shapesig, sndobj *signal );
sndobj *repeathw(double maxlength, sndobj *repetitions, sndobj *grouplength, sndobj *signal);
sndobj *voiceact( float threshold, double gap, sndobj *cmp, sndobj *signal );
sndobj *wlinresample( double bottom, sndobj *factor, sndobj *signal );
sndobj *wasyncinterleave(double bottom, int ninputs, ... );
#define RESYNC_OUT	0
#define RESYNC_IN	1
#define RESYNC_TRIG	2
#define RESYNC_TIME	3
#define RESYNC_FREQ	4
sndobj *resync(int id, int mode, sndobj *signal, sndobj *factor);
void printresyncreg();


/*-----------------------------------------------------------------------------
    sndfilt.c
-----------------------------------------------------------------------------*/

sndobj *avg( sndobj *freq, sndobj *fade, sndobj *signal );
sndobj *gaussavg( int quality, sndobj *freq, sndobj *fade, sndobj *signal );
sndobj *rms( sndobj *freq, sndobj *signal );
extern float canonicalbands[];
sndobj *bands( float *tab, sndobj *fmult, sndobj *signal );
sndobj *graphiceq( float *tab, sndobj *fmult, sndobj *signal );
sndobj *filter( float *forwardfeed, float *backfeed, sndobj *signal );
sndobj *freqfir(float *response, int order, sndobj *signal);
sndobj *formant( sndobj *freq, sndobj *bandwidth, sndobj *signal );
sndobj *shelve( sndobj *freq, sndobj *gain, sndobj *signal );
sndobj *_allpass( int delay, float coeff, sndobj *signal );
sndobj *_foallpass( float coeff, sndobj *signal );
sndobj *_comb( int delay, float backfeed, sndobj *signal );
sndobj *_fcomb( int delay, float back2delay, float back2delayp1, sndobj *signal );
sndobj *_12_0pole( int iszero, float a0, float ab1, float ab2, sndobj *signal );
sndobj *lowpass( sndobj *cutoff, sndobj *signal );
sndobj *highpass( sndobj *cutoff, sndobj *signal );
sndobj *fastrms( sndobj *cutoff, sndobj *signal );
sndobj *median( int notmedianmode, sndobj *freq, sndobj *signal );
sndobj *dcblock( sndobj *in );
sndobj *peak( sndobj *sharpness, sndobj *freq, sndobj *signal );
sndobj *peaknotch(sndobj *freq, sndobj *gain, sndobj *bandwidth, sndobj *signal);
sndobj *boost(sndobj *freq, sndobj *gain, sndobj *bandwidth, sndobj *signal);
sndobj *rmpeak( float peakminheight, sndobj *signal );


/*-----------------------------------------------------------------------------
    sndeff.c
-----------------------------------------------------------------------------*/

sndobj *freeverb( sndobj *roomsize, sndobj *damping, sndobj *stereomix,
		    sndobj *revmix, sndobj *signal );
sndobj *_lbcf( int delay, float feedback, float lowcoeff, sndobj *signal );
sndobj *_fvap( int bufsize, float feedback, sndobj *signal );
sndobj *nrev( float length, float lowpass, float feedback, sndobj *signal );
sndobj *josrev( sndobj *lowpassf, sndobj *in );
sndobj *fdnrev( int *delays, double t60low, double t60high, sndobj *signal );
sndobj *boxrev(float *del_refl, float cutoff, float random, sndobj *signal);
sndobj *convrev( double revtime, double refltime, double directtime, double direct2noise, double refl2noise, sndobj *mixratio, sndobj *signal );
sndobj *revkernel( double revtime, double refltime, double directtime, double direct2noise, double refl2noise, int stereo );
sndobj *clip( float quench, sndobj *max, sndobj *signal );
sndobj *softclip( sndobj *max, sndobj *signal );
sndobj *crossoverdist( sndobj *ampl, sndobj *smooth, sndobj *signal );
sndobj *lee( sndobj *freq, sndobj *noiselvl, sndobj *signal );
sndobj *flanger( float coeff, float backfeed, float maxdelay, sndobj *del, sndobj *signal );
sndobj *phaser( int nfreq, float directgain, sndobj *basefreq, sndobj *freqfactor, sndobj *signal );
sndobj *leslie( sndobj *freq, sndobj *strength, sndobj *signal );
sndobj *chirp( float chirptime, sndobj *signal );
sndobj *spike( int nsam, sndobj *maxamp, sndobj *signal );
sndobj *fbmdistort(sndobj *lacunarity, sndobj *H, sndobj *cutoff, sndobj *source, sndobj *reset);
sndobj *avghw(int nmax, sndobj *n, sndobj *signal);


/*-----------------------------------------------------------------------------
    sndsource.c
-----------------------------------------------------------------------------*/

sndobj *harmonic( float *tab, sndobj *freq, sndobj *ampl );
sndobj *sine( float phase_2pi, sndobj *freq, sndobj *ampl );
sndobj *rect( float phase, sndobj *freq, sndobj *ampl, sndobj *ratio );
sndobj *saw( float phase, sndobj *freq, sndobj *ampl, sndobj *ratio );
sndobj *file(const char *name, sndobj *freq, ...);
sndobj *oneshot( const char *filename, sndobj *trigger, sndobj *speed );
sndobj *rapidfire( int nmax, const char *filename, sndobj *speed, sndobj *trigger );
#define MWTYPE_LINEAR   1
#define MWTYPE_WRAP     2
#define MWTYPE_FREQ     3
#define MWTYPE_MASK     0xFF
#define MWFLAG_CONTINUE 0x100
#define MWFLAGS_MASK    0xFFFFFF00
sndobj *multiwave( int dimension, const int *indextypes, 
	const char **ascfiles, const float *positions, sndobj *freq, ... );
sndobj *flatrand(float max);
sndobj *gaussrand(float sigma);
sndobj *fbm(sndobj *lacunarity, sndobj *H, sndobj *cutoff, sndobj *reset);
sndobj *linnoise( sndobj *slope, sndobj *nastyness );
sndobj *walk(float steprms);
sndobj *fwalk(sndobj *speed);
sndobj *slide( sndobj *minfreq, sndobj *maxfreq, sndobj *lacunarity, sndobj *exponent, sndobj *meandev );
#ifndef NO_OSS
sndobj *mike( int linein, float level );
#endif


/*-----------------------------------------------------------------------------
    sndsink.c
-----------------------------------------------------------------------------*/

sndobj *writeout(const char *filename, sndobj *in);
sndobj *print( FILE *stream, sndobj *in );
sndobj *tprint( FILE *stream, sndobj *in );
sndobj *sndstat( float min, float max, float binsize, sndobj *signal );
sndobj *beancounter( sndobj *signal );
#ifndef NO_OSS
sndobj *oss( float volume, sndobj *signal );
#endif
#ifndef NO_ALSA
sndobj *alsa(const char *device, sndobj *signal);
#endif
sndobj *meantime( sndobj *signal );


/*-----------------------------------------------------------------------------
    sndcron.c
-----------------------------------------------------------------------------*/

#define RAMP_L  0
#define RAMP_E  1
#define RAMP_C  2
sndobj *ramp( double offset, double *list );
sndobj *at( double offset, double *list );
sndobj *cron( double offset, double interval, float value, double *exceptions );
sndobj *exttrigger( float *data );
sndobj *expodecay( sndobj *decaylength, sndobj *trigger );
sndobj *idexp( sndobj *decaylength, sndobj *trigger );
sndobj *halfgauss( sndobj *halflength, sndobj *trigger );
sndobj *adsr( double attack, double decay, double sustain, double release, 
		float sustainfrac, sndobj *trigger );
sndobj *flipflop( sndobj *trigger, sndobj *input );
sndobj *holdtrigger( sndobj *trigger, sndobj *duration );
sndobj *holdnon0(int mode, sndobj *signal);
sndobj *timeout( double delay, sndobj *trigger, sndobj *signal );
#define SCHED_RR	0
sndobj *scheduler(int mode, int nclients, sndobj *trigger);


/*-----------------------------------------------------------------------------
    sndscore.c
-----------------------------------------------------------------------------*/

typedef struct {
  int       maxchord, jazzynotes, accentsinchord, enforcechordlength, silentfreq;
  float     mfvolume, fortefactor, transposefactor;
  double    len, stacclen, tenutolen;   //  relative to nominal length
  float     cue, stacccue, legatocue, legstacccue, endcue;
  char      *accentstring[PLENTY];  //  must include !...!
  float     accentcue[PLENTY];      //  < 0 means add to normal/stacc. ... cue
}
abcplayopts;

sndobj *abc( const char *filename, int voice, abcplayopts *opts );

sndobj *xyz(const char *filename, const char *voice);


/*-----------------------------------------------------------------------------
    sndmodel.c
-----------------------------------------------------------------------------*/

sndobj *karplusfilt(float loss, float blend, sndobj *signal);
sndobj *karplusbypass(double len, float ratio, double exponent, sndobj *signal);
sndobj *stringfeedback(float superinit, float si_loopscale, float si_initscale,
    sndobj *cue, sndobj *delay, sndobj *feed, sndobj *init);


/*-----------------------------------------------------------------------------
    sndwt.c
-----------------------------------------------------------------------------*/

sndobj *wthaarnn(double bottom, sndobj *signal);
sndobj *iwthaarnn(double bottom, sndobj *signal);
sndobj *wthaar(double bottom, sndobj *signal);
sndobj *iwthaar(double bottom, sndobj *signal);
sndobj *wtdaub( int n, double bottom, sndobj *signal );
sndobj *iwtdaub( int n, double bottom, sndobj *signal );
sndobj *wtdaubla( int n, double bottom, sndobj *signal );
sndobj *iwtdaubla( int n, double bottom, sndobj *signal );
sndobj *wtdaubbl( int n, double bottom, sndobj *signal );
sndobj *iwtdaubbl( int n, double bottom, sndobj *signal );
sndobj *wtcoif( int n, double bottom, sndobj *signal );
sndobj *iwtcoif( int n, double bottom, sndobj *signal );
sndobj *wtany( double bottom, float *scalingfunc, int llap, sndobj *signal );
sndobj *iwtany( double bottom, float *scalingfunc, int llap, sndobj *signal );
#define WDELAY_CENTRE	-1
sndobj *wdelay( double bottom, int n, int llap, sndobj *signal );
sndobj *iwdelay( double bottom, int n, int llap, sndobj *signal );
sndobj *wprint(FILE *out, double bottom, sndobj *signal);
sndobj *wtrunc(double bottom, sndobj *freq, sndobj *crossover, sndobj *signal);
sndobj *wthresh( sndobj *threshold, sndobj *wtsignal );
#define WSHIFT_SELECT	0
#define WSHIFT_SPREAD	1
#define WSHIFT_SUCC	2
sndobj *wshift( double bottom, int mode, int nshifts, sndobj *signal );
sndobj *wspread( double bottom, int mode, sndobj *spreadup, sndobj *spreaddown, sndobj *signal );
sndobj *wlinear(double bottom, int mode, int dim, float *matrix, sndobj *signal);
sndobj *wuninterleave(double bottom, sndobj *signal);
sndobj *winterleave(double bottom, sndobj *signal);
#define WSELECT_INDEX		0
#define WSELECT_SORTMAX		1
#define WSELECT_SORTMEAN	2
#define WSELECT_SORTMIN		3
#define WSELECT_SORTRMS		4
#define WSELECT_SAMERMS		0x100
sndobj *wselect( double bottom, sndobj *mode, sndobj *mask, sndobj *signal);
#define WFREQ_MODEMASK	0xFF
#define WFREQ_CONST	0
#define WFREQ_NOISE	1
#define WFREQ_ALIGNMASK	0x300
#define WFREQ_CENTRE	0x000
#define WFREQ_TOP	0x100
#define WFREQ_BOTTOM	0x200
sndobj *wfreq(double bottom, int mode, sndobj *freq, sndobj *mask);


#endif  // #ifndef __sndsys_h


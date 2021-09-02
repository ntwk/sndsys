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
#include <limits.h>
#include <errno.h>
#include <ctype.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/soundcard.h>
#include <unistd.h>
#include <signal.h>

#include "sndsys.h"


typedef struct _sharedbuf {
  void  *buf;
  int   usercount;
  char  *name;
  struct _sharedbuf *next;
}
sharedbuf;

struct loop {
  const char *loopname;
  int  loopid;
};

float complainundefined(sndobj *me, sndobj *caller, int innr);

static sndobj *firstinchain= NULL;

static sharedbuf *firstshared= NULL;

static sndobj nonexistent= { complainundefined, dontskip, {}, {}, {}, 0, 0, "nonexistent", "undefined", {}, NULL, NULL };

double starttime= 0.0, stoptime;

/*=============================================================================
    float complainundefined(sndobj *, sndobj *caller, int innr)

Outputs "undefined input" warning and returns 0.0.
=============================================================================*/

float complainundefined(sndobj *me, sndobj *caller, int innr)
{
  if( !caller )
    fprintf( stderr, "Warning: Undefined input used by unknown caller (NULL).\n" );
  else
    fprintf( stderr, "Warning: Undefined #%d input used by sndobj `%s', type `%s'.\n", innr, caller->name, caller->type );
  return 0.0;
}


/*=============================================================================
    void skipinputs(sndobj *me, sndobj *, int)

Default skip function.  Just skips all inputs.
=============================================================================*/

void skipinputs(sndobj *me, sndobj *caller, int innr)
{
  int in;
  
  for( in= 0; in< me->nin; ++in )
    INSKIP(in);
}


/*=============================================================================
    void dontskip(sndobj *me, sndobj *caller, int innr)

Skip function which calls the calc function instead of skipping.
=============================================================================*/

void dontskip(sndobj *me, sndobj *caller, int innr)
{
  me->calc(me, caller, innr);
}


/*=============================================================================
    sndobj *newsndo()

sndobj "constructor".
=============================================================================*/

sndobj *newsndo(float (*calc)(sndobj *,sndobj *, int), const char *name, 
		const char *type, int nch, int nin, ...)
{
  sndobj *p;
  va_list in;
  int i;
  
  if( nch <= 0 ) {
    fprintf(stderr, "Warning: newsndo: number of channels <= 0.  Name `%s', type `%s'.  Returning NULL.", name, type);
    return NULL;
  }
  if( nin < 0 )
    fprintf(stderr, "Warning: newsndo: number of inputs < 0.");
  p= (sndobj *)malloc(sizeof(sndobj));
  p->calc= calc;
  p->skip= skipinputs;
  p->nextinchain= firstinchain;
  firstinchain= p;
  if(nin>PLENTY)
    nin=PLENTY;
  p->nin= nin;
  p->nch= nch;
  p->name= name;
  p->type= type;
  p->exit= NULL;
  va_start(in, nin);
  for( i= 0; i< PLENTY; ++i ) {
    if( i<nin ) {
      p->in[i]= va_arg(in, sndobj*);
      p->inch[i]= p->in[i]->ch;
    }
    else {
      p->in[i]= &nonexistent;
      p->inch[i]= NULL;
    }
    p->ch[i]= 0.0;
    p->private[i]= NULL;
  }
  va_end(in);
  return p;
}


/*=============================================================================
    void addinput(sndobj *obj, sndobj *input)

Add another input to an object.  Returns 1 if successful, 0 if there were
already PLENTY of inputs.
=============================================================================*/

int addinput(sndobj *obj, sndobj *input)
{
  if( obj->nin >= PLENTY )
    return 0;
  obj->in[obj->nin]= input;
  obj->inch[obj->nin]= input->ch;
  ++obj->nin;
  return 1;
}


/*=============================================================================
    void dumpsndo(sndobj *p, int dumpwhat)

Print out information about the sound object.
=============================================================================*/

void dumpfloatarray( float *ar, int size );

void dumpsndo(sndobj *p, int dumpwhat)
{
  int in, previousin;
  
  if( !p ) {
    printf( "Dump of zero sndobj pointer attempted.\n");
    return;
  }
  printf("Dump of `%s', type `%s', address %p.\n", p->name, p->type, p );
  printf("  %d inputs, %d output channels.\n", p->nin, p->nch );
  if( dumpwhat & DUMP_INALL ) {
    previousin= PLENTY+1;
    for( in= 0; in< PLENTY; ++in )
      if( in < p->nin || p->in[in] != &nonexistent ) {
	if( previousin < in-1 )
	  printf("  ...\n");
	if( dumpwhat & DUMP_IN ) {
	  printf("  Input % 2d: ", in );
	  if( p->in[in] && (dumpwhat & DUMP_INNAME) )
	    printf("`%s', type `%s', address %p.\n", p->in[in]->name, 
				      p->in[in]->type, p->in[in] );
	  else
	    printf("%p\n", p->in[in] );
	}
	if( (dumpwhat & DUMP_INCH) )  {
	  printf("    Input channels #%d at %p", in, p->inch[in] );
	  if( p->inch[in] ) {
	    printf(": ");
	    dumpfloatarray( p->inch[in], PLENTY );
	    printf("\n");
	  }
	}
	if( p->in[in] && (dumpwhat & DUMP_INCH) && p->in[in]->ch!=p->inch[in] )
	{
	  printf("    Channels of input %d at %p: ", in, p->inch[in] );
	  dumpfloatarray( p->in[in]->ch, PLENTY );
	  printf("\n");
	}
	previousin= in;
      }
  }
  if( dumpwhat & DUMP_OUTCH ) {
    printf("  Output channels: ");
    dumpfloatarray( p->ch, PLENTY );
    printf("\n");
  }
  if( dumpwhat & DUMP_PRIV ) {
    printf("  Private pointers: ");
    for( in= 0; in< PLENTY-1; ++in )
      printf( "%p, ", p->private[in] );
    printf("%p\n", p->private[in] );
  }
}

void dumpfloatarray( float *ar, int size )
{
  float currval;
  int ind, samecount;

  for( ind= 0; ind< size; ) {
    currval= ar[ind];
    for( samecount= 1, ++ind; ind< size && ar[ind]==currval;  ++samecount, ++ind );
    if( samecount > 2 )
      printf("%g (%d times)", currval, samecount );
    else {
      if( samecount==2 )
	printf("%g, ", currval );
      printf("%g", currval );
    }
    if( ind<size )
	printf(", ");
  }
}


/*=============================================================================
    void disableskip()

Disables the skipping mechanism by replacing all skipping functions with 
dontskip.
=============================================================================*/

void disableskip()
{
  sndobj *scan;

  for( scan= firstinchain; scan; scan= scan->nextinchain )
    scan->skip= dontskip;
}


/*=============================================================================
    void closeloops(void)

This function removes the loop and defloop objects and wires the inputs of the
remaining objects so as to create the loops.  makefanouts() has to be called
after this function to take care of the resulting forks in the data stream.
=============================================================================*/

void closeloops(void)
{
  sndobj *def, *use, *next, *client;
  struct loop *defd, *used;
  int i;
  
  def= firstinchain;
  while( def ) {
    if( strcmp("defloop", def->type) ) {
      def= def->nextinchain;
      continue;
    }
    defd= (struct loop *)def->private[0];
    use= firstinchain;
    while( use ) {
      used= (struct loop *)use->private[0];
      if( use!=def && !strcmp("defloop", use->type) &&  
	    used->loopid == defd->loopid && 
	    !strcmp(used->loopname, defd->loopname) ) {
	fprintf( stderr, "Error: Source for loop `%s', id %d, doubly defined:"
	    " `%s', type `%s', address %p, and `%s', type `%s', address %p.\n",
	    defd->loopname, defd->loopid, def->in[0]->name, def->in[0]->type,
	    def->in[0], use->in[0]->name, use->in[0]->type, use->in[0] );
	exit(-1);
      }
      if( strcmp("loop", use->type) || used->loopid != defd->loopid ||
	    strcmp(used->loopname, defd->loopname) ) {
	use= use->nextinchain;
	continue;
      }
      if( use->nch != def->nch )
	fprintf(stderr, "Warning: Mismatching width in loop `%s', id %d: "
		"%d channels in source, %d channels used.\n", use->name, 
		*(int*)use->private[0], def->nch, use->nch );
      for( client= firstinchain; client; client= client->nextinchain )
	for( i= 0; i< client->nin; ++i )
	  if( client->in[i] == use ) {
	    if( client == def ) {
	      fprintf(stderr, "Error: Loop `%s', id %d, has no sndobj's except"
		  "loop and defloop.\n", defd->loopname, defd->loopid);
	      exit(-1);
	    }
	    client->in[i]= def->in[0];
	    client->inch[i]= def->inch[0];
	  }
      next= use->nextinchain;
      killsndo(use);
      use= next;
    }
    for( client= firstinchain; client; client= client->nextinchain )
      for( i= 0; i< client->nin; ++i )
	if( client->in[i] == def ) {
	  client->in[i]= def->in[0];
	  client->inch[i]= def->inch[0];
	}
    next= def->nextinchain;
    killsndo(def);
    def= next;
  }
  // Check for undefined loops left over:
  for( use= firstinchain; use; use= use->nextinchain )
    if( !strcmp("loop", use->type) ) {
      used= (struct loop *)use->private[0];
      fprintf(stderr, "Error: Loop `%s', id %d, without source.\n", 
			used->loopname, used->loopid);
      exit(-1);
    }
}


/*=============================================================================
    void prune(sndobj *top)

Remove objects unconnected to top-level object.
=============================================================================*/

void prune(sndobj *top)
{
  sndobj *p, *read, *write, **prevnext;
  int i;
  
  for( p= firstinchain, prevnext= &firstinchain; p; 
	    prevnext= &(p->nextinchain), p= p->nextinchain )
    if( p==top )
    {
      *prevnext= p->nextinchain;
      top->nextinchain= 0;
      break;
    }

  write= top;
  for( read= top; read; read= read->nextinchain )
    for( i= 0; i< read->nin; ++i )
      for( p= firstinchain, prevnext= &firstinchain; p; 
		prevnext= &(p->nextinchain), p= p->nextinchain )
	if( p==read->in[i] )
	{
	  *prevnext= p->nextinchain;
	  write->nextinchain= read->in[i];
	  write= write->nextinchain;
	  write->nextinchain= NULL;
	  break;
	}

  if( firstinchain )
  {
    fprintf( stderr, "Warning: Unconnected objects found.  Objects to be removed:\n");
    dumpchain();
    killemall();
  }
  firstinchain= top;
}


/*=============================================================================
    void makefanouts(void)

Insert fanout buffers whereever an object serves as several inputs to one or
several other objects.
=============================================================================*/

void fanout( int nclients, sndobj **clients, int *innrs );
// for implementation see below under buf()

void makefanouts(void)
{
  sndobj *p, *input, **clientlist;
  int *innrlist;
  int i, n, clientnr;
  
  for( input= firstinchain; input; input= input->nextinchain)
  {
    n= 0;
    for( p= firstinchain; p; p= p->nextinchain)
      for( i= 0; i< p->nin; ++i )
	if( input==p->in[i] )
	  ++n;
    if( n>1 && strcmp("fanout", input->type) )
    {
      if( !strcmp("buf", input->type) ) {
	fprintf( stderr, "Error: buffer `%s' fans out!\n", input->name );
	exit(-1);
      }
      clientlist= news(sndobj *, n);
      innrlist= news(int, n);
      clientnr= 0;
//    printf("Fanout for `%s', type `%s'\n", input->name, input->type );
      for( p= firstinchain; p; p= p->nextinchain)
	for( i= 0; i< p->nin; ++i )
	  if( input==p->in[i] ) {
	    clientlist[clientnr] = p;
	    innrlist[clientnr]= i;
//          printf("client %d: `%s', type `%s', input %d\n", clientnr, p->name, p->type, i );
	    ++clientnr;
	  }
      fanout(n, clientlist, innrlist);
      free(innrlist);
      free(clientlist);
    }
  }
}


/*=============================================================================
    void killsndo(sndobj *victim)

Destroy object and its privately allocated structures after calling its exit()
function.
=============================================================================*/

void killsndo(sndobj *victim)
{
  sndobj *p;
  int i;
  
  if( firstinchain==victim )
    firstinchain= victim->nextinchain;
  else
    for( p= firstinchain; p; p= p->nextinchain )
      if( p->nextinchain==victim ) {
	p->nextinchain= victim->nextinchain;
	break;
      }
  if( victim->exit )
    victim->exit(victim);
  for( i=0; i<PLENTY; ++i )
    if( !freeshared( victim->private[i] ) )
      free(victim->private[i]);
  free(victim);
}


/*=============================================================================
    void killemall(void)

Destroy all objects and their privately allocated structures.  The exit()
functions of all sndthread objects are called first to terminate all other
threads.  Then the exit() functions of all other objects are called in an
undefined order before the objects are destroyed.
=============================================================================*/

void killemall(void)
{
  sndobj *p, *nextp;
  int i;
  
  // First terminate threads by calling exit() functions of all sndthread
  // objects
  for( p= firstinchain; p; p= p->nextinchain )
    if( !strcmp(p->type, "sndthread") && p->exit )
      p->exit(p);
  for( p= firstinchain; p; p= p->nextinchain )
    if( strcmp(p->type, "sndthread") && p->exit )
      p->exit(p);
  p= firstinchain; 
  while( p )
  {
//fprintf( stderr, "Destroying object `%s', type `%s'.\n", p->name, p->type);
//fprintf( stderr, "freeing private pointers\n");
    for( i=0; i<PLENTY; ++i ) {
//    fprintf( stderr, "number %d: %p\n", i, p->private[i] );
      if( !freeshared( p->private[i] ) )
	free(p->private[i]);
    }
    nextp= p->nextinchain;
//fprintf( stderr, "freeing sndobj\n");
    free(p);
    p= nextp;
  }
  firstinchain= NULL;
}


/*=============================================================================
    void sndexecute(double t_start, double t_end, sndobj *top)

Loops over evaluation of the object top, throwing the results away.  t_start
and t_end are assigned to the global variables starttime and stoptime.  Their
difference gives the number of samples to evaluate, while the absolute
starttime can influence score and scheduling objects which use it to skip
events and cues if it is > 0.  Likewise, stoptime can be reset by an object to
stop evaluation prematurely; stoptime = -DBL_MAX will cause sndexecute to abort
immediately.

Before the start of the evalutation, closeloops(), prune() and makefanouts()
are called, saving the caller these chores.  The computation time is noted, and
after the evaluation is completed, the time, number of samples and samples
computed per second is output.  Since the wallclock time, not cpu time is used,
this serves only for illustration.  Before returning, all objects are deleted
with killemall().
=============================================================================*/

void sndexecute(double t_start, double t_end, sndobj *top)
{
  struct timeval calct;
  double sndtime, sampersec;
  long n, calcs, calcus;

  closeloops();
  prune(top);
  makefanouts();
// dumptree(top);
//  printresyncreg();
  starttime= t_start;
  stoptime= t_end;
  gettimeofday(&calct, NULL);
  calcs= calct.tv_sec;
  calcus= calct.tv_usec;
  for( sndtime= t_start + 0.5/(double)SAMPRATE, n= 0; sndtime< stoptime; 
       sndtime += 1.0/(double)SAMPRATE, ++n )
    top->calc(top, NULL, 0);
  gettimeofday(&calct, NULL);
  calcs= calct.tv_sec - calcs;
  calcus= calct.tv_usec - calcus;
  if( calcus<0 ) {
    --calcs;
    calcus += 1000000;
  }
  sampersec= floor((double)n/((double)calcs+1e-6*(double)calcus) + 0.5);
  printf("Computation finished.  Took");
  if( calcs>=3600 ) {
    printf(" %ld hour%s", calcs/3600, calcs >= 7200? "s" : "");
    calcs %= 3600;
    if( calcs || calcus ) {
      printf(" %ld minute%s", calcs/60, calcs >= 120? "s" : "" );
    }
  }
  else if( calcs>=60 ) {
    printf(" %ld minute%s", calcs/60, calcs >= 120? "s" : "" );
    calcs %= 60;
    if( calcs )
      printf(" %ld second%s", calcs, calcs > 1? "s" : "" );
  }
  else
    printf(" %.3g seconds", (double)(calcs%60) + 1e-6*calcus );
  printf(" for %ld samples (%g seconds); %ld samples per second; %.2g times audio speed.\n", n, (double)n/(double)SAMPRATE, (long)sampersec, sampersec/SAMPRATE);
  killemall();
}


/*=============================================================================
    void dumpchain(void)
    void dumptree(sndobj *top)
=============================================================================*/

void dumpchain(void)
{
  sndobj *p;
  
  for( p= firstinchain; p; p= p->nextinchain )
    printf( "`%s' (%s) :%d\n", p->name, p->type, p->nch );
}


void dumptree(sndobj *top)
{
  sndobj *p, **track;
  short *path;
  int n, i, depth, rdepth;
  
  for( p= firstinchain, n= 0; p; p= p->nextinchain, ++n );
  track= (sndobj **)malloc((n+1) * sizeof(sndobj*));
  path= (short *)malloc((n+1) * sizeof(short));
  depth= 0;
  track[0]= top;
  path[0]= 0;
  while( depth>=0 )
  {
    for( i= 0; i<depth; ++i )
      printf( "  " );
    p= track[depth];
    for( rdepth= depth-1; rdepth>=0; --rdepth )
      if( p==track[rdepth] )
	break;
    if( rdepth >= 0 ) {
      printf( "`%s' (%s) :%d (recursive %d-loop)\n", p->name, p->type, p->nch, depth-rdepth );
      --depth;
    }
    else {
      printf( "`%s' (%s) :%d\n", p->name, p->type, p->nch );
    }
    while( depth>=0 && path[depth] >= track[depth]->nin )
      --depth;
    if( depth>=0 )
    {
      track[depth+1]= track[depth]->in[path[depth]++];
      path[depth+1]= 0;
      ++depth;
    }
  }
  free(path);
  free(track);
}


/*=============================================================================
    defloop "name" (id) <source>
    loop "name" (id) (nch)

Define or reference a loop.  The loop is identified with the "name" and (id)
you pass to both functions.  Every loop can have several references (loop
objects) but only one definition (defloop object).  <source> is the sndobj the
data from which will be looped around.  (nch) should be set to the number of
channels of <source>.  It is used only for providing a channel number to
objects of which loop() is an input and which want to know about their inputs'
number of channels.

    sndobj *defloop( const char *name, int id, sndobj *source )
    sndobj *loop( const char *name, int id, int nch )
    float calcloop( sndobj *me, sndobj *, int)
=============================================================================*/

sndobj *loop( const char *name, int id, int nch )
{
  sndobj *p;
  struct loop *d;

  p= newsndo( calcloop, name, "loop", nch, 0 );
  p->private[0]= d= new(struct loop);
  d->loopname= name;
  d->loopid= id;
  return p;
}

sndobj *defloop( const char *name, int id, sndobj *source )
{
  sndobj *p;
  struct loop *d;

  p= newsndo( calcloop, name, "defloop", source->nch, 1, source );
  p->private[0]= d= new(struct loop);
  d->loopname= name;
  d->loopid= id;
  return p;
}

float calcloop( sndobj *me, sndobj *caller, int innr)
{
  fprintf(stderr, "Error: %s object `%s', id %d, called by `%s', type `%s' "
	"(input %d).  Use closeloops() to close loops before running!\n", 
	me->type, me->name, *(int*)me->private[0], 
	caller->name, caller->type, innr );
  exit(-1);
  return 0.0;
}


/*=============================================================================
    void buf( float , float , float , int nch, sndobj *client, int innr )
    float calcbuf( sndobj *me, sndobj *caller, int innr )  
    float intercalcbuf( sndobj *me, sndobj *caller, int innr, int index, float frac )  
    void updatebase( sndobj *me, struct bufread *rincr, int skip )
    void skipbuf( sndobj *me, sndobj *caller, int innr)
    void incrbuf( sndobj *me, sndobj *caller, int innr, int inc )
    void getbufptr( sndobj *me, sndobj *caller, int innr, int index, float **ptr, int *nsamples, int *nch )
    void fanout( int nclients, sndobj **clients, int *innrs )
    void set_loopclients( sndobj *fanout )
=============================================================================*/

#define BUF_MINSIZE     256
struct bufwrite {
  float     *buf, *base, *top;
  long      size, writeoff, maxoff;
  int       nch, haveloop, writelock;
  char      name[NAMESTRLEN];
};
struct bufread {
  sndobj    *client;
  long      minoff, off0, maxoff;
  float     ch[PLENTY];
  int       innr, bufclient, loopclient;
};

float calcbuf( sndobj *me, sndobj *caller, int innr );
void updatebase( sndobj *me, struct bufread *rincr, int skip );
void skipbuf( sndobj *me, sndobj *caller, int innr );
void set_loopclients( sndobj *fanout );

void buf( float backtime, float forwtime, float jitter, int nch, sndobj *client, int innr )
{
  sndobj *signal, *p;
  struct bufwrite *w;
  struct bufread *r;
  
  if( forwtime < 0.0 )
    forwtime= 0.0;
  if( backtime < 0.0 )
    backtime= 0.0;
  if( !forwtime && !backtime )
    printf( "Warning: buffer for zero time interval created by `%s', type %s for input %d.\n", client->name, client->type, innr);
  signal= client->in[innr];
  w= new(struct bufwrite);
  w->writelock= 0;
  w->nch= nch;
  snprintf( w->name, NAMESTRLEN, "+%g/-%g [%s]", forwtime+jitter, backtime+jitter, signal->name );
  w->maxoff= ((int)ceil(forwtime*SAMPRATE) + (int)ceil(backtime*SAMPRATE) + 
		2*(int)ceil(jitter*SAMPRATE)) * nch;
  w->size= w->maxoff + nch;
  if( w->size < BUF_MINSIZE*nch )
    w->size= BUF_MINSIZE*nch;           // avoid too many wraparounds
  p= newsndo( calcbuf, w->name, "buf", nch, 1, signal );
  p->skip= skipbuf;
  p->private[0]= w;
  p->private[1]= w->buf= w->base= (float*)calloc( w->size, sizeof(float) );
  w->top= w->buf + w->size;
  p->private[2]= r= news(struct bufread, 2);
  r->minoff= 0;
  w->writeoff= r->off0= nch * (int)ceil(backtime*SAMPRATE);
  r->maxoff= r->off0 + nch * (int)ceil(forwtime*SAMPRATE);
  r->client= client;
  r->innr= innr;
  r->bufclient= 1;
  r->loopclient= 0;
  r[1].client= NULL;
  client->in[innr]= p;
  client->inch[innr]= p->ch;   //  Just by way of initialisation
//printf("created buf for `%s', type `%s', client `%s', type `%s', %d channels, with wmaxoff %ld, forwtime %g, rmaxoff %ld\n", signal->name, signal->type, client->name, client->type, w->nch, w->maxoff, forwtime, r->maxoff );
}


float calcbuf( sndobj *me, sndobj *caller, int innr )
{
  return intercalcbuf(me, caller, innr, 0, 0.0);
}


float intercalcbuf( sndobj *me, sndobj *caller, int innr, int index, float frac )  
// frac must be in [0,1]
{
  struct bufwrite *w;
  struct bufread *r;
  float *write, *read, *prev;
  long offset;
  int calcvals, ch;

//fprintf(stderr, "intercalcbuf for `%s', type `%s', client `%s', type `%s', input %d, index %d\n", me->in[0]->name, me->in[0]->type, caller->name, caller->type, innr, index);
  w= (struct bufwrite *)me->private[0];
  r= (struct bufread *)me->private[2];
  while( r->client && (r->client!=caller || r->innr!=innr) )
    ++r;
  offset= r->off0 + index * w->nch;
//fprintf(stderr, "off %ld, rmax %ld, wmax %ld, size %ld\n", offset, r->maxoff, w->maxoff, w->size);
  if( offset < r->minoff || offset < 0 ) {
    caller->inch[innr]= me->ch;  // zero
    return 0.0;
  }
  if( frac )
    offset += w->nch;
  if( offset > r->maxoff || offset > w->maxoff ) {
    caller->inch[innr]= me->ch;   // zero
    return 0.0;
  }
  if( offset >= w->writeoff && !w->writelock ) {
    w->writelock= 1;
    write= w->base + w->writeoff;
    if( write >= w->top )
      write -= w->size;
    for( calcvals= (int)(offset-w->writeoff+1); calcvals> 0; calcvals-=w->nch )
    {
      *write++ = INCALC(0);
      for( ch= 1; ch < w->nch; ++ch )
	*write++ = INPUTC(0, ch);
      if( write >= w->top )
	write -= w->size;
    }
    w->writeoff= offset + w->nch;
    w->writelock= 0;
  }
  read= w->base + offset;
  if( read >= w->top )
    read -= w->size;
  if( frac ) {
    prev= w->base + offset - 1;
    if( prev >= w->top )
      prev -= w->size;
    for( ch= 0; ch< me->nch; ++ch )
      r->ch[ch]= prev[ch] + frac*(read[ch]-prev[ch]);
    caller->inch[innr]= r->ch;
    return r->ch[0];
    //  No increment is ever necessary for buffered clients, so we return here
  }
  caller->inch[innr]= read;
  if( !r->bufclient ) {
    r->off0 += w->nch;
    r->minoff += w->nch;
    r->maxoff += w->nch;
    updatebase(me, r, 0);
  }
  return *read;
}


// This function serves to call the input's skip function when it is certain
// that past values cannot be needed any more by any of the clients.

void updatebase( sndobj *me, struct bufread *rincr, int skip )
{
  struct bufwrite *w;
  struct bufread *rlist, *r;
  float *write;
  long tailminoff;
  int ch;

//fprintf(stderr, "updatebase for `%s', type `%s', client `%s', type `%s', input %d\n", me->in[0]->name, me->in[0]->type, rincr->client->name, rincr->client->type, rincr->innr);
  w= (struct bufwrite *)me->private[0];
  if( w->writelock )
    return;
  rlist= (struct bufread *)me->private[2];
  tailminoff= rincr->minoff;
  if( tailminoff <= 0 )
    return;
//  fprintf(stderr, "min offsets: ");
//  for( r= rlist; r->client; ++r ) {
//    if( r->loopclient )
//	  fprintf(stderr, "L");
//    if( r==rincr )
//      fprintf(stderr, "*");
//    fprintf(stderr, "%ld, ", r->minoff );
//  }
//  fprintf(stderr, "\n");
  for( r= rlist; r->client; ++r )
    if( r->minoff < tailminoff ) {
//  rincr is the one whose off0 (and minoff) field has been incremented.  So we
//  know that if any client lags behind rincr and we have no loops, no updating
//  is necessary.
      if( !w->haveloop )
	return;
      tailminoff= r->minoff;
    }

  if( w->writeoff < tailminoff )
  {
    if( skip && !w->writelock )
    {
      long skipvals= tailminoff - w->writeoff;
      w->writelock= 1;
      for( ; skipvals> 0;  skipvals -= w->nch )
	INSKIP(0);
      w->writelock= 0;
    }
    else
      tailminoff= w->writeoff;
    w->writeoff= 0;
  }
  else
    w->writeoff -= tailminoff;
  w->base += tailminoff;
  if( w->base >= w->top )
    w->base -= w->size;
  for( r= rlist; r->client; ++r ) {
    r->minoff -= tailminoff;
    r->off0 -= tailminoff;
    r->maxoff -= tailminoff;
  }
}


void skipbuf( sndobj *me, sndobj *caller, int innr )
{
  struct bufwrite *w;
  struct bufread *r;
  
  w= (struct bufwrite *)me->private[0];
  r= (struct bufread *)me->private[2];
  while( r->client && (r->client!=caller || r->innr!=innr) )
    ++r;
  r->off0 += w->nch;
  r->minoff += w->nch;
  r->maxoff += w->nch;
  updatebase(me, r, 1);
}


void incrbuf( sndobj *me, sndobj *caller, int innr, int incr )
{
  struct bufwrite *w;
  struct bufread *r;
  
  if( incr <= 0 ) {
    if( incr )
      fprintf(stderr, "incrbuf: Warning: `%s\', type `%s\', attempted to decrement buffer index - ignored.\n", caller->name, caller->type );
    return;
  }
  w= (struct bufwrite *)me->private[0];
  r= (struct bufread *)me->private[2];
  while( r->client && (r->client!=caller || r->innr!=innr) )
    ++r;
  if( !r->client ) {
    fprintf(stderr, "incrbuf: Warning: called by sndobj `%s\', type `%s\', which is not a client - ignored.\n", caller->name, caller->type );
    return;
  }
  r->off0 += incr * w->nch;
  r->minoff += incr * w->nch;
  r->maxoff += incr * w->nch;
  if( w->haveloop ) {
    intercalcbuf(me, caller, innr, r->minoff-r->off0, 0.0f);
    updatebase(me, r, 0);
  }
  else
    updatebase(me, r, 1);
}


void getbufptr( sndobj *me, sndobj *caller, int innr, int index, float **ptr, int *nsamples, int *nch )
{
  struct bufwrite *w;
  struct bufread *r;
  float *result;
  long offset;

  w= (struct bufwrite *)me->private[0];
  r= (struct bufread *)me->private[2];
  while( r->client && (r->client!=caller || r->innr!=innr) )
    ++r;
  offset= r->off0 + index * w->nch;
  if( offset < r->minoff || offset < 0 || offset > r->maxoff || offset >= w->writeoff ) {
    *ptr= NULL;
    *nsamples= 0;
    *nch= w->nch;
    return;
  }
  result= w->base + offset;
  if( result >= w->top )
    result -= w->size;
  *ptr= result;
  if( result >= w->base )
    *nsamples= (w->top - result) / w->nch;
  else
    *nsamples= (w->writeoff - offset) / w->nch;
  *nch= w->nch;
}


void fanout( int nclients, sndobj **clients, int *innrs )
{
  sndobj *signal, *p;
  struct bufwrite *w, *otherw;
  struct bufread *r, *otherr;
  long maxbackind, maxforwind;
  int i;
  
  signal= (*clients)->in[*innrs];
  if( !strcmp(signal->type, "buf") )
    signal= signal->in[0];

  //  Now we determine the needed buffer size from the fields writeoff and
  //  maxoff of the bufwrite structures, which unlike the information in the
  //  bufread structs include the jitter.  Non-buffered clients are allowed up
  //  to HORIZON jitter, but no additional forward or backward tolerance.
  maxforwind= maxbackind= (int)ceil(HORIZON*SAMPRATE);
  for( i= 0; i< nclients; ++i )
    if( !strcmp( clients[i]->type, "buf" ) ) {
      otherw= (struct bufwrite *)clients[i]->private[0];
      if( maxbackind < otherw->writeoff/otherw->nch )
        maxbackind= otherw->writeoff/otherw->nch;
      if( maxforwind < (otherw->maxoff-otherw->writeoff)/otherw->nch )
        maxforwind= (otherw->maxoff-otherw->writeoff)/otherw->nch;
    }
  w= new(struct bufwrite);
  w->writelock= 0;
  w->haveloop= 0;
  w->nch= signal->nch;
  w->writeoff= maxbackind * w->nch;
  w->maxoff= (maxbackind + maxforwind + 1) * w->nch;
  // Additional 1 sample tolerance since clients execute in turn and therefore
  // may be one sample out of sync.
  w->size= w->maxoff + w->nch;
  // Additional space for one sample since one value has to be stored even if
  // maxbackind and maxforwind are 0
  if( w->size < BUF_MINSIZE*w->nch )
    w->size= BUF_MINSIZE*w->nch;           // avoid too many wraparounds
  snprintf( w->name, NAMESTRLEN, "+%g/-%g [%s]", (double)maxforwind/(double)SAMPRATE, (double)maxbackind/(double)SAMPRATE, signal->name );
  p= newsndo( calcbuf, w->name, "fanout", signal->nch, 1, signal );
  p->skip= skipbuf;
  p->private[0]= w;
  p->private[1]= w->buf= w->base= (float*)calloc( w->size, sizeof(float) );
  w->top= w->buf + w->size;
  p->private[2]= r= news(struct bufread, nclients+1);
  for( i= 0; i< nclients; ++i )
    if( !strcmp( clients[i]->type, "buf" ) ) {
      otherw= (struct bufwrite *)clients[i]->private[0];
      otherr= (struct bufread *)clients[i]->private[2];
      r[i].client= otherr->client;
      r[i].innr= otherr->innr;
      r[i].bufclient= 1;
      r[i].loopclient= 0;
      r[i].off0= w->writeoff;
      r[i].minoff= r[i].off0 - 
		    (otherr->off0 - otherr->minoff) / otherw->nch * w->nch;
      r[i].maxoff= r[i].off0 +
		    (otherr->maxoff - otherr->off0) / otherw->nch * w->nch;
      r[i].client->in[r[i].innr]= p;
      r[i].client->inch[r[i].innr]= p->ch;   //  Just by way of initialisation
      killsndo( clients[i] );
    }
    else {
      r[i].client= clients[i];
      r[i].innr= innrs[i];
      r[i].bufclient= 0;
      r[i].loopclient= 0;
      r[i].off0= w->writeoff;
      r[i].minoff= r[i].off0;
      r[i].maxoff= r[i].off0;
      clients[i]->in[innrs[i]]= p;
      clients[i]->inch[innrs[i]]= p->ch;   //  Just by way of initialisation
    }
  
  r[i].client= NULL;
  set_loopclients(p);
  if( w->haveloop )
    p->skip= dontskip;
//printf("created fanout for `%s', type `%s', %d channels, with wmaxoff %ld, rmaxoff %ld, maxbackind %ld, maxforwind %ld\n", signal->name, signal->type, w->nch, w->maxoff, r->maxoff, maxbackind, maxforwind );
}


void set_loopclients( sndobj *fanout )
{
  sndobj *p, **track;
  struct bufwrite *w;
  struct bufread *r;
  short *path;
  int n, depth, rdepth;
  
  for( p= firstinchain, n= 0; p; p= p->nextinchain, ++n );
  track= (sndobj **)malloc((n+1) * sizeof(sndobj*));
  path= (short *)malloc((n+1) * sizeof(short));
  w= (struct bufwrite *)fanout->private[0];
  w->haveloop= 0;
  depth= 0;
  track[0]= fanout;
  path[0]= 0;
  while( depth>=0 )
  {
    p= track[depth];
    for( rdepth= depth-1; rdepth>=0; --rdepth )
      if( p==track[rdepth] )   // search for loops, whether or not they contain
	break;		       // our fanout
    if( rdepth >= 0 )
      --depth;	  // backtrack so we don't go round in circles in the while loop
    if( !rdepth ) { // loop containing fanout -> previous sndobj is loop client
      r= (struct bufread *)fanout->private[2];
      while( r->client && (r->client!=track[depth] || r->innr!=path[depth]-1) )
	++r;
      r->loopclient= 1;
      w->haveloop= 1;
    }
    while( depth>=0 && path[depth] >= track[depth]->nin )
      --depth;
    if( depth>=0 )
    {
      track[depth+1]= track[depth]->in[path[depth]++];
      path[depth+1]= 0;
      ++depth;
    }
  }
  free(path);
  free(track);
}


/*=============================================================================
    void registershared( void *buf, char *name )

Registers a shared buffer.  Both buf and name have to be dynamically allocated.
They are taken care of by the shared buffer functions and should not be freed
by hand.  buf has to be added to the private pointer list of the calling
sndobj.
=============================================================================*/

void registershared( void *buf, char *name )
{
  sharedbuf *end, *newshared, **prevnext;
  
//fprintf(stderr, "registering %s\n", name );
  prevnext= &firstshared;
  for( end= firstshared; end; prevnext= &end->next, end= end->next );
  *prevnext= newshared= new(sharedbuf);
  newshared->buf= buf;
  newshared->usercount= 1;
  newshared->name= name;
  newshared->next= NULL;
}


/*=============================================================================
    void *findshared( const char *name, int registerme )

Finds a shared buffer by name and returns a pointer to the buffer (not to the
sharedbuf structure).  If the name is not found, NULL is returned.  If the name
is found and registerme is set, the user count is incremented.  The calling
sndobj must then add the buffer pointer to its private pointer list.
=============================================================================*/

void *findshared( const char *name, int registerme )
{
  sharedbuf *search;
  
//fprintf(stderr, "finding shared %s\n", name );
  for( search= firstshared; search; search= search->next )
    if( !strcmp(search->name, name) ) {
      if( registerme )
	++search->usercount;
//    fprintf(stderr, "found\n");
      return search->buf;
    }
//fprintf(stderr, "not found\n");
  return NULL;
}


/*=============================================================================
    int freeshared( void *buf )

Decrements the user count in the shared buffer corresponding to buf and frees
the buffer, name and sharedbuf structure if this was the last user.  Returns 1
if successfull, 0 if the buffer was not found.
=============================================================================*/

int freeshared( void *buf )
{
  sharedbuf *search, **prevnext;
  
  if( !buf )
    return 0;
  prevnext= &firstshared;
  for( search= firstshared; search; prevnext= &search->next, search= search->next )
    if( search->buf==buf ) {
//    fprintf(stderr, "freeing shared %s\n", search->name );
      if( !--search->usercount ) {
	*prevnext= search->next;
//      fprintf(stderr, "last user deregistered, really freeing now\n");
	free( search->buf );
	free( search->name );
	free( search );
      }
      return 1;
    }
  return 0;
}


/*=============================================================================
    float *sinelookup(long size)

Creates a sine lookup table of the specified length.  The table has to be freed
by the caller.
=============================================================================*/

float *sinelookup( long size )
{
  float *write, *tab;
  double phase, inc;
  long i;
  
  tab= malloc( size * sizeof(float) );
  inc= 2.0*M_PI/(double)size;
  phase= 0.0;
  write= tab;
  for( i= size; i> 0; --i, phase += inc )
    *write++ = (float)sin(phase);
  return tab;
}


/*=============================================================================
    rngdesc *makerng(void)
    float flatrng( rngdesc *rng )       equally distributed in [-1,1]
    float gaussrng( rngdesc *rng )      mean= 0, sigma= 1
=============================================================================*/

void initknuthrng( knuthrngstate *rngdata );
float knuthrand( knuthrngstate *rngdata );

rngdesc *makerng(void)
{
  static unsigned callcount= 0;

  rngdesc *rng;
  unsigned seed;

  // Both time() and clock() can return the same value if multiple RNGs are 
  // created without delay.  So we add the number of times this has been
  // called.
  seed= (time(NULL)^clock()) + callcount;
  ++callcount;
  rng= new(rngdesc);
  rng->re_buf.state= NULL;
  if( initstate_r( seed, rng->state, RNG_ARRAYSIZE, &rng->re_buf ) )
    fprintf(stderr, "Warning: Random number generator initialisation: initstate_r returned error.\n" );
  if( setstate_r(rng->state, &rng->re_buf) )
    fprintf(stderr, "Warning: Random number generator initialisation: setstate_r returned error.\n" );
  return rng;
}

float flatrng( rngdesc *rng )
{
  int intrand;

  random_r(&rng->re_buf, &intrand);
  return 2.0f*(float)intrand/(float)RAND_MAX - 1.0f;
}

float gaussrng( rngdesc *rng )
{
  float r;
  int intrand;

  random_r(&rng->re_buf, &intrand);
  r = (float)intrand;
  random_r(&rng->re_buf, &intrand);
  r += (float)intrand;
  random_r(&rng->re_buf, &intrand);
  r += (float)intrand;
  random_r(&rng->re_buf, &intrand);
  r += (float)intrand;
  random_r(&rng->re_buf, &intrand);
  r += (float)intrand;
  /* The sum of n equally-distributed random variables tends towards a 
     gaussian distribution for large n; here: n=5  */
  r = 1.54919333848296675407 * r / (float)RAND_MAX - 3.87298334620741688517;
  /* normalisation to mean=0, sigma=1: sqrt(12/5) * r / RAND_MAX - sqrt(3*n) */
  return r;
}

#define KRNG_MOD  1000000000

void initknuthrng( knuthrngstate *rngdata )
{
  static unsigned callcount= 0;

  int initval, previnitval, ind, count;

  // Both time() and clock() can return the same value if multiple RNGs are 
  // created without delay.  So we add the number of times this has been
  // called.
  previnitval= abs(161803398 - ((time(NULL)^clock()) + callcount++));
  previnitval %= KRNG_MOD;
  rngdata->state[55]= previnitval;
  initval= 1;
  for( count= 1; count< 55; ++count ) {
    ind= (21*ind) % 55;
    rngdata->state[ind]= initval;
    initval= previnitval - initval;
    if( initval< 0 )
      initval += KRNG_MOD;
    previnitval= rngdata->state[ind];
  }
  for( count= 0; count< 4; ++count )
    for( ind= 1; ind< 56; ++ind ) {
      rngdata->state[ind] -= rngdata->state[1+(ind+30)%55];
      if( rngdata->state[ind] < 0 )
	rngdata->state[ind] += KRNG_MOD;
    }
  rngdata->inext= 0;
  rngdata->inextp= 31;
}

float knuthrand( knuthrngstate *rngdata )
{
  int difference;

  if( ++rngdata->inext >= 56 )
    rngdata->inext= 1;
  if( ++rngdata->inextp >= 56 )
    rngdata->inextp= 1;
  difference= rngdata->state[rngdata->inext] - rngdata->state[rngdata->inextp];
  if( difference < 0 )
    difference += KRNG_MOD;
  rngdata->state[rngdata->inext]= difference;
  return (float)difference*(2.0f/KRNG_MOD) - 1.0f;
}




/*=============================================================================
    float lininterpol( float *buf, long size, float pos )
=============================================================================*/

float lininterpol( float *buf, long size, float fpos )
{
  float *read, *next;
  
  read= buf + (long)floorf(fpos);
  while( read >= buf + size )
    read -= size;
  while( read < buf )
    read += size;
  fpos -= floorf(fpos);
  next= read + 1;
  if( next >= buf+size )
    next -= size;
  return *read + fpos * (*next - *read);
}



/*=============================================================================
    int getprime( int lowlimit )

Returns prime >= lowlimit.
=============================================================================*/

static int primetable[]= { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 0 };

int getprime( int lowlimit )
{
  int *tab;
  int num, div, divlimit;
  
  for( tab= primetable; *tab; ++tab )
    if( *tab>=lowlimit )
      return *tab;
  for( num= lowlimit | 1; 13; num += 2 )
  {
    divlimit= (int)floor(sqrt(num));
    for( tab= primetable; *tab && *tab<divlimit && num%*tab!=0; ++tab );
    if( *tab ) {
      if( *tab > divlimit )
	return num;
      if( num%*tab==0 )
	continue;
    }
    for( div= tab[-1]; div <= divlimit && num%div!=0; div += 2 );
    if( div>divlimit )
      return num;
  }
}


/*=============================================================================
    long binfindf( float *sorted, long size, float val )

Finds the position in a sorted array (in ascending order) before which a new
value must be inserted.  The return value is in [0, size].  If the array
already contains a value equal to val, its position is returned (if there are
several, the exact position is undefined).  size must be at least 1.
=============================================================================*/

long binfindf( float *sorted, long size, float val )
{
  float *search, *left, *right;
  
  left= sorted;
  if( val <= *left )
    return 0;
  right= sorted+size-1;
  if( val > *right )
    return size;
  else if( val == *right )
    return size-1;
  search= sorted + size/2;
  while( right-left > 1 )
    if( val > *search ) {
      left= search;
      search += (right-left+1)/2;
    }
    else if( val < *search ) {
      right= search;
      search -= (right-left+1)/2;
    }
    else
      return search-sorted;
  return right-sorted;
}


/*=============================================================================
    int dspdesc()

Returns the file descriptor of the DSP device after opening it if necessary.
Unavailable if NO_OSS has been defined in sndsys.h.
=============================================================================*/

#ifndef NO_OSS

void dspatexit(void);
void dsponinterrupt(int signr);

static void (*dsp_oldinthandler)(int);

int dspdesc(int restore)
{
  static struct {
    int valid, dspd, old_pcm, old_recsrc, old_mic, old_line, old_igain;
  } dsp = { 0 };
  int ioarg;
  signed short testend;
  
  if( restore )
  {
    if( !dsp.valid )
      return 0;
    ioctl( dsp.dspd, SNDCTL_DSP_SYNC );     // flush, wait and reset
    ioctl( dsp.dspd, SOUND_MIXER_WRITE_PCM, &dsp.old_pcm );
    ioctl( dsp.dspd, SOUND_MIXER_WRITE_RECSRC, &dsp.old_recsrc );
    ioctl( dsp.dspd, SOUND_MIXER_WRITE_MIC, &dsp.old_mic );
    ioctl( dsp.dspd, SOUND_MIXER_WRITE_LINE, &dsp.old_line );
    ioctl( dsp.dspd, SOUND_MIXER_WRITE_IGAIN, &dsp.old_igain );
    dsp.valid= 0;   // no need to do it twice
    return 0;
  }
  else if( !dsp.valid ) 
  {
    dsp.dspd= open("/dev/dsp", O_RDWR );
    ioctl( dsp.dspd, SNDCTL_DSP_RESET );
    ioarg= APF_CPUINTENS;
    ioctl( dsp.dspd, SNDCTL_DSP_PROFILE, &ioarg );
    testend= 256;
    if( *(char*)&testend )
      ioarg= AFMT_S16_BE;
    else
      ioarg= AFMT_S16_LE;
    ioctl( dsp.dspd, SNDCTL_DSP_SETFMT, &ioarg );
    ioarg= 2;
    ioctl( dsp.dspd, SNDCTL_DSP_STEREO, &ioarg );
    ioarg= SAMPRATE;
    ioctl( dsp.dspd, SNDCTL_DSP_SPEED, &ioarg );
    if( ioarg != SAMPRATE )
      fprintf(stderr, "Warning: OSS DSP sampling rate is %d not %d.\n", ioarg, SAMPRATE );
    ioctl( dsp.dspd, SOUND_MIXER_READ_PCM, &dsp.old_pcm );
    ioctl( dsp.dspd, SOUND_MIXER_READ_RECSRC, &dsp.old_recsrc );
    ioctl( dsp.dspd, SOUND_MIXER_READ_MIC, &dsp.old_mic );
    ioctl( dsp.dspd, SOUND_MIXER_READ_LINE, &dsp.old_line );
    ioctl( dsp.dspd, SOUND_MIXER_READ_IGAIN, &dsp.old_igain );
    dsp.valid= 1;
    atexit( dspatexit );
    dsp_oldinthandler= signal( SIGINT, dsponinterrupt );
    if( dsp_oldinthandler==SIG_ERR || dsp_oldinthandler==SIG_DFL || 
	dsp_oldinthandler==SIG_IGN )
      dsp_oldinthandler= NULL;
  }
  return dsp.dspd;
}

void dspatexit(void)
{
  dspdesc(1);
}

void dsponinterrupt(int signr)
{
  dspdesc(1);
  if( dsp_oldinthandler )
    dsp_oldinthandler(signr);
  else
    exit(-1);
}

#endif



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
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "sndsys.h"

/*=============================================================================
    abc "file" (voice) *playopts*

Parses (part of) an .abc file.  The output channels are: Cue, volume and
frequency/ies.  Besides the .abc file name, this object takes the number of the
voice and a pointer to an abcplayopts structure as an argument.  This structure
contains the following variables: maxchord is the maximal number of notes in a
chord and determines the number of output channels of this part (maxchord+2).
If jazzynotes is !=0, pairs of notes with a single broken rhythm symbol (< or >)
in between get a length ratio of 2:1 instead of 3:1, which is common in Jazz.
accentsinchord determines whether accents should stand before the square
brackets enclosing a chord or within.  The abc standard gives an example in
which they are inside, but common practice seems to differ.  Since "accents"
include dynamics marks, I also find it more sensible to put them outside.  If
enforcechordlengths!=0, an error message is output if the lengths of different
notes in the same chord differ.  mfvolume is the volume to output for
mezzoforte, fortefactor the factor by which the volume of forte is larger (and
of piano is smaller).  Fortissimo and similar dynamical notations will be
larger or smaller by the appropriate number of fortefactors.  len, stacclen,
tenutolen are the lengths of ordinary, staccato and tenuto notes relative to
their nominal length.  The values cue, stacccue etc. specify the values which
are output on the cue channel (# 0) for the respective notes.  The cue values
are output for one sample at the start of each note.  If endcue is !=0, it is
output at the end of each tone (using the net ratheer than nominal note
length).

abc standard compliance: Pretty sloppy ;).  The only header fields interpreted
are T, L, Q, K and M.  In multi-voice files there MUST be EXACTLY one line
"V:#", where # is the number specifying the voice (starting at 1 or 0 as you
choose).  Parsing starts at the following line and ends at the next V: line or
at the end of the file.  Single-voice files are parsed from the first line
which does not have a ':' as the second character and is not empty or a
comment.

Some things that may deviate from the standard or be unexpected: The > or <
characters for broken rythms must immediately follow the preceding note.
Legato parentheses must follow afterwards.  All notes in a chord must be the
same length.  If not, the length of the first note counts.  If a chord is the
first note of a broken rythm, it is legal (with this parser) to add the > or <
characters after the bracket ] closing the chord.  The tuplets (5, (7 and (9
are always taken to be so many notes in the time of two.  In the abc standard,
that is only correct in the case of a 4/4 metre.  There is an additional accent
"!glissando!" which signifies a glissando on the note which it precedes (!) to
the pitch of the following one.  More than four piano marks ("!ppppp!" etc.)
set the volume to 0, which can be used for crescendi from nowhere.  More than
four forte marks amount to the same as four.  For (re)definition of accent
symbols, "u:" should be used, which defines how they are played, not "U:",
which determines how they are typeset.  The keys HP and Hp for highland bagpipe
music are not supported ;).  Last but not least, one should remember that in
abc, an empty line means the end of the tune and therefore stops the parser.
(However, this parser does not recognise lines that contain just spaces or tabs
as the end of the tune.)

See also: xyz
=============================================================================*/

struct abc {
  char      *buf;
  const char *pos, *title;
  abcplayopts *opts;
  int       voice;
  double    unitnote, wholelength, unitlength, barlength;
  double    currlength, accentfactor, tupletfactor;
  double    time, noteend, restend, lastbar;
  int       legato, tie, prevlegato, tupletcount, crescdim, iscresc, glissando;
  double    currvol, volinc, currfreq[PLENTY], freqfact[PLENTY];
  signed char accidentals[7], globaccs[7];
  const char *userdef_upper[19], *userdef_lower[16], *userdef_tilde;
};

float calcabc(sndobj *me, sndobj *, int);
int abcoutputnextnote( sndobj *me, struct abc *d );
int abcparsefraction( double *var, const char *start, int mode );
int abcparsetempo( struct abc *d, const char *tempo );
int abcparsekey( struct abc *d, const char *key );
int abcparseuserdef( struct abc *d, const char *read, const char **end );
int abcparsetuplet( struct abc *d, const char *read, const char **end );
float abcparseaccents( struct abc *d, const char *read, const char **end );
int cmpabcaccstr( const char *str1, const char *str2 );
void abcskipaccents( const char **pos );
double abcparsenote( struct abc *d, const char *read, const char **end, int freqonly );
void abcsetfreqfact( struct abc *d, const char *read );
void abcsetvolinc( struct abc *d, const char *read );

sndobj *abc( const char *filename, int voice, abcplayopts *opts )
{
  sndobj *p;
  struct abc *d;
  FILE *fp;
  const char *ext, *read, *endofheader;
  char *sharedname, *write, *title;
  long filesize;
  int count, retcode;
  
  for( ext= filename; *ext; ++ext );
  if( (ext -= 4) < filename || strcmp(ext, ".abc") ) {
    fprintf(stderr, "Warning: abc: File `%s' doesn't have extension abc. Returning 0.\n", filename);
    return c(0);
  }
  if( opts->maxchord < 1 ) {
    fprintf( stderr, "Warning: abc: Max. number of notes in chord must be at least 1.  Set to 1.\n" );
    opts->maxchord= 1;
  }
  else if( opts->maxchord > PLENTY-2 ) {
    fprintf( stderr, "Warning: abc: Max. number of notes in chord exceeded available output channels.  Set to Maximum (%d).\n", (int)(PLENTY-2) );
    opts->maxchord= PLENTY-2;
  }
  p= newsndo( calcabc, filename, "abc", 2+opts->maxchord, 0 );
  p->private[0]= d= new( struct abc );
  d->opts= opts;
  d->voice= voice;
  d->currlength= d->accentfactor= d->tupletfactor= 0.0;
  d->time= d->noteend= d->restend= d->lastbar= 0.0;
  d->currvol= opts->mfvolume;
  d->volinc= 0.0;
  for( count= 0; count< p->nch; ++count )
    d->freqfact[count]= 0.0;
  d->legato= d->tie= d->prevlegato= d->tupletcount= d->crescdim= d->glissando= 0;
  for( count= 0; count< 7; ++count )
    d->accidentals[count]= d->globaccs[count]= 0;
  for( count= 0; count< 16; ++count )
    d->userdef_upper[count]= d->userdef_lower[count]= NULL;
  d->userdef_upper[16]= d->userdef_upper[17]= d->userdef_upper[18]= 
    d->userdef_tilde= NULL;
  //  default accent symbols:
  d->userdef_upper['T'-'H']= "!trill!";
  d->userdef_upper['H'-'H']= "!fermata!";
  d->userdef_upper['L'-'H']= "!emphasis!";
  d->userdef_upper['M'-'H']= "!lowermordent!";
  d->userdef_upper['P'-'H']= "!uppermordent!";
  d->userdef_upper['S'-'H']= "!segno!";
  d->userdef_upper['O'-'H']= "!coda!";
  d->userdef_lower['u'-'h']= "!upbow!";
  d->userdef_lower['v'-'h']= "!downbow!";
  sharedname= news(char, strlen(filename)+1);
  strcpy(sharedname, filename );
  if( (d->buf= findshared(sharedname, 1))!=NULL )
  {
    free(sharedname);
    p->private[1]= d->buf;
    filesize= strlen(d->buf);
  }
  else {
    fp= fopen( filename, "r" );
    if( !fp ) {
      fprintf(stderr, "Warning: abc: Could not open file `%s'.  Returning 0.\n", filename );
      return c(0);
    }
    fseek( fp, 0L, SEEK_END );
    filesize= ftell( fp );
    fseek( fp, 0L, SEEK_SET );
    p->private[1]= d->buf= news( char, filesize + 3L );
    fread( d->buf, filesize, 1L, fp );
    if( ferror(fp) ) {
      fprintf(stderr, "Error: abc: Error reading data in `%s'. Aborting.\n", filename );
      exit(1);
    }
    fclose(fp);
    d->buf[filesize]= d->buf[filesize+1]= d->buf[filesize+2]= 0;
    registershared( d->buf, sharedname );
  }
  
  if( (read= strstr( d->buf, "\nK:" ))!=NULL ) {
    endofheader= read;
    retcode= abcparsekey( d, read+3 );
    if( retcode==0 )
      fprintf(stderr, "Warning: abc: Key could not be parsed in file `%s'.\n", filename);
    else if( retcode<0 )
      fprintf(stderr, "Warning: abc: Some explicit accidentals following the key could not be parsed in file `%s'.\n", filename);
  }
  else {
    fprintf(stderr, "Warning: abc: No K: header field found in file `%s'. Assuming C major. End of header unknown.\n", filename);
    endofheader= d->buf + filesize;
  }
  
  d->barlength= 0.0;
  if( (read= strstr( d->buf, "\nM:" ))!=NULL && read< endofheader )
    if( !abcparsefraction( &d->barlength, read+3, 0 ) ) {
      if( read[3]=='C' )
	d->barlength= 1.0;
      else
	fprintf(stderr, "Warning: abc: Could not parse metre in M: header field in file `%s'.\n", filename);
    }

  d->unitnote= 0.125;
  if( (read= strstr( d->buf, "\nL:" ))!=NULL && read< endofheader ) {
    if( !abcparsefraction( &d->unitnote, read+3, 0 ) )
      fprintf(stderr, "Warning: abc: Could not read unit note length in L: header field in file `%s'.\n", filename );
  }
  else if( d->barlength ) {
    if( d->barlength >= 0.75 )      d->unitnote= 0.125;
    else                            d->unitnote= 0.0625;
  }
  
  d->wholelength= 0.5 / d->unitnote;
  if( (read= strstr( d->buf, "\nQ:" ))!=NULL && read< endofheader )
    abcparsetempo( d, read+3 );
  d->unitlength= d->unitnote*d->wholelength;
  d->barlength *= d->wholelength;
  
  read= strstr( d->buf, "\nu:" ); 
  while( read && read< endofheader )
  {
    abcparseuserdef( d, read+3, &read );
    read= strstr( read, "\nu:" );
  }
  
  p->private[2]= title= news(char, NAMESTRLEN);
  if( (read= strstr( d->buf, "\nT:" ))!=NULL ) {
    read += 3;
    write= title;
    for( count= NAMESTRLEN-10; *read && *read!='\n' && *read!='%' && count> 0; --count )
      *write++= *read++;
    --write;
    while( isspace(*write) && write>=title )
      --write;
    ++write;
    *write++ = ' ';
  }
  else {
    count= strlen(filename);
    if( count<=NAMESTRLEN-10 ) {
      strcpy( title, filename );
      write= title+count;
    }
    else {
      title[0]= title[1]= title[2]= '.';
      for( read= filename+count-NAMESTRLEN+13, write= title+3, 
		count= NAMESTRLEN-13; count> 0; --count )
	*write++ = *read++;
    }
    *write++ = ' ';
  }

  d->pos= NULL;
  read= strstr( d->buf, "\nV:" );
  while( read ) {
    if( atoi( read+3 ) == d->voice && read >= endofheader) {
      read+=3;
      while( *read && *read!='\n' )
	++read;
      ++read;
      d->pos= read;
      snprintf(write, 10, "V:%d", d->voice);
      break;
    }
    read= strstr( read+3, "\nV:" );
  }
  if( !d->pos ) {
    if( *endofheader )      //  endofheader is well-defined
      d->pos= endofheader;
    else
      d->pos= d->buf;
    write[-1]= 0;
  }
  p->name= d->title= title;
  
  return p;
}



float calcabc(sndobj *me, sndobj *caller, int innr)
{
  struct abc *d;
  const char *oldpos;
  int count;
  
  OUTPUT(0, 0.0);   //  reset cue
  d= (struct abc *)me->private[0];
  if( d->noteend && d->time < d->noteend ) {
    if( d->crescdim )
      OUTPUT(1, d->currvol);
    if( d->glissando )
      for( count= 2; count< me->nch; ++count )
	me->ch[count] *= d->freqfact[count-2];
  }
  else {
    if( d->noteend ) {
      if( !d->legato && !d->tie )
	OUTPUT(0, d->opts->endcue);
      OUTPUT(1, 0.0);
      d->noteend= 0.0;
    }
    if( !d->restend || d->time >= d->restend ) {
      d->restend= 0.0;
      oldpos= d->pos;
      while( abcoutputnextnote( me, d ) && d->pos==oldpos ) {
	fprintf(stderr, "Warning: abc file `%s': Could not parse note at position %d. Trying next character.\n", d->title, (int)(d->pos-d->buf) );
      }
    }
  }
  if( d->crescdim )
    d->currvol += d->volinc;
  d->time += 1.0/SAMPRATE;
  CALCRETURN;
}


/*
Parses the next abc note at d->pos and sets d and the output values in me
accordingly.  If the end of the tune is reached, 0 is returned, otherwise 1.  */

#define BARTOLERANCE    1e-3

int abcoutputnextnote( sndobj *me, struct abc *d )
{
  const char *read, *endparse;
  double lengthfactor, newfreq;
  int count, retcode;
  
  read= d->pos;
  if( *read=='>' ) {
    lengthfactor= 0.5;
    while( *++read == '>' )
      lengthfactor *= 0.5;
    if( lengthfactor==0.5 && d->opts->jazzynotes )
      lengthfactor= 2.0/3.0;
  }
  else if( *read=='<' ) {
    lengthfactor= 0.5;
    while( *++read == '<' )
      lengthfactor *= 0.5;
    if( lengthfactor==0.5 && d->opts->jazzynotes )
      lengthfactor= 4.0/3.0;
    else
      lengthfactor= 2.0 - lengthfactor;
  }
  else
    lengthfactor= 1.0;
  while( !*read || isspace(*read) || *read==')' || *read==']' )
    if( !*read || (*read=='\n' && read[1]=='\n') ) {
      d->pos= read;
      OUTPUT(1, 0.0);
      d->noteend= 0.0;
      d->restend= DBL_MAX;
      return 0;
    }
    else ++read;

  while( 13 )       //  parse anything which may occur between notes
  {
    if( *read=='|' || read[1]=='|' || (*read==':' && read[1]==':') )
    {       //  bar symbols
      for( count= 0; count< 7; ++count )
	d->accidentals[count]= d->globaccs[count];
      if( d->barlength && d->lastbar ) {
	if( fabs(d->time - d->lastbar - d->barlength) > BARTOLERANCE )
	  fprintf(stderr, "Warning: abc file `%s': Length of bar before position %d deviates from metre by %g whole notes.\n", d->title, (int)(read-d->buf), (d->time - d->lastbar - d->barlength)/d->wholelength );
      }
      d->lastbar= d->time;
      if( read[1]=='|' )
	if( read[2]==']' || read[2]==':' || isdigit(read[2]) )  read+=3;
	else    read+=2;
      else
	if( read[1]==':' || read[1]==']' || isdigit(read[1]) )  read+=2;
	else    ++read;
    }
    else if( *read=='[' && isdigit(read[1]) ) {
      ++read;
      while( isdigit(*read) )
	++read;
    }
    else if( *read=='[' && (read[1]<'H' || read[1]>'Z') && 
			    (read[1]<'h' || read[1]>'w') )  //  start of chord
      break;
    //  "header" fields:
    else if( (read[-1]=='\n' && read[1]==':') || (*read=='[' && read[2]==':') )
    {
      if( *read=='[' )
	++read;
      switch( *read ) {
	case 'K':retcode= abcparsekey( d, read+2 );
		if( retcode==0 )
		  fprintf(stderr, "Warning: abc file `%s': Could not parse key"
		    " change at position %d.\n", d->title, (int)(read-d->buf) );
		else if( retcode<0 )
		  fprintf(stderr, "Warning: abc file `%s': Could not parse "
		    "explicit accidentals in key change at position %d.\n", 
			d->title, (int)(read-d->buf) );
	    break;
	case 'L':if( !abcparsefraction( &d->unitnote, read+2, 0 ) )
		    fprintf(stderr, "Warning: abc file `%s': Could not parse "
			"change of default note length at position %d.\n", 
			d->title, (int)(read-d->buf) );
		else
		    d->unitlength= d->unitnote*d->wholelength;
	    break;
	case 'M':d->barlength /= d->wholelength;
		if( read[2] == 'C' )
		  d->barlength= 1.0;
		else
		  if( !abcparsefraction( &d->barlength, read+2, 0 ) )
		    fprintf(stderr, "Warning: abc file `%s': Could not parse "
		  "metre change at position %d.", d->title, (int)(read-d->buf));
		d->barlength *= d->wholelength;
	    break;
	case 'Q':if( !abcparsetempo( d, read+2 ) )
		    fprintf(stderr, "Warning: abc file `%s': Could not parse "
	      "tempo change at position %d.\n", d->title, (int)(read-d->buf) );
	    break;
	case 'u':abcparseuserdef( d, read+2, &read );
	  break;
	default:
	    break;
      }
      if( read[-1]=='[' )
	while( *read && *read!=']' )
	  ++read;
      else 
	while( *read && *read!='\n' )
	  ++read;
      ++read;
    }
    else if( *read=='{' ) {     //  grace notes are not played
      while( *read && *read !='}' )
	++read;
      ++read;
    }
    else if( *read=='%' ) {     //  comment to end of line
      while( *read && *read !='\n' )
	++read;
    }
    else if( *read=='\\' ) {    //  line-linking backslash
      ++read;
      while( *read && *read !='\n' ) {
	if( *read=='%' )
	  while( *read && *read!='\n' )
	    ++read;
	else  {
	  if( !isspace(*read) )
	    fprintf(stderr, "Warning: abc file `%s': non-whitespace characters"
		" after line-linking backslash at position %d.\n", d->title, 
		(int)(read-d->buf) );
	  ++read;
	}
      }
    }
    else
      break;
    while( !*read || isspace(*read) )
      if( !*read || (*read=='\n' && read[1]=='\n') ) {
	d->pos= read;
	OUTPUT(1, 0.0);
	d->noteend= 0.0;
	d->restend= DBL_MAX;
	return 0;
      }
      else ++read;
  }
  
  if( *read=='"' ) {        //  Guitar chords are ignored for now
    ++read;
    while( *read && *read!='"' )
      ++read;
  }
  
  d->prevlegato= d->legato || d->tie;
  
  while( *read=='(' )
    if( isdigit(read[1]) ) {
      if( d->tupletcount )
	fprintf(stderr, "Warning: abc file `%s': Nested tuplets found at "
		"position %d. Not supported.\n", d->title, (int)(read-d->buf) );
      ++read;
      if( !abcparsetuplet( d, read, &endparse ) )
	fprintf(stderr, "Warning: abc file `%s': Error parsing "
		"tuplet at position %d. Tuplet ignored.\n", 
		d->title, (int)(read-d->buf) );
      read= endparse;
    }
    else {
      ++d->legato;
      ++read;
    }
  
  d->tie= 0;
  d->glissando= 0;
  d->currlength= 0.0;
  if( !d->opts->accentsinchord )
    OUTPUT(0, abcparseaccents( d, read, &read ));
  if( *read=='[' )      //  chord
  {
    ++read;
    if( d->opts->accentsinchord )
      OUTPUT(0, abcparseaccents( d, read, &read ));
    for( count= 0; count< d->opts->maxchord && *read!=']'; ++count ) {
      if( *read=='.' || *read=='!' || *read=='~' || 
	    (*read>='h' && *read<='w') || (*read<='H' && *read>='Z') ) {
	if( d->opts->accentsinchord )
	  fprintf(stderr, "Warning: abc file `%s': Different accents for different chord notes are ignored (position %d).\n", d->title, (int)(read-d->buf) );
	else
	  fprintf(stderr, "Warning: abc file `%s': No accents within chord allowed with this options setting (position %d).\n", d->title, (int)(read-d->buf) );
	abcskipaccents( &read );
      }
      newfreq= abcparsenote( d, read, &endparse, 0 );
      if( newfreq > 1.0 || !d->opts->silentfreq )
        d->currfreq[count]= newfreq;
      OUTPUT(count+2, d->currfreq[count]);
      read= endparse;
    }
    for( ; count< d->opts->maxchord; ++count ) {
      d->currfreq[count]= 0.0;
      OUTPUT(2+count, 0.0);
    }
    while( *read && *read!=']' )
      ++read;
    ++read;
    if( *read=='<' ) {
      lengthfactor= 0.5;
      while( *++read == '<' )
	lengthfactor *= 0.5;
      if( lengthfactor==0.5 && d->opts->jazzynotes )
	lengthfactor= 2.0/3.0;
    }
    else if( *read=='>' ) {
      lengthfactor= 0.5;
      while( *++read == '>' )
	lengthfactor *= 0.5;
      if( lengthfactor==0.5 && d->opts->jazzynotes )
	lengthfactor= 4.0/3.0;
      else
	lengthfactor= 2.0 - lengthfactor;
    }
  }
  else {
    if( d->opts->accentsinchord )
      OUTPUT(0, abcparseaccents( d, read, &read ));
    newfreq= abcparsenote( d, read, &endparse, 0 );
    if( newfreq > 1.0 || !d->opts->silentfreq )
      d->currfreq[0]= newfreq;
    OUTPUT(2, d->currfreq[0]);
    read= endparse;
    for( count= 1; count< d->opts->maxchord; ++count ) {
      d->currfreq[count]= 0.0;
      OUTPUT(2+count, 0.0);
    }
  }

  while( *read==')' ) {
    if( d->legato>0 )
      --d->legato;
    ++read;
  }
  
  d->currlength *= lengthfactor;
  if( d->tupletcount ) {
    d->currlength *= d->tupletfactor;
    --d->tupletcount;
  }
  if( d->legato || d->tie )
    d->accentfactor= 1.0;
  if( d->accentfactor>0.0 ) {
    d->noteend= d->time + d->currlength*d->accentfactor;
    OUTPUT(1, d->currvol);
  }
  else {
    d->noteend= 0.0;
    OUTPUT(1, 0.0);
  }
  if( d->accentfactor<1.0 )
    d->restend= d->time + d->currlength;
  else
    d->restend= 0.0;

  if( d->glissando )
    abcsetfreqfact( d, read );
  if( d->crescdim && !d->volinc )
    abcsetvolinc( d, read );

  d->pos= read;
  while( read[-1]==')' )
    --read;
  if( read[-1]==']' )
    --read;
  if( read[-1]=='<' )  {
    while( read[-1]=='<' )
      --read;
    d->pos= read;
  }
  else if( read[-1]=='>' )  {
    while( read[-1]=='>' )
      --read;
    d->pos= read;
  }
  
//fprintf(stderr, "abcoutputnextnote: cue %g, vol %g, freq %g; length %g\n", me->ch[0], me->ch[1], me->ch[2], d->currlength);
//fprintf(stderr, "abcoutputnextnote: cue %g, vol %g, freq %g, %g, %g; length %g\n", me->ch[0], me->ch[1], me->ch[2], me->ch[3], me->ch[4], d->currlength);
//fprintf(stderr, "position: %.10s\n", d->pos );
  return !!*read;
}



/*
Parses a fraction starting at <start> which is of the type:
[white space] decimal digits [white space] '/' [white space] decimal digits
and writes the result to <var>.  Returns 1 if successful; otherwise no 
assignment is made and 0 is returned.  If <ignorenumerator> is !=0, the result 
is 1/denominator.
*/
int abcparsefraction( double *var, const char *start, int ignorenumerator )
{
  const char *endparse;
  long num, den;
  
  while( *start!='\n' && isspace(*start) )
    ++start;
  if( ignorenumerator ) {
    while( isdigit(*start) )
      ++start;
    num= 1L;
  }
  else {
    num= strtol( start, (char **)&endparse, 10 );
    if( endparse==start )
      return 0;
    start= endparse;
  }
  while( *start!='\n' && isspace(*start) )
    ++start;
  if( *start!='/' )
    return 0;
  ++start;
  while( *start!='\n' && isspace(*start) )
    ++start;
  den= strtol( start, (char **)&endparse, 10 );
  if( endparse==start )
    return 0;
  *var= (double)num/(double)den;
  return 1;
}


/*
Parses an abc tempo field of the form <tempo of default note length> or <note
length num.> "/" <note length denom.> "=" <tempo>.  If it is of the first form,
d->unitnote has to be defined already.  If successful, d->wholelength and
d->unitlength are assigned proper values.  Returns 1 on success, 0 otherwise.  */

int abcparsetempo( struct abc *d, const char *tempo )
{
  char *endparse;
  long notenum, noteden, notetempo;
  
  notenum= strtol( tempo, &endparse, 10 );
  if( endparse==tempo )
    return 0;
  tempo= endparse;
  if( isspace(*tempo) || *tempo=='%' ) {
    notetempo= notenum;
    d->wholelength= 60.0 / (notetempo * d->unitnote);
    d->unitlength= 60.0 / (double)notetempo;
    return 1;
  }
  if( *tempo!='/' )
    return 0;
  ++tempo;
  noteden= strtol( tempo, &endparse, 10 );
  if( endparse==tempo )
    return 0;
  tempo= endparse;
  if( *tempo!='=' )
    return 0;
  ++tempo;
  notetempo= strtol( tempo, &endparse, 10 );
  if( endparse==tempo )
    return 0;
  d->wholelength= 60.0 / ((double)notetempo*(double)notenum/(double)noteden);
  d->unitlength= d->unitnote*d->wholelength;
  return 1;
}


/*
Parses key in abc notation.  must point to an abc structure.  Its globaccs and
accidentals fields are filled with values between -2 (double flat) and 2
(double sharp).  The return value is 1 if everything went fine, 0 if the key
could not be parsed, or -1 if not all of the explicit additional accidentals
after the key could be parsed.  In the former case no assignment is made.  */

int abcparsekey( struct abc *d, const char *key )
{
  // number of sharps (>0)/flats (<0) for major scale starting [index]
  // semitones above C
  static int baseaccs[12]= { 0, -5, 2, -3, 4, -1, 6, 1, -4, 3, -2, 5 };
  // order in which accidentals are set.  0 is C, 1 is D, ...
  static int sharpsorder[6]= { 3, 0, 4, 1, 5, 2 };
  static int flatsorder[6]= { 6, 2, 5, 1, 4, 0 };
  
  int base, fsharp, count, naccs, acc, retcode;
  
  while( *key!='\n' && isspace(*key) )
    ++key;
  if( *key<'A' || *key>'G' )
    return 0;
  base= *key-'C';
  if( base==6 )
    base= -1;       // H is B
  base *= 2;
  if( base<0 )  ++base;
  else if( base>4 ) --base;
  ++key;
  fsharp= 0;
  if( *key=='#' ) {
    ++base;
    if( base==6 )
      fsharp= 1;
  }
  else if( *key=='b' )  --base;
  while( *key!='\n' && isspace(*key) )
    ++key;
  switch( *key )    // adjust base to get equivalent major key
  {
    case 'a': case 'A':
	    if( (key[1]=='e' || key[1]=='E') && (key[2]=='o' || key[2]=='O') )
		    base += 3;
	    else    return 0;
	break;
    case 'd': case 'D':
	    if( (key[1]=='o' || key[1]=='O') && (key[2]=='r' || key[2]=='R') )
		    base -= 2;
	    else    return 0;
	break;
    case 'i': case 'I':
	    if( (key[1]!='o' && key[1]!='O') || (key[2]!='n' && key[2]!='N') )
		return 0;
	break;
    case 'l': case 'L':
	    if( (key[1]=='y' || key[1]=='Y') && (key[2]=='d' || key[2]=='D') )
		    base -= 5;
	    else if( (key[1]=='o' || key[1]=='O') && 
		    (key[2]=='k' || key[2]=='K') )          base += 1;
	    else    return 0;
	break;
    case 'm': case 'M':
	    if( !isalpha(key[1]) || ((key[1]=='i' || key[1]=='I') 
		&& (key[2]=='n' || key[2]=='N')) )          base += 3;
	    else if( (key[1]=='i' || key[1]=='I') && 
			(key[2]=='x' || key[2]=='X') )      base -= 7;
	    else if( (key[1]!='a' && key[1]!='A') || 
			(key[2]!='j' && key[2]!='J') )      return 0;
	break;
    case 'p': case 'P':
	    if( (key[1]=='h' || key[1]=='H') && (key[2]=='r' || key[2]=='R') )
		    base -= 4;
	    else    return 0;
	break;
    case '_': case '=': case '^': case '\n': case '%': case ']':
	break;
    default:return 0;
  }
  while( base<0 )
    base += 12;
  while( base>11 )
    base -= 12;
  if( base==6 && !fsharp )
    naccs= -6;
  else
    naccs= baseaccs[base];
  if( naccs>=0 ) {
    for( count= 0; count< naccs; ++count )
      d->globaccs[sharpsorder[count]]= 1;
  }
  else {
    naccs= -naccs;
    for( count= 0; count< naccs; ++count )
      d->globaccs[flatsorder[count]]= -1;
  }
  
  // now for explicit additional accidentals
  retcode= 1;
  while( isalpha(*key) )
    ++key;
  while( *key!='\n' && isspace(*key) )
    ++key;
  while( *key=='_' || *key=='^' || *key=='=' )
  {
    if( *key=='=' )     acc= 0;
    else if( *key=='_' ) {
      if( key[1]=='_' ) { acc= -2;  ++key; }
      else              acc= -1;
    }
    else {
      if( key[1]=='^' ) { acc= 2;  ++key; }
      else              acc= 1;
    }
    ++key;
    base= tolower(*key) - 'c';
    if( base<-2 || base>4 )
      retcode= -1;
    else {
      if( base<0 )
	base += 7;
      d->globaccs[base]= acc;
    }
    ++key;
    while( *key!='\n' && isspace(*key) )
      ++key;
  }
  
  for( count= 0; count< 7; ++count )
    d->accidentals[count]= d->globaccs[count];
  
  return retcode;
}


/*
Parses the contents of the u: header field which redefines the characters
corresponding to certain accents.  <read> must point to the first character
after "u:", <end> is a pointer to which the position of the first character of
the following line is written, and <d> points to an abc structure to which the
definition is written (components userdef_...).  The return value is 1 on
success, 0 otherwise.  Error messages are already output by this function.  */

int abcparseuserdef( struct abc *d, const char *read, const char **end )
{
  const char *accstr;
  char accchar;

  while( isspace(*read) && *read!='\n' )
    ++read;
  accchar= *read++;
  if( accchar!='~' && (accchar<'H' || accchar>'Z') && 
			(accchar<'h' || accchar>'w') )      {
    fprintf(stderr, "Warning: abc file `%s': Character in accent definition at position %d cannot be redefined. Definition ignored.\n", d->title, (int)(read-d->buf) );
    while( *read && *read!='\n' )
      ++read;
    if( end )
      *end= read;
    return 0;
  }
  while( isspace(*read) && *read!='\n' )
    ++read;
  if( *read++!='=' ) {
    fprintf(stderr, "Warning: abc file `%s': Could not parse accent definition: Equals sign expected at position %d.\n", d->title, (int)(read-d->buf) );
    while( *read && *read!='\n' )
      ++read;
    if( end )
      *end= read;
    return 0;
  }
  while( isspace(*read) && *read!='\n' )
    ++read;
  if( *read!='!' ) {
    fprintf(stderr, "Warning: abc file `%s': Accent string in accent definition at position %d does not start with an `!'. Definition ignored.\n", d->title, (int)(read-d->buf) );
    while( *read && *read!='\n' )
      ++read;
    if( end )
      *end= read;
    return 0;
  }
  if( !cmpabcaccstr(read, "!none!") || !cmpabcaccstr(read, "!nil!") )
    accstr= NULL;
  else
    accstr= read;
  if( accchar=='~' )
    d->userdef_tilde= accstr;
  else if( isupper(accchar) )
    d->userdef_upper[accchar-'H']= accstr;
  else
    d->userdef_lower[accchar-'h']= accstr;
  while( *read && *read!='\n' )
    ++read;
  if( end )
    *end= read;
  return 1;
}


/*
Parses an abc expression for n-tuplets.  <d> is a pointer to an abc structure,
the components tupletfactor and tupletcount of which are set.  <read> must
point to the first digit of the tuplet expression (not the parenthesis!).
Returns 1 if successful, 0 on an unspecified error, and -1 if an error message
has already been output.  In the last two cases no assignment is made.  */

int abcparsetuplet( struct abc *d, const char *read, const char **end )
{
  //  default tuplet note length numerator depending on denominator
  //  (5, 7, 9 conforming to standard only for 4/4 metre)
  static long defnum[10]= { 0, 0, 3, 2, 3, 2, 2, 2, 3, 2 };
  
  char *endparse;
  long tupletn, tupletnum, tupletden;
  
  tupletden= strtol( read, &endparse, 10 );
  if( tupletden<2 ) {
    fprintf(stderr, "Warning: abc file `%s': Denominator of tuplet smaller "
		"than 2 or unreadable at position %d. Tuplet ignored.\n", 
		d->title, (int)(read-d->buf) );
    if( end )
      *end= endparse;
    return -1;
  }
  read= endparse;
  if( *read!=':' )
  {
    if( tupletden<10 ) {
      d->tupletcount= tupletden;
      d->tupletfactor= (double)defnum[tupletden]/(double)tupletden;
      if( end )
	*end= read;
      return 1;
    }
    else {
      fprintf(stderr, "Warning: abc file `%s': No default length for %ld>9 "
		"tuplet at position %d. Tuplet ignored.\n", 
		d->title, tupletden, (int)(read-d->buf) );
      if( end )
	*end= read;
      return -1;
    }
  }
  ++read;
  if( *read==':' ) {    //  numerator not given explicitly
    if( tupletden<10 )
      tupletnum= defnum[tupletden];
    else {
      fprintf(stderr, "Warning: abc file `%s': No default length for %ld>9 "
		"tuplet at position %d. Tuplet ignored.\n", 
		d->title, tupletden, (int)(read-d->buf) );
      while( *read==':' || isdigit(*read)  )
	++read;
      if( end )
	*end= read;
      return -1;
    }
  }
  else if( isdigit(*read) ) {
    tupletnum= strtol( read, (char **)&read, 10 );
  }
  else {
    while( *read==':' || isdigit(*read)  )
      ++read;
    if( end )
      *end= read;
    return 0;
  }
  if( *read!=':' ) {
    if( end )
      *end= read;
    return 0;
  }
  ++read;
  if( !isdigit(*read) )
    tupletn= tupletden;
  else {
    tupletn= strtol( read, &endparse, 10 );
    read= endparse;
  }
  d->tupletcount= tupletn;
  d->tupletfactor= (double)tupletnum/(double)tupletden;
  if( end )
    *end= read;
  return 1;
}


/*
Parses the accents preceding an abc note.  <d> points to an abc structure which
may be modified if relevant accents are found (components currvol, crescdim and
glissando); <read> points to the start of the accent strings; and <end> is a
pointer to which the position of the first character which does not belong to
an accent is written.  Return value is the cue value determined from d->opts
and the relevant accents.  

As an extension of the abc standard, "!glissando!" is recognised as signifying
a glissando on the current note to the pitch of the next.  Besides, more than
four piano marks ("!ppppp!" etc.) set the volume to 0.  More than four forte
marks amount to the same as four.  */

float abcparseaccents( struct abc *d, const char *read, const char **end )
{
  const char *accstr, *scan;
  float cue, cueadd;
  int count;
  
  cue= d->prevlegato? d->opts->legatocue : d->opts->cue;
  cueadd= 0.0;
  d->accentfactor= d->opts->len;
  while( *read=='!' || *read=='.' || *read=='~' ||
	(*read>='h' && *read<='w') || (*read>='H' && *read<='Z') )
  {
    if( *read=='.' ) {
      cue= d->prevlegato? d->opts->legstacccue : d->opts->stacccue;
      if( !d->prevlegato )
	d->accentfactor= d->opts->stacclen;
      ++read;
      continue;
    }
    if( *read=='!' )
      accstr= read;
    else if( *read=='~' )
      accstr= d->userdef_tilde;
    else if( isupper(*read) )
      accstr= d->userdef_upper[*read-'H'];
    else
      accstr= d->userdef_lower[*read-'h'];
    
    if( accstr )
    {
      for( scan= accstr+2; *scan==accstr[1]; ++scan );
      if( *scan=='!' && (accstr[1]=='p' || accstr[1]=='f') )
      {
	d->currvol= d->opts->mfvolume;
	count= scan-accstr-1;
	if( accstr[1]=='f' ) {
	  if( count>=5 )
	    count= 4;
	  for( ; count> 0; --count )
	    d->currvol *= d->opts->fortefactor;
	}
	else {
	  if( count>=5 )
	    d->currvol= 0;
	  else
	    for( count= scan-accstr-1; count> 0; --count )
	      d->currvol /= d->opts->fortefactor;
	}
      }
      else if( !cmpabcaccstr(accstr, "!mf!") ) {
	d->currvol= d->opts->mfvolume;
      }
      else if( !cmpabcaccstr(accstr, "!crescendo(!") ) {
	d->crescdim= 1;
	d->volinc= 0.0;
	d->iscresc= 1;
      }
      else if( !cmpabcaccstr(accstr, "!diminuendo(!") ) {
	d->crescdim= 1;
	d->volinc= 0.0;
	d->iscresc= 0;
      }
      else if( !cmpabcaccstr(accstr, "!crescendo)!") ||
	      !cmpabcaccstr(accstr, "!diminuendo)!") ) {
	d->crescdim= 0;
      }
      else if( !cmpabcaccstr(accstr, "!glissando!") ) {
	d->glissando= 1;
      }
      else if( !cmpabcaccstr(accstr, "!tenuto!") ) {
	d->accentfactor= d->opts->tenutolen;
      }
      
      for( count= 0; count < PLENTY; ++count )
	if( d->opts->accentstring[count] && 
	      !cmpabcaccstr(accstr, d->opts->accentstring[count]) )  {
	  if( d->opts->accentcue[count]>= 0 )
	    cue= d->opts->accentcue[count];
	  else
	    cueadd += -d->opts->accentcue[count];
	}
    }
    else
      fprintf(stderr, "Warning: abc file `%s': Undefined accent character `%c' encountered at position %d. Ignored.\n", d->title, *read, (int)(read-d->buf) );
    if( *read=='!' ) {
	for( ++read; *read && *read!='!'; ++read );
	++read;
    }
    else
	++read;
  }
  if( end )
    *end= read;

  return cue + cueadd;
}


/*
Compares two accent strings enclosed in exclamation marks.  The return value is
as for strcmp.  */

int cmpabcaccstr( const char *str1, const char *str2 )
{
  do {
    if( *str1!=*str2 )
	return *str2 - *str1;
    ++str1;
    ++str2;
  }
  while( *str1 && *str1!='!' );
  
  if( *str1!=*str2 )
    return *str2 - *str1;
  else
    return 0;
}


/*
Skips the accents after *pos.  */

void abcskipaccents( const char **pos )
{
  while( 13 )
  {
    if( **pos=='!' )
      for( ++*pos; **pos!='!'; ++*pos );
    else if( **pos=='.' || **pos=='~' || 
	    (**pos>='h' && **pos<='w') || (**pos<='H' && **pos>='Z') )
      ++*pos;
    else
      return;
  }
}


/*
Parses the accidental, note and note length in an abc file.  <read> should
point to the position of the note in the abc file, <end> will be assigned the
position of the first character not belonging to this note.  <d> must point to
an abc structure.  Unless <freqonly> is !=0, the component d->currlength is
assigned the right value, d->tie may be set, the current accidentals may be
adjusted and d->accentfactor is set to 0 in the case of rests.  The frequency
is returned.  All notes in a chord have to have the same length (probably
deviating from the abc standard).  If this note is the second or later one in a
chord, d->currlength should contain the length of the previous note so this can
be checked.  If the note lengths in a chord differ, and if
d->opts->enforcechordlength is set, an error message is output.  */

double abcparsenote( struct abc *d, const char *read, const char **end, int freqonly )
{
  static int scalesemitones[]= { 0, 2, -9, -7, -5, -4, -2 };
  
  const char *firstdotchar;
  double freq, dotlengthfactor, notelength;
  long lengthden, lengthnum;
  int accidental, semitones, noteindex;
  
  accidental= 3;
  if( *read=='=' || *read=='_' || *read=='^' )
  {
    if( *read=='=' )
      accidental= 0;
    else if( *read=='_' )
    {
      accidental= -1;
      if( read[1]=='_' ) {
	accidental= -2;
	++read;
      }
    }
    else if( *read=='^' )
    {
      accidental= 1;
      if( read[1]=='^' ) {
	accidental= 2;
	++read;
      }
    }
    ++read;
  }
  
  if( *read>='A' && *read<='G' ) {
    freq= 110.0;
    semitones= scalesemitones[*read-'A'];
    noteindex= *read-'C';
    if( noteindex< 0 )
      noteindex += 7;
    ++read;
    if( *read==',' )
      for( ; *read==','; ++read )
        freq /= 2.0;
  }
  else if( *read>='a' && *read<='g' ) {
    freq= 220.0;
    semitones= scalesemitones[*read-'a'];
    noteindex= *read-'c';
    if( noteindex< 0 )
      noteindex += 7;
    ++read;
    if( *read=='\'' )
      for( ; *read=='\''; ++read )
        freq *= 2.0;
  }
  else {
    if( *read!='z' )
      fprintf( stderr, "Warning: abc file `%s': Illegal note symbol `%c' at position %d. Interpreted as rest.\n", d->title, *read, (int)(read-d->buf) ); 
    freq= 1.0;
    noteindex= -1;
    if( !freqonly )
      d->accentfactor= 0.0;
    ++read;
  }

  if( noteindex>=0 ) {
    if( accidental!=3 ) {
      if( !freqonly )
	d->accidentals[noteindex]= accidental;
    }
    else
      accidental= d->accidentals[noteindex];
    semitones += accidental;
    for( ; semitones>0; --semitones )
      freq *= SEMITONE;
    for( ; semitones<0; ++semitones )
      freq /= SEMITONE;
    freq *= d->opts->transposefactor;
  }

  if( isdigit(*read) )
    lengthnum= strtol( read, (char**)&read, 10 );
  else
    lengthnum= 1;
  
  if( *read=='/' ) {
    if( isdigit(read[1]) )
      lengthden= strtol( read+1, (char**)&read, 10 );
    else
      for( ++read, lengthden= 2; *read=='/'; ++read )
	lengthden *= 2;
  }
  else
    lengthden= 1;
  notelength= d->unitlength*(double)lengthnum/(double)lengthden;
  
  if( *read=='-' ) {
    if( !freqonly )
      d->tie= 1;
    ++read;
  }
  else if( *read=='>' || *read=='<' ) {
    firstdotchar= read;
    dotlengthfactor= 0.5;
    for( ++read; *read==*firstdotchar; ++read )
      dotlengthfactor *= 0.5;
    if( dotlengthfactor==0.5 && d->opts->jazzynotes )
      dotlengthfactor= 2.0/3.0;
    if( *firstdotchar=='>' )
      dotlengthfactor= 2.0 - dotlengthfactor;
    notelength *= dotlengthfactor;
  }
  
  if( end )
    *end= read;
  
  if( !freqonly ) {
    if( d->currlength ) {
      if( d->opts->enforcechordlength && d->currlength!=notelength )
	fprintf(stderr, "Warning: abc file `%s': Differing lengths in chord before position %d. First note wins.\n", d->title, (int)(read-d->buf) );
    }
    else
      d->currlength= notelength;
  }
  
  return freq;
}


/*
Parses the following note or chord (but not rest) and stores the frequency
multipliers for glissando in the component d->freqfact of the abc structure
<d>.  <read> has to be the position in the abc file after the current note.  If
there was no next note, d->glissando is set to 0.  */

void abcsetfreqfact( struct abc *d, const char *read )
{
  signed char saveaccs[7];
  int count, chord= 0;
  
  for( count= 0; count< 7; ++count )
    saveaccs[count]= d->accidentals[count];
  
  while( (*read<'a' || *read>'g') && (*read<'A' || *read>'G')
	    && *read!='=' && *read!='_' && *read!='^' )
  {
    if( *read=='!' ) {
      ++read;
      while( *read && *read!='!' )
	++read;
    }
    else if( *read=='"' ) {
      ++read;
      while( *read && *read!='"' )
	++read;
    }
    else if( *read=='[' ) {
      if( read[1]=='|' ) {
	for( count= 0; count< 7; ++count )
	  d->accidentals[count]= d->globaccs[count];
	++read;
      }
      else if( read[2]==':' )
	while( *read && *read!=']' )
	  ++read;
      else if( !isdigit(read[1]) )
	chord= 1;
    }
    else if( !*read || (*read=='\n' && read[1]=='\n') ) {
      d->glissando= 0;
      return;
    }
    if( *read=='|' || read[1]=='|' || (*read==':' && read[1]==':') )
      for( count= 0; count< 7; ++count )
	d->accidentals[count]= d->globaccs[count];
    ++read;
  }
  
  d->freqfact[0]= abcparsenote( d, read, &read, 1 );
  if( chord ) {
    for( count= 1; count< d->opts->maxchord; ++count )  {
      while( (*read<'a' || *read>'g') && (*read<'A' || *read>'G')
	    && *read!='=' && *read!='_' && *read!='^' && *read!=']' &&
	    *read && (*read!='\n' || read[1]!='\n') )
	++read;
      if( !*read || (*read=='\n' && read[1]=='\n') )
	break;
      if( *read!=']' )
	d->freqfact[count]= abcparsenote( d, read, &read, 1 );
      else
	break;
    }
    for( ; count< d->opts->maxchord; ++count )
      d->freqfact[count]= 0.0;
  }
  else
    for( count= 1; count< d->opts->maxchord; ++count )
      d->freqfact[count]= 0.0;

  for( count= 0; count< 7; ++count )
    d->accidentals[count]= saveaccs[count];

  for( count= 0; count< d->opts->maxchord && d->currfreq[count]
		  && d->freqfact[count]; ++count ) {
    d->freqfact[count]= exp( 1.0/(d->currlength*SAMPRATE)
			  *log(d->freqfact[count]/d->currfreq[count]) );
  }
  for( ; count< d->opts->maxchord; ++count )
    d->freqfact[count]= d->freqfact[count-1];
}


/*
Skips through the part of the song following <read> until the end of the
current crescendo or diminuendo (d->iscresc decides which).  d->volinc is set
appropriately.  If any dynamics mark, start of a new (de)crescendo, or a
mismatching end (of dim. instead of cresc. or vice versa) is found, an error
message is output and the (de)crescendo is concluded.  The first dynamics mark
between the end-of-(de)crescendo mark and the following note is used as the
volume to be reached by the end of the crescendo, unless it is softer than the
original volume for crescendo or louder for diminuendo or equal to it.  If no
dynamics mark is found, the difference in dynamics is assumed to be two factors
of d->opts->fortefactor.  The fact that the first dynamics mark is taken makes
it possible to switch to a different volume immediately afterwards, for
instance: "!crescendo)!!ff!!p! ..."  */

void abcsetvolinc( struct abc *d, const char *read )
{
  const char *scan;
  double timespan, newvol;
  long lengthnum, lengthden;
  int count, chord;
  
  chord= 0;
  newvol= 0.0;
  timespan= d->currlength;
  while( 13 )
  {
    if( !*read || (*read=='\n' && read[1]=='\n') )
      break;
    if( *read=='!' ) {
      if( !cmpabcaccstr( read, "!crescendo)!" ) ) {
	if( !d->iscresc )
	  fprintf( stderr, "Warning: abc file `%s': Encountered end of crescendo while searching for end of diminuendo at position %d.\n", d->title, (int)(read-d->buf) );
	break;
      }
      if( !cmpabcaccstr( read, "!diminuendo)!" ) ) {
	if( d->iscresc )
	  fprintf( stderr, "Warning: abc file `%s': Encountered end of diminuendo while searching for end of crescendo at position %d.\n", d->title, (int)(read-d->buf) );
	break;
      }
      if( !cmpabcaccstr( read, "!crescendo(!" ) || 
	    !cmpabcaccstr( read, "!diminuendo(" ) ) {
	fprintf( stderr, "Warning: abc file `%s': Encountered start of new (de)crescendo while searching for end of previous one at position %d. Ending first (de)crescendo here.\n", d->title, (int)(read-d->buf) );
	break;
      }
      for( scan= read+2; *scan==read[1]; ++scan );
      if( (*scan=='!' && (read[1]=='f' || read[1]=='p')) ||
	    !cmpabcaccstr(read, "!mf!") ) {
	fprintf( stderr, "Warning: abc file `%s': Encountered dynamics mark at position %d while searching for end of (de)crescendo. Ending (de)crescendo here.\n", d->title, (int)(read-d->buf) );
	break;
      }
      ++read;
      while( *read && *read!='!' )
	++read;
    }
    else if( *read=='"' ) {
      ++read;
      while( *read && *read!='"' )
	++read;
    }
    else if( *read=='[' ) {
      if( read[2]==':' ) {
	if( read[1]!='|' )
	  while( *read && *read!='\n' && *read!=']' )
	    ++read;
      }
      else if( !isdigit(read[1]) )
	chord= 1;
    }
    else if( *read=='%' )
      while( *read && *read!='\n' )
	++read;
    else if( (*read>='a' && *read<='g') || (*read>='A' && *read<='G') || *read=='z' )
    {
      ++read;
      if( isdigit(*read) )
	lengthnum= strtol( read, (char**)&read, 10 );
      else
	lengthnum= 1;
      
      if( *read=='/' ) {
	if( isdigit(read[1]) )
	  lengthden= strtol( read+1, (char**)&read, 10 );
	else
	  for( ++read, lengthden= 2; *read=='/'; ++read )
	    lengthden *= 2;
      }
      else
	lengthden= 1;
      timespan += d->unitlength*(double)lengthnum/(double)lengthden;
      if( chord ) {
	while( *read && *read!=']' )
	  ++read;
	chord= 0;
      }
      else
	--read;
    }
    ++read;
  }
  
  while( (*read<'a' || *read>'g') && (*read<'A' || *read>'G') && *read!='z'
	&& *read && (*read!='\n' || read[1]!='\n') )
  {
    if( *read=='!' )
    {
      if( !cmpabcaccstr( read, "!mf!" ) ) {
	newvol= d->opts->mfvolume;
	break;
      }
      for( scan= read+2; *scan==read[1]; ++scan );
      if( *scan=='!' && (read[1]=='p' || read[1]=='f') )
      {
	newvol= d->opts->mfvolume;
	count= scan-read-1;
	if( read[1]=='f' ) {
	  if( count>=5 )
	    count= 4;
	  for( ; count> 0; --count )
	    newvol *= d->opts->fortefactor;
	}
	else {
	  if( count>=5 )
	    newvol= 0;
	  else
	    for( ; count> 0; --count )
	      newvol /= d->opts->fortefactor;
	}
	break;
      }
      ++read;
      while( *read && *read!='!' )
	++read;
    }
    else if( *read=='%' )
      while( *read && *read!='\n' )
	++read;
    ++read;
  }
  
  if( (d->iscresc && newvol<=d->currvol) || (!d->iscresc && newvol>=d->currvol) )
    newvol= 0;
  if( newvol )
    d->volinc= (newvol-d->currvol)/(timespan*SAMPRATE);
  else if( d->iscresc )
    d->volinc= (d->opts->fortefactor*d->opts->fortefactor-1.0)*d->currvol /
		    (timespan*SAMPRATE);
  else
    d->volinc= (1.0/(d->opts->fortefactor*d->opts->fortefactor)-1.0) *
		    d->currvol / (timespan*SAMPRATE);
}



/*=============================================================================
    xyz "filename" "voice"

An alternative score format based on comma-separated lists.  The file contains
several voices which consist in a name (not containing white space) followed by
a block of numbers enclosed in braces {}.  The numbers are organised in rows
(which need not be equivalent to text lines).  Each row contains the same
number of values separated by commas, which will be output in the corresponding
output channels of this object.  Most rows are followed by a semicolon.  The
end of a bar is denoted by a '|' instead.  The first row can contain for each
channel: an optional keyword denoting special columns, a prefactor, the word
"mute" and the mute value.  The first row is followed by a colon.  The
prefactor is applied to all the numbers in the given column.  The keyword gives
the type of the data in the column.  The following types are available: "abc"
allows to give frequencies in abc-like notation (see below); the prefactor is
then a transposition factor.  "length" denotes the column giving the gross note
length and may be followed by two numbers: the unit length (metre) in seconds
and the bar length in units.  The keyword "net" denotes the column giving the
net note length as a fraction of the gross length.  The keyword "cue" causes
the channel to be reset to 0 after the given value has been output for one
sample.

If the word "mute" and a number follows after the prefactor, this
number is output on this channel when the voice is silent (that is, wenn there
is a 'z' in all abc columns or when a note has ended and before the next
starts).  The keyword "length" has to be present in one of the columns.  If the
bar length is given, the bars are checked to be of equal length.  Keywords must
precede the corresponding prefactors.  Unspecified prefactors default to 1, but
the commas between them have to be present.  If "mute" is present, the
prefactor has to be given explicitly.  Unspecified values in the data default
to the previous value (or zero if there was no precedent).  Values may also be
given as fractions of real numbers in the notation <numerator> '/' 
<denominator>.

The abc type columns can contain frequencies in a notation similar to abc.  The
',' denoting the double bass octave in abc is replaced by a '.' since the comma
is the column delimiter.  The prime giving the higher octaves is replaced by a
'`' to keep the preprocessor happy (see below).  The key is always C major,
accidentals have to be given for every note, and note lengths have to be given
separately in the length column.  Numerical values are also allowed in the abc
column.

The length and net length ratio columns are treated specially.  They are not 
output as data channels.

xyz runs the C preprocessor over the file before parsing it.  This means that
one can use preprocessor macros and C-style comments in the file.  It also 
makes it necessary to replace the prime in abc syntax by a '`' (see above)
because the preprocessor will otherwise treat it as the start of a character
constant.

See also: abc, at
=============================================================================*/

typedef struct {
  double    val;
  int	    flags;
}
xyzdata;

struct xyzparser {
  double    prefactor[PLENTY], current[PLENTY];
  float	    muteval[PLENTY];
  int	    type[PLENTY];
  const char *filename, *voice;
  char	    *buf, *read;
  int	    ncolumns, ndata;
  long	    line;
  xyzdata   *data, *write, *end;
  double    unitnote, bar, barsum;
};

struct xyz {
  double    target[PLENTY], slope[PLENTY], time[PLENTY];
  double    note;
  xyzdata   *buf, *pos, *end;
  int	    startup, startnote, firstval;
};

void xyzskipspace(struct xyzparser *pa);
void xyzskipcolumn(struct xyzparser *pa);
void xyzskiprow(struct xyzparser *pa);
void xyzparserow(struct xyz *d, struct xyzparser *pa);
double xyzparseabc(struct xyzparser *pa);
double xyzparsefraction(struct xyzparser *pa);
float calcxyz(sndobj *me, sndobj *caller, int innr);

#define XYZTYPE_NUMERIC	0
#define XYZTYPE_ABC	1
#define XYZTYPE_LENGTH	2
#define XYZTYPE_NET	3
#define XYZTYPE_CUE	4

#define XYZFLAGS_CUE	    1
#define XYZFLAGS_LINSTART   2
#define XYZFLAGS_LINEND	    4
#define XYZFLAGS_MUTED	    8

sndobj *xyz(const char *filename, const char *voice)
{
  sndobj *p;
  struct xyz *d;
  struct xyzparser pa;
  FILE *file;
  char *syscmd, *tmpfile, *read2, *objname;
  long length, nrows;
  int count, ncols, lengthtypes, lengthcol;
  char endc;

  length= strlen(filename);
  for( read2= (char*)filename+length-1; read2> filename && *read2!='/'; --read2 );
  if( *read2=='/' )
    ++read2;
  syscmd= malloc(2*length + 20);
  snprintf(syscmd, 2*length + 20, "cpp %s /tmp/%s.cpp", filename, read2);
  if( system(syscmd) ) {
    fprintf(stderr, "Error: xyz: Error while executing preprocessor! Aborting.\n");
    exit(1);
  }
  tmpfile= syscmd;
  snprintf(tmpfile, 2*length + 20, "/tmp/%s.cpp", read2);
  file= fopen(tmpfile, "r");
  free(syscmd);
  fseek(file, 0L, SEEK_END);
  length= ftell(file);
  fseek(file, 0L, SEEK_SET);
  pa.buf= malloc(length+2);
  fread(pa.buf, 1L, length, file);
  if( ferror(file) ) {
    fprintf(stderr, "Error: xyz: Error reading data in `%s'. Aborting.\n", filename );
    exit(1);
  }
  fclose(file);
  pa.buf[length]= pa.buf[length+1]= 0;
  pa.bar= 0.0;
  pa.barsum= 0.0;
  pa.voice= voice;
  pa.filename= filename;
  pa.read= pa.buf;
  pa.line= 0;
  for( pa.read= pa.buf; *pa.read; ++pa.read )
  {
    if( *pa.read=='\n' )
      ++pa.line;
    else if( *pa.read=='#' ) {
      ++pa.read;
      pa.line= strtol(pa.read, &pa.read, 10);
      while( *pa.read && *pa.read!='\n' )
        ++pa.read;
    }
    else if( *pa.read=='{' ) {
      for( read2= pa.read - 1; read2> pa.buf && isspace(*read2); --read2 );
      endc= read2[1];
      read2[1]= 0;
      while( read2> pa.buf && !isspace(*read2) )
        --read2;
      if( isspace(*read2) )
        ++read2;
      if( !strcmp(voice, read2) )
        break;
      for( read2= pa.read - 1; read2> pa.buf && isspace(*read2); --read2 );
      read2[1]= endc;
    }
  }
  if( !*pa.read ) {
    fprintf(stderr, "Error: xyz: cannot find voice `%s' in file `%s'.  Aborting.\n", voice, filename);
    exit(1);
  }
  ++pa.read;
  nrows= 0;
  ncols= 1;
  read2= pa.read;
  while( 13 )
  {
    if( *read2=='}' || !*read2 ) {
      ++nrows;
      break;
    }
    else if( *read2==';' || *read2=='|' || *read2==':' )
      ++nrows;
    else if( *read2==',' )
      if( !nrows )
        ++ncols;
    ++read2;
  }
  nrows += 2;
  if( !*read2 )
    fprintf(stderr, "Warning: xyz: Voice `%s' ends with end of file `%s'.\n", voice, filename);
  *read2= 0;
  //  Now parse the voice we want
  objname= news(char, NAMESTRLEN);
  if( strlen(filename)+strlen(voice) > NAMESTRLEN-3 ) {
    const char *filefrom, *voicefrom;
    voicefrom= voice + strlen(voice) - (int)(NAMESTRLEN * 0.7)-3;
    if( voicefrom< voice )
      voicefrom= voice;
    filefrom= filename + strlen(filename) - (NAMESTRLEN - (int)(NAMESTRLEN * 0.7) - 6);
    snprintf(objname, NAMESTRLEN, "...%s: ...%s", filefrom, voicefrom);
  }
  else
    snprintf(objname, NAMESTRLEN, "%s: %s", filename, voice);
  p= newsndo(calcxyz, objname, "xyz", 1, 0);
  p->skip= dontskip;
  p->private[0]= d= new(struct xyz);
  p->private[1]= d->buf= news(xyzdata, 2*nrows*ncols);
  d->end= d->buf + nrows*ncols;
  p->private[2]= objname;
  xyzskipspace(&pa);
  for( read2= pa.read; *read2 && *read2!=':'; ++read2 );
  if( !*read2 ) {
    fprintf(stderr, "Warning: xyz: voice `%s' in file `%s' has no explicit prefactors and types (line %ld).  Assuming first column is length.\n", pa.voice, pa.filename, pa.line);
    for( count= 0; count< PLENTY; ++count ) {
      pa.prefactor[count]= 1.0;
      pa.type[count]= XYZTYPE_NUMERIC;
      pa.muteval[count]= MAGIC;
    }
    pa.type[0]= XYZTYPE_LENGTH;
    for( read2= pa.read; *read2 && *read2!=';' && *read2 !='|'; ++read2 )
      if( *read2==',' )
        ++p->nch;
    p->nch -= 2;
    if( !p->nch ) {
      fprintf(stderr, "Warning: xyz: No output channels remain in voice `%s' in file `%s', line %ld.  Returning zero.\n", pa.voice, pa.filename, pa.line);
      return c(0);
    }
  }
  else {
    count= 0;
    lengthtypes= 0;
    p->nch= 0;
    do {
      xyzskipspace(&pa);
      if( isalpha(*pa.read) ) {
        if( !strncmp(pa.read, "abc", 3L) && !isalnum(pa.read[3]) ) {
	  pa.type[count]= XYZTYPE_ABC;
	  pa.read += 3;
	}
        else if( !strncmp(pa.read, "cue", 3L) && !isalnum(pa.read[3]) ) {
	  pa.type[count]= XYZTYPE_CUE;
	  pa.read += 3;
	}
        else if( !strncmp(pa.read, "length", 6L) && !isalnum(pa.read[6]) ) {
	  if( (lengthtypes&1)!=0 ) {
	    fprintf(stderr, "Warning: xyz: \"length\" given as type more than once in column %d, line %ld in voice `%s' of `%s'.  Ignored.\n", count, pa.line, pa.voice, pa.filename);
	    pa.type[count]= XYZTYPE_NUMERIC;
	  }
	  else {
	    pa.type[count]= XYZTYPE_LENGTH;
	    lengthtypes |= 1;
	  }
	  pa.read += 6;
	}
        else if( !strncmp(pa.read, "net", 3L) && !isalnum(pa.read[3]) ) {
	  if( (lengthtypes&2)!=0 ) {
	    fprintf(stderr, "Warning: xyz: \"net\" given as type more than once in column %d, line %ld in voice `%s' of `%s'.  Ignored.\n", count, pa.line, pa.voice, pa.filename);
	    pa.type[count]= XYZTYPE_NUMERIC;
	  }
	  else {
	    pa.type[count]= XYZTYPE_NET;
	    lengthtypes |= 2;
	  }
	  pa.read += 3;
	}
	else {
	  fprintf(stderr, "Warning: xyz: Unknown type `");
	  while( isalnum(*pa.read) )
	    fprintf(stderr, "%c", (int)*pa.read++);
	  fprintf(stderr, "' in column %d of voice `%s' in `%s', line %ld.  Substituting default.\n", count, pa.voice, pa.filename, pa.line );
	  pa.type[count]= XYZTYPE_NUMERIC;
	}
	xyzskipspace(&pa);
      }
      else
        pa.type[count]= XYZTYPE_NUMERIC;

      if( *pa.read=='-' || *pa.read=='.' || isdigit(*pa.read) ) {
        pa.prefactor[count]= xyzparsefraction(&pa);
	xyzskipspace(&pa);
	if( pa.type[count]==XYZTYPE_LENGTH && (isdigit(*pa.read) || *pa.read=='.') ) {
  	  pa.bar= strtod(pa.read, &pa.read) * pa.prefactor[count];
	  xyzskipspace(&pa);
	}
      }
      else
        pa.prefactor[count]= 1.0;

      if( !strncmp(pa.read, "mute", 4L) && !isalnum(pa.read[4]) ) {
        pa.read += 4;
	xyzskipspace(&pa);
	if( *pa.read=='-' || *pa.read=='.' || isdigit(*pa.read) ) {
	  pa.muteval[count]= strtod(pa.read, &pa.read) * pa.prefactor[count];
	  xyzskipspace(&pa);
	}
	else
	  pa.muteval[count]= 0.0;
	if( pa.type[count]==XYZTYPE_LENGTH || pa.type[count]==XYZTYPE_NET ) {
	  fprintf(stderr, "Warning: xyz: Column %d, line %ld of voice `%s' in `%s': length and net type columns cannot be muted.  Ignoring mute value.\n", count, pa.line, pa.voice, pa.filename);
	  pa.muteval[count]= MAGIC;
	}
      }
      else
        pa.muteval[count]= MAGIC;   // no muting

      if( *pa.read!=',' && *pa.read!=':' ) {
        fprintf(stderr, "Warning: xyz: Surplus type or illegal character found at column %d, line %ld of voice `%s' in `%s'.\n", count, pa.line, pa.voice, pa.filename);
	xyzskipcolumn(&pa);
      }
      if( pa.type[count]!=XYZTYPE_LENGTH && pa.type[count]!=XYZTYPE_NET )
        ++p->nch;
      ++count;
    }
    while( *pa.read++!=':' && count< PLENTY );
    if( pa.read[-1]!=':' ) {
      fprintf(stderr, "Warning: xyz: Too many columns (> %d) in type row in voice `%s' in `%s', line %ld.  Ignoring surplus columns.\n", (int)PLENTY, pa.voice, pa.filename, pa.line);
      xyzskiprow(&pa);
      if( *pa.read )
        ++pa.read;
    }
    if( !p->nch ) {
      fprintf(stderr, "Warning: xyz: No output channels remain in voice `%s' in file `%s', line %ld.  Returning zero.\n", pa.voice, pa.filename, pa.line);
      return c(0);
    }
    pa.ncolumns= count;
    xyzskipspace(&pa);
  }
  for( lengthcol= 0; lengthcol< PLENTY; ++lengthcol )
    if( pa.type[lengthcol]==XYZTYPE_LENGTH )
      break;
  if( lengthcol==PLENTY ) {
    fprintf(stderr, "Warning: xyz: No length column found in voice `%s' in `%s' (line %ld).  Choosing first column.\n", pa.voice, pa.filename, pa.line);
    pa.type[0]= XYZTYPE_LENGTH;
    lengthcol= 0;
  }
  pa.unitnote= pa.prefactor[lengthcol];
  for( count= 0; count< pa.ncolumns; ++count )
    if( pa.type[count]==XYZTYPE_LENGTH )
      pa.current[count]= pa.unitnote;
    else if( pa.type[count]==XYZTYPE_NET )
      pa.current[count]= 1.0;
    else
      pa.current[count]= 0.0;

  for( count= 0; count< PLENTY; ++count ) {
    d->target[count]= 0.0;
    d->slope[count]= 0.0;
    d->time[count]= 0.0;
  }
  d->pos= d->buf;
  pa.data= pa.write= d->buf;
  pa.end= d->end;
  pa.ndata= p->nch + 1;

  while( *pa.read )
    xyzparserow(d, &pa);
  if( pa.write + p->nch+1 <= pa.end ) {    //  should be true
    pa.write->val= DBL_MAX;
    pa.write->flags= 0;
    ++pa.write;
    for( count= 0; count< pa.ncolumns; ++count ) {
      if( pa.type[count]==XYZTYPE_LENGTH || pa.type[count]==XYZTYPE_NET )
	continue;
      if( pa.muteval[count]!=MAGIC )
        pa.write->val= pa.muteval[count];
      else
        pa.write->val= pa.current[count];
      if( pa.type[count]==XYZTYPE_CUE ) {
        pa.write->flags= XYZFLAGS_CUE;
	if( pa.write - pa.ndata < pa.data || 
			pa.write->val == pa.write[-pa.ndata].val )
	  pa.write->val= 0.0;
      }
      else
        pa.write->flags= 0;
      ++pa.write;
    }
  }
  
  free(pa.buf);
  d->startup= 1;
  return p;
}


//  Skips white space, keeping line number up to date
void xyzskipspace(struct xyzparser *pa)
{
  while( isspace(*pa->read) )
    if( *pa->read=='\n' ) {
      ++pa->read;
      if( *pa->read=='#' ) {
        ++pa->read;
        pa->line= strtol(pa->read, &pa->read, 10);  // number of following line
        while( *pa->read && *pa->read!='\n' )
          ++pa->read;
	++pa->read;
      }
      else
        ++pa->line;
    }
    else
      ++pa->read;
}


//  Skips until end of column (but does not skip over end-of-column character)
void xyzskipcolumn(struct xyzparser *pa)
{
  while( *pa->read && *pa->read!=',' && *pa->read!=';' && *pa->read!=':' && *pa->read!='|' )
    if( *pa->read=='\n' ) {
      ++pa->read;
      if( *pa->read=='#' ) {
        ++pa->read;
        pa->line= strtol(pa->read, &pa->read, 10);  // number of following line
        while( *pa->read && *pa->read!='\n' )
          ++pa->read;
	++pa->read;
      }
      else
        ++pa->line;
    }
    else
      ++pa->read;
}


//  Skips until end of row (but does not skip over end-of-row character)
void xyzskiprow(struct xyzparser *pa)
{
  while( *pa->read && *pa->read!=':' && *pa->read!=';' && *pa->read!='|' )
    if( *pa->read=='\n' ) {
      ++pa->read;
      if( *pa->read=='#' ) {
        ++pa->read;
        pa->line= strtol(pa->read, &pa->read, 10);  // number of following line
        while( *pa->read && *pa->read!='\n' )
          ++pa->read;
	++pa->read;
      }
      else
        ++pa->line;
    }
    else
      ++pa->read;
}


void xyzparserow(struct xyz *d, struct xyzparser *pa)
{
  int linmode[PLENTY];
  double note, mutenote;
  int col, allrest, prev_endcue;
  
  for( col= 0; col< pa->ncolumns && *pa->read; ++col )
  {
    xyzskipspace(pa);
    linmode[col]= 0;
    if( *pa->read==')' ) {
      linmode[col] |= XYZFLAGS_LINEND;
      ++pa->read;
    }
    if( *pa->read=='(' ) {
      linmode[col] |= XYZFLAGS_LINSTART;
      ++pa->read;
    }
    if( *pa->read==',' || *pa->read=='|' || *pa->read==';' || !*pa->read )
      ;
    else if( pa->type[col]==XYZTYPE_ABC && *pa->read!='-' && *pa->read!='.' && !isdigit(*pa->read) )
    {
      pa->current[col]= xyzparseabc(pa) * pa->prefactor[col];
      if( pa->current[col]< 0.0 ) {
        fprintf( stderr, "Warning: xyz: Illegal note symbol `%c' in column %d, line %ld in voice `%s' in `%s'. Interpreted as rest.\n", *pa->read, col, pa->line, pa->voice, pa->filename ); 
	pa->current[col]= pa->muteval[col];
      }
      else if( pa->current[col] == 0.0 && pa->muteval[col]!=MAGIC )
	pa->current[col]= pa->muteval[col];
    }
    else
      pa->current[col]= xyzparsefraction(pa) * pa->prefactor[col];
    xyzskipcolumn(pa);
    if( pa->type[col]==XYZTYPE_LENGTH )
      note= pa->current[col];
    else if( pa->type[col]==XYZTYPE_NET ) {
      mutenote= 1.0 - pa->current[col];
      if( mutenote< 0.0 )
        mutenote= 0.0;
      else if( mutenote> 1.0 )
        mutenote= 1.0;
    }
    if( *pa->read=='|' || *pa->read==';' || !*pa->read )
      break;
    ++pa->read;
  }
  mutenote *= note;
  pa->barsum += note;
  xyzskipspace(pa);
  if( *pa->read=='|' ) {
    if( fabs(pa->barsum - pa->bar)> 1e-5 )
      fprintf(stderr, "Warning: xyz: Bar ending on line %ld of voice `%s' in `%s' deviates from standard length by %g unit notes (%g seconds).\n", pa->line, pa->voice, pa->filename, (pa->barsum-pa->bar)/pa->unitnote, pa->barsum-pa->bar );
    pa->barsum= 0.0;
  }
  if( *pa->read ) {
    ++pa->read;
    xyzskipspace(pa);
  }
  
  allrest= 0;
  for( col= 0; col< pa->ncolumns; ++col )
    if( pa->type[col]==XYZTYPE_ABC )
      if( pa->current[col]==pa->muteval[col] )
        allrest= 1;
      else {
        allrest= 0;
	break;
      }

  if( pa->write >= pa->end )
    return;

  if( allrest )
    mutenote= note;

  prev_endcue= pa->write<=pa->data || 
		(pa->write[-pa->ndata].flags & XYZFLAGS_MUTED)!=0;

  note -= mutenote;
  pa->write->val= note;
  pa->write->flags= 0;
  ++pa->write;
  for( col= 0; col< pa->ncolumns && pa->write< pa->end; ++col ) {
    if( pa->type[col]==XYZTYPE_LENGTH || pa->type[col]==XYZTYPE_NET )
	continue;
    pa->write->val= pa->current[col];
    pa->write->flags= linmode[col];
    if( pa->type[col]==XYZTYPE_CUE ) {
      pa->write->flags |= XYZFLAGS_CUE;
      if( note==0 )
	pa->write->val= 0.0;
    }
    ++pa->write;
  }

  if( mutenote> 0.0 ) {
    pa->write->val= mutenote;
    pa->write->flags= XYZFLAGS_MUTED;
    ++pa->write;
    for( col= 0; col< pa->ncolumns && pa->write< pa->end; ++col ) {
      if( pa->type[col]==XYZTYPE_LENGTH || pa->type[col]==XYZTYPE_NET )
  	continue;
      if( pa->muteval[col]!=MAGIC )
  	pa->write->val= pa->muteval[col];
      else
  	pa->write->val= pa->current[col];
      pa->write->flags= XYZFLAGS_MUTED;
      if( pa->type[col]==XYZTYPE_CUE ) {
        pa->write->flags |= XYZFLAGS_CUE;
	if( !note && prev_endcue )
	  pa->write->val= 0.0;
      }
      ++pa->write;
    }
  }
}


//  Returns frequency of abc note or 0 for rest.  -1 means illegal note symbol.
double xyzparseabc(struct xyzparser *pa)
{
  static double notefreqfact[]= { 1.0, 1.12246204830937298142, 
    0.59460355750136053336, 0.66741992708501718241, 0.74915353843834074940, 
    0.79370052598409973737, 0.89089871814033930474 };

  double accfactor, note;
  
  xyzskipspace(pa);
  accfactor= 1.0;
  if( *pa->read=='=' )
    ++pa->read;
  else if( *pa->read=='^' ) {
    accfactor= SEMITONE;
    ++pa->read;
  }
  else if( *pa->read=='_' ) {
    accfactor= 1.0/SEMITONE;
    ++pa->read;
  }
  
  if( *pa->read>='A' && *pa->read<='G' ) {
    note= 110.0 * accfactor * notefreqfact[*pa->read-'A'];
    ++pa->read;
    if( *pa->read=='.' )
      for( ; *pa->read=='.'; ++pa->read )
        note /= 2.0;
    return note;
  }
  else if( *pa->read>='a' && *pa->read<='g' ) {
    note= 220.0 * accfactor * notefreqfact[*pa->read-'a'];
    ++pa->read;
    if( *pa->read=='`' )
      for( ; *pa->read=='`'; ++pa->read )
        note *= 2.0;
    return note;
  }
  else if( *pa->read=='z' )
    return 0.0;
  else
    return -1.0;
}


double xyzparsefraction(struct xyzparser *pa)
{
  double denom, num;
  
  xyzskipspace(pa);
  if( *pa->read=='/' )
    num= 1.0;
  else
    num= strtod(pa->read, &pa->read);
  if( *pa->read!='/' )
    return num;
  ++pa->read;
  if( !*pa->read )
    return num;
  denom= strtod(pa->read, &pa->read);
  if( denom==0.0 )
    denom= 2.0;
  return num/denom;
}


float calcxyz(sndobj *me, sndobj *caller, int innr)
{
  struct xyz *d;
  xyzdata *searchend;
  double lintime;
  int ch;
  
  d= (struct xyz *)me->private[0];
  if( d->startup ) {
    //  Skip until starttime:
    double time= 0;
    while( 13 ) {
      time += d->pos->val;
      if( time >= starttime )
  	break;
      d->pos += me->nch+1;
    }
    d->pos -= me->nch + 1;
    d->note= -(time - starttime);
    d->startnote= 1;
    d->startup= 0;
  }
  if( d->note> 0.0 ) {
    if( d->firstval ) {
      d->firstval= 0;
      FORCH
        if( (d->pos[ch+1].flags & XYZFLAGS_CUE)!=0 ) {
	  OUTPUT(ch, 0.0);
	}
    }
    FORCH
      if( d->time[ch] > 0.0 ) {
        if( (d->pos[ch+1].flags & XYZFLAGS_MUTED)==0 )
	  OUTPUT(ch, d->target[ch] - d->time[ch]*d->slope[ch]);
	d->time[ch] -= 1.0/SAMPRATE;
      }
  }
  else
    while( d->pos + me->nch+1 < d->end && d->note <= 0.0 )
    {
      d->pos += me->nch+1;
      d->firstval= 1;
      if( d->startnote ) {
        d->note= -d->note;
	d->startnote= 0;
      }
      else
        d->note= d->pos[0].val;
#ifdef XYZ_DEBUG
      printf("xyz outputs ");
#endif
      FORCH  {
  	if( d->time[ch] > 0.0 && (d->pos[ch+1].flags & XYZFLAGS_MUTED)==0  ) {
  	  OUTPUT(ch, d->target[ch] - d->time[ch]*d->slope[ch]);
  	  d->time[ch] -= 1.0/SAMPRATE;
  	}
  	else
	  OUTPUT(ch, d->pos[ch+1].val);
#ifdef XYZ_DEBUG
	printf("%g  ", (double)me->ch[ch]);
#endif
	if( (d->pos[ch+1].flags & XYZFLAGS_LINSTART) != 0 ) {
	  lintime= d->note;
	  for( searchend= d->pos+me->nch+ch+2; searchend< d->end; searchend+=me->nch+1 ) {
	    if( (searchend->flags & XYZFLAGS_LINEND) != 0 )
	      break;
	    lintime += searchend[-ch-1].val;
	  }
	  if( searchend< d->end ) {
  	    d->target[ch]= searchend->val;
  	    d->slope[ch]= (searchend->val - d->pos[ch+1].val)/lintime;
  	    d->time[ch]= lintime - 1.0/SAMPRATE;
	  }
	  else
	    d->time[ch]= 0.0;
	}
      }
#ifdef XYZ_DEBUG
      printf(" length %g\n", (double)d->note);
#endif
    }

  d->note -= 1.0/(double)SAMPRATE;

  CALCRETURN;
}



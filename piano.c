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

#ifdef NO_X
#error piano needs X11!  Make sure X is available and comment "#define NO_X" in sndsys.h.
#endif

#ifdef NO_OSS
#error piano needs Open Sound System!  Make sure OSS is available and comment "#define NO_OSS" in sndsys.h.
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/soundcard.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/keysym.h>

#include "sndsys.h"

//  X Bitmap for keyboard image:
#include "piano.xbm"

const int keys[]= { XK_z, XK_s, XK_x, XK_d, XK_c, XK_v, XK_g, XK_b, XK_h, XK_n, XK_j, XK_m, XK_comma, XK_l, XK_period, XK_semicolon, XK_slash, XK_q, XK_2, XK_w, XK_3, XK_e, XK_4, XK_r, XK_t, XK_6, XK_y, XK_7, XK_u, XK_i, XK_9, XK_o, XK_0, XK_p, XK_minus, XK_bracketleft, XK_bracketright, -1 };

const int shiftkeys[]= { XK_Z, XK_S, XK_X, XK_D, XK_C, XK_V, XK_G, XK_B, XK_H, XK_N, XK_J, XK_M, XK_less, XK_L, XK_greater, XK_colon, XK_question, XK_Q, XK_at, XK_W, XK_numbersign, XK_E, XK_dollar, XK_R, XK_T, XK_asciicircum, XK_Y, XK_ampersand, XK_U, XK_I, XK_parenleft, XK_O, XK_parenright, XK_P, XK_underscore, XK_braceleft, XK_braceright, -1 };

#define PIANO_VOL       0.7
#define PIANO_AMPL      0.5
#define PIANO_DECAY     0.5
#define PIANO_THRESHOLD 0.001

typedef struct {
  float     freq, trigger;
  int       freqflag, gate, shiftgate;
  sndobj    *generator, *amplitude;
}
pianokey;

int totalkeys;
pianokey *piano;
pianokey *keylookup[256];
int prerun;
sndobj *dummymike;
sndobj *top;

void eventloop(int argc, char **argv);
void setup( void );
float addpiano(sndobj *me, sndobj *caller, int innr);
void tunechromatic( void );

int main(int argc, char *argv[])
{
  setup();
  eventloop(argc, argv);
  killemall();
  free(piano);
  return 0;
}


void eventloop(int argc, char **argv)
{
  Display *display;
  Window win;
  XTextProperty name;
  XSizeHints *size_hints;
  XWMHints *wm_hints;
  XClassHint *class_hints;
  XImage *keyboardimg;
  GC graphcontext;
  XEvent event;
  KeySym key;
  Atom protocol[2], message;
  pianokey *hitkey;
  char *namestr;
  int screen, nevents, quit, count;
  
  display= XOpenDisplay(NULL);
  if( !display ) {
    fprintf(stderr, "piano: Error: could not connect to X server!\n");
    exit(1);
  }
  screen= DefaultScreen( display );
  win= XCreateSimpleWindow( display, RootWindow(display, screen),
	10,10,piano_width,piano_height,
	4, BlackPixel(display, screen), WhitePixel(display, screen) );
  size_hints= XAllocSizeHints();
  size_hints->flags = PSize | PMinSize;
  size_hints->min_width = piano_width;
  size_hints->min_height = piano_height;
  wm_hints= XAllocWMHints();
  wm_hints->flags = StateHint | InputHint;
  wm_hints->initial_state= NormalState;
  wm_hints->input= True;
  class_hints= XAllocClassHint();
  class_hints->res_name= "piano";
  class_hints->res_class= "piano";
  namestr= "Piano";
  XStringListToTextProperty( &namestr, 1, &name );
  XSetWMProperties(display, win, &name, &name, argv, argc, size_hints, wm_hints, class_hints );
  XSelectInput(display, win, KeyPressMask|KeyReleaseMask|FocusChangeMask|ExposureMask );
  protocol[0]= XInternAtom(display, "WM_PROTOCOLS", False );
  protocol[1]= XInternAtom(display, "WM_DELETE_WINDOW", False );
  XSetWMProtocols( display, win, protocol, 2 );
  XMapWindow(display, win);
  keyboardimg= XCreateImage( display, DefaultVisual(display, screen), 1, XYBitmap, 0, (char*)piano_bits, piano_width, piano_height, 8, 0);
  graphcontext= XCreateGC(display, win, 0, NULL);
  XSetBackground(display, graphcontext, WhitePixel(display, screen));
  XSetForeground(display, graphcontext, BlackPixel(display, screen));
  quit= 0;
  do {
    for( nevents= XPending(display); nevents> 0; --nevents ) {
      XNextEvent( display, &event );
      if( event.type==KeyPress ) {
	key= XLookupKeysym( &event.xkey, 0 );
	if( key==XK_Escape )
	  quit= 1;
	if( (unsigned)key< 256 && NULL!=(hitkey=keylookup[key&0xFF]) ) {
	  hitkey->gate= 1;
	  hitkey->trigger= PIANO_AMPL;
	  hitkey->shiftgate= event.xkey.state & (ShiftMask|LockMask);
	}
	else if( key==XK_Shift_L || key==XK_Shift_R )
	  for( count= 0; count< totalkeys; ++count )
	    piano[count].shiftgate |= piano[count].gate;
	else if( key==XK_Caps_Lock || key==XK_Shift_Lock ) {
	  if( (event.xkey.state & LockMask)==0 )
	    for( count= 0; count< totalkeys; ++count )
	      piano[count].shiftgate |= piano[count].gate;
	  else if( (event.xkey.state & ShiftMask)==0 )
	    for( count= 0; count< totalkeys; ++count )
	      piano[count].shiftgate= 0;
	}
      }
      else if( event.type==KeyRelease ) {
	key= XLookupKeysym( &event.xkey, 0 );
	if( (unsigned)key< 256 && NULL!=(hitkey=keylookup[key&0xFF]) )
	  hitkey->gate= 0;
	else if( key==XK_Shift_L || key==XK_Shift_R ) {
	  if( (event.xkey.state & LockMask)==0 )
	    for( count= 0; count< totalkeys; ++count )
	      piano[count].shiftgate= 0;
	}
      }
      else if( event.type==FocusOut )
	XAutoRepeatOn(display);
      else if( event.type==FocusIn )
	XAutoRepeatOff(display);
      else if( event.type==ClientMessage ) {
	if( event.xclient.message_type==protocol[0] ) {
	  if( event.xclient.format==8 )
	    message= event.xclient.data.b[0];
	  else if( event.xclient.format==16 )
	    message= event.xclient.data.s[0];
	  else if( event.xclient.format==32 )
	    message= event.xclient.data.l[0];
	  else message= protocol[1] - 1;
	  if( message==protocol[1] )
	    quit= 1;
	}
      }
      else if( event.type==Expose && event.xexpose.count==0 )
      	XPutImage(display, win, graphcontext, keyboardimg, 0,0, 0,0, piano_width, piano_height);
    }
    top->calc(top, NULL, 0);
  }
  while( !quit );
  XAutoRepeatOn(display);
  XFreeGC(display, graphcontext);
  XFree(size_hints);
  XFree(wm_hints);
  XFree(class_hints);
  XCloseDisplay(display);
}


void setup( void )
{
  int count, keynr, fillbuf;
  
  for( totalkeys= 0; keys[totalkeys]>=0; ++totalkeys );
  for( ; shiftkeys[totalkeys]>=0; ++totalkeys );
  piano= news( pianokey, totalkeys );
  for( count= 0; count< totalkeys; ++count ) {
    piano[count].freq= 0.0;
    piano[count].trigger= 0.0;
    piano[count].freqflag= 0;
    piano[count].gate= 0;
    piano[count].shiftgate= 0;
    piano[count].generator= file("pianostring.asc", extc( 1, &piano[count].freq, &piano[count].freqflag ));
    piano[count].amplitude= expodecay( c(PIANO_DECAY), exttrigger( &piano[count].trigger ) );
//    piano[count].generator= sine( 0.0, extc( 1, &piano[count].freq, &piano[count].freqflag ), c(1.0) );
  }
//  makefanouts();  not necessary here
  for( count= 0; count< 256; ++count )
    keylookup[count]= NULL;
  for( keynr= 0; keys[keynr]>=0; ++keynr )
    keylookup[0xFF&keys[keynr]]= piano+keynr;
  for( keynr= 0; shiftkeys[keynr]>=0; ++keynr )
    keylookup[0xFF&shiftkeys[keynr]]= piano+keynr;
  top= oss( PIANO_VOL, lowpass(c(500), newsndo( addpiano, "addpiano", "addpiano", 1, 0) ) );
  tunechromatic();
  dummymike= mike( 0, -1.0 );   // needed with my soundcard to reduce latency
  prerun= 1;        //  prerun to fill output buffer
  ioctl( dspdesc(0), SNDCTL_DSP_GETBLKSIZE, &fillbuf );
  for( fillbuf /= 4; fillbuf> 0; --fillbuf )
    top->calc(top, NULL, 0);
  prerun= 0;
}

float addpiano(sndobj *me, sndobj *caller, int innr)
{
  float outval, ampl;
  int count;
  
  if( prerun ) {
    OUTPUT(0, 0.0);
    CALCRETURN;
  }
  for( count= 0, outval= 0.0; count< totalkeys; ++count )
    if( piano[count].gate || piano[count].shiftgate ) {
      ampl= piano[count].amplitude->calc( piano[count].amplitude, me, 0 );
      outval += ampl * piano[count].generator->calc( piano[count].generator, me, 0 );
      if( ampl < PIANO_THRESHOLD )
	piano[count].gate= piano[count].shiftgate= 0;
    }
  OUTPUT(0, outval );
  dummymike->calc(dummymike, NULL, 0);  // needed with my soundcard to reduce latency
  CALCRETURN;
}

void tunechromatic( void )
{
  float afreq= 110.0, multiplier= pow( 2, 3.0/12.0 ), factor= pow(2, 1.0/12.0);
  int count;
  
  for( count= 0; count< totalkeys; ++count ) {
    if( multiplier> 1.99 ) {
      afreq *= 2;
      multiplier= 1.0;
    }
    piano[count].freq= afreq * multiplier;
    piano[count].freqflag= 1;
    multiplier *= factor;
  }
}





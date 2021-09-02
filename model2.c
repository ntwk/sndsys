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


#include "sndsys.h"

abcplayopts abcopts= {
    1, 1, 0, 1, 1,
    1.0, 1.2, 2.0,
    1, 1, 1,
    0.8, 1.0, 0.0, 0.8, 0.0, 
    {  },
    {  }
};

double trigger0[]= { 0.0, 1, END };
double gate[]= { 0.0, 1, RAMP_C, 1e-3, 0.0, END };

sndobj *plucked(sndobj *score);

sndobj *plucked(sndobj *score)
{
  static int loopnr= 0;
  sndobj *cue, *vol, *damping, *freq, *model;
  
  vol= c1(score);
  cue= mul(2, c0(score), vol);
  damping= c(1);//conditional(vol, c(1e-4), c(1), c(0));
  freq= c2(score);//mul(2, c2(score), damping);
  model= defloop("plucked", loopnr, delay(0.1,arithmetic("/", c(1.0), freq), 
	    add( 2, mul(2, holdtrigger(cue, c(1e-3)), linnoise(c(5000),c(1))),
		 _12_0pole( 1, 0.5, 0.5, 0.0, mul(2, damping,
		loop("plucked", loopnr, 1 ))))));
  ++loopnr;
//  return switchboard("0 1 2 3", model, cue, vol, freq);
  return model;
}




int main(int argc, char *argv[])
{
  sndobj *output;
  
  output= linear(0, 0.7, plucked(abc("melody2.abc", 2, &abcopts)));
  sndexecute( 0, 15, writeout("model2.wav", output));
}


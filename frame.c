/*
    This file is part of sndsys, a digital signal processing system.
    Copyright (c) 2007 Volker Schatz (noise at volkerschatz dot com).

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


int main(int argc, char *argv[])
{
  sndobj *output;
  double duration;
  
  output= c(0);
  duration= 2.0;
  sndexecute( 0, duration, writeout("frame.wav", output));
  return 0;
}


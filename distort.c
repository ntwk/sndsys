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

#ifdef NO_OSS
#error distort needs Open Sound System!  Make sure OSS is available and comment "#define NO_OSS" in sndsys.h.
#endif

#include <signal.h>
#include "sndsys.h"

int main(int argc, char *argv[])
{
  sndobj *top, *pickup;
  
  pickup= linear(0, 10,  mike(0, 1.0));
//  pickup= shelve( c(15000), c(0.01), linear(0, 50,  mike(0, 1.0)));
//  pickup= avg( c(15000), c(0), linear(0, 50,  mike(0, 1.0)));
  top= oss( 1, clip( 0.0, linear(0, 0.8, fastrms( c(300), pickup)), pickup ) );
  closeloops();
  prune(top);
  makefanouts();
  dumptree(top);
  while( 1 )
    top->calc(top, NULL, 0);
  return 0;
}


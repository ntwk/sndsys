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

double trigger0[]= { 0.0, 1, END };

int main(int argc, char *argv[])
{
  sndobj *model;
  
  model= defloop("karplus-strong", 1, cdelay( 1.0/440, 
	      add( 2, at( 0.0, trigger0 ), _12_0pole( 1, 0.5, 0.5, 0.0, 
		  loop("karplus-strong", 1, 1 )))));
  sndexecute( 0, 5, writeout("model.wav", model));
}


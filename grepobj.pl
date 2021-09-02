#!/usr/bin/perl

#    This file is part of sndsys, a digital signal processing system.
#    Copyright (c) 2004-06 Volker Schatz (noise at volkerschatz dot com).
#
#    sndsys is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    sndsys is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with sndsys; if not, see <http://www.gnu.org/licenses/>.

use strict;

my $objname;
my $nargs= @ARGV;
my $regex= $ARGV[0];

if ( $nargs!=1 || $ARGV[0] eq "-h" || $ARGV[0] eq "--help" ) {
  print "usage: grepobj <regular expression>\n";
  print "Searches snd*.c files for sndobj's the documentation of which contains <regular expression>.\n";
  exit 0;
}

while( my $filename= glob("snd*.c") )
{
  next if $filename eq "sndsys.c";	# not a file which defines sndobj's
  open( file, "<", $filename);
  while( <file> )
  {
    if( /^\/\*==========/ ) {
      if( defined( $_= <file> ) && ! /=======\*\// && /^\s+(\w+)\W+/ )
      {
	$objname= $1;
	while( defined( $_=<file> ) && ! /=======\*\// ) {
	  print "$objname: $_" if /$regex/;
	}
      }
      else {
	while( defined( $_=<file> ) && ! /=======\*\// ) { }
      }
      if( defined($_) )
	{ redo; }	  # so we can still scan the line after ====*/
    }
  }
}


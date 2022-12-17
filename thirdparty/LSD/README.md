LSD - Line Segment Detector
===========================

Version 1.5 - December 3, 2010

by Rafael Grompone von Gioi <grompone@gmail.com>


Introduction
------------

LSD is an implementation of the Line Segment Detector on digital
images described in the paper:

  "LSD: A Fast Line Segment Detector with a False Detection Control"
  by Rafael Grompone von Gioi, Jeremie Jakubowicz, Jean-Michel Morel,
  and Gregory Randall, IEEE Transactions on Pattern Analysis and
  Machine Intelligence, vol. 32, no. 4, pp. 722-732, April, 2010.

and in more details in the CMLA Technical Report:

  "LSD: A Line Segment Detector, Technical Report",
  by Rafael Grompone von Gioi, Jeremie Jakubowicz, Jean-Michel Morel,
  Gregory Randall, CMLA, ENS Cachan, 2010.

The version implemented here includes some further improvements
described on the LSD page at www.ipol.im. That same page includes more
information, including this code and an online demo version:

  http://www.ipol.im/pub/algo/gjmr_line_segment_detector


Files
-----

README.txt          - This file.
COPYING             - GNU AFFERO GENERAL PUBLIC LICENSE Version 3.
Makefile            - Compilation instructions for 'make'.
lsd.c               - LSD module ANSI C code.
lsd.h               - LSD module ANSI C header.
lsd_cmd.c           - LSD command line interface, ANSI C code.
lsd_call_example.c  - Minimal example of calling LSD from a C language program.
chairs.pgm          - Test image in PGM format.
chairs.lsd.txt      - Expected result for 'chairs.pgm' image as an ASCII file.
chairs.lsd.eps      - Expected result for 'chairs.pgm' image as an EPS file.
doc                 - Html code documentation.
doxygen.config      - doxygen configuration file for documentation generation.


Compiling
---------

LSD is an ANSI C Language program and can be used as a module
to be called from a C language program or as an independent
command.

In the distribution is included a Makefile file with instructions
to build the command lines program 'lsd', as well as minimal
example program on how to call LSD from C code.

To build both programs, a C compiler (called with 'cc') must be
installed on your system, as well as the program 'make'.
LSD only uses the standard C library so it should compile
in any ANSI C Language environment. In particular, it should
compile in an Unix like system.

The compiling instruction is just

  make

from the directory where the source codes and the Makefile are located.

To verify a correct compilation you can apply LSD to the test
image 'chairs.pgm' and compare the result to the provided ones.

An explicit example of how to compile a program using LSD as a module
is provided. The compilation line for 'lsd_call_example.c' is just

  cc -lm -o lsd_call_example lsd_call_example.c lsd.c


Running LSD Command
-------------------

The simplest LSD command execution is just

  lsd

or

  ./lsd

if the command is not in the path. That should print LSD version
and the command line interface, including the available options.
The only input image format handled by LSD is PGM, in its two
versions, ASCII and Binary. A useful execution would be:

  lsd chairs.pgm chairs.result.txt

That should give the result as an ASCII file 'chairs.result.txt'
with the coordinates each line segment detected as a line in
the file like the following:

  159.232890 134.369601 160.325338 105.613616 2.735466 

which means that a line segment starting at point (159.232890,134.369601)
and ending at point (160.325338 105.613616) and of width 2.735466
was detected. The unit is the pixel and the origin of coordinates
is the center of pixel (0,0).

For easier visualization of the result, the LSD command can also
give the output in EPS or SVG file formats. For example,

  lsd -P chairs.result.eps chairs.pgm chairs.result.txt

will, in addition to the ASCII output file, produce the EPS file
'chairs.result.eps'.

To see the full options, execute LSD command without parameters,
as in './lsd'.

Optional arguments should always appear before the needed arguments
input and output. For example, the following line is wrong:

  lsd chairs.pgm -s 0.5 chairs.result.txt   -> WRONG!!

and should be

  lsd -s 0.5 chairs.pgm chairs.result.txt

If the name of an input file is just - (one dash), then that
file will be read from the standard input. Analogously, if the
name of an output file is just - (one dash), then that file
will be written to the standard output. For example,

  lsd - -

will work as a filter, taking the input from standard input and
giving the output to standard output.


Code Documentation
------------------

There is a HTML documentation of the code on the directory 'doc'. The
entry point is the file 'doc/index.html' that should be opened with a
web browser. The documentation was automatically generated from the
source code files using the Doxygen documentation system, see
http://www.stack.nl/~dimitri/doxygen/.


Copyright and License
---------------------

Copyright 2007-2010 rafael grompone von gioi (grompone@gmail.com)

LSD is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

LSD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


Thanks
------

I would be grateful to receive any comment, especially about errors,
bugs, or strange results.

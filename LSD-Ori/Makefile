# -----------------------------------------------------------------------------
#
# LSD - Line Segment Detector on digital images
#
# Copyright 2007-2010 rafael grompone von gioi (grompone@gmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# -----------------------------------------------------------------------------

all: lsd lsd_call_example

lsd: lsd.c lsd.h lsd_cmd.c
	cc -O3 -lm -o lsd lsd_cmd.c lsd.c

lsd_call_example: lsd.c lsd.h lsd_call_example.c
	cc -lm -o lsd_call_example lsd_call_example.c lsd.c

strict: lsd.c lsd.h lsd_cmd.c
	cc -fPIC -ansi -Wall -Wextra -Werror -lm -o lsd_strict lsd.c lsd_cmd.c

doc: lsd.c lsd.h doxygen.config
	doxygen doxygen.config

clean:
	rm -f lsd lsd_call_example lsd_strict

cleandoc:
	rm -rf doc

# -----------------------------------------------------------------------------

#   bam2td: HOMER tag directory creator from BAM/SAM/CRAM
#   Copyright (C) 2026  Benjamin Jean-Marie Tremblay
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Build layout mirrors quaqc:
#   ./libs/htslib   full htslib source tree -> libs/htslib/libhts.a
#   ./libs/zlib     full zlib source tree   -> libs/zlib/libz.a
#   ./src/htslib    header-only vendored copy for -Isrc and the editor
#
# Default `make` produces a release build that links the static archives
# under libs/. Use hts_dyn=1 and/or z_dyn=1 to link the system libraries
# instead, e.g. `make release hts_dyn=1 z_dyn=1` to skip the libs/ build.

.PHONY: test debug release release-full clean mostlyclean libz libhts install uninstall

CC      ?=cc
CFLAGS  +=-std=gnu99 -Isrc
LDLIBS  +=-lm -lpthread
ZDIR    ?=./libs/zlib
HTSDIR  ?=./libs/htslib
PREFIX  ?=/usr/local
BINDIR  ?=bin
MANDIR  ?=share/man/man1
TESTDIR  =./tests

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer -DDEBUG
debug: bam2td

release: CFLAGS+=-DNDEBUG -O3 -flto
release: bam2td

release-full: libz libhts release

ZOPS=
ZLIB=

HTSOPS=
HTSLIB=

ifneq ($(native),)
	CFLAGS+=-march=native
	HTSOPS+=CFLAGS=-march=native
endif

ifeq ($(with_lzma),)
	HTSOPS+=--disable-lzma
else
	LDLIBS+=-llzma
endif

ifeq ($(with_bz2),)
	HTSOPS+=--disable-bz2
else
	LDLIBS+=-lbz2
endif

ifeq ($(with_curl),)
	HTSOPS+=--disable-libcurl
else
	LDLIBS+=-lcurl
ifneq ($(shell uname -s),Darwin)
	LDLIBS+=-lssl -lcrypto
endif
endif

ifeq ($(with_deflate),)
	HTSOPS+=--without-libdeflate
endif

ifeq ($(z_dyn),)
	HTSOPS+=LDFLAGS=-L$(ZDIR) CPPFLAGS=-I$(ZDIR)
	ZLIB=$(ZDIR)/libz.a
else
	LDLIBS+=-lz
endif

ifeq ($(hts_dyn),)
	HTSLIB=$(HTSDIR)/libhts.a
else
	LDLIBS+=-lhts
endif

libz/libz.a:
	(cd $(ZDIR) && ./configure --prefix=./ --static)
	$(MAKE) -C $(ZDIR)

libz: libz/libz.a

libhts/libhts.a:
	(cd $(HTSDIR) && ./configure $(HTSOPS))
	$(MAKE) -C $(HTSDIR) lib-static

libhts: libhts/libhts.a

clean/libz:
	$(MAKE) -C $(ZDIR) clean

clean/libhts:
	$(MAKE) -C $(HTSDIR) clean

clean/bam2td:
	-rm -f ./src/*.o
	-rm -f ./bam2td

clean: clean/libz clean/libhts clean/bam2td

mostlyclean/bam2td:
	-rm -f ./src/*.o

mostlyclean: clean/libz clean/libhts mostlyclean/bam2td

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $^ -o $@

objects := $(patsubst %.c,%.o,$(wildcard src/*.c))

bam2td: $(objects)
	$(CC) $(CFLAGS) $(LDFLAGS) $(objects) -o $@ $(HTSLIB) $(ZLIB) $(LDLIBS)

test: bam2td
	bash $(TESTDIR)/run.sh

install: bam2td
	install -p ./bam2td $(PREFIX)/$(BINDIR)

uninstall:
	-rm -f $(PREFIX)/$(BINDIR)/bam2td

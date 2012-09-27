CC = g++
OPTN = -O3
EXPORT = irbis-1.7.4
PROGRAMS=map_single best_match transf getsegm getwind getmuf trim irbis extend indexing correlation

################# ATTENTION THIS IS SET UP TO BE DP'S HOME DIRECTORY TO SAVE SPACE  ######################
################# RENAME THIS VARIABLE IF YOU HAVE TO USER CALLED DP IN HOME FOLDER ######################

#DBDIR=/home/dp/ # Please use this if you are a ma.fbb.msu.ru user
DBDIR=~/db/	# Use this path otherwise

CHAIN    = $(DBDIR)chain/
SEQ      = $(DBDIR)genome/
SOURCE   = $(DBDIR)download/
METADATA = ~/db/metadata/
OUTDIR   = ~/db/output/

.PHONY:	all install clean

all::	install

install:	${PROGRAMS} database download metacalc

export:
		mkdir $(EXPORT)/
		cp *.c $(EXPORT)/               
		cp *.h $(EXPORT)/               
		cp _* $(EXPORT)/
		cp *.pm $(EXPORT)/
		rm -f $(EXPORT)/setup.pm
		cp *.dat $(EXPORT)/
		cp -f *.mk $(EXPORT)/
		cp -f *.cfg $(EXPORT)/
		cp -f *.nwc $(EXPORT)/
		cp special $(EXPORT)/
		cp README $(EXPORT)/            
		cp VERSION $(EXPORT)/
		cp LICENCE $(EXPORT)/
		cp makefile $(EXPORT)/
		tar -cf $(EXPORT).tar $(EXPORT)/        
		gzip $(EXPORT).tar                      
		rm -f -r $(EXPORT)/
		mv $(EXPORT).tar.gz ..

setup.pm:	_setup
		./_setup setup.pm

database download: _makedatabase config.dat utils.pm setup.pm
		./_makedatabase

metacalc:	_makemetacalc config.dat utils.pm setup.pm
		./_makemetacalc

genutils.o:	genutils.c genutils.h
		$(CC) $(OPTN) -c genutils.c

map_single:	map_single.c genutils.o
		$(CC) $(OPTN) -o map_single map_single.c genutils.o

best_match:	best_match.c genutils.o
		$(CC) $(OPTN) -o best_match best_match.c genutils.o

getsegm:	getsegm.c genutils.o
		$(CC) $(OPTN) -o getsegm getsegm.c genutils.o

getwind:	getwind.c genutils.o
		$(CC) $(OPTN) -o getwind getwind.c genutils.o

transf:		transf.c genutils.o
		$(CC) $(OPTN) -o transf transf.c genutils.o

getmuf:		getmuf.c genutils.o
		$(CC) $(OPTN) -o getmuf getmuf.c genutils.o

trim:		trim.c dictionary.h genutils.h subset.h orderedset.h conservation.h genutils.o
		$(CC) $(OPTN) genutils.o trim.c -o trim

irbis: 	irbis.c dictionary.h genutils.h subset.h orderedset.h conservation.h genutils.o relation.h
		$(CC) $(OPTN) genutils.o irbis.c -o irbis

extend:		extend.c
		$(CC) $(OPTN) extend.c -o extend

indexing:	indexing.c
		$(CC) indexing.c -o indexing

correlation:	correlation.c tdistr.c
		g++ correlation.c genutils.o -o correlation

muscle:
		wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
		tar -xf muscle3.8.31_i86linux64.tar.gz
		rm muscle3.8.31_i86linux64.tar.gz
		mv muscle3.8.31_i86linux64 muscle

clean:
		rm -f genutils.o map_single transf getsegm getwind best_match getmuf indexing
		rm -f $(EXPORT).tar.gz database metacalc download setup.pm
		rm -f trim irbis extend correlation tmp*

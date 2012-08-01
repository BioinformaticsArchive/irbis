//	Copyright 2011,2012 Dmitri Pervouchine (dp@crg.eu)
//	This file is a part of the IRBIS package.
//	IRBIS package is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//	
//	IRBIS package is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with IRBIS package.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/dir.h>
#include "genutils.h"

int main(int argc, char* argv[]) {
    char outfilename[MAXBUFFLENGTH]="";
    char alnfilename[MAXBUFFLENGTH]="";
    char dbfilename[MAXBUFFLENGTH]="";
    char idxfilename[MAXBUFFLENGTH];
    char dbxfilename[MAXBUFFLENGTH];

    struct direct **files;
    char *ptr;
    long offset;
    long seqlen;
    int output_type=OUTPUT_TBMAF;
    long intronic_window = 150;
    long exonic_window = 0;

    char filename[MAXBUFFLENGTH];
    char name[MAXBUFFLENGTH];

    char chr_name[MAXBUFFLENGTH];
    char buff[MAXBUFFLENGTH];
    char longbuff[MAXLONGBUFFLENGTH];
    char longbuffm[MAXLONGBUFFLENGTH];

    FILE *idxfile;
    FILE *dbxfile;
    FILE *alnfile;
    FILE *outfile;
   
    int b,i,j,q,k,n,a,m,s;
    char c;
    char *pc;
    long x,y;

    long** pos;
    char** str;
    char** typ;
    char*** ids;
    long p,l;

    int record_count[MAXCHR];
    int record_idx[MAXCHR];
    int cis=0;
    int tolerant=0;

    char format[][64] = {"%*s %*i %*c %s %li %c %*s %s %c", "%s %li %c %*s %s %c"};

    if(argc==1) {
	fprintf(stderr,"Get sequence windows from custom compressed FASTA repository\n");
        fprintf(stderr,"Last update by (dp) on Oct 21, 2011\n");
        fprintf(stderr," -i ALN file name\n -d database file name (idx, dbx)\n -o output file name\n");
        fprintf(stderr," -we exonic window [0] -wi intronic window [150] -c cis [NO] -t treat S/E sites similarly to D/A sites\n");
	fprintf(stderr," -v suppress verbose output [NO] -fasta FASTA output -tab TAB delimited output -tabmaf TAB for MAF files [DEFAULT]\n");
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
	pc = argv[i];
	if(*pc != '-') continue;
        if(*(pc+1) == 'i') {
	    sscanf(argv[++i], "%s", &alnfilename[0]);
	}
	if(*(pc+1) == 'd') {
	    sscanf(argv[++i], "%s", &dbfilename[0]);
	}
        if(*(pc+1) == 'o') {
            sscanf(argv[++i], "%s", &outfilename[0]);
        }
        if(strcmp(pc+1,"we")==0) {
            sscanf(argv[++i], "%li", &exonic_window);
        }
        if(strcmp(pc+1,"wi")==0) {
            sscanf(argv[++i], "%li", &intronic_window);
        }
        if(*(pc+1) == 'v') {
	    verbose=0;
	}
        if(*(pc+1) == 'c') {
            cis=1;
        }
        if(*(pc+1) == 't') {
            tolerant=1;
        }
	if(strcmp(pc+1,"fasta")==0) {
	    output_type = OUTPUT_FASTA;
        }
        if(strcmp(pc+1,"tab")==0) {
            output_type = OUTPUT_TABSQ;
        }
        if(strcmp(pc+1,"tabmaf")==0) {
            output_type = OUTPUT_TBMAF;
        }
    }

    if(outfilename[0]==0) {
        fprintf(stderr,"No output file privided, exiting\n");
        exit(1);
    }
    outfile = fopen(outfilename,"w");
    if(outfile == NULL) {
        fprintf(stderr,"Can't open output file, exiting\n");
        exit(1);
    }

    alnfile = fopen(alnfilename,"r");
    if(alnfile == NULL) {
	fprintf(stderr,"Can't open BED input file, exiting\n");
        exit(1);
    }
    if(verbose) fprintf(stderr,"Reading ALN file pass 1");
    while(fgets(buff,MAXBUFFLENGTH,alnfile)) {
        if(strlen(buff)<2) break;
      	sscanf(buff,format[cis], &chr_name[0], &x, &c, &name[0], &c);
      	n = assign_code(chr_name);
      	record_count[n]++;
    }
    fclose(alnfile);

    pos = (long**) malloc(sizeof(long*)*(N_CHR_NAMES+1));
    str = (char**) malloc(sizeof(char*)*(N_CHR_NAMES+1));
    typ = (char**) malloc(sizeof(char*)*(N_CHR_NAMES+1));
    ids = (char***)malloc(sizeof(char**)*(N_CHR_NAMES+1));

    for(i=0;i<N_CHR_NAMES;i++) {
    	if(record_count[i]>0) {
	    pos[i] = (long*) malloc(sizeof(long)*(record_count[i]+1));
	    str[i] = (char*) malloc(sizeof(char)*(record_count[i]+1));
            typ[i] = (char*) malloc(sizeof(char)*(record_count[i]+1));
	    ids[i] = (char**) malloc(sizeof(char*)*(record_count[i]+1));
	    record_idx[i]=0;
	}
    }

    alnfile = fopen(alnfilename,"r");
    if(verbose) fprintf(stderr,"\nReading BED file pass 2");
    while(fgets(buff,MAXBUFFLENGTH,alnfile)) {
        if(strlen(buff)<2) break;
	sscanf(buff,format[cis],&chr_name[0], &x, &c, &name[0], &c);
	i = get_chr_code(chr_name);
	j = record_idx[i];
	sscanf(buff,format[cis],&chr_name[0], &pos[i][j], &str[i][j], &name[0], &typ[i][j]);
	ids[i][j] = (char*) malloc(sizeof(char)*(strlen(name)+1));
	strcpy(ids[i][j],name);
	record_idx[i]++;
    }
    if(verbose) fprintf(stderr,"\n");

    if(verbose) fprintf(stderr,"Reading database...");
    strcpy(idxfilename,dbfilename);
    strcat(idxfilename,".idx");
    strcpy(dbxfilename,dbfilename);
    strcat(dbxfilename,".dbx");
    idxfile = fopen(idxfilename,"r");
    dbxfile = fopen(dbxfilename,"r");
    if(idxfile == NULL || dbxfile == NULL) {
        if(verbose) fprintf(stderr,"Can't open database, exiting\n");
        exit(1);
    }

    offset = 0;
    while(fgets(buff,MAXBUFFLENGTH,idxfile)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s" , &name[0]);
	while(fgets(buff,MAXBUFFLENGTH,idxfile)) {
            if(strlen(buff)<2) break;
	    sscanf(buff,"%s %li" , &chr_name[0], &seqlen);
	    i = get_chr_code(chr_name);
	    for(k=0;k<record_count[i];k++) {
               if(pos[i][k]>seqlen) {
                    fprintf(stderr,"[WARNING: %li out of range %s:0-%li, ignored]\n", pos[i][k],chr_name,seqlen);
                    continue;
                }
		l = exonic_window + intronic_window;
		if(typ[i][k]=='D' || typ[i][k]=='A' || tolerant) {
		    if(str[i][k]=='+') {
		    	p = pos[i][k] - 1 - ((typ[i][k]=='D' || typ[i][k]=='E') ? exonic_window : intronic_window);
		    }
		    else {
			p = pos[i][k] - ((typ[i][k]=='D' ||  typ[i][k]=='E') ? intronic_window : exonic_window);
		    }
		    fget_segment(longbuff, dbxfile, offset, p, l);
		    if(str[i][k]=='-') {
		    	rev1(longbuff);
		    }
               	    if(is_all_n(longbuff)) continue;
		    if(output_type==OUTPUT_FASTA) 
			fprintf(outfile,">%s\n%s\n",ids[i][k], longbuff);
                    if(output_type==OUTPUT_TABSQ) 
			fprintf(outfile,"%s\t%s\n",ids[i][k], longbuff);
                    if(output_type==OUTPUT_TBMAF)
			fprintf(outfile,"%s\t%s\t%li\t%li\t%c\t%li\t%s\n",ids[i][k], chr_name, (str[i][k]=='+' ? p + 1 : seqlen - (p + l)), l, str[i][k], seqlen, longbuff);
		}
	    }
            offset+= (seqlen % 8 == 0) ? seqlen/8 : (seqlen/8 + 1);
	}
    }
    fclose(outfile);
    fclose(idxfile);
    fclose(dbxfile);
    if(verbose) fprintf(stderr,"\n");
    timestamp_report();
    exit(0);
}

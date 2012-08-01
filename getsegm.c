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
    char bedfilename[MAXBUFFLENGTH]="";
    char dbfilename[MAXBUFFLENGTH]="";
    char idxfilename[MAXBUFFLENGTH];
    char dbxfilename[MAXBUFFLENGTH];

    struct direct **files;
    char *ptr;
    long offset;
    long seqlen;
    int output_type=OUTPUT_TBMAF;
    long length_limit = (MAXLONGBUFFLENGTH-20);
    long half_length_limit = (MAXLONGBUFFLENGTH-20)/2;

    char filename[MAXBUFFLENGTH];
    char name[MAXBUFFLENGTH];

    char chr_name[MAXBUFFLENGTH];
    char buff[MAXBUFFLENGTH];
    char longbuff[MAXLONGBUFFLENGTH];

    FILE *idxfile;
    FILE *dbxfile;
    FILE *bedfile;
    FILE *outfile;
   
    int b,i,j,q,k,n,a,m,s;
    char c;
    char *pc;
    long x,y;

    long** beg;
    long** end;
    char** str;
    char*** ids;
    long pos,p,l;
    long margin=0;

    int record_count[MAXCHR];
    int record_idx[MAXCHR];
    int cis=0;

    if(argc==1) {
	fprintf(stderr,"Get sequences from custom compressed FASTA repository\n");
        fprintf(stderr,"Revised (dp) on Nov 8, 2011\n");
        fprintf(stderr," -i BED file name\n -d database file name (idx, dbx)\n -o output file name\n");
        fprintf(stderr," -v suppress verbose output [NO] -l [length limit] -m [margin]\n");
	fprintf(stderr," -fasta FASTA output -tab TAB delimited output -tabmaf TAB for MAF files [DEFAULT]\n");
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
	pc = argv[i];
	if(*pc != '-') continue;
        if(*(pc+1) == 'i') {
	    sscanf(argv[++i], "%s", &bedfilename[0]);
	}
	if(*(pc+1) == 'd') {
	    sscanf(argv[++i], "%s", &dbfilename[0]);
	}
        if(*(pc+1) == 'o') {
            sscanf(argv[++i], "%s", &outfilename[0]);
        }
        if(*(pc+1) == 'l') {
            sscanf(argv[++i], "%li", &length_limit);
        }
        if(*(pc+1) == 'v') {
	    verbose=0;
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
	if(*(pc+1) == 'm') {
            sscanf(argv[++i], "%li", &margin);
        }

    }
    half_length_limit = length_limit/2;

    if(outfilename[0]==0) {
        fprintf(stderr,"No output file privided, exiting\n");
        exit(1);
    }
    outfile = fopen(outfilename,"w");
    if(outfile == NULL) {
        fprintf(stderr,"Can't open output file, exiting\n");
        exit(1);
    }

    bedfile = fopen(bedfilename,"r");
    if(bedfile == NULL) {
	fprintf(stderr,"Can't open BED input file, exiting\n");
        exit(1);
    }
    if(verbose) fprintf(stderr,"Reading BED file pass 1");
    while(fgets(buff,MAXBUFFLENGTH,bedfile)) {
        if(strlen(buff)<2) break;
      	sscanf(buff,"%s %*i %*i %*i %*s %*c" , &chr_name[0]);
      	n = assign_code(chr_name);
      	record_count[n]++;
    }
    fclose(bedfile);

    beg = (long**) malloc(sizeof(long*)*(N_CHR_NAMES+1));
    end = (long**) malloc(sizeof(long*)*(N_CHR_NAMES+1));
    str = (char**) malloc(sizeof(char*)*(N_CHR_NAMES+1));
    ids = (char***)malloc(sizeof(char**)*(N_CHR_NAMES+1));

    for(i=0;i<N_CHR_NAMES;i++) {
    	if(record_count[i]>0) {
	    beg[i] = (long*) malloc(sizeof(long)*(record_count[i]+1));
	    end[i] = (long*) malloc(sizeof(long)*(record_count[i]+1));
	    str[i] = (char*) malloc(sizeof(char)*(record_count[i]+1));
	    ids[i] = (char**) malloc(sizeof(char*)*(record_count[i]+1));
	    record_idx[i]=0;
	}
    }

    bedfile = fopen(bedfilename,"r");
    if(verbose) fprintf(stderr,"\nReading BED file pass 2");
    while(fgets(buff,MAXBUFFLENGTH,bedfile)) {
        if(strlen(buff)<2) break;
	sscanf(buff,"%s %*i %*i %*i %*s %*c" , &chr_name[0]);
	i = get_chr_code(chr_name);
	j = record_idx[i];
	sscanf(buff,"%*s %li %li %*s %s %c" ,&beg[i][j], &end[i][j], &name[0], &str[i][j]);
        if(beg[i][j]>end[i][j]) {
            pos=beg[i][j];beg[i][j]=end[i][j];end[i][j]=pos;
        }
	ids[i][j] = (char*)malloc(sizeof(char)*(strlen(name)+1));
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
		if(end[i][k]>seqlen) {
		    fprintf(stderr,"[WARNING: %li out of range %s:0-%li, ignored]\n", end[i][k],chr_name,seqlen);
		    continue;
		}
		l = end[i][k]-beg[i][k];
		if(l==0) continue;
		beg[i][k]-=margin;
		l+=2*margin;
//		if(l>MAXLONGBUFFLENGTH-2) l = MAXLONGBUFFLENGTH-2;
		p = beg[i][k] - 1 + (str[i][k]=='-' ? 1 : 0);
		if(l<length_limit) {
		    fget_segment(longbuff, dbxfile, offset, p, l);
		}
		else {
		    fget_segment(longbuff, dbxfile, offset, p, half_length_limit);
		    strcat(longbuff + half_length_limit,".....");
                    fget_segment(longbuff+half_length_limit+5, dbxfile, offset, p+l-half_length_limit, half_length_limit);
		}
		if(str[i][k]=='-') {
		    rev1(longbuff);
		}
		if(is_all_n(longbuff)) continue;
		m = strlen(longbuff)-1;
		if(output_type==OUTPUT_FASTA) 
		    fprintf(outfile,">%s\n%s\n",ids[i][k], longbuff);
                if(output_type==OUTPUT_TABSQ) 
		    fprintf(outfile,"%s\t%s\n",ids[i][k], longbuff);
		if(output_type==OUTPUT_TBMAF) 
		    fprintf(outfile,"%s\t%s\t%li\t%li\t%c\t%li\t%s\n",ids[i][k], chr_name, (str[i][k]=='+' ? beg[i][k] : seqlen - end[i][k]), l, str[i][k], seqlen, longbuff);
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

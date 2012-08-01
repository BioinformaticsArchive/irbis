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

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/dir.h>
#include "genutils.h"

#define MAXKEYS 500000

int main(int argc, char** argv) {
    char outputfilename[MAXBUFFLENGTH]="";
    char chr[MAXBUFFLENGTH]; 
    char* filename[MAXSPECIES];

    char buff[MAXLONGBUFFLENGTH];
    char seq[MAXLONGBUFFLENGTH];

    int  key[MAXKEYS];

    FILE *inpfile;
    FILE *outputfile;

    char *pc;
    int	x, max_key;
    int	i, j, nk, n_species;

    long pos, len, total;
    char strand;
    int *ind;
    char **muf;

    if(argc==1) {
	fprintf(stderr, "MUF builder\n");
	fprintf(stderr, "Usage: %s -i [input_file_1] [input_file_2] ... [input_file_n] -o [output file]\n STDIN: cps file for subsetting\n",argv[0]);
	exit(1);
    }

    timestamp_set();
    n_species=0;
    for(i=1; i<argc; i++) {
        pc = argv[i];
        if(*pc != '-') continue;
        if(*(pc+1) == 'o') {
            sscanf(argv[++i], "%s", &outputfilename[0]);
        }
        if(*(pc+1) == 'i') {
	    for(;i+1<argc;i++) {
		pc = argv[i+1];
		if(*pc == '-') break;
		filename[n_species]    = (char*)malloc(sizeof(char)*MAXBUFFLENGTH);
		strcpy(filename[n_species++],pc);
	    }
        }
	if(*(pc+1) == 'v') {
	    verbose=0;
	}
    }

    outputfile = fopen(outputfilename,"w");
    if(outputfile == NULL) {
        fprintf(stderr,"Can't open output file, exiting\n");
        exit(1);
    }

    max_key=0;
    nk=0;
    while(fgets(buff,MAXBUFFLENGTH,stdin)) {
        sscanf(buff,"%*s %*i %*c %*i %i",&x);
        if(x>max_key) max_key=x;
        key[nk++]=x;
        if(nk>=MAXKEYS-1) {fprintf(stderr,"[WARNING: too many keys to hold, ignoring the rest]\n");break;}
    }

    ind = (int*)malloc(sizeof(int)*(max_key + 16));
    muf = (char**)malloc(sizeof(char*)*(max_key + 16));
    if(ind==NULL) {
        fprintf(stderr,"Can't create index (%i), exiting\n",max_key);
        exit(1);
    }

    for(i=0;i<=max_key;i++) ind[i]=0;
    for(j=0;j<nk;j++) {
	ind[key[j]]=1;
	muf[key[j]]=(char*)malloc(sizeof(char)*MAXLONGBUFFLENGTH*MAXSPECIES);
	if(muf[key[j]]==NULL) {
	    fprintf(stderr,"Can't allocate enough memory, exiting\n");
	    exit(1);
	}
	muf[key[j]][0]=0;
	if(verbose) fprintf(stderr,"[%i]",key[j]);
    }
    if(verbose) fprintf(stderr,"\n");

    for(i=0;i<n_species;i++) {
	inpfile = fopen(filename[i], "r");
	if(verbose) fprintf(stderr,"[%s, ", filename[i]);
	if(inpfile==NULL) {
	    fprintf(stderr,"Cannot open %s, exiting",filename[i]);
            exit(1);
        }
	for(j=strlen(filename[i])-1;filename[i][j]!='.' && j>=0;j--);
	filename[i][j]=0;
	for(j--;filename[i][j]!='/' && j>=0;j--);
        if(verbose) fprintf(stderr,"tag=%s", filename[i] + j + 1);
	while(fgets(buff, MAXLONGBUFFLENGTH, inpfile)) {
	    x=0;
	    sscanf(buff,"%i %s %li %li %c %li %s", &x, &chr[0], &pos, &len, &strand, &total, &seq[0]);
	    if(x==0) {fprintf(stderr,"!");continue; }
	    if(x<=max_key) {
		if(ind[x]) {
                    sprintf(buff,"s %s.%s %li %li %c %li %s\n", filename[i] + j + 1, chr, pos, len, strand, total, seq);	
		    strcat(muf[x],buff);
		}
 	    }
	}
	if(verbose) fprintf(stderr,"]\n");
    }

    for(i=0;i<=max_key;i++) {
	if(ind[i]) if(strlen(muf[i])>1) fprintf(outputfile,"a id=%i\n%s\n", i, muf[i]);	
    }

    fclose(outputfile);

    timestamp_report();
    exit(0);
}

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

#define MAXBUFFLENGTH 100000

int main(int argc, char* argv[]) {
    char buff[MAXBUFFLENGTH];
    char dirname[MAXBUFFLENGTH]="";
    char inpfilename[MAXBUFFLENGTH]="";
    char outfilename[MAXBUFFLENGTH]="";
    char extention[MAXBUFFLENGTH]="";
    char prefix[MAXBUFFLENGTH]="";
    char suffix[MAXBUFFLENGTH]="";
    char singlefilename[MAXBUFFLENGTH]="";
    long pos;

    FILE *inpfile, *outfile;
    struct direct **files;
    char *pc;
    int i, n, key;
    char *ptr;

    if(argc==1) {
	fprintf(stderr,"This script index for single sequence (sus/suw) files\n");
        fprintf(stderr,"Usage: %s -dir dirname -file filename -extention extention -suffix suffix -prefix prefix\n",argv[0]);
	fprintf(stderr,"If filename unspecified then every file in dirname/ with extention extention is processed\n");
	fprintf(stderr,"The output is of the form dirname/<prefix><filename><suffix>\n");
        exit(1);
    }

    for(i=1;i<argc;i++) {
	pc = argv[i];
        if(*pc != '-') continue;
        if(strcmp(pc+1,"file")==0)	sscanf(argv[++i], "%s",  &singlefilename[0]);
	if(strcmp(pc+1,"dir")==0)	sscanf(argv[++i], "%s",  &dirname[0]);
	if(strcmp(pc+1,"extension")==0)	sscanf(argv[++i], "%s",  &extention[0]);
	if(strcmp(pc+1,"prefix")==0) 	sscanf(argv[++i], "%s",  &prefix[0]);
	if(strcmp(pc+1,"suffix")==0)    sscanf(argv[++i], "%s",  &suffix[0]);
    }

    if(singlefilename[0]) {
	strcpy(inpfilename, dirname);
	strcat(inpfilename, singlefilename);
	inpfile = fopen(inpfilename, "r");
	if(inpfile == NULL) exit(1);
	strcpy(outfilename, dirname);
	strcat(outfilename, prefix);
	strcat(outfilename, singlefilename);
	strcat(outfilename, suffix);
        outfile = fopen(outfilename, "w");
        fprintf(stderr, "[%s", outfilename);
        if(outfile == NULL) exit(1);
	pos=0;
        while(fgets(buff,MAXBUFFLENGTH,inpfile)) {
            sscanf(buff,"%i",&key);
            fprintf(outfile,"%i\t%li\n", key, pos);
            pos = ftell(inpfile);
        }
        fclose(inpfile);
        fclose(outfile);
        fprintf(stderr, "]\n");
	exit(0);
    }

    fprintf(stderr, "Scanning %s for *.%s\n",dirname,extention); 
    n=scandir(dirname,&files,0,alphasort);
    for(i=0;i<n;i++) {
	ptr = rindex(files[i]->d_name,'.');
        if(ptr != NULL && strcmp(ptr, extention)==0) {
	    strcpy(inpfilename, dirname);
	    strcat(inpfilename, files[i]->d_name);
	    inpfile = fopen(inpfilename, "r");
	    if(inpfile == NULL) continue;
	    strcpy(outfilename, dirname);
	    strcat(outfilename, prefix);
	    strcat(outfilename, files[i]->d_name);
	    strcat(outfilename, suffix);
	    outfile = fopen(outfilename, "w");
	    if(outfile == NULL) continue;
	    fprintf(stderr, "[%s",outfilename);
	    pos=0;
	    while(fgets(buff,MAXBUFFLENGTH,inpfile)) {
 		sscanf(buff,"%i",&key);
		fprintf(outfile,"%i\t%li\n", key, pos);
		pos = ftell(inpfile);
	    }
	    fclose(inpfile);
	    fclose(outfile);
	    fprintf(stderr, "]\n");
	}
    }
    exit(0);
}

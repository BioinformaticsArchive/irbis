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

#include "genutils.h"

#define MAXALN 2000000

FILE *cpsfile;
FILE *outfile;
FILE *chainfile;

int main(int argc, char* argv[]) {
    char cpsfilename[MAXBUFFLENGTH];
    char chainfilename[MAXBUFFLENGTH];
    char outfilename[MAXBUFFLENGTH]="";
 
    int marginlength = 0;

    int MAXREC; 

    char c;
    char *pc;
    int x;

    char buff[MAXBUFFLENGTH+1];
    char aux[MAXBUFFLENGTH+1];

    long score;
    int start1,end1,len1,start2,end2,len2;
    char strand1, strand2, chr1[MAXBUFFLENGTH], chr2[MAXBUFFLENGTH];

    int *size, *dq, *dt;
    int a,b,k,i,j,s,m;

    int *position;
    char *strand;
    int *ids;
    int *idg;
    char *type;

    int chridx[MAXCHR+1];
    int chroff[MAXCHR+1];

    char resstr;
    int  rescrd;


    if(argc==1) {
	fprintf(stderr,"Finds matches of the given set of sites (CPS file) in the BLASTZ chain alignment (CHAIN file)\n");
        fprintf(stderr,"Last update by (dp) on Sep 21, 2011\n");
	fprintf(stderr,"Keys:\n -i CPS file (remember to sort by position in ascending order)\n -d CHAIN alignment file\n -o output file\n");
 	fprintf(stderr," -m margin length [0]\n -v suppress verbose output [NO]\n");
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
	pc = argv[i];
	if(*pc != '-') continue;
        if(*(pc+1) == 'i') {
	   sscanf(argv[++i], "%s", &cpsfilename[0]);
	}
	if(*(pc+1) == 'd') {
	   sscanf(argv[++i], "%s", &chainfilename[0]);
	}
        if(*(pc+1) == 'o') {
           sscanf(argv[++i], "%s", &outfilename[0]);
        }
        if(*(pc+1) == 'm') {
           sscanf(argv[++i], "%i", &marginlength);
        }
        if(*(pc+1) == 'v') {
	   verbose=0;
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

    for(i=0;i<MAXCHR;i++) {
      	chridx[i]=chroff[i]=0;
    }

    MAXREC = 0;
    cpsfile= fopen(cpsfilename,"r");
    if(cpsfile==NULL) {
	fprintf(stderr,"Can't access CPS file. Exiting\n");
	exit(1);
    }

    if(verbose) fprintf(stderr,"Reading CPS input pass 1");

    while(!feof(cpsfile)) {
      	buff[0]=0;
      	fgets(buff,MAXBUFFLENGTH,cpsfile);
      	if(strlen(buff)<2) break;
      	sscanf(buff,"%s" , aux);
      	chridx[assign_code(aux)]++;
	MAXREC++;
    }
    fclose(cpsfile);

    for(s=i=0;i<MAXCHR;i++) {
        x = chridx[i];
	chridx[i] =s;
	s+=x;
    }
    chridx[i] = s;

    position   = (int*)  malloc(sizeof(int)*(s+4));
    strand     = (char*) malloc(sizeof(char)*(s+4));
    type       = (char*) malloc(sizeof(char)*(s+4));
    ids	       = (int*)  malloc(sizeof(int)*(s+4));
    idg        = (int*)  malloc(sizeof(int)*(s+4));

    if(position==NULL || strand==NULL || type==NULL || ids==NULL || idg==NULL) {
        fprintf(stderr,"Not enough memory. Terminated\n");
        exit(1);
    }

    cpsfile= fopen(cpsfilename,"r");
    if(verbose) fprintf(stderr,", records = %i\nReading CPS input pass 2",MAXREC);
    while(!feof(cpsfile)) {
      	buff[0]=0;
      	fgets(buff,MAXBUFFLENGTH,cpsfile);
        if(strlen(buff)<2) break;
        sscanf(buff,"%s" , aux);
        i = assign_code(aux);
        m = chridx[i]+chroff[i];
        sscanf(buff,"%*s %i %c %i %i %c" , position+m, strand+m, idg+m, ids+m, type+m);
        chroff[i]++;
    }
    fclose(cpsfile);

    if(verbose) fprintf(stderr,"\nSorting segments");

/*
    for(i=0;i<MAXCHR;i++) {
	quickSort_ic(position,strand,chridx[i],chridx[i+1]-1);
    }
*/
	

    for(i=0;i<MAXCHR;i++) {
        k=1;
        while(k) {
            k=0;
            for(j=chridx[i];j<chridx[i+1]-1;j++) {
                if(position[j]>position[j+1]) {
                    k=1;
                    swapi(position+j,position+j+1);
                    swapc(strand+j,strand+j+1);
                    swapc(type+j,type+j+1);
		    swapi(ids+j,ids+j+1);
                    swapi(idg+j,idg+j+1);
                }
            }
        }
    }

    if(verbose) fprintf(stderr," done\nProcessing chains");

/**********************************************************************************************/
    size = (int*) malloc(sizeof(int)*MAXALN);
    dq   = (int*) malloc(sizeof(int)*MAXALN);
    dt   = (int*) malloc(sizeof(int)*MAXALN);

    if(size ==0 || dq ==0 || dt==0) {
        fprintf(stderr,"Not enough memory for such long chains. Terminated\n");
        exit(1);
    }


/**********************************************************************************************/

    chainfile = fopen(chainfilename,"r");
    while(!feof(chainfile)) {
     	buff[0]=0;
     	fgets(buff,MAXBUFFLENGTH,chainfile);
     	if(strlen(buff)<2) break;
     	buff[5]=0;
     	if(strcmp(buff,"chain")==0) {
       	    sscanf(buff+6,"%li %s %i %c %i %i %s %i %c %i %i",&score, &chr1[0], &len1, &strand1, &start1, &end1, &chr2[0], &len2, &strand2, &start2, &end2);
	    k=0;
	    while(!feof(chainfile)) {
	    	buff[0]=0;
	    	fgets(buff,MAXBUFFLENGTH,chainfile);
            	if(strlen(buff)<2) break;
	    	sscanf(buff,"%i %i %i",&size[k],&dt[k],&dq[k]);
	    	k++;
	    	if(k>MAXALN) {
		    fprintf(stderr,"Chain length exceeded. Terminating");
		    exit(1);
	    	}
	    }

	    x = get_chr_code(chr1);
	    if(x<0) continue;

            a=start1;b=start2;
            j=0;

	    for(i=chridx[x];i<chridx[x+1] && position[i]<start1;i++);
	    for(;i<chridx[x+1]&& position[i]<end1;i++) {
	    	while(position[i]>a+size[j]+dt[j] && j<k){
		    a+=size[j]+dt[j];
		    b+=size[j]+dq[j];
		    j++;
	    	}
	        if(j>=k) break;
	        if(position[i]-a > marginlength && a+size[j]-position[i] >= marginlength) {
                    if(strand1==strand2) {
                    	resstr = strand[i];
                    	rescrd = position[i] - a + b;
                    }
                    else {
                    	resstr = (strand[i]=='+') ? '-' : '+';
                    	rescrd = len2 - (position[i] - a + b - 1) ;
                    }
		    fprintf(outfile,"%s\t%i\t%c\t%s\t%i\t%c\t%i\t%i\t%c\t%li\n",chr1, position[i], strand[i],chr2,rescrd, resstr, idg[i], ids[i],type[i],score);
	    	}
	    }
     	}
    }
    if(verbose) fprintf(stderr," done\n");
    fclose(chainfile);
    fclose(outfile);
    timestamp_report();
    exit(0);
}

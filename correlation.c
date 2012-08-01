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
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "tdistr.c"

#define MAXBUFFLENGTH 10000
#define NAMELENGTH 64
#define MAXCOL 100
#define MAXROW 300000

typedef struct {
    int ncols;
    int nrows;
    char** header;
    char** rowname;
    double** data;
    int** flag;
} data;

int spearman = 0;
int mindf = 3;
int symmetric = 0;
double  alpha = 0.05;
int texoutput = 0;
char texline[MAXBUFFLENGTH];

void compute_corr(double *x, int *fx, double *y, int *fy, char **h, int n, double *res, double *pv, int *df);

data* read_file(FILE *f);

int main(int argc, char* argv[]) {
    char leftfilename[MAXBUFFLENGTH];
    char rightfilename[MAXBUFFLENGTH];
    char outputfilename[MAXBUFFLENGTH];

    FILE *leftfile;
    FILE *rightfile;
    FILE *outputfile;

    int i,j,n,m,k;
    char *pc;

    data* left;
    data* right;

    FILE*   t_dist;
    double r,p;
    int df;
    int flag = 0;

    outputfilename[0]=0;

    if(argc==1) {
	fprintf(stderr, "Compute correlation between two tables\n -l left partner or stdin\n -r right partner or stdin\n -a significance level [0.05]\n");
	fprintf(stderr, " -s use Spearman correlation [no]\n -d min degrees of fredom [3]\n -u Upper triangular (ignore indices row>column)\n");
	exit(1);
    }

    for(i=1;i<argc;i++) {
        pc = argv[i];
        if(*pc == '-') {
            if(strcmp(pc+1,"l") == 0) sscanf(argv[++i], "%s",  &leftfilename[0]);
	    if(strcmp(pc+1,"r") == 0) sscanf(argv[++i], "%s",  &rightfilename[0]);
	    if(strcmp(pc+1,"o") == 0) sscanf(argv[++i], "%s",  &outputfilename[0]);
	    if(strcmp(pc+1,"a") == 0) sscanf(argv[++i], "%lf", &alpha);
	    if(strcmp(pc+1,"d") == 0) sscanf(argv[++i], "%i",  &mindf);
	    if(strcmp(pc+1,"s") == 0) spearman  = 1;
	    if(strcmp(pc+1,"u") == 0) symmetric = 1;
	    if(strcmp(pc+1,"t") == 0) texoutput = 1;
        }
    }

    leftfile = fopen(leftfilename, "r");
    if(leftfile==NULL) {
	fprintf(stderr, "[WARNING: left input set to STDIN]\n");
	leftfile=stdin;
    }

    rightfile = fopen(rightfilename,"r");
    if(rightfile==NULL) {
	fprintf(stderr, "[WARNING: right input set to STDIN]\n");
	rightfile=stdin;
    }

    if(leftfile==stdin && rightfile==stdin) {
	fprintf(stderr, "Both inputs cannot be STDIN, exiting\n");
	exit(1);
    }

    outputfile = (outputfilename[0]==0) ? NULL : fopen(outputfilename, "w");
    if(outputfile==NULL) {
	fprintf(stderr, "[WARNING: output set to STDOUT]\n");
	outputfile=stdout;
    }

    left  = read_file(leftfile);
    fprintf(stderr, "[left=%i]\n", left->nrows);

    right = read_file(rightfile);
    fprintf(stderr, "[right=%i]\n", right->nrows);

    if(left->ncols!=right->ncols) {
	fprintf(stderr, "Number of columns are not equal, exiting\n");
        exit(1);
    }

    flag = 0;
    for(i=0;i<left->ncols;i++) {
	if(strcmp(left->header[i],right->header[i])!=0) flag=1;
    }
    fprintf(stderr, "[dim=%i]\n",left->ncols);

    if(flag) {
	fprintf(stderr, "Column names are not the same, exiting\n");
        exit(1);
    }

    k=0;
    for(i=0;i<left->nrows;i++) {
        for(j=0;j<right->nrows;j++) {
            if(symmetric && i>j) continue;
            compute_corr(left->data[i], left->flag[i], right->data[j], right->flag[j], left->header, left->ncols, &r, &p, &df);
            if(p<alpha && df>=mindf) {
                fprintf(outputfile, "%s\t%s\t%lf\t%lf\t%i\t%s\n",left->rowname[i], right->rowname[j], r, p, df,texline);
                k++;
            }
        }
    }
    fprintf(stderr, "[output=%i]\n", k);

}

void replace_by_ranks(double *x, int n) {
    int i, j, k;
    int *r;
    double *s;
    double y;
    int flag=1;

    r = (int*) malloc(sizeof(int)*(n+1));
    s = (double*) malloc(sizeof(double)*(n+1));

    for(i=0;i<n;i++) r[i]=i;
    while(flag) {
        flag=0;
        for(i=0;i<n-1;i++) {
            if(x[i]>x[i+1]) {
                y = x[i]; x[i] = x[i+1]; x[i+1] = y;
                j = r[i]; r[i] = r[i+1]; r[i+1] = j;
                flag = 1;
            }
        }
    }
    for(i=0;i<n;i=j) {
        for(j=i; x[i]==x[j] && j<n; j++);
        for(k=i; k<j; k++) s[r[k]]=(double)(i + j - 1)/2;

    }
    for(i=0;i<n;i++) {
        x[i]=s[i];
    }

    free(r);
    free(s);
}


data* read_file(FILE *f) {
    char buff[MAXBUFFLENGTH];
    int i, j, n, m;
    data *ptr = (data*)malloc(sizeof(data));
    ptr->ncols = ptr->nrows = 0;
    ptr->header = (char**)malloc(sizeof(char*)*MAXCOL);

    fgets(buff, MAXBUFFLENGTH, f);
    n = strlen(buff);
    for(i=j=0;i<n;i++) {
        if(buff[i]<=32 && buff[i+1]>32) j = i + 1;
        if(buff[i]>32 && buff[i+1]<=32) {
            buff[i+1]=0;
            ptr->header[ptr->ncols] = (char*)malloc(sizeof(char)*(i-j+2));
            strcpy(ptr->header[ptr->ncols++], &buff[j]);
        }
    }

    ptr->rowname = (char**)malloc(sizeof(char*)*(MAXROW + 1));
    ptr->data    = (double**)malloc(sizeof(double*)*(MAXROW + 1));
    ptr->flag    = (int**)malloc(sizeof(int*)*(MAXROW + 1));

    if(ptr->rowname==NULL || ptr->data==NULL || ptr->flag==NULL) {
        fprintf(stderr, "Not enough memory, exiting\n");
        exit(1);
    }

    while(fgets(buff, MAXBUFFLENGTH, f)) {
        if(strlen(buff)<2) break;
        ptr->rowname[ptr->nrows] =   (char*) malloc(sizeof(char)*NAMELENGTH);
        ptr->data[ptr->nrows]    = (double*) malloc(sizeof(double)*(ptr->ncols+1));
        ptr->flag[ptr->nrows]    =    (int*) malloc(sizeof(int)*(ptr->ncols+1));
        n = strlen(buff);
        m = 0;
        for(i=0;i<n;i=j) {
	    for(;i<n && buff[i]<=32;i++);
            for(j=i;j<n && buff[j]>32;j++);
            buff[j]=0;
            if(m==0) {
                sscanf(&buff[i], "%s",  &ptr->rowname[ptr->nrows][0]);
            }
            else {
                if(strcmp(&buff[i],(char*)"NA")==0) {
                    ptr->data[ptr->nrows][m-1] = INFTY;
                    ptr->flag[ptr->nrows][m-1] = 0;
                }
                else {
                    sscanf(&buff[i], "%lf", &ptr->data[ptr->nrows][m-1]);
                    ptr->flag[ptr->nrows][m-1] = 1;
                }
            }
            m++;
            j++;
        }
        if(spearman) replace_by_ranks(ptr->data[ptr->nrows], m);
        ptr->nrows++;
	if(ptr->nrows>=MAXROW) {
	    fprintf(stderr, "More than %i lines; skipping the rest",MAXROW);
	    break;
	}
    }
/*
    for(i=0;i<ptr->ncols;i++) {
	fprintf(stderr,"\t%s",ptr->header[i]);
    }
    fprintf(stderr,"\n");

    for(j=0;j<ptr->nrows;j++) {
	fprintf(stderr,"%s",ptr->rowname[j]);
	for(i=0;i<ptr->ncols;i++) {
	    if(ptr->flag[j][i]) fprintf(stderr,"\t%2.2lf", ptr->data[j][i]); else fprintf(stderr,"\tNA");
	}
	fprintf(stderr,"\n");
    }
*/
    return(ptr);
}

void compute_corr(double *x, int *fx, double *y, int *fy, char **h, int n, double *res, double *pv, int *df) {
    char buff[MAXBUFFLENGTH];
    int i, m;
    double sx, sy, sx2, sy2, sxy;
    double r, t;

    sx = sy = sx2 = sy2 = sxy = 0;
    m = 0;
    texline[0]=0;

    for(i=0;i<n;i++) {
        if(fx[i]>0 && fy[i]>0) {
            sx+=x[i];
            sy+=y[i];
            sx2+=x[i]*x[i];
            sy2+=y[i]*y[i];
            sxy+=x[i]*y[i];
            m++;
	    if(texoutput) {
		sprintf(&buff[0], "(%1.3lf,%1.3lf)[%s]", x[i], y[i], h[i]);
		strcat(texline,buff);
	    }
        }
    }
    *res = 0;
    *pv  = 1;
    *df  = 0;
    if(m<2) return;

    *df  = m-2;

    sx2 = sx2 - sx*sx/m;
    sy2 = sy2 - sy*sy/m;
    sxy = sxy - sx*sy/m;

    if(sx2<=0 || sy2<=0) return;

    r = sxy/sqrt(sx2*sy2);
    t = abs(r)<1 ? r*sqrt(m-2)/sqrt(1-r*r) : INFTY;

    *res = r;
    *pv  = eval_t(t,m-2);
}

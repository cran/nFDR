/*////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////*/
/*                                                    */
/*                 C Program for                      */
/*            Nonparametric FDR Estimate              */
/*         Using the Bernstein Polynomials            */
/*                                                    */
/*////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////////////////////////*/
/*                                                                                */
/*  Reference:                                                                    */
/*    Zhong Guan, Baolin Wu and Hongyu Zhao,                                      */
/*          Nonparametric estimator of false discovery rate based on              */
/*              Bernstein polynomials                                             */
/*                                                                                */
/*////////////////////////////////////////////////////////////////////////////////*/
/*#include <conio.h> */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define RNMX (1.0-EPS)
#define MAXIT 10000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define DBL_MAX         1.7976931348623157e+308
#define DBL_MAX_EXP     1024
/*long RX=1,RY=1,RZ=1; */
/*#define M_LN2     0.693147180559945309417232121458 */   /* ln(2) */
#define PI 3.141592653589793238462643383279502884197169399375
#define R_D_Lval(p) (lower_tail ? (p) : (1 - (p)))  /*  p  */

#define R_D_val(x)  (log_p  ? log(x) : (x))     /*  x  in pF(x,..) */

#define R_DT_val(x) R_D_val(R_D_Lval(x))        /*  x  in pF */
/*#define min(i, j) (((i) < (j)) ? (i) : (j))  */
/*############################################################  */
/*static double dsqrarg;                                        */
/*#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)*/

/*static double dmaxarg1,dmaxarg2;                               */
/*#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\ */
/*        (dmaxarg1) : (dmaxarg2))          */


double fmax2(double x, double y){
    return (x < y) ? y : x;
}

double fmin2(double x, double y){
    return (x < y) ? x : y;
}
int imin2(int x, int y){
    return (x < y) ? x : y;
}
int imax2(int x, int y){
    return (x < y) ? y : x;
}
int fcompare (const void * a, const void * b){
    if(*(double*)a - *(double*)b>0) return 1;
    else if(*(double*)a - *(double*)b<0) return -1;
    else return 0;
}
void ErrorMessage(char error_text[]){
        fprintf(stderr,"Run-Time Error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}
char *cvector(int nl,int nh){
        char *v;

        v=(char *)malloc((unsigned) (nh-nl+1)*sizeof(char));
        if (!v) ErrorMessage("allocation failure in cvector()");
        return v-nl;
}

void free_cvector(char *v, int nl, int nh){
        free((char*) (v+nl));
}

int *ivector(int nl,int nh){
        int *v;

        v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
        if (!v) ErrorMessage("allocation failure in ivector()");
        return v-nl;
}

void free_ivector(int *v, int nl, int nh){
        free((char*) (v+nl));
}

int **imatrix(int nrl,int nrh,int ncl,int nch){
        int i, **m;

        m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
        if (!m) ErrorMessage("allocation failure 1 in imatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++){
            m[i]=ivector(ncl, nch);
         }
        return m;
}
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch){
        int i;

        for(i=nrh;i>=nrl;i--) free_ivector(m[i], ncl, nch);
}
int ***iarray(int ml, int mh, int nl, int nh, int pl, int ph){
    int i, ***a;

        a=(int ***) malloc((unsigned) (mh-ml+1)*sizeof(int **));
        if (!a) ErrorMessage("allocation failure in iarray()");
        a -= ml;

        for(i=ml;i<=mh;i++) {
            a[i]=imatrix(nl, nh, pl, ph);
        }
        return a;
}
void free_iarray(int ***a, int ml, int mh, int nl, int nh, int pl, int ph){
    int i;

    for(i=mh;i>=ml;i--) free_imatrix(a[i], nl, nh, pl, ph);
}
double *dvector(int nl,int nh){
        double *v;

        v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
        if (!v) ErrorMessage("allocation failure in dvector()");
        return v-nl;
}

void free_dvector(double *v, int nl, int nh){
        free((char*) (v+nl));
}
double **dmatrix(int nrl,int nrh,int ncl,int nch){
        int i;
        double **m;

        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m) ErrorMessage("allocation failure 1 in dmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i]) ErrorMessage("allocation failure 2 in dmatrix()");
                m[i] -= ncl;
        }
        return m;
}
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch){
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}

double ***darray(int ml, int mh, int nl, int nh, int pl, int ph){
    int i;
    double ***a;

        a=(double ***) malloc((unsigned) (mh-ml+1)*sizeof(double**));
        if (!a) ErrorMessage("allocation failure in darray()");
        a -= ml;

        for(i=ml;i<=mh;i++) {
            a[i]=dmatrix(nl, nh, pl, ph);
        }
        return a;
}
void free_darray(double ***a, int ml, int mh, int nl, int nh, int pl, int ph){
    int i;
    for(i=mh;i>=ml;i--) free_dmatrix(a[i], nl, nh, pl, ph);
}
double Mean(double *x, int s){
    int i;
    double sum;
    sum = 0.0;
    for (i=1;i<=s;i++) sum += x[i];
    sum = sum / s;
    return sum;
}
double Median(double *x, int s)
{
    int i;
    double *vec;
    double mdn;
    vec = dvector(0,s);

    for(i=0; i<=s-1; i++) *(vec + i) = x[i+1];
    qsort(vec, s, sizeof(double), fcompare);
/*    quicksort(vec, s);   */
    if ( s % 2 == 0) mdn = 0.5*(vec[ (int) (s/2) -1] + vec[(int) (s/2)]);
    else mdn = vec[(int) ((s-1)/2)];

    free_dvector(vec, 0,s);
    return mdn;

}

double StDev(double *x, int s)
{
    int i;
    double mean, sum;
    if(s<=1) ErrorMessage("Sample size must be greater that 1 when calculating Stdev!");
    mean = 0;
    sum=0;
    for(i=1;i<=s;i++) mean+=x[i];
    mean = mean / (double) s;
    for(i=1;i<=s;i++) sum += (x[i]-mean)*(x[i]-mean);
    sum = sqrt(sum / (double) (s-1));
    return sum;
}
double betai(double a, double b, double x)
{
    double betacf(double a, double b, double x);
    double gammln(double xx);
    void ErrorMessage(char error_text[]);
    double bt;
    if (x < 0.0 || x > 1.0) ErrorMessage("Bad x in Routine betai");
    if (x == 0.0 || x == 1.0) bt=0.0;
    else bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
    if (x < (a+1.0)/(a+b+2.0))
        return bt*betacf(a,b,x)/a;
    else
        return 1.0-bt*betacf(b,a,1.0-x)/b;
}
double betacf(double a, double b, double x)
{
    void ErrorMessage(char error_text[]);
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++){
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    if (m > MAXIT) ErrorMessage("a or b too big, or MAXIT too small in betacf");
    return h;
}

double gammln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}


/* Function order Returns an integer vector containing  */
/* the permutation that will sort the input into ascending order. */

void forder(double *vec, int leng, int *ord)
{
    int i, j;
    double *temp,*z;
    z=dvector(0,leng-1);
    temp=dvector(0,leng-1);
    for(i=0;i<leng;i++) {
        z[i]=vec[i];
        temp[i]=z[i];
    }
    qsort(z, leng, sizeof(double), fcompare);
    for(i=0;i<leng;i++) {
        for(j=leng-1;j>=0;j--){
            if(z[i]==temp[j]){
                ord[i]=j+1;
            }
        }
        temp[ord[i]-1]=z[0]-1000.0;
    }
    free_dvector(z,0,leng-1);
    free_dvector(temp,0,leng-1);
}


int compare (const void * a, const void * b){
    return ( *(int*)a - *(int*)b );
}


void iorder(int *vec, int leng, int *ord){
    int i, j,*temp,*z;
        z=ivector(0,leng-1);
        temp=ivector(0,leng-1);
    for(i=0;i<leng;i++){
        z[i]=vec[i];
        temp[i]=z[i];
    }
    qsort(z, leng, sizeof(int), compare);
    for(i=0;i<leng;i++){
        for(j=leng-1;j>=0;j--){
            if(z[i]==temp[j]){
                ord[i]=j+1;
            }
        }
        temp[ord[i]-1]=z[0]-100;
    }
    free_ivector(z,0,leng-1);
    free_ivector(temp,0,leng-1);
}
double Rank(double y, double *x, int leng){
    int i,k;
    double rnk=0.0,*z;
    z=dvector(0,leng);
    for(i=0;i<leng;i++)
        z[i]=x[i];
    qsort(z,leng, sizeof(double), fcompare);
    k=0;
    for(i=0;i<leng;i++){
        if(z[i]==y){
            k+=1;
            rnk=i;
        }
    }
    free_dvector(z,0,leng);
    rnk=rnk-(double)((k-3.0)/2.0);
    return(rnk);
}
int choose(int n, int m, int order)
{
    int nCm = 0;
    if(n<=0 || m<0) ErrorMessage("n<=0 or m<0!");
    if(order != 0) nCm = (int) (exp(gammln((double) (n + 1)) - gammln((double) (n - m + 1))));
    else nCm = (int) (exp(gammln((double) (n + 1)) - gammln((double) (n - m + 1))-gammln((double) (m + 1))));
    return(nCm);
}

void ecdf(double *x, double *y, int s, int t, double *edf)
/* x = original data                                 */
/* y = vector of points at which ecdf are evaluated  */
/* s = sample size, length of x                      */
/* t = length of y                                   */
{
    int i, j;
/*  qsort(x, s, sizeof(double), fcompare);  */
    for(j=0;j<t;j++){
        edf[j] = 0.;
        for(i=0;i<s;i++){
            edf[j] += (x[i]<=y[j]);
        }
        edf[j]/=(double)s;
    }
}
double Bernstein_polynom(int j, int k, double t)
/*  Calculate Bernstein ploynomials b(j; k, t)*/
/* bp = Bernstein polynomial at t             */
{
    double bp;
    if(j==0) bp = pow(1.0-t,k);
    else if(j==k) bp = pow(t,k);
    else bp = betai(j, k-j+1.0, t)-betai(j+1.0, k-j, t);
/*  if(t<0.0 || t>1.0) printf("Bernstein_polynom: t = %f\n", t); */
    return bp;
}
void Bf(double *f, int s, double *x, int t, double *bf)
/*  f = sample value of f                            */
/*  s = sample size, length of f                     */
/*  x = vector of points at which Bf is evaluated    */
/*  t = length of x                                  */
/* bf = Bernstein polynomial approx of f at x        */
{
    int i, j;
    for(j = 0; j<t; j++){
        bf[j] = f[0]*pow(1.-x[j],s-1);
        for(i = 1; i<s-1; i++){
            bf[j] += f[i]*(betai(i, s-i, x[j])-betai(i+1., s-i-1., x[j]));
            if(x[j]<0.0 || x[j]>1.0) printf("Bf: x[j] = %f\n", x[j]);
        }
        bf[j] += f[s-1]*pow(x[j],s-1);
    }
}
double hk(int r, int k){
/* r = integers  */
/* k = int       */
    int i, j;
    double hr, mtemp = 0.0;
    if(r>k) ErrorMessage("r is greater than k in hk");
    if(fmin2(r,k)<=0) ErrorMessage("Zero r or k value in hk");
    hr = 0.0;
    for(j = 0; j<k; j++){
        mtemp = 0.0;
        for(i=0; i<r; i++) mtemp+=Bernstein_polynom(j, k-1, 1.- (double) (i+1.)/(double)k);
        mtemp /=(double) r;
        hr += mtemp*mtemp;
    }
    return hr;
}



double mse_pi0_approx(int m, int r, int k, double *ftilde, double *fp1, double *fp2){
    int i, j;
    double R12, PI0, mse, bbar;
    double temp0 = 0., temp2 = 0., temp3 = 0., bias = 0.;
    PI0 = 0.;
    for(i=k-r; i<k; i++) PI0 += ftilde[i];
    PI0/= (double) r;
    R12 = .0;
/*printf("PI0=%f\n", PI0);  */
    for(j=0; j<k; j++){
        bbar = 0.0;
        for(i=1; i<=r; i++) bbar+=Bernstein_polynom(j, k-1, 1.0-(double)i/(double)k);
        bbar/=(double)r;
        temp0 = 0.;
        R12 += bbar*fabs(fp1[j])*fabs(1.-1./(double)k-(double)j/(double)(k-1.));
        temp2 += bbar*fabs(fp1[j]);
        temp3 += bbar*fabs((double)j/(double)(k-1.0)*fp1[j]);
    }
    bias = ((temp2/2.0+temp3)/(double)k+fabs(fp2[0]))/(double)k+R12;
/*printf("R12=%f\n", R12); */
    mse = bias*bias+PI0*(double)k*hk(r,k)/(double)m;
/*printf("mse=%f\n", mse); */
    return mse;
}
double mse_pi0_est(int m, int r, int k, double *ftilde, double *Ftilde){
    int i, j;
    double PI0, mse, bbar, f1hat;
    PI0 = 0.;
    f1hat = k*(Ftilde[k]-Ftilde[k-1]);
/*printf("PI0=%f\n", PI0);   */
    mse = 0.;
    for(j=0; j<k-1; j++){
        bbar = 0.0;
        for(i=1; i<=r; i++) bbar+=Bernstein_polynom(j, k-1, 1.0-(double)i/(double)k);
        bbar /=(double)r;
        mse += bbar*(Ftilde[j+2]-Ftilde[j+1]);
    }
    mse = (k*mse-f1hat)*(k*mse-f1hat)+f1hat*k*hk(r,k)/(double)m;
/*printf("mse=%f\n",mse);    */
    return mse;
}

/*  Choosing optimal r and k  */

/*void est_rk(double *y, int* m, int* r0, int* k0, int* K, int rk[2], int* Trial_r, int* Trial_k, */
/*               bool approx = true, bool Smooth = true, bool Tracking = true){                   */
void est_rk(double *y, int* m, int* r0, int* k0, int* K, int *rk, int* Trial_r, int* Trial_k,
                 int* approx, int* Smooth){
    int i, j, s, t, jj, tr = 0, *rr, *ord;
    int r_opt = *r0, k_opt = *k0, tk = 0, k = *k0;
    double *MSE, min_mse = 100., *tt, *ftilde,*Ftilde, *fp1, *fp2;
    double *edf, *f;

    jj = 0;
/*    printf("r = %d, k=%d\n", rk[0], rk[1]); */
    while(k>3 && tk<=*Trial_k){
        MSE = dvector(0, k+1);
        jj = jj+1;
        rr = ivector(0, k+1);
        i = imin2(r_opt, k-1);
        j = 0; tr = 0;
        tt = dvector(0,k);
        edf = dvector(0,k);
        f = dvector(0,k-1);
        Ftilde = dvector(0,k);
        ftilde = dvector(0,k);
        for(t=0; t<=k; t++) tt[t] = (double)t/(double)k;
        ecdf(y, tt, *m, k+1, edf);
        Bf(edf, k+1, tt, k+1, Ftilde);
        if(*Smooth>0) for(t=0; t<k; t++) f[t] = k*(Ftilde[t+1]-Ftilde[t]);
        else for(t=0; t<k; t++) f[t] = k*(edf[t+1]-edf[t]);
        Bf(f, k, tt, k+1, ftilde);
        fp1 = dvector(0,k);
        fp2 = dvector(0,k);
        if(*approx>0){
            for(t=0; t<k; t++) {
                fp1[t] = 0.0;
                fp2[t] = 0.0;
                for(s=0; s<k-1; s++){
                    fp1[t] += (ftilde[s+1]-ftilde[s])*Bernstein_polynom(s, k-2, (double)t/(double)(k-1.));
                    fp2[t] += (ftilde[s+1]-ftilde[s])*Bernstein_polynom(s, k-2, 1.-(double)(t+1)/(double)k);
                }
                fp1[t] *= (double) (k-1.0);
                fp2[t] *= (double) (k-1.0);
            }
        }
        while(i>1 && tr<=*Trial_r){
            if(*approx>0) MSE[j] = mse_pi0_approx(*m, i, k, ftilde, fp1, fp2);
            else MSE[j] = mse_pi0_est(*m, i, k, ftilde, Ftilde);
            rr[j] = i;
            i = i-1;
            if(j>1 && MSE[j]>MSE[j-1]) tr = tr+1;
            j = j+1;
        }
        i = r_opt+1;
        tr = 0;
        while(i<k && tr<=*Trial_r){
            if(*approx>0) MSE[j] = mse_pi0_approx(*m, i, k, ftilde, fp1, fp2);
            else MSE[j] = mse_pi0_est(*m, i, k, ftilde, Ftilde);
            rr[j] = i;
            i = i+1;
            if(MSE[j]>MSE[j-1]) tr = tr+1;
            j = j+1;
        }
        ord = ivector(0, j-2);
        forder(MSE, j-1, ord);
        if(MSE[ord[0]-1] <= min_mse){
            r_opt = rr[ord[0]-1];
            k_opt = k;
            min_mse = MSE[ord[0]-1];
            tk = 0;
        }
        else tk = tk+1;
        /*if(Tracking)  */
        /*printf("k = %d, j=%d, r_opt = %d, k_opt = %d, Min(MSE) = %f,
            MSE = %f\n", k, j, r_opt, k_opt, min_mse, MSE[ord[0]-1]);  */
        free_dvector(MSE, 0, k+1);
        free_ivector(rr, 0, k+1);
        free_ivector(ord, 0, j-2);
        free_dvector(ftilde, 0,k);
        free_dvector(tt,0,k);
        /*if(*approx>0){  */
            free_dvector(fp1, 0,k);
            free_dvector(fp2, 0,k);
        /*}   */
        free_dvector(Ftilde, 0, k);
        free_dvector(edf, 0,k);
        free_dvector(f, 0,k-1);
        k = k-1;
    }
    k = *k0+1; tk = 0;
    while(k<*K && tk<=*Trial_k){
        MSE = dvector(0, k+1);
        jj = jj+1;
        rr = ivector(0, k+1);
        i = imin2(r_opt, k-1);
        j = 0; tr = 0;
        tt = dvector(0,k);
        edf = dvector(0,k);
        f = dvector(0,k-1);
        Ftilde = dvector(0,k);
        ftilde = dvector(0,k);
        for(t=0;t<=k;t++) tt[t] = (double)t/(double)k;
        ecdf(y, tt, *m, k+1, edf);
        Bf(edf, k+1, tt, k+1, Ftilde);
        if(*Smooth>0) for(t=0;t<k;t++) f[t] = k*(Ftilde[t+1]-Ftilde[t]);
        else for(t=0;t<k;t++) f[t] = k*(edf[t+1]-edf[t]);
        Bf(f, k, tt, k+1, ftilde);
            fp1 = dvector(0,k);
            fp2 = dvector(0,k);
        if(*approx>0){
            for(t=0;t<k;t++) {
                fp1[t] = 0.0;
                fp2[t] = 0.0;
                for(s=0; s<k-1;s++){
                    fp1[t] += (ftilde[s+1]-ftilde[s])*Bernstein_polynom(s, k-2, (double)t/(double)(k-1.));
                    fp2[t] += (ftilde[s+1]-ftilde[s])*Bernstein_polynom(s, k-2, 1.-(double)(t+1)/(double)k);
                }
                fp1[t] *= (double) (k-1.0);
                fp2[t] *= (double) (k-1.0);
            }
        }
        while(i>1 && tr<=*Trial_r){
            if(*approx>0) MSE[j] = mse_pi0_approx(*m, i, k, ftilde, fp1, fp2);
            else MSE[j] = mse_pi0_est(*m, i, k, ftilde, Ftilde);
            rr[j] = i;
            i = i-1;
            if(j>1 && MSE[j]>MSE[j-1]) tr = tr+1;
            j = j+1;
        }
        i = r_opt+1;
        tr = 0;
        while(i<k && tr<=*Trial_r){
            if(*approx>0) MSE[j] = mse_pi0_approx(*m, i, k, ftilde, fp1, fp2);
            else MSE[j] = mse_pi0_est(*m, i, k, ftilde, Ftilde);
            rr[j] = i;
            i = i+1;
            if(MSE[j]>MSE[j-1]) tr = tr+1;
            j = j+1;
        }
        ord = ivector(0, j-2);
        forder(MSE, j-1, ord);
        if(MSE[ord[0]-1]<=min_mse){
            r_opt = rr[ord[0]-1];
            k_opt = k;
            min_mse = MSE[ord[0]-1];
            tk = 0;
        }
        else tk = tk+1;
        /*if(Tracking) */
        /*  printf("k = %d, j=%d, r_opt = %d, k_opt = %d, Min(MSE) = %f,
                MSE = %f\n", k, j, r_opt, k_opt, min_mse,  MSE[ord[0]-1]);  */
        free_dvector(ftilde, 0,k);
        free_dvector(MSE, 0, k+1);
        free_ivector(rr, 0, k+1);
        free_ivector(ord, 0, j-2);
        free_dvector(tt,0,k);
        /*if(*approx>0){   */
            free_dvector(fp1, 0,k);
            free_dvector(fp2, 0,k);
        /*}        */
        free_dvector(Ftilde, 0, k);
        free_dvector(edf, 0,k);
        free_dvector(f, 0,k-1);
        k = k+1;
    }
    rk[0] = r_opt;
    rk[1] = k_opt;
    printf("r = %d, k=%d\n", rk[0], rk[1]);
}

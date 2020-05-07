/*******************************************************************************************

Copyright (c) 2019 Neil Cornish

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_spline.h>

#include "Header.h"
#include "Constants.h"

double cosw(double t, double Tint);

#define PI 3.141592653589793

// gcc -o segmentAET segmentAET.c -lm -lgsl

int main(int argc, char *argv[])
{
  int i, j, k, N, M, is;
  int Ns, Ng;
  double fac;
  double t, f, x, Tend;
  char filename[1024];
  char command[1024];
  char strn[300];
  double alpha, s1, s2, Tcut, t_tuke;
  double *A, *E, *T, *Times;
  double *AC, *EC;
  double *SA, *SE, *ST, *SN, *SNS;
  double *X, *Y, *Z;
  double SAE, SXYZ;

  FILE *in;
  FILE *out;

    N = NFFT;

    double Tseg;
    int Nseg, segs, seg;
    
    t_tuke = 1.0e5;
    
    segs = 16;
    
    Nseg = N/segs;
    
    Tseg = Tobs/(double)(segs);
    
      Times = (double*)malloc(sizeof(double)*(N));
      X = (double*)malloc(sizeof(double)*(N));
      Y = (double*)malloc(sizeof(double)*(N));
      Z = (double*)malloc(sizeof(double)*(N));
      A = (double*)malloc(sizeof(double)*(N));
      E = (double*)malloc(sizeof(double)*(N));
      T = (double*)malloc(sizeof(double)*(N));
      SN = (double*)malloc(sizeof(double)*(N));
    
         // Read in the data
         in = fopen("LDC1-1_MBHB_v2_TD.txt","r");
        // in = fopen("LDC1-1_MBHB_v2_TD_noiseless.txt","r");
         fgets(strn, 100, in);  // strip comment line
         printf("\n %s\n", strn);
         for (i = 0; i < N; ++i)
         {
            fscanf(in,"%lf%lf%lf%lf", &Times[i], &X[i], &Y[i], &Z[i]);
         }
          fclose(in);
    
    alpha = (2.0*t_tuke/Tseg);
    
    fac = Tseg/(double)(Nseg/2);   // Fourier scaling
    
   for (j = 0; j < segs; ++j)
   {
       
     for (i = 0; i < Nseg; ++i)
      {
       k = i+j*Nseg;
       A[i] = (2.0*X[k]-Y[k]-Z[k])/3.0;
       E[i] = (Z[k]-Y[k])/sqrt(3.0);
       T[i] = (X[k]+Y[k]+Z[k])/3.0;
      }
       
          tukey(A, alpha, Nseg);
          tukey(E, alpha, Nseg);
          tukey(T, alpha, Nseg);
       
          gsl_fft_real_radix2_transform(A, 1, Nseg);
          gsl_fft_real_radix2_transform(E, 1, Nseg);
          gsl_fft_real_radix2_transform(T, 1, Nseg);
       
       for(i=0; i< Nseg; i++)
       {
           A[i] *= fac;
           E[i] *= fac;
           T[i] *= fac;
       }
       
       sprintf(command, "AET_seg%d_f.dat", j);
       out = fopen(command,"w");
       SN[0] = 1.0;
       SAE = 1.0;
       for(i=0; i< Nseg; i++)
       {
           f = (double)(i)/Tseg;
           
           if(i > 0 && i < Nseg/2)
           {
               instrument_noise(f, &SAE, &SXYZ);
               SN[i] = SAE;
           }

           fprintf(out,"%e %.12e %.12e %.12e %.12e\n", f, A[i], E[i], T[i], SAE);
       }
       fclose(out);
       
   }
    
     alpha = (2.0*t_tuke/Tobs);
     fac = Tobs/(double)(N/2);   // Fourier scaling
    
    for (i = 0; i < N; ++i)
    {
        A[i] = (2.0*X[i]-Y[i]-Z[i])/3.0;
        E[i] = (Z[i]-Y[i])/sqrt(3.0);
        T[i] = (X[i]+Y[i]+Z[i])/3.0;
    }
    
       tukey(A, alpha, N);
       tukey(E, alpha, N);
       tukey(T, alpha, N);
    
    // use X, Y, Z to hold a copy
    for (i = 0; i < N; ++i)
    {
        X[i] = A[i];
        Y[i] = E[i];
        Z[i] = T[i];
    }
    
    // if we want to cut the data off
    Tend = 2.492000e+07;
    Tcut = Tend+t_rise/1.5;
    
    j = (int)((Tcut-t_rise)/dt)-2;
    k = (int)((Tcut)/dt)+2;
    if(k > NFFT) k = NFFT;
    
      for (i = j; i < k; ++i)
       {
           t = (double)(i)*dt;
           x = 0.5*(1.0-cos(PI*(t-Tcut)/t_rise));
           if(t < Tcut-t_rise) x = 1.0;
           if(t > Tcut) x = 0.0;
           X[i] *= x;
           Y[i] *= x;
           Z[i] *= x;
          // printf("%e %f\n", t, x);
       }
    
    for (i = k; i < NFFT; ++i)
     {
         X[i] = 0.0;
         Y[i] = 0.0;
         Z[i] = 0.0;
     }
    
    
    out = fopen("cut.dat","w");
    for (i = 0; i < NFFT; ++i)
    {
        fprintf(out,"%e %e %e\n", (double)(i)*dt, X[i], Y[i]);
    }
    fclose(out);
    
    out = fopen("full.dat","w");
    for (i = 0; i < NFFT; ++i)
    {
        fprintf(out,"%e %e %e\n", (double)(i)*dt, A[i], E[i]);
    }
    fclose(out);
    
    
    
          gsl_fft_real_radix2_transform(X, 1, N);
          gsl_fft_real_radix2_transform(Y, 1, N);
          gsl_fft_real_radix2_transform(Z, 1, N);
    
    
       gsl_fft_real_radix2_transform(A, 1, N);
       gsl_fft_real_radix2_transform(E, 1, N);
       gsl_fft_real_radix2_transform(T, 1, N);
    
    for(i=0; i< N; i++)
    {
        A[i] *= fac;
        E[i] *= fac;
        T[i] *= fac;
        X[i] *= fac;
        Y[i] *= fac;
        Z[i] *= fac;
    }
    
    out = fopen("AET_f.dat","w");
    SN[0] = 1.0;
    SAE = 1.0;
    
    for(i=0; i< N; i++)
    {
        f = (double)(i)/Tobs;
        
        if(i > 0 && i < N/2)
        {
            instrument_noise(f, &SAE, &SXYZ);
            SN[i] = SAE;
        }

        fprintf(out,"%e %.12e %.12e %.12e %.12e\n", f, A[i], E[i], T[i], SAE);
    }
    fclose(out);
    
    out = fopen("AET_f_cut.dat","w");
    SN[0] = 1.0;
    SAE = 1.0;
    for(i=0; i< N; i++)
    {
        f = (double)(i)/Tobs;
        
        if(i > 0 && i < N/2)
        {
            instrument_noise(f, &SAE, &SXYZ);
            SN[i] = SAE;
        }
        
       

        fprintf(out,"%e %.12e %.12e %.12e %.12e\n", f, X[i], Y[i], Z[i], SAE);
    }
    fclose(out);
                  
    x = 0.0;
    for(i=0; i< N/2; i++)
    {
        x += 4.0*(X[i]*X[i]+X[N-i]*X[N-i]+Y[i]*Y[i]+Y[N-i]*Y[N-i])/SN[i];
    }
    printf("cut SNR = %f\n", sqrt(x/Tobs));
    
    j = (int)(Tobs*2.0e-4);
    k = (int)(Tobs*6.0e-4);
    
    
      out = fopen("ratio.dat","w");
       for(i=j; i< k; i++)
       {
           f = (double)(i)/Tobs;
           
           fprintf(out,"%e %.12e %.12e\n", f, sqrt((X[i]*X[i]+X[N-i]*X[N-i])/(A[i]*A[i]+A[N-i]*A[N-i])),
               sqrt((Y[i]*Y[i]+Y[N-i]*Y[N-i])/(E[i]*E[i]+E[N-i]*E[N-i])));
       }
       fclose(out);
     
    return 0;

}

void instrument_noise(double f, double *SAE, double *SXYZ)
{
    //Power spectral density of the detector noise and transfer frequency
    double Sn, red, confusion_noise;
    double Sloc, fonfs;
    double f1, f2;
    double A1, A2, slope, LC;
    FILE *outfile;
    
    fonfs = f/fstar;
    
    LC = 16.0*fonfs*fonfs;
    
     red = 16.0*((1.0e-4/f)*(1.0e-4/f));
   // red = 0.0;
    
    
    // Calculate the power spectral density of the detector noise at the given frequency
    
    *SAE = LC*16.0/3.0*pow(sin(fonfs),2.0)*( (2.0+cos(fonfs))*(Sps) + 2.0*(3.0+2.0*cos(fonfs)+cos(2.0*fonfs))*(Sacc/pow(2.0*PI*f,4.0)*(1.0+red)) ) / pow(2.0*Larm,2.0);
    
    *SXYZ = LC*4.0*pow(sin(fonfs),2.0)*( 4.0*(Sps) + 8.0*(1.0+pow(cos(fonfs),2.0))*(Sacc/pow(2.0*PI*f,4.0)*(1.0+red)) ) / pow(2.0*Larm,2.0);
    
}



void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2)
{
    /* Butterworth bandpass filter
     n = filter order 4,8,12,...
     s = sampling frequency
     f1 = upper half power frequency
     f2 = lower half power frequency  */
    
    if(n % 4){ printf("Order must be 4,8,12,16,...\n"); return;}
    
    int i, j;
    double a = cos(PI*(f1+f2)/s)/cos(PI*(f1-f2)/s);
    double a2 = a*a;
    double b = tan(PI*(f1-f2)/s);
    double b2 = b*b;
    double r;
    
    n = n/4;
    double *A = (double *)malloc(n*sizeof(double));
    double *d1 = (double *)malloc(n*sizeof(double));
    double *d2 = (double *)malloc(n*sizeof(double));
    double *d3 = (double *)malloc(n*sizeof(double));
    double *d4 = (double *)malloc(n*sizeof(double));
    double *w0 = (double *)malloc(n*sizeof(double));
    double *w1 = (double *)malloc(n*sizeof(double));
    double *w2 = (double *)malloc(n*sizeof(double));
    double *w3 = (double *)malloc(n*sizeof(double));
    double *w4 = (double *)malloc(n*sizeof(double));
    double x;
    
    for(i=0; i<n; ++i)
    {
        r = sin(PI*(2.0*(double)i+1.0)/(4.0*(double)n));
        s = b2 + 2.0*b*r + 1.0;
        A[i] = b2/s;
        d1[i] = 4.0*a*(1.0+b*r)/s;
        d2[i] = 2.0*(b2-2.0*a2-1.0)/s;
        d3[i] = 4.0*a*(1.0-b*r)/s;
        d4[i] = -(b2 - 2.0*b*r + 1.0)/s;
        w0[i] = 0.0;
        w1[i] = 0.0;
        w2[i] = 0.0;
        w3[i] = 0.0;
        w4[i] = 0.0;
    }
    
    for(j=0; j< M; ++j)
    {
        if(fwrv == 1) x = in[j];
        if(fwrv == -1) x = in[M-j-1];
        for(i=0; i<n; ++i)
        {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i]+ d3[i]*w3[i]+ d4[i]*w4[i] + x;
            x = A[i]*(w0[i] - 2.0*w2[i] + w4[i]);
            w4[i] = w3[i];
            w3[i] = w2[i];
            w2[i] = w1[i];
            w1[i] = w0[i];
        }
        if(fwrv == 1) out[j] = x;
        if(fwrv == -1) out[M-j-1] = x;
    }
    
    free(A);
    free(d1);
    free(d2);
    free(d3);
    free(d4);
    free(w0);
    free(w1);
    free(w2);
    free(w3);
    free(w4);
    
    return;
}

double tukeyt(double t, double tr, double Tint)
{
    double filter;
    
    filter = 1.0;
    if(t < tr) filter = 0.5*(1.0+cos(PI*(t/tr-1.0)));
    if(t > (Tint-tr)) filter = 0.5*(1.0+cos(PI*((Tint-t)/tr-1.0 )));
    if(t < 0.0) filter = 0.0;
    if(t > Tint) filter = 0.0;
    
    return filter;
    
}

double cosw(double t, double Tint)
{
    double filter;
    
    filter = sin(PI*(t/Tint));
    if(t < 0) filter = 0.0;
    if(t > Tint) filter = 0.0;
    
    return filter;
    
}


void tukey(double *data, double alpha, int N)
{
    int i, imin, imax;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i < imin) filter = 0.5*(1.0+cos(PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(PI*( (double)(i)/(double)(imin)-2.0/alpha+1.0 )));
        data[i] *= filter;
    }
    
}

void tukey_scale(double *s1, double *s2, double alpha, int N)
{
    int i, imin, imax;
    double x1, x2;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    x1 = 0.0;
    x2 = 0.0;
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i < imin) filter = 0.5*(1.0+cos(PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(PI*( (double)(i)/(double)(imin)-2.0/alpha+1.0 )));
        x1 += filter;
        x2 += filter*filter;
    }
    x1 /= (double)(N);
    x2 /= (double)(N);
    
    *s1 = x1;
    *s2 = x2;
    
}



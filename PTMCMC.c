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
#include <stdlib.h>
#include <math.h>
#include "Constants.h"
#include "IMRPhenomD.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <time.h>

#include "Declarations.h"

#ifndef _OPENMP
#define omp ignore
#endif

// OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o  PTMCMC PTMCMC.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

// Linux
// gcc -std=gnu99 -fopenmp -w -o PTMCMC PTMCMC.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm




//##############################################
//MT modifications

gsl_rng **rvec;
//##############################################

int main(int argc,char **argv)
{

  double f, fdot, theta, phi, A, iota, psi, phase;
  char Gfile[50];
  double *params, *pnew;
  double *AS, *ES;
  double *AQ, *EQ;
  double AR, AI, ER, EI;
  double fonfs, Sn, Sm, Acut;
  double Aar, Aai;
  double x, y, z;
  long M, N, q;
  long i, j, k, cnt, mult;
  double Mc, fstart, fstop, fr;
  double SNR;
  double HH, HD, HDQ, DD, Match, MX, ts, ps;
  double HHA, HDA, DDA, HHE, HDE, DDE, MA, ME;
  double Tend;
    
  double *SN;
    
  double *AC, *EC, *TC;
    
    
    clock_t start, end;
    double cpu_time_used;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
  FILE* in;
  FILE* out;
    
    //##############################################
    //open MP modifications
    omp_set_num_threads(NC);
    rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (NC+1));
    for(i = 0 ; i<= NC; i++){
        rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(rvec[i] , i);
    }
    //##############################################
    
    
    params = (double*)malloc(sizeof(double)* (NP));
    SN = (double*)malloc(sizeof(double)* (NFFT));  // goes out further in frequency than we need.
    AC = (double*)malloc(sizeof(double)* (NFFT));
    EC = (double*)malloc(sizeof(double)* (NFFT));
    TC = (double*)malloc(sizeof(double)* (NFFT));
    

    
    if(nflag == 1)
    {
    // start using output from  search
    in = fopen("skymax.dat","r");
    fscanf(in,"%lf", &x);
    for(i=0; i< NP; i++) fscanf(in,"%lf", &params[i]);
    fclose(in);
    }
    else  // no noise
    {
    // start using true values
    in = fopen("truths.dat","r");
    for(i=0; i< NP; i++) fscanf(in,"%lf", &params[i]);
    fclose(in);
    SNRFast(params);
    }
    
    
     // Read in FFTed LDC data and A,E PSD from segmentAET.c
    
      if(cut == 1)
      {
         in = fopen("AET_f_cut.dat","r");
         for(i=0; i< NFFT; i++)
         {
             fscanf(in,"%lf%lf%lf%lf%lf\n", &f, &AC[i], &EC[i], &TC[i], &SN[i]);
         }
         fclose(in);
         Tend = 2.492000e+07;
      }
     else
     {
         in = fopen("AET_f.dat","r");
         for(i=0; i< NFFT; i++)
         {
             fscanf(in,"%lf%lf%lf%lf%lf\n", &f, &AC[i], &EC[i], &TC[i], &SN[i]);
         }
         fclose(in);
         Tend = Tobs;
     }
    
    //FisherPlot(0, Tend, params);
    
    //return 1;
    
    MCMC(params, Tend, NFFT, AC, EC, SN);
    
    //###############################################
    //MT modification
   	for(i =0 ;i<= NC; i++){
        gsl_rng_free(rvec[i]);
    }
    free(rvec);
    //###############################################
    
    free(params);
    free(SN);
    free(AC);
    free(EC);
    free(TC);
    

    return 1;
    
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


double fourier_nwip(double *a, double *b, double *Sn, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-i;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product/(Sn[i]);
    }
    
    return(4.0*arg);
    
}

void fourier_nwip_dual_time(double *abt, double *aA, double *bA, double *aE, double *bE, double *Sn, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    double scale;
    
    scale = 2.0*(double)(n)/Tobs;
    
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-i;
        ReA = aA[j]; ImA = aA[k];
        ReB = bA[j]; ImB = bA[k];
        abt[j] = scale*(ReA*ReB + ImA*ImB)/Sn[i];
        abt[k] = scale*(ImA*ReB - ReA*ImB)/Sn[i];
        ReA = aE[j]; ImA = aE[k];
        ReB = bE[j]; ImB = bE[k];
        abt[j] += scale*(ReA*ReB + ImA*ImB)/Sn[i];
        abt[k] += scale*(ImA*ReB - ReA*ImB)/Sn[i];
    }
    
    abt[0]=0.0;
    abt[n/2]=0.0;
    
    gsl_fft_halfcomplex_radix2_inverse(abt, 1, n);
    
}

void fourier_nwip_time(double *abt, double *a, double *b, double *Sn, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    double scale;
    
    scale = 2.0*(double)(n)/Tobs;
    
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-i;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        abt[j] = scale*(ReA*ReB + ImA*ImB)/Sn[i];
        abt[k] = scale*(ImA*ReB - ReA*ImB)/Sn[i];
    }
    
    abt[0]=0.0;
    abt[n/2]=0.0;
    
    gsl_fft_halfcomplex_radix2_inverse(abt, 1, n);
    
}


// merger time at guiding center of detector
double Tmerger(double *params, double t)
{
    
    /*   Indicies   */
    int i,j, k, n, m, a, M;
    
    /*   Gravitational Wave basis vectors   */
    double *kv;
    
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;

    /*   Dot products   */
    double kdotx;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double xa, ya, za;
    
    /*   Miscellaneous  */
    double xm, td;
    
    /*   Allocating Arrays   */
    
    kv = double_vector(4);
    x = double_vector(4); y = double_vector(4); z = double_vector(4);

    phi = params[8];   // EclipticLongitude
    //Calculate cos and sin of sky position
    costh = params[7];   sinth = sqrt(1.0-costh*costh);
    cosph = cos(phi);     sinph = sin(phi);

    kv[1] = -sinth*cosph;
    kv[2] = -sinth*sinph;
    kv[3] = -costh;
    
    spacecraft(t, x, y, z);
        
        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
    kdotx = (xa*kv[1]+ya*kv[2]+za*kv[3])/clight;
    
    td = t+kdotx;
    
    free_double_vector(kv);
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    
    return(td);
}

void RAantenna(double *params, int NF, double *TF, double *FF, double *xi, double *FpAR, double *FpAI, double *FcAR, double *FcAI,
               double *FpER, double *FpEI, double *FcER, double *FcEI)
{
    
    /*   Indicies   */
    int i,j, k, n, m, a, M;
    
    /*   Gravitational Wave basis vectors   */
    double *u,*v,*kv;
    
    /*   Polarization basis tensors   */
    double **eplus, **ecross;
    
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;
    double *r12, *r13, *r21, *r23, *r31, *r32;
    double *r10, *r20, *r30;
    double *vx, *vy, *vz;
    
    double q1, q2, q3, q4;
    
    /*   Dot products   */
    double kdotx;
    
    /*   Convenient quantities   */
    double **dplus, **dcross;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double t, xa, ya, za;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx;
    
    double delt;
    
    double fpx, fcx, fpy, fcy, fpz, fcz;
    
    double **TR, **TI, **kdr;
    
    double *kdg;
    
    double fr;
    
    /*   Allocating Arrays   */
    
    u = double_vector(4); v = double_vector(4); kv = double_vector(4);
    
    eplus  = double_matrix(4,4); ecross = double_matrix(4,4);
    
    dplus  = double_matrix(4,4); dcross = double_matrix(4,4);
    
    TR  = double_matrix(4,4); TI = double_matrix(4,4);
    
    kdr = double_matrix(4,4);
    
    kdg = double_vector(4);
    
    x = double_vector(4); y = double_vector(4); z = double_vector(4);
    
    r12 = double_vector(4); r21 = double_vector(4); r31 = double_vector(4);
    r13 = double_vector(4); r23 = double_vector(4); r32 = double_vector(4);
    r10 = double_vector(4); r20 = double_vector(4); r30 = double_vector(4);
    
    phi = params[8];   // EclipticLongitude
    psi = params[9];   // polarization
    //Calculate cos and sin of sky position, inclination, polarization
    costh = params[7];   sinth = sqrt(1.0-costh*costh);
    cosph = cos(phi);     sinph = sin(phi);
    cosps = cos(2.*psi);  sinps = sin(2.*psi);
    
    
    /*   Tensor basis  */
    v[1] =  -costh*cosph;
    v[2] =  -costh*sinph;
    v[3] = sinth;
    
    u[1] =  sinph;
    u[2] = -cosph;
    u[3] =  0.;
    
    
    kv[1] = -sinth*cosph;
    kv[2] = -sinth*sinph;
    kv[3] = -costh;
    
    
    
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }
    
    /*   Main Loop   */
    for(n=0; n< NF; n++)
    {
        // Barycenter time
        t = TF[n];
        fr = FF[n]/(2.0*fstar);
        
        spacecraft(t, x, y, z);
        
        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
        kdotx = (xa*kv[1]+ya*kv[2]+za*kv[3])/clight;
        
        // detector time and frequency
        xi[n]  = t - kdotx;
        
        //Unit separation vector from spacecraft i to j
        r12[1] = (x[2] - x[1])/Larm;   r13[1] = (x[3] - x[1])/Larm;   r23[1] = (x[3] - x[2])/Larm;
        r12[2] = (y[2] - y[1])/Larm;   r13[2] = (y[3] - y[1])/Larm;   r23[2] = (y[3] - y[2])/Larm;
        r12[3] = (z[2] - z[1])/Larm;   r13[3] = (z[3] - z[1])/Larm;   r23[3] = (z[3] - z[2])/Larm;
        
        // These are not unit vectors. Just pulling out the Larm scaling
        r10[1] = (xa-x[1])/Larm;   r10[2] = (ya-y[1])/Larm;  r10[3] = (za-z[1])/Larm;
        r20[1] = (xa-x[2])/Larm;   r20[2] = (ya-y[2])/Larm;  r20[3] = (za-z[2])/Larm;
        r30[1] = (xa-x[3])/Larm;   r30[2] = (ya-y[3])/Larm;  r30[3] = (za-z[3])/Larm;
        
        kdr[1][2] = 0.0;
        for(k=1; k<=3; k++) kdr[1][2] += kv[k]*r12[k];
        kdr[1][3] = 0.0;
        for(k=1; k<=3; k++) kdr[1][3] += kv[k]*r13[k];
        kdr[2][3] = 0.0;
        for(k=1; k<=3; k++) kdr[2][3] += kv[k]*r23[k];
        
        kdr[2][1] = -kdr[1][2];  kdr[3][1] = -kdr[1][3];  kdr[3][2] = -kdr[2][3];
        
        kdg[1] = 0.0;
        for(k=1; k<=3; k++) kdg[1] += kv[k]*r10[k];
        kdg[2] = 0.0;
        for(k=1; k<=3; k++) kdg[2] += kv[k]*r20[k];
        kdg[3] = 0.0;
        for(k=1; k<=3; k++) kdg[3] += kv[k]*r30[k];
        
        //Make use of symmetry
        for(i=1; i<=3; i++)
        {
            r21[i] = -r12[i];
            r31[i] = -r13[i];
            r32[i] = -r23[i];
        }
        
        
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                q1 = fr*(1.0-kdr[i][j]);
                q2 = fr*(1.0+kdr[i][j]);
                q3 = -fr*(3.0+kdr[i][j]-2.0*kdg[i]);
                q4 = -fr*(1.0+kdr[i][j]-2.0*kdg[i]);
                q1 = (sin(q1)/q1);
                q2 = (sin(q2)/q2);
                TR[i][j] = 0.5*(q1*cos(q3)+q2*cos(q4));   // goes to 1 when f/fstat small
                TI[i][j] = 0.5*(q1*sin(q3)+q2*sin(q4));   // goes to 0 when f/fstat small
            }
        }
        

        
        dplus[1][2] = dplus[1][3] = dplus[2][1] = dplus[2][3] = dplus[3][1] = dplus[3][2] = 0.;
        dcross[1][2] = dcross[1][3] = dcross[2][1] = dcross[2][3] = dcross[3][1] = dcross[3][2] = 0.;
        //Convenient quantities d+ & dx
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                dplus[1][2]  += r12[i]*r12[j]*eplus[i][j];   dcross[1][2] += r12[i]*r12[j]*ecross[i][j];
                dplus[2][3]  += r23[i]*r23[j]*eplus[i][j];   dcross[2][3] += r23[i]*r23[j]*ecross[i][j];
                dplus[1][3]  += r13[i]*r13[j]*eplus[i][j];   dcross[1][3] += r13[i]*r13[j]*ecross[i][j];
            }
        }
        
        dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
        dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
        dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
        
        fpx = -0.5*( (dplus[1][2]*cosps+dcross[1][2]*sinps)*TR[1][2] - (dplus[1][3]*cosps+dcross[1][3]*sinps)*TR[1][3] );
        fcx = -0.5*( (-dplus[1][2]*sinps+dcross[1][2]*cosps)*TR[1][2] - (-dplus[1][3]*sinps + dcross[1][3]*cosps)*TR[1][3] );
        
        fpy = -0.5*( (dplus[2][3]*cosps+dcross[2][3]*sinps)*TR[2][3] - (dplus[2][1]*cosps+dcross[2][1]*sinps)*TR[2][1] );
        fcy = -0.5*( (-dplus[2][3]*sinps+dcross[2][3]*cosps)*TR[2][3] - (-dplus[2][1]*sinps + dcross[2][1]*cosps)*TR[2][1] );
        
        fpz = -0.5*( (dplus[3][1]*cosps+dcross[3][1]*sinps)*TR[3][1] - (dplus[3][2]*cosps+dcross[3][2]*sinps)*TR[3][2] );
        fcz = -0.5*( (-dplus[3][1]*sinps+dcross[3][1]*cosps)*TR[3][1] - (-dplus[3][2]*sinps + dcross[3][2]*cosps)*TR[3][2] );

        FpAR[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAR[n] = (2.0*fcx-fcy-fcz)/3.0;
        
        FpER[n] = (fpz-fpy)/sq3;
        FcER[n] = (fcz-fcy)/sq3;
                   
        fpx = -0.5*( (dplus[1][2]*cosps+dcross[1][2]*sinps)*TI[1][2] - (dplus[1][3]*cosps+dcross[1][3]*sinps)*TI[1][3] );
        fcx = -0.5*( (-dplus[1][2]*sinps+dcross[1][2]*cosps)*TI[1][2] - (-dplus[1][3]*sinps + dcross[1][3]*cosps)*TI[1][3] );
                                         
        fpy = -0.5*( (dplus[2][3]*cosps+dcross[2][3]*sinps)*TI[2][3] - (dplus[2][1]*cosps+dcross[2][1]*sinps)*TI[2][1] );
        fcy = -0.5*( (-dplus[2][3]*sinps+dcross[2][3]*cosps)*TI[2][3] - (-dplus[2][1]*sinps + dcross[2][1]*cosps)*TI[2][1] );
                                                               
        fpz = -0.5*( (dplus[3][1]*cosps+dcross[3][1]*sinps)*TI[3][1] - (dplus[3][2]*cosps+dcross[3][2]*sinps)*TI[3][2] );
        fcz = -0.5*( (-dplus[3][1]*sinps+dcross[3][1]*cosps)*TI[3][1] - (-dplus[3][2]*sinps + dcross[3][2]*cosps)*TI[3][2] );
                                                                                     
        FpAI[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAI[n] = (2.0*fcx-fcy-fcz)/3.0;
                                                                                     
        FpEI[n] = (fpz-fpy)/sq3;
        FcEI[n] = (fcz-fcy)/sq3;
        
        
    }

    
    free_double_vector(u); free_double_vector(v); free_double_vector(kv);
    
    free_double_matrix(eplus,4); free_double_matrix(ecross,4);
    
    free_double_matrix(dplus,4); free_double_matrix(dcross,4);
    
    free_double_matrix(TR,4); free_double_matrix(TI,4);
    
    free_double_matrix(kdr,4);
    
    free_double_vector(kdg);
    
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    
    free_double_vector(r12); free_double_vector(r21); free_double_vector(r31);
    free_double_vector(r13); free_double_vector(r23); free_double_vector(r32);
    free_double_vector(r10); free_double_vector(r20); free_double_vector(r30);
    
    return;
}

void RAfilters(double *params, int NF, double *TF, double *FF, double *xi, double *FpAR, double *FpAI, double *FcAR, double *FcAI,
               double *FpER, double *FpEI, double *FcER, double *FcEI)
{
    
    /*   Indicies   */
    int i,j, k, n, m, a, M;
    
    /*   Gravitational Wave basis vectors   */
    double *u,*v,*kv;
    
    /*   Polarization basis tensors   */
    double **eplus, **ecross;
    
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;
    double *r12, *r13, *r21, *r23, *r31, *r32;
    double *r10, *r20, *r30;
    double *vx, *vy, *vz;
    
    double q1, q2, q3, q4;
    
    /*   Dot products   */
    double kdotx;
    
    /*   Convenient quantities   */
    double **dplus, **dcross;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double t, xa, ya, za;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx;
    
    double delt;
    
    double fpx, fcx, fpy, fcy, fpz, fcz;
    
    double **TR, **TI, **kdr;
    
    double *kdg;
    
    double fr;
    
    /*   Allocating Arrays   */
    
    u = double_vector(4); v = double_vector(4); kv = double_vector(4);
    
    eplus  = double_matrix(4,4); ecross = double_matrix(4,4);
    
    dplus  = double_matrix(4,4); dcross = double_matrix(4,4);
    
    TR  = double_matrix(4,4); TI = double_matrix(4,4);
    
    kdr = double_matrix(4,4);
    
    kdg = double_vector(4);
    
    x = double_vector(4); y = double_vector(4); z = double_vector(4);
    
    r12 = double_vector(4); r21 = double_vector(4); r31 = double_vector(4);
    r13 = double_vector(4); r23 = double_vector(4); r32 = double_vector(4);
    r10 = double_vector(4); r20 = double_vector(4); r30 = double_vector(4);

    
    phi = params[8];   // EclipticLongitude
    //Calculate cos and sin of sky position, inclination, polarization
    costh = params[7];   sinth = sqrt(1.0-costh*costh);
    cosph = cos(phi);     sinph = sin(phi);
    
    
    /*   Tensor basis  */
    v[1] =  -costh*cosph;
    v[2] =  -costh*sinph;
    v[3] = sinth;
    
    u[1] =  sinph;
    u[2] = -cosph;
    u[3] =  0.;
    
    
    kv[1] = -sinth*cosph;
    kv[2] = -sinth*sinph;
    kv[3] = -costh;
    
    
    
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }
    
    /*   Main Loop   */
    for(n=0; n< NF; n++)
    {
        // Barycenter time
        t = TF[n];
        fr = FF[n]/(2.0*fstar);
        
        spacecraft(t, x, y, z);
        
        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
        kdotx = (xa*kv[1]+ya*kv[2]+za*kv[3])/clight;
        
        // detector time and frequency
        xi[n]  = t - kdotx;
        
        //Unit separation vector from spacecraft i to j
        r12[1] = (x[2] - x[1])/Larm;   r13[1] = (x[3] - x[1])/Larm;   r23[1] = (x[3] - x[2])/Larm;
        r12[2] = (y[2] - y[1])/Larm;   r13[2] = (y[3] - y[1])/Larm;   r23[2] = (y[3] - y[2])/Larm;
        r12[3] = (z[2] - z[1])/Larm;   r13[3] = (z[3] - z[1])/Larm;   r23[3] = (z[3] - z[2])/Larm;
        
        // These are not unit vectors. Just pulling out the Larm scaling
        r10[1] = (xa-x[1])/Larm;   r10[2] = (ya-y[1])/Larm;  r10[3] = (za-z[1])/Larm;
        r20[1] = (xa-x[2])/Larm;   r20[2] = (ya-y[2])/Larm;  r20[3] = (za-z[2])/Larm;
        r30[1] = (xa-x[3])/Larm;   r30[2] = (ya-y[3])/Larm;  r30[3] = (za-z[3])/Larm;
        
        kdr[1][2] = 0.0;
        for(k=1; k<=3; k++) kdr[1][2] += kv[k]*r12[k];
        kdr[1][3] = 0.0;
        for(k=1; k<=3; k++) kdr[1][3] += kv[k]*r13[k];
        kdr[2][3] = 0.0;
        for(k=1; k<=3; k++) kdr[2][3] += kv[k]*r23[k];
        
        kdr[2][1] = -kdr[1][2];  kdr[3][1] = -kdr[1][3];  kdr[3][2] = -kdr[2][3];
        
        kdg[1] = 0.0;
        for(k=1; k<=3; k++) kdg[1] += kv[k]*r10[k];
        kdg[2] = 0.0;
        for(k=1; k<=3; k++) kdg[2] += kv[k]*r20[k];
        kdg[3] = 0.0;
        for(k=1; k<=3; k++) kdg[3] += kv[k]*r30[k];
        
        //Make use of symmetry
        for(i=1; i<=3; i++)
        {
            r21[i] = -r12[i];
            r31[i] = -r13[i];
            r32[i] = -r23[i];
        }
        
        
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                q1 = fr*(1.0-kdr[i][j]);
                q2 = fr*(1.0+kdr[i][j]);
                q3 = -fr*(3.0+kdr[i][j]-2.0*kdg[i]);
                q4 = -fr*(1.0+kdr[i][j]-2.0*kdg[i]);
                q1 = (sin(q1)/q1);
                q2 = (sin(q2)/q2);
                TR[i][j] = 0.5*(q1*cos(q3)+q2*cos(q4));   // goes to 1 when f/fstat small
                TI[i][j] = 0.5*(q1*sin(q3)+q2*sin(q4));   // goes to 0 when f/fstat small
            }
        }
        
        
        
        dplus[1][2] = dplus[1][3] = dplus[2][1] = dplus[2][3] = dplus[3][1] = dplus[3][2] = 0.;
        dcross[1][2] = dcross[1][3] = dcross[2][1] = dcross[2][3] = dcross[3][1] = dcross[3][2] = 0.;
        //Convenient quantities d+ & dx
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                dplus[1][2]  += r12[i]*r12[j]*eplus[i][j];   dcross[1][2] += r12[i]*r12[j]*ecross[i][j];
                dplus[2][3]  += r23[i]*r23[j]*eplus[i][j];   dcross[2][3] += r23[i]*r23[j]*ecross[i][j];
                dplus[1][3]  += r13[i]*r13[j]*eplus[i][j];   dcross[1][3] += r13[i]*r13[j]*ecross[i][j];
            }
        }
        
        dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
        dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
        dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
        
        fpx = -0.5*( dplus[1][2]*TR[1][2] - dplus[1][3]*TR[1][3] );
        fcx = -0.5*( dcross[1][2]*TR[1][2] - dcross[1][3]*TR[1][3] );
        
        fpy = -0.5*( dplus[2][3]*TR[2][3] - dplus[2][1]*TR[2][1] );
        fcy = -0.5*( dcross[2][3]*TR[2][3] - dcross[2][1]*TR[2][1] );
        
        fpz = -0.5*( dplus[3][1]*TR[3][1] - dplus[3][2]*TR[3][2] );
        fcz = -0.5*( dcross[3][1]*TR[3][1] - dcross[3][2]*TR[3][2] );
        
        FpAR[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAR[n] = (2.0*fcx-fcy-fcz)/3.0;
        
        FpER[n] = (fpz-fpy)/sq3;
        FcER[n] = (fcz-fcy)/sq3;
        
        fpx = -0.5*( dplus[1][2]*TI[1][2] - dplus[1][3]*TI[1][3] );
        fcx = -0.5*( dcross[1][2]*TI[1][2] - dcross[1][3]*TI[1][3] );
        
        fpy = -0.5*( dplus[2][3]*TI[2][3] - dplus[2][1]*TI[2][1] );
        fcy = -0.5*( dcross[2][3]*TI[2][3] - dcross[2][1]*TI[2][1] );
        
        fpz = -0.5*( dplus[3][1]*TI[3][1] - dplus[3][2]*TI[3][2] );
        fcz = -0.5*( dcross[3][1]*TI[3][1] - dcross[3][2]*TI[3][2] );
        
        FpAI[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAI[n] = (2.0*fcx-fcy-fcz)/3.0;
        
        FpEI[n] = (fpz-fpy)/sq3;
        FcEI[n] = (fcz-fcy)/sq3;
        
        
    }
    
    free_double_vector(u); free_double_vector(v); free_double_vector(kv);
    
    free_double_matrix(eplus,4); free_double_matrix(ecross,4);
    
    free_double_matrix(dplus,4); free_double_matrix(dcross,4);
    
    free_double_matrix(TR,4); free_double_matrix(TI,4);
    
    free_double_matrix(kdr,4);
    
    free_double_vector(kdg);
    
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    
    free_double_vector(r12); free_double_vector(r21); free_double_vector(r31);
    free_double_vector(r13); free_double_vector(r23); free_double_vector(r32);
    free_double_vector(r10); free_double_vector(r20); free_double_vector(r30);
    
    return;
}


// The fast response code works, but is sometimes not accurate enough at early times (low frequency). Phase interpolation is good, but does not
// appear to be good enough sometimes.
void ResponseFast(int ll, double Tend, double *params, long N, double *AS, double *ES)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    double af, fr, df, DT, fac, deltaF, f, fmax, fmin, x;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess, kxm;
    double m1_SI, m2_SI, distance, tc, phic;
    int i, NF;
    double A, P;
    double px, fnew;
    
    int NFmax = 5000;
    
    double *AF, *PF, *FF, *TF;
    
    
    FILE *out;
    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FF = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
    
    SetUp(ll, params, NFmax, &NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    
    
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    
    /*
     out = fopen("inAP.dat","w");
     for (i=0; i< NF; i++)
     {
     fprintf(out,"%e %e %e %e\n", FF[i], TF[i], PF[i], AF[i]);
     }
     fclose(out);
     */
    
     /*
     out = fopen("fastAPr.dat","w");
     for (i=0; i< NF; i++)
     {
     fprintf(out,"%e %e %e %e %e\n", FF[i], AAmp[i], APhase[i], EAmp[i], EPhase[i]);
     }
     fclose(out);
      */
    
    gsl_interp_accel *PAacc = gsl_interp_accel_alloc();
    gsl_spline *PAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PAspline, FF, APhase, NF);
    
    gsl_interp_accel *PEacc = gsl_interp_accel_alloc();
    gsl_spline *PEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PEspline, FF, EPhase, NF);
    
    gsl_interp_accel *AAacc = gsl_interp_accel_alloc();
    gsl_spline *AAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(AAspline, FF, AAmp, NF);
    
    gsl_interp_accel *AEacc = gsl_interp_accel_alloc();
    gsl_spline *AEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(AEspline, FF, EAmp, NF);
    
    
    
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    
    AS[0] = 0.0;
    ES[0] = 0.0;
    AS[N/2] = 0.0;
    ES[N/2] = 0.0;
    
    // printf("%e %e\n", FF[0], FF[NF-1]);
    
   // out = fopen("fastAP_full.dat","w");
    
    for (i=1; i< N/2; i++)
    {
        f = (double)(i)/Tobs;
        AS[i] = 0.0;
        ES[i] = 0.0;
        AS[N-i] = 0.0;
        ES[N-i] = 0.0;
        
        if(f > FF[0] && f < FF[NF-1])
        {
            px = 2.0*PI*f*(Tobs-tc + dt/2.0)-2.0*phic;
            P = gsl_spline_eval (PAspline, f, PAacc);
            A = gsl_spline_eval (AAspline, f, AAacc);
            AS[i] = A*cos(P+px);
            AS[N-i] = A*sin(P+px);
            //fprintf(out,"%e %e %e ", f, A, P);
            P = gsl_spline_eval (PEspline, f, PEacc);
            A = gsl_spline_eval (AEspline, f, AEacc);
            //fprintf(out,"%e %e\n", A, P);
            ES[i] = A*cos(P+px);
            ES[N-i] = A*sin(P+px);
        }
        
    }
    
    //fclose(out);
    
    free(TF);
    free(FF);
    free(PF);
    free(AF);

    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    
    gsl_spline_free(PAspline);
    gsl_spline_free(PEspline);
    gsl_spline_free(AAspline);
    gsl_spline_free(AEspline);
    gsl_interp_accel_free(PAacc);
    gsl_interp_accel_free(PEacc);
    gsl_interp_accel_free(AAacc);
    gsl_interp_accel_free(AEacc);
    
    
}

double FofT(int ll, double *params, double tref)
{
    
    double m1, m2, m1_SI, m2_SI, chi1, chi2, tc;
    double Mtot, eta, Mc, af, fr;
    double Amp, Phase, distance;
    double fref;
    double fnew, tf;
    int i;
    
    if(ll == 0)
    {
    m1 = params[0];
    m2 = params[1];
    distance = params[6]*1.0e9*PC_SI; // distance
    }
    else
    {
    m1 = exp(params[0]);
    m2 = exp(params[1]);
    distance = exp(params[6])*1.0e9*PC_SI; // distance
    }
    
    m1_SI = m1*MSUN_SI;
    m2_SI = m2*MSUN_SI;
    chi1 = params[2];
    chi2 = params[3];
    tc = params[5];
    
    Mtot = (m1+m2)*TSUN;
    eta = m1*m2/((m1+m2)*(m1+m2));
    
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0)*TSUN;
    
    af = FinalSpin0815(eta, chi1, chi2);
    fr = fring(eta, chi1, chi2, af)/Mtot;

    
    // Here we space the frequency array to give approximately equal spacing in time
    // The dynamic frequency spacing is capped to be between 1e-6 and 1e-4 Hz
    
    if(tref > tc)
    {
        fref = 2.0*fr;
    }
    else
    {
    
    // guess at f of t
    fref = 0.9*pow( (pow(Mc,5.0/3.0)*(tc-tref)/5.0) ,-3.0/8.0)/(8.0*PI);
    
    // find the frequency at t= tref.
    i = 0;
    do
    {
        getfreq(&fnew, &tf, &Amp, &Phase, tref, fref, 0.0, PDfref, m1_SI, m2_SI, chi1, chi2, distance, tc);
        fref = fnew;
        i++;
    }while(i < 10 && fabs(tf-tref) > 1.0 && fref == fref);
    
    // nan catcher
    if(fref != fref) fref = 1.0/Tobs;
    if(fref < 0.0) fref = 1.0/Tobs;
        
    }
    
    return(fref);
    
}


void StartStop(int ll, double *params, double Tseg, double tstart, double tstop, double *fstart, double *fstop, double *frg)
{
    
    double m1, m2, m1_SI, m2_SI, chi1, chi2, tc;
    double Mtot, eta, Mc, af, fr;
    double Amp, Phase, distance;
    double fmin, fmax;
    double fnew, tf;
    int i;
    
    if(ll == 0)
    {
    m1 = params[0];
    m2 = params[1];
    distance = params[6]*1.0e9*PC_SI; // distance
    }
    else
    {
    m1 = exp(params[0]);
    m2 = exp(params[1]);
    distance = exp(params[6])*1.0e9*PC_SI; // distance
    }
    
    m1_SI = m1*MSUN_SI;
    m2_SI = m2*MSUN_SI;
    chi1 = params[2];
    chi2 = params[3];
    tc = params[5];
    
    Mtot = (m1+m2)*TSUN;
    eta = m1*m2/((m1+m2)*(m1+m2));
    
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0)*TSUN;
    
    af = FinalSpin0815(eta, chi1, chi2);
    fr = fring(eta, chi1, chi2, af)/Mtot;
    
    if(fr < 0.0) fr = 1.0/Mtot;  // should not happen ....
    
    *frg = fr;
    
    // Here we space the frequency array to give approximately equal spacing in time
    // The dynamic frequency spacing is capped to be between 1e-6 and 1e-4 Hz
    
    // guess at fmin (whre the signal starts at t=0
    fmin = 0.9*pow( (pow(Mc,5.0/3.0)*(tc-tstart)/5.0) ,-3.0/8.0)/(8.0*PI);
    if(fmin < 1.0/Tobs) fmin = 1.0/Tobs;
    
    // find the frequency at t= tstart.
    i = 0;
    do
    {
        getfreq(&fnew, &tf, &Amp, &Phase, tstart, fmin, 0.0, PDfref, m1_SI, m2_SI, chi1, chi2, distance, tc);
        if(fnew < 1.0/Tobs) fnew = 1.0/Tobs;
        fmin = fnew;
        i++;
    }while(i < 10 && fabs(tf-tstart) > 1.0 && fmin == fmin);
    
    // nan catcher
    if(fmin != fmin) fmin = 1.0/Tobs;
    if(fmin < 0.0) fmin = 1.0/Tobs;
    
    fmax = 2.0*fr;

    i = 0;
    if(tc > tstop)
    {
        fmax = 0.9*pow( (pow(Mc,5.0/3.0)*(tc-tstop)/5.0) ,-3.0/8.0)/(8.0*PI);
        if(fmax < fmin) fmax = fmin+1.0/Tobs;
        // find the frequency at t= tstop.
        do
        {
            getfreq(&fnew, &tf, &Amp, &Phase, tstop, fmax, 0.0, PDfref, m1_SI, m2_SI, chi1, chi2, distance, tc);
            if(fnew < fmin) fnew = fmin+1.0/Tobs;
            fmax = fnew;
            i++;
        }while(i < 10 && fabs(tf-tstop) > 1.0);
    }
    
    // nan catcher
    if(fmax != fmax) fmax = 2.0*fr;
    if(fmax < fmin) fmax = 2.0*fmin;
    
    *fstart = fmin;
    *fstop = fmax;
    
}

void ResponseFreq(int ll, double Tend, double *params, long N, double *AS, double *ES)
{
    
    /*   Indicies   */
    int i,j, k, n, m, a, M, nn, nmin, nmax;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double xi, t;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx;
    
    double Amp, Phase, fonfs, f, x;
    
    double Aprime, Pprime, fi, fend;
    
    double HC, HS, hp, hc;
    
    double m1, m2, chi1, chi2, phic, tc, distance, Mtot, eta, fr, af;
    
    double *ta, *xia, *FF;
    
    double Fp, Fc, kdotx, delt, fmin, fmax;
    
    double XR, XI, YR, YI, ZR, ZI;
    
    double fstart, fstop, t_tuke, Tcut;
    
    int nfmin, nfmax, nf;
    
    int NA;
    
    clock_t start, end;
    double cpu_time_used;
    
    double *FpAR, *FpAI, *FcAR, *FcAI;
    double *FpER, *FpEI, *FcER, *FcEI;
    
    FILE *out;
    
    t_tuke = 1.0e5;
    
    // Have to generate full signal to get the merger phase correct
    // There are probably ways to work around this
    StartStop(ll, params, Tobs, 0.0, Tobs, &fstart, &fstop, &fr);
 

    if(ll == 0)  // linear
    {
        m1 = params[0];    // Mass1
        m2 = params[1];    // Mass2
        distance = params[6]*1.0e9*PC_SI; // distance
    }
    else  // log
    {
        m1 = exp(params[0]);    // Mass1
        m2 = exp(params[1]);    // Mass2
        distance = exp(params[6])*1.0e9*PC_SI; // distance
    }


    chi1 = params[2];  // Spin1
    chi2 = params[3];  // Spin2
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    cosi = params[10];  // inclination
    
    Mtot = (m1+m2)*TSUN;
    eta = m1*m2/((m1+m2)*(m1+m2));
    
    Aplus = 0.5*(1.+cosi*cosi);
    Across = -cosi;
    
    AmpPhaseFDWaveform *ap = NULL;
    double m1_SI, m2_SI, deltaF;
    double fRef_in=PDfref;
    double *AF, *TF;
    int ret, flag1, flag2;
    
    m1_SI =  m1*MSUN_SI;
    m2_SI =  m2*MSUN_SI;
    
    AF = (double*)malloc(sizeof(double)* (N/2));
    TF = (double*)malloc(sizeof(double)* (N/2));

    
    nfmin = (int)(fstart*Tobs);
    nfmax = (int)(fstop*Tobs);
    if(nfmax > N/2) nfmax = N/2;
    nf = nfmax-nfmin;
    
    fmin = (double)(nfmin)/Tobs;
    fmax = (double)(nfmax)/Tobs;
    
    deltaF = 1.0/Tobs;
    
    ta = (double*)malloc(sizeof(double)* (nf));
    xia = (double*)malloc(sizeof(double)* (nf));
    FpAR = (double*)malloc(sizeof(double)* (nf));
    FcAR = (double*)malloc(sizeof(double)* (nf));
    FpER = (double*)malloc(sizeof(double)* (nf));
    FcER = (double*)malloc(sizeof(double)* (nf));
    FpAI = (double*)malloc(sizeof(double)* (nf));
    FcAI = (double*)malloc(sizeof(double)* (nf));
    FpEI = (double*)malloc(sizeof(double)* (nf));
    FcEI = (double*)malloc(sizeof(double)* (nf));
    FF = (double*)malloc(sizeof(double)* (nf));
    
    RealVector *freq;
    freq = CreateRealVector(nf);
    for (i=0; i< nf; i++) freq->data[i] = fmin+(double)(i)*deltaF;
    for (i=0; i< nf; i++) FF[i] = freq->data[i];
    
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(
                                          &ap,      /**< [out] FD waveform */
                                          freq, /**< Input: frequencies (Hz) on which to evaluate h22 FD - will be copied in the output AmpPhaseFDWaveform. Frequencies exceeding max freq covered by PhenomD will be given 0 amplitude and phase. */
                                          0.0,                  /**< Orbital phase at fRef (rad) */
                                          fRef_in,               /**< reference frequency (Hz) */
                                          m1_SI,                 /**< Mass of companion 1 (kg) */
                                          m2_SI,                 /**< Mass of companion 2 (kg) */
                                          chi1,                  /**< Aligned-spin parameter of companion 1 */
                                          chi2,                  /**< Aligned-spin parameter of companion 2 */
                                          distance               /**< Distance of source (m) */
                                          );
    
    
    // compute the Frequency time series
    timearray(params, freq, nf, TF, ap);
    
    RAantenna(params, nf, TF, FF, xia, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
    for (i=0; i< nf; i++)
    {
        AF[i] =  h22fac*ap->amp[i];
        fonfs = freq->data[i]/fstar;
        AF[i] *= (8.0*fonfs*sin(fonfs));   // conversion to fractional frequency and leading order TDI transfer function
    }
    
    nmin = nfmin;
    nmax = nfmax;
    
    // find the usable frequency range
    flag1 = 0;
    for (i=0; i< nf; i++)
    {
        if(TF[i] > 0.0 && flag1 == 0)
        {
            flag1 = 1;
            nmin = i+nfmin;
        }
    }
    

    for(n=0; n< N; n++)
    {
        AS[n] = 0.0;
        ES[n] = 0.0;
    }
        
    Tcut = Tend+t_rise/1.5;
    
    
    /*   Main Loop   */
    
    nn = 0;
    for(n=nmin; n< nmax; n++)
    {
        // Barycenter time and frequency
        
        // The frequency and time arrays start at nfmin
        m = n-nfmin;
        
        if(m > -1 && m < nf)
        {
            t = TF[m];
            f = FF[m];
            xi = xia[m];
            
            x = 1.0;
            
            // Tukey filter to match what is done in time domain
            if(t < t_tuke) x = 0.5*(1.0+cos(PI*(t/t_tuke-1.0)));
            
            // taper
            if(t > Tcut-t_rise && t < Tcut)
            {
            x = 0.5*(1.0-cos(PI*(t-Tcut)/t_rise));
            }
            if(t > Tcut) x = 0.0;
            
            kdotx = t-xi;
            
            Amp = x*AF[m];
            Phase = ap->phase[m]+2.0*phic;
            
            HC = Amp*cos(2.0*PI*f*(Tobs-tc+dt/2.0-kdotx)-Phase);
            HS = Amp*sin(2.0*PI*f*(Tobs-tc+dt/2.0-kdotx)-Phase);
            

            AS[n] = FpAR[m]*Aplus*HC - FpAI[m]*Aplus*HS - FcAR[m]*Across*HS - FcAI[m]*Across*HC;
            AS[N-n] = FpAI[m]*Aplus*HC + FpAR[m]*Aplus*HS - FcAI[m]*Across*HS + FcAR[m]*Across*HC;
            
            ES[n] = FpER[m]*Aplus*HC - FpEI[m]*Aplus*HS - FcER[m]*Across*HS - FcEI[m]*Across*HC;
            ES[N-n] = FpEI[m]*Aplus*HC + FpER[m]*Aplus*HS - FcEI[m]*Across*HS + FcER[m]*Across*HC;
           
            
        }
        
        
    }
    
    
    /*   Deallocate Arrays   */
        
        DestroyAmpPhaseFDWaveform(ap);
        DestroyRealVector(freq);

    
    free(AF);
    free(TF);
    
    free(ta);
    free(xia);
    free(FcAR);
    free(FpAR);
    free(FcER);
    free(FpER);
    free(FcAI);
    free(FpAI);
    free(FcEI);
    free(FpEI);
    free(FF);
    
    return;
}



double nh(int ll, double Tend, double *params, double *ATR, double *ATI, double *ETR, double *ETI, int NF, double *FF, double *AP, double *EP, double *FN, double *ASD)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    int i, n, nx, flag;
    double logL;
    double *AF, *PF, *TF;
    double x, y, u, v, z, px, phic, tc, f, xx, yy, t;
    double kxm, ca, ce, df;
    int nn;
    double yold, zold;
    double *dA, *dE;
    double *AR, *ER;
    double *AI, *EI;
    double tau, iold;
    double *cA, *cE;
    double *sA, *sE;
    int NM;
    double *Fs, *As, *Es;
    double *Ac, *Ec;
    double fr;
    
    FILE * out;

    Ac = (double*)malloc(sizeof(double)* (MF));
    Ec = (double*)malloc(sizeof(double)* (MF));
    As = (double*)malloc(sizeof(double)* (MF));
    Es = (double*)malloc(sizeof(double)* (MF));

    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    
    AR = (double*)malloc(sizeof(double)* (NF));
    AI = (double*)malloc(sizeof(double)* (NF));
    ER = (double*)malloc(sizeof(double)* (NF));
    EI = (double*)malloc(sizeof(double)* (NF));
    
    cA = (double*)malloc(sizeof(double)* (NF));
    cE = (double*)malloc(sizeof(double)* (NF));
    sA = (double*)malloc(sizeof(double)* (NF));
    sE = (double*)malloc(sizeof(double)* (NF));
    
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    flag = 0;
    
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    for(n=0; n< NF; n++)
    {
        px = TPI*FF[n]*(Tobs-tc+ dt/2.0)-2.0*phic;
        APhase[n] += px;
        EPhase[n] += px;
    }
    
    nx = NF/4;
    
    
    for(n=0; n< NF; n++)
    {
        cA[n] = AAmp[n]*cos(APhase[n]-AP[n])/ASD[n];
        cE[n] = EAmp[n]*cos(EPhase[n]-EP[n])/ASD[n];
        sA[n] = AAmp[n]*sin(APhase[n]-AP[n])/ASD[n];
        sE[n] = EAmp[n]*sin(EPhase[n]-EP[n])/ASD[n];
    }
    
    
    
    // count the number of zero crossings
    nn = 0;
    x = cA[0];
    y = cE[0];
    for(n=1; n< NF; n++)
    {
        if(x*cA[n] < 0.0 || y*cE[n] < 0) nn++;
        x = cA[n];
        y = cE[n];
    }
    
    z = -1.0e6;

    
    if(nn < 40)
    {
    
    gsl_interp_accel *Aacc = gsl_interp_accel_alloc();
    gsl_spline *Aspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_interp_accel *Eacc = gsl_interp_accel_alloc();
    gsl_spline *Espline = gsl_spline_alloc (gsl_interp_cspline, NF);
    

    gsl_spline_init(Aspline, FF, cA, NF);
    gsl_spline_init(Espline, FF, cE, NF);

    
    for(i=0; i< MF; i++)
    {
        f = FN[i];
        Ac[i] = 0.0;
        Ec[i] = 0.0;
        
       if(f >= FF[0] && f <= FF[NF-1])
        {
         Ac[i] = gsl_spline_eval (Aspline, f, Aacc);
         Ec[i] = gsl_spline_eval (Espline, f, Eacc);
        }
    
    }
        
        gsl_spline_init(Aspline, FF, sA, NF);
        gsl_spline_init(Espline, FF, sE, NF);
        
        for(i=0; i< MF; i++)
        {
            f = FN[i];
            As[i] = 0.0;
            Es[i] = 0.0;
            
            if(f >= FF[0] && f <= FF[NF-1])
            {
                As[i] = gsl_spline_eval (Aspline, f, Aacc);
                Es[i] = gsl_spline_eval (Espline, f, Eacc);
            }
            
            
        }
    
        
    gsl_spline_free(Aspline);
    gsl_interp_accel_free(Aacc);
    gsl_spline_free(Espline);
    gsl_interp_accel_free(Eacc);
        
    // The tapering has very little impact. Spectral leakage not a problem
        
     /*
      fr = 1.0e-4;
      out = fopen("cossin.dat","w");
     for(i=0; i< MF; i++)
       {
          f = FN[i];
           if(f >= FF[0] && f <= FF[0]+fr)
           {
               x = 0.5*(1.0-cos(PI*(f-FF[0])/fr));
               As[i] *= x;
               Es[i] *= x;
               Ac[i] *= x;
               Ec[i] *= x;
           }
           
           fprintf(out,"%e %e %e %e %e\n", f, Ac[i], Ec[i], As[i], Es[i]);
           
       }
    fclose(out);
      */
    
    gsl_fft_real_radix2_transform(As, 1, MF);
    gsl_fft_real_radix2_transform(Es, 1, MF);
    gsl_fft_real_radix2_transform(Ac, 1, MF);
    gsl_fft_real_radix2_transform(Ec, 1, MF);
    
    z = 0.0;
        
    //out = fopen("hetint.dat","w");
        
    for(i=0; i< MF/2; i++)
    {
        if(i > 0)
        {
            z += Ac[i]*ATR[i]+Ac[MF-i]*ATR[MF-i]+Ec[i]*ETR[i]+Ec[MF-i]*ETR[MF-i]
                +(As[MF-i]*ATI[MF-i]+Es[MF-i]*ETI[MF-i]+As[i]*ATI[i]+Es[i]*ETI[i]);
        }
        else
        {
            z += 0.5*(Ac[i]*ATR[i]+Ec[i]*ETR[i]+(As[i]*ATI[i]+Es[i]*ETI[i]));
        }
        
       // fprintf(out, "%e %e\n", (double)(i)*dt, z*8.0/(Tobs*(double)(MF)));
        
    }
        
        z *= 8.0/(Tobs*(double)(MF));
        
        // fclose(out);
        
        
    }
    
   
    /*
    out = fopen("het.dat","w");
    
    i = 0;
    fprintf(out, "%e %e %e %e %e %e %e %e %e\n", (double)(i)*dt, Ac[i], 0.0, As[i], 0.0, ATR[i], 0.0, ATI[i], 0.0);
    
    for(i=1; i< MF/2; i++)
    {
        fprintf(out, "%e %e %e %e %e %e %e %e %e\n", (double)(i)*dt, Ac[i],  Ac[MF-i], As[i], As[MF-i], ATR[i], ATR[MF-i], ATI[i], ATI[MF-i]);
    }
    
    fclose(out);
     
    */
    
    free(cA);
    free(cE);
    free(sA);
    free(sE);
    
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(TF);
    free(AF);
    free(PF);
    free(AR);
    free(AI);
    free(ER);
    free(EI);
    
    free(Ac);
    free(Ec);
    free(As);
    free(Es);
    
    return(z);
    
    
}

double ndh(int ll, double Tend, double *params, double *ATR, double *ATI, double *ETR, double *ETI, int NF, double *FF, double *AA, double *EE, double *AP, double *EP, double *FN, double *ASD)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    int i, n, nx, flag;
    double logL;
    double *AF, *PF, *TF;
    double x, y, u, v, z, px, phic, tc, f, t;
    double kxm, ca, ce, df;
    int nn;
    double yold, zold;
    double *dA, *dE;
    double *AR, *ER;
    double *AI, *EI;
    double tau, iold;
    double *cA, *cE, *sA, *sE;
    int NM;
    double *Fs, *As, *Es, *Ac, *Ec;
    
    FILE * out;
    
    Ac = (double*)malloc(sizeof(double)* (MF));
    Ec = (double*)malloc(sizeof(double)* (MF));
    As = (double*)malloc(sizeof(double)* (MF));
    Es = (double*)malloc(sizeof(double)* (MF));
    
    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    
    AR = (double*)malloc(sizeof(double)* (NF));
    AI = (double*)malloc(sizeof(double)* (NF));
    ER = (double*)malloc(sizeof(double)* (NF));
    EI = (double*)malloc(sizeof(double)* (NF));
    
    cA = (double*)malloc(sizeof(double)* (NF));
    cE = (double*)malloc(sizeof(double)* (NF));
    sA = (double*)malloc(sizeof(double)* (NF));
    sE = (double*)malloc(sizeof(double)* (NF));
    
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    flag = 0;
    
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    for(n=0; n< NF; n++)
    {
        px = TPI*FF[n]*(Tobs-tc+ dt/2.0)-2.0*phic;
        APhase[n] += px;
        EPhase[n] += px;
    }
    
    nx = NF/4;
    

    for(n=0; n< NF; n++)
    {
        cA[n] = (AA[n]-AAmp[n]*cos(APhase[n]-AP[n]))/ASD[n];
        cE[n] = (EE[n]-EAmp[n]*cos(EPhase[n]-EP[n]))/ASD[n];
        sA[n] = -AAmp[n]*sin(APhase[n]-AP[n])/ASD[n];
        sE[n] = -EAmp[n]*sin(EPhase[n]-EP[n])/ASD[n];
    }
    
    // count the number of zero crossings
    nn = 0;
    x = cA[0];
    y = cE[0];
    for(n=1; n< NF; n++)
    {
        if(x*cA[n] < 0.0 || y*cE[n] < 0) nn++;
        x = cA[n];
        y = cE[n];
    }
    
    z = 1.0e6;
    
    
    if(nn < 40)
    {
        
        gsl_interp_accel *Aacc = gsl_interp_accel_alloc();
        gsl_spline *Aspline = gsl_spline_alloc (gsl_interp_cspline, NF);
        gsl_interp_accel *Eacc = gsl_interp_accel_alloc();
        gsl_spline *Espline = gsl_spline_alloc (gsl_interp_cspline, NF);
        
        
        gsl_spline_init(Aspline, FF, cA, NF);
        gsl_spline_init(Espline, FF, cE, NF);
        
        
        for(i=0; i< MF; i++)
        {
            f = FN[i];
            Ac[i] = 0.0;
            Ec[i] = 0.0;
            
            if(f >= FF[0] && f <= FF[NF-1])
            {
                Ac[i] = gsl_spline_eval (Aspline, f, Aacc);
                Ec[i] = gsl_spline_eval (Espline, f, Eacc);
            }
            
            
        }
        
        gsl_spline_init(Aspline, FF, sA, NF);
        gsl_spline_init(Espline, FF, sE, NF);
        
        
        for(i=0; i< MF; i++)
        {
            f = FN[i];
            As[i] = 0.0;
            Es[i] = 0.0;
            
            if(f >= FF[0] && f <= FF[NF-1])
            {
                As[i] = gsl_spline_eval (Aspline, f, Aacc);
                Es[i] = gsl_spline_eval (Espline, f, Eacc);
            }
        }

        
        
        gsl_spline_free(Aspline);
        gsl_interp_accel_free(Aacc);
        gsl_spline_free(Espline);
        gsl_interp_accel_free(Eacc);
        
        gsl_fft_real_radix2_transform(Ac, 1, MF);
        gsl_fft_real_radix2_transform(Ec, 1, MF);
        gsl_fft_real_radix2_transform(As, 1, MF);
        gsl_fft_real_radix2_transform(Es, 1, MF);
        
        z = 0.0;
        
        for(i=0; i< MF/2; i++)
        {
            if(i > 0)
            {
                z += Ac[i]*ATR[i]+Ac[MF-i]*ATR[MF-i]+Ec[i]*ETR[i]+Ec[MF-i]*ETR[MF-i]
                +(As[MF-i]*ATI[MF-i]+Es[MF-i]*ETI[MF-i]+As[i]*ATI[i]+Es[i]*ETI[i]);
            }
            else
            {
                z += 0.5*(Ac[i]*ATR[i]+Ec[i]*ETR[i]+(As[i]*ATI[i]+Es[i]*ETI[i]));
            }
            
        }
        
        z *= 8.0/(Tobs*(double)(MF));
        
    }
    

    
    
    
    //printf("%d %e\n", nn, z);
    
    free(cA);
    free(cE);
    free(sA);
    free(sE);
    
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(TF);
    free(AF);
    free(PF);
    free(AR);
    free(AI);
    free(ER);
    free(EI);
    
    free(As);
    free(Es);
    free(Ac);
    free(Ec);
    
    return(z);
    
    
}

void pbt_shift(double *corr, double *corrf, double *data1, double *data2, double *Sn, int n)
{
    int nb2, i, l, k, j;
    int imax, imin;
    
    nb2 = n/2;
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        
        corr[l]    = (data1[l]*data2[l] + data1[k]*data2[k])/Sn[i];
        corr[k]    = (data1[k]*data2[l] - data1[l]*data2[k])/Sn[i];
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
    }
    
    corr[0] = 0.0;
    corrf[0] = 0.0;
    corr[nb2] = 0.0;
    corrf[nb2] = 0.0;
    
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, n);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, n);
    
    
}


double log_likelihood_max_dual(int ll, double *A, double *E, double *params, double *SN, int N)
{
    int i, j, k, NMAX;
    int ii, jj, Nend;
    double sum;
    double HH, LL;
    double logL;
    double logLfs;
    double fmax;
    double HA, HE, LD, x;
    double normA, normE, deltH, pshiftA, pshiftE;
    double *AS, *ES;
    double *AC, *AF;
    double *EC, *EF;
    double T1, T2, a, b, c, dtx;
   
    //Chi^2 = (d-h|d-h)/Sn
    
    AS = double_vector(N); ES = double_vector(N);
    AC=double_vector(N);  AF=double_vector(N);
    EC=double_vector(N);  EF=double_vector(N);
    
    ResponseFreq(ll, -1.0, params, N, AS, ES);
    
    HH = fourier_nwip(AS, AS, SN, N)+fourier_nwip(ES, ES, SN, N);
    
    pbt_shift(AC, AF, A, AS, SN, N);
    pbt_shift(EC, EF, E, ES, SN, N);
    
    // only allow time shifts up to 500 seconds
    Nend = (int)(1000.0/dt);
    
    for(i = 0; i < Nend/2; i++) AS[i+N/2] = sqrt(AC[i]*AC[i]+AF[i]*AF[i]) + sqrt(EC[i]*EC[i]+EF[i]*EF[i]);
    for(i = -Nend/2; i < 0; i++) AS[i+N/2] = sqrt(AC[N+i]*AC[N+i]+AF[N+i]*AF[N+i]) + sqrt(EC[N+i]*EC[N+i]+EF[N+i]*EF[N+i]);
    
    x = 0;
    for (i = -Nend/2; i < Nend/2; ++i)
    {
        if(AS[i+N/2] > x)
        {
            x = AS[i+N/2];
            k = i;
        }
        
        printf("%d %e\n", i, AS[i+N/2]);
    }
    
    
    T1 = sqrt(AS[k+1+N/2]);
    T2 = sqrt(AS[k-1+N/2]);
    
    c = sqrt(AS[k+N/2]);
    b = (T1-T2)/(2.0*dt);
    a = (T1-b*dt-c)/(dt*dt);
    
    dtx = -b/(2.0*a);
    
    HA = 2.0*(double)(N)*(AS[k+N/2]);
    
    printf("Max = %f  delay = %f %f\n", x, (double)(k)*dt, dtx);
    
    deltH = (double)(k)*dt+dtx;
    
    if(k < 0) k += N;
    
    // Inverse FFTs in fourier_nwip are un-normlaized
    HH /= Tobs;
    HA /= Tobs;

    pshiftA = 0.5*(atan2(AF[k],AC[k])+atan2(EF[k],EC[k]));
    
    printf("%e %f\n", HH, sqrt(HH));
    
    normA = HA/HH;
  
    if(ll==0)
    {
     params[6] /= normA;
    }
    else
    {
     params[6] -= log(normA);
    }
    params[5] += deltH;
    params[4] -= pshiftA/2.0;
    
    printf("\n dphase %f dtime %f scale %f\n", pshiftA, deltH, normA);

    logL = (HA*HA)/(2.0*HH);
    
    free_double_vector(AC);  free_double_vector(AF);
    free_double_vector(EC);  free_double_vector(EF);

    free_double_vector(AS);
    free_double_vector(ES);
    
    return logL;
}

double SNRFast(double *params)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    double af, fr, df, DT, fac, deltaF, f, fmax, fmin, x, y, t, t10;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess, kxm;
    double m1_SI, m2_SI, distance, tc, phic;
    int i, NF, NW, flag;
    double A, P;
    double SNR;
    double px, fnew;
    double *Aint, *Eint, *SAE;
    double *SNRSQA, *SNRSQE;
    
    FILE *out;
    
    int NFmax = 5000;
    
    double *AF, *PF, *FF, *TF;
    
    //clock_t start, end;
    //double cpu_time_used;
    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FF = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
                         
    SetUp(0, params, NFmax, &NF, FF, TF, PF, AF);
                         
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
                         
    Extrinsic(params, Tobs, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    
    Aint = (double*)malloc(sizeof(double)* (NF));
    Eint = (double*)malloc(sizeof(double)* (NF));
    SAE = (double*)malloc(sizeof(double)* (NF));
    

   // out = fopen("intg.dat","w");
    for (i=0; i< NF; i++)
    {
        instrument_noise(FF[i], &SAE[i], &x);
        Aint[i] = 4.0*AAmp[i]*AAmp[i]/SAE[i];
        Eint[i] = 4.0*EAmp[i]*EAmp[i]/SAE[i];
       // fprintf(out,"%e %e %e %e %e\n", TF[i], FF[i], Aint[i], Eint[i], SAE[i]);
    }
    //fclose(out);
    
    gsl_interp_accel *PAacc = gsl_interp_accel_alloc();
    gsl_spline *PAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PAspline, FF, APhase, NF);
    
    gsl_interp_accel *PEacc = gsl_interp_accel_alloc();
    gsl_spline *PEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PEspline, FF, EPhase, NF);
    
    gsl_interp_accel *AAacc = gsl_interp_accel_alloc();
    gsl_spline *AAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(AAspline, FF, AAmp, NF);
    
    gsl_interp_accel *AEacc = gsl_interp_accel_alloc();
    gsl_spline *AEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(AEspline, FF, EAmp, NF);
    
    gsl_interp_accel *IAacc = gsl_interp_accel_alloc();
    gsl_spline *IAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(IAspline, FF, Aint, NF);
    
    gsl_interp_accel *IEacc = gsl_interp_accel_alloc();
    gsl_spline *IEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(IEspline, FF, Eint, NF);
    
    SNRSQA = (double*)malloc(sizeof(double)* (NF));
    SNRSQE = (double*)malloc(sizeof(double)* (NF));
    
    SNRSQA[0] = 0.0;
    SNRSQE[0] = 0.0;
    for (i=1; i< NF; i++)
    {
     SNRSQA[i] = SNRSQA[i-1]+gsl_spline_eval_integ(IAspline, FF[i-1], FF[i], IAacc);
     SNRSQE[i] = SNRSQE[i-1]+gsl_spline_eval_integ(IEspline, FF[i-1], FF[i], IEacc);
    }
    
    double SA, SE, SAT, SET, tmax;
    
    SA = SNRSQA[NF-1];
    SE = SNRSQE[NF-1];
    
    SNR = sqrt(SA+SE);
    
   // printf("%f %f %e\n", SA, SE, SNRSQA[NF-1]+SNRSQE[NF-1]);
    
    

    gsl_interp_accel *SAacc = gsl_interp_accel_alloc();
    gsl_spline *SAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(SAspline, TF, SNRSQA, NF);
    
    gsl_interp_accel *SEacc = gsl_interp_accel_alloc();
    gsl_spline *SEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(SEspline, TF, SNRSQE, NF);
    
    
    NW = (int)(Tobs/1000.0);
    tmax = TF[NF-1];
    out = fopen("SNRofT.dat","w");
    
    flag = 0;
    
    for (i=0; i< NW; i++)
    {
        t = (double)(i)*1000.0;
        SAT = SA;  // make sure we get to the final value
        SET = SE;
        if(t < tmax)
        {
            SAT = gsl_spline_eval (SAspline, t, SAacc);
            SET = gsl_spline_eval (SEspline, t, SEacc);
        }
        
        x = sqrt(SAT+SET);
        
        if(x > 10.0 && flag == 0)
        {
            flag = 1;
            y = t;
        }
        
        fprintf(out,"%e %e %e %e\n", t, x, sqrt(SAT), sqrt(SET));
    }
    fclose(out);
    
    flag = 0;
    
    for (i=-1000; i< 1000; i++)
    {
        t = y+(double)(i);
        
        SAT = SA;  // make sure we get to the final value
        SET = SE;
        if(t < tmax)
        {
            SAT = gsl_spline_eval (SAspline, t, SAacc);
            SET = gsl_spline_eval (SEspline, t, SEacc);
        }
        
        x = sqrt(SAT+SET);
        
        if(x > 10.0 && flag == 0)
        {
            flag = 1;
            t10 = t;
        }
        
    }
 
    
    printf("SNR=10 at t = %e, tc-t = %e\n", t10, params[5]-t10);
     
    
    
    
    free(SNRSQA);
    free(SNRSQE);
    free(SAE);
    free(Aint);
    free(Eint);
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(FF);
    free(TF);
    free(AF);
    free(PF);
    
    gsl_spline_free(IAspline);
    gsl_spline_free(IEspline);
    gsl_spline_free(PAspline);
    gsl_spline_free(PEspline);
    gsl_spline_free(AAspline);
    gsl_spline_free(AEspline);
    gsl_interp_accel_free(IAacc);
    gsl_interp_accel_free(IEacc);
    gsl_interp_accel_free(PAacc);
    gsl_interp_accel_free(PEacc);
    gsl_interp_accel_free(AAacc);
    gsl_interp_accel_free(AEacc);
    
    return(SNR);
    
}

double Likelihood(int ll, double Tend, double *params, long N, double *AD, double *ED, double *SN)
{
    double *AS, *ES;
    double HH, HD;
    double logL;
    
    AS = (double*)malloc(sizeof(double)* (N));
    ES = (double*)malloc(sizeof(double)* (N));
    
    ResponseFreq(ll, Tend, params, N, AS, ES);
    
    HH = (fourier_nwip(AS, AS, SN, NFFT)+fourier_nwip(ES, ES, SN, NFFT))/Tobs;
    HD = (fourier_nwip(AS, AD, SN, NFFT)+fourier_nwip(ES, ED, SN, NFFT))/Tobs;
    
    logL = HD-0.5*HH;
    
    free(AS);
    free(ES);
    
    return(logL);
    
}

double LikelihoodFstat(int ll, double Tend, double *params, double tm, int NF, double *FF, double *AA, double *EA, double *AP, double *EP, double *SN)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    int i, n, flag;
    double logL;
    double *AF, *PF, *TF, *IF;
    double x, y, z, px, phic, tc, tx;
    double *pnew;

    pnew = (double*)malloc(sizeof(double)* (6));
    
    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, NF, FF, TF, PF, AF);
    
    // time maximize
    if(tm > 0.0)
    {
     tx = Tmerger(params,params[5]);
        params[5] += (tm-tx);
    }
    
    FstatRA(ll, Tend, params, pnew, NF, FF, TF, PF, AF, AA, EA, AP, EP, SN);
    
    logL = pnew[0];
    
    // [0] ln Mass1  [1] ln Mass2  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln distance
    // [7] cos EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] cos inclination
    
    params[10] = pnew[1];
    params[6] += log(pnew[2]);
    params[9] = pnew[3];
    params[4] = pnew[4];

    free(pnew);
    free(TF);
    free(AF);
    free(PF);
    
    return(logL);
    
}

double LikelihoodDeltaMaxT(int ll, double Tend, double *params, int NF, double *FF, double *AA, double *EA, double *AP, double *EP, double *SN)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    int i, n, flag;
    double logL;
    double *AF, *PF, *TF, *IC, *IS, *IA;
    double *ICF, *ISF, *ICFF, *ISFF;
    double x, y, z, px, phic, tc, t;
    double kxm, ca, ce;
    double C, S, A, dphase;
    double CF, CFF, SF, SFF;
    int nn;
    double yold, zold;
    double T1, T2, T3, T4, T5, a, b, c;
    double deltT, dtx;
    
    FILE *out;

    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    IC = (double*)malloc(sizeof(double)* (NF));
    IS = (double*)malloc(sizeof(double)* (NF));
    ICF = (double*)malloc(sizeof(double)* (NF));
    ISF = (double*)malloc(sizeof(double)* (NF));
    ICFF = (double*)malloc(sizeof(double)* (NF));
    ISFF = (double*)malloc(sizeof(double)* (NF));
    IA = (double*)malloc(sizeof(double)* (NF));
    
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    for(n=0; n< NF; n++)
    {
      px = TPI*FF[n]*(Tobs-tc +dt/2.0)-2.0*phic;
       // printf("%d %e\n", n, APhase[n]);
      APhase[n] += px;
      EPhase[n] += px;
    }
    
    deltT = 10.0;
    
    yold = 1.0;
    zold = 1.0;
    nn = 0;
    flag = 0;
    n = 0;
    IA[n] = 2.0*(AAmp[n]*AAmp[n]+EAmp[n]*EAmp[n])/SN[n];
    IC[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n])+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]))/SN[n];
    ICF[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n]+TPI*FF[n]*deltT)+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]+TPI*FF[n]*deltT))/SN[n];
    ICFF[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n]-TPI*FF[n]*deltT)+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]-TPI*FF[n]*deltT))/SN[n];
    IS[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n])+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]))/SN[n];
    ISF[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n]+TPI*FF[n]*deltT)+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]+TPI*FF[n]*deltT))/SN[n];
    ISFF[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n]-TPI*FF[n]*deltT)+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]-TPI*FF[n]*deltT))/SN[n];
    
    for(n=1; n< NF; n++)
    {
        ca = cos(AP[n]-APhase[n]);
        ce = cos(EP[n]-EPhase[n]);
        
        IA[n] = 2.0*(AAmp[n]*AAmp[n]+EAmp[n]*EAmp[n])/SN[n];
        IC[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n])+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]))/SN[n];
        ICF[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n]+TPI*FF[n]*deltT)+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]+TPI*FF[n]*deltT))/SN[n];
        ICFF[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n]-TPI*FF[n]*deltT)+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]-TPI*FF[n]*deltT))/SN[n];
        IS[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n])+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]))/SN[n];
        ISF[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n]+TPI*FF[n]*deltT)+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]+TPI*FF[n]*deltT))/SN[n];
        ISFF[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n]-TPI*FF[n]*deltT)+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]-TPI*FF[n]*deltT))/SN[n];
        

       // fprintf(out,"%e %e %e %e %e\n", FF[n], (AA[n]*AA[n])/SN[n], (AAmp[n]*AAmp[n])/SN[n], (AA[n]*AAmp[n]*ca)/SN[n], ca);

        x = (FF[n]-FF[n-1])*IA[n];  // contribution to the integrand
        if(x > 10.0)  // only care about phase differences in regions that contribute signifcantly
        {
            y = ca;
            z = ce;
            if(y*yold < 0.0)
            {
                nn++;
                yold = y;
            }
            if(z*zold < 0.0)
            {
                nn++;
                zold = z;
            }
        }
    }
    
   // fclose(out);
    
    if(nn > 20) flag = 1;
    
    logL = 0.0;
    
    if(flag == 0)
    {

        
    gsl_interp_accel *Iacc = gsl_interp_accel_alloc();
    gsl_spline *Ispline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(Ispline, FF, IC, NF);
    C =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, ICF, NF);
    CF =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, ICFF, NF);
    CFF =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, IS, NF);
    S =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, ISF, NF);
    SF =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, ISFF, NF);
    SFF =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, IA, NF);
    A =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);

        
        T1 = sqrt(CF*CF+SF*SF);
        T2 = sqrt(CFF*CFF+SFF*SFF);
        
        c = sqrt(C*C+S*S);
        b = (T1-T2)/(2.0*deltT);
        a = (T1-b*deltT-c)/(deltT*deltT);
        
        dtx = -b/(2.0*a);
        
    /*
    out = fopen("quad.dat","w");
        
    fprintf(out,"%e %e\n", -deltT, sqrt(CFF*CFF+SFF*SFF));
    fprintf(out,"%e %e\n", 0.0, sqrt(C*C+S*S));
    fprintf(out,"%e %e\n", deltT, sqrt(CF*CF+SF*SF));
        
    fclose(out);
    
    out = fopen("quadT.dat","w");
    for(n=-50; n <=50 ; n++)
    {
      t = (double)n*deltT/25.0;
      fprintf(out,"%e %e\n", t, a*t*t+b*t+c);
    }
   fclose(out); */
        
        for(n=0; n< NF; n++)
        {
            IC[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n]+TPI*FF[n]*dtx)+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]+TPI*FF[n]*dtx))/SN[n];
            IS[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n]+TPI*FF[n]*dtx)+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]+TPI*FF[n]*dtx))/SN[n];
        }
        
        gsl_spline_init(Ispline, FF, IC, NF);
        C =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
        gsl_spline_init(Ispline, FF, IS, NF);
        S =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
        gsl_spline_free(Ispline);
        gsl_interp_accel_free(Iacc);
        
        
    dphase = atan2(S,C);
        
      //  printf("%f %f\n", dtx, dphase/2.0);
        
        
    x = sqrt(C*C+S*S);
        
    y = x/(2.0*A);
        
    logL  = -y*y*A + y*x;
        
    if(logL > 0.0)
    {
    //printf("%d %e %f %f\n", nn, logL, dphase, y);
        
    params[4] -= dphase/2.0;
    params[5] += dtx;
      
      if(ll = 1)
      {
       params[6] -= log(y);
      }
      else
      {
      params[6] /= y;
      }
        
      }
    else
    {
     logL = 0.0;
    }
        
    }
    

    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(TF);
    free(AF);
    free(PF);
    free(IC);
    free(IS);
    free(IA);
    free(ICF);
    free(ISF);
    free(ICFF);
    free(ISFF);
    
    return(logL);

}



// Maximizes over phase and distance
double LikelihoodDeltaMax(int ll, double Tend, double *params, int NF, double *FF, double *AA, double *EA, double *AP, double *EP, double *SN)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    int i, n, flag;
    double logL;
    double *AF, *PF, *TF, *IC, *IS, *IA;
    double x, y, z, px, phic, tc;
    double kxm, ca, ce;
    double C, S, A, dphase;
    int nn;
    double yold, zold;
    
    FILE *out;

    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    IC = (double*)malloc(sizeof(double)* (NF));
    IS = (double*)malloc(sizeof(double)* (NF));
    IA = (double*)malloc(sizeof(double)* (NF));
    
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    /*
    out = fopen("fastAPl.dat","w");
    for (i=0; i< NF; i++)
    {
        fprintf(out,"%e %e %e %e %e\n", FF[i], AAmp[i], APhase[i], EAmp[i], EPhase[i]);
    }
    fclose(out);
    */
    
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    for(n=0; n< NF; n++)
    {
      px = TPI*FF[n]*(Tobs-tc +dt/2.0)-2.0*phic;
       // printf("%d %e\n", n, APhase[n]);
      APhase[n] += px;
      EPhase[n] += px;
    }
    
    //out = fopen("fast.dat","w");
    
    yold = 1.0;
    zold = 1.0;
    nn = 0;
    flag = 0;
    n = 0;
    IC[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n])+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]))/SN[n];
    IS[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n])+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]))/SN[n];
    IA[n] = 2.0*(AAmp[n]*AAmp[n]+EAmp[n]*EAmp[n])/SN[n];
    for(n=1; n< NF; n++)
    {
        ca = cos(AP[n]-APhase[n]);
        ce = cos(EP[n]-EPhase[n]);
        
        IC[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n])+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]))/SN[n];
        IS[n] = 4.0*(AA[n]*AAmp[n]*sin(AP[n]-APhase[n])+EA[n]*EAmp[n]*sin(EP[n]-EPhase[n]))/SN[n];
        IA[n] = 2.0*(AAmp[n]*AAmp[n]+EAmp[n]*EAmp[n])/SN[n];

       // fprintf(out,"%e %e %e %e %e\n", FF[n], (AA[n]*AA[n])/SN[n], (AAmp[n]*AAmp[n])/SN[n], (AA[n]*AAmp[n]*ca)/SN[n], ca);

        x = (FF[n]-FF[n-1])*IA[n];  // contribution to the integrand
        if(x > 10.0)  // only care about phase differences in regions that contribute signifcantly
        {
            y = ca;
            z = ce;
            if(y*yold < 0.0)
            {
                nn++;
                yold = y;
            }
            if(z*zold < 0.0)
            {
                nn++;
                zold = z;
            }
        }
    }
    
   // fclose(out);
    
    if(nn > 20) flag = 1;
    
    logL = 0.0;
    
    if(flag == 0)
    {

        
    gsl_interp_accel *Iacc = gsl_interp_accel_alloc();
    gsl_spline *Ispline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(Ispline, FF, IC, NF);
    C =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, IS, NF);
    S =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, IA, NF);
    A =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_free(Ispline);
    gsl_interp_accel_free(Iacc);
        
    dphase = atan2(S,C);
        
    x = sqrt(C*C+S*S);
        
    y = x/(2.0*A);
        
    logL  = -y*y*A + y*x;
        
    if(logL > 0.0)
    {
    //printf("%d %e %f %f\n", nn, logL, dphase, y);
        
    params[4] -= dphase/2.0;
      
      if(ll = 1)
      {
       params[6] -= log(y);
      }
      else
      {
      params[6] /= y;
      }
        
      }
    else
    {
     logL = 0.0;
    }
        
    }


    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(TF);
    free(AF);
    free(PF);
    free(IC);
    free(IS);
    free(IA);
    
    return(logL);

}



double LikelihoodDelta(int ll, double Tend, double *params, int NF, double *FF, double *AA, double *EA, double *AP, double *EP, double *SN)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    int i, n, flag;
    double logL;
    double *AF, *PF, *TF, *IF1, *IF2;
    double x, y, z, px, phic, tc, zz;
    double kxm, ca, ce, Lc;
    int nn, nt, ns;
    double yold, zold;
    
    FILE *out;

    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    IF1 = (double*)malloc(sizeof(double)* (NF));
    IF2 = (double*)malloc(sizeof(double)* (NF));
    
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    /*
    out = fopen("fastAPl.dat","w");
    for (i=0; i< NF; i++)
    {
        fprintf(out,"%e %e %e %e %e\n", FF[i], AAmp[i], APhase[i], EAmp[i], EPhase[i]);
    }
    fclose(out);
    */
    
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    for(n=0; n< NF; n++)
    {
      px = TPI*FF[n]*(Tobs-tc +dt/2.0)-2.0*phic;
       // printf("%d %e\n", n, APhase[n]);
      APhase[n] += px;
      EPhase[n] += px;
    }
    
    //out = fopen("fast.dat","w");
    
    
    nn = 0;
    nt = 0;
    flag = 0;
    n = 0;
    IF1[n] = 4.0*(AA[n]*AAmp[n]*cos(AP[n]-APhase[n])+EA[n]*EAmp[n]*cos(EP[n]-EPhase[n]))/SN[n];
    IF2[n] = -2.0*(AAmp[n]*AAmp[n]+EAmp[n]*EAmp[n])/SN[n];
    yold = cos(AP[n]-APhase[n]);
    zold = cos(AP[n]-APhase[n]);
    
    zz = 0.0;
    
    for(n=1; n< NF; n++)
    {
        ca = cos(AP[n]-APhase[n]);
        ce = cos(EP[n]-EPhase[n]);
        
        IF1[n] = 4.0*(AA[n]*AAmp[n]*ca+EA[n]*EAmp[n]*ce)/SN[n];
        IF2[n] = -2.0*(AAmp[n]*AAmp[n]+EAmp[n]*EAmp[n])/SN[n];
        
        //fprintf(out,"%e %e %e %e %e\n", FF[n], (AA[n]*AA[n])/SN[n], (AAmp[n]*AAmp[n])/SN[n], (AA[n]*AAmp[n]*ca)/SN[n], ca);
        
        /*
        y = ca;
        if(y*yold < 0.0)
        {
            nt++;
            yold = y;
        }
        */

        x = (FF[n]-FF[n-1])*(IF1[n]+IF2[n]);  // contribution to the integrand
        
        if(fabs(x) > 2.0)  // only care about phase differences in regions that contribute signifcantly
        {
            z = ca;
            if(z*zold < 0.0)
            {
                nn++;
                zold = z;
            }
        }
        
    }
    
    //fclose(out);
    
    logL = 0.0;
    
    if(nn < 10)
    {

    gsl_interp_accel *Iacc = gsl_interp_accel_alloc();
    gsl_spline *Ispline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(Ispline, FF, IF1, NF);
    logL =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_init(Ispline, FF, IF2, NF);
    logL +=  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
    gsl_spline_free(Ispline);
    gsl_interp_accel_free(Iacc);
        
    }
    
  //  printf("%d %d %e %e %e\n", nt, nn, zz, x, logL);
    

    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(TF);
    free(AF);
    free(PF);
    free(IF1);
    free(IF2);
    
    return(logL);

}


double tvol(double *params)
{
    double *paramX;
    double **Fish;
    double **Stencil;
    double epsilon=1.0e-5;
    double tv;
    int i, j;
    
    
    paramX = double_vector(NP);

    
    Fish = double_matrix(2,2);
    Stencil = double_matrix(3,3);
    
    for (i=0; i< NP; i++) paramX[i] = params[i];
    
    for (i=-1; i<= 1; i++)
    {
        paramX[7] = params[7]+(double)(i)*epsilon;
        for (j=-1; j<= 1; j++)
        {
            paramX[8] = params[8]+(double)(j)*epsilon;
            Stencil[i+1][j+1] = Tmerger(paramX,paramX[5]);
        }
    }
    
    Fish[0][0] = -(Stencil[2][1]-2.0*Stencil[1][1]+Stencil[0][1])/(epsilon*epsilon);
    Fish[1][1] = -(Stencil[1][2]-2.0*Stencil[1][1]+Stencil[1][0])/(epsilon*epsilon);
    Fish[0][1] = -(Stencil[2][2]-Stencil[2][0]-Stencil[0][2]+Stencil[0][0])/(epsilon*epsilon);
    Fish[1][0] = Fish[0][1];
    
    //printf("%e %e %e\n", Fish[0][0], Fish[1][1], Fish[0][1]);
    
    tv = 1.0/sqrt(fabs(Fish[0][0]*Fish[1][1]-Fish[0][1]*Fish[1][0]));
    
    free_double_matrix(Fish,2);
    free_double_matrix(Stencil,3);
    free_double_vector(paramX);
    
    return(tv);
    
}

void InChl(int ll, double *params, double **Fisher, double **iChl)
{
    int i, j, k, kk;
    double **CovI;
    double *scale;
    
    scale = (double*)malloc(sizeof(double)* (NI));
    
    CovI = double_matrix(NI,NI);

    Inverse(Fisher, CovI, NI);

   // Fisher uses log derivatives on m1, m2,  Want to switch to linear derivatives since using
   // uniform prior in m1, m2

    for(j = 0; j < NI; j++) scale[j] = 1.0;
    
    if(ll==0)
    {
    scale[0] = params[0];
    scale[1] = params[1];
    }
    else
    {
    scale[0] = exp(params[0]);
    scale[1] = exp(params[1]);
    }

    for(j = 0; j < NI; j++)
      {
       for(k=0; k < NI; k++)
       {
           CovI[j][k] *= scale[j]*scale[k];
       }
    }

    cholesky(CovI, iChl, NI);
  
    free(scale);
    free_double_matrix(CovI,NI);

    
}

void Ext_In(int ll, double *params, double **Fisher, double **eChl, double **iChl)
{
    int i, j, k, kk;
    double **FishE;
    double **CovE, **CovI;
    double *scale;
    
    kk = NP-NE;
    
    scale = (double*)malloc(sizeof(double)* (NP));
    
    FishE = double_matrix(NE,NE);
    CovE = double_matrix(NE,NE);
    CovI = double_matrix(NI,NI);
    
       for(j = 0; j < NE; j++)
         {
          for(k=0; k < NE; k++)
          {
              FishE[j][k] = Fisher[j+kk][k+kk];
          }
         }

    Inverse(Fisher, CovI, NI);
    Inverse(FishE, CovE, NE);

   // Fisher uses log derivatives on m1, m2, DL. Want to switch to linear derivatives since using
   // uniform prior in m1, m2 and DL

    for(j = 0; j < NP; j++) scale[j] = 1.0;
    
    if(ll==0)
    {
    scale[0] = params[0];
    scale[1] = params[1];
    scale[6] = params[6];
    }
    else
    {
    scale[0] = exp(params[0]);
    scale[1] = exp(params[1]);
    scale[6] = exp(params[6]);
    }

    for(j = 0; j < NI; j++)
      {
       for(k=0; k < NI; k++)
       {
           CovI[j][k] *= scale[j]*scale[k];
       }
    }
    
    for(j = 0; j < NE; j++)
      {
       for(k=0; k < NE; k++)
       {
           CovE[j][k] *= scale[j+kk]*scale[k+kk];
       }
    }


    cholesky(CovI, iChl, NI);
    cholesky(CovE, eChl, NE);
    
    free(scale);
    free_double_matrix(CovI,NI);
    free_double_matrix(CovE,NE);
    free_double_matrix(FishE,NE);
    
}


void FisherPlot(int ll, double Tend, double *params)
{
    int i, j, k;
    
    double **Fisher;
    double **Cov, **Chl;
    double *zv, *scale, *pref;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);

    FILE *out;
    
    Fisher = double_matrix(NP,NP);
    Cov = double_matrix(NP,NP);
    Chl = double_matrix(NP,NP);
    zv = double_vector(NP);
    
    pref = (double*)malloc(sizeof(double)* (NP));
    scale = (double*)malloc(sizeof(double)* (NP));
    
    for (i=0; i< NP; i++) pref[i] = params[i];
    
    if(ll == 1)
    {
    pref[0] = exp(params[0]);
    pref[1] = exp(params[1]);
    pref[6] = exp(params[6]);
    }
    
    FisherFast(ll, Tend, params, Fisher);
    
    printf("\n");
    for(j = 0; j < NP; j++)
     {
      for(k=0; k < NP; k++)
      {
          printf("%e ", Fisher[j][k]);
      }
         printf("\n");
    }
    printf("\n");
    
    Inverse(Fisher, Cov, NP);

  // Fisher uses log derivatives on m1, m2, DL. Want to switch to linear derivatives since using
 // uniform prior in m1, m2 and DL

     for(j = 0; j < NP; j++) scale[j] = 1.0;
     scale[0] = pref[0];
     scale[1] = pref[1];
     scale[6] = pref[6];
 
     for(j = 0; j < NP; j++)
       {
        for(k=0; k < NP; k++)
        {
            Cov[j][k] *= scale[j]*scale[k];
        }
      }
 
 
      cholesky(Cov, Chl, NP);
 
 
      out = fopen("fisher.dat","w");
 
      i = 0;
      k = 0;
       do
       {
 
        for(j = 0; j < NP; j++) zv[j] = gsl_ran_gaussian(r,1.0);
     
        for(j = 0; j < NP; j++) params[j] = pref[j];
 
       for(j = 0; j < NP; j++)
        {
         for(k=0; k < NP; k++)
         {
         params[j] += Chl[j][k]*zv[k];
         }
       }
     
     // maximum 100% error in distance and spins, cosines in physical range
           if(params[0] > 5.0e3 && params[1] > 5.0e3 && fabs(params[2]) < 1.0 && fabs(params[3]) < 1.0 && params[6] > 0.1 && params[6] < 1000.0 && fabs(params[7]) < 1.0 && fabs(params[10]) < 1.0)
     {
         
         while(params[4] < 0.0) params[4] += PI;
         while(params[4] > PI) params[4] -= PI;
         while(params[9] < 0.0) params[9] += PI;
         while(params[9] > PI) params[9] -= PI;
         while(params[8] < 0.0) params[8] += TPI;
         while(params[8] > TPI) params[8] -= TPI;
     
            fprintf(out, "%d %e ", k, 0.0);
            for(j = 0; j < NP; j++) fprintf(out, "%.15e ", params[j]);
            fprintf(out,"\n");
         
         k++;
  
     }
 
       }while(i < 10000000 && k < 10000);
    
    fclose(out);
    
    for (i=0; i< NP; i++) params[i] = pref[i];
     
     if(ll == 1)
     {
     params[0] = exp(pref[0]);
     params[1] = exp(pref[1]);
     params[6] = exp(pref[6]);
     }
    
    

}


void MCMC(double *params, double Tend, long N, double *AD, double *ED, double *SN)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    double kxm, *logLx, logLy, SNR, x, y, z, Lref, fc, fe;
    double *nhx;
    double lm1, lm2, chi1, chi2;
    double ctheta, phi;
    double m1, m2, DL;
    int i, j, k, km, n, q, mc, NF;
    int NFmax = 5000;
    double *AF, *PF, *FF, *TF, *SAE, *FS, *ASD;
    double *FN, *ASN;
    int NH=1000;
    int MP;
    int M = 1000000;
    int *who;
    double *heat;
    double **paramx, **paramy;
    double *pnew;
    double *pref, *scale;
    double ***Fisher;
    double ***history;
    double *min, *max;
    double **ejump, ***evec, **diag;
    double **ejumpI, ***evecI;
    int ll;
    double alpha, beta, tau, f, df;
    double phic, tc, px, tm;
    int *m;
    int *scount, *sacc, hold;
    int mcount, macc;
    int **av, **cv;
    
    double *ATR, *ATI, *ETR, *ETI;
    double *AS, *ES;
    double *ASR, *ESR;
    
    double **Cov, **Chl;
    double *zv;
    
    double ***iChl, ***eChl;
    
    clock_t start, end;
    double cpu_time_used;
    
    FILE *in;
    FILE *chain;
    FILE *out;
    FILE *levels;
    FILE *swaps;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    MP = 6;
    
    eChl = double_tensor(NC,NE,NE);
    iChl = double_tensor(NC,NI,NI);
    
    pnew = (double*)malloc(sizeof(double)* (NP));
    pref = (double*)malloc(sizeof(double)* (NP));
    scale = (double*)malloc(sizeof(double)* (NP));
    
    for (i=0; i< NP; i++) pref[i] = params[i];

    // switch to log for masses and distance
    
    ll = 1;
    params[0] = log(params[0]);
    params[1] = log(params[1]);
    params[6] = log(params[6]);
    
    while(params[4] > PI)   params[4] -= PI;
    while(params[4] < 0.0)  params[4] += PI;
    while(params[8] > TPI)  params[8] -= TPI;
    while(params[8] < 0.0)  params[8] += TPI;
    while(params[9] > PI)   params[9] -= PI;
    while(params[9] < 0.0)  params[9] += PI;
    
    if(params[1] > params[0])  // catch if m2 > m1 and flip
    {
        lm1 = params[1];
        chi1 = params[3];
        lm2 = params[0];
        chi2 = params[2];
        params[0] = lm1;
        params[1] = lm2;
        params[2] = chi1;
        params[3] = chi2;
    }

    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FS = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
    
    SetUp(ll, params, NFmax, &NF, FS, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    SAE = (double*)malloc(sizeof(double)* (NF));
    ASD = (double*)malloc(sizeof(double)* (NF));
    ASN = (double*)malloc(sizeof(double)* (NFFT/2));
    
    FN = (double*)malloc(sizeof(double)* (MF));

    sacc = int_vector(NC);
    scount = int_vector(NC);
    who = int_vector(NC);
    heat = double_vector(NC);
    logLx = double_vector(NC);
    nhx = double_vector(NC);
    paramx = double_matrix(NC,NP);
    paramy = double_matrix(NC,NP);
    history = double_tensor(NC,NH,NP);
    Fisher = double_tensor(NC,NP,NP);
    Cov = double_matrix(NP,NP);
    Chl = double_matrix(NP,NP);
    zv = double_vector(NP);
    ejump = double_matrix(NC,NP);
    evec = double_tensor(NC,NP,NP);
    ejumpI = double_matrix(NC,NP);
    evecI = double_tensor(NC,NP,NP);
    diag = double_matrix(NC,NP);

    Extrinsic(params, Tend, NF, FS, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    FF = (double*)malloc(sizeof(double)* (NF+1));
    for(n=0; n< NF; n++) FF[n] = FS[n];
    
    free(FS);
    free(TF);
    free(PF);
    free(AF);
    
    /*
    out = fopen("ft.dat", "w");
    for(n=0; n< NF; n++)
       {
           fprintf(out,"%e %e\n", TF[n], FF[n]);
       }
    fclose(out);
    */
    
    // Have to add in the merger time and phase since they are extrinsic
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    for(n=0; n< NF; n++)
    {
        px = TPI*FF[n]*(Tobs-tc + dt/2.0)-2.0*phic;
        //printf("%d %e\n", n, APhase[n]);
        APhase[n] += px;
        EPhase[n] += px;
    }
    
    // use additional element in FF array to store maximum frequency of the reference signal
    FF[NF] = FF[NF-1];
    if(Tend < tc)
    {
     FF[NF] = FofT(ll, params, Tend);
    }
    
    
    for (i=0; i< NF; i++)
    {
      instrument_noise(FF[i], &SAE[i], &x);
      ASD[i] = sqrt(SAE[i]);
    }
    
    df = (1.0/dt)/(double)(MF);
    for (i=0; i< MF; i++)
    {
       FN[i] = (double)(i)*df;
    }
    
    ATR = (double*)malloc(sizeof(double)* (MF));
    ATI = (double*)malloc(sizeof(double)* (MF));
    ETR = (double*)malloc(sizeof(double)* (MF));
    ETI = (double*)malloc(sizeof(double)* (MF));
    
    ASR = (double*)malloc(sizeof(double)* (NFFT));
    ESR = (double*)malloc(sizeof(double)* (NFFT));
    AS = (double*)malloc(sizeof(double)* (NFFT));
    ES = (double*)malloc(sizeof(double)* (NFFT));
    
    double fr;
    
    fr = 1.0e-5;
    
    ResponseFreq(ll, Tend, params, NFFT, ASR, ESR);
    
     // subtract the reference signal
     for (i=0; i< NFFT; i++)
     {
         AD[i] -= ASR[i];
         ED[i] -= ESR[i];
     }
    
    // whiten the resdiual
    AD[0] = 0.0;
    ED[0] = 0.0;
    for (i=1; i< NFFT/2; i++)
    {
        ASN[i] = sqrt(SN[i]);
        AD[i] /= ASN[i];
        ED[i] /= ASN[i];
        AD[NFFT-i] /= ASN[i];
        ED[NFFT-i] /= ASN[i];
    }
    
    Hetrodyne(ll, Tend, params, NFFT, ATR, ATI, ETR, ETI, AD, ED, SN);
    
    Lref = (fourier_nwip(ASR, ASR, SN, NFFT)/Tobs + fourier_nwip(ESR, ESR, SN, NFFT)/Tobs)/2.0;
    x = nh(ll, Tend, params, ATR, ATI, ETR, ETI, NF, FF, APhase, EPhase, FN, ASD);
    y = fourier_nwip(AD, ASR, ASN, NFFT)/Tobs+ fourier_nwip(ED, ESR, ASN, NFFT)/Tobs;
    printf("%f %f %f\n", Lref, x, y);
    
    //return;
    
    nhx[0] = 0.0;
    if(nflag == 1)
    {
     nhx[0] = -ndh(ll, Tend, params, ATR, ATI, ETR, ETI, NF, FF, AAmp, EAmp, APhase, EPhase, FN, ASD);
    }
    
    logLx[0] = 0.0;
    if(lhold == 0) logLx[0] = LikelihoodDelta(ll, Tend, params, NF, FF, AAmp, EAmp, APhase, EPhase, SAE) + nhx[0];
    
    
    start = clock();
    for(i=0; i<1000; i++) x = LikelihoodDelta(ll, Tend, params, NF, FF, AAmp, EAmp, APhase, EPhase, SAE);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("likelihood calculation took %f seconds\n", cpu_time_used/1000.0);
    
    
    start = clock();
    for(i=0; i<1000; i++) x = ndh(ll, Tend, params, ATR, ATI, ETR, ETI, NF, FF, AAmp, EAmp, APhase, EPhase, FN, ASD);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("(n|h) calculation took %f seconds\n", cpu_time_used/1000.0);
    
    
    
    
    printf("%f %f\n", logLx[0], nhx[0]);
    for (i=1; i< NC; i++)
    {
    logLx[i] = logLx[0];
    nhx[i] = nhx[0];
    }
    
    
    for (i=0; i< NC; i++) who[i] = i;
    heat[0] = 1.0;
    // Hot chain should have an effective SNR of ~5
    SNR = sqrt(logLx[0]*2.0);
    x = pow((SNR/5.0),1.0/(double)(NC-1));
    for (i=1; i< NC; i++) heat[i] = heat[i-1]*x;
    printf("SNR %f increment %f max heat %f SNReff = %f\n", SNR, x, heat[NC-1], SNR/heat[NC-1]);
    
    //return;
   
    // prior boundaries
    max = (double*)malloc(sizeof(double)* (NP));
    min = (double*)malloc(sizeof(double)* (NP));
    
    // [0] ln Mass1  [1] ln Mass2  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln distance
    // [7] cos EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] cos inclination
    
    max[0] = log(1.0e8);
    max[1] = log(1.0e8);
    max[2] = 0.999;
    max[3] = 0.999;
    max[4] = PI;
    max[5] = 2.0*Tobs;
    max[6] = log(400.0);
    max[7] = 1.0;
    max[8] = 2.0*PI;
    max[9] = PI;
    max[10] = 1.0;
    
    
    min[0] = log(5.0e4);
    min[1] = log(5.0e4);
    min[2] = -0.999;
    min[3] = -0.999;
    min[4] = 0.0;
    min[5] = 0.0;
    min[6] = log(0.1);
    min[7] = -1.0;
    min[8] = 0.0;
    min[9] = 0.0;
    min[10] = -1.0;

    
    for (i=0; i< NC; i++)
    {
        for (j=0; j< NP; j++) paramx[i][j] = params[j];
        
        for (k=0; k< NH; k++)
        {
            for (j=0; j< NP; j++) history[i][k][j] = params[j];
        }
    }
    

    
    FisherFast(ll, Tend, paramx[0], Fisher[0]);
    
    for (i=1; i< NC; i++)
    {
        for (k=0; k< NP; k++)
        {
            
            for (j=0; j< NP; j++) Fisher[i][k][j] = Fisher[0][k][j];
        }
        
    }
    
    for (i=0; i< NC; i++)
    {
        for (k=0; k< NP; k++) diag[i][k] = 1.0/sqrt(Fisher[i][k][k]);
    }
    
    FisherEvec(Fisher[0], ejump[0], evec[0], NP);
    FisherEvec(Fisher[0], ejumpI[0], evecI[0], 6);
    

    

    for (i=1; i< NC; i++)
    {
        
        for (k=0; k< NP; k++)
        {
            ejump[i][k] = ejump[0][k];
            for (j=0; j< NP; j++) evec[i][k][j] = evec[0][k][j];
        }
        
        for (k=0; k< NI; k++)
        {
            ejumpI[i][k] = ejumpI[0][k];
            for (j=0; j< NI; j++) evecI[i][k][j] = evecI[0][k][j];
        }
        
    }

    av = int_matrix(4,NC);
    cv = int_matrix(4,NC);
    
    for(j = 0; j < 4; j++)
    {
        for(k=0; k < NC; k++)
        {
            av[j][k] = 0;
            cv[j][k] = 1;
        }
    }
    
    m = int_vector(NC);
    
    for(k=0; k < NC; k++)
    {
        m[k] = 0;
        sacc[k] = 0;
        scount[k] = 0;
    }
    
    mcount = 1;
    macc = 0;
    
    // wind the randoms
    for(i = 0; i < 3; i++)
    {
    for(k=0; k < NC; k++)
    {
      x = gsl_ran_gaussian(r,1.0);
      x = gsl_ran_gaussian(rvec[k],1.0);
    }
    }
    
    
    swaps = fopen("swaps.dat","w");
    levels = fopen("levels.dat","w");
    chain = fopen("chain.dat","w");
    
    for(mc = 1; mc <= M; mc++)
    {
        
        if(mc%5000==0  && lhold == 0)
        {
            // update the Fisher matrices
            #pragma omp parallel for
            for(k = 0; k < NC; k++)
            {
                FisherFast(ll, Tend, paramx[k], Fisher[k]);
                FisherEvec(Fisher[k], ejump[k], evec[k], NP);
                FisherEvec(Fisher[k], ejumpI[k], evecI[k], NI);
            }
            
            for(k = 0; k < NC; k++)
            {
                for(i = 0; i < NP; i++) diag[k][i] = 1.0/sqrt(Fisher[k][i][i]);
            }
            
        }
        
        alpha = gsl_rng_uniform(r);
        
        if((NC > 1) && (alpha < 0.2))  // decide if we are doing a MCMC update of all the chains or a PT swap
        {
            
            // chain swap
           
            alpha = (double)(NC-1)*gsl_rng_uniform(r);
            j = (int)(alpha);
            beta = exp((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
            alpha = gsl_rng_uniform(r);
            if(beta > alpha)
            {
                hold = who[j];
                who[j] = who[j+1];
                who[j+1] = hold;
                sacc[j]++;
            }
            
             scount[j]++;
            
        }
        else      // MCMC update
        {
            
            mcount++;
            
            for(j = 0; j < NC; j++)
            {
                for(i = 0; i < NP; i++)
                {
                    paramy[j][i] = paramx[j][i];
                }
            }

        #pragma omp parallel for
        for(k=0; k < NC; k++)
        {
          update(k, ll, Tend, logLx, nhx, paramx, paramy, min, max, who, heat, history, SAE, ejump, evec, ejumpI, evecI,
                 diag, iChl, tm, NH, NF, FF, AAmp, EAmp, APhase, EPhase, cv, av, ATR, ATI, ETR, ETI, FN, ASD, rvec[k]);
         }
            
            // add to the history file
            if(mc%10 == 0)
            {
                for(k=0; k < NC; k++)
                {
                    q = who[k];
                    i = m[k]%NH;
                    // the history file is kept for each temperature
                    for(j=0; j<NP; j++) history[k][i][j] = paramx[q][j];
                    m[k]++;
                }
            }
            
        }
        

        if(mc > 0 && mc%20 == 0)
        {
                q = who[0];
                m1 = exp(paramx[q][0]);
                m2 = exp(paramx[q][1]);
                DL = exp(paramx[q][6]);
                
                fprintf(chain,"%d %.12e %.12e %.12e ", mc, logLx[q], m1, m2);
                for(i = 2; i < NP; i++)
                 {
                    if(i == 6)
                     {
                     fprintf(chain, "%.12e ", DL);
                     }
                     else
                     {
                      fprintf(chain, "%.15e ", paramx[q][i]);
                     }
                 }
                fprintf(chain,"%d\n", q);
            
            fprintf(swaps, "%d ", mc);
            for(k = 0; k < NC-1; k++)
            {
                fprintf(swaps, "%f ", (double)(sacc[k])/(double)(scount[k]));
            }
            fprintf(swaps, "\n");
            
            for(k = 0; k < NC; k++)
            {
                fprintf(levels, "%.12e ", logLx[who[k]]);
            }
            fprintf(levels, "\n");
            
            
        }
        
        if(mc%100 == 0)
        {
            

            q = who[0];
            m1 = exp(paramx[q][0]);
            m2 = exp(paramx[q][1]);
            
            
            printf("%d %e %e %e %f %f %f %f\n", mc, logLx[q], m1, m2,
            
                   (double)(sacc[0])/(double)(scount[0]),
                   (double)(av[0][q])/(double)(cv[0][q]),
                   (double)(av[1][q])/(double)(cv[1][q]),
                   (double)(av[2][q])/(double)(cv[2][q]));
        }

        
    }
    
    fclose(chain);
    fclose(levels);
    fclose(swaps);
    

    free_int_matrix(av,4);
    free_int_matrix(cv,4);
    free_int_vector(m);
    free_int_vector(who);
    free_int_vector(sacc);
    free_int_vector(scount);
    free_double_vector(heat);
    free_double_vector(logLx);
    free_double_vector(nhx);
    free_double_matrix(paramx,NC);
    free_double_matrix(paramy,NC);
    free_double_tensor(history,NC,NH);
    free_double_tensor(Fisher,NC,NP);
    free_double_matrix(ejump,NC);
    free_double_tensor(evec,NC,NP);
    free_double_matrix(ejumpI,NC);
    free_double_tensor(evecI,NC,NP);
    free_double_matrix(diag,NC);
    
    free_double_tensor(iChl,NC,NI);
    free_double_tensor(eChl,NC,NE);
    
    free(FN);
    free(ASN);

    
    free(SAE);
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(FF);

    
    
}



void update(int k, int ll, double Tend, double *logLx, double *nhx, double **paramx, double **paramy, double *min, double *max, int *who, double *heat, double ***history, double *SN, double **ejump, double ***evec, double **ejumpI, double ***evecI, double **diag, double ***iChl, double tm, int NH, int NF, double *FF, double *AA, double *EA, double *AP, double *EP, int **cv, int **av, double *ATR, double *ATI, double *ETR, double *ETI, double *FN, double *ASD, gsl_rng *r)
{
    int q, i, j;
    int fstat;
    double qx, qy, a, b, c, tx, nhy;
    double alpha, beta, DL, pDx, pDy, H;
    double logLy, eta, leta, pMcMx, pMcMy, lglx;
    double w, v;
    double lm1, lm2, chi1, chi2;
    int typ, flag, id1, id2;
    int Ndraw, Nden;
    long itheta, iphi, iDL, iiota, ipsi, iphic, itc, label;
    double dtheta, dphi, dDL, diota, dpsi, dphic, dtc;
    double cth, phi;
    double size, x, y, z;
    double Lf;
    double *jump, *zv;
    double *u;
    
    // set to 1 to use Fstat likelihood
    fstat = 0;
    
    u = double_vector(5);
    
    jump = double_vector(NI);
    zv = double_vector(NI);

    a = 0.5;
    b = 0.2;
    
    q = who[k];
    
    qx = qy = 0.0;    // log proposal densities
    
    alpha = gsl_rng_uniform(r);
    
    if(alpha > a) // fisher jump
    {
        typ = 0;
        
        alpha = gsl_rng_uniform(r);
        
        if(alpha > 0.8)  // all parameters
        {
        
        // pick an eigendirection to jump in
        beta = gsl_rng_uniform(r);
        // we know that 0, 1, 2 are bad
        i = 3+(int)(beta*(NP-3));
        //i = (int)(beta*(NP));
        // draw the jump size
        beta = sqrt(heat[k])*ejump[q][i]*gsl_ran_gaussian(r,1.0);
        for(j = 0; j < NP; j++) paramy[q][j] = paramx[q][j]+beta*evec[q][i][j];
            
            
        }
        else if (alpha > 0.2)  // just parameters in the phase
        {
            
            for(j = 0; j < NP; j++) paramy[q][j] = paramx[q][j];
            
            
          
          beta = gsl_rng_uniform(r);
          i = (int)(beta*(NI));
          beta = sqrt(heat[k])*ejumpI[q][i]*gsl_ran_gaussian(r,1.0);
          for(j = 0; j < NI; j++) paramy[q][j] += beta*evecI[q][i][j];
            
            flag = 0;
            for(j = 0; j < 4; j++)
            {
                if(paramy[q][j] > max[j] || paramy[q][j] < min[j]) flag = 1;
            }
            
            if(i==0 && flag ==0)
            {
                    int *pmap;
                    double **Fish;
                    double *params;
                    double **Svec;
                    double *Sval;
                
                if(paramy[q][1] > paramy[q][0])  // catch if m2 > m1 and flip
                {
                    lm1 = paramy[q][1];
                    chi1 = paramy[q][3];
                    lm2 = paramy[q][0];
                    chi2 = paramy[q][2];
                    paramy[q][0] = lm1;
                    paramy[q][1] = lm2;
                    paramy[q][2] = chi1;
                    paramy[q][3] = chi2;
                }
                    
                    Svec = double_matrix(3,3);
                    Sval = double_vector(3);
                    
                    pmap = (int*)malloc(sizeof(int)* (NP));
                    Fish = double_matrix(3,3);
                    params = (double*)malloc(sizeof(double)* (NP));
                    
                    for(j = 0; j < NP; j++) pmap[j] = -1;
                    // Here we only vary tc, phic and DL
                    pmap[4] = 0;
                    pmap[5] = 1;
                    pmap[6] = 2;
                
                    for(j = 0; j < NP; j++) params[j] = paramx[q][j];
        
                    lglx = LikelihoodDeltaMaxT(ll, Tend, params, NF, FF, AA, EA, AP, EP, SN);
                
                     FisherSub(ll, Tend, pmap, params, Fish);
                 
                 y = logLx[q]-lglx;
                 
                 qx = 0.0;
                 for (i=0; i< NP; i++)
                 {
                     if(pmap[i] > -1)
                     {
                       for (j=0; j< NP; j++)
                        {
                         if(pmap[j] > -1) qx -= 0.5*Fish[pmap[i]][pmap[j]]*(paramx[q][i]-params[i])*(paramx[q][j]-params[j]);
                        }
                     }
                 }
                 
                 //printf("%f %f\n", y, qx);
                 
                 x = det(Fish,3)/(heat[k]*heat[k]*heat[k]);
                 
                 qx /= heat[k];
                 
                // printf("%f %f\n", qx, 0.5*log(x));
                 
                 qx += 0.5*log(x);
                 
                 
                 for(j = 0; j < NP; j++) params[j] = paramy[q][j];
                 
                 Lf = LikelihoodDeltaMaxT(ll, Tend, params, NF, FF, AA, EA, AP, EP, SN);
                 
                 FisherSub(ll, Tend, pmap, params, Fish);
                 
                 FisherEvec(Fish, Sval, Svec, 3);
                 
                 for(j = 0; j < NP; j++) paramy[q][j] = params[j];
                 
                 //w = LikelihoodDelta(ll, Tend, paramy[q], NF, FF, AA, EA, AP, EP, SN);
                 
                 
                 // pick an eigendirection to jump in
                 beta = gsl_rng_uniform(r);
                 i = (int)(beta*3);
                 // draw the jump size
                 beta = sqrt(heat[k])*Sval[i]*gsl_ran_gaussian(r,1.0);
                 for(j = 0; j < NP; j++)
                 {
                   if(pmap[j] > -1) paramy[q][j] = params[j]+beta*Svec[i][pmap[j]];
                 }
                 
                 //v = LikelihoodDelta(ll, Tend, paramy[q], NF, FF, AA, EA, AP, EP, SN);
                 
                 qy = 0.0;
                 for (i=0; i< NP; i++)
                 {
                     if(pmap[i] > -1)
                     {
                         for (j=0; j< NP; j++)
                         {
                             if(pmap[j] > -1) qy -= 0.5*Fish[pmap[i]][pmap[j]]*(paramy[q][i]-params[i])*(paramy[q][j]-params[j]);
                         }
                     }
                 }
                 
                 //printf("%f %f\n", v-w, qy);
                 
                 x = det(Fish,3)/(heat[k]*heat[k]*heat[k]);
                 
                 qy /= heat[k];
                 
                 // printf("%f %f\n", qy, 0.5*log(x));
                 
                 qy += 0.5*log(x);
                
                
                free_double_matrix(Fish,3);
                free(pmap);
                free(params);
                free_double_matrix(Svec,3);
                free_double_vector(Sval);
                
            }
            
        }
        else  // no correlation jump
        {
        x = sqrt(heat[k]/(double)(NP));
        for(j = 0; j < NP; j++) paramy[q][j] = paramx[q][j]+x*diag[q][j]*gsl_ran_gaussian(r,1.0);
        }
         
        
        tx = -1.0;
        
    }
    else if(alpha > b)// differential evolution
    {
        typ = 1;
        // the history file is kept for each temperature
        
        de_jump(paramx[q], paramy[q], history[k], NH, NP, r);
        
        tx = Tmerger(paramx[q],paramx[q][5]);
        
    }
    else  // big sky
    {
        typ = 2;
        
         if(fstat == 1)
         {
             
             for(j = 0; j < NP; j++) paramy[q][j] = paramx[q][j];
             beta = gsl_rng_uniform(r);
             if(beta > 0.7)
             {
             paramy[q][7] = -1.0+2.0*gsl_rng_uniform(r);
             paramy[q][8] = TPI*gsl_rng_uniform(r);
             }
             else if (beta > 0.3)
             {
                 paramy[q][7] = paramx[q][7] + gsl_ran_gaussian(r,0.1);
                 paramy[q][8] = paramx[q][8] + gsl_ran_gaussian(r,0.1);
             }
             else
             {
                 paramy[q][7] = paramx[q][7] + gsl_ran_gaussian(r,0.05);
                 paramy[q][8] = paramx[q][8] + gsl_ran_gaussian(r,0.05);
             }
             
         }
         else
         {
             
             beta = gsl_rng_uniform(r);
             if(beta > 0.7)
             {
                 cth = -1.0+2.0*gsl_rng_uniform(r);
                 phi = TPI*gsl_rng_uniform(r);
             }
             else if (beta > 0.3)
             {
                 cth = paramx[q][7] + gsl_ran_gaussian(r,0.1);
                 phi = paramx[q][8] + gsl_ran_gaussian(r,0.1);
             }
             else
             {
                 cth = paramx[q][7] + gsl_ran_gaussian(r,0.05);
                 phi = paramx[q][8] + gsl_ran_gaussian(r,0.05);
             }
             
             if(phi < 0.0) phi += TPI;
             if(phi > TPI) phi -= TPI;
             

         if(fabs(cth) < 1.0)
          {
        
        int *pmap;
        double **Fish;
        double *params;
        double **Svec;
        double *Sval;
        
    
        Svec = double_matrix(4,4);
        Sval = double_vector(4);
        
        pmap = (int*)malloc(sizeof(int)* (NP));
        Fish = double_matrix(4,4);
        params = (double*)malloc(sizeof(double)* (NP));
        
        pmap[0] = pmap[1] = pmap[2] = pmap[3] = -1;
        pmap[5] = pmap[7] = pmap[8] = -1;
        pmap[4] = 0;
        pmap[6] = 1;
        pmap[9] = 2;
        pmap[10] = 3;
        
        
        for(j = 0; j < NP; j++) params[j] = paramx[q][j];
        tx = -1.0;  // time offset will be at correct value, no need to maximize
        
        lglx = LikelihoodFstat(ll, Tend, params, tx, NF, FF, AA, EA, AP, EP, SN);
        
        // The signal is invariant under these shifts. Find out which one gives
        // the smallest parameter deltas
        j = 0;
        u[0] = fabs(params[4]-paramx[q][4]) + fabs(params[9]-paramx[q][9]);
        u[1] = fabs(params[4]+PI/2-paramx[q][4]) + fabs(params[9]+PI/2-paramx[q][9]);
        u[2] = fabs(params[4]+PI/2-paramx[q][4]) + fabs(params[9]-PI/2-paramx[q][9]);
        u[3] = fabs(params[4]-PI/2-paramx[q][4]) + fabs(params[9]+PI/2-paramx[q][9]);
        u[4] = fabs(params[4]-PI/2-paramx[q][4]) + fabs(params[9]-PI/2-paramx[q][9]);
        
        j = 0;
        y = u[0];
        for (i=1; i< 5; i++)
        {
           if(u[i] < y)
            {
                y = u[i];
                j = i;
            }
        }
        
        if(j == 1)
        {
            params[4] += PI/2;
            params[9] += PI/2;
        }
        
        if(j == 2)
        {
            params[4] += PI/2;
            params[9] -= PI/2;
        }
        
        if(j == 3)
        {
            params[4] -= PI/2;
            params[9] += PI/2;
        }
        
        if(j == 4)
        {
            params[4] -= PI/2;
            params[9] -= PI/2;
        }

        
       // printf("%f %f\n", params[4]-paramx[q][4], params[9]-paramx[q][9]);
        
        FisherSub(ll, Tend, pmap, params, Fish);
        
        //y = logLx[q]-lglx;
        
        qx = 0.0;
        for (i=0; i< NP; i++)
        {
            if(pmap[i] > -1)
            {
              for (j=0; j< NP; j++)
               {
                if(pmap[j] > -1) qx -= 0.5*Fish[pmap[i]][pmap[j]]*(paramx[q][i]-params[i])*(paramx[q][j]-params[j]);
               }
            }
        }
        
        //printf("%f %f\n", y, qx);
        
        x = det(Fish,4)/(heat[k]*heat[k]*heat[k]*heat[k]);
        
        
        
        qx /= heat[k];
        
       // printf("%f %f\n", qx, 0.5*log(x));
        
        qx += 0.5*log(x);
        
        
        for(j = 0; j < NP; j++) params[j] = paramx[q][j];
        params[7] = cth;
        params[8] = phi;
        
        tx = Tmerger(paramx[q],paramx[q][5]);
        Lf = LikelihoodFstat(ll, Tend, params, tx, NF, FF, AA, EA, AP, EP, SN);
        
        FisherSub(ll, Tend, pmap, params, Fish);
        
        FisherEvec(Fish, Sval, Svec, 4);
        
        for(j = 0; j < NP; j++) paramy[q][j] = params[j];
        
        //w = LikelihoodDelta(ll, Tend, paramy[q], NF, FF, AA, EA, AP, EP, SN);
        
        
        // pick an eigendirection to jump in
        beta = gsl_rng_uniform(r);
        i = (int)(beta*4);
        // draw the jump size
        beta = sqrt(heat[k])*Sval[i]*gsl_ran_gaussian(r,1.0);
        for(j = 0; j < NP; j++)
        {
          if(pmap[j] > -1) paramy[q][j] = params[j]+beta*Svec[i][pmap[j]];
        }
        
        //v = LikelihoodDelta(ll, Tend, paramy[q], NF, FF, AA, EA, AP, EP, SN);
        
        qy = 0.0;
        for (i=0; i< NP; i++)
        {
            if(pmap[i] > -1)
            {
                for (j=0; j< NP; j++)
                {
                    if(pmap[j] > -1) qy -= 0.5*Fish[pmap[i]][pmap[j]]*(paramy[q][i]-params[i])*(paramy[q][j]-params[j]);
                }
            }
        }
        
        //printf("%f %f\n", v-w, qy);
        
        x = det(Fish,4)/(heat[k]*heat[k]*heat[k]*heat[k]);
        
        qy /= heat[k];
        
        // printf("%f %f\n", qy, 0.5*log(x));
        
        qy += 0.5*log(x);
        
       // printf("%f %f\n", qx, qy);
        
       // printf("%f %f %f %f %f ", qx, qy, lglx, logLx[q], Lf);
        
        
        // Need to include a Jacobian that accounts for the deteministic time mapping.
        // The Tmerger function is used to re-set tc so that the merger occurs at the
        // same time in the detector. Depending on the sky location, the time range dt
        // mapped to by a unit dcostheta dphi around the reference location will be different.
        
        // Test were unclear on wether this kind of term was needed
        
        // tvol is the time volume surrounding theta, phi
        //qy += log(tvol(paramy[q]));
        //qx += log(tvol(paramx[q]));
        
        
        // To cover the full range, apply a pi shift to 2psi and 2phic
        beta = gsl_rng_uniform(r);
        if(beta > 0.5)
        {
            paramy[q][4] += PI/2.0;
            paramy[q][9] += PI/2.0;
        }
              
        
        free_double_matrix(Fish,4);
        free(pmap);
        free(params);
        free_double_matrix(Svec,4);
        free_double_vector(Sval);
              
          }
             
         }
        
         
         // Note that testing this proposal using logL=const requires is not very informative since the
         // Ftstat likelihood mostly gets skipped since the intrinsic parameters are not close to the true.

        
        
    }
    
    // [0] ln Mass1  [1] ln Mass2  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln distance
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination

    
    if(paramy[q][1] > paramy[q][0])  // catch if m2 > m1 and flip
    {
        lm1 = paramy[q][1];
        chi1 = paramy[q][3];
        lm2 = paramy[q][0];
        chi2 = paramy[q][2];
        paramy[q][0] = lm1;
        paramy[q][1] = lm2;
        paramy[q][2] = chi1;
        paramy[q][3] = chi2;
    }
    
    
    cv[typ][k]++;
    
    // [0] ln Mass1  [1] ln Mass2  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln distance
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    

    // re-map angular parameters to their proper range
    if(paramy[q][4] > PI)   paramy[q][4] -= PI;
    if(paramy[q][4] < 0.0)  paramy[q][4] += PI;
    if(paramy[q][8] > TPI)  paramy[q][8] -= TPI;
    if(paramy[q][8] < 0.0)  paramy[q][8] += TPI;
    if(paramy[q][9] > PI)   paramy[q][9] -= PI;
    if(paramy[q][9] < 0.0)  paramy[q][9] += PI;

    
    // check proposed values are in prior range
    flag = 0;
    for(i = 0; i < NP; i++)
    {
        if(paramy[q][i] > max[i] || paramy[q][i] < min[i]) flag = 1;
    }
    
    if(flag == 0)
    {
        
        // Jacobian that makes the prior flat in m1, m2.
        // Jacobian is m1*m2, but we combine the probablities as logs
        pMcMy = paramy[q][0]+paramy[q][1];
        pMcMx = paramx[q][0]+paramx[q][1];
        
        logLy = 0.0;
        nhy = 0.0;
        
        if(lhold == 0)
        {
            if(fstat == 1)
            {
            logLy = LikelihoodFstat(ll, Tend, paramy[q], tx, NF, FF, AA, EA, AP, EP, SN);
        
                    // To cover the full range, apply a pi shift to 2psi and 2phic
                beta = gsl_rng_uniform(r);
                if(beta > 0.5)
                {
                    paramy[q][4] += PI/2.0;
                    paramy[q][9] += PI/2.0;
                }

                if(paramy[q][4] > PI)   paramy[q][4] -= PI;
                if(paramy[q][4] < 0.0)  paramy[q][4] += PI;
                if(paramy[q][9] > PI)   paramy[q][9] -= PI;
                if(paramy[q][9] < 0.0)  paramy[q][9] += PI;
                
                /*
                x = LikelihoodDelta(ll, Tend, paramy[q], NF, FF, AA, EA, AP, EP, SN);
                printf("%f %f %f\n", logLy-x, logLy, logLx);
                */
                
            }
            else
            {
              // this is the (n|dh) term in the likelihood
             if(nflag == 1)
             {
                  nhy = -ndh(ll, Tend, paramy[q], ATR, ATI, ETR, ETI, NF, FF, AA, EA, AP, EP, FN, ASD);
                 if(fabs(nhy) > (double)(NP)) nhy = -1.0e8; // for this to happen something must have gone wrong with the splines
                //nhy = nh(ll, Tend, paramy[q], ATR, ATI, ETR, ETI, NF, FF, AP, EP, FN, ASD);
                logLy = LikelihoodDelta(ll, Tend, paramy[q], NF, FF, AA, EA, AP, EP, SN) + nhy;
             }
             else
             {
                logLy = LikelihoodDelta(ll, Tend, paramy[q], NF, FF, AA, EA, AP, EP, SN);
             }
               // if(typ==2) printf("%f\n", logLy);
            }
            
        }
        
    
    // variable in MCMC is x=logD, so p(x) = dD/dx p(D) = D p(D) = D^3
    
    //pDx = 3.0*paramx[q][6];   // uniform in volume prior
    //pDy = 3.0*paramy[q][6];   // uniform in volume prior
    
    pDx = paramx[q][6];   // uniform in distance prior
    pDy = paramy[q][6];   // uniform in distance prior
    
    
    H = (logLy-logLx[q])/heat[k] + pMcMy + pDy - qy - pDx - pMcMx + qx;
    
    
    alpha = log(gsl_rng_uniform(r));

    
    if(H > alpha)
    {
        // copy over new state if accepted
        logLx[q] = logLy;
        nhx[q] = nhy;
        for(i = 0; i < NP; i++) paramx[q][i] = paramy[q][i];
        av[typ][k]++;
    }
        
    }
    
    free_double_vector(u);
    free_double_vector(jump);
    free_double_vector(zv);
    
}

void de_jump(double *paramsx, double *paramsy, double **history, int m, int d, gsl_rng *r)
{
    int i, j, k;
    double alpha, beta;
    
    // pick two points from the history
    i = (int)((double)(m)*gsl_rng_uniform(r));
    k = 0;
    do
    {
        j = (int)((double)(m)*gsl_rng_uniform(r));
        k++;
    } while(i==j);
    
   // printf("%d\n", k);
    
    alpha = 1.0;
    beta = gsl_rng_uniform(r);
    if(beta < 0.9) alpha = gsl_ran_gaussian(r,0.5);
    
    beta = gsl_rng_uniform(r);
    if(beta > 0.5)
    {  // all parameters
    for(k=0; k< d; k++) paramsy[k] = paramsx[k]+alpha*(history[i][k]-history[j][k]);
    }
    else
    { // just the intrinsic (plus phase)
      for(k=0; k< NI; k++) paramsy[k] = paramsx[k]+alpha*(history[i][k]-history[j][k]);
      for(k=NI; k< d; k++) paramsy[k] = paramsx[k];
    }
    
}

void FisherEvec(double **fish, double *ej, double **ev, int d)
{
    int i, j, ecc, sc;
    double x, maxc;
 
    ecc = 0;
    for (i = 0 ; i < d ; i++) if(fabs(fish[i][i]) < 1.0e-16) ecc = 1;
    
    if(ecc == 0)
    {
        
        gsl_matrix *m = gsl_matrix_alloc (d, d);
        
        for (i = 0 ; i < d ; i++)
        {
            for (j = 0 ; j < d ; j++)
            {
                gsl_matrix_set(m, i, j, fish[i][j]);
            }
        }
        
        
        gsl_vector *eval = gsl_vector_alloc (d);
        gsl_matrix *evec = gsl_matrix_alloc (d, d);
        
        gsl_eigen_symmv_workspace * w =
        gsl_eigen_symmv_alloc (d);
        
        ecc = gsl_eigen_symmv (m, eval, evec, w);
        
        gsl_eigen_symmv_free (w);
        
        
        sc = gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
        
        for (i = 0; i < d; i++)
        {
            ej[i] = gsl_vector_get (eval, i);
            
           // printf("eigenvalue = %g\n", ej[i]);
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = gsl_matrix_get(evec, j, i);
               // printf("%f ", ev[i][j]);
            }
           // printf("\n");
            
        }
        
        for (i = 0; i < d; i++)
        {
            // make sure no eigenvalue is too small
            //if(ej[i] < 100.0) ej[i] = 100.0;
            // turn into 1-sigma jump amplitudes
            ej[i] = 1.0/sqrt(fabs(ej[i]));
            //printf("jump %d = %g\n", i, ej[i]);
        }
        
        gsl_matrix_free (m);
        gsl_vector_free (eval);
        gsl_matrix_free (evec);
        
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            ej[i] = 10000.0;
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
                if(i==j) ev[i][j] = 1.0;
            }
        }
        
    }
    
    return;
    
}

void FisherEvecSVD(double **fish, double *ej, double **ev, int d)
{
    int i, j, ecc, sc;
    double x, maxc;
 
    ecc = 0;
    for (i = 0 ; i < d ; i++) if(fabs(fish[i][i]) < 1.0e-16) ecc = 1;
    
    if(ecc == 0)
    {
        
        gsl_matrix *m = gsl_matrix_alloc (d, d);
        
        for (i = 0 ; i < d ; i++)
        {
            for (j = 0 ; j < d ; j++)
            {
                gsl_matrix_set(m, i, j, fish[i][j]);
            }
        }
        
        
        gsl_vector *S = gsl_vector_alloc (d);
        gsl_matrix *V = gsl_matrix_alloc (d, d);
        
        gsl_linalg_SV_decomp_jacobi (m, V, S);
        
        
        for (i = 0; i < d; i++)
        {
            ej[i] = gsl_vector_get (S, i);
            
           // printf("eigenvalue = %g\n", ej[i]);
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = gsl_matrix_get(m, j, i);
               // printf("%f ", ev[i][j]);
            }
           // printf("\n");
            
        }
        
        for (i = 0; i < d; i++)
        {
            // make sure no eigenvalue is too small
            //if(ej[i] < 100.0) ej[i] = 100.0;
            // turn into 1-sigma jump amplitudes
            //ej[i] = 1.0/sqrt(fabs(ej[i]));
            //printf("jump %d = %g\n", i, ej[i]);
        }
        
        gsl_matrix_free (m);
        gsl_vector_free (S);
        gsl_matrix_free (V);
        
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            ej[i] = 10000.0;
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
                if(i==j) ev[i][j] = 1.0;
            }
        }
        
    }
    
    return;
    
}

void FisherEvecSplit(double **fish, double *ej, double **ev, int d)
{
    int i, j, ecc, sc;
    double x, maxc;
    int NIN = 6;
    int NEX = 7;
    double **cov;
    
    
    ecc = 0;
    for (i = 0 ; i < d ; i++) if(fabs(fish[i][i]) < 1.0e-16) ecc = 1;
    
    
    
    if(ecc == 0)
    {
        
        for (i = 0; i < d; i++)
        {
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
            }
        }
        
        cov = double_matrix(d,d);
        
        Inverse(fish,cov,d);
        
        // start with intrinsic
        
        gsl_matrix *m = gsl_matrix_alloc (NIN, NIN);
        
        for (i = 0 ; i < NIN ; i++)
        {
            for (j = 0 ; j < NIN ; j++)
            {
                gsl_matrix_set(m, i, j, cov[i][j]);
            }
        }
        
        
        gsl_vector *eval = gsl_vector_alloc (NIN);
        gsl_matrix *evec = gsl_matrix_alloc (NIN, NIN);
        
        gsl_eigen_symmv_workspace * w =
        gsl_eigen_symmv_alloc (NIN);
        
        ecc = gsl_eigen_symmv (m, eval, evec, w);
        
        gsl_eigen_symmv_free (w);
        
        
        sc = gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
        
        for (i = 0; i < NIN-1; i++)
        {
            ej[i] = gsl_vector_get (eval, i);
            
            // printf("eigenvalue = %g\n", ej[i]);
            for (j = 0 ; j < NIN ; j++)
            {
                ev[i][j] = gsl_matrix_get(evec, j, i);
                // printf("%f ", ev[i][j]);
            }
            // printf("\n");
            
        }
        
        for (i = 0; i < NIN-1; i++)
        {
            // make sure no eigenvalue is too small
            //if(ej[i] < 100.0) ej[i] = 100.0;
            // turn into 1-sigma jump amplitudes
            ej[i] = sqrt(fabs(ej[i]));
            //printf("jump %d = %g\n", i, ej[i]);
        }
        
        gsl_matrix_free (m);
        gsl_vector_free (eval);
        gsl_matrix_free (evec);
        
        // now the extrinsic
        
        gsl_matrix *mx = gsl_matrix_alloc (NEX, NEX);
        
        for (i = 0 ; i < NEX ; i++)
        {
            for (j = 0 ; j < NEX ; j++)
            {
                gsl_matrix_set(mx, i, j, cov[i+(NP-NEX)][j+(NP-NEX)]);
            }
        }
        
        
        gsl_vector *evalx = gsl_vector_alloc (NEX);
        gsl_matrix *evecx = gsl_matrix_alloc (NEX, NEX);
        
        gsl_eigen_symmv_workspace * wx =
        gsl_eigen_symmv_alloc (NEX);
        
        ecc = gsl_eigen_symmv (mx, evalx, evecx, wx);
        
        gsl_eigen_symmv_free (wx);
        
        
        sc = gsl_eigen_symmv_sort (evalx, evecx, GSL_EIGEN_SORT_ABS_ASC);
        
        for (i = 0; i < NEX-1; i++)
        {
            ej[i+NIN-1] = gsl_vector_get (evalx, i);
            
            // printf("eigenvalue = %g\n", ej[i]);
            for (j = 0 ; j < NEX ; j++)
            {
                ev[i+NIN-1][j+NP-NEX] = gsl_matrix_get(evecx, j, i);
                // printf("%f ", ev[i][j]);
            }
            // printf("\n");
            
        }
        
        for (i = 0; i < NEX; i++)
        {
            // make sure no eigenvalue is too small
            //if(ej[i] < 100.0) ej[i] = 100.0;
            // turn into 1-sigma jump amplitudes
            ej[i+NIN-1] = sqrt(fabs(ej[i]));
            //printf("jump %d = %g\n", i, ej[i]);
        }
        
        gsl_matrix_free (mx);
        gsl_vector_free (evalx);
        gsl_matrix_free (evecx);

        free_double_matrix(cov,d);
        
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            ej[i] = 10000.0;
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
                if(i==j) ev[i][j] = 1.0;
            }
        }
        
    }
    
    return;
    
}



void Inverse(double **M, double **IM, int d)
{
    int i, j;
    int s;
    double x, maxc;
    
    gsl_matrix *m = gsl_matrix_alloc (d, d);
    
    for (i = 0 ; i < d ; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            gsl_matrix_set(m, i, j, M[i][j]);
        }
    }
    
    gsl_permutation *p = gsl_permutation_alloc(d);
    
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);
    
    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(d, d);
    
    gsl_linalg_LU_invert(m, p, inv);
    
    gsl_permutation_free(p);
    
    
    for (i = 0; i < d; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            IM[i][j] = gsl_matrix_get(inv, j, i);
        }
    }
    
    
    gsl_matrix_free (inv);
    gsl_matrix_free (m);
    
    return;
    
}



void FisherDirect(int ll, double Tend, double *params, double **Fisher, double *SN)
{
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW, lx;
    double A, P;
    double px, fnew, kxm;
    double *epsilon;
    double **AH, **AP, **AM;
    double **EH, **EP, **EM;
    
    int kmin, kmax, flag;

    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    AH = double_matrix(NP,NFFT);
    AP = double_matrix(NP,NFFT);
    AM = double_matrix(NP,NFFT);
    EH = double_matrix(NP,NFFT);
    EP = double_matrix(NP,NFFT);
    EM = double_matrix(NP,NFFT);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-4;
    epsilon[5] = 1.0e-2;
    epsilon[6] = 1.0e-4;
    epsilon[7] = 1.0e-4;
    epsilon[8] = 1.0e-4;
    epsilon[9] = 1.0e-4;
    epsilon[10] = 1.0e-4;
    
    for (i=0; i< NP; i++)
    {
        
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        
        if(ll==0)
        {
            if(i==0 || i==1 || i==6) // log derivatives
            {
            x = exp(epsilon[i]);
            paramsP[i] *= x;
            paramsM[i] /= x;
            }
            else
            {
                paramsP[i] += epsilon[i];
                paramsM[i] -= epsilon[i];
            }
            
        }
        else
        {
            paramsP[i] += epsilon[i];
            paramsM[i] -= epsilon[i];
        }

        ResponseFreq(ll, Tend, paramsP, NFFT, AP[i], EP[i]);
        ResponseFreq(ll, Tend, paramsM, NFFT, AM[i], EM[i]);

        for (j=0; j< NFFT; j++) AH[i][j] = (AP[i][j]-AM[i][j])/(2.0*epsilon[i]);
        for (j=0; j< NFFT; j++) EH[i][j] = (EP[i][j]-EM[i][j])/(2.0*epsilon[i]);
        
    }

       free_double_matrix(AP,NP);
       free_double_matrix(AM,NP);
       free_double_matrix(EP,NP);
       free_double_matrix(EM,NP);
    
       free(epsilon);
       free(paramsP);
       free(paramsM);

   
    for (i=0; i< NP; i++)
    {
        for (j=i; j< NP; j++)
        {
  
            Fisher[i][j] =  fourier_nwip(AH[i], AH[j], SN, NFFT)/Tobs + fourier_nwip(EH[i], EH[j], SN, NFFT)/Tobs;

        }
        
    }
    
    for (i=0; i< NP; i++)
    {
        for (j=i+1; j< NP; j++)
        {
            Fisher[j][i] =  Fisher[i][j];
        }
    }
    
    free_double_matrix(AH,NP);
    free_double_matrix(EH,NP);
    
    
}

void FisherDirectShift(int ll, double Tend, double *params, double **Fisher, double *SN)
{
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW, lx;
    double A, P;
    double px, fnew, kxm;
    double *epsilon;
    double **AH, **AP, **AM;
    double **EH, **EP, **EM;
    
    int kmin, kmax, flag;

    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    AH = double_matrix(NP,NFFT);
    AP = double_matrix(NP,NFFT);
    AM = double_matrix(NP,NFFT);
    EH = double_matrix(NP,NFFT);
    EP = double_matrix(NP,NFFT);
    EM = double_matrix(NP,NFFT);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-4;
    epsilon[5] = 1.0e-2;
    epsilon[6] = 1.0e-4;
    epsilon[7] = 1.0e-4;
    epsilon[8] = 1.0e-4;
    epsilon[9] = 1.0e-4;
    epsilon[10] = 1.0e-4;
    
    for (i=0; i< NP; i++)
    {
        
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        
        if(ll==0)
        {
            if(i==0 || i==1 || i==6) // log derivatives
            {
            x = exp(epsilon[i]);
            paramsP[i] *= x;
            paramsM[i] /= x;
            }
            else
            {
                paramsP[i] += epsilon[i];
                paramsM[i] -= epsilon[i];
            }
            
        }
        else
        {
            paramsP[i] += epsilon[i];
            paramsM[i] -= epsilon[i];
        }

        ResponseFreq(ll, Tend, paramsP, NFFT, AP[i], EP[i]);
        ResponseFreq(ll, Tend, paramsM, NFFT, AM[i], EM[i]);

        for (j=0; j< NFFT; j++) AH[i][j] = (AP[i][j]-AM[i][j])/(2.0*epsilon[i]);
        for (j=0; j< NFFT; j++) EH[i][j] = (EP[i][j]-EM[i][j])/(2.0*epsilon[i]);
        
    }

       free_double_matrix(AP,NP);
       free_double_matrix(AM,NP);
       free_double_matrix(EP,NP);
       free_double_matrix(EM,NP);
    
       free(epsilon);
       free(paramsP);
       free(paramsM);

   
    for (i=0; i< NP; i++)
    {
        for (j=i; j< NP; j++)
        {
  
            Fisher[i][j] =  fourier_nwip(AH[i], AH[j], SN, NFFT)/Tobs + fourier_nwip(EH[i], EH[j], SN, NFFT)/Tobs;

        }
        
    }
    
    for (i=0; i< NP; i++)
    {
        for (j=i+1; j< NP; j++)
        {
            Fisher[j][i] =  Fisher[i][j];
        }
    }
    
    free_double_matrix(AH,NP);
    free_double_matrix(EH,NP);
    
    
}

void FisherFast(int ll, double Tend, double *params, double **Fisher)
{
    double *AAmp, *EAmp, *APhase, *EPhase, *TF, *FF, *PF, *AF;
    double *PFref, *AFref, *TFref;
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW, lx;
    double A, P;
    double px, fnew, kxm;
    double *Aint, *Eint, *SAE;
    double *AphaseP, *AphaseM, *AampP, *AampM;
    double *EphaseP, *EphaseM, *EampP, *EampM;
    double **DphaseA, **DampA;
    double **DphaseE, **DampE;
    double *epsilon;
    
    int NFmax = 5000;
    
    int kmin, kmax, flag;

    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FF = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
    
    // Note: This call sets the FF array, which then gets held and used even when the parameters change a little
    SetUp(ll, params, NFmax, &NF, FF, TF, PF, AF);
    
    TFref = (double*)malloc(sizeof(double)* (NF));
    PFref = (double*)malloc(sizeof(double)* (NF));
    AFref = (double*)malloc(sizeof(double)* (NF));
    
    for (i=0; i< NF; i++)
     {
      TFref[i] = TF[i];
      PFref[i] = PF[i];
      AFref[i] = AF[i];
      //printf("%d %e\n", i, FF[i]);
     }
    
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-4;
    epsilon[5] = 1.0e-1;
    epsilon[6] = 1.0e-4;
    epsilon[7] = 1.0e-4;
    epsilon[8] = 1.0e-4;
    epsilon[9] = 1.0e-4;
    epsilon[10] = 1.0e-4;
    
    AphaseP = (double*)malloc(sizeof(double)* (NF));
    AphaseM = (double*)malloc(sizeof(double)* (NF));
    
    AampP = (double*)malloc(sizeof(double)* (NF));
    AampM = (double*)malloc(sizeof(double)* (NF));
    
    EphaseP = (double*)malloc(sizeof(double)* (NF));
    EphaseM = (double*)malloc(sizeof(double)* (NF));
    
    EampP = (double*)malloc(sizeof(double)* (NF));
    EampM = (double*)malloc(sizeof(double)* (NF));
    
    DphaseA = double_matrix(NP, NF);
    DampA = double_matrix(NP, NF);
    
    DphaseE = double_matrix(NP, NF);
    DampE = double_matrix(NP, NF);
 
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    
    Aint = (double*)malloc(sizeof(double)* (NF));
    Eint = (double*)malloc(sizeof(double)* (NF));
    SAE = (double*)malloc(sizeof(double)* (NF));
    
    // reference amplitude and phase
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    // masses and spins
    for (i=0; i< 4; i++)
    {
        
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        
        if(ll==0)
        {
            if(i==0 || i==1) // log derivatives
            {
            x = exp(epsilon[i]);
            paramsP[i] *= x;
            paramsM[i] /= x;
            }
            else
            {
                paramsP[i] += epsilon[i];
                paramsM[i] -= epsilon[i];
            }
            
        }
        else
        {
            paramsP[i] += epsilon[i];
            paramsM[i] -= epsilon[i];
        }

            Intrinsic(ll,paramsP, NF, FF, TF, PF, AF);  // intrinsic
            Extrinsic(paramsP, Tend, NF, FF, TF, PF, AF, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
                
            Intrinsic(ll,paramsM, NF, FF, TF, PF, AF); // intrinsic
            Extrinsic(paramsM, Tend, NF, FF, TF, PF, AF, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic

            
            for (k=0; k< NF; k++)
            {
                DphaseA[i][k] = (AphaseP[k]-AphaseM[k]);
                DphaseE[i][k] = (EphaseP[k]-EphaseM[k]);
                DampA[i][k] = (AampP[k]-AampM[k]);
                DampE[i][k] = (EampP[k]-EampM[k]);
            }
        
    }

    

   for (i=4; i< 7; i++)
    {
    
      if(i == 4)  // Put in coallesence phase manually. Note this is orbital phase, so times two in GW phase
       {
           
        x = -4.0*epsilon[4];
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = x;
            DphaseE[i][k] = x;
            DampA[i][k] = 0.0;
            DampE[i][k] = 0.0;
        }
        
      }
     else if(i == 5)  // tc is not included in the phase subroutine
      {
        for (k=0; k< NF; k++)
        {
            x = -4.0*PI*FF[k]*epsilon[5];
            DphaseA[i][k] = x;
            DphaseE[i][k] = x;
            DampA[i][k] = 0.0;
            DampE[i][k] = 0.0;
        }
        
       }
      else if(i == 6)  // ln DL derivative
       {
           for (k=0; k< NF; k++)
           {
               DphaseA[i][k] = 0;
               DphaseE[i][k] = 0;
               DampA[i][k] = -2.0*epsilon[6]*AAmp[k];
               DampE[i][k] = -2.0*epsilon[6]*EAmp[k];
           }

       }

  }

    for (i=7; i< NP; i++)
    {
    
       for (j=0; j< NP; j++)
        {
        paramsP[j] = params[j];
        paramsM[j] = params[j];
        }
    
        paramsP[i] += epsilon[i];
        paramsM[i] -= epsilon[i];

        Extrinsic(paramsP, Tend, NF, FF, TFref, PFref, AFref, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
        Extrinsic(paramsM, Tend, NF, FF, TFref, PFref, AFref, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic
        
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = (AphaseP[k]-AphaseM[k]);
            DphaseE[i][k] = (EphaseP[k]-EphaseM[k]);
            DampA[i][k] = (AampP[k]-AampM[k]);
            DampE[i][k] = (EampP[k]-EampM[k]);
        }
        
    }

    for (k=0; k< NF; k++) instrument_noise(FF[k], &SAE[k], &x);


    gsl_interp_accel *IAacc = gsl_interp_accel_alloc();
    gsl_spline *IAspline = gsl_spline_alloc (gsl_interp_cspline, NF);


    gsl_interp_accel *IEacc = gsl_interp_accel_alloc();
    gsl_spline *IEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    
    
   
    for (i=0; i< NP; i++)
    {
        for (j=i; j< NP; j++)
        {
  
            for (k=0; k< NF; k++)
            {
            Aint[k] = 4.0*(DampA[i][k]*DampA[j][k]+AAmp[k]*AAmp[k]*DphaseA[i][k]*DphaseA[j][k])/SAE[k];
            Eint[k] = 4.0*(DampE[i][k]*DampE[j][k]+EAmp[k]*EAmp[k]*DphaseE[i][k]*DphaseE[j][k])/SAE[k];
            }
            
            gsl_spline_init(IAspline, FF, Aint, NF);
            gsl_spline_init(IEspline, FF, Eint, NF);
            
            Fisher[i][j] = gsl_spline_eval_integ(IAspline, FF[0], FF[NF-1], IAacc)+gsl_spline_eval_integ(IEspline, FF[0], FF[NF-1], IEacc);
            Fisher[i][j] /= (4.0*epsilon[i]*epsilon[j]);
        }
        
    }
    
    for (i=0; i< NP; i++)
    {
        for (j=i+1; j< NP; j++)
        {
            Fisher[j][i] =  Fisher[i][j];
        }
    }
    
    
    free_double_matrix(DphaseA,NP);
    free_double_matrix(DampA,NP);
    free_double_matrix(DphaseE,NP);
    free_double_matrix(DampE,NP);
 
    free(epsilon);
    free(SAE);
    free(Aint);
    free(Eint);
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(FF);
    free(TF);
    free(AF);
    free(PF);
    free(TFref);
    free(AFref);
    free(PFref);
    
    free(paramsP);
    free(paramsM);
    
    free(AphaseP);
    free(AphaseM);
    free(AampP);
    free(AampM);
    
    free(EphaseP);
    free(EphaseM);
    free(EampP);
    free(EampM);

    gsl_spline_free(IAspline);
    gsl_spline_free(IEspline);

    gsl_interp_accel_free(IAacc);
    gsl_interp_accel_free(IEacc);
 
    
    
}

void FisherSub(int ll, double Tend, int *pmap, double *params, double **Fisher)
{
    double *AAmp, *EAmp, *APhase, *EPhase, *TF, *FF, *PF, *AF;
    double *PFref, *AFref, *TFref;
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW;
    double A, P;
    double px, fnew, kxm;
    double *Aint, *Eint, *SAE;
    double *AphaseP, *AphaseM, *AampP, *AampM;
    double *EphaseP, *EphaseM, *EampP, *EampM;
    double **DphaseA, **DampA;
    double **DphaseE, **DampE;
    double *epsilon;
    
    int NFmax = 5000;
    
    int kmin, kmax, flag;
    
    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FF = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
    
    // Note: This call sets the FF array, which then gets held and used even when the parameters change a little
    SetUp(ll, params, NFmax, &NF, FF, TF, PF, AF);
    
    TFref = (double*)malloc(sizeof(double)* (NF));
    PFref = (double*)malloc(sizeof(double)* (NF));
    AFref = (double*)malloc(sizeof(double)* (NF));
    
    for (i=0; i< NF; i++)
    {
        TFref[i] = TF[i];
        PFref[i] = PF[i];
        AFref[i] = AF[i];
    }
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-4;
    epsilon[5] = 1.0e-1;
    epsilon[6] = 1.0e-4;
    epsilon[7] = 1.0e-4;
    epsilon[8] = 1.0e-4;
    epsilon[9] = 1.0e-4;
    epsilon[10] = 1.0e-4;
    
    AphaseP = (double*)malloc(sizeof(double)* (NF));
    AphaseM = (double*)malloc(sizeof(double)* (NF));
    
    AampP = (double*)malloc(sizeof(double)* (NF));
    AampM = (double*)malloc(sizeof(double)* (NF));
    
    EphaseP = (double*)malloc(sizeof(double)* (NF));
    EphaseM = (double*)malloc(sizeof(double)* (NF));
    
    EampP = (double*)malloc(sizeof(double)* (NF));
    EampM = (double*)malloc(sizeof(double)* (NF));
    
    DphaseA = double_matrix(NP, NF);
    DampA = double_matrix(NP, NF);
    
    DphaseE = double_matrix(NP, NF);
    DampE = double_matrix(NP, NF);
    
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    
    Aint = (double*)malloc(sizeof(double)* (NF));
    Eint = (double*)malloc(sizeof(double)* (NF));
    SAE = (double*)malloc(sizeof(double)* (NF));
    
    // reference amplitude and phase
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    // masses and spins
    for (i=0; i< 4; i++)
    {
        if(pmap[i] > -1)
        {
        
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        if(ll==0)
        {
            if(i==0 || i==1) // log derivatives
            {
                x = exp(epsilon[i]);
                paramsP[i] *= x;
                paramsM[i] /= x;
            }
            else
            {
                paramsP[i] += epsilon[i];
                paramsM[i] -= epsilon[i];
            }
            
        }
        else
        {
            paramsP[i] += epsilon[i];
            paramsM[i] -= epsilon[i];
        }
        
        
        
        Intrinsic(ll,paramsP, NF, FF, TF, PF, AF);  // intrinsic
        Extrinsic(paramsP, Tend, NF, FF, TF, PF, AF, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
        
        Intrinsic(ll,paramsM, NF, FF, TF, PF, AF); // intrinsic
        Extrinsic(paramsM, Tend, NF, FF, TF, PF, AF, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic
        
        
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = (AphaseP[k]-AphaseM[k]);
            DphaseE[i][k] = (EphaseP[k]-EphaseM[k]);
            DampA[i][k] = (AampP[k]-AampM[k]);
            DampE[i][k] = (EampP[k]-EampM[k]);
        }
            
       }
        
    }
    
    
    
    for (i=4; i< 7; i++)
    {
        if(pmap[i] > -1)
        {
            
        if(i == 4)  // Put in coallesence phase manually. Note this is orbital phase, so times two in GW phase
        {
            
            x = -4.0*epsilon[4];
            for (k=0; k< NF; k++)
            {
                DphaseA[i][k] = x;
                DphaseE[i][k] = x;
                DampA[i][k] = 0.0;
                DampE[i][k] = 0.0;
            }
            
        }
        else if(i == 5)  // tc is not included in the phase subroutine
        {
            for (k=0; k< NF; k++)
            {
                x = -4.0*PI*FF[k]*epsilon[5];
                DphaseA[i][k] = x;
                DphaseE[i][k] = x;
                DampA[i][k] = 0.0;
                DampE[i][k] = 0.0;
            }
            
        }
        else if(i == 6)  // ln DL derivative
        {
            for (k=0; k< NF; k++)
            {
                DphaseA[i][k] = 0;
                DphaseE[i][k] = 0;
                DampA[i][k] = -2.0*epsilon[6]*AAmp[k];
                DampE[i][k] = -2.0*epsilon[6]*EAmp[k];
            }
            
        }
        
        }
    }
    
    for (i=7; i< NP; i++)
    {
        if(pmap[i] > -1)
        {
            
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        paramsP[i] += epsilon[i];
        paramsM[i] -= epsilon[i];
        
        Extrinsic(paramsP, Tend, NF, FF, TFref, PFref, AFref, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
        Extrinsic(paramsM, Tend, NF, FF, TFref, PFref, AFref, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic
        
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = (AphaseP[k]-AphaseM[k]);
            DphaseE[i][k] = (EphaseP[k]-EphaseM[k]);
            DampA[i][k] = (AampP[k]-AampM[k]);
            DampE[i][k] = (EampP[k]-EampM[k]);
        }
            
        }
        
    }
    
    for (k=0; k< NF; k++) instrument_noise(FF[k], &SAE[k], &x);
    
    
    gsl_interp_accel *IAacc = gsl_interp_accel_alloc();
    gsl_spline *IAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    
    
    gsl_interp_accel *IEacc = gsl_interp_accel_alloc();
    gsl_spline *IEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    
    
    
    for (i=0; i< NP; i++)
    {
        if(pmap[i] > -1)
        {
            
        for (j=i; j< NP; j++)
        {
            if(pmap[j] > -1)
            {
                
            for (k=0; k< NF; k++)
            {
                Aint[k] = 4.0*(DampA[i][k]*DampA[j][k]+AAmp[k]*AAmp[k]*DphaseA[i][k]*DphaseA[j][k])/SAE[k];
                Eint[k] = 4.0*(DampE[i][k]*DampE[j][k]+EAmp[k]*EAmp[k]*DphaseE[i][k]*DphaseE[j][k])/SAE[k];
            }
            
            gsl_spline_init(IAspline, FF, Aint, NF);
            gsl_spline_init(IEspline, FF, Eint, NF);
            
            Fisher[pmap[i]][pmap[j]] = gsl_spline_eval_integ(IAspline, FF[0], FF[NF-1], IAacc)+gsl_spline_eval_integ(IEspline, FF[0], FF[NF-1], IEacc);
            Fisher[pmap[i]][pmap[j]] /= (4.0*epsilon[i]*epsilon[j]);
                
            }
        }
            
        }
        
    }
    
    for (i=0; i< NP; i++)
    {
        if(pmap[i] > -1)
        {
            
        for (j=i+1; j< NP; j++)
        {
            if(pmap[j] > -1)
            {
            Fisher[pmap[j]][pmap[i]] =  Fisher[pmap[i]][pmap[j]];
            }
        }
            
        }
    }
    
    free_double_matrix(DphaseA,NP);
    free_double_matrix(DampA,NP);
    free_double_matrix(DphaseE,NP);
    free_double_matrix(DampE,NP);
    
    free(epsilon);
    free(SAE);
    free(Aint);
    free(Eint);
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(FF);
    free(TF);
    free(AF);
    free(PF);
    free(TFref);
    free(AFref);
    free(PFref);
    
    free(paramsP);
    free(paramsM);
    
    free(AphaseP);
    free(AphaseM);
    free(AampP);
    free(AampM);
    
    free(EphaseP);
    free(EphaseM);
    free(EampP);
    free(EampM);
    
    gsl_spline_free(IAspline);
    gsl_spline_free(IEspline);
    
    gsl_interp_accel_free(IAacc);
    gsl_interp_accel_free(IEacc);

    
}



void Intrinsic(int ll, double *params, int NF, double *FF, double *TF, double *PF, double *AF)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    RealVector *freq;
    
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double m1_SI, m2_SI, distance, tc, phic;
    
    int i, ii, ret, flag;
    double fonfs, t, told, tx;
    
    double fRef_in=PDfref;
    
    if(ll == 0)  // linear
    {
    m1 = params[0];    // Mass1
    m2 = params[1];    // Mass2
    distance = params[6]*1.0e9*PC_SI; // distance
    }
    else  // log
    {
     m1 = exp(params[0]);    // Mass1
     m2 = exp(params[1]);    // Mass2
     distance = exp(params[6])*1.0e9*PC_SI; // distance
    }
    
    chi1 = params[2];  // Spin1
    chi2 = params[3];  // Spin2
    m1_SI =  m1*MSUN_SI;
    m2_SI =  m2*MSUN_SI;
    tc = params[5];    // merger time
    
    Mtot = (m1+m2)*TSUN;
    eta = m1*m2/((m1+m2)*(m1+m2));
    
    freq = CreateRealVector(NF);
    
    for (i=0; i< NF; i++) freq->data[i] = FF[i];
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(
                                          &ap,      /**< [out] FD waveform */
                                          freq,    /**< Input: frequencies (Hz) on which to evaluate h22 FD */
                                          0.0,                  /**< Orbital phase at fRef (rad) */
                                          fRef_in,               /**< reference frequency (Hz) */
                                          m1_SI,                 /**< Mass of companion 1 (kg) */
                                          m2_SI,                 /**< Mass of companion 2 (kg) */
                                          chi1,                  /**< Aligned-spin parameter of companion 1 */
                                          chi2,                  /**< Aligned-spin parameter of companion 2 */
                                          distance               /**< Distance of source (m) */
                                          );
    
    
    flag = 0;
    told = ap->time[0]+tc;
    for (i=0; i< NF; i++)
    {
        PF[i] = ap->phase[i];
        
        AF[i] =  h22fac*ap->amp[i];
        fonfs = freq->data[i]/fstar;
        AF[i] *= (8.0*fonfs*sin(fonfs));   // conversion to fractional frequency and leading order TDI transfer function
        
        t = ap->time[i]+tc;
        if(t < told && flag == 0)
        {
            flag = 1;
            ii = i-1;
            tx = told;
        }
        TF[i] = t;
        if(t < -Tpad) TF[i] = -Tpad;
        if(flag == 1) TF[i] = tx +(double)(i-ii)*Mtot;
        told = t;
        //printf("%d %e %e\n", i, FF[i], TF[i]);
    }
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);

    
}


void SetUp(int ll, double *params, int NFmax, int *NFS, double *FF, double *TF, double *PF, double *AF)
{
    double af, fr, df, DT, fac, deltaF, f, fmax, fmin, x;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic, told;
    int i, ii, NF;
    double A, P;
    double px, fnew, fonfs;
    double t, tx;
    double dfmin, dfmax;
    
    double fRef_in=PDfref;
    
    int ret, flag, flag1, flag2;
    
    // This subroutine sets up the frequency sample array, the time-frequency map and finds the PhenomD amplitude and phase
    
    // NFmax is the size of the holder arrays. NFS is the actual size.
    
    if(ll == 0)  // linear
    {
        m1 = params[0];    // Mass1
        m2 = params[1];    // Mass2
        distance = params[6]*1.0e9*PC_SI; // distance
    }
    else  // log
    {
        m1 = exp(params[0]);    // Mass1
        m2 = exp(params[1]);    // Mass2
        distance = exp(params[6])*1.0e9*PC_SI; // distance
    }
    
    tc = params[5];    // merger time
    
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0)*TSUN;
    
    StartStop(ll, params, Tobs, 0.0, Tobs, &fmin, &fmax, &fr);
    
    //printf("%e %e %e\n", fmin, fmax, 1.0/Mtot);
    
    // this can happen when tc is really small and the masses are large
    if(fmax < fmin)
    {
        fmin = 0.5*fmax;
    }
    
    dfmin = 1.0/Tobs;
    dfmax = fr/100.0;
    DT = 3.0e5;
    
    fac = DT*pow(8.0*PI, 8.0/3.0)*3.0/40.0*pow(Mc,5.0/3.0);
    
    f = fmin;
    NF = 1;
    do
    {
        df = fac*pow(f,11.0/3.0);
        //printf("%e %e %e %e\n", f, df, dfmin, dfmax);
        if(df < dfmin) df = dfmin;
        if(df > dfmax) df = dfmax;
        f += df;
        NF++;
         //printf("%d %e\n", NF, f);
    }while(f < fmax);
    
    //printf("%d %e %e\n", NF, fmin, fmax);
    
     // Need to catch is NF > NFmax
    
    f = fmin;
    FF[0] = fmin;
    for (i=1; i< NF; i++)
    {
        df = fac*pow(f,11.0/3.0);
        if(df < dfmin) df = dfmin;
        if(df > dfmax) df = dfmax;
        f += df;
        FF[i] = f;
    }
    
    if(NF < 4)
    {
        NF = 4;
        df = (fmax-fmin)/3.0;
        FF[0] = fmin;
        for (i=1; i< NF; i++) FF[i] = fmin+df*(double)(i);
    }
    
    
    //for (i=0; i< NF; i++) printf("%d %e\n", i, FF[i]);
    
    Intrinsic(ll, params, NF, FF, TF, PF, AF);

    *NFS = NF;
    
}


void Hetrodyne(int ll, double Tend, double *params, long N, double *ATR, double *ATI, double *ETR, double *ETI, double *AD, double *ED, double *SN)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    double af, fr, df, DT, fac, deltaF, t, f, fmin;
    double x, y, u, v;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess, kxm, Tcut;
    double m1_SI, m2_SI, distance, tc, phic, fd, fe, ft;
    double cp, sp;
    int i, n, NF;
    double A, P;
    double px, fnew;
    double tau;
    
    int NFmax = 5000;
    
    double *AR, *ER;
    double *AI, *EI;
    
    double *AF, *PF, *FF, *TF;
    
    
    FILE *out;
    
    AR = (double*)malloc(sizeof(double)* (N));
    ER = (double*)malloc(sizeof(double)* (N));
    AI = (double*)malloc(sizeof(double)* (N));
    EI = (double*)malloc(sizeof(double)* (N));
    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FF = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
    
    SetUp(ll, params, NFmax, &NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    gsl_interp_accel *Tacc = gsl_interp_accel_alloc();
    gsl_spline *Tspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(Tspline, FF, TF, NF);
    
    gsl_interp_accel *PAacc = gsl_interp_accel_alloc();
    gsl_spline *PAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PAspline, FF, APhase, NF);
    
    gsl_interp_accel *PEacc = gsl_interp_accel_alloc();
    gsl_spline *PEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PEspline, FF, EPhase, NF);
    
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    
    for (i=0; i< N; i++)
    {
        AR[i] = 0.0;
        ER[i] = 0.0;
        AI[i] = 0.0;
        EI[i] = 0.0;
    }

    for (i=1; i< N/2; i++)
    {
        f = (double)(i)/Tobs;
        
        P = 0.0;
        
        if(f >= FF[0] && f <= FF[NF-1])
        {
            
            px = 2.0*PI*f*(Tobs-tc + dt/2.0)-2.0*phic;
            
            P = gsl_spline_eval (PAspline, f, PAacc);
            
            cp = cos(P+px);
            sp = sin(P+px);
            
            AR[i] = (AD[i]*cp+AD[N-i]*sp);
            AI[i] = (-AD[i]*sp+AD[N-i]*cp);
            
            P = gsl_spline_eval (PEspline, f, PEacc);
            
            cp = cos(P+px);
            sp = sin(P+px);
            
            ER[i] = (ED[i]*cp+ED[N-i]*sp);
            EI[i] = (-ED[i]*sp+ED[N-i]*cp);

        }
        
    }
    
    gsl_fft_real_radix2_transform(AR, 1, N);
    gsl_fft_real_radix2_transform(ER, 1, N);
    gsl_fft_real_radix2_transform(AI, 1, N);
    gsl_fft_real_radix2_transform(EI, 1, N);
    
    ATR[0] = AR[0];
    ETR[0] = ER[0];
    ATI[0] = AI[0];
    ETI[0] = EI[0];
    
    for (i=1; i< MF/2; i++)
    {
        ATR[i] = AR[i];
        ATR[MF-i] = AR[N-i];
        ETR[i] = ER[i];
        ETR[MF-i] = ER[N-i];
        ATI[i] = AI[i];
        ATI[MF-i] = AI[N-i];
        ETI[i] = EI[i];
        ETI[MF-i] = EI[N-i];
    }
    
    
    gsl_spline_free(PAspline);
    gsl_spline_free(PEspline);
    gsl_interp_accel_free(PAacc);
    gsl_interp_accel_free(PEacc);


    free(AR);
    free(ER);
    free(AI);
    free(EI);
    
    free(TF);
    free(FF);
    free(PF);
    free(AF);
    
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    
    
}





void Extrinsic(double *params, double Tend, int NF, double *FF, double *TF, double *PF, double *AF, double *AAmp, double *EAmp, double *APhase, double *EPhase, double *kxm)
{
    
    /*   Indicies   */
    int i,j, k, n, m;
    
    /*   Time and distance variables   */
    double *xi;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx, x;
    
    double *ta, *xia;
    double *FpAR, *FcAR, *FpER, *FcER;
    double *FpAI, *FcAI, *FpEI, *FcEI;
    
    double t, f, kdotx, Amp, Phase, RR, II;
    
    int NA, flag;
    
    double iota, cosi;
    double Aplus, Across;
    
    double Tcut;
    
    clock_t start, end;
    double cpu_time_used;
    
    FILE *out;
    
    
    xi = (double*)malloc(sizeof(double)* (NF));
    FpAR = (double*)malloc(sizeof(double)* (NF));
    FcAR = (double*)malloc(sizeof(double)* (NF));
    FpER = (double*)malloc(sizeof(double)* (NF));
    FcER = (double*)malloc(sizeof(double)* (NF));
    FpAI = (double*)malloc(sizeof(double)* (NF));
    FcAI = (double*)malloc(sizeof(double)* (NF));
    FpEI = (double*)malloc(sizeof(double)* (NF));
    FcEI = (double*)malloc(sizeof(double)* (NF));

    
    RAantenna(params, NF, TF, FF, xi, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
    for(n=0; n< NF; n++)
    {
        AAmp[n] = 0.0;
        EAmp[n] = 0.0;
        APhase[n] = 0.0;
        EPhase[n] = 0.0;
    }
    
    cosi = params[10];  // cos of inclination
    
    Aplus = 0.5*(1.+cosi*cosi);
    Across = -cosi;
    
    // Merger kdotx
    *kxm = (TF[NF-1]-xi[NF-1]);
    
    Tcut = Tend+t_rise/1.5;

    //out = fopen("map_fast.dat","w");
    for(n=0; n< NF; n++)
    {
        // Barycenter time and frequency
        t = TF[n];
        f = FF[n];
        
        x = 1.0;
        if(t > Tcut-t_rise && t < Tcut)
        {
        x = 0.5*(1.0-cos(PI*(t-Tcut)/t_rise));
        }
        if(t > Tcut) x = 0.0;
        
        kdotx = t-xi[n];
        
        Amp = x*AF[n];
        
        //  - FF[n]/fstar  Derivation says this is needed in the phase. Doesn't seem to be.
        
        Phase = -2.0*PI*f*kdotx-PF[n];
        
        RR = FpAR[n]*Aplus - FcAI[n]*Across;
        II = FcAR[n]*Across + FpAI[n]*Aplus;
        
        AAmp[n] = Amp*sqrt(RR*RR+II*II);
        APhase[n] = Phase+atan2(II,RR);
        
        RR = FpER[n]*Aplus - FcEI[n]*Across;
        II = FcER[n]*Across + FpEI[n]*Aplus;
        
        EAmp[n] = Amp*sqrt(RR*RR+II*II);
        EPhase[n] = Phase+atan2(II,RR);
        
        //printf("%d %e %e\n", n, t, f);
        
        //fprintf(out,"%d %e %e %e %e %e\n", n, t, f, kdotx, Amp, Phase);
        
    }
    //fclose(out);
    
    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("response loop took %f seconds\n", cpu_time_used);
    
    
    /*   Deallocate Arrays   */
    
    free(xi);
    
    free(FcAR);
    free(FpAR);
    free(FcER);
    free(FpER);
    free(FcAI);
    free(FpAI);
    free(FcEI);
    free(FpEI);
    
    return;
}



// Takes an exisitng intrinsic waveform
void FstatRA(int ll, double Tend, double *params, double *pnew, int NF, double *FF, double *TF, double *PF, double *AF, double *AAmp, double *EAmp, double *APhase, double *EPhase, double *SN)
{
    
    /*   Indicies   */
    int i,j, k, n, m, imax;
    
    /*   Time and distance variables   */
    double *xi;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx;
    
    double *ta, *xia, *PX;
    
    double *FpAR, *FcAR, *FpER, *FcER;
    double *FpAI, *FcAI, *FpEI, *FcEI;
    
    double **AS, **ES;
    
    double *AC, *EC;
    
    double **MM, **MI;
    
    double ***MMI;
    
    double *NV, *aV, **NVI;
    
    double *filtAR, *filtER;
    double *filtAI, *filtEI;
    
    double AR, AI, ER, EI;
    
    double t, f, kdotx, Amp, Phase, Fp, Fc;
    
    double A, P, px, tc, tcs, pxi, pxj;
    
    double x, y, u, v, xold, yold;
    
    double cp, sp;
    
    int NA, flag;
    
    double cx, sx, logL, logLmax;
    
    int nn;
    
    double iota, cosi;
    double Ap, Ac, scale;
    
    double a1, a2, a3, a4;
    double cpsi, spsi;
    double psi, phic;
    double cpc, spc;
    
    FILE *out;
    
    
    clock_t start, end;
    double cpu_time_used;
    
    filtAR = (double*)malloc(sizeof(double)* (4));
    filtER = (double*)malloc(sizeof(double)* (4));
    filtAI = (double*)malloc(sizeof(double)* (4));
    filtEI = (double*)malloc(sizeof(double)* (4));
    
    xi = (double*)malloc(sizeof(double)* (NF));
    PX = (double*)malloc(sizeof(double)* (NF));
    
    FpAR = (double*)malloc(sizeof(double)* (NF));
    FcAR = (double*)malloc(sizeof(double)* (NF));
    FpER = (double*)malloc(sizeof(double)* (NF));
    FcER = (double*)malloc(sizeof(double)* (NF));
    
    FpAI = (double*)malloc(sizeof(double)* (NF));
    FcAI = (double*)malloc(sizeof(double)* (NF));
    FpEI = (double*)malloc(sizeof(double)* (NF));
    FcEI = (double*)malloc(sizeof(double)* (NF));
    
    RAfilters(params, NF, TF, FF, xi, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
    //out = fopen("map_fast.dat","w");
    for(n=0; n< NF; n++)
    {
        // Barycenter time and frequency
        t = TF[n];
        f = FF[n];
        
        kdotx = t-xi[n];
        
        Amp = AF[n];
        
        // Detector phase (sans time and phase shift)
        PX[n]  = -2.0*PI*f*kdotx-PF[n];
        
        FpAR[n] *= Amp;
        FcAR[n] *= Amp;
        
        FpER[n] *= Amp;
        FcER[n] *= Amp;
        
        FpAI[n] *= Amp;
        FcAI[n] *= Amp;
        
        FpEI[n] *= Amp;
        FcEI[n] *= Amp;
        
    }
    
    tc = params[5];    // merger time
    
    NVI = double_matrix(4,NF);
    MMI = double_tensor(4,4,NF);
    
    nn = 0;
    
    px = PX[0] + TPI*FF[0]*(Tobs-tc + dt/2.0);
    xold = cos(px-APhase[0]);
    
    for(n=0; n< NF; n++)
    {
        px = PX[n] + TPI*FF[n]*(Tobs-tc + dt/2.0);
        
            cp = cos(px);
            sp = sin(px);
        
            AR = AAmp[n]*cos(APhase[n]);
            AI = AAmp[n]*sin(APhase[n]);
            ER = EAmp[n]*cos(EPhase[n]);
            EI = EAmp[n]*sin(EPhase[n]);
        
            // counts the windings
            x = cos(px-APhase[n]);
            if(x*xold < 0.0) nn++;
            xold = x;
        
           filtAR[0] = (FpAR[n]*cp-FpAI[n]*sp);
           filtER[0] = (FpER[n]*cp-FpEI[n]*sp);
        
           filtAR[1] = -(FcAR[n]*cp-FcAI[n]*sp);
           filtER[1] = -(FcER[n]*cp-FcEI[n]*sp);
        
           filtAR[2] = -(FpAI[n]*cp+FpAR[n]*sp);
           filtER[2] = -(FpEI[n]*cp+FpER[n]*sp);
        
           filtAR[3] = (FcAI[n]*cp+FcAR[n]*sp);
           filtER[3] = (FcEI[n]*cp+FcER[n]*sp);
        
           filtAI[0] = (FpAR[n]*sp+FpAI[n]*cp);
           filtEI[0] = (FpER[n]*sp+FpEI[n]*cp);
        
           filtAI[1] = -(FcAR[n]*sp+FcAI[n]*cp);
           filtEI[1] = -(FcER[n]*sp+FcEI[n]*cp);
        
           filtAI[2] = (-FpAI[n]*sp+FpAR[n]*cp);
           filtEI[2] = (-FpEI[n]*sp+FpER[n]*cp);
        
           filtAI[3] = -(-FcAI[n]*sp+FcAR[n]*cp);
           filtEI[3] = -(-FcEI[n]*sp+FcER[n]*cp);

        

          for (i=0; i< 4; i++)
           {
               NVI[i][n] = 4.0*(filtAR[i]*AR+filtAI[i]*AI+filtER[i]*ER+filtEI[i]*EI)/SN[n];
            for (j=i; j< 4; j++)
            {
                MMI[i][j][n] = 4.0*(filtAR[i]*filtAR[j]+filtAI[i]*filtAI[j]+filtER[i]*filtER[j]+filtEI[i]*filtEI[j])/SN[n];
            }
           }
        
    }
    
    if(nn < 20)
    {
        
        NV = (double*)malloc(sizeof(double)* (4));
        aV = (double*)malloc(sizeof(double)* (4));
        MM = double_matrix(4,4);
        MI = double_matrix(4,4);
        
        //printf("\n");
        for (i=0; i< 4; i++)
        {
            gsl_interp_accel *Iacc = gsl_interp_accel_alloc();
            gsl_spline *Ispline = gsl_spline_alloc (gsl_interp_cspline, NF);
            gsl_spline_init(Ispline, FF, NVI[i], NF);
            NV[i] =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
            gsl_spline_free(Ispline);
            gsl_interp_accel_free(Iacc);
           // printf("%e ", NV[i]);
        }
        //printf("\n");
        
       // printf("\n");
        for (i=0; i< 4; i++)
        {
            for (j=i; j< 4; j++)
            {
                gsl_interp_accel *Iacc = gsl_interp_accel_alloc();
                gsl_spline *Ispline = gsl_spline_alloc (gsl_interp_cspline, NF);
                gsl_spline_init(Ispline, FF, MMI[i][j], NF);
                MM[i][j] =  gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
               // printf("%e ", MM[i][j]);
                gsl_spline_free(Ispline);
                gsl_interp_accel_free(Iacc);
            }
            //printf("\n");
        }
       // printf("\n");
        
        
        
        for (i=0; i< 4; i++)
        {
            for (j=i; j< 4; j++)
            {
                MM[j][i] = MM[i][j];
            }
        }
        
        // a small stabilizer
        for (i=0; i< 4; i++) MM[i][i] += 0.1;
        
        //x = det(MM,4);
        
        //printf("%e\n", x);
        
        //printf("%e %e %e %e\n", MM[0][0], MM[1][0], MM[2][3], MM[1][2]);
        
        Inverse(MM, MI, 4);
        
        logL = 0.0;
        
        for (i=0; i< 4; i++)
        {
            aV[i] = 0.0;
            for (j=0; j< 4; j++)
            {
                aV[i] += MI[i][j]*NV[j];
                logL += 0.5*MI[i][j]*NV[i]*NV[j];
            }
        }
        
        //printf("%e %f\n", logL, sqrt(2.0*logL));
        
        //printf("%e %e %e %e\n", aV[0], aV[1], aV[2], aV[3]);
        
        x = (aV[0]+aV[3]);
        x = x*x;
        y = (aV[1]-aV[2]);
        y = y*y;
        u = (aV[0]-aV[3]);
        u = u*u;
        v = (aV[1]+aV[2]);
        v = v*v;
        
        Ap = sqrt(x+y)+sqrt(u+v);
        Ac = sqrt(x+y)-sqrt(u+v);
        A = Ap + sqrt(Ap*Ap-Ac*Ac);
        
        //  + + - +   need -pi/2 on both
        //  + + + +   ?
        //   + + + -  ?
        //  + + - - ?
        //  + - + +
        //  + - + -   // phic = psi_true, psi = phic_true + pi/2
        //  + - - +
        //  + - - -
        
        //  - + + +  ?
        //  - + - + ?
        //  - + - -  ?
        //  - + + -
        //  - - + +
        //  - - + -  // phic = phic_true + pi/2, psi = psi_true
        //  - - - +
        //  - - - - ?
        
        x = atan2((aV[1]-aV[2]),(aV[0]+aV[3]));
        
        y = atan2(-(aV[1]+aV[2]),(aV[0]-aV[3]));
        
        psi = 0.25*(y-x);
        while(psi < 0.0) psi += PI;
        while(psi > PI) psi -= PI;
        
        
        phic = 0.25*(x+y);
        while(phic < 0.0) phic += PI;
        while(phic > PI) phic -= PI;
        
        cosi = Ac/A;
        
        scale = 2.0/A;
        
       // printf("%f %f %f %f %f %f %f %f %f\n", scale, x, y, phic, psi, cosi, params[4], params[9], params[10]);
        
        pnew[0] = logL;
        pnew[1] = cosi;
        pnew[2] = scale;
        pnew[3] = psi;
        pnew[4] = phic;
        
        
        free(NV);
        free(aV);
        free_double_matrix(MM,4);
        free_double_matrix(MI,4);
        
    }
    else
    {
        pnew[0] = 0.0;
        pnew[1] = params[10];
        pnew[2] = 1.0;
        pnew[3] = params[9];
        pnew[4] = params[4];
    }
    
    /*   Deallocate Arrays   */
    

    
    
    free_double_matrix(NVI,4);
    free_double_tensor(MMI,4,4);
    
    free(xi);
    free(PX);

    free(filtAR);
    free(filtER);
    free(filtAI);
    free(filtEI);

    free(FcAR);
    free(FpAR);
    free(FcER);
    free(FpER);
    
    free(FcAI);
    free(FpAI);
    free(FcEI);
    free(FpEI);
    
    
    
    return;
}

// Uses the full data (slow)
void FstatFull(int ll, double Tend, double *params, double *pnew, int N, double *AC, double *EC, double *SN)
{
    
    /*   Indicies   */
    int i,j, k, n, m, imax;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx;
    
    double *pfilt;
    
    double **filtA, **filtE;
    
    double **MM, **MI;
    
    double *NV, *aV;
    
    double AR, AI, ER, EI;
    
    double t, f, kdotx, Amp, Phase, Fp, Fc;
    
    double A, P, px, tc, tcs, pxi, pxj;
    
    double x, y, u, v, xold, yold;
    
    double cp, sp;
    
    int NA, flag;
    
    double cx, sx, logL, logLmax;
    
    int nn;
    
    double iota, cosi;
    double Ap, Ac, scale;
    
    double a1, a2, a3, a4;
    double cpsi, spsi;
    double psi, phic;
    double cpc, spc;
    
    FILE *out;
    
    
    clock_t start, end;
    double cpu_time_used;
    
    pfilt = (double*)malloc(sizeof(double)* (NP));
    
    filtA = double_matrix(4,N);
    filtE = double_matrix(4,N);
    
    NV = (double*)malloc(sizeof(double)* (4));
    aV = (double*)malloc(sizeof(double)* (4));
    MM = double_matrix(4,4);
    MI = double_matrix(4,4);
    
    for (i=0; i< NP; i++) pfilt[i] = params[i];
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
       // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    pfilt[10] = 0.0;
    
    pfilt[4] = 0.0;
    pfilt[9] = 0.0;
    ResponseFreq(ll, Tend, pfilt, N, filtA[0], filtE[0]);
    
    pfilt[4] = 0.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(ll, Tend, pfilt, N, filtA[1], filtE[1]);
    
    pfilt[4] = PI/4.0;
    pfilt[9] = 0.0;
    ResponseFreq(ll, Tend, pfilt, N, filtA[2], filtE[2]);
    
    pfilt[4] = PI/4.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(ll, Tend, pfilt, N, filtA[3], filtE[3]);
    
    for (i=0; i< N; i++)
    {
        filtA[1][i] *= -1.0;
        filtE[1][i] *= -1.0;
        filtA[2][i] *= -1.0;
        filtE[2][i] *= -1.0;
    }
    
   //printf("\n");
        for (i=0; i< 4; i++)
        {
            NV[i] = 2.0*(fourier_nwip(filtA[i], AC, SN, N) + fourier_nwip(filtE[i], EC, SN, N))/Tobs;
            for (j=i; j< 4; j++)
            {
                MM[i][j] =  4.0*(fourier_nwip(filtA[i], filtA[j], SN, N) + fourier_nwip(filtE[i], filtE[j], SN, N))/Tobs;
                //printf("%e ", MM[i][j]);
            }
            //printf("\n");
        }
        
        for (i=0; i< 4; i++)
        {
            for (j=i; j< 4; j++)
            {
                MM[j][i] = MM[i][j];
            }
        }
        
        // a small stabilizer
        for (i=0; i< 4; i++) MM[i][i] += 0.1;
        
        //x = det(MM,4);
        
        //printf("%e\n", x);
        
       // printf("%e %e %e %e\n", MM[0][0], MM[1][1], MM[2][2], MM[1][2]);
       //printf("%e %e %e %e\n", NV[0], NV[1], NV[2], NV[3]);
        
        Inverse(MM, MI, 4);
        
        logL = 0.0;
        
        for (i=0; i< 4; i++)
        {
            aV[i] = 0.0;
            for (j=0; j< 4; j++)
            {
                aV[i] += MI[i][j]*NV[j];
                logL += 0.5*MI[i][j]*NV[i]*NV[j];
            }
        }
        
        printf("%e %f\n", logL, sqrt(2.0*logL));
        
        //printf("%e %e %e %e\n", aV[0], aV[1], aV[2], aV[3]);
        
        x = (aV[0]+aV[3]);
        x = x*x;
        y = (aV[1]-aV[2]);
        y = y*y;
        u = (aV[0]-aV[3]);
        u = u*u;
        v = (aV[1]+aV[2]);
        v = v*v;
    
    x = (aV[0]+aV[3]);
    x = x*x;
    y = (aV[1]-aV[2]);
    y = y*y;
    u = (aV[0]-aV[3]);
    u = u*u;
    v = (aV[1]+aV[2]);
    v = v*v;
        
        Ap = sqrt(x+y)+sqrt(u+v);
        Ac = sqrt(x+y)-sqrt(u+v);
        A = Ap + sqrt(Ap*Ap-Ac*Ac);
        
        //  + + - +   need -pi/2 on both
        //  + + + +   ?
        //   + + + -  ?
        //  + + - - ?
        //  + - + +
        //  + - + -   // phic = psi_true, psi = phic_true + pi/2
        //  + - - +
        //  + - - -
        
        //  - + + +  ?
        //  - + - + ?
        //  - + - -  ?
        //  - + + -
        //  - - + +
        //  - - + -  // phic = phic_true + pi/2, psi = psi_true
        //  - - - +
        //  - - - - ?
        
        x = atan2((aV[1]-aV[2]),(aV[0]+aV[3]));
        
        y = atan2(-(aV[1]+aV[2]),(aV[0]-aV[3]));
        
        psi = 0.25*(y-x);
        while(psi < 0.0) psi += PI;
        while(psi > PI) psi -= PI;
        
        
        phic = 0.25*(x+y);
        while(phic < 0.0) phic += PI;
        while(phic > PI) phic -= PI;
        
        cosi = Ac/A;
        
        scale = 2.0/A;
        
       printf("%f %f %f %f %f %f %f %f %f\n", scale, x, y, phic, psi, cosi, params[4], params[9], params[10]);
    
       // maxed on data without noise get
       // scale = 0.996465 phic = 2.837624 psi =  2.930106 cosi = 0.339382
        // maxed on data with noise get
        // scale = 0.992009 phic = 2.835180 psi = 2.927718 cosi = 0.337941
    
        // reference values
        // 1 2.848705 2.937147 0.339386
        
        pnew[0] = logL;
        pnew[1] = cosi;
        pnew[2] = scale;
        pnew[3] = psi;
        pnew[4] = phic;
        
        
 

    /*   Deallocate Arrays   */
    
     free(NV);
     free(aV);
     free_double_matrix(MM,4);
     free_double_matrix(MI,4);
     free_double_matrix(filtA,4);
     free_double_matrix(filtE,4);
    
    
    
    return;
}



void cholesky(double **A, double **C, int N)
{
    int *IPIV;
    int LWORK = N*N;
    int INFO;
    int i, j;
    double dx, dy;
    

    
    gsl_matrix *m = gsl_matrix_alloc (N, N);
    gsl_vector *s = gsl_vector_alloc (N);
    
    for (i = 0 ; i < N ; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            gsl_matrix_set(m, i, j, A[i][j]);
        }
    }
    
   /* gsl_linalg_cholesky_decomp1(m);
    for (i = 0; i < N; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            C[i][j] = gsl_matrix_get(m, i, j);
        }
    }
    */
    
    gsl_linalg_cholesky_decomp2(m, s);
    
    
    for (i = 0; i < N; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            C[i][j] = gsl_matrix_get(m, i, j)/gsl_vector_get(s, i);
        }
    }
    
    for (i = 0; i < N; i++)
    {
        for (j = i+1 ; j < N ; j++)
        {
            C[i][j] = 0.0;
        }
    }


    gsl_matrix_free(m);
    gsl_vector_free(s);
    
    return;
    
    
}



double det(double **A, int N)
{
    int *IPIV;
    int LWORK = N*N;
    int INFO;
    int i, j;
    double dx, dy;
    int s;
    
    gsl_permutation *p = gsl_permutation_alloc(N);
    
    gsl_matrix *m = gsl_matrix_alloc (N, N);
    
    for (i = 0 ; i < N ; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            gsl_matrix_set(m, i, j, A[i][j]);
        }
    }
    
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);
    
    dx = 1.0;
    for (i = 0; i < N; i++) dx *= gsl_matrix_get (m, i, i);
    dx = fabs(dx);
    
    //returns the absolute value of the determinant.
    
    gsl_permutation_free(p);
    gsl_matrix_free(m);
    
    
    return dx;
    
    
}



void getfreq(double *fnew, double *tf, double *Amp, double *Phase, double t, double fguess, double phic, double fRef_in, double m1_SI, double m2_SI, double chi1, double chi2, double distance, double tc)
{
    AmpPhaseFDWaveform *ap = NULL;
    double ep, u, v, tnew, x;
    double delT, delF, dtdf, fonfs;
    int ret;
    double M_sec;
    RealVector *freq;
    
    M_sec = (m1_SI+m2_SI) * MTSUN_SI/MSUN_SI;
    
    if(fguess < 1.0/Tobs) fguess = 1.0/Tobs;
    
    ep = 1.0e-6/M_sec;
    v = (4.0*PI*ep);
    u = (2.0*PI*ep*ep);
    
    if(fguess-ep < 0.0) ep = 0.5*fguess;
    
    freq = CreateRealVector((3));
    
    freq->data[0] = fguess-ep;
    freq->data[1] = fguess;
    freq->data[2] = fguess+ep;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap,freq,phic,fRef_in,m1_SI, m2_SI, chi1, chi2, distance);
    
    tnew = (ap->phase[2]-ap->phase[0])/v +tc;
    
    dtdf = (ap->phase[2]+ap->phase[0]-2.0*ap->phase[1])/u;
    
    delT = t-tnew;
    
    delF = delT/dtdf;
    
    *fnew = fguess + delF;
    
    *tf = tnew;
    
    x = h22fac*ap->amp[1];
    fonfs = fguess/fstar;
    x *= 8.0*fonfs*sin(fonfs);   // conversion to fractional frequency and leading order TDI transfer function

    *Amp = x;
    
    *Phase = ap->phase[1];
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
}


void timearray(double *params, RealVector *freq, long N, double *TF, AmpPhaseFDWaveform *ap)
{
    
    int flag;
    int i, j;
    double v;
    double tc,  deltaF, fmax;
    
    tc = params[5];    // merger time

    v = (4.0*PI*deltaF);
    
    for (i = 1; i < N-1; ++i) TF[i] = ( ap->phase[i+1]- ap->phase[i-1])/(2.0*PI*(freq->data[i+1]-freq->data[i-1])) + tc;
    
    TF[0] = TF[1];
    TF[N-1] = TF[N-2];
    
    j = N-1;
    flag = 0;
    for (i = 0; i < N-1; ++i)
    {
        // catch where time turns over
        if(TF[i+1] < TF[i] && flag == 0)
        {
            j = i;
            flag = 1;
        }
        // don't allow time to go too far into the past
        if(TF[i] < -Tpad) TF[i] = -Tpad;
    }
    
    // freeze time at turn over
    for (i = j; i < N; ++i)
    {
        TF[i] = TF[j];
    }
    
}


	
/*************************************************************************************/
/*                                                                                   */
/*                                    Subroutines                                    */
/*                                                                                   */
/*************************************************************************************/ 


/*************************************************************************/
/*        Rigid approximation position of each LISA spacecraft           */
/*************************************************************************/
void spacecraft(double t, double *x, double *y, double *z)
{

  double alpha;
  double beta1, beta2, beta3;
  double sa, sb, ca, cb;
 
  alpha = 2.*pi*fm*t + kappa0;

  beta1 = 0. + lambda0;
  beta2 = 2.*pi/3. + lambda0;
  beta3 = 4.*pi/3. + lambda0;

  sa = sin(alpha);
  ca = cos(alpha);


  sb = sin(beta1);
  cb = cos(beta1);
  x[1] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[1] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[1] = -sq3*AU*ec*(ca*cb + sa*sb);

 
  sb = sin(beta2);
  cb = cos(beta2);
  x[2] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[2] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[2] = -sq3*AU*ec*(ca*cb + sa*sb);

  sb = sin(beta3);
  cb = cos(beta3);
  x[3] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[3] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[3] = -sq3*AU*ec*(ca*cb + sa*sb);
  
}



int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
}

double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
    free(v);
}

double **double_matrix(int N, int M)
{
    int i;
    double **m = malloc( (N+1) * sizeof(double *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(double));
    }
    
    return m;
}

void free_double_matrix(double **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_double_vector(m[i]);
    free(m);
}


int **int_matrix(int N, int M)
{
    int i;
    int **m = malloc( (N+1) * sizeof(int *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(int));
    }
    
    return m;
}

void free_int_matrix(int **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_int_vector(m[i]);
    free(m);
}

double ***double_tensor(int N, int M, int L)
{
    int i,j;
    
    double ***t = malloc( (N+1) * sizeof(double **));
    for(i=0; i<N+1; i++)
    {
        t[i] = malloc( (M+1) * sizeof(double *));
        for(j=0; j<M+1; j++)
        {
            t[i][j] = malloc( (L+1) * sizeof(double));
        }
    }
    
    return t;
}

void free_double_tensor(double ***t, int N, int M)
{
    int i;
    
    for(i=0; i<N+1; i++) free_double_matrix(t[i],M);
    
    free(t);
}

void ang2pix_ring( const long nside, double theta, double phi, long *ipix) {
    /*
     c=======================================================================
     c     gives the pixel number ipix (RING)
     c     corresponding to angles theta and phi
     c=======================================================================
     */
    
    int nl2, nl4, ncap, npix, jp, jm, ipix1;
    double  z, za, tt, tp, tmp;
    int ir, ip, kshift;
    
    double piover2 = 0.5*M_PI;
    double twopi=2.0*M_PI;
    double z0=2.0/3.0;
    long ns_max=8192;
    
    if( nside<1 || nside>ns_max ) {
        fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
        exit(0);
    }
    
    if( theta<0. || theta>PI) {
        fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
        exit(0);
    }
    
    z = cos(theta);
    za = fabs(z);
    if( phi >= twopi)  phi = phi - twopi;
    if (phi < 0.)     phi = phi + twopi;
    tt = phi / piover2;//  ! in [0,4)
    
    nl2 = 2*nside;
    nl4 = 4*nside;
    ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
    npix  = 12*nside*nside;
    
    if( za <= z0 ) {
        
        jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
        jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
        
        ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
        kshift = 0;
        if (fmod(ir,2)==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise
        
        ip = (int)floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1;// ! in {1,4n}
        if( ip>nl4 ) ip = ip - nl4;
        
        ipix1 = ncap + nl4*(ir-1) + ip ;
    }
    else {
        
        tp = tt - floor(tt);//      !MOD(tt,1.d0)
        tmp = sqrt( 3.*(1. - za) );
        
        jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
        jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index
        
        ir = jp + jm + 1;//        ! ring number counted from the closest pole
        ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
        if( ip>4*ir ) ip = ip - 4*ir;
        
        ipix1 = 2*ir*(ir-1) + ip;
        if( z<=0. ) {
            ipix1 = npix - 2*ir*(ir+1) + ip;
        }
    }
    *ipix = ipix1 - 1;// ! in {0, npix-1}
    
}

void pix2ang_ring( long nside, long ipix, double *theta, double *phi) {
    /*
     c=======================================================================
     c     gives theta and phi corresponding to pixel ipix (RING)
     c     for a parameter nside
     c=======================================================================
     */
    
    int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
    double  fact1, fact2, fodd, hip, fihip;
    //      PARAMETER (pi     = 3.1415926535897932384626434d0)
    //      parameter (ns_max = 8192) ! 2^13 : largest nside available
    
    int ns_max=8192;
    
    if( nside<1 || nside>ns_max ) {
        fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
        exit(0);
    }
    npix = 12*nside*nside;      // ! total number of points
    if( ipix<0 || ipix>npix-1 ) {
        fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
        exit(0);
    }
    
    ipix1 = ipix + 1; // in {1, npix}
    nl2 = 2*nside;
    nl4 = 4*nside;
    ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
    fact1 = 1.5*nside;
    fact2 = 3.0*nside*nside;
    
    if( ipix1 <= ncap ) {  //! North Polar cap -------------
        
        hip   = ipix1/2.;
        fihip = floor(hip);
        iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
        iphi  = ipix1 - 2*iring*(iring - 1);
        
        *theta = acos( 1. - iring*iring / fact2 );
        *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
    }
    else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
        
        ip    = ipix1 - ncap - 1;
        iring = (int)floor( ip / nl4 ) + nside;// ! counted from North pole
        iphi  = (int)fmod(ip,nl4) + 1;
        
        fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
        *theta = acos( (nl2 - iring) / fact1 );
        *phi   = (1.*iphi - fodd) * PI /(2.*nside);
    }
    else {//! South Polar cap -----------------------------------
        
        ip    = npix - ipix1 + 1;
        hip   = ip/2.;
        /* bug corrige floor instead of 1.* */
        fihip = floor(hip);
        iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
        iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
        
        *theta = acos( -1. + iring*iring / fact2 );
        *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
    }
}


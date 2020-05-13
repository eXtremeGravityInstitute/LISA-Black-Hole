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


#include "IMRPhenomD.h"

#include <ctype.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include "Constants_search.h"
#include "Detector.h"
#include "Declarations_search.h"

#define NR_END 1
#define FREE_ARG char*

#ifndef _OPENMP
#define omp ignore
#endif

// OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o search search.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

// Linux
// gcc -std=gnu99 -fopenmp -w -o search search.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

//##############################################
//MT modifications

gsl_rng **rvec;
//##############################################

int main(int argc, char *argv[])
{
    double phi0;
    double deltaF;
    double m1_SI;
    double m2_SI;
    double chi1;
    double chi2;
    double distance;
    double tmerger;
    double zred;
    double Mtot, Mc;
    double f, cp, sp, p, tshift;
    double cv, sv, v, A, dtm, tlength, t0;
    double logL;
    int ret, seg;
    double *params, *pref, *zv;
    double *pmax;
    double HH, HD, DD;
    double thetaL, phiL, betaE, lambdaE, psi, iota;
    double fstart, fstop;
    double Tzero, beta;
    int ll, k;
    char command[1024];
   
    
    FILE *in;
    FILE *out;
    
    double x, dt, Tobs;
    int i, j, N;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    //##############################################
    //open MP modifications
    omp_set_num_threads(NC);
    rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (NC));
    for(i = 0 ; i< NC; i++){
        rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(rvec[i] , i);
    }
    //##############################################
    
    if(argc<2)
    {
        printf("./search segment\n");
        printf("segment numbers run from 0 to 15\n");
        return 0;
    }
    
    seg = atoi(argv[1]);
    dt = 10.0;
    N = 262144;
    Tobs = (double)(N)*dt;
    // the merger for LDC Radler is in segment 9.
    Tzero = (double)(seg)*Tobs;

    pmax = double_vector(NP);
    params = (double*)malloc(sizeof(double)*NP);
    
    double *AC, *EC, *TC, *SN;
    AC = (double*)malloc(sizeof(double)* (N));
    EC = (double*)malloc(sizeof(double)* (N));
    TC = (double*)malloc(sizeof(double)* (N));
    SN = (double*)malloc(sizeof(double)* (N));  // goes out further in frequency than we need.
    
    // Read in FFTed LDC data and A,E PSD from segmentAET.c
    sprintf(command, "AET_seg%d_f.dat", seg);
    in = fopen(command,"r");
    for(i=0; i< N; i++)
    {
        fscanf(in,"%lf%lf%lf%lf%lf\n", &f, &AC[i], &EC[i], &TC[i], &SN[i]);
    }
    fclose(in);
    
    search(pmax, AC, EC, SN, Tobs, Tzero, N);
    
    printf("\n Sky search \n");

    searchsky(pmax, AC, EC, SN, Tobs, Tzero, N);
    
    return 1;
    
    
}


void search(double *pmax, double *AC, double *EC, double *SN, double Tobs, double Tzero, int N)
{
    int i, j, k, mc, hold;
    double x, logLmax;
    double **paramx;
    double *logLx;
    double *min, *max;
    int *scount, *sacc;
    int **av, **cv;
    int ll=1;
    int *who;
    double *heat;
    int q;
    double m1, m2, DL;
    double alpha, beta;
    
    FILE *chain;
    FILE *levels;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    sacc = int_vector(NC);
    scount = int_vector(NC);
    who = int_vector(NC);
    max = double_vector(NP);
    min = double_vector(NP);
    heat = double_vector(NC);
    logLx = double_vector(NC);
    paramx = double_matrix(NC,NP);
    av = int_matrix(4,NC);
    cv = int_matrix(4,NC);
    
    for (i=0; i< NC; i++) who[i] = i;
    heat[0] = 1.0;
    for (i=1; i< NC; i++) heat[i] = heat[i-1]*1.5;
    
    for(j = 0; j < 4; j++)
    {
        for(k=0; k < NC; k++)
        {
            av[j][k] = 0;
            cv[j][k] = 1;
        }
    }
    
    for(k=0; k < NC; k++)
    {
        sacc[k] = 0;
        scount[k] = 1;
    }
    
    max[0] = log(1.0e8);
    max[1] = log(1.0e8);
    max[2] = 0.999;
    max[3] = 0.999;
    max[4] = PI;
    max[5] = Tzero+2.0*Tobs;  // merger could be up to 2 segments away
    max[6] = log(400.0);
    max[7] = 10.0;
    max[8] = 2.0*PI;
    max[9] = PI;
    max[10] = 1.0;
   
    min[0] = log(5.0e4);
    min[1] = log(5.0e4);
    min[2] = -0.999;
    min[3] = -0.999;
    min[4] = 0.0;
    min[5] = Tzero+1000.0;
    min[6] = log(0.1);
    min[7] = 0.1;
    min[8] = 0.0;
    min[9] = 0.0;
    min[10] = -1.0;
    

    for (i=0; i< NC; i++)
       {
           for (j=0; j< NP; j++) paramx[i][j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
           logLx[i] = log_likelihood_max_dual(ll, AC, EC, paramx[i], SN, N, Tobs, Tzero);
           
       }
    
    chain = fopen("search.dat","w");
    levels = fopen("likes.dat","w");
    
    logLmax = -1.0e10;
    
    for(mc = 1; mc < MS; mc++)
    {
    
     alpha = gsl_rng_uniform(r);
    
      if((NC > 1) && (alpha < 0.3))  // decide if we are doing a MCMC update of all the chains or a PT swap
       {
        // chain swap
       
        
        alpha = (double)(NC-1)*gsl_rng_uniform(r);
        j = (int)(alpha);
        beta = ((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
        alpha = log(gsl_rng_uniform(r));
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
           #pragma omp parallel for
           for(k=0; k < NC; k++)
            {
             update(mc, k, ll, logLx, paramx, min, max, who, heat, av, cv, AC, EC, SN, Tobs, Tzero, N, rvec[k]);
            }
       }
        
        for(k=0; k < NC; k++)
        {
          if(logLx[k] > logLmax)
          {
              logLmax = logLx[k];
              for(i = 0; i < NP; i++) pmax[i] = paramx[k][i];
          }
        }
        
        if(mc%100 == 0)
        {
            
            // clone the best solutioin into the hottest chain
            x = -1.0e4;
            for(k=0; k < NC; k++)
            {
              if(logLx[k] > x)
              {
                  x = logLx[k];
                  j = k;
              }
            }
            
            q = who[NC-1];
            logLx[q] = logLx[j];
            for(i = 0; i < NP; i++) paramx[q][i] = paramx[j][i];
            
        }
        
        if(mc%10 == 0)
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
            
              fprintf(levels, "%d ", mc);
               for(k = 0; k < NC; k++)
                {
                fprintf(levels, "%.12e ", logLx[who[k]]);
                }
                 fprintf(levels, "\n");
          }
        
              if(mc%10 == 0)
                {
                    q = who[0];
                    m1 = exp(paramx[q][0]);
                    m2 = exp(paramx[q][1]);
                    printf("%d %e %e %e %e %f %f %f %f\n", mc, logLx[q], paramx[q][5], m1, m2,
                    (double)(sacc[0])/(double)(scount[0]),
                    (double)(av[0][q])/(double)(cv[0][q]),
                    (double)(av[1][q])/(double)(cv[1][q]),
                    (double)(av[2][q])/(double)(cv[2][q]));
                }
        
     }
    
    fclose(chain);
    fclose(levels);
    
    chain = fopen("searchmax.dat","w");
     m1 = exp(pmax[0]);
     m2 = exp(pmax[1]);
     DL = exp(pmax[6]);
    fprintf(chain,"%.12e %.12e %.12e ", logLmax, m1, m2);
    for(i = 2; i < NP; i++)
     {
     if(i == 6)
      {
      fprintf(chain, "%.12e ", DL);
      }
      else
      {
       fprintf(chain, "%.15e ", pmax[i]);
      }
    }
    fclose(chain);

    free_int_vector(sacc);
    free_int_vector(scount);
    free_int_vector(who);
    free_double_vector(max);
    free_double_vector(min);
    free_double_vector(heat);
    free_double_vector(logLx);
    free_double_matrix(paramx,NC);
    free_int_matrix(av,4);
    free_int_matrix(cv,4);
   
}

void update(int mc, int k, int ll, double *logLx, double **paramx, double *min, double *max, int *who, double *heat, int **av, int **cv, double *AC, double *EC, double *SN, double Tobs, double Tzero, int N, gsl_rng *r)
{
    int q, i, j;
    int typ, flag;
    double lm1, lm2;
    double alpha, beta;
    double logLy, H, x;
    double *paramy;
    double **Fisher;
    double **evec;
    double *eval;
    double a, b;
    double **Chl, **Cov;
    double *zv;
    
    paramy = double_vector(NP);
    
    q = who[k];
    
    for(i = 0; i < NP; i++) paramy[i] = paramx[q][i];
    
    alpha = gsl_rng_uniform(r);
    
    if(mc < MS/4)
    {
        a = 0.5;
        b = 0.3;
    }
    else if (mc < MS/2)
    {
        a = 0.8;
        b = 0.4;
    }
    else
    {
        a = 0.9;
        b = 0.2;
    }
    
    if(alpha > a ) // uniform draw
    {
        typ = 0;
        
        
        for (j=0; j< NP; j++) paramy[j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
        
        // keep tc at current value
        beta = gsl_rng_uniform(r);
        if(beta > 0.1) paramy[5] = paramx[q][5];
        
    }
    else if (alpha > b && logLx[q] > 20.0)// fisher jump
    {
        typ = 1;
        
        
        Fisher = double_matrix(NV,NV);
        evec = double_matrix(NV,NV);
        eval = double_vector(NV);
        
        FisherMax(ll, paramx[q], Fisher, Tobs, Tzero);
        FisherEvec(Fisher, eval, evec, NV);
        
        beta = gsl_rng_uniform(r);
        i = (int)(beta*(NV));
        beta = sqrt(heat[k])*eval[i]*gsl_ran_gaussian(r,1.0);
        for(j = 0; j < NV; j++) paramy[j] += beta*evec[i][j];
        
        free_double_matrix(Fisher,NV);
        free_double_matrix(evec,NV);
        free_double_vector(eval);
        
    }
    else
    {
        typ = 2;
        
        beta = gsl_rng_uniform(r);
        
        if(beta > 0.5)
        {
        paramy[0] += gsl_ran_gaussian(r,0.1);
        paramy[1] += gsl_ran_gaussian(r,0.1);
        paramy[2] += gsl_ran_gaussian(r,0.1);
        paramy[3] += gsl_ran_gaussian(r,0.1);
        }
        else
        {
            paramy[0] += gsl_ran_gaussian(r,0.01);
            paramy[1] += gsl_ran_gaussian(r,0.01);
            paramy[2] += gsl_ran_gaussian(r,0.05);
            paramy[3] += gsl_ran_gaussian(r,0.05);
        }
        
    }
    
    lm1 = paramy[0];
    lm2 = paramy[1];
    if(lm2 > lm1)
    {
        paramy[0] = lm2;
        paramy[1] = lm1;
    }
    
    cv[typ][k]++;
    
    
    for(i = 0; i < NV; i++)
    {
        if(paramy[i] > max[i] || paramy[i] < min[i])
        {
            paramy[i] = min[i]+(max[i]-min[i])*gsl_rng_uniform(r);
        }
    }

    
    
    logLy = log_likelihood_max_dual(ll, AC, EC, paramy, SN, N, Tobs, Tzero);
    
    H = (logLy-logLx[q])/heat[k];
     
    alpha = log(gsl_rng_uniform(r));
        
     if(H > alpha)
     {
         // copy over new state if accepted
         logLx[q] = logLy;
         for(i = 0; i < NP; i++) paramx[q][i] = paramy[i];
         av[typ][k]++;
     }
    
    free_double_vector(paramy);
    
}

// maps tc so that merger time at detector is held fixed
double Tmap(double *params, double tdet)
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
    
    spacecraft(tdet, x, y, z);
        
        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
    kdotx = (xa*kv[1]+ya*kv[2]+za*kv[3])/clight;
    
    free_double_vector(kv);
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    
    return kdotx;
}

void searchsky(double *params, double *AC, double *EC, double *SN, double Tobs, double Tzero, int N)
{
    int i, j, k, mc, hold;
    double x, logLmax, ts;
    double **paramx;
    double *pmax;
    double *logLx;
    double *min, *max;
    int *scount, *sacc;
    int **av, **cv;
    int ll=1;
    int *who;
    double *heat;
    int q;
    double m1, m2, DL;
    double alpha, beta;
    
    FILE *chain;
    FILE *levels;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    sacc = int_vector(NC);
    scount = int_vector(NC);
    who = int_vector(NC);
    max = double_vector(NP);
    min = double_vector(NP);
    pmax = double_vector(NP);
    heat = double_vector(NC);
    logLx = double_vector(NC);
    paramx = double_matrix(NC,NP);
    av = int_matrix(4,NC);
    cv = int_matrix(4,NC);

    
    for (i=0; i< NC; i++) who[i] = i;
    heat[0] = 1.0;
    for (i=1; i< NC; i++) heat[i] = heat[i-1]*1.5;
    
    for(j = 0; j < 4; j++)
    {
        for(k=0; k < NC; k++)
        {
            av[j][k] = 0;
            cv[j][k] = 1;
        }
    }
    
    for(k=0; k < NC; k++)
    {
        sacc[k] = 0;
        scount[k] = 1;
    }
    
    max[0] = log(1.0e9);
    max[1] = log(1.0e9);
    max[2] = 0.999;
    max[3] = 0.999;
    max[4] = PI;
    max[5] = params[5]+600.0;  // within light travel time
    max[6] = log(1000.0);
    max[7] = 1.0;
    max[8] = 2.0*PI;
    max[9] = PI;
    max[10] = 1.0;
   
    min[0] = log(1.0e3);
    min[1] = log(1.0e3);
    min[2] = -0.999;
    min[3] = -0.999;
    min[4] = 0.0;
    min[5] = params[5]-600.0;
    min[6] = log(0.1);
    min[7] = -1.0;
    min[8] = 0.0;
    min[9] = 0.0;
    min[10] = -1.0;
    
    for (i=0; i< NC; i++)
    {
        for (j=0; j< NP; j++) paramx[i][j] = params[j];
    }
    
    for (i=0; i< NC; i++)
       {
           do
           {
           for (j=7; j< 9; j++) paramx[i][j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
           ts = Tmap(paramx[i], params[5]);
           paramx[i][5] = params[5]-ts;
           //logLx[i] = likelihoodFstat(ll, paramx[i], N, AC, EC, SN, Tobs, Tzero);
           logLx[i] = likelihoodFstatTmax(ll, paramx[i], N, AC, EC, SN, Tobs, Tzero);
           }while(paramx[i][5] > max[5] || paramx[i][5] < min[5]);
       }
    
    chain = fopen("searchsky.dat","w");
    levels = fopen("likesky.dat","w");
    
    logLmax = -1.0e10;
    
    for(mc = 1; mc < MSS; mc++)
    {
    
     alpha = gsl_rng_uniform(r);
    
      if((NC > 1) && (alpha < 0.3))  // decide if we are doing a MCMC update of all the chains or a PT swap
       {
        // chain swap
       
        
        alpha = (double)(NC-1)*gsl_rng_uniform(r);
        j = (int)(alpha);
        beta = ((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
        alpha = log(gsl_rng_uniform(r));
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
           #pragma omp parallel for
           for(k=0; k < NC; k++)
            {
             updatesky(mc, k, ll, logLx, paramx, min, max, who, heat, av, cv, AC, EC, SN, Tobs, Tzero, N, rvec[k]);
            }
       }
        
        for(k=0; k < NC; k++)
        {
          if(logLx[k] > logLmax)
          {
              logLmax = logLx[k];
              for(i = 0; i < NP; i++) pmax[i] = paramx[k][i];
          }
        }
        
        if(mc%100 == 0)
        {
            
            // clone the best solution into the hottest chain
            x = -1.0e4;
            for(k=0; k < NC; k++)
            {
              if(logLx[k] > x)
              {
                  x = logLx[k];
                  j = k;
              }
            }
            
            q = who[NC-1];
            logLx[q] = logLx[j];
            for(i = 0; i < NP; i++) paramx[q][i] = paramx[j][i];
            
        }
        
        if(mc%10 == 0)
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
            
              fprintf(levels, "%d ", mc);
               for(k = 0; k < NC; k++)
                {
                fprintf(levels, "%.12e ", logLx[who[k]]);
                }
                 fprintf(levels, "\n");
          }
        
              if(mc%10 == 0)
                {
                    q = who[0];
                    printf("%d %e %e %e %e %f %f %f %f\n", mc, logLx[q], paramx[q][5], paramx[q][7], paramx[q][8],
                    (double)(sacc[0])/(double)(scount[0]),
                    (double)(av[0][q])/(double)(cv[0][q]),
                    (double)(av[1][q])/(double)(cv[1][q]),
                    (double)(av[2][q])/(double)(cv[2][q]));
                }
        
     }
    
    fclose(chain);
    fclose(levels);
    
    chain = fopen("skymax.dat","w");
     m1 = exp(pmax[0]);
     m2 = exp(pmax[1]);
     DL = exp(pmax[6]);
    fprintf(chain,"%.12e %.12e %.12e ", logLmax, m1, m2);
    for(i = 2; i < NP; i++)
     {
     if(i == 6)
      {
      fprintf(chain, "%.12e ", DL);
      }
      else
      {
       fprintf(chain, "%.15e ", pmax[i]);
      }
    }
    fclose(chain);

    free_int_vector(sacc);
    free_int_vector(scount);
    free_int_vector(who);
    free_double_vector(max);
    free_double_vector(min);
    free_double_vector(pmax);
    free_double_vector(heat);
    free_double_vector(logLx);
    free_double_matrix(paramx,NC);
    free_int_matrix(av,4);
    free_int_matrix(cv,4);
   
}

void updatesky(int mc, int k, int ll, double *logLx, double **paramx, double *min, double *max, int *who, double *heat, int **av, int **cv, double *AC, double *EC, double *SN, double Tobs, double Tzero, int N, gsl_rng *r)
{
    int q, i, j;
    int typ, flag;
    double lm1, lm2;
    double alpha, beta;
    double tx, ty;
    double logLy, H, x;
    double *paramy;
    double **Fisher;
    double **evec;
    double *eval;
    double a, b;
    double pDx, pDy;
    double **Chl, **Cov;
    double *zv;
    
    paramy = double_vector(NP);
    
    q = who[k];
    
    for(i = 0; i < NP; i++) paramy[i] = paramx[q][i];
    
    alpha = gsl_rng_uniform(r);
    
    if(mc < MSS/4)
    {
        a = 0.6;
        b = 0.3;
    }
    else if (mc < MSS/2)
    {
        a = 0.8;
        b = 0.4;
    }
    else
    {
        a = 0.9;
        b = 0.2;
    }
    
    if(alpha > a ) // uniform draw on sky
    {
        typ = 0;
        
        
        for (j=7; j< 9; j++) paramy[j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
    
        // time delay for current sky location
        tx = Tmap(paramx[q], paramx[q][5]);
        // time delay for new current sky location
        ty = Tmap(paramy, paramy[5]);
        
        paramy[5] = paramx[q][5]+tx-ty;
        
        
    }
    else if (alpha > b) // small jump in time and sky
    {
        typ = 1;
        
        beta = gsl_rng_uniform(r);
        if(beta > 0.5)
        {
        paramy[7] += gsl_ran_gaussian(r,0.01);
        paramy[8] += gsl_ran_gaussian(r,0.02);
        paramy[5] += gsl_ran_gaussian(r,1.0);
        }
        else
        {
         paramy[7] += gsl_ran_gaussian(r,0.05);
         paramy[8] += gsl_ran_gaussian(r,0.1);
         paramy[5] += gsl_ran_gaussian(r,5.0);
        }
 
        
    }
    else  // small jump in  sky
    {
        typ = 2;
        
        if(beta > 0.5)
        {
        paramy[7] += gsl_ran_gaussian(r,0.01);
        paramy[8] += gsl_ran_gaussian(r,0.02);
        }
        else
        {
         paramy[7] += gsl_ran_gaussian(r,0.05);
         paramy[8] += gsl_ran_gaussian(r,0.1);
        }
        
        // time delay for current sky location
        tx = Tmap(paramx[q], paramx[q][5]);
        // time delay for new current sky location
        ty = Tmap(paramy, paramy[5]);
        
        paramy[5] = paramx[q][5]+tx-ty;
        
    }
    
    if(paramy[7] > 1.0)  paramy[7] = 0.99;
    if(paramy[7] < -1.0)  paramy[7] = -0.99;
    if(paramy[8] > TPI)  paramy[8] -= TPI;
    if(paramy[8] < 0.0)  paramy[8] += TPI;
    

    cv[typ][k]++;
    
    if(paramy[5] < min[5] || paramy[5] > max[5])
    {
        logLy = -1.0e6;
    }
    else
    {
     logLy = likelihoodFstat(ll, paramy, N, AC, EC, SN, Tobs, Tzero);
    }
    
    pDx = paramx[q][6];   // uniform in distance prior
    pDy = paramy[6];   // uniform in distance prior
    
    H = (logLy-logLx[q])/heat[k] + pDy - pDx;
     
    alpha = log(gsl_rng_uniform(r));
        
     if(H > alpha)
     {
         // copy over new state if accepted
         logLx[q] = logLy;
         for(i = 0; i < NP; i++) paramx[q][i] = paramy[i];
         av[typ][k]++;
     }
    
    free_double_vector(paramy);
    
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
        
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = gsl_matrix_get(evec, j, i);
            }
            
        }
        
        for (i = 0; i < d; i++)
        {
            // turn into 1-sigma jump amplitudes
            ej[i] = 1.0/sqrt(fabs(ej[i]));
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

double fourier_nwip_Tshift(double delT, double Tobs, double *a, double *b, double *Sn, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    double f, c, s;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-i;
        
        f = (double)(i)/Tobs;
        c = cos(TPI*f*delT);
        s = sin(TPI*f*delT);
        
        ReA = a[j]*c+a[k]*s;
        ImA = -a[j]*s+a[k]*c;
        
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product/(Sn[i]);
    }
    
    return(4.0*arg);
    
}

void StartStop(int ll, double *params, double Tobs, double tstart, double tstop, double *fstart, double *fstop, double *frg)
{
    
    double m1, m2, m1_SI, m2_SI, chi1, chi2, tc;
    double Mtot, eta, Mc, af, fr;
    double fmin, fmax;
    double fnew, tf;
    int i;
    
    if(ll == 0)
    {
    m1 = params[0];
    m2 = params[1];
    }
    else
    {
    m1 = exp(params[0]);
    m2 = exp(params[1]);
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
        getfreq(&fnew, &tf, tstart, fmin, m1_SI, m2_SI, chi1, chi2, tc);
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
            getfreq(&fnew, &tf, tstop, fmax, m1_SI, m2_SI, chi1, chi2, tc);
            if(fnew < fmin) fnew = fmin+1.0/Tobs;
            fmax = fnew;
            i++;
        }while(i < 10 && fabs(tf-tstop) > 1.0);
    }
    
    if(fmax < fmin) fmax = 2.0*fmin;
    
    *fstart = fmin;
    *fstop = fmax;
    
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
    
    // Calculate the power spectral density of the detector noise at the given frequency
    
    *SAE = LC*16.0/3.0*pow(sin(fonfs),2.0)*( (2.0+cos(fonfs))*(Sps) + 2.0*(3.0+2.0*cos(fonfs)+cos(2.0*fonfs))*(Sacc/pow(2.0*PI*f,4.0)*(1.0+red)) ) / pow(2.0*L,2.0);
    
    *SXYZ = LC*4.0*pow(sin(fonfs),2.0)*( 4.0*(Sps) + 8.0*(1.0+pow(cos(fonfs),2.0))*(Sacc/pow(2.0*PI*f,4.0)*(1.0+red)) ) / pow(2.0*L,2.0);
    
}


void PDwave(int ll, double *wavef, double *params, int N, double Tobs, double Tzero)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    double phi0, mc, q,  m1_SI, m2_SI, chi1, chi2, f_min, f_max, distance;
    int ret, flag, i, k;
    double p, cp, sp, Amp;
    double f, x, y, deltaF, ts, fr;
    double tc, fstart, fstop;
    int imin, imax;
    
    double mt, eta, dm;
    
    RealVector *freq;
    
    chi1 = params[2];
    chi2 = params[3];
    phi0 = params[4];
    tc = params[5];
    ts = Tzero+Tobs-tc;
    if(ll == 0)
    {
        distance = params[6]*1.0e9*PC_SI;
        m1_SI = params[0]*MSUN_SI;
        m2_SI = params[1]*MSUN_SI;
    }
    else
    {
        distance = exp(params[6])*1.0e9*PC_SI;
        m1_SI = exp(params[0])*MSUN_SI;
        m2_SI = exp(params[1])*MSUN_SI;
    }
    
    StartStop(ll, params, Tobs, Tzero, Tzero+Tobs, &fstart, &fstop, &fr);
    
    imin = (int)(fstart*Tobs);
    imax = (int)(fstop*Tobs)+1;
    
    if(imin < 0) imin = 0;
    if(imax < imin) imax = imin+1;
    if(imax > N/2) imax = N/2;
    
    freq = CreateRealVector((imax-imin));
    for (i=0; i< (imax-imin); i++) freq->data[i] = (double)(i+imin)/Tobs;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, 0.0, PDfref, m1_SI, m2_SI, chi1, chi2, distance);
    
    for (i=0; i< N; i++) wavef[i] = 0.0;
    
    for (i=imin; i< imax; i++)
    {
        k = i-imin;
        f = freq->data[k];
        x = f/fstar;
        Amp = 8.0*x*sin(x)*h22fac*ap->amp[k];
        p = TPI*f*ts-ap->phase[k]-2.0*phi0;
        wavef[i] = Amp*cos(p);
        wavef[N-i] = Amp*sin(p);
    }
   
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
}

void FisherMax(int ll, double *params, double **Fisher, double Tobs, double Tzero)
{
    double *AAmp, *EAmp, *APhase, *EPhase, *TF, *FF, *PF, *AF;
    double *PFref, *AFref, *TFref;
    double af, fr, df, DT, fac, deltaF, f, fmax, fmin, x, y, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW, lx;
    double A, P;
    double px, fnew, kxm;
    double *Xint, *SAE;
    double *phaseP, *phaseM, *ampP, *ampM;
    double **Dphase, **Damp;
    double *epsilon;

    
    int NFmax = 10000;
    
    int kmin, kmax, flag;
    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FF = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
    
    // Note: This call sets the FF array, which then gets held and used even when the parameters change a little
    SetUp(ll, params, NFmax, &NF, FF, TF, PF, AF, Tobs, Tzero);
    
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
    // [7] E/A amplitude [8] PhaseA-PhaseE
    
    // combines A annd E channel amplitudes
    fac = (1.0+params[7]*params[7]);
    
    // will take log derivatives for m1, m2
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-4;
    epsilon[5] = 1.0e-1;
    
    phaseP = (double*)malloc(sizeof(double)* (NF));
    phaseM = (double*)malloc(sizeof(double)* (NF));
    
    ampP = (double*)malloc(sizeof(double)* (NF));
    ampM = (double*)malloc(sizeof(double)* (NF));
    
    Dphase = double_matrix(NV, NF);
    Damp = double_matrix(NV, NF);
    
    Xint = (double*)malloc(sizeof(double)* (NF));
    SAE = (double*)malloc(sizeof(double)* (NF));
    
    
    // masses and spins
    for (i=0; i< NV; i++)
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


        
            if(i == 4)  // Put in coallesence phase manually. Note this is orbital phase, so times two in GW phase
              {
                  
               x = -4.0*epsilon[4];
               for (k=0; k< NF; k++)
               {
                   Dphase[i][k] = x;
                   Damp[i][k] = 0.0;
               }
               
             }
            else if(i == 5)  // tc is not included in the phase subroutine
             {
               for (k=0; k< NF; k++)
               {
                   x = -4.0*PI*FF[k]*epsilon[5];
                   Dphase[i][k] = x;
                   Damp[i][k] = 0.0;
               }
             }
            else
            {
                
            Intrinsic(ll,paramsP, NF, FF, TF, phaseP, ampP, Tobs, Tzero);
            Intrinsic(ll,paramsM, NF, FF, TF, phaseM, ampM, Tobs, Tzero);
            
            for (k=0; k< NF; k++)
            {
                Dphase[i][k] = (phaseP[k]-phaseM[k]);
                Damp[i][k] = (ampP[k]-ampM[k]);
            }
                
            }
        
    }

    

    for (k=0; k< NF; k++) instrument_noise(FF[k], &SAE[k], &x);

    gsl_interp_accel *Iacc = gsl_interp_accel_alloc();
    gsl_spline *Ispline = gsl_spline_alloc (gsl_interp_cspline, NF);

   
    for (i=0; i< NV; i++)
    {
        for (j=i; j< NV; j++)
        {
  
            for (k=0; k< NF; k++)
            {
            Xint[k] = 4.0*fac*(Damp[i][k]*Damp[j][k]+AFref[k]*AFref[k]*Dphase[i][k]*Dphase[j][k])/SAE[k];
            }
            
            gsl_spline_init(Ispline, FF, Xint, NF);
    
            Fisher[i][j] = gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
            
            Fisher[i][j] /= (4.0*epsilon[i]*epsilon[j]);
        }
        
    }
    
    for (i=0; i< NV; i++)
    {
        for (j=i+1; j< NV; j++)
        {
            Fisher[j][i] =  Fisher[i][j];
        }
    }
    
    free_double_matrix(Dphase,NV);
    free_double_matrix(Damp,NV);
 
    free(epsilon);
    free(SAE);
    free(Xint);
    free(FF);
    free(TF);
    free(AF);
    free(PF);
    free(TFref);
    free(AFref);
    free(PFref);
    
    free(paramsP);
    free(paramsM);
    
    free(phaseP);
    free(phaseM);
    free(ampP);
    free(ampM);

    gsl_spline_free(Ispline);
    gsl_interp_accel_free(Iacc);
   
    
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




void SetUp(int ll, double *params, int NFmax, int *NFS, double *FF, double *TF, double *PF, double *AF, double Tobs, double Tzero)
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
    
    Mtot = (m1+m2)*TSUN;
    eta = m1*m2/((m1+m2)*(m1+m2));
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0)*TSUN;
    
    StartStop(ll, params, Tobs, Tzero, Tzero+Tobs, &fmin, &fmax, &fr);
    
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
        if(df < dfmin) df = dfmin;
        if(df > dfmax) df = dfmax;
        f += df;
        NF++;
    }while(f < fmax);
    
    //printf("%d %e\n", NF, log10(Mtot/TSUN));
    
     // Need to catch is NF > NFmax
    
    if(NF > NFmax)
    {
        NF = NFmax;
        df = (fmax-fmin)/(double)(NFmax-1);
        for (i=0; i< NF; i++)
        {
            FF[i] = fmin + df*(double)(i);
        }
    }
    else
    {
    
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
    }
    
    if(NF < 4)
    {
        NF = 4;
        df = (fmax-fmin)/3.0;
        FF[0] = fmin;
        for (i=1; i< NF; i++) FF[i] = fmin+df*(double)(i);
    }
    
    Intrinsic(ll, params, NF, FF, TF, PF, AF, Tobs, Tzero);

    *NFS = NF;
    
}

void Intrinsic(int ll, double *params, int NF, double *FF, double *TF, double *PF, double *AF, double Tobs, double Tzero)
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
    told = Tzero+ap->time[0]+tc;
    for (i=0; i< NF; i++)
    {
        PF[i] = ap->phase[i];
        
        AF[i] =  h22fac*ap->amp[i];
        fonfs = freq->data[i]/fstar;
        AF[i] *= (8.0*fonfs*sin(fonfs));   // conversion to fractional frequency and leading order TDI transfer function
        
        t = Tzero+ap->time[i]+tc;
        if(t < told && flag == 0)
        {
            flag = 1;
            ii = i-1;
            tx = told;
        }
        TF[i] = t;
        if(flag == 1) TF[i] = tx +(double)(i-ii)*Mtot;
        told = t;
    }
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);

    
}



double log_likelihood_max_dual(int ll, double *A, double *E, double *params, double *SN, int N, double Tobs, double Tzero)
{
    int i, j, k, NMAX;
    int ii, jj, Nend;
    double sum;
    double HH, LL;
    double logL;
    double logLfs;
    double fmax;
    double HA, HE, LD, dt, x;
    double normA, normE, deltH, pshiftA, pshiftE;
    double *AS, *ES;
    double *AC, *AF;
    double *EC, *EF;
    double *h;
    
    //Chi^2 = (d-h|d-h)/Sn
    
    dt = Tobs/(double)(N);
    
    h = (double*)malloc(sizeof(double)*(N));
    
    PDwave(ll, h, params, N, Tobs, Tzero);
    
    AS = double_vector(N); ES = double_vector(N);
    AC=double_vector(N);  AF=double_vector(N);
    EC=double_vector(N);  EF=double_vector(N);
    
    pbt_shift(AC, AF, A, h, SN, N);
    pbt_shift(EC, EF, E, h, SN, N);
    
    // only allow time shifts up to +/- Tobs/8 (otherwise frequency range can be way off)
    Nend = N/8;
    
    for(i = 0; i < Nend/2; i++) AS[i+N/2] = sqrt(AC[i]*AC[i]+AF[i]*AF[i]);
    for(i = -Nend/2; i < 0; i++) AS[i+N/2] = sqrt(AC[N+i]*AC[N+i]+AF[N+i]*AF[N+i]);
    
    for(i = 0; i < Nend/2; i++) ES[i+N/2] = sqrt(EC[i]*EC[i]+EF[i]*EF[i]);
    for(i = -Nend/2; i < 0; i++) ES[i+N/2] = sqrt(EC[N+i]*EC[N+i]+EF[N+i]*EF[N+i]);
    
    x = 0;
    for (i = -Nend/2; i < Nend/2; ++i)
    {
        if((AS[i+N/2]+ES[i+N/2]) > x)
        {
            x = AS[i+N/2]+ES[i+N/2];
            k = i;
        }
    }
    
    HA = 2.0*(double)(N)*(AS[k+N/2]);
    HE = 2.0*(double)(N)*(ES[k+N/2]);
    
    deltH = dt*(double)(k);
    
    if(k < 0) k += N;
    
    HH = fourier_nwip(h, h, SN, N);
    
    // Inverse FFTs in fourier_nwip are un-normlaized
    HH /= Tobs;
    HA /= Tobs;
    HE /= Tobs;

    pshiftA = atan2(AF[k],AC[k]);
    pshiftE = atan2(EF[k],EC[k]);
    normA = HA/HH;
    normE = HE/HH;
    
    // reference to channel A
    if(ll==0)
    {
     params[6] /= normA;
    }
    else
    {
     params[6] -= log(normA);
    }
    params[5] += deltH;
    params[4] -= pshiftA/2.0;  //
    
    // [7] E/A amplitude [8] PhaseA-PhaseE
    
    params[7] = normE/normA;
    params[8] = 0.5*(pshiftA-pshiftE);
    
    logL = (HA*HA+HE*HE)/(2.0*HH);
    
    free_double_vector(AC);  free_double_vector(AF);
    free_double_vector(EC);  free_double_vector(EF);

    free_double_vector(AS);
    free_double_vector(ES);
    
    free(h);

    
    return logL;
}

double log_likelihood_max(int ll, double *H, double *params, double *SN, int N, double Tobs, double Tzero)
{
    int i, j, k, NMAX;
    int ii, jj;
    double sum;
    double HH, LL;
    double logL;
    double logLfs;
    double fmax;
    double HD, LD, dt, x;
    double normH, deltH, pshiftH;
    double *HS;
    double *HC, *HF;
    double *h;
    
    //Chi^2 = (d-h|d-h)/Sn
    
    dt = Tobs/(double)(N);
    
    h = (double*)malloc(sizeof(double)*(N));
    
    PDwave(ll, h, params, N, Tobs, Tzero);
    
    HS = double_vector(N);
    HC=double_vector(N);  HF=double_vector(N);
    
    pbt_shift(HC, HF, H, h, SN, N);
    
    for(i = 0; i < N/2; i++) HS[i+N/2] = sqrt(HC[i]*HC[i]+HF[i]*HF[i]);
    for(i = -N/2; i < 0; i++) HS[i+N/2] = sqrt(HC[N+i]*HC[N+i]+HF[N+i]*HF[N+i]);
    
    x = 0;
    for (i = -N/2; i < N/2; ++i)
    {
        if(HS[i+N/2] > x)
        {
            x = HS[i+N/2];
            k = i;
        }
    }
    
    HD = 2.0*(double)(N)*HS[k+N/2];
    
    deltH = dt*(double)(k);
    
    if(k < 0) k += N;
    
    pshiftH = atan2(HF[k],HC[k]);
    
    HH = fourier_nwip(h, h, SN, N);
    
    // Inverse FFTs in fourier_nwip are un-normlaized
    HH /= Tobs;
    HD /= Tobs;

    normH = HD/HH;

    logL = (HD*HD/HH)/2.0;
    
    // Shift the waverform to match
    
    
    if(ll==0)
    {
     params[6] /= normH;
    }
    else
    {
     params[6] -= log(normH);
    }
    params[5] += deltH;
    params[4] -= pshiftH/2.0;  // PhenomD uses orbital phase, while here we have GW phase
    
    free_double_vector(HC);  free_double_vector(HF);

    free_double_vector(HS);
    
    free(h);

    
    return logL;
}

double likelihoodFstat(int ll, double *params, int N, double *AC, double *EC, double *SN, double Tobs, double Tzero)
{
    
    /*   Indicies   */
    int i,j, k, n, m, imax;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx;
    
    double *pfilt;
    
    double **filtA, **filtE;
    
    double **MM, **MI;
    
    double *KV, *aV;
    
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
    
    pfilt = (double*)malloc(sizeof(double)* (NP));
    
    filtA = double_matrix(4,N);
    filtE = double_matrix(4,N);
    
    KV = (double*)malloc(sizeof(double)* (4));
    aV = (double*)malloc(sizeof(double)* (4));
    MM = double_matrix(4,4);
    MI = double_matrix(4,4);
    
    for (i=0; i< NP; i++) pfilt[i] = params[i];
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
       // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    pfilt[10] = 0.0;
    
    pfilt[4] = 0.0;
    pfilt[9] = 0.0;
    ResponseFreq(ll, pfilt, N, filtA[0], filtE[0], Tobs, Tzero);
    
    pfilt[4] = PI/2.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(ll, pfilt, N, filtA[1], filtE[1], Tobs, Tzero);
    
    pfilt[4] = 3.0*PI/4.0;
    pfilt[9] = 0.0;
    ResponseFreq(ll, pfilt, N, filtA[2], filtE[2], Tobs, Tzero);
    
    pfilt[4] = PI/4.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(ll, pfilt, N, filtA[3], filtE[3], Tobs, Tzero);
    
    
        for (i=0; i< 4; i++)
        {
            KV[i] = 2.0*(fourier_nwip(filtA[i], AC, SN, N) + fourier_nwip(filtE[i], EC, SN, N))/Tobs;
            for (j=i; j< 4; j++)
            {
                MM[i][j] =  4.0*(fourier_nwip(filtA[i], filtA[j], SN, N) + fourier_nwip(filtE[i], filtE[j], SN, N))/Tobs;
               
            }
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
        
        Inverse(MM, MI, 4);
        
        logL = 0.0;
        
        for (i=0; i< 4; i++)
        {
            aV[i] = 0.0;
            for (j=0; j< 4; j++)
            {
                aV[i] += MI[i][j]*KV[j];
                logL += 0.5*MI[i][j]*KV[i]*KV[j];
            }
        }
    
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
        

      if(ll == 0)
      {
        params[6] *= scale;
      }
      else
      {
       params[6] += log(scale);
      }
    
       params[10] = cosi;
       params[9] = psi;
       params[4] = phic;

    /*   Deallocate Arrays   */
    
     free(KV);
     free(aV);
     free_double_matrix(MM,4);
     free_double_matrix(MI,4);
     free_double_matrix(filtA,4);
     free_double_matrix(filtE,4);
    
    
    return logL;
}

// Fstat likelihood that allows for a small (roughly +/- 30 s) time maximization
double likelihoodFstatTmax(int ll, double *params, int N, double *AC, double *EC, double *SN, double Tobs, double Tzero)
{
    
    /*   Indicies   */
    int i,j, k, n, m, imax;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx;
    
    double *pfilt;
    
    double **filtA, **filtE;
    
    double **MM, **MI;
    
    double *KV, *aV;
    
    double *KP, *KM;
    
    double AR, AI, ER, EI;
    
    double t, f, kdotx, Amp, Phase, Fp, Fc;
    
    double A, P, px, tc, tcs, pxi, pxj;
    
    double x, y, u, v, xold, yold;
    
    double cp, sp;
    
    int NA, flag;
    
    double cx, sx, logL, logLmax;
    
    double delT = 15.0;
    
    int nn;
    
    double iota, cosi;
    double Ap, Ac, scale;
    
    double a1, a2, a3, a4;
    double cpsi, spsi;
    double psi, phic;
    double cpc, spc;
    
    double T0, T1, T2, a, b, c, dtx;
    
    FILE *out;
    
    pfilt = (double*)malloc(sizeof(double)* (NP));
    
    filtA = double_matrix(4,N);
    filtE = double_matrix(4,N);
    
    KV = (double*)malloc(sizeof(double)* (4));
    KP = (double*)malloc(sizeof(double)* (4));
    KM = (double*)malloc(sizeof(double)* (4));
    
    aV = (double*)malloc(sizeof(double)* (4));
    MM = double_matrix(4,4);
    MI = double_matrix(4,4);
    
    for (i=0; i< NP; i++) pfilt[i] = params[i];
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
       // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    pfilt[10] = 0.0;
    
    pfilt[4] = 0.0;
    pfilt[9] = 0.0;
    ResponseFreq(ll, pfilt, N, filtA[0], filtE[0], Tobs, Tzero);
    
    pfilt[4] = PI/2.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(ll, pfilt, N, filtA[1], filtE[1], Tobs, Tzero);
    
    pfilt[4] = 3.0*PI/4.0;
    pfilt[9] = 0.0;
    ResponseFreq(ll, pfilt, N, filtA[2], filtE[2], Tobs, Tzero);
    
    pfilt[4] = PI/4.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(ll, pfilt, N, filtA[3], filtE[3], Tobs, Tzero);
    
   //printf("\n");
        for (i=0; i< 4; i++)
        {
            KV[i] = 2.0*(fourier_nwip(filtA[i], AC, SN, N) + fourier_nwip(filtE[i], EC, SN, N))/Tobs;
            KP[i] = 2.0*(fourier_nwip_Tshift(delT, Tobs, filtA[i], AC, SN, N) + fourier_nwip_Tshift(delT, Tobs, filtE[i], EC, SN, N))/Tobs;
            KM[i] = 2.0*(fourier_nwip_Tshift(-delT, Tobs, filtA[i], AC, SN, N) + fourier_nwip_Tshift(-delT, Tobs, filtE[i], EC, SN, N))/Tobs;
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
    
    
    T0 = 0.0;
    T1 = 0.0;
    T2 = 0.0;
    
    for (i=0; i< 4; i++)
    {
        for (j=0; j< 4; j++)
        {
            T0 += 0.5*MI[i][j]*KV[i]*KV[j];
            T1 += 0.5*MI[i][j]*KP[i]*KP[j];
            T2 += 0.5*MI[i][j]*KM[i]*KM[j];
        }
    }
    
    c = T0;
    b = (T1-T2)/(2.0*delT);
    a = (T1-b*delT-c)/(delT*delT);
    
    dtx = -b/(2.0*a);
                        
   // printf("%f\n", dtx);
    
    // This method is not well define for shifts of more than a few delT
    if(dtx < -2.0*delT) dtx = -2.0*delT;
    if(dtx > 2.0*delT) dtx = 2.0*delT;
    
    params[5] += dtx;
    
    for (i=0; i< 4; i++)
    {
        KV[i] = 2.0*(fourier_nwip_Tshift(dtx, Tobs, filtA[i], AC, SN, N) + fourier_nwip_Tshift(dtx, Tobs, filtE[i], EC, SN, N))/Tobs;
    }
        
        logL = 0.0;
        
        for (i=0; i< 4; i++)
        {
            aV[i] = 0.0;
            for (j=0; j< 4; j++)
            {
                aV[i] += MI[i][j]*KV[j];
                logL += 0.5*MI[i][j]*KV[i]*KV[j];
            }
        }
        
       //printf("%f %f %f %f\n", logL, T0, T1, T2);
        
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

      if(ll == 0)
      {
        params[6] *= scale;
      }
      else
      {
       params[6] += log(scale);
      }
    
       params[10] = cosi;
       params[9] = psi;
       params[4] = phic;

    /*   Deallocate Arrays   */
    
     free(KV);
     free(KP);
     free(KM);
     free(aV);
     free_double_matrix(MM,4);
     free_double_matrix(MI,4);
     free_double_matrix(filtA,4);
     free_double_matrix(filtE,4);
    
    
    return logL;
}


void ResponseFreq(int ll, double *params, long N, double *AS, double *ES, double Tobs, double Tzero)
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
    
    double Amp, Phase, fonfs, f;
    
    double Aprime, Pprime, fi;
    
    double HC, HS, hp, hc;
    
    double m1, m2, chi1, chi2, phic, tc, distance, Mtot, eta, fr, af;
    
    double *ta, *xia, *FF;
    
    double Fp, Fc, kdotx, delt, fmin, fmax, dt;
    
    double XR, XI, YR, YI, ZR, ZI;
    
    int nfmin, nfmax, nf;
    
    int NA;
    
    double fstart, fstop;
    
 
    double *FpAR, *FpAI, *FcAR, *FcAI;
    double *FpER, *FpEI, *FcER, *FcEI;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    
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
    
    StartStop(ll, params, Tobs, Tzero, Tzero+Tobs, &fstart, &fstop, &fr);
    
    //printf("%e %e\n", fstart, fstop);
    
    nfmin = (int)(fstart*Tobs);
    nfmax = (int)(fstop*Tobs);
    if(nfmax > N/2) nfmax = N/2;
    nf = nfmax-nfmin;
    
    fmin = (double)(nfmin)/Tobs;
    fmax = (double)(nfmax)/Tobs;
    
    deltaF = 1.0/Tobs;
    
    //printf("%e %e %d %ld\n", fmin, fmax, nf, N/2);
    
    
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
    
    // start = clock();
    
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
    
    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("PhenomD call took %f seconds\n", cpu_time_used);
    
    // start = clock();
    
    
    // compute the Frequency time series
    timearray(params, freq, nf, TF, ap);
    
    RAantenna(params, nf, TF, FF, xia, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
    /*
     out = fopen("check.dat","w");
     for (i=0; i< nf; i++)
     {
     fprintf(out,"%d %e %e %e\n", i, freq->data[i], TF[i], ap->amp[i]);
     }
     fclose(out);
     */
    
    
    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    // printf("time array call took %f seconds\n", cpu_time_used);
    
    // start = clock();
    
    for (i=0; i< nf; i++)
    {
        AF[i] =  h22fac*ap->amp[i];
        fonfs = freq->data[i]/fstar;
        AF[i] *= (8.0*fonfs*sin(fonfs));   // conversion to fractional frequency and leading order TDI transfer function
    }
    
    
    for(n=0; n< N; n++)
    {
        AS[n] = 0.0;
        ES[n] = 0.0;
    }
    
    /*   Main Loop   */
    
    // start = clock();
    
    nn = 0;
    for(n=nfmin; n< nfmax; n++)
    {
        // Barycenter time and frequency
        
        // The frequency and time arrays start at nfmin
        m = n-nfmin;
        
        if(m > -1 && m < nf)
        {
            t = TF[m];
            f = FF[m];
            xi = xia[m];
            
            kdotx = t-xi;
            
            Amp = AF[m];
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
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
    return;
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



void pbt_shift(double *corr, double *corrf, double *data1, double *data2, double *Sn, int n)
{
    int nb2, i, l, k, j;
    int imax, imin;
    
    nb2 = n/2;
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        
        corr[l]	= (data1[l]*data2[l] + data1[k]*data2[k])/Sn[i];
        corr[k]	= (data1[k]*data2[l] - data1[l]*data2[k])/Sn[i];
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


void getfreq(double *fnew, double *tf, double t, double fguess, double m1_SI, double m2_SI, double chi1, double chi2, double tc)
{
    AmpPhaseFDWaveform *ap = NULL;
    double ep, u, v, tnew, x;
    double delT, delF, dtdf, fonfs;
    int ret;
    double M_sec;
    RealVector *freq;
    
    M_sec = (m1_SI+m2_SI) * MTSUN_SI/MSUN_SI;
    
    ep = 1.0e-6/M_sec;
    v = (4.0*PI*ep);
    u = (2.0*PI*ep*ep);
    
    freq = CreateRealVector((3));
    
    freq->data[0] = fguess-ep;
    freq->data[1] = fguess;
    freq->data[2] = fguess+ep;
    
    if(freq->data[0] < 0.0)
    {
        freq->data[0] = fguess;
        freq->data[1] = fguess+ep;
        freq->data[2] = fguess+2.0*ep;
    }
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, 0.0, PDfref, m1_SI, m2_SI, chi1, chi2, PC_SI);
    
    tnew = (ap->phase[2]-ap->phase[0])/v +tc;
    
    dtdf = (ap->phase[2]+ap->phase[0]-2.0*ap->phase[1])/u;
    
    delT = t-tnew;
    
    delF = delT/dtdf;
    
    *fnew = fguess + delF;
    
    *tf = tnew;
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
}



void bwlf(double *in, double *out, int fwrv, int M, int n, double s, double f)
{
    /* Butterworth bandpass filter
     n = filter order 2,4,6,8,...
     s = sampling frequency
     f = half power frequency
     */
    
    if(n % 2){ printf("Order must be 2,4,6,8,...\n"); return;}
    
    int i, j;
    n = n/2;
    double a = tan(PI*f/s);
    double a2 = a*a;
    double r;
    double *A = (double *)malloc(n*sizeof(double));
    double *d1 = (double *)malloc(n*sizeof(double));
    double *d2 = (double *)malloc(n*sizeof(double));
    double *w0 = (double *)calloc(n, sizeof(double));
    double *w1 = (double *)calloc(n, sizeof(double));
    double *w2 = (double *)calloc(n, sizeof(double));
    double x;
    
    for(i=0; i<n; ++i)
    {
        r = sin(PI*(2.0*i+1.0)/(4.0*n));
        s = a2 + 2.0*a*r + 1.0;
        A[i] = a2/s;
        d1[i] = 2.0*(1-a2)/s;
        d2[i] = -(a2 - 2.0*a*r + 1.0)/s;
        w0[i] = 0.0;
        w1[i] = 0.0;
        w2[i] = 0.0;
    }
    
    for(j=0; j< M; ++j)
    {
        if(fwrv == 1) x = in[j];
        if(fwrv == -1) x = in[M-j-1];
        for(i=0; i<n; ++i)
        {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i] + x;
            x = A[i]*(w0[i] + 2.0*w1[i] + w2[i]);
            w2[i] = w1[i];
            w1[i] = w0[i];
        }
        if(fwrv == 1) out[j] = x;
        if(fwrv == -1) out[M-j-1] = x;
    }
    
    free(A);
    free(d1);
    free(d2);
    free(w0);
    free(w1);
    free(w2);
    
    return;
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

double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
    free(v);
}

int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
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

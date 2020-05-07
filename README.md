# LISA-Black-Hole
Codes to find and characterize massive black hole binary signals

The codes in this release can be used to analyze data from the Radler release of the LISA  Data Challenge. With some minor modifications the codes can also be used to study the parameter estimation performance of LISA on catalogs of massive black holes.

The first thing to do is download the Radler data: The file you need is here
https://lisa-ldc.lal.in2p3.fr/media/uploads/LDC1-1_MBHB_v2_TD_9gc2s16.hdf5

To see the parameters used to generate the data you can use some of the LDC python tools to query the hdf5 file:

python3 LISAh5_display.py LDC1-1_MBHB_v2_TD_9gc2s16.hdf5

The next step is to dump the data to ASCII so that our codes can read it in. The command is

python3 LISAh5_tdi2ascii.py LDC1-1_MBHB_v2_TD_9gc2s16.hdf5 LDC1-1_MBHB_v2_TD.txt

The next code to run is segementAET.c, which is compiled using 

gcc -o segmentAET segmentAET.c -lm -lgsl

and run by typing

./segmentAET

This code reads in the time domain data, forms up the A,E,T TDI combinations, applies a Tukey window, then Fourier transforms the data. The segmentAET code outputs the full data set, in addition to transforms of the data broken up into 16 time chunks. These chunks are used in the search phase.

Next comes the search. The search code is compiled with 

OSX
clang -Xpreprocessor -fopenmp -lomp -w -o search search.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

Linux
gcc -std=gnu99 -fopenmp -w -o search search.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

The search can be run on any of the 16 data segments. The Radler source merges in segment 9 (the count starts at zero), so the most interesting case to run is

./search 9

As the code runs it prints information to the screen. During the initial search the output is

iteration# logL t_c m_1 m_2 swap_acceptance uniform_acceptance fisher_acceptance normal_acceptance

After the initial search, the code moves on to find the sky location. During this stage the screen output is

iteration# logL t_c cos_theta phi swap_acceptance uniform_acceptance time-sky_acceptance sky_acceptance

Once both stages of the search are complete (about 10 minutes on a laptop), the search code outputs a file named skymax.dat containing the maximum likelihood solution for all the source parameters. This file serves as the starting point for the parameter estimation code.

The parameter estimation code is compiled via

OSX
clang -Xpreprocessor -fopenmp -lomp -w -o  PTMCMC PTMCMC.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

Linux
gcc -std=gnu99 -fopenmp -w -o PTMCMC PTMCMC.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

To run the PTMCMC type

./PTMCMC

Options for the analysis are currently set in the header file Constants.h. A more flexible solution would be to pass arguments at the command line. The main parameters you might want to change are the nflag and cut flags. When nflag is set equal to 0 the noise free likelihood is used. This is useful for running on catalogs. Setting nflag to 1 causes the full likelihood to be used. The cut flag selects between using the full Radler data set (cut = 0) and the truncated data set (cut = 1) with the data truncated in time prior to merger and ringdown.

As the PTMCMC codes runs it prints a summary to screen. The columns are:

iteration# logL m_1 m_2 swap_acceptance fisher_acceptance DE_acceptance sky_acceptance

The default setting is for the PTMCMC code to use NC=16 parallel chains and to perform 10^6 iterations. So 1.6e7 likelihood evaluations. There is also additional overhead in updating the Fisher matrices, computing the more expensive F-statistic etc. The total run time on a 2.9 GHz quad-core laptop is about 5 hours. 

